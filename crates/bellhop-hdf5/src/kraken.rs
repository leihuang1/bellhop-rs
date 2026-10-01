//! Sequential legacy KRAKEN/KRAKENC runs with atomic HDF5 schema-v1 output.
//! This schema is independent of BELLHOP schema v3; no solver result types are unified.

use std::fmt;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use hdf5::{File, Group, H5Type};
use kraken::{Case, DiagnosticReport, ModeSolver, SimulationResult, SourceGeometry};
use sha2::{Digest, Sha256};

use super::{hdf5_error, write_scalar_attribute, write_string_attribute};

pub const SCHEMA_VERSION: u32 = 1;
pub const DEFAULT_MAX_OUTPUT_BYTES: u64 = 256 * 1024 * 1024;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RunSummary {
    pub frequency_count: usize,
    pub mode_count: u64,
    pub pressure_count: u64,
}

#[derive(Debug)]
pub enum RunError {
    Input(String),
    Simulation {
        frequency_index: usize,
        frequency_hz: f64,
        report: DiagnosticReport,
    },
    Output(String),
}

impl fmt::Display for RunError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input(message) => formatter.write_str(message),
            Self::Simulation {
                frequency_index,
                frequency_hz,
                report,
            } => {
                write!(
                    formatter,
                    "frequency[{frequency_index}] ({frequency_hz} Hz): {report}"
                )
            }
            Self::Output(message) => write!(formatter, "error[KR0402]: {message}"),
        }
    }
}

impl std::error::Error for RunError {}

/// Solve one supported legacy pair and publish a complete, versioned HDF5 result.
///
/// Exact input snapshots are parsed and hashed once. Frequencies are solved/written
/// sequentially; no collection of all frequency results is retained. Dataset payload
/// and flushed file size are bounded by `max_output_bytes` (default 256 MiB in the CLI).
/// Existing per-frequency numerical work limits remain unchanged.
///
/// Publication is atomic: without `overwrite`, a hard link refuses an existing
/// destination, including one created during the run. With it, rename replaces the
/// old output only after every frequency succeeds. Inputs cannot be destinations.
/// A pre-existing `<output>.tmp` is never truncated or removed.
///
/// # Errors
///
/// Returns input diagnostics, numerical diagnostics, or a quota/filesystem/HDF5 error.
/// No incomplete result is published on failure.
#[allow(clippy::too_many_lines)]
pub fn run_legacy(
    env_path: &Path,
    flp_path: &Path,
    output_path: &Path,
    solver: ModeSolver,
    overwrite: bool,
    max_output_bytes: u64,
) -> Result<RunSummary, RunError> {
    if max_output_bytes == 0 {
        return Err(RunError::Output("max_output_bytes must be positive".into()));
    }
    let env_source = read_source(env_path)?;
    let flp_source = read_source(flp_path)?;
    let bottom_table = kraken::legacy::bottom_table_extension(&env_source, env_path, solver)
        .map_err(|report| RunError::Input(report.to_string()))?
        .map(|extension| {
            let path = env_path.with_extension(extension);
            let source = read_source(&path)?;
            Ok::<_, RunError>((extension, path, source))
        })
        .transpose()?;
    let surface_table = kraken::legacy::surface_table_extension(&env_source, env_path, solver)
        .map_err(|report| RunError::Input(report.to_string()))?
        .map(|extension| {
            let path = env_path.with_extension(extension);
            let source = read_source(&path)?;
            Ok::<_, RunError>((extension, path, source))
        })
        .transpose()?;
    let cases = kraken::legacy::load_frequency_cases_with_boundary_tables(
        &env_source,
        &flp_source,
        env_path,
        flp_path,
        solver,
        surface_table.as_ref().map(|(_, _, source)| source.as_str()),
        bottom_table.as_ref().map(|(_, _, source)| source.as_str()),
    )
    .map_err(|report| RunError::Input(report.to_string()))?;
    let mut inputs = vec![
        ("env", env_path, env_source.as_str()),
        ("flp", flp_path, flp_source.as_str()),
    ];
    for table in [&surface_table, &bottom_table].into_iter().flatten() {
        inputs.push((table.0, table.1.as_path(), table.2.as_str()));
    }
    protect_inputs(
        output_path,
        &inputs.iter().map(|(_, path, _)| *path).collect::<Vec<_>>(),
    )
    .map_err(RunError::Output)?;
    match fs::symlink_metadata(output_path) {
        Ok(_) if !overwrite => {
            return Err(RunError::Output(format!(
                "output already exists: {}; pass --overwrite to replace it",
                output_path.display()
            )));
        }
        Err(error) if error.kind() != std::io::ErrorKind::NotFound => {
            return Err(RunError::Output(error.to_string()));
        }
        _ => {}
    }
    // Admit the known FIELD payload before doing numerical work. Modes are charged
    // as they are written, because their counts are not known before solving.
    let pressure_bytes = cases.iter().try_fold(0_u64, |total, case| {
        total.checked_add(
            (case.source_depths_m.len() as u64)
                .checked_mul(case.receiver_depths_m.len() as u64)?
                .checked_mul(case.receiver_ranges_m.len() as u64)?
                .checked_mul(8)?,
        )
    });
    if pressure_bytes.is_none_or(|bytes| bytes > max_output_bytes) {
        return Err(RunError::Output(format!(
            "FIELD payload exceeds {max_output_bytes} bytes"
        )));
    }

    let mut temporary = output_path.as_os_str().to_os_string();
    temporary.push(".tmp");
    let temporary = PathBuf::from(temporary);
    // HDF5's exclusive create is the reservation; do not check then truncate.
    let file =
        File::create_excl(&temporary).map_err(|error| RunError::Output(hdf5_error(error)))?;
    let cleanup = TemporaryOutput(temporary);
    (|| {
        // Own the handle in this inner scope: every early return drops it (and
        // child groups) before the outer guard removes the scratch path.
        let file = file;
        let mut budget = Budget {
            payload: 0,
            maximum: max_output_bytes,
        };
        write_header(&file, &cases, &inputs, &mut budget).map_err(RunError::Output)?;
        check_file(&file, &cleanup.0, max_output_bytes).map_err(RunError::Output)?;
        let frequencies = file
            .create_group("frequencies")
            .map_err(|e| RunError::Output(hdf5_error(e)))?;
        let mut summary = RunSummary {
            frequency_count: cases.len(),
            ..RunSummary::default()
        };
        for (index, case) in cases.iter().enumerate() {
            let result = kraken::solve(case).map_err(|report| RunError::Simulation {
                frequency_index: index,
                frequency_hz: case.frequency_hz,
                report,
            })?;
            let group = frequencies
                .create_group(&index.to_string())
                .map_err(|e| RunError::Output(hdf5_error(e)))?;
            write_frequency(&group, case, &result, &mut budget).map_err(RunError::Output)?;
            summary.mode_count += result.modes.modes.len() as u64;
            summary.pressure_count += result.field.pressure.len() as u64;
            check_file(&file, &cleanup.0, max_output_bytes).map_err(RunError::Output)?;
        }
        drop(frequencies);
        write_scalar_attribute(&file, "mode_count", &summary.mode_count)
            .map_err(RunError::Output)?;
        write_scalar_attribute(&file, "pressure_count", &summary.pressure_count)
            .map_err(RunError::Output)?;
        check_file(&file, &cleanup.0, max_output_bytes).map_err(RunError::Output)?;
        file.close().map_err(|e| RunError::Output(hdf5_error(e)))?;
        check_file_size(&cleanup.0, max_output_bytes).map_err(RunError::Output)?;
        fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&cleanup.0)
            .and_then(|file| file.sync_all())
            .map_err(|e| RunError::Output(e.to_string()))?;
        if overwrite {
            fs::rename(&cleanup.0, output_path)
        } else {
            fs::hard_link(&cleanup.0, output_path)
        }
        .map_err(|e| {
            RunError::Output(format!("unable to install {}: {e}", output_path.display()))
        })?;
        Ok(summary)
    })()
}

fn read_source(path: &Path) -> Result<String, RunError> {
    let mut source = String::new();
    fs::File::open(path)
        .and_then(|file| {
            file.take(kraken::legacy::MAX_INPUT_BYTES + 1)
                .read_to_string(&mut source)
        })
        .map_err(|e| {
            RunError::Input(format!(
                "{}:1:1: error[KR0001]: unable to read input: {e} (input)",
                path.display()
            ))
        })?;
    Ok(source)
}

fn protect_inputs(output: &Path, inputs: &[&Path]) -> Result<(), String> {
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let filename = output.file_name().ok_or("output must name a file")?;
    let destination = fs::canonicalize(parent)
        .map_err(|e| e.to_string())?
        .join(filename);
    let existing = fs::canonicalize(output).ok();
    for input in inputs {
        let input = fs::canonicalize(input).map_err(|e| e.to_string())?;
        if destination == input || existing.as_ref() == Some(&input) {
            return Err("output must not replace an input file or its symlink alias".into());
        }
    }
    Ok(())
}

struct TemporaryOutput(PathBuf);

impl Drop for TemporaryOutput {
    fn drop(&mut self) {
        if let Err(error) = fs::remove_file(&self.0)
            && error.kind() != std::io::ErrorKind::NotFound
        {
            eprintln!(
                "warning[KR0402]: unable to remove scratch {}: {error}",
                self.0.display()
            );
        }
    }
}

struct Budget {
    payload: u64,
    maximum: u64,
}

impl Budget {
    fn reserve(&mut self, bytes: u64) -> Result<(), String> {
        self.payload = self
            .payload
            .checked_add(bytes)
            .filter(|&n| n <= self.maximum)
            .ok_or_else(|| format!("dataset payload exceeds {} bytes", self.maximum))?;
        Ok(())
    }
}

fn check_file(file: &File, path: &Path, maximum: u64) -> Result<(), String> {
    // ponytail: physical quota checked after header/frequency flushes; a bounded HDF5 VFD
    // is needed for a strict in-write disk quota, not for sequential bounded results.
    file.flush().map_err(hdf5_error)?;
    check_file_size(path, maximum)
}

fn check_file_size(path: &Path, maximum: u64) -> Result<(), String> {
    if fs::metadata(path).map_err(|e| e.to_string())?.len() > maximum {
        return Err(format!("HDF5 file exceeds {maximum} bytes"));
    }
    Ok(())
}

fn write_header(
    file: &File,
    cases: &[Case],
    inputs: &[(&str, &Path, &str)],
    budget: &mut Budget,
) -> Result<(), String> {
    write_scalar_attribute(file, "schema_version", &SCHEMA_VERSION)?;
    write_string_attribute(file, "schema_name", "kraken")?;
    write_string_attribute(
        file,
        "implementation",
        concat!("bellhop-rs kraken ", env!("CARGO_PKG_VERSION")),
    )?;
    write_string_attribute(
        file,
        "compatibility_reference",
        "Acoustics Toolbox v2023.5 (475108519289c6fb488b58980c644ea14eccc604)",
    )?;
    write_string_attribute(file, "title", &cases[0].title)?;
    write_string_attribute(
        file,
        "solver",
        match cases[0].mode_solver {
            ModeSolver::Kraken => "kraken",
            ModeSolver::Krakenc => "krakenc",
        },
    )?;
    write_string_attribute(
        file,
        "coordinate_convention",
        "range origin at source; depth positive downward",
    )?;
    write_scalar_attribute(file, "frequency_count", &(cases.len() as u64))?;
    write_scalar_attribute(file, "max_output_bytes", &budget.maximum)?;
    let frequencies: Vec<_> = cases.iter().map(|case| case.frequency_hz).collect();
    dataset(
        file,
        "frequency_hz",
        &frequencies,
        &[frequencies.len()],
        "Hz",
        budget,
    )?;
    let group = file.create_group("inputs").map_err(hdf5_error)?;
    for (name, path, source) in inputs {
        let input = group.create_group(name).map_err(hdf5_error)?;
        write_string_attribute(&input, "filename", &path.to_string_lossy())?;
        write_scalar_attribute(&input, "size_bytes", &(source.len() as u64))?;
        write_string_attribute(
            &input,
            "sha256",
            &format!("{:x}", Sha256::digest(source.as_bytes())),
        )?;
    }
    Ok(())
}

fn write_frequency(
    group: &Group,
    case: &Case,
    result: &SimulationResult,
    budget: &mut Budget,
) -> Result<(), String> {
    write_scalar_attribute(group, "frequency_hz", &result.modes.frequency_hz)?;
    write_string_attribute(
        group,
        "surface_boundary",
        match &case.surface_boundary {
            kraken::SurfaceBoundary::Vacuum => "V",
            kraken::SurfaceBoundary::Rigid => "R",
            kraken::SurfaceBoundary::FluidHalfSpace => "A",
            kraken::SurfaceBoundary::Reflection(_) => "F",
            kraken::SurfaceBoundary::Impedance { .. } => "P",
        },
    )?;
    write_string_attribute(
        group,
        "bottom_boundary",
        match &case.bottom_boundary {
            kraken::BottomBoundary::Vacuum => "V",
            kraken::BottomBoundary::FluidHalfSpace => "A",
            kraken::BottomBoundary::Rigid => "R",
            kraken::BottomBoundary::Reflection(_) => "F",
            kraken::BottomBoundary::Impedance { .. } => "P",
        },
    )?;
    write_scalar_attribute(
        group,
        "mesh_reference_frequency_hz",
        &case
            .mesh_reference_frequency_hz
            .unwrap_or(case.frequency_hz),
    )?;
    write_scalar_attribute(group, "requested_mesh_points", &(case.mesh_points as u64))?;
    write_scalar_attribute(group, "max_range_m", &case.max_range_m)?;
    write_scalar_attribute(group, "field_mode_limit", &(case.mode_limit as u64))?;
    write_string_attribute(
        group,
        "source_geometry",
        match case.source_geometry {
            SourceGeometry::Line => "line",
            SourceGeometry::Point => "point",
        },
    )?;

    write_scalar_attribute(
        group,
        "finite_fluid_layer_count",
        &((1 + case.additional_fluid_layers.len()) as u64),
    )?;
    let media = group.create_group("media").map_err(hdf5_error)?;
    let mut top_depth = 0.0;
    for index in 0..=case.additional_fluid_layers.len() {
        let (bottom_depth, density, mesh_points) = if index == 0 {
            (
                case.water_depth_m,
                case.water_density_g_cm3,
                case.mesh_points,
            )
        } else {
            let layer = &case.additional_fluid_layers[index - 1];
            (layer.bottom_depth_m, layer.density_g_cm3, layer.mesh_points)
        };
        let layer = media.create_group(&index.to_string()).map_err(hdf5_error)?;
        write_scalar_attribute(&layer, "top_depth_m", &top_depth)?;
        write_scalar_attribute(&layer, "bottom_depth_m", &bottom_depth)?;
        write_scalar_attribute(&layer, "density_g_cm3", &density)?;
        write_scalar_attribute(&layer, "requested_mesh_points", &(mesh_points as u64))?;
        top_depth = bottom_depth;
    }

    write_modes(group, result, budget)?;
    write_field(group, result, budget)
}

fn write_modes(
    group: &Group,
    result: &SimulationResult,
    budget: &mut Budget,
) -> Result<(), String> {
    let modes = group.create_group("modes").map_err(hdf5_error)?;
    let count = result.modes.modes.len();
    let depths = &result.modes.sampled_depths_m;
    write_string_attribute(&modes, "eigenfunction_axis_order", "mode,sample_depth")?;
    write_string_attribute(
        &modes,
        "normalization",
        "pinned AT density-weighted normalization; fluid half-space contribution when present; arbitrary unit phase",
    )?;
    dataset(
        &modes,
        "sample_depth_m",
        depths,
        &[depths.len()],
        "m",
        budget,
    )?;
    for (name, unit, values) in [
        (
            "horizontal_wavenumber_real",
            "rad/m",
            result
                .modes
                .modes
                .iter()
                .map(|m| m.horizontal_wavenumber_rad_per_m.re)
                .collect::<Vec<_>>(),
        ),
        (
            "horizontal_wavenumber_imaginary",
            "rad/m",
            result
                .modes
                .modes
                .iter()
                .map(|m| m.horizontal_wavenumber_rad_per_m.im)
                .collect(),
        ),
        (
            "phase_speed_mps",
            "m/s",
            result
                .modes
                .modes
                .iter()
                .map(|m| m.phase_speed_mps)
                .collect(),
        ),
        (
            "group_speed_mps",
            "m/s",
            result
                .modes
                .modes
                .iter()
                .map(|m| m.group_speed_mps)
                .collect(),
        ),
        (
            "attenuation_nepers_per_m",
            "neper/m",
            result
                .modes
                .modes
                .iter()
                .map(|m| m.attenuation_nepers_per_m)
                .collect(),
        ),
    ] {
        dataset(&modes, name, &values, &[count], unit, budget)?;
    }
    for imaginary in [false, true] {
        let values: Vec<_> = result
            .modes
            .modes
            .iter()
            .flat_map(|m| &m.eigenfunction)
            .map(|value| if imaginary { value.im } else { value.re })
            .collect();
        dataset(
            &modes,
            if imaginary {
                "eigenfunction_imaginary"
            } else {
                "eigenfunction_real"
            },
            &values,
            &[count, depths.len()],
            "reference_normalized",
            budget,
        )?;
    }
    Ok(())
}

#[allow(clippy::cast_possible_truncation)]
fn write_field(
    group: &Group,
    result: &SimulationResult,
    budget: &mut Budget,
) -> Result<(), String> {
    let field = group.create_group("field").map_err(hdf5_error)?;
    write_string_attribute(
        &field,
        "pressure_axis_order",
        "source_depth,receiver_depth,receiver_range",
    )?;
    for (name, values) in [
        ("source_depth_m", &result.field.source_depths_m),
        ("receiver_depth_m", &result.field.receiver_depths_m),
        ("receiver_range_m", &result.field.receiver_ranges_m),
        ("receiver_offset_m", &result.field.receiver_offsets_m),
    ] {
        dataset(&field, name, values, &[values.len()], "m", budget)?;
    }
    let shape = [
        result.field.source_depths_m.len(),
        result.field.receiver_depths_m.len(),
        result.field.receiver_ranges_m.len(),
    ];
    for imaginary in [false, true] {
        let values: Vec<_> = result
            .field
            .pressure
            .iter()
            .map(|value| (if imaginary { value.im } else { value.re }) as f32)
            .collect();
        dataset(
            &field,
            if imaginary {
                "pressure_imaginary"
            } else {
                "pressure_real"
            },
            &values,
            &shape,
            "1",
            budget,
        )?;
    }
    Ok(())
}

fn dataset<T: H5Type + Copy + Into<f64>>(
    group: &Group,
    name: &str,
    values: &[T],
    shape: &[usize],
    unit: &str,
    budget: &mut Budget,
) -> Result<(), String> {
    if values.iter().any(|&value| !value.into().is_finite()) {
        return Err(format!("non-finite output in {name}"));
    }
    budget.reserve(
        (values.len() as u64)
            .checked_mul(std::mem::size_of::<T>() as u64)
            .ok_or("dataset byte count overflow")?,
    )?;
    let data = group
        .new_dataset::<T>()
        .shape(shape)
        .create(name)
        .map_err(hdf5_error)?;
    data.write_raw(values).map_err(hdf5_error)?;
    write_string_attribute(&data, "unit", unit)
}
