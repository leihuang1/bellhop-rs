use std::fs;
use std::path::Path;
use std::str::FromStr;

use bellhop::solver::{RayTermination, SimulationResult};
use hdf5::types::VarLenUnicode;
use hdf5::{File, Group, H5Type};
use sha2::{Digest, Sha256};

const SCHEMA_VERSION: u32 = 1;

pub(crate) fn write_hdf5(
    path: &Path,
    input_path: &Path,
    result: &SimulationResult,
) -> Result<(), String> {
    let input =
        fs::read(input_path).map_err(|error| format!("unable to read input metadata: {error}"))?;
    let input_size =
        u64::try_from(input.len()).map_err(|_| "input file size does not fit in u64".to_owned())?;
    let input_sha256 = format!("{:x}", Sha256::digest(&input));
    let input_filename = input_path
        .file_name()
        .map_or_else(String::new, |name| name.to_string_lossy().into_owned());

    let file = File::create(path).map_err(hdf5_error)?;
    write_scalar_attribute(&file, "schema_version", &SCHEMA_VERSION)?;
    write_string_attribute(
        &file,
        "implementation",
        concat!("bellhop-rs ", env!("CARGO_PKG_VERSION")),
    )?;
    write_string_attribute(
        &file,
        "compatibility_reference",
        "Acoustics Toolbox v2023.5 (475108519289c6fb488b58980c644ea14eccc604)",
    )?;
    write_string_attribute(&file, "input_filename", &input_filename)?;
    write_scalar_attribute(&file, "input_size_bytes", &input_size)?;
    write_string_attribute(&file, "input_sha256", &input_sha256)?;
    write_string_attribute(
        &file,
        "coordinate_convention",
        "range origin at source; depth positive downward; launch angle positive downward",
    )?;

    let rays = file.create_group("rays").map_err(hdf5_error)?;
    write_rays(&rays, result)?;
    file.flush().map_err(hdf5_error)
}

#[allow(clippy::too_many_lines)]
fn write_rays(group: &Group, result: &SimulationResult) -> Result<(), String> {
    let source_depth_m: Vec<f32> = result
        .sources
        .iter()
        .map(|source| source.source_depth_m)
        .collect();
    let ray_count: usize = result.sources.iter().map(|source| source.rays.len()).sum();
    let point_count: usize = result
        .sources
        .iter()
        .flat_map(|source| &source.rays)
        .map(|ray| ray.points.len())
        .sum();
    let mut ray_source_index = Vec::with_capacity(ray_count);
    let mut launch_angle_degrees = Vec::with_capacity(ray_count);
    let mut point_offset = Vec::with_capacity(ray_count + 1);
    let mut top_bounces = Vec::with_capacity(ray_count);
    let mut bottom_bounces = Vec::with_capacity(ray_count);
    let mut termination = Vec::with_capacity(ray_count);
    let mut range_m = Vec::with_capacity(point_count);
    let mut depth_m = Vec::with_capacity(point_count);
    let mut travel_time_s = Vec::with_capacity(point_count);
    let mut attenuation_time_s = Vec::with_capacity(point_count);
    let mut amplitude = Vec::with_capacity(point_count);
    let mut phase_radians = Vec::with_capacity(point_count);
    point_offset.push(0_u64);

    for (source_index, source) in result.sources.iter().enumerate() {
        let source_index = u32::try_from(source_index)
            .map_err(|_| "source index does not fit in u32".to_owned())?;
        for ray in &source.rays {
            ray_source_index.push(source_index);
            launch_angle_degrees.push(ray.launch_angle_degrees);
            top_bounces.push(ray.top_bounces);
            bottom_bounces.push(ray.bottom_bounces);
            termination.push(termination_code(ray.termination));
            for point in &ray.points {
                range_m.push(point.range_m);
                depth_m.push(point.depth_m);
                travel_time_s.push(point.travel_time_s);
                attenuation_time_s.push(point.attenuation_time_s);
                amplitude.push(point.amplitude);
                phase_radians.push(point.phase_radians);
            }
            point_offset.push(
                u64::try_from(range_m.len())
                    .map_err(|_| "ray point offset does not fit in u64".to_owned())?,
            );
        }
    }

    write_dataset(group, "source_depth_m", &source_depth_m, "m")?;
    write_dataset(group, "ray_source_index", &ray_source_index, "index")?;
    write_dataset(
        group,
        "launch_angle_degrees",
        &launch_angle_degrees,
        "degree",
    )?;
    write_dataset(group, "point_offset", &point_offset, "index")?;
    write_dataset(group, "top_bounces", &top_bounces, "count")?;
    write_dataset(group, "bottom_bounces", &bottom_bounces, "count")?;
    write_dataset(group, "termination", &termination, "enum")?;
    write_dataset(group, "range_m", &range_m, "m")?;
    write_dataset(group, "depth_m", &depth_m, "m")?;
    write_dataset(group, "travel_time_s", &travel_time_s, "s")?;
    write_dataset(group, "attenuation_time_s", &attenuation_time_s, "s")?;
    write_dataset(group, "amplitude", &amplitude, "1")?;
    write_dataset(group, "phase_radians", &phase_radians, "rad")?;
    write_string_attribute(
        group,
        "termination_codes",
        "0=exited_trace_box,1=lost_energy,2=escaped_boundary,3=source_outside_boundaries,4=step_limit",
    )?;
    Ok(())
}

fn termination_code(termination: RayTermination) -> u8 {
    match termination {
        RayTermination::ExitedTraceBox => 0,
        RayTermination::LostEnergy => 1,
        RayTermination::EscapedBoundary => 2,
        RayTermination::SourceOutsideBoundaries => 3,
        RayTermination::StepLimit => 4,
    }
}

fn write_dataset<T: H5Type>(
    group: &Group,
    name: &str,
    values: &[T],
    unit: &str,
) -> Result<(), String> {
    let dataset = group
        .new_dataset_builder()
        .with_data(values)
        .create(name)
        .map_err(hdf5_error)?;
    write_string_attribute(&dataset, "unit", unit)
}

fn write_scalar_attribute<T: H5Type>(
    parent: &hdf5::Location,
    name: &str,
    value: &T,
) -> Result<(), String> {
    parent
        .new_attr::<T>()
        .create(name)
        .and_then(|attribute| attribute.write_scalar(value))
        .map_err(hdf5_error)
}

fn write_string_attribute(parent: &hdf5::Location, name: &str, value: &str) -> Result<(), String> {
    let value = VarLenUnicode::from_str(value)
        .map_err(|error| format!("invalid metadata string: {error}"))?;
    parent
        .new_attr::<VarLenUnicode>()
        .create(name)
        .and_then(|attribute| attribute.write_scalar(&value))
        .map_err(hdf5_error)
}

#[allow(clippy::needless_pass_by_value)]
fn hdf5_error(error: hdf5::Error) -> String {
    format!("HDF5 output failed: {error}")
}

#[cfg(test)]
mod tests {
    use std::fs;

    use bellhop::solver::{
        RayPoint, RayTermination, RayTrajectory, SimulationResult, SourceRaySet,
    };
    use hdf5::File;

    use super::write_hdf5;

    #[test]
    fn writes_versioned_flattened_ray_schema() {
        let directory =
            std::env::temp_dir().join(format!("bellhop-hdf5-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        let input = directory.join("case.env");
        let output = directory.join("case.h5");
        fs::write(&input, "test input").unwrap();
        let result = SimulationResult {
            sources: vec![SourceRaySet {
                source_depth_m: 10.0,
                rays: vec![RayTrajectory {
                    launch_angle_degrees: 5.0,
                    points: vec![
                        RayPoint {
                            range_m: 0.0,
                            depth_m: 10.0,
                            travel_time_s: 0.0,
                            attenuation_time_s: 0.0,
                            amplitude: 1.0,
                            phase_radians: 0.0,
                        },
                        RayPoint {
                            range_m: 100.0,
                            depth_m: 20.0,
                            travel_time_s: 0.1,
                            attenuation_time_s: -0.01,
                            amplitude: 0.8,
                            phase_radians: 1.0,
                        },
                    ],
                    top_bounces: 0,
                    bottom_bounces: 1,
                    termination: RayTermination::ExitedTraceBox,
                }],
            }],
        };

        write_hdf5(&output, &input, &result).unwrap();
        let file = File::open(&output).unwrap();
        assert_eq!(
            file.attr("schema_version")
                .unwrap()
                .read_scalar::<u32>()
                .unwrap(),
            1
        );
        assert_eq!(
            file.dataset("rays/point_offset")
                .unwrap()
                .read_raw::<u64>()
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            file.dataset("rays/depth_m")
                .unwrap()
                .read_raw::<f64>()
                .unwrap(),
            vec![10.0, 20.0]
        );
        drop(file);
        fs::remove_dir_all(directory).unwrap();
    }
}
