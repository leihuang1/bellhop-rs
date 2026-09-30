// Serialization round-trips deliberately require exact floating-point equality.
#![allow(clippy::float_cmp)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use hdf5::types::VarLenUnicode;
use hdf5::{File, Group, H5Type};
use kraken::{Case, ModeSolver, SimulationResult};
use sha2::{Digest, Sha256};

fn directory(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("kraken-cli-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&path);
    fs::create_dir_all(&path).unwrap();
    path
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../kraken/tests/fixtures")
        .join(name)
}

fn run(env: &Path, output: &Path, solver: &str, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_kraken"))
        .arg("run")
        .arg(env)
        .args(["--solver", solver, "--output"])
        .arg(output)
        .args(extra)
        .output()
        .unwrap()
}

fn attribute(location: &hdf5::Location, name: &str) -> String {
    location
        .attr(name)
        .unwrap()
        .read_scalar::<VarLenUnicode>()
        .unwrap()
        .as_str()
        .to_owned()
}

fn data<T: H5Type>(group: &Group, name: &str, shape: &[usize], unit: &str) -> Vec<T> {
    let dataset = group.dataset(name).unwrap();
    assert_eq!(dataset.shape(), shape, "{name}");
    assert!(dataset.dtype().unwrap().is::<T>(), "{name} datatype");
    assert_eq!(attribute(&dataset, "unit"), unit, "{name}");
    dataset.read_raw().unwrap()
}

#[test]
fn cli_round_trips_single_and_multifrequency_results() {
    let root = directory("round-trip");
    for (name, solver) in [
        ("Pekeris", "kraken"),
        ("PekerisFiltered", "kraken"),
        ("PekerisHardBoth", "kraken"),
        ("MunkBottomLoss", "kraken"),
        ("PekerisComplexBlank", "krakenc"),
        ("PekerisComplexRefined", "krakenc"),
        ("MunkLeakyPartialC", "krakenc"),
        ("PekerisBroadband", "kraken"),
        ("PekerisComplexBroadband", "krakenc"),
    ] {
        let env = fixture(name).with_extension("env");
        let output = root.join(name).with_extension("h5");
        let process = run(&env, &output, solver, &[]);
        assert!(
            process.status.success(),
            "{}",
            String::from_utf8_lossy(&process.stderr)
        );
        assert_product(&output, &env, &env.with_extension("flp"), solver);
        assert!(!output.with_extension("h5.tmp").exists());
    }
    // Duplicate frequencies must have distinct index groups, never Hz-named keys.
    let env = root.join("duplicate.env");
    fs::write(
        &env,
        fs::read_to_string(fixture("PekerisBroadband.env"))
            .unwrap()
            .replace("75.0 50.0 62.5 /", "75.0 75.0 50.0 /"),
    )
    .unwrap();
    fs::copy(fixture("PekerisBroadband.flp"), env.with_extension("flp")).unwrap();
    let output = root.join("duplicate.h5");
    let process = run(&env, &output, "kraken", &[]);
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    assert_product(&output, &env, &env.with_extension("flp"), "kraken");
    fs::remove_dir_all(root).unwrap();
}

#[allow(clippy::too_many_lines)]
fn assert_product(output: &Path, env: &Path, flp: &Path, solver: &str) {
    let engine = if solver == "krakenc" {
        ModeSolver::Krakenc
    } else {
        ModeSolver::Kraken
    };
    let cases = kraken::legacy::load_frequency_cases(env, flp, engine).unwrap();
    let file = File::open(output).unwrap();
    assert_eq!(
        file.attr("schema_version")
            .unwrap()
            .read_scalar::<u32>()
            .unwrap(),
        1
    );
    assert_eq!(attribute(&file, "schema_name"), "kraken");
    assert_eq!(attribute(&file, "solver"), solver);
    assert_eq!(attribute(&file, "title"), cases[0].title);
    assert!(
        attribute(&file, "compatibility_reference")
            .contains("475108519289c6fb488b58980c644ea14eccc604")
    );
    assert_eq!(
        file.attr("frequency_count")
            .unwrap()
            .read_scalar::<u64>()
            .unwrap(),
        cases.len() as u64
    );
    assert_eq!(
        data::<f64>(&file, "frequency_hz", &[cases.len()], "Hz"),
        cases
            .iter()
            .map(|case| case.frequency_hz)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        file.group("frequencies")
            .unwrap()
            .member_names()
            .unwrap()
            .len(),
        cases.len()
    );
    for (name, path) in [("env", env), ("flp", flp)] {
        let bytes = fs::read(path).unwrap();
        let group = file.group(&format!("inputs/{name}")).unwrap();
        assert_eq!(attribute(&group, "filename"), path.to_string_lossy());
        assert_eq!(
            group
                .attr("size_bytes")
                .unwrap()
                .read_scalar::<u64>()
                .unwrap(),
            bytes.len() as u64
        );
        assert_eq!(
            attribute(&group, "sha256"),
            format!("{:x}", Sha256::digest(bytes))
        );
    }
    let mut mode_count = 0_u64;
    let mut pressure_count = 0_u64;
    for (index, case) in cases.iter().enumerate() {
        let expected = kraken::solve(case).unwrap();
        let group = file.group(&format!("frequencies/{index}")).unwrap();
        assert_eq!(
            group
                .attr("frequency_hz")
                .unwrap()
                .read_scalar::<f64>()
                .unwrap(),
            case.frequency_hz
        );
        assert_eq!(
            group
                .attr("mesh_reference_frequency_hz")
                .unwrap()
                .read_scalar::<f64>()
                .unwrap(),
            case.mesh_reference_frequency_hz
                .unwrap_or(case.frequency_hz)
        );
        assert_eq!(
            group
                .attr("requested_mesh_points")
                .unwrap()
                .read_scalar::<u64>()
                .unwrap(),
            case.mesh_points as u64
        );
        assert_eq!(
            group
                .attr("max_range_m")
                .unwrap()
                .read_scalar::<f64>()
                .unwrap(),
            case.max_range_m
        );
        assert_eq!(
            group
                .attr("field_mode_limit")
                .unwrap()
                .read_scalar::<u64>()
                .unwrap(),
            case.mode_limit as u64
        );
        assert_eq!(
            attribute(&group, "source_geometry"),
            match case.source_geometry {
                kraken::SourceGeometry::Line => "line",
                kraken::SourceGeometry::Point => "point",
            }
        );
        assert_modes(&group.group("modes").unwrap(), &expected);
        assert_field(&group.group("field").unwrap(), case, &expected);
        mode_count += expected.modes.modes.len() as u64;
        pressure_count += expected.field.pressure.len() as u64;
    }
    assert_eq!(
        file.attr("mode_count")
            .unwrap()
            .read_scalar::<u64>()
            .unwrap(),
        mode_count
    );
    assert_eq!(
        file.attr("pressure_count")
            .unwrap()
            .read_scalar::<u64>()
            .unwrap(),
        pressure_count
    );
}

fn assert_modes(group: &Group, expected: &SimulationResult) {
    let depths = &expected.modes.sampled_depths_m;
    let count = expected.modes.modes.len();
    assert_eq!(
        attribute(group, "eigenfunction_axis_order"),
        "mode,sample_depth"
    );
    assert_eq!(
        data::<f64>(group, "sample_depth_m", &[depths.len()], "m"),
        *depths
    );
    for (name, unit, values) in [
        (
            "horizontal_wavenumber_real",
            "rad/m",
            expected
                .modes
                .modes
                .iter()
                .map(|m| m.horizontal_wavenumber_rad_per_m.re)
                .collect::<Vec<_>>(),
        ),
        (
            "horizontal_wavenumber_imaginary",
            "rad/m",
            expected
                .modes
                .modes
                .iter()
                .map(|m| m.horizontal_wavenumber_rad_per_m.im)
                .collect(),
        ),
        (
            "phase_speed_mps",
            "m/s",
            expected
                .modes
                .modes
                .iter()
                .map(|m| m.phase_speed_mps)
                .collect(),
        ),
        (
            "group_speed_mps",
            "m/s",
            expected
                .modes
                .modes
                .iter()
                .map(|m| m.group_speed_mps)
                .collect(),
        ),
        (
            "attenuation_nepers_per_m",
            "neper/m",
            expected
                .modes
                .modes
                .iter()
                .map(|m| m.attenuation_nepers_per_m)
                .collect(),
        ),
    ] {
        assert_eq!(data::<f64>(group, name, &[count], unit), values, "{name}");
    }
    for (name, imaginary) in [
        ("eigenfunction_real", false),
        ("eigenfunction_imaginary", true),
    ] {
        let values: Vec<_> = expected
            .modes
            .modes
            .iter()
            .flat_map(|mode| &mode.eigenfunction)
            .map(|value| if imaginary { value.im } else { value.re })
            .collect();
        assert_eq!(
            data::<f64>(group, name, &[count, depths.len()], "reference_normalized"),
            values
        );
    }
}

fn assert_field(group: &Group, case: &Case, expected: &SimulationResult) {
    assert_eq!(
        attribute(group, "pressure_axis_order"),
        "source_depth,receiver_depth,receiver_range"
    );
    for (name, values) in [
        ("source_depth_m", &case.source_depths_m),
        ("receiver_depth_m", &case.receiver_depths_m),
        ("receiver_range_m", &case.receiver_ranges_m),
        ("receiver_offset_m", &case.receiver_offsets_m),
    ] {
        assert_eq!(data::<f64>(group, name, &[values.len()], "m"), *values);
    }
    let shape = [
        case.source_depths_m.len(),
        case.receiver_depths_m.len(),
        case.receiver_ranges_m.len(),
    ];
    for (name, imaginary) in [("pressure_real", false), ("pressure_imaginary", true)] {
        let values: Vec<_> = expected
            .field
            .pressure
            .iter()
            .map(|p| if imaginary { p.im } else { p.re })
            .collect();
        let stored: Vec<_> = data::<f32>(group, name, &shape, "1")
            .into_iter()
            .map(f64::from)
            .collect();
        assert_eq!(stored, values, "{name}");
    }
}

fn assert_failure(process: &Output, code: i32, output: &Path, old: &[u8]) {
    assert_eq!(
        process.status.code(),
        Some(code),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    assert_eq!(fs::read(output).unwrap(), old);
    assert!(!output.with_extension("h5.tmp").exists());
}

#[test]
fn cli_preserves_old_outputs_and_inputs_on_failure() {
    let root = directory("failure");
    let env = root.join("case.env");
    let flp = env.with_extension("flp");
    let source = fs::read_to_string(fixture("PekerisBroadband.env")).unwrap();
    fs::write(&env, &source).unwrap();
    fs::copy(fixture("PekerisBroadband.flp"), &flp).unwrap();
    let output = root.join("result.h5");
    let old = b"previous complete result";
    fs::write(&output, old).unwrap();
    assert_failure(&run(&env, &output, "kraken", &[]), 4, &output, old);
    assert_failure(
        &run(
            &env,
            &output,
            "kraken",
            &["--overwrite", "--max-output-bytes", "100"],
        ),
        4,
        &output,
        old,
    );
    assert_failure(
        &run(
            &env,
            &output,
            "kraken",
            &["--overwrite", "--max-output-bytes", "4096"],
        ),
        4,
        &output,
        old,
    );

    // First frequency succeeds; 5 Hz has no trapped mode. No partial publication.
    fs::write(&env, source.replace("75.0 50.0 62.5 /", "75.0 5.0 62.5 /")).unwrap();
    let process = run(&env, &output, "kraken", &["--overwrite"]);
    assert_failure(&process, 3, &output, old);
    assert!(String::from_utf8_lossy(&process.stderr).contains("frequency[1] (5 Hz)"));
    fs::write(&env, &source).unwrap();
    for input in [&env, &flp] {
        let bytes = fs::read(input).unwrap();
        let process = run(&env, input, "kraken", &["--overwrite"]);
        assert_eq!(process.status.code(), Some(4));
        assert_eq!(fs::read(input).unwrap(), bytes);
    }
    let temporary = output.with_extension("h5.tmp");
    fs::write(&temporary, b"someone else's scratch file").unwrap();
    let process = run(&env, &output, "kraken", &["--overwrite"]);
    assert_eq!(process.status.code(), Some(4));
    assert_eq!(
        fs::read(&temporary).unwrap(),
        b"someone else's scratch file"
    );
    assert_eq!(fs::read(&output).unwrap(), old);
    fs::remove_file(temporary).unwrap();

    // A quota just below the complete file size also catches cumulative/metadata growth.
    let process = run(&env, &output, "kraken", &["--overwrite"]);
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    let complete = fs::read(&output).unwrap();
    let quota = (complete.len() - 1).to_string();
    assert_failure(
        &run(
            &env,
            &output,
            "kraken",
            &["--overwrite", "--max-output-bytes", &quota],
        ),
        4,
        &output,
        &complete,
    );
    fs::write(&env, "unsupported input").unwrap();
    assert_failure(
        &run(&env, &output, "kraken", &["--overwrite"]),
        2,
        &output,
        &complete,
    );
    fs::write(&env, [0xff]).unwrap();
    assert_failure(
        &run(&env, &output, "kraken", &["--overwrite"]),
        2,
        &output,
        &complete,
    );
    fs::write(&env, vec![b' '; 1_048_577]).unwrap();
    assert_failure(
        &run(&env, &output, "kraken", &["--overwrite"]),
        2,
        &output,
        &complete,
    );
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn explicit_field_input_and_default_output_work() {
    let root = directory("explicit-flp");
    let env = fixture("PekerisBroadband.env");
    let flp = fixture("PekerisBroadband.flp");
    let process = Command::new(env!("CARGO_BIN_EXE_kraken"))
        .current_dir(&root)
        .arg("run")
        .arg(&env)
        .arg("--flp")
        .arg(&flp)
        .output()
        .unwrap();
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    assert_product(&root.join("PekerisBroadband.h5"), &env, &flp, "kraken");
    fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn symlink_input_alias_cannot_be_overwritten() {
    let root = directory("symlink");
    let env = fixture("Pekeris.env");
    let output = root.join("alias.h5");
    std::os::unix::fs::symlink(&env, &output).unwrap();
    let old = fs::read(&env).unwrap();
    let process = run(&env, &output, "kraken", &["--overwrite"]);
    assert_eq!(process.status.code(), Some(4));
    assert_eq!(fs::read(&env).unwrap(), old);
    assert!(output.is_symlink());
    fs::remove_dir_all(root).unwrap();
}
