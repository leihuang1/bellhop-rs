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
        ("MunkLeakyPartialP", "krakenc"),
        ("MunkLeakyPartialS", "krakenc"),
        ("PekerisComplexSpline3", "krakenc"),
        ("MunkAnalyticComplex", "krakenc"),
        ("MunkLeakyPchipBroadband", "krakenc"),
        ("TabRefBrcN", "krakenc"),
        ("TabRefBrcC", "krakenc"),
        ("TabRefIrcN", "krakenc"),
        ("TabRefIrcC", "krakenc"),
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

#[test]
#[ignore = "requires unmodified upstream BroadBand/MunkK input snapshots"]
fn original_complex_spline_broadband_retains_the_work_limit_and_output() {
    let env = PathBuf::from(std::env::var_os("KRAKEN_DIFFERENTIAL_ENV").unwrap());
    let cases =
        kraken::legacy::load_frequency_cases(&env, env.with_extension("flp"), ModeSolver::Krakenc)
            .unwrap();
    assert_eq!(cases.len(), 2);
    assert_eq!(cases[1].frequency_hz, 500.0);
    assert_eq!(cases[1].interpolation, kraken::Interpolation::Spline);
    let root = directory("complex-spline-budget");
    let output = root.join("previous.h5");
    fs::write(&output, b"previous output").unwrap();
    // 50 Hz is written first; 500 Hz then exceeds the unchanged 300M root budget.
    let process = run(&env, &output, "krakenc", &["--overwrite"]);
    assert_failure(&process, 3, &output, b"previous output");
    let message = String::from_utf8_lossy(&process.stderr);
    assert!(
        message.contains("frequency[1] (500 Hz)")
            && message.contains("complex root work limit exceeded")
    );
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_round_trips_water_material_results() {
    let root = directory("water-materials");
    for name in [
        "WaterLossN",
        "WaterLossC",
        "WaterLossP",
        "WaterLossS",
        "WaterLossUnitN",
        "WaterLossUnitM",
        "WaterLossUnitF",
        "WaterLossUnitQ",
        "WaterLossUnitL",
        "WaterLossPower",
        "WaterLossThorp",
        "WaterLossFg",
        "WaterLossBio",
        "WaterLossLeaky",
    ] {
        for engine in ["kraken", "krakenc"] {
            if name == "WaterLossLeaky" && engine == "kraken" {
                continue;
            }
            let env = fixture(name).with_extension("env");
            let output = root.join(format!("{name}-{engine}.h5"));
            let process = run(&env, &output, engine, &[]);
            assert!(
                process.status.success(),
                "{}",
                String::from_utf8_lossy(&process.stderr)
            );
            assert_product(&output, &env, &env.with_extension("flp"), engine);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_round_trips_smooth_boundaries_and_top_tables() {
    let root = directory("fluid-boundaries");
    for name in [
        "FluidBoundaryVV",
        "FluidBoundaryVR",
        "FluidBoundaryVA",
        "FluidBoundaryRV",
        "FluidBoundaryRR",
        "FluidBoundaryRA",
        "FluidBoundaryAV",
        "FluidBoundaryAR",
        "FluidBoundaryAA",
        "FluidRigidPlaneLoss",
        "FluidBoundaryAir",
        "FluidTrcN",
        "FluidTrcC",
        "FluidTrcRigid",
    ] {
        for engine in ["kraken", "krakenc"] {
            if engine == "kraken" && (name == "FluidBoundaryAir" || name.starts_with("FluidTrc")) {
                continue;
            }
            let env = fixture(name).with_extension("env");
            let output = root.join(format!("{name}-{engine}.h5"));
            let process = run(&env, &output, engine, &[]);
            assert!(
                process.status.success(),
                "{}",
                String::from_utf8_lossy(&process.stderr)
            );
            assert_product(&output, &env, &env.with_extension("flp"), engine);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_round_trips_layered_fluids() {
    let root = directory("layered-fluids");
    for name in [
        "LayeredFluidN",
        "LayeredFluidC",
        "LayeredFluidP",
        "LayeredFluidS",
        "LayeredBoundaryVV",
        "LayeredBoundaryVR",
        "LayeredBoundaryRV",
        "LayeredBoundaryRR",
        "LayeredBoundaryRA",
        "LayeredBoundaryAV",
        "LayeredBoundaryAR",
        "LayeredBoundaryAA",
        "LayeredFluidThree",
        "LayeredFluidThreeWide",
        "LayeredFluidPlane",
        "LayeredFluidPower",
        "LayeredFluidBio",
        "LayeredFluidLeaky",
        "LayeredDoubleRefined",
        "LayeredNormalization",
        "LayeredFluidFractional",
    ] {
        for engine in ["kraken", "krakenc"] {
            if name == "LayeredFluidLeaky" && engine == "kraken" {
                continue;
            }
            let env = fixture(name).with_extension("env");
            let output = root.join(format!("{name}-{engine}.h5"));
            let process = run(&env, &output, engine, &[]);
            assert!(
                process.status.success(),
                "{}",
                String::from_utf8_lossy(&process.stderr)
            );
            assert_product(&output, &env, &env.with_extension("flp"), engine);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_round_trips_finite_elastic_layers() {
    let root = directory("finite-elastic-products");
    for name in [
        "FiniteElasticBothC",
        "FiniteElasticBothN",
        "FiniteElasticBothP",
        "FiniteElasticBothS",
        "FiniteElasticBottomC",
        "FiniteElasticBottomN",
        "FiniteElasticBottomP",
        "FiniteElasticBottomS",
        "FiniteElasticPower",
        "FiniteElasticRigid",
        "FiniteElasticShearOnly",
        "FiniteElasticStack",
        "FiniteElasticTopC",
        "FiniteElasticTopN",
        "FiniteElasticTopP",
        "FiniteElasticTopS",
        "FiniteElasticVacuum",
        "FiniteSingleIceC",
        "FiniteSingleIceP",
        "FiniteSingleIceS",
        "FiniteSingleSedimentC",
        "FiniteSingleSedimentP",
        "FiniteSingleSedimentS",
        "OriginalElasticIce",
        "OriginalElasticSediment",
    ] {
        for solver in ["kraken", "krakenc"] {
            if solver == "kraken" && name.starts_with("FiniteElastic") {
                continue;
            }
            let env = fixture(&format!("{name}.env"));
            let output = root.join(format!("{name}-{solver}.h5"));
            let process = run(&env, &output, solver, &[]);
            assert!(
                process.status.success(),
                "{}",
                String::from_utf8_lossy(&process.stderr)
            );
            assert_product(&output, &env, &env.with_extension("flp"), solver);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn finite_elastic_failures_preserve_output_and_remove_scratch() {
    let root = directory("finite-elastic-failures");
    let output = root.join("old.h5");
    let old = b"old finite elastic output";
    fs::write(&output, old).unwrap();
    let process = run(
        &fixture("FiniteElasticBothN.env"),
        &output,
        "kraken",
        &["--overwrite"],
    );
    assert_failure(&process, 2, &output, old);
    assert!(String::from_utf8_lossy(&process.stderr).contains("multi-fluid elastic secant parity"));
    let env = root.join("bad.env");
    fs::write(
        &env,
        fs::read_to_string(fixture("OriginalElasticSediment.env"))
            .unwrap()
            .replace("1300.0  2000.0", "2500.0  3000.0"),
    )
    .unwrap();
    fs::copy(
        fixture("OriginalElasticSediment.flp"),
        env.with_extension("flp"),
    )
    .unwrap();
    assert_failure(
        &run(&env, &output, "kraken", &["--overwrite"]),
        3,
        &output,
        old,
    );
    assert_failure(
        &run(
            &fixture("OriginalElasticIce.env"),
            &output,
            "krakenc",
            &["--overwrite", "--max-output-bytes", "1"],
        ),
        4,
        &output,
        old,
    );
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn cli_round_trips_elastic_half_spaces() {
    let root = directory("elastic-half-spaces");
    for name in [
        "ElasticHalfBottomN",
        "ElasticHalfBottomC",
        "ElasticHalfBottomP",
        "ElasticHalfBottomS",
        "ElasticHalfTopN",
        "ElasticHalfTopC",
        "ElasticHalfTopP",
        "ElasticHalfTopS",
        "ElasticHalfBothN",
        "ElasticHalfBothC",
        "ElasticHalfBothP",
        "ElasticHalfBothS",
        "ElasticHalfLeaky",
        "ElasticHalfPower",
        "ElasticHalfShearOnly",
        "OriginalElasticScholte",
        "OriginalElasticNormal",
        "OriginalElasticFlused",
    ] {
        for engine in ["kraken", "krakenc"] {
            if engine == "kraken"
                && (name.starts_with("ElasticHalfTop")
                    || name.starts_with("ElasticHalfBoth")
                    || name == "ElasticHalfLeaky")
            {
                continue;
            }
            let env = fixture(name).with_extension("env");
            let output = root.join(format!("{name}-{engine}.h5"));
            let process = run(&env, &output, engine, &[]);
            assert!(
                process.status.success(),
                "{}",
                String::from_utf8_lossy(&process.stderr)
            );
            assert_product(&output, &env, &env.with_extension("flp"), engine);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn finite_elastic_metadata_checks_the_last_frequency() {
    let root = directory("finite-elastic-metadata-corruption");
    let env = fixture("FiniteElasticPower.env");
    let output = root.join("power.h5");
    let process = run(&env, &output, "krakenc", &[]);
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    // Retain one RW handle; parallel children must not force a native lock upgrade.
    let file = File::open_rw(&output).unwrap();
    assert_product(&output, &env, &env.with_extension("flp"), "krakenc");
    file.group("frequencies/3/elastic_media/bottom/0")
        .unwrap()
        .attr("shear_attenuation_db_per_wavelength")
        .unwrap()
        .write_scalar(&0.25_f64)
        .unwrap();
    file.flush().unwrap();
    assert!(
        std::panic::catch_unwind(|| assert_product(
            &output,
            &env,
            &env.with_extension("flp"),
            "krakenc"
        ))
        .is_err()
    );
    drop(file);
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn elastic_metadata_checks_the_last_frequency() {
    let root = directory("elastic-metadata-corruption");
    let env = fixture("ElasticHalfPower.env");
    let output = root.join("power.h5");
    let process = run(&env, &output, "krakenc", &[]);
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    // Avoid a read-to-write reopen: parallel CLI children can inherit native read handles.
    let file = File::open_rw(&output).unwrap();
    assert_product(&output, &env, &env.with_extension("flp"), "krakenc");
    file.group("frequencies/3")
        .unwrap()
        .attr("bottom_shear_sound_speed_mps")
        .unwrap()
        .write_scalar(&1999.0_f64)
        .unwrap();
    file.flush().unwrap();
    assert!(
        std::panic::catch_unwind(|| assert_product(
            &output,
            &env,
            &env.with_extension("flp"),
            "krakenc"
        ))
        .is_err()
    );
    drop(file);
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn elastic_failures_preserve_output_and_remove_scratch() {
    let root = directory("elastic-failures");
    let output = root.join("previous.h5");
    fs::write(&output, b"old elastic output").unwrap();
    let process = run(
        &fixture("ElasticHalfTopN.env"),
        &output,
        "kraken",
        &["--overwrite"],
    );
    assert_failure(&process, 2, &output, b"old elastic output");
    assert!(String::from_utf8_lossy(&process.stderr).contains("elastic top requires KRAKENC"));

    let env = root.join("outside.env");
    fs::write(
        &env,
        fs::read_to_string(fixture("ElasticHalfBottomN.env"))
            .unwrap()
            .replace("1400.0 1800.0", "2500.0 3000.0"),
    )
    .unwrap();
    fs::copy(fixture("ElasticHalfBottomN.flp"), env.with_extension("flp")).unwrap();
    let process = run(&env, &output, "kraken", &["--overwrite"]);
    assert_failure(&process, 3, &output, b"old elastic output");
    assert!(String::from_utf8_lossy(&process.stderr).contains("phase-speed limits"));

    let process = run(
        &fixture("ElasticHalfBothN.env"),
        &output,
        "krakenc",
        &["--overwrite", "--max-output-bytes", "1"],
    );
    assert_failure(&process, 4, &output, b"old elastic output");
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn layered_failures_preserve_output_and_remove_scratch() {
    let root = directory("layered-failures");
    let output = root.join("previous.h5");
    let original = fixture("OriginalLayeredDouble.env");
    for engine in ["kraken", "krakenc"] {
        fs::write(&output, b"old layered output").unwrap();
        let process = run(&original, &output, engine, &["--overwrite"]);
        assert_failure(&process, 3, &output, b"old layered output");
        assert!(String::from_utf8_lossy(&process.stderr).contains("mode count changed"));
    }
    let env = root.join("later.env");
    fs::write(
        &env,
        fs::read_to_string(fixture("LayeredFluidPower.env"))
            .unwrap()
            .replace("75.0 50.0 62.5 50.0 /", "75.0 50.0 62.5 7500.0 /"),
    )
    .unwrap();
    fs::copy(fixture("LayeredFluidPower.flp"), env.with_extension("flp")).unwrap();
    let process = run(&env, &output, "krakenc", &["--overwrite"]);
    assert_failure(&process, 3, &output, b"old layered output");
    let message = String::from_utf8_lossy(&process.stderr);
    assert!(
        message.contains("frequency[3] (7500 Hz)")
            && message.contains("complex root work limit exceeded"),
        "{message}"
    );
    let env = fixture("LayeredFluidN.env");
    let process = run(
        &env,
        &output,
        "krakenc",
        &["--overwrite", "--max-output-bytes", "1"],
    );
    assert_failure(&process, 4, &output, b"old layered output");
    fs::remove_dir_all(root).unwrap();
}

fn assert_finite_elastic_metadata(group: &hdf5::Group, case: &kraken::Case) {
    let count = case.top_elastic_layers.len() + case.bottom_elastic_layers.len();
    assert_eq!(
        group
            .attr("finite_elastic_layer_count")
            .unwrap()
            .read_scalar::<u64>()
            .unwrap(),
        count as u64
    );
    if count == 0 {
        return;
    }
    for (side, layers, mut top) in [
        ("top", &case.top_elastic_layers, 0.0),
        (
            "bottom",
            &case.bottom_elastic_layers,
            case.fluid_bottom_depth_m(),
        ),
    ] {
        let media = group.group(&format!("elastic_media/{side}")).unwrap();
        assert_eq!(media.member_names().unwrap().len(), layers.len());
        for (index, layer) in layers.iter().enumerate() {
            let material = media.group(&index.to_string()).unwrap();
            assert_eq!(attribute(&material, "material"), "elastic");
            assert_eq!(
                attribute(&material, "attenuation_model"),
                if case.mode_solver == ModeSolver::Kraken {
                    "reference_real_stiffness"
                } else {
                    "complex"
                }
            );
            assert_eq!(
                material
                    .attr("requested_mesh_points")
                    .unwrap()
                    .read_scalar::<u64>()
                    .unwrap(),
                layer.mesh_points as u64
            );
            for (name, expected) in [
                ("top_depth_m", top),
                ("bottom_depth_m", layer.bottom_depth_m),
                (
                    "compressional_sound_speed_mps",
                    layer.compressional_sound_speed_mps,
                ),
                ("shear_sound_speed_mps", layer.shear_sound_speed_mps),
                ("density_g_cm3", layer.density_g_cm3),
                (
                    "compressional_attenuation_db_per_wavelength",
                    layer.compressional_attenuation_db_per_wavelength,
                ),
                (
                    "shear_attenuation_db_per_wavelength",
                    layer.shear_attenuation_db_per_wavelength,
                ),
            ] {
                let attr = material.attr(name).unwrap();
                assert!(attr.dtype().unwrap().is::<f64>());
                assert_eq!(
                    attr.read_scalar::<f64>().unwrap().to_bits(),
                    expected.to_bits()
                );
            }
            top = layer.bottom_depth_m;
        }
    }
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
    let mut inputs = vec![("env", env.to_path_buf()), ("flp", flp.to_path_buf())];
    if matches!(
        cases[0].surface_boundary,
        kraken::SurfaceBoundary::Reflection(_)
    ) {
        inputs.push(("trc", env.with_extension("trc")));
    }
    match &cases[0].bottom_boundary {
        kraken::BottomBoundary::Reflection(_) => inputs.push(("brc", env.with_extension("brc"))),
        kraken::BottomBoundary::Impedance { .. } => inputs.push(("irc", env.with_extension("irc"))),
        _ => {}
    }
    assert_eq!(
        file.group("inputs").unwrap().member_names().unwrap().len(),
        inputs.len()
    );
    for (name, path) in inputs {
        let bytes = fs::read(&path).unwrap();
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
        assert_eq!(
            attribute(&group, "surface_boundary"),
            match &case.surface_boundary {
                kraken::SurfaceBoundary::Vacuum => "V",
                kraken::SurfaceBoundary::Rigid => "R",
                kraken::SurfaceBoundary::FluidHalfSpace
                | kraken::SurfaceBoundary::ElasticHalfSpace { .. } => "A",
                kraken::SurfaceBoundary::Reflection(_) => "F",
                kraken::SurfaceBoundary::Impedance { .. } => "P",
            }
        );
        assert_eq!(
            attribute(&group, "bottom_boundary"),
            match &case.bottom_boundary {
                kraken::BottomBoundary::Vacuum => "V",
                kraken::BottomBoundary::FluidHalfSpace
                | kraken::BottomBoundary::ElasticHalfSpace { .. } => "A",
                kraken::BottomBoundary::Rigid => "R",
                kraken::BottomBoundary::Reflection(_) => "F",
                kraken::BottomBoundary::Impedance { .. } => "P",
            }
        );
        for (name, boundary, cp, density, loss) in [
            (
                "surface",
                &case.surface_boundary,
                case.surface_sound_speed_mps,
                case.surface_density_g_cm3,
                case.surface_attenuation_db_per_wavelength,
            ),
            (
                "bottom",
                &case.bottom_boundary,
                case.bottom_sound_speed_mps,
                case.bottom_density_g_cm3,
                case.bottom_attenuation_db_per_wavelength,
            ),
        ] {
            if let kraken::Boundary::ElasticHalfSpace {
                shear_sound_speed_mps,
                shear_attenuation_db_per_wavelength,
            } = boundary
            {
                assert_eq!(
                    attribute(&group, &format!("{name}_half_space_material")),
                    "elastic"
                );
                assert_eq!(
                    attribute(&group, &format!("{name}_elastic_attenuation_model")),
                    if engine == ModeSolver::Kraken {
                        "reference_real"
                    } else {
                        "complex"
                    }
                );
                for (name, value) in [
                    (format!("{name}_sound_speed_mps"), cp),
                    (format!("{name}_density_g_cm3"), density),
                    (format!("{name}_attenuation_db_per_wavelength"), loss),
                    (
                        format!("{name}_shear_sound_speed_mps"),
                        *shear_sound_speed_mps,
                    ),
                    (
                        format!("{name}_shear_attenuation_db_per_wavelength"),
                        *shear_attenuation_db_per_wavelength,
                    ),
                ] {
                    assert_eq!(
                        group.attr(&name).unwrap().read_scalar::<f64>().unwrap(),
                        value
                    );
                }
            }
        }
        assert_finite_elastic_metadata(&group, case);
        let layer_count = 1 + case.additional_fluid_layers.len();
        assert_eq!(
            group
                .attr("finite_fluid_layer_count")
                .unwrap()
                .read_scalar::<u64>()
                .unwrap(),
            layer_count as u64
        );
        let media = group.group("media").unwrap();
        assert_eq!(media.member_names().unwrap().len(), layer_count);
        let mut top_depth = case.fluid_top_depth_m();
        for index in 0..layer_count {
            let (bottom_depth, density, mesh) = if index == 0 {
                (
                    case.water_depth_m,
                    case.water_density_g_cm3,
                    case.mesh_points,
                )
            } else {
                let layer = &case.additional_fluid_layers[index - 1];
                (layer.bottom_depth_m, layer.density_g_cm3, layer.mesh_points)
            };
            let layer = media.group(&index.to_string()).unwrap();
            assert_eq!(
                layer
                    .attr("top_depth_m")
                    .unwrap()
                    .read_scalar::<f64>()
                    .unwrap(),
                top_depth
            );
            assert_eq!(
                layer
                    .attr("bottom_depth_m")
                    .unwrap()
                    .read_scalar::<f64>()
                    .unwrap(),
                bottom_depth
            );
            assert_eq!(
                layer
                    .attr("density_g_cm3")
                    .unwrap()
                    .read_scalar::<f64>()
                    .unwrap(),
                density
            );
            assert_eq!(
                layer
                    .attr("requested_mesh_points")
                    .unwrap()
                    .read_scalar::<u64>()
                    .unwrap(),
                mesh as u64
            );
            top_depth = bottom_depth;
        }
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

#[test]
fn cli_protects_consumed_tables_and_preserves_outputs_on_table_errors() {
    let root = directory("tables");
    for (name, extension) in [
        ("TabRefBrcC", "brc"),
        ("TabRefIrcC", "irc"),
        ("FluidTrcC", "trc"),
    ] {
        let env = root.join(name).with_extension("env");
        let table = env.with_extension(extension);
        for ext in ["env", "flp", extension] {
            fs::copy(fixture(name).with_extension(ext), env.with_extension(ext)).unwrap();
        }
        let bytes = fs::read(&table).unwrap();
        assert_failure(
            &run(&env, &table, "krakenc", &["--overwrite"]),
            4,
            &table,
            &bytes,
        );
        #[cfg(unix)]
        {
            let alias = root.join(format!("{name}-alias.h5"));
            std::os::unix::fs::symlink(&table, &alias).unwrap();
            assert_failure(
                &run(&env, &alias, "krakenc", &["--overwrite"]),
                4,
                &alias,
                &bytes,
            );
        }
        let output = root.join(name).with_extension("h5");
        let old = b"old complete output";
        fs::write(&output, old).unwrap();
        for bad in [b"0\n".to_vec(), vec![0xff], vec![b' '; 1_048_577]] {
            fs::write(&table, bad).unwrap();
            assert_failure(
                &run(&env, &output, "krakenc", &["--overwrite"]),
                2,
                &output,
                old,
            );
        }
        fs::remove_file(&table).unwrap();
        assert_failure(
            &run(&env, &output, "krakenc", &["--overwrite"]),
            2,
            &output,
            old,
        );
        if matches!(extension, "brc" | "trc") {
            // Valid finite rows, singular F impedance: failure after HDF5 creation.
            fs::write(&table, "2\n0.0 1.0 0.0\n90.0 1.0 0.0\n").unwrap();
            assert_failure(
                &run(&env, &output, "krakenc", &["--overwrite"]),
                3,
                &output,
                old,
            );
        }
        fs::write(&table, &bytes).unwrap();
        assert!(
            run(&env, &output, "krakenc", &["--overwrite"])
                .status
                .success()
        );
        assert_product(&output, &env, &env.with_extension("flp"), "krakenc");
    }
    fs::remove_dir_all(root).unwrap();
}

fn assert_failure(process: &Output, code: i32, output: &Path, old: &[u8]) {
    assert_eq!(
        process.status.code(),
        Some(code),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    assert_eq!(fs::read(output).unwrap(), old);
    let mut scratch = output.as_os_str().to_os_string();
    scratch.push(".tmp");
    assert!(!PathBuf::from(scratch).exists());
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
