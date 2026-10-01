use std::fs;
use std::path::{Path, PathBuf};

use kraken::{BottomBoundary, Case, Interpolation, ModeSolver, SurfaceBoundary, legacy};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn table_snapshots_do_not_reread_files_and_old_snapshot_api_requires_the_resource() {
    for (name, extension) in [
        ("TabRefBrcN", "brc"),
        ("TabRefIrcC", "irc"),
        ("FluidTrcC", "trc"),
    ] {
        let root = fixture(name);
        let env = fs::read_to_string(root.with_extension("env")).unwrap();
        let flp = fs::read_to_string(root.with_extension("flp")).unwrap();
        let table = fs::read_to_string(root.with_extension(extension)).unwrap();
        let env_path = Path::new("not-on-disk.env");
        let flp_path = Path::new("not-on-disk.flp");
        assert_eq!(
            if extension == "trc" {
                legacy::surface_table_extension(&env, env_path, ModeSolver::Krakenc)
            } else {
                legacy::bottom_table_extension(&env, env_path, ModeSolver::Krakenc)
            }
            .unwrap(),
            Some(extension)
        );
        let from_sources = legacy::load_frequency_cases_with_boundary_tables(
            &env,
            &flp,
            env_path,
            flp_path,
            ModeSolver::Krakenc,
            (extension == "trc").then_some(table.as_str()),
            (extension != "trc").then_some(table.as_str()),
        )
        .unwrap();
        let from_files =
            legacy::load_complex_case(root.with_extension("env"), root.with_extension("flp"))
                .unwrap();
        assert_eq!(from_sources, vec![from_files]);
        let report = legacy::load_frequency_cases_from_sources(
            &env,
            &flp,
            env_path,
            flp_path,
            ModeSolver::Krakenc,
        )
        .unwrap_err();
        assert_eq!(
            report.diagnostics()[0].path,
            env_path.with_extension(extension)
        );
        assert!(report.to_string().contains("snapshot is missing"));
        assert!(legacy::load_case(root.with_extension("env"), root.with_extension("flp")).is_err());
    }
}

#[test]
fn surface_table_limits_are_validated_at_the_public_case_boundary() {
    let root = fixture("FluidTrcN");
    let input = legacy::load_complex_case(root.with_extension("env"), root.with_extension("flp"))
        .unwrap()
        .into_definition();
    for change in [
        "solver",
        "interpolation",
        "loss",
        "refinement",
        "broadband",
        "evanescent",
        "material",
        "topP",
        "bottomTable",
        "count",
    ] {
        let mut invalid = input.clone();
        match change {
            "solver" => invalid.mode_solver = ModeSolver::Kraken,
            "interpolation" => invalid.interpolation = Interpolation::Spline,
            "loss" => {
                invalid.water_attenuation_db_per_wavelength =
                    vec![0.1; invalid.sound_speed_profile.len()];
            }
            "refinement" => invalid.max_range_m = 1.0,
            "broadband" => invalid.mesh_reference_frequency_hz = Some(50.0),
            "evanescent" => invalid.c_low_mps = 1400.0,
            "material" => invalid.surface_density_g_cm3 = 1.0,
            "topP" => {
                invalid.surface_boundary = SurfaceBoundary::Impedance {
                    frequency_hz: 50.0,
                    points: Vec::new(),
                }
            }
            "bottomTable" => {
                invalid.bottom_boundary = invalid.surface_boundary.clone();
                invalid.bottom_sound_speed_mps = 0.0;
                invalid.bottom_density_g_cm3 = 0.0;
                invalid.bottom_attenuation_db_per_wavelength = 0.0;
            }
            _ => {
                if let SurfaceBoundary::Reflection(p) = &mut invalid.surface_boundary {
                    p.truncate(1);
                }
            }
        }
        assert!(Case::from_definition(invalid).is_err(), "{change}");
    }
    let env = fs::read_to_string(root.with_extension("env")).unwrap();
    let flp = fs::read_to_string(root.with_extension("flp")).unwrap();
    let old = legacy::load_frequency_cases_with_bottom_table(
        &env,
        &flp,
        &root.with_extension("env"),
        &root.with_extension("flp"),
        ModeSolver::Krakenc,
        None,
    )
    .unwrap_err();
    assert_eq!(old.diagnostics()[0].path, root.with_extension("trc"));
}

#[test]
fn table_combinations_are_validated_at_the_public_case_boundary() {
    for name in ["TabRefBrcC", "TabRefIrcN"] {
        let root = fixture(name);
        let input =
            legacy::load_complex_case(root.with_extension("env"), root.with_extension("flp"))
                .unwrap()
                .into_definition();
        for change in [
            "solver",
            "interpolation",
            "surface",
            "refinement",
            "broadband",
            "material",
        ] {
            let mut invalid = input.clone();
            match change {
                "solver" => invalid.mode_solver = ModeSolver::Kraken,
                "interpolation" => invalid.interpolation = Interpolation::Spline,
                "surface" => invalid.surface_boundary = SurfaceBoundary::Rigid,
                "refinement" => invalid.max_range_m = 1.0,
                "broadband" => invalid.mesh_reference_frequency_hz = Some(50.0),
                _ => invalid.bottom_density_g_cm3 = 1.0,
            }
            assert!(Case::from_definition(invalid).is_err(), "{name} {change}");
        }
        for change in ["count", "order", "finite", "scale", "frequency"] {
            let mut invalid = input.clone();
            match &mut invalid.bottom_boundary {
                BottomBoundary::Reflection(points) => match change {
                    "count" => points.truncate(1),
                    "order" => points[1].angle_degrees = points[0].angle_degrees,
                    "finite" => points[1].phase_radians = f64::NAN,
                    "scale" => points[0].magnitude = -1.0,
                    _ => continue,
                },
                BottomBoundary::Impedance {
                    frequency_hz,
                    points,
                } => match change {
                    "count" => points.truncate(1),
                    "order" => points[1].wavenumber_squared = points[0].wavenumber_squared,
                    "finite" => points[1].f.im = f64::INFINITY,
                    "scale" => {
                        points.truncate(2);
                        points[1].power = 101;
                    }
                    _ => *frequency_hz = 51.0,
                },
                _ => unreachable!(),
            }
            assert!(Case::from_definition(invalid).is_err(), "{name} {change}");
        }
    }
}

#[test]
fn parsers_bound_and_locate_bad_tables() {
    for (name, extension) in [("TabRefBrcC", "brc"), ("TabRefIrcC", "irc")] {
        let root = fixture(name);
        let env = fs::read_to_string(root.with_extension("env")).unwrap();
        let flp = fs::read_to_string(root.with_extension("flp")).unwrap();
        let table = fs::read_to_string(root.with_extension(extension)).unwrap();
        let env_path = root.with_extension("env");
        let flp_path = root.with_extension("flp");
        for bad in [
            String::new(),
            "100001\n".to_owned(),
            table.clone() + "extra\n",
            " ".repeat(1_048_577),
        ] {
            let report = legacy::load_frequency_cases_with_bottom_table(
                &env,
                &flp,
                &env_path,
                &flp_path,
                ModeSolver::Krakenc,
                Some(&bad),
            )
            .unwrap_err();
            assert_eq!(report.diagnostics()[0].path, root.with_extension(extension));
        }
    }
    let root = fixture("TabRefIrcC");
    let env = fs::read_to_string(root.with_extension("env")).unwrap();
    let flp = fs::read_to_string(root.with_extension("flp")).unwrap();
    let table = fs::read_to_string(root.with_extension("irc")).unwrap();
    let bad = table.replace(" 50.0\n3", " 51.0\n3");
    let report = legacy::load_frequency_cases_with_bottom_table(
        &env,
        &flp,
        &root.with_extension("env"),
        &root.with_extension("flp"),
        ModeSolver::Krakenc,
        Some(&bad),
    )
    .unwrap_err();
    assert_eq!(report.diagnostics()[0].path, root.with_extension("irc"));
    assert!(report.to_string().contains("frequency"));
    // Fixed-width rows cannot be rewritten as whitespace-delimited records.
    let bad = table
        .lines()
        .enumerate()
        .map(|(i, line)| {
            if i < 2 {
                line.to_owned()
            } else {
                line.split_whitespace().collect::<Vec<_>>().join(" ")
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        legacy::load_frequency_cases_with_bottom_table(
            &env,
            &flp,
            &root.with_extension("env"),
            &root.with_extension("flp"),
            ModeSolver::Krakenc,
            Some(&bad)
        )
        .is_err()
    );
}

#[test]
fn field_endpoint_extension_is_bounded_and_stays_in_water() {
    let root = fixture("TabRefBrcC");
    let mut input =
        legacy::load_complex_case(root.with_extension("env"), root.with_extension("flp"))
            .unwrap()
            .into_definition();
    input.mode_sample_depths_m = vec![1.0, 2.0, 98.0, 99.0];
    input.source_depths_m = vec![0.0, 100.0];
    assert!(Case::from_definition(input.clone()).is_ok());
    input.frequency_hz = 1501.0; // one-metre extension exceeds 1500/f
    assert!(Case::from_definition(input.clone()).is_err());
    input.frequency_hz = 50.0;
    input.source_depths_m[0] = -0.01;
    assert!(Case::from_definition(input.clone()).is_err());
    input.source_depths_m = vec![50.0];
    input.receiver_depths_m = vec![50.0];
    input.receiver_offsets_m = vec![0.0];
    input.mode_sample_depths_m = vec![50.0];
    assert!(Case::from_definition(input.clone()).is_ok());
    input.receiver_depths_m = vec![51.0]; // a singleton cannot extrapolate
    assert!(Case::from_definition(input).is_err());
}
