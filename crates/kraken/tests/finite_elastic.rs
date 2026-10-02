use kraken::{Boundary, Case, ModeSolver, legacy, solve};
use std::path::{Path, PathBuf};

fn original(name: &str, solver: ModeSolver) -> Case {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    legacy::load_frequency_cases(
        root.with_extension("env"),
        root.with_extension("flp"),
        solver,
    )
    .unwrap()
    .remove(0)
}

#[test]
fn finite_elastic_layers_keep_the_absolute_fluid_interval() {
    let ice = original("OriginalElasticIce", ModeSolver::Krakenc);
    assert_eq!(
        (
            ice.fluid_top_depth_m(),
            ice.fluid_bottom_depth_m(),
            ice.total_depth_m()
        ),
        (30.0, 5000.0, 5000.0)
    );
    assert_eq!(
        ice.sound_speed_profile[0].depth_m.to_bits(),
        30.0_f64.to_bits()
    );
    assert_eq!(ice.top_elastic_layers.len(), 1);
    assert_eq!(ice.mesh_points, 500);
    let mut input = ice.into_definition();
    input.source_depths_m[0] = 10.0;
    assert!(
        Case::from_definition(input)
            .unwrap_err()
            .diagnostics()
            .iter()
            .any(|d| d.field == "source_depths_m")
    );
    let sediment = original("OriginalElasticSediment", ModeSolver::Kraken);
    assert_eq!(
        (
            sediment.fluid_top_depth_m(),
            sediment.fluid_bottom_depth_m(),
            sediment.total_depth_m()
        ),
        (0.0, 5000.0, 5100.0)
    );
    let mut input = sediment.into_definition();
    input.mode_sample_depths_m[0] = 5100.0;
    assert!(
        Case::from_definition(input)
            .unwrap_err()
            .diagnostics()
            .iter()
            .any(|d| d.field == "mode_sample_depths_m")
    );
}

#[test]
fn finite_elastic_materials_and_unvalidated_combinations_are_rejected() {
    let base = original("OriginalElasticIce", ModeSolver::Krakenc).into_definition();
    for change in ["cp", "cs", "density", "p_loss", "s_loss", "bulk", "depth"] {
        let mut input = base.clone();
        let layer = &mut input.top_elastic_layers[0];
        match change {
            "cp" => layer.compressional_sound_speed_mps = f64::NAN,
            "cs" => layer.shear_sound_speed_mps = f64::INFINITY,
            "density" => layer.density_g_cm3 = 0.0,
            "p_loss" => layer.compressional_attenuation_db_per_wavelength = -1.0,
            "s_loss" => layer.shear_attenuation_db_per_wavelength = f64::NAN,
            "bulk" => layer.compressional_sound_speed_mps = 1500.0,
            "depth" => layer.bottom_depth_m = f64::INFINITY,
            _ => unreachable!(),
        }
        assert!(
            Case::from_definition(input)
                .unwrap_err()
                .diagnostics()
                .iter()
                .any(|d| d.field.starts_with("top_elastic_layers[0]")),
            "{change}"
        );
    }
    let mut input = original("OriginalElasticSediment", ModeSolver::Krakenc).into_definition();
    input.bottom_boundary = Boundary::FluidHalfSpace;
    assert!(
        Case::from_definition(input)
            .unwrap_err()
            .to_string()
            .contains("compound state is undefined")
    );
    let env = include_str!("fixtures/OriginalElasticIce.env")
        .replace("30.0 3000.0 1400.0 1.0", "30.0 3100.0 1400.0 1.0");
    assert!(
        legacy::load_frequency_cases_from_sources(
            &env,
            include_str!("fixtures/OriginalElasticIce.flp"),
            Path::new("ice.env"),
            Path::new("ice.flp"),
            ModeSolver::Krakenc
        )
        .unwrap_err()
        .to_string()
        .contains("homogeneous elastic")
    );
    let mut input = base;
    input.top_elastic_layers[0].mesh_points = 1_000_000;
    assert!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap_err()
            .to_string()
            .contains("mesh limit")
    );
}

#[test]
fn finite_loss_and_mesh_scaling_preserve_frequency_order() {
    let cases = legacy::load_frequency_cases_from_sources(
        include_str!("fixtures/FiniteElasticPower.env"),
        include_str!("fixtures/FiniteElasticPower.flp"),
        Path::new("power.env"),
        Path::new("power.flp"),
        ModeSolver::Krakenc,
    )
    .unwrap();
    assert_eq!(
        cases.iter().map(|c| c.frequency_hz).collect::<Vec<_>>(),
        [75.0, 50.0, 62.5, 50.0]
    );
    for (case, factor) in cases
        .iter()
        .zip([1.1_f64.powf(0.25), 1.0, 1.1_f64.powf(0.25), 1.0])
    {
        assert!(
            (case.bottom_elastic_layers[0].shear_attenuation_db_per_wavelength - 0.3 * factor)
                .abs()
                < 1e-14
        );
    }
    assert_eq!(cases[1], cases[3]);
    let mut input = cases[0].clone().into_definition();
    input.mode_solver = ModeSolver::Kraken;
    assert!(
        Case::from_definition(input)
            .unwrap_err()
            .to_string()
            .contains("multi-fluid elastic secant parity")
    );
}

#[test]
fn real_finite_stiffness_is_not_the_lossless_half_space_approximation() {
    let case = original("OriginalElasticIce", ModeSolver::Kraken);
    let with_loss = solve(&case).unwrap();
    assert!(
        with_loss
            .modes
            .modes
            .iter()
            .all(|m| m.attenuation_nepers_per_m == 0.0)
    );
    let mut input = case.into_definition();
    input.top_elastic_layers[0].compressional_attenuation_db_per_wavelength = 0.0;
    input.top_elastic_layers[0].shear_attenuation_db_per_wavelength = 0.0;
    let lossless = solve(&Case::from_definition(input).unwrap()).unwrap();
    assert!(
        with_loss
            .modes
            .modes
            .iter()
            .zip(lossless.modes.modes)
            .any(
                |(a, b)| (a.horizontal_wavenumber_rad_per_m - b.horizontal_wavenumber_rad_per_m)
                    .norm()
                    > 1e-12
            )
    );
    let complex = solve(&original("OriginalElasticIce", ModeSolver::Krakenc)).unwrap();
    assert!(
        complex
            .modes
            .modes
            .iter()
            .all(|m| m.attenuation_nepers_per_m > 0.0)
    );
}
