//! Characterization through the unchanged Case/load/solve interfaces.
use kraken::{Case, Interpolation, ModeSolver, legacy, solve};
use std::path::Path;

#[test]
fn profile_diagnostics_keep_first_and_later_layer_rules_and_order() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let base = legacy::load_complex_case(
        root.join("LayeredFluidN.env"),
        root.join("LayeredFluidN.flp"),
    )
    .unwrap()
    .into_definition();
    for first in [true, false] {
        for kind in [
            "empty",
            "zero",
            "nan_speed",
            "depth",
            "short_loss",
            "gain",
            "nan_loss",
            "excess_loss",
            "spline",
        ] {
            let mut input = base.clone();
            let (points, losses) = if first {
                (
                    &mut input.sound_speed_profile,
                    &mut input.water_attenuation_db_per_wavelength,
                )
            } else {
                let layer = &mut input.additional_fluid_layers[0];
                (
                    &mut layer.sound_speed_profile,
                    &mut layer.attenuation_db_per_wavelength,
                )
            };
            match kind {
                "empty" => points.clear(),
                "zero" => points[1].sound_speed_mps = 0.0,
                "nan_speed" => points[1].sound_speed_mps = f64::NAN,
                "depth" => points[1].depth_m = points[0].depth_m,
                "short_loss" => {
                    losses.pop();
                }
                "gain" => losses[1] = -0.1,
                "nan_loss" => losses[1] = f64::NAN,
                "excess_loss" => losses[1] = 55.0,
                "spline" => {
                    *losses = vec![0.02, 0.0, 0.0];
                    input.interpolation = Interpolation::Spline;
                }
                _ => unreachable!(),
            }
            let profile_field = if first {
                "sound_speed_profile"
            } else {
                "additional_fluid_layers[0].sound_speed_profile"
            };
            let loss_field = if first {
                "water_attenuation_db_per_wavelength"
            } else {
                "additional_fluid_layers[0].attenuation_db_per_wavelength"
            };
            let expected_fields = match kind {
                "empty" => vec![profile_field, loss_field],
                "nan_speed" if !first => vec![profile_field, profile_field],
                "nan_loss" if !first => vec![loss_field, profile_field],
                "spline" if !first => vec![loss_field, loss_field],
                "zero" | "nan_speed" | "depth" => vec![profile_field],
                _ => vec![loss_field],
            };
            let report = Case::from_definition(input).unwrap_err();
            assert_eq!(
                report
                    .diagnostics()
                    .iter()
                    .map(|d| d.field.as_str())
                    .collect::<Vec<_>>(),
                expected_fields,
                "{first} {kind}: {report}"
            );
            for d in report.diagnostics() {
                assert_eq!(
                    (d.code, d.path.as_path(), d.line, d.column),
                    ("KR0201", Path::new("<case>"), 1, 1)
                );
            }
            if kind == "spline" {
                assert!(report.diagnostics().iter().all(
                    |d| d.message == "interpolated complex speed requires 0 <= Im(c) <= Re(c)"
                ));
            }
        }
    }
}

#[test]
fn each_frequency_converts_its_own_raw_materials_not_the_previous_case() {
    for (env, flp) in [
        (
            include_str!("fixtures/WaterLossPower.env"),
            include_str!("fixtures/WaterLossPower.flp"),
        ),
        (
            include_str!("fixtures/LayeredFluidPower.env"),
            include_str!("fixtures/LayeredFluidPower.flp"),
        ),
        (
            include_str!("fixtures/ElasticHalfPower.env"),
            include_str!("fixtures/ElasticHalfPower.flp"),
        ),
        (
            include_str!("fixtures/FiniteElasticPower.env"),
            include_str!("fixtures/FiniteElasticPower.flp"),
        ),
    ] {
        let load = |source: &str| {
            legacy::load_frequency_cases_from_sources(
                source,
                flp,
                Path::new("material.env"),
                Path::new("material.flp"),
                ModeSolver::Krakenc,
            )
            .unwrap()
        };
        let cases = load(env);
        for case in cases {
            let single = env.replace(
                "4\n75.0 50.0 62.5 50.0 /",
                &format!("1\n{} /", case.frequency_hz),
            );
            assert_ne!(single, env);
            assert_eq!(load(&single), [case]);
        }
    }
}

#[test]
fn refinement_keeps_shapes_and_group_speeds_from_the_first_mesh() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    for solver in [ModeSolver::Kraken, ModeSolver::Krakenc] {
        let refined = legacy::load_frequency_cases(
            root.join("PekerisRefined.env"),
            root.join("PekerisRefined.flp"),
            solver,
        )
        .unwrap()
        .remove(0);
        let mut first = refined.clone().into_definition();
        first.max_range_m = 0.0;
        let first = solve(&Case::from_definition(first).unwrap()).unwrap();
        let refined = solve(&refined).unwrap();
        assert_eq!(refined.modes.modes.len(), first.modes.modes.len());
        for (a, b) in refined.modes.modes.iter().zip(&first.modes.modes) {
            assert_eq!(a.eigenfunction, b.eigenfunction);
            assert_eq!(a.group_speed_mps.to_bits(), b.group_speed_mps.to_bits());
            assert_ne!(
                a.horizontal_wavenumber_rad_per_m,
                b.horizontal_wavenumber_rad_per_m
            );
        }
    }
}
