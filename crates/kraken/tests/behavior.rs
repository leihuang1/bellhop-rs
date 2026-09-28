use std::f64::consts::PI;
use std::path::Path;

use kraken::{Case, CaseDefinition, legacy::load_case, solve};
use num_complex::Complex64;

fn definition() -> CaseDefinition {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/Pekeris");
    load_case(root.with_extension("env"), root.with_extension("flp"))
        .unwrap()
        .into_definition()
}

#[test]
fn validation_collects_errors_and_bounds_allocations() {
    let mut input = definition();
    input.frequency_hz = f64::NAN;
    input.bottom_density_g_cm3 = -1.0;
    input.source_depths_m = vec![101.0];
    input.receiver_offsets_m.clear();
    let report = Case::from_definition(input).unwrap_err();
    for field in [
        "frequency_hz",
        "bottom_density_g_cm3",
        "source_depths_m",
        "receiver_offsets_m",
    ] {
        assert!(
            report.diagnostics().iter().any(|d| d.field == field),
            "{report}"
        );
    }

    let mut input = definition();
    input.source_depths_m = vec![75.0; 100_000];
    input.receiver_depths_m = vec![75.0; 100_000];
    input.receiver_offsets_m = vec![0.0; 100_000];
    assert!(
        Case::from_definition(input)
            .unwrap_err()
            .diagnostics()
            .iter()
            .any(|d| d.field == "field_grid")
    );

    for bad in [f64::NAN, f64::INFINITY, -1.0, 101.0] {
        let mut input = definition();
        input.mode_sample_depths_m[1] = bad;
        assert!(Case::from_definition(input).is_err());
    }
    let mut input = definition();
    input.receiver_ranges_m = vec![500.0, 500.0];
    assert!(Case::from_definition(input).is_err());
}

#[test]
fn phase_speed_selection_and_field_mode_limit_are_independent() {
    let input = definition();
    let baseline = solve(&Case::from_definition(input.clone()).unwrap()).unwrap();
    let mut selected = input.clone();
    selected.c_low_mps = baseline.modes.modes[0]
        .phase_speed_mps
        .midpoint(baseline.modes.modes[1].phase_speed_mps);
    selected.c_high_mps = baseline.modes.modes[1]
        .phase_speed_mps
        .midpoint(baseline.modes.modes[2].phase_speed_mps);
    let result = solve(&Case::from_definition(selected).unwrap()).unwrap();
    assert_eq!(result.modes.modes.len(), 1);
    assert!(
        (result.modes.modes[0].horizontal_wavenumber_rad_per_m
            - baseline.modes.modes[1].horizontal_wavenumber_rad_per_m)
            .norm()
            < 1e-13
    );

    let mut limited = input;
    limited.mode_limit = 1;
    let result = solve(&Case::from_definition(limited).unwrap()).unwrap();
    assert_eq!(result.modes, baseline.modes);
    assert!((result.field.pressure[0] - baseline.field.pressure[0]).norm() > 1e-3);
    let mode = &result.modes.modes[0];
    let k = mode.horizontal_wavenumber_rad_per_m;
    let expected = Complex64::new(0.0, 1.0)
        * (2.0 * PI).sqrt()
        * Complex64::from_polar(1.0, PI / 4.0)
        * mode.eigenfunction[1]
        * mode.eigenfunction[0]
        / k
        * (Complex64::new(0.0, -500.0) * k).exp();
    // FIELD rounds intermediate modal products and the sum to single precision.
    assert!((result.field.pressure[0] - expected).norm() < 2e-7);
}

#[test]
fn offsets_interpolation_and_multiple_sources_preserve_field_layout() {
    let mut input = definition();
    input.source_depths_m = vec![25.0, 50.0, 75.0];
    input.receiver_depths_m = vec![25.0, 50.0, 75.0];
    input.receiver_offsets_m = vec![0.0; 3];
    let result = solve(&Case::from_definition(input.clone()).unwrap()).unwrap();
    // 50 m is halfway between the stored 25 and 75 m samples, for either end.
    for range in 0..3 {
        let p = &result.field.pressure;
        assert!((p[3 + range] - (p[range] + p[6 + range]) * 0.5).norm() < 2e-7);
        assert!((p[9 + range] - (p[range] + p[18 + range]) * 0.5).norm() < 2e-7);
    }
    input.receiver_offsets_m = vec![10.0; 3];
    let shifted = solve(&Case::from_definition(input.clone()).unwrap()).unwrap();
    input.receiver_offsets_m = vec![0.0; 3];
    for range in &mut input.receiver_ranges_m {
        *range += 10.0;
    }
    let moved = solve(&Case::from_definition(input).unwrap()).unwrap();
    // Separate offset/range exponentials have different single-precision rounding.
    for (shifted, moved) in shifted.field.pressure.iter().zip(moved.field.pressure) {
        assert!((*shifted - moved).norm() < 2e-7);
    }
}

#[test]
fn trapped_mode_cutoff_and_group_speed() {
    let mut input = definition();
    let cutoff_frequency = 1.0
        / (4.0
            * input.water_depth_m
            * (input.sound_speed_profile[0].sound_speed_mps.recip().powi(2)
                - input.bottom_sound_speed_mps.recip().powi(2))
            .sqrt());
    input.frequency_hz = cutoff_frequency * (1.0 - 1e-6);
    assert_eq!(
        solve(&Case::from_definition(input.clone()).unwrap())
            .unwrap_err()
            .diagnostics()[0]
            .code,
        "KR0301"
    );
    // KRAKEN's xMin=1.00001*omega²/cHigh² excludes modes arbitrarily close to cutoff.
    input.frequency_hz = cutoff_frequency * 1.05;
    assert_eq!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap()
            .modes
            .modes
            .len(),
        1
    );

    let mut input = definition();
    input.max_range_m = 0.0; // group velocity is computed on the first mesh, not extrapolated.
    let baseline = solve(&Case::from_definition(input.clone()).unwrap()).unwrap();
    let mut lower = input.clone();
    let mut upper = input;
    lower.frequency_hz -= 0.001;
    upper.frequency_hz += 0.001;
    let lower = solve(&Case::from_definition(lower).unwrap()).unwrap();
    let upper = solve(&Case::from_definition(upper).unwrap()).unwrap();
    for (index, mode) in baseline.modes.modes.iter().enumerate() {
        let dk = upper.modes.modes[index].horizontal_wavenumber_rad_per_m.re
            - lower.modes.modes[index].horizontal_wavenumber_rad_per_m.re;
        let numerical_group_speed = 2.0 * PI * 0.002 / dk;
        assert!((mode.group_speed_mps - numerical_group_speed).abs() < 1e-5);
    }
}

#[test]
fn modal_roots_do_not_depend_on_an_absolute_q_tolerance() {
    let mut input = definition();
    let baseline = solve(&Case::from_definition(input.clone()).unwrap()).unwrap();
    let scale = 1e12;
    input.frequency_hz /= scale;
    input.water_depth_m *= scale;
    for point in &mut input.sound_speed_profile {
        point.depth_m *= scale;
    }
    input.max_range_m *= scale;
    for values in [
        &mut input.mode_sample_depths_m,
        &mut input.source_depths_m,
        &mut input.receiver_depths_m,
        &mut input.receiver_ranges_m,
        &mut input.receiver_offsets_m,
    ] {
        for value in values {
            *value *= scale;
        }
    }
    let scaled = solve(&Case::from_definition(input).unwrap()).unwrap();
    assert_eq!(baseline.modes.modes.len(), scaled.modes.modes.len());
    for (baseline, scaled) in baseline.modes.modes.iter().zip(scaled.modes.modes) {
        assert!(
            (baseline.horizontal_wavenumber_rad_per_m
                - scaled.horizontal_wavenumber_rad_per_m * scale)
                .norm()
                < 1e-12
        );
    }
}

#[test]
fn solver_limits_and_numeric_overflow_are_errors() {
    let mut input = definition();
    input.frequency_hz = 1_000_000.0;
    input.mesh_points = 1_000_000;
    assert_eq!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap_err()
            .diagnostics()[0]
            .field,
        "mode_count"
    );

    let mut input = definition();
    input.frequency_hz = f64::MAX;
    assert_eq!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap_err()
            .diagnostics()[0]
            .code,
        "KR0302"
    );

    let mut input = definition();
    input.frequency_hz = 1000.0; // Many modes, but a small mode-sample grid.
    input.source_depths_m = vec![75.0; 100_000]; // 900,000 pressure samples: allowed by Case.
    assert_eq!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap_err()
            .diagnostics()[0]
            .field,
        "field_grid"
    );

    let mut input = definition();
    input.frequency_hz = 1000.0;
    input.mode_sample_depths_m = (0..100_000).map(|i| f64::from(i) * 0.001).collect();
    assert_eq!(
        solve(&Case::from_definition(input).unwrap())
            .unwrap_err()
            .diagnostics()[0]
            .field,
        "mode_sample_depths_m"
    );
}
