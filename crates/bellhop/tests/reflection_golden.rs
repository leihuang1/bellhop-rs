use std::f64::consts::PI;
use std::path::Path;

use bellhop::legacy::load_case;
use bellhop::solver::{SimulationLimits, run};

fn first_reflected_point(case_name: &str) -> bellhop::solver::RayPoint {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/golden")
        .join(case_name);
    let case = load_case(&path).unwrap().value;
    let result = run(&case, SimulationLimits::default()).unwrap();
    result.sources[0].rays[0]
        .points
        .iter()
        .copied()
        .find(|point| point.amplitude < 0.99)
        .expect("ray must have a lossy bottom reflection")
}

#[test]
fn acousto_elastic_reflection_matches_v2023_5_amplitude_and_phase() {
    let reflected = first_reflected_point("ElasticReflection.env");
    assert!((reflected.amplitude - 0.932_224_853_630_714_3).abs() < 2.0e-14);
    assert!((reflected.phase_radians - 0.768_934_880_079_570_1).abs() < 2.0e-14);
}

#[test]
fn grain_size_reflection_matches_v2023_5_amplitude_and_phase() {
    let reflected = first_reflected_point("GrainReflection.env");
    assert!((reflected.amplitude - 0.878_176_295_538_360_7).abs() < 2.0e-14);
    assert!((reflected.phase_radians - 0.396_522_146_363_542_27).abs() < 2.0e-14);
}

#[test]
fn precalculated_internal_reflection_matches_bounce_v2023_5() {
    let reflected = first_reflected_point("InternalReflection.env");

    // BOUNCE v2023.5's .brc output at 27.858558377213921 degrees.
    let expected_magnitude = 0.251_283_475_629_951_57;
    let expected_phase = 2.015_890_066_844_079_f64 * PI / 180.0;
    assert!((reflected.amplitude - expected_magnitude).abs() < 2.0e-6);
    assert!((reflected.phase_radians - expected_phase).abs() < 2.0e-6);
}
