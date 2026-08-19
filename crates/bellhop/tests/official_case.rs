use std::path::{Path, PathBuf};

use bellhop::legacy::load_case;
use bellhop::model::BoundaryInterpolation;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn loads_official_range_dependent_case() {
    let case = load_case(&fixture("Gulf_ray_rd.env")).unwrap().value;
    let sound_speed = case.range_dependent_sound_speed.unwrap();
    assert_eq!(sound_speed.ranges_m.len(), 9);
    assert_eq!(sound_speed.depths_m.len(), 10);
    assert_eq!(sound_speed.speeds_mps.len(), 10);
    assert_eq!(sound_speed.speeds_mps[0].len(), 9);

    let bathymetry = case.bathymetry.unwrap();
    assert_eq!(bathymetry.points.len(), 8);
    assert_eq!(
        bathymetry.interpolation,
        BoundaryInterpolation::PiecewiseLinear
    );
}

#[test]
#[allow(clippy::float_cmp)]
fn loads_official_long_boundary_format() {
    let case = load_case(&fixture("PekerisRDB.env")).unwrap().value;
    let bathymetry = case.bathymetry.unwrap();
    assert_eq!(bathymetry.points.len(), 3);
    let material = bathymetry.points[0].material.as_ref().unwrap();
    assert_eq!(material.compressional_speed_mps, 1700.0);
    assert_eq!(material.density_g_cm3, 1.2);
}

#[test]
fn loads_official_source_beam_pattern() {
    let case = load_case(&fixture("shaded.env")).unwrap().value;
    assert!(case.environment.run.has_source_beam_pattern);
    let pattern = case.source_beam_pattern.unwrap();
    assert_eq!(pattern.points.len(), 37);
    assert!((pattern.points[0].amplitude - 10.0_f64.sqrt()).abs() < 1.0e-14);
}
