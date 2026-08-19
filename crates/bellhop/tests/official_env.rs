use std::path::PathBuf;

use bellhop::legacy::load_env;
use bellhop::model::{RunKind, SspInterpolation, VolumeAttenuation};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn parses_representative_official_2d_environments() {
    for name in [
        "calibB.env",
        "calibB_Cerveny.env",
        "MunkB_OneBeam.env",
        "MunkB_ray.env",
        "Gulf_ray_rd.env",
        "free_FGB.env",
        "sbcx_Arr_bin.env",
        "omni.env",
    ] {
        let path = fixture(name);
        if let Err(report) = load_env(&path) {
            panic!("failed to parse {}:\n{report}", path.display());
        }
    }
}

#[test]
#[allow(clippy::float_cmp)]
fn preserves_reference_specific_semantics() {
    let one_beam = load_env(&fixture("MunkB_OneBeam.env")).unwrap().value;
    assert_eq!(
        one_beam.top_options.interpolation,
        SspInterpolation::CubicSpline
    );
    assert_eq!(one_beam.run.kind, RunKind::Coherent);
    assert_eq!(one_beam.trace.selected_launch_angle, Some(80));
    assert_eq!(one_beam.trace.launch_angles_degrees.len(), 100);

    let range_dependent = load_env(&fixture("Gulf_ray_rd.env")).unwrap().value;
    assert_eq!(
        range_dependent.top_options.interpolation,
        SspInterpolation::Quadrilateral
    );
    assert!(range_dependent.bottom_boundary.has_shape_file);

    let volume = load_env(&fixture("free_FGB.env")).unwrap().value;
    assert!(matches!(
        volume.top_options.volume_attenuation,
        VolumeAttenuation::FrancoisGarrison { .. }
    ));

    let binary_arrivals = load_env(&fixture("sbcx_Arr_bin.env")).unwrap().value;
    assert_eq!(binary_arrivals.run.kind, RunKind::Arrivals);

    // A slash-only half-space record retains the medium values from the
    // preceding SSP record in the reference Fortran implementation.
    let free_space = load_env(&fixture("omni.env")).unwrap().value;
    let bellhop::model::BoundaryCondition::AcoustoElastic(bottom) =
        free_space.bottom_boundary.condition
    else {
        panic!("expected an acousto-elastic bottom");
    };
    assert_eq!(bottom.depth_m, 10_000.0);
    assert_eq!(bottom.compressional_speed_mps, 1500.0);
    assert!(free_space.positions.receiver_ranges_m[0] < 0.0);
}
