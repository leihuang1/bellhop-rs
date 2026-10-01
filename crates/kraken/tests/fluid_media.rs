use std::path::Path;

use kraken::{Case, ModeSolver, legacy::load_frequency_cases_from_sources, solve};

fn cases(env: &str) -> Result<Vec<Case>, kraken::DiagnosticReport> {
    load_frequency_cases_from_sources(
        env,
        include_str!("fixtures/Pekeris.flp"),
        Path::new("water.env"),
        Path::new("water.flp"),
        ModeSolver::Krakenc,
    )
}

#[test]
fn biological_loss_is_sampled_at_ssp_nodes_but_not_in_half_spaces() {
    let case = cases(include_str!("fixtures/WaterLossBio.env"))
        .unwrap()
        .remove(0);
    assert!(
        case.water_attenuation_db_per_wavelength
            .iter()
            .all(|&a| a > 0.0)
    );
    // UpdateHSLoss passes HUGE(depth), not the interface depth, to CRCI.
    assert_eq!(
        case.bottom_attenuation_db_per_wavelength.to_bits(),
        0.0_f64.to_bits()
    );
    assert!(
        solve(&case)
            .unwrap()
            .modes
            .modes
            .iter()
            .all(|m| m.attenuation_nepers_per_m > 0.0)
    );
}

#[test]
fn material_loss_preserves_frequency_order_and_validates_public_and_legacy_inputs() {
    let power = cases(include_str!("fixtures/WaterLossPower.env")).unwrap();
    assert_eq!(
        power.iter().map(|c| c.frequency_hz).collect::<Vec<_>>(),
        [75.0, 50.0, 62.5, 50.0]
    );
    assert_eq!(power[1], power[3]);
    let base = power[1].water_attenuation_db_per_wavelength[0];
    let below = (50.0_f64 / 50.0).powf(1.5) / 50.0;
    let above = (75.0_f64 / 50.0) * (60.0_f64 / 50.0).sqrt() / 75.0;
    assert!((power[0].water_attenuation_db_per_wavelength[0] / base - above / below).abs() < 1e-14);
    assert!(
        (power[0].bottom_attenuation_db_per_wavelength
            / power[1].bottom_attenuation_db_per_wavelength
            - (70.0_f64 / 50.0).powf(-0.25))
        .abs()
            < 1e-14
    );

    let base = power[1].clone().into_definition();
    for loss in [
        vec![0.1],
        vec![0.1; 100_001],
        vec![f64::NAN; 3],
        vec![f64::INFINITY; 3],
        vec![-0.1; 3],
        vec![55.0; 3],
    ] {
        let mut invalid = base.clone();
        invalid.water_attenuation_db_per_wavelength = loss;
        assert!(
            Case::from_definition(invalid)
                .unwrap_err()
                .diagnostics()
                .iter()
                .any(|d| d.field == "water_attenuation_db_per_wavelength")
        );
    }
    let env = include_str!("fixtures/WaterLossBio.env");
    for invalid in [
        env.replacen("\n2\n", "\n201\n", 1),
        env.replace("500.0 2.5 0.25", "500.0 0.0 0.25"),
        env.replace("600.0 3.0 0.15", "-600.0 3.0 0.15"),
        env.replace("43.0 120.0", "120.0 43.0"),
        env.replace("0.0 80.0", "NaN 80.0"),
    ] {
        let report = cases(&invalid).unwrap_err();
        assert!(
            report
                .diagnostics()
                .iter()
                .all(|d| d.path == Path::new("water.env") && d.line > 1)
        );
    }
    for invalid in [
        include_str!("fixtures/WaterLossFg.env").replace("10.0 35.0", "-273.0 35.0"),
        include_str!("fixtures/WaterLossFg.env").replace("10.0 35.0", "10.0 -35.0"),
        include_str!("fixtures/WaterLossPower.env").replace("1.5 60.0", "1.5 0.0"),
        include_str!("fixtures/WaterLossC.env")
            .replace("43.0 1530.0 0.0 1.0", "43.0 1530.0 0.0 1.1"),
    ] {
        assert!(cases(&invalid).is_err());
    }

    let mut undershoot = base;
    undershoot.interpolation = kraken::Interpolation::Spline;
    undershoot.sound_speed_profile[1].depth_m = 1.0;
    undershoot.water_attenuation_db_per_wavelength = vec![0.0, 0.0, 30.0];
    assert!(
        Case::from_definition(undershoot)
            .unwrap_err()
            .diagnostics()
            .iter()
            .any(|d| d.field == "water_attenuation_db_per_wavelength")
    );

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let mut analytic = kraken::legacy::load_complex_case(
        root.join("MunkAnalyticComplex.env"),
        root.join("MunkAnalyticComplex.flp"),
    )
    .unwrap()
    .into_definition();
    analytic.water_attenuation_db_per_wavelength = vec![0.1];
    assert!(Case::from_definition(analytic).is_err());
    let mut table =
        kraken::legacy::load_complex_case(root.join("TabRefBrcN.env"), root.join("TabRefBrcN.flp"))
            .unwrap()
            .into_definition();
    table.water_attenuation_db_per_wavelength = vec![0.1; table.sound_speed_profile.len()];
    assert!(
        Case::from_definition(table)
            .unwrap_err()
            .diagnostics()
            .iter()
            .any(|d| d.field == "bottom_boundary")
    );
}
