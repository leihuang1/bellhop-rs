use std::fs;
use std::path::{Path, PathBuf};

use bellhop::legacy::load_case;
use bellhop::model::BeamComponent;
use bellhop::solver::{SimulationLimits, run};

const PRESSURE_TOLERANCE: f32 = 5.0e-8;

#[test]
fn field_modes_match_v2023_5_golden_pressures() {
    for stem in [
        "Field_G",
        "Field_B_incoherent",
        "Field_B_semi",
        "Field_S",
        "Field_CervenyC",
        "Field_CervenyR",
    ] {
        let case = load_case(&fixture(&format!("{stem}.env"))).unwrap().value;
        let result = run(&case, SimulationLimits::default()).unwrap();
        let actual = &result.field_sources[0].samples;
        let expected = read_pressure(&fixture(&format!("{stem}.csv")));
        assert_eq!(actual.len(), expected.len(), "{stem}");

        for (sample_index, (sample, (real, imaginary))) in actual.iter().zip(expected).enumerate() {
            assert!(
                (sample.pressure.re - real).abs() <= PRESSURE_TOLERANCE,
                "{stem} sample {sample_index} real: actual={}, expected={real}",
                sample.pressure.re
            );
            assert!(
                (sample.pressure.im - imaginary).abs() <= PRESSURE_TOLERANCE,
                "{stem} sample {sample_index} imaginary: actual={}, expected={imaginary}",
                sample.pressure.im
            );
        }
    }
}

#[test]
fn cerveny_velocity_components_match_v2023_5_samples() {
    let mut case = load_case(&fixture("Field_CervenyR.env")).unwrap().value;
    let options = case.environment.trace.cerveny.as_mut().unwrap();
    options.component = BeamComponent::Vertical;
    let vertical = run(&case, SimulationLimits::default())
        .unwrap()
        .field_sources[0]
        .samples[16]
        .pressure;
    assert!((vertical.re - 1.022_244_3e-3).abs() <= PRESSURE_TOLERANCE);
    assert!((vertical.im - 2.047_689e-3).abs() <= PRESSURE_TOLERANCE);

    case.environment.trace.cerveny.as_mut().unwrap().component = BeamComponent::Horizontal;
    let horizontal = run(&case, SimulationLimits::default())
        .unwrap()
        .field_sources[0]
        .samples[16]
        .pressure;
    assert!((horizontal.re - 6.562_056_5e-3).abs() <= PRESSURE_TOLERANCE);
    assert!((horizontal.im + 1.197_408_8e-2).abs() <= PRESSURE_TOLERANCE);
}

#[test]
fn field_grid_obeys_the_configured_resource_limit() {
    let case = load_case(&fixture("Field_G.env")).unwrap().value;
    let error = run(
        &case,
        SimulationLimits {
            max_field_cells: 32,
            ..SimulationLimits::default()
        },
    )
    .unwrap_err();

    assert_eq!(error.diagnostics()[0].code, "BH0303");
}

fn read_pressure(path: &Path) -> Vec<(f32, f32)> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let values: Vec<&str> = line.split(',').collect();
            (values[2].parse().unwrap(), values[3].parse().unwrap())
        })
        .collect()
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/golden")
        .join(name)
}
