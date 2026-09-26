use std::path::PathBuf;

use kraken::{legacy::load_case, solve};
use num_complex::Complex64;

// Extracted from the pinned Fortran .prt, .mod, and .shd outputs.
const MODE_K: [f64; 3] = [0.207_648_652_2, 0.202_104_615_5, 0.192_450_526_0];
const MODE_GROUP_SPEED: [f64; 3] = [1490.38, 1459.18, 1406.66];
const MODE_SHAPES: [[f64; 3]; 3] = [
    [0.083_455_719_05, 0.117_315_866_05, 0.055_813_901_13],
    [0.130_341_440_44, -0.110_319_383_44, -0.099_323_667_59],
    [0.116_119_377_32, -0.011_325_890_20, 0.124_988_913_54],
];
const FIELD: [[f64; 2]; 9] = [
    [0.100_606_732_07, -0.273_737_311_36],
    [-0.099_672_563_37, -0.050_616_275_52],
    [-0.205_574_810_50, 0.057_186_506_69],
    [0.060_523_387_04, 0.012_440_915_22],
    [-0.038_632_467_39, 0.291_528_135_54],
    [0.117_217_414_08, 0.205_000_877_38],
    [-0.005_802_171_30, 0.060_673_333_70],
    [-0.001_534_514_13, 0.219_899_296_76],
    [0.103_178_843_86, 0.101_857_766_51],
];

#[test]
fn pekeris_modes_and_field_match_pinned_reference() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/Pekeris");
    let case = load_case(fixture.with_extension("env"), fixture.with_extension("flp")).unwrap();
    let result = solve(&case).unwrap();

    assert_eq!(result.modes.sampled_depths_m, [25.0, 75.0, 99.0]);
    assert_eq!(result.modes.modes.len(), MODE_K.len());
    for (index, mode) in result.modes.modes.iter().enumerate() {
        assert!((mode.horizontal_wavenumber_rad_per_m.re - MODE_K[index]).abs() < 5e-10);
        assert!(mode.horizontal_wavenumber_rad_per_m.im.abs() < 1e-12);
        assert!((mode.group_speed_mps - MODE_GROUP_SPEED[index]).abs() < 0.005);
        for (actual, expected) in mode.eigenfunction.iter().zip(MODE_SHAPES[index]) {
            assert!((*actual - Complex64::new(expected, 0.0)).norm() < 1e-6);
        }
    }

    assert_eq!(result.field.receiver_offsets_m, [0.0; 3]);
    assert_eq!(result.field.pressure.len(), FIELD.len());
    for (actual, [real, imaginary]) in result.field.pressure.iter().zip(FIELD) {
        assert!((*actual - Complex64::new(real, imaginary)).norm() < 2e-6);
    }
}
