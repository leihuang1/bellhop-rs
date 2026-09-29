//! Test-only reader for pinned little-endian, single-frequency fluid .mod/.shd files.
//! These binary formats are reference artifacts, not a production output contract.
use std::fs;
use std::path::{Path, PathBuf};

use kraken::{
    Case, ModeSet, PressureField,
    legacy::{load_case, load_complex_case},
    solve, solve_complex_modes,
};
use num_complex::Complex64;

const CASES: &[&str] = &[
    "Pekeris",
    "PekerisFiltered",
    "PekerisDense",
    "PekerisDenseLoss",
    "PekerisRefined",
    "PekerisSpline3",
    "PekerisRigid",
    "PekerisRigidLoss",
    "PekerisHard",
    "PekerisHardBoth",
    "PekerisRigidPlane",
    "MunkLossless",
    "MunkBottomLoss",
    "MunkAnalytic",
    "SductTrapped",
    "SductPchip",
    "SductSpline",
];

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

#[test]
fn fluid_modes_and_field_match_pinned_goldens() {
    for name in CASES {
        compare(
            &fixtures().join(name).with_extension("env"),
            &fixtures().join("golden").join(name),
        );
    }
}

#[test]
fn complex_fluid_modes_match_pinned_fortran() {
    for (name, golden) in [
        ("PekerisComplex", "PekerisComplex"),
        ("PekerisComplexBlank", "PekerisComplex"),
        ("PekerisComplexSlow", "PekerisComplexSlow"),
        ("PekerisComplexCLow", "PekerisComplexCLow"),
        ("PekerisComplexGradient", "PekerisComplexGradient"),
        (
            "PekerisComplexReverseGradient",
            "PekerisComplexReverseGradient",
        ),
        ("MunkLeakyPartial", "MunkLeakyPartial"),
    ] {
        compare_complex(
            &fixtures().join(format!("{name}.env")),
            &fixtures().join("golden").join(golden),
        );
    }
}

#[test]
#[ignore = "requires a pinned external Fortran reference run"]
fn complex_fluid_matches_fresh_reference() {
    let name = std::env::var("KRAKEN_COMPLEX_CASE").expect("KRAKEN_COMPLEX_CASE is required");
    let reference = PathBuf::from(
        std::env::var_os("KRAKEN_COMPLEX_REFERENCE_ROOT")
            .expect("KRAKEN_COMPLEX_REFERENCE_ROOT is required"),
    );
    compare_complex(&fixtures().join(format!("{name}.env")), &reference);
}

fn compare_complex(env: &Path, reference: &Path) {
    let case = load_complex_case(env, fixtures().join("Pekeris.flp")).unwrap();
    let actual = solve_complex_modes(&case).unwrap();
    let errors = compare_modes(
        &case,
        &actual,
        &Records::read(&reference.with_extension("mod")),
        &fs::read_to_string(reference.with_extension("prt")).unwrap(),
    );
    eprintln!(
        "{}: {} KRAKENC modes, errors {errors:?}",
        env.display(),
        actual.modes.len()
    );
}

#[test]
#[ignore = "requires a pinned external Fortran reference run"]
fn fluid_matches_fresh_reference() {
    let env =
        std::env::var_os("KRAKEN_DIFFERENTIAL_ENV").expect("KRAKEN_DIFFERENTIAL_ENV is required");
    let root =
        std::env::var_os("KRAKEN_DIFFERENTIAL_ROOT").expect("KRAKEN_DIFFERENTIAL_ROOT is required");
    compare(Path::new(&env), Path::new(&root));
}

fn compare(env: &Path, reference: &Path) {
    let case = load_case(env, env.with_extension("flp")).unwrap();
    let actual = solve(&case).unwrap();
    let modes = Records::read(&reference.with_extension("mod"));
    let field = Records::read(&reference.with_extension("shd"));
    let printed = fs::read_to_string(reference.with_extension("prt")).unwrap();
    let (k_error, loss_error, binary_loss_error, shape_error) =
        compare_modes(&case, &actual.modes, &modes, &printed);
    let pressure_error = compare_field(&case, &actual.field, &field);
    eprintln!(
        "{}: {} modes, {} pressure samples; max |dk|={k_error:e}, |dalpha(.prt)|={loss_error:e}, |dalpha(.mod)|={binary_loss_error:e}, |dphi|={shape_error:e}, |dp|={pressure_error:e}",
        env.display(),
        actual.modes.modes.len(),
        actual.field.pressure.len()
    );
}

fn close(actual: f64, expected: f64, tolerance: f64, field: &str) -> f64 {
    let error = (actual - expected).abs();
    assert!(
        actual.is_finite() && expected.is_finite() && error <= tolerance,
        "{field}: actual {actual:e}, reference {expected:e}, error {error:e} > {tolerance:e}"
    );
    error
}

fn complex_close(actual: Complex64, expected: Complex64, tolerance: f64, field: &str) -> f64 {
    close((actual - expected).norm(), 0.0, tolerance, field)
}

#[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
fn compare_modes(
    case: &Case,
    actual: &ModeSet,
    file: &Records,
    printed: &str,
) -> (f64, f64, f64, f64) {
    let header = file.record(0);
    assert_eq!(count(header, 84), 1, "mode frequency count");
    assert_eq!(count(header, 88), 1, "medium count");
    assert_eq!(
        count(header, 92),
        actual.sampled_depths_m.len(),
        "mode depth count"
    );
    assert_eq!(
        count(header, 96),
        actual.sampled_depths_m.len(),
        "fluid mode shape size"
    );
    assert_eq!(&file.record(1)[4..12], b"ACOUSTIC", "fluid material");
    close(
        double(file.record(3), 0),
        case.frequency_hz,
        1e-12,
        "mode frequency",
    );
    for (index, &depth) in actual.sampled_depths_m.iter().enumerate() {
        close(
            f64::from(depth as f32),
            single(file.record(4), 4 * index),
            0.0,
            "mode depth",
        );
    }
    let mode_count = count(file.record(5), 0);
    assert_eq!(mode_count, actual.modes.len(), "mode count");
    let modes_per_record = file.record_bytes / 8;
    assert_eq!(
        file.len(),
        7 + mode_count + mode_count.div_ceil(modes_per_record),
        "mode file record count"
    );

    // .prt preserves extrapolated k to ten decimal places; .mod stores complex32.
    let (_, table) = printed
        .split_once("Group Speed\n")
        .expect("modal print table");
    let rows: Vec<_> = table
        .lines()
        .skip(1)
        .take_while(|line| {
            line.split_whitespace()
                .next()
                .is_some_and(|word| word.parse::<usize>().is_ok())
        })
        .collect();
    let print_stride = (mode_count / 30).max(1);
    assert_eq!(
        rows.len(),
        mode_count.div_ceil(print_stride),
        "printed mode count"
    );
    let mut max_k = 0.0_f64;
    let mut max_loss = 0.0_f64;
    let mut max_binary_loss = 0.0_f64;
    let mut max_shape = 0.0_f64;
    for (row_index, row) in rows.into_iter().enumerate() {
        let index = row_index * print_stride;
        let mode = &actual.modes[index];
        let columns: Vec<_> = row.split_whitespace().collect();
        assert_eq!(columns.len(), 5, "modal print columns");
        assert_eq!(columns[0].parse::<usize>().unwrap(), index + 1);
        let values: Vec<f64> = columns[1..]
            .iter()
            .map(|value| value.parse().unwrap())
            .collect();
        max_k = max_k.max(close(
            mode.horizontal_wavenumber_rad_per_m.re,
            values[0],
            5e-10,
            "wavenumber (.prt)",
        ));
        // Fortran prints signed Im(k) at G10.2; the public attenuation is -Im(k).
        // Check to half of the last printed decimal place, not an arbitrary loose tolerance.
        let attenuation_tolerance = if mode.attenuation_nepers_per_m == 0.0 {
            1e-12
        } else {
            let (mantissa, exponent) = columns[2].split_once(['E', 'e']).unwrap();
            let digits = mantissa.split_once('.').unwrap().1.len();
            0.5 * 10_f64.powi(exponent.parse::<i32>().unwrap() - i32::try_from(digits).unwrap())
        };
        max_loss = max_loss.max(close(
            -mode.attenuation_nepers_per_m,
            values[1],
            attenuation_tolerance.max(1e-12),
            "attenuation (.prt)",
        ));
        close(mode.phase_speed_mps, values[2], 5e-6, "phase speed");
        close(mode.group_speed_mps, values[3], 0.005, "group speed");
    }
    for (index, mode) in actual.modes.iter().enumerate() {
        let k_record = file.record(7 + mode_count + index / modes_per_record);
        let stored = complex(k_record, 8 * (index % modes_per_record));
        close(
            f64::from(mode.horizontal_wavenumber_rad_per_m.re as f32),
            stored.re,
            0.0,
            &format!("real wavenumber (.mod) mode {}", index + 1),
        );
        // Loss depends on the first-mesh root; a different root rounding can
        // shift small imaginary parts even when the extrapolated real k agrees.
        max_binary_loss = max_binary_loss.max(close(
            f64::from(mode.horizontal_wavenumber_rad_per_m.im as f32),
            stored.im,
            if case.bottom_attenuation_db_per_wavelength > 0.0 {
                1e-10
            } else {
                0.0
            },
            "imaginary wavenumber (.mod)",
        ));
        let shape = file.record(7 + index);
        assert_eq!(mode.eigenfunction.len(), actual.sampled_depths_m.len());
        let reference: Vec<_> = (0..mode.eigenfunction.len())
            .map(|i| complex(shape, 8 * i))
            .collect();
        // Align only the arbitrary unit phase, never rescale amplitude to hide normalization errors.
        let overlap: Complex64 = reference
            .iter()
            .zip(&mode.eigenfunction)
            .map(|(a, b)| a.conj() * b)
            .sum();
        assert!(
            overlap.norm().is_finite() && overlap.norm() > 0.0,
            "mode phase cannot be aligned"
        );
        let phase = overlap / overlap.norm();
        for (&value, expected) in mode.eigenfunction.iter().zip(reference) {
            max_shape = max_shape.max(complex_close(value, expected * phase, 1e-6, "mode shape"));
        }
    }
    (max_k, max_loss, max_binary_loss, max_shape)
}

#[allow(clippy::cast_possible_truncation)]
fn compare_field(case: &Case, actual: &PressureField, file: &Records) -> f64 {
    assert_eq!(&file.record(1)[..10], b"          ", "FIELD plot type");
    let header = file.record(2);
    for offset in [0, 4, 8, 12] {
        assert_eq!(count(header, offset), 1, "frequency/bearing/x/y count");
    }
    assert_eq!(count(header, 16), actual.source_depths_m.len());
    assert_eq!(count(header, 20), actual.receiver_depths_m.len());
    assert_eq!(count(header, 24), actual.receiver_ranges_m.len());
    // FIELD leaves header freq0 unset; freqVec is the authoritative frequency.
    close(
        double(file.record(3), 0),
        case.frequency_hz,
        1e-12,
        "field frequency",
    );
    for (record, depths) in [(7, &actual.source_depths_m), (8, &actual.receiver_depths_m)] {
        for (index, &depth) in depths.iter().enumerate() {
            close(
                f64::from(depth as f32),
                single(file.record(record), 4 * index),
                0.0,
                "field depth",
            );
        }
    }
    for (index, &range) in actual.receiver_ranges_m.iter().enumerate() {
        close(
            range,
            double(file.record(9), 8 * index),
            1e-9,
            "field range",
        );
    }
    assert_eq!(actual.receiver_offsets_m, case.receiver_offsets_m);
    let rows = actual.source_depths_m.len() * actual.receiver_depths_m.len();
    let ranges = actual.receiver_ranges_m.len();
    assert_eq!(file.len(), 10 + rows, "field record count");
    assert_eq!(actual.pressure.len(), rows * ranges);
    let mut maximum = 0.0_f64;
    for (index, &pressure) in actual.pressure.iter().enumerate() {
        let expected = complex(file.record(10 + index / ranges), 8 * (index % ranges));
        maximum = maximum.max(complex_close(
            pressure,
            expected,
            2e-6,
            &format!("pressure[{index}]"),
        ));
    }
    maximum
}

#[derive(Clone)]
struct Records {
    bytes: Vec<u8>,
    record_bytes: usize,
}

impl Records {
    fn read(path: &Path) -> Self {
        let bytes = fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let words = count(&bytes, 0);
        assert!((25..=1_000_000).contains(&words), "invalid record size");
        let record_bytes = words * 4;
        assert_eq!(bytes.len() % record_bytes, 0, "truncated reference file");
        Self {
            bytes,
            record_bytes,
        }
    }

    fn len(&self) -> usize {
        self.bytes.len() / self.record_bytes
    }

    fn record(&self, index: usize) -> &[u8] {
        self.bytes
            .get(index * self.record_bytes..(index + 1) * self.record_bytes)
            .expect("missing reference record")
    }
}

fn count(bytes: &[u8], offset: usize) -> usize {
    usize::try_from(i32::from_le_bytes(
        bytes[offset..offset + 4].try_into().unwrap(),
    ))
    .expect("negative count")
}

fn single(bytes: &[u8], offset: usize) -> f64 {
    f64::from(f32::from_le_bytes(
        bytes[offset..offset + 4].try_into().unwrap(),
    ))
}

fn double(bytes: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn complex(bytes: &[u8], offset: usize) -> Complex64 {
    Complex64::new(single(bytes, offset), single(bytes, offset + 4))
}

#[test]
fn comparator_detects_corruption_but_accepts_eigenvector_sign_changes() {
    let root = fixtures().join("Pekeris");
    let case = load_case(root.with_extension("env"), root.with_extension("flp")).unwrap();
    let result = solve(&case).unwrap();
    let root = fixtures().join("golden/Pekeris");
    let modes = Records::read(&root.with_extension("mod"));
    let printed = fs::read_to_string(root.with_extension("prt")).unwrap();
    let field = Records::read(&root.with_extension("shd"));

    let mut flipped = modes.clone();
    for i in 0..case.mode_sample_depths_m.len() {
        let offset = 7 * flipped.record_bytes + 8 * i;
        let value = f32::from_le_bytes(flipped.bytes[offset..offset + 4].try_into().unwrap());
        flipped.bytes[offset..offset + 4].copy_from_slice(&(-value).to_le_bytes());
    }
    compare_modes(&case, &result.modes, &flipped, &printed);
    let lossy = fixtures().join("MunkBottomLoss");
    let lossy_case = load_case(lossy.with_extension("env"), lossy.with_extension("flp")).unwrap();
    let lossy_modes = Records::read(&fixtures().join("golden/MunkBottomLoss.mod"));
    let lossy_printed = fs::read_to_string(fixtures().join("golden/MunkBottomLoss.prt")).unwrap();
    let lossy_result = solve(&lossy_case).unwrap();
    let mut corrupted_loss = lossy_modes.clone();
    let imag = (7 + lossy_result.modes.modes.len()) * corrupted_loss.record_bytes + 4;
    corrupted_loss.bytes[imag..imag + 4].copy_from_slice(&1e-7_f32.to_le_bytes());
    assert!(
        std::panic::catch_unwind(|| compare_modes(
            &lossy_case,
            &lossy_result.modes,
            &corrupted_loss,
            &lossy_printed
        ))
        .is_err()
    );

    for value in [1.0_f32, f32::NAN] {
        for record in [4, 7, 7 + result.modes.modes.len()] {
            let mut corrupted = modes.clone();
            let offset = record * corrupted.record_bytes;
            corrupted.bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
            assert!(
                std::panic::catch_unwind(|| compare_modes(
                    &case,
                    &result.modes,
                    &corrupted,
                    &printed
                ))
                .is_err()
            );
        }
        let mut corrupted = field.clone();
        let offset = 10 * corrupted.record_bytes;
        corrupted.bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
        assert!(
            std::panic::catch_unwind(|| compare_field(&case, &result.field, &corrupted)).is_err()
        );
    }
}
