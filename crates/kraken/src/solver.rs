use crate::{
    Case, Diagnostic, DiagnosticReport, MAX_FIELD_SAMPLES, ModeSet, PressureField,
    SimulationResult, SourceGeometry,
};
use num_complex::{Complex32, Complex64};
use std::path::Path;

// ponytail: 250m contributions cover original sductK's 217m; revisit with measured throughput.
const MAX_FIELD_WORK: usize = 250_000_000;

pub(super) fn solve(case: &Case) -> Result<SimulationResult, DiagnosticReport> {
    let modes = match case.mode_solver {
        crate::ModeSolver::Kraken => crate::modes::solve(case)?,
        crate::ModeSolver::Krakenc => crate::complex_modes::solve(case)?,
    };
    let field = synthesize_field(case, &modes)?;
    Ok(SimulationResult { modes, field })
}

fn synthesize_field(case: &Case, mode_set: &ModeSet) -> Result<PressureField, DiagnosticReport> {
    let mode_count = mode_set.modes.len().min(case.mode_limit);
    let sample_count = case
        .source_depths_m
        .len()
        .checked_mul(case.receiver_depths_m.len())
        .and_then(|count| count.checked_mul(case.receiver_ranges_m.len()))
        .ok_or_else(|| error("KR0302", "field grid size overflow", "field_grid"))?;
    if sample_count > MAX_FIELD_SAMPLES
        || sample_count
            .checked_mul(mode_count)
            .is_none_or(|work| work > MAX_FIELD_WORK)
    {
        return Err(error(
            "KR0302",
            format!("field calculation exceeds the {MAX_FIELD_WORK} modal-operation limit"),
            "field_grid",
        ));
    }
    // EvaluateMod.f90: single-precision modal arithmetic and separate range/offset phases.
    #[allow(clippy::approx_constant)]
    let field_pi = 3.141_592_6_f32;
    let factor = Complex32::new(0.0, 1.0)
        * (2.0 * field_pi).sqrt()
        * Complex32::from_polar(1.0, field_pi * 0.25);
    let mut pressure = Vec::with_capacity(sample_count);
    for &source_depth in &case.source_depths_m {
        let source_shapes: Vec<_> = mode_set
            .modes
            .iter()
            .take(mode_count)
            .map(|mode| {
                sample_shape(
                    &mode_set.sampled_depths_m,
                    &mode.eigenfunction,
                    source_depth,
                )
            })
            .collect();
        for (receiver_index, &receiver_depth) in case.receiver_depths_m.iter().enumerate() {
            let receiver_offset = case.receiver_offsets_m[receiver_index];
            let receiver_shapes: Vec<_> = mode_set
                .modes
                .iter()
                .take(mode_count)
                .map(|mode| {
                    sample_shape(
                        &mode_set.sampled_depths_m,
                        &mode.eigenfunction,
                        receiver_depth,
                    )
                })
                .collect();
            for &range in &case.receiver_ranges_m {
                let mut value = Complex32::new(0.0, 0.0);
                for ((mode, source_shape), receiver_shape) in mode_set
                    .modes
                    .iter()
                    .take(mode_count)
                    .zip(&source_shapes)
                    .zip(&receiver_shapes)
                {
                    let k = single(mode.horizontal_wavenumber_rad_per_m);
                    let amplitude = factor * single(*source_shape)
                        / if case.source_geometry == SourceGeometry::Point {
                            k.sqrt()
                        } else {
                            k
                        };
                    let ik = double(Complex32::new(0.0, -1.0) * k);
                    let offset_shape = single(
                        double(amplitude * single(*receiver_shape)) * (ik * receiver_offset).exp(),
                    );
                    value += offset_shape * single((ik * range).exp());
                }
                if case.source_geometry == SourceGeometry::Point && range + receiver_offset > 0.0 {
                    #[allow(clippy::cast_possible_truncation)]
                    let spreading = (range + receiver_offset).sqrt() as f32;
                    value /= spreading;
                }
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(error(
                        "KR0302",
                        "field pressure is not finite",
                        "field_grid",
                    ));
                }
                pressure.push(double(value));
            }
        }
    }
    Ok(PressureField {
        source_depths_m: case.source_depths_m.clone(),
        receiver_depths_m: case.receiver_depths_m.clone(),
        receiver_ranges_m: case.receiver_ranges_m.clone(),
        receiver_offsets_m: case.receiver_offsets_m.clone(),
        pressure,
    })
}

fn double(value: Complex32) -> Complex64 {
    Complex64::new(f64::from(value.re), f64::from(value.im))
}

#[allow(clippy::cast_possible_truncation)]
fn single(value: Complex64) -> Complex32 {
    Complex32::new(value.re as f32, value.im as f32)
}

fn sample_shape(depths: &[f64], values: &[Complex64], depth: f64) -> Complex64 {
    if depth <= depths[0] {
        return values[0];
    }
    let index = depths.partition_point(|sample| *sample < depth);
    if index == depths.len() {
        return values[values.len() - 1];
    }
    let weight = (depth - depths[index - 1]) / (depths[index] - depths[index - 1]);
    values[index - 1] + (values[index] - values[index - 1]) * weight
}

pub(super) fn error(
    code: &'static str,
    message: impl Into<String>,
    field: &str,
) -> DiagnosticReport {
    DiagnosticReport::one(Diagnostic::new(
        code,
        message,
        field,
        Path::new("<case>"),
        1,
        1,
    ))
}

#[cfg(test)]
#[test]
fn field_work_is_bounded_before_modal_products() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/Pekeris");
    let mut input =
        crate::legacy::load_case(root.with_extension("env"), root.with_extension("flp"))
            .unwrap()
            .into_definition();
    input.source_depths_m = vec![75.0; 100_000];
    input.receiver_depths_m = vec![75.0; 5];
    input.receiver_offsets_m = vec![0.0; 5];
    input.receiver_ranges_m = vec![500.0, 1000.0];
    let case = Case::from_definition(input).unwrap();
    let mode = crate::NormalMode {
        horizontal_wavenumber_rad_per_m: Complex64::new(1.0, 0.0),
        phase_speed_mps: 1.0,
        group_speed_mps: 1.0,
        attenuation_nepers_per_m: 0.0,
        eigenfunction: vec![],
    };
    let modes = ModeSet {
        frequency_hz: case.frequency_hz,
        sampled_depths_m: case.mode_sample_depths_m.clone(),
        modes: vec![mode; 251],
    };
    assert_eq!(
        synthesize_field(&case, &modes).unwrap_err().diagnostics()[0].field,
        "field_grid"
    );
}
