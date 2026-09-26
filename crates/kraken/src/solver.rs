use std::f64::consts::PI;
use std::path::Path;

use num_complex::Complex64;

use crate::{
    Case, Diagnostic, DiagnosticReport, MAX_FIELD_SAMPLES, ModeSet, NormalMode, PressureField,
    SimulationResult,
};

const MAX_MODES: u32 = 20_000;
const MAX_MODE_SHAPE_VALUES: usize = 5_000_000;
// ponytail: fixed work ceiling keeps eager fields bounded; raise when a real case needs it.
const MAX_FIELD_WORK: usize = 50_000_000;
const ROOT_ITERATIONS: usize = 100;

// ponytail: closed-form homogeneous Pekeris only; replace when layered profiles land.
pub(super) fn solve(case: &Case) -> Result<SimulationResult, DiagnosticReport> {
    let definition = &case.0;
    let omega = 2.0 * PI * definition.frequency_hz;
    let water_k = omega / definition.water_sound_speed_mps;
    let bottom_k = omega / definition.bottom_sound_speed_mps;
    let cutoff_q_squared = water_k.mul_add(water_k, -bottom_k * bottom_k);
    let q_cutoff = cutoff_q_squared.sqrt();
    if !omega.is_finite() || !q_cutoff.is_finite() || q_cutoff <= 0.0 {
        return Err(error(
            "KR0302",
            "frequency and waveguide values exceed the supported numeric range",
            "frequency_hz",
        ));
    }
    let slow_phase_speed = definition.c_low_mps.max(definition.water_sound_speed_mps);
    let fast_phase_speed = definition.c_high_mps.min(definition.bottom_sound_speed_mps);
    if slow_phase_speed >= fast_phase_speed {
        return Err(error(
            "KR0301",
            "phase-speed limits contain no trapped Pekeris modes",
            "phase_speed_limits",
        ));
    }

    let q_min = q_for_phase_speed(omega, definition.water_sound_speed_mps, slow_phase_speed);
    let q_max =
        q_for_phase_speed(omega, definition.water_sound_speed_mps, fast_phase_speed).min(q_cutoff);
    if !q_min.is_finite() || !q_max.is_finite() {
        return Err(error(
            "KR0302",
            "phase-speed limits exceed the supported numeric range",
            "phase_speed_limits",
        ));
    }
    let roots = find_roots(
        q_min,
        q_max,
        q_cutoff,
        definition.water_depth_m,
        definition.water_density_g_cm3 / definition.bottom_density_g_cm3,
    )?;
    if roots.is_empty() {
        return Err(error(
            "KR0301",
            "no trapped modes fall inside the phase-speed limits",
            "phase_speed_limits",
        ));
    }
    if roots
        .len()
        .checked_mul(definition.mode_sample_depths_m.len())
        .is_none_or(|count| count > MAX_MODE_SHAPE_VALUES)
    {
        return Err(error(
            "KR0302",
            format!("mode shapes exceed the {MAX_MODE_SHAPE_VALUES} sample limit"),
            "mode_sample_depths_m",
        ));
    }

    let modes = roots
        .into_iter()
        .map(|q| build_mode(case, q, omega))
        .collect::<Vec<_>>();
    if modes.iter().any(|mode| {
        !mode.horizontal_wavenumber_rad_per_m.re.is_finite()
            || !mode.horizontal_wavenumber_rad_per_m.im.is_finite()
            || !mode.phase_speed_mps.is_finite()
            || !mode.group_speed_mps.is_finite()
            || !mode.attenuation_nepers_per_m.is_finite()
            || mode
                .eigenfunction
                .iter()
                .any(|value| !value.re.is_finite() || !value.im.is_finite())
    }) {
        return Err(error(
            "KR0302",
            "mode calculation exceeded the supported numeric range",
            "modes",
        ));
    }
    let mode_set = ModeSet {
        frequency_hz: definition.frequency_hz,
        sampled_depths_m: definition.mode_sample_depths_m.clone(),
        modes,
    };
    let field = synthesize_field(case, &mode_set)?;
    Ok(SimulationResult {
        modes: mode_set,
        field,
    })
}

fn q_for_phase_speed(omega: f64, water_speed: f64, phase_speed: f64) -> f64 {
    ((omega / water_speed).powi(2) - (omega / phase_speed).powi(2))
        .max(0.0)
        .sqrt()
}

fn find_roots(
    q_min: f64,
    q_max: f64,
    q_cutoff: f64,
    water_depth: f64,
    density_ratio: f64,
) -> Result<Vec<f64>, DiagnosticReport> {
    let mut roots = Vec::new();
    for interval in 0..=MAX_MODES {
        let start = (f64::from(interval) + 0.5) * PI / water_depth;
        if start >= q_max {
            break;
        }
        if interval == MAX_MODES {
            return Err(error(
                "KR0302",
                format!("mode count exceeds the {MAX_MODES} mode limit"),
                "mode_count",
            ));
        }
        let end = (f64::from(interval) + 1.0) * PI / water_depth;
        let lower = start.max(q_min);
        let upper = end.min(q_max);
        if lower >= upper {
            continue;
        }
        let lower_value = dispersion(lower, q_cutoff, water_depth, density_ratio);
        let upper_value = dispersion(upper, q_cutoff, water_depth, density_ratio);
        if lower_value == 0.0 {
            roots.push(lower);
        } else if upper_value == 0.0 {
            let cutoff_margin = 16.0 * f64::EPSILON * q_cutoff.max(1.0);
            if q_cutoff - upper > cutoff_margin {
                roots.push(upper);
            }
        } else if lower_value.is_sign_positive() != upper_value.is_sign_positive() {
            roots.push(bisect_root(
                lower,
                upper,
                lower_value,
                q_cutoff,
                water_depth,
                density_ratio,
            ));
        }
    }
    Ok(roots)
}

fn dispersion(q: f64, q_cutoff: f64, water_depth: f64, density_ratio: f64) -> f64 {
    let gamma = (q_cutoff.mul_add(q_cutoff, -q * q)).max(0.0).sqrt();
    let qh = q * water_depth;
    qh.cos() + density_ratio * gamma / q * qh.sin()
}

fn bisect_root(
    mut lower: f64,
    mut upper: f64,
    mut lower_value: f64,
    q_cutoff: f64,
    water_depth: f64,
    density_ratio: f64,
) -> f64 {
    for _ in 0..ROOT_ITERATIONS {
        let middle = lower + (upper - lower) * 0.5;
        let middle_value = dispersion(middle, q_cutoff, water_depth, density_ratio);
        if middle_value == 0.0 || upper - lower <= 4.0 * f64::EPSILON * middle.abs().max(1.0) {
            return middle;
        }
        if lower_value.is_sign_positive() == middle_value.is_sign_positive() {
            lower = middle;
            lower_value = middle_value;
        } else {
            upper = middle;
        }
    }
    lower + (upper - lower) * 0.5
}

fn build_mode(case: &Case, q: f64, omega: f64) -> NormalMode {
    let definition = &case.0;
    let k_real = ((omega / definition.water_sound_speed_mps).powi(2) - q * q).sqrt();
    let k = Complex64::new(k_real, 0.0);
    let gamma = (k_real * k_real - (omega / definition.bottom_sound_speed_mps).powi(2)).sqrt();
    let water_depth = definition.water_depth_m;
    let interface_shape = (q * water_depth).sin();
    let water_norm = (water_depth * 0.5 - (2.0 * q * water_depth).sin() / (4.0 * q))
        / definition.water_density_g_cm3;
    let bottom_norm = interface_shape.powi(2) / (2.0 * gamma * definition.bottom_density_g_cm3);
    let normalization = (water_norm + bottom_norm).sqrt().recip();
    let eigenfunction = definition
        .mode_sample_depths_m
        .iter()
        .map(|depth| Complex64::new(normalization * (q * depth).sin(), 0.0))
        .collect();
    let group_integral = water_norm / definition.water_sound_speed_mps.powi(2)
        + bottom_norm / definition.bottom_sound_speed_mps.powi(2);

    NormalMode {
        horizontal_wavenumber_rad_per_m: k,
        phase_speed_mps: omega / k_real,
        group_speed_mps: k_real * (water_norm + bottom_norm) / (omega * group_integral),
        attenuation_nepers_per_m: 0.0,
        eigenfunction,
    }
}

fn synthesize_field(case: &Case, mode_set: &ModeSet) -> Result<PressureField, DiagnosticReport> {
    let definition = &case.0;
    let mode_count = mode_set.modes.len().min(definition.mode_limit);
    let sample_count = definition
        .source_depths_m
        .len()
        .checked_mul(definition.receiver_depths_m.len())
        .and_then(|count| count.checked_mul(definition.receiver_ranges_m.len()))
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

    let factor =
        Complex64::new(0.0, 1.0) * (2.0 * PI).sqrt() * Complex64::from_polar(1.0, PI * 0.25);
    let mut pressure = Vec::with_capacity(sample_count);
    for source_depth in &definition.source_depths_m {
        let source_shapes: Vec<_> = mode_set
            .modes
            .iter()
            .take(mode_count)
            .map(|mode| {
                sample_shape(
                    &mode_set.sampled_depths_m,
                    &mode.eigenfunction,
                    *source_depth,
                )
            })
            .collect();
        for (receiver_index, receiver_depth) in definition.receiver_depths_m.iter().enumerate() {
            let receiver_offset = definition.receiver_offsets_m[receiver_index];
            let receiver_shapes: Vec<_> = mode_set
                .modes
                .iter()
                .take(mode_count)
                .map(|mode| {
                    sample_shape(
                        &mode_set.sampled_depths_m,
                        &mode.eigenfunction,
                        *receiver_depth,
                    )
                })
                .collect();
            for range in &definition.receiver_ranges_m {
                let mut value = Complex64::new(0.0, 0.0);
                for ((mode, source_shape), receiver_shape) in mode_set
                    .modes
                    .iter()
                    .take(mode_count)
                    .zip(&source_shapes)
                    .zip(&receiver_shapes)
                {
                    let amplitude = factor * *source_shape * *receiver_shape
                        / mode.horizontal_wavenumber_rad_per_m;
                    let phase_argument = Complex64::new(0.0, -(range + receiver_offset))
                        * mode.horizontal_wavenumber_rad_per_m;
                    if !phase_argument.re.is_finite() || !phase_argument.im.is_finite() {
                        return Err(error(
                            "KR0302",
                            "field phase exceeds the supported numeric range",
                            "field_grid",
                        ));
                    }
                    value += amplitude * phase_argument.exp();
                }
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(error(
                        "KR0302",
                        "field pressure is not finite",
                        "field_grid",
                    ));
                }
                pressure.push(value);
            }
        }
    }

    Ok(PressureField {
        source_depths_m: definition.source_depths_m.clone(),
        receiver_depths_m: definition.receiver_depths_m.clone(),
        receiver_ranges_m: definition.receiver_ranges_m.clone(),
        receiver_offsets_m: definition.receiver_offsets_m.clone(),
        pressure,
    })
}

fn sample_shape(depths: &[f64], values: &[Complex64], depth: f64) -> Complex64 {
    if depth <= depths[0] {
        return values[0];
    }
    for index in 1..depths.len() {
        if depth <= depths[index] {
            let weight = (depth - depths[index - 1]) / (depths[index] - depths[index - 1]);
            return values[index - 1] + (values[index] - values[index - 1]) * weight;
        }
    }
    values[values.len() - 1]
}

fn error(code: &'static str, message: impl Into<String>, field: &str) -> DiagnosticReport {
    DiagnosticReport::one(Diagnostic::new(
        code,
        message,
        field,
        Path::new("<case>"),
        1,
        1,
    ))
}
