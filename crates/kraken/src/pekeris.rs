//! Independent analytical oracle for homogeneous-fluid tests, not the production solver.
use crate::{Case, ModeSet, NormalMode};
use num_complex::Complex64;
use std::f64::consts::PI;

#[allow(clippy::float_cmp)]
pub(super) fn solve(case: &Case) -> ModeSet {
    let c = case.sound_speed_profile[0].sound_speed_mps;
    assert!(
        case.sound_speed_profile
            .iter()
            .all(|point| point.sound_speed_mps == c)
    );
    let omega = 2.0 * PI * case.frequency_hz;
    let water_k = omega / c;
    let bottom_k = omega / case.bottom_sound_speed_mps;
    let cutoff = (water_k * water_k - bottom_k * bottom_k).sqrt();
    let ratio = case.water_density_g_cm3 / case.bottom_density_g_cm3;
    let dispersion = |q: f64| {
        (q * case.water_depth_m).cos()
            + ratio * (cutoff * cutoff - q * q).max(0.0).sqrt() / q * (q * case.water_depth_m).sin()
    };
    let mut modes = Vec::new();
    for n in 0..20_000 {
        let mut lower = (f64::from(n) + 0.5) * PI / case.water_depth_m;
        let mut upper = ((f64::from(n) + 1.0) * PI / case.water_depth_m).min(cutoff);
        if lower >= upper {
            break;
        }
        let mut left = dispersion(lower);
        if left.is_sign_positive() == dispersion(upper).is_sign_positive() {
            continue;
        }
        for _ in 0..100 {
            let middle = lower.midpoint(upper);
            if middle == lower || middle == upper {
                break;
            }
            let value = dispersion(middle);
            if value.is_sign_positive() == left.is_sign_positive() {
                lower = middle;
                left = value;
            } else {
                upper = middle;
            }
        }
        let q = lower.midpoint(upper);
        let k = (water_k * water_k - q * q).sqrt();
        let phase_speed = omega / k;
        if phase_speed < case.c_low_mps || phase_speed > case.c_high_mps {
            continue;
        }
        let gamma = (k * k - bottom_k * bottom_k).sqrt();
        let water_norm = (case.water_depth_m * 0.5
            - (2.0 * q * case.water_depth_m).sin() / (4.0 * q))
            / case.water_density_g_cm3;
        let bottom_norm =
            (q * case.water_depth_m).sin().powi(2) / (2.0 * gamma * case.bottom_density_g_cm3);
        let normalization = (water_norm + bottom_norm).sqrt().recip();
        let slow = water_norm / c.powi(2) + bottom_norm / case.bottom_sound_speed_mps.powi(2);
        modes.push(NormalMode {
            horizontal_wavenumber_rad_per_m: Complex64::new(k, 0.0),
            phase_speed_mps: phase_speed,
            group_speed_mps: k * (water_norm + bottom_norm) / (omega * slow),
            attenuation_nepers_per_m: 0.0,
            eigenfunction: case
                .mode_sample_depths_m
                .iter()
                .map(|z| Complex64::new(normalization * (q * z).sin(), 0.0))
                .collect(),
        });
    }
    ModeSet {
        frequency_hz: case.frequency_hz,
        sampled_depths_m: case.mode_sample_depths_m.clone(),
        modes,
    }
}
