// Adapted from Acoustics Toolbox v2023.5 misc/AttenMod.f90,
// Copyright (C) 2009 Michael B. Porter. GPL-3.0-or-later; see LICENSE.
//! Legacy loss conversion. Public cases use dB/wavelength at their solve frequency.
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub(crate) enum VolumeLoss {
    None,
    Thorp,
    FrancoisGarrison([f64; 4]), // temperature, salinity, pH, mean depth
    Biological(Vec<[f64; 5]>),  // top, bottom, resonance frequency, Q, a0
}

#[allow(clippy::excessive_precision)] // preserve the reference's f32 divisor literal
pub(crate) fn db_per_wavelength(
    unit: u8,
    volume: &VolumeLoss,
    value: f64,
    depth: f64,
    speed: f64,
    frequency: f64,
    power_law: (f64, f64, f64), // reference frequency, beta, transition frequency
) -> f64 {
    // Keep the existing W path bit-for-bit, including its zero-loss cases.
    if unit == b'W' && matches!(volume, VolumeLoss::None) {
        return value;
    }
    let omega = 2.0 * PI * frequency;
    let mut nepers = match unit {
        b'N' => value,
        b'M' => value / 8.685_889_6,
        b'm' => {
            let (reference, beta, transition) = power_law;
            let scaling = if frequency < transition {
                (frequency / reference).powf(beta)
            } else {
                (frequency / reference) * (transition / reference).powf(beta - 1.0)
            };
            value / 8.685_889_6 * scaling
        }
        b'F' => value * frequency / 8_685.889_6,
        b'Q' if value != 0.0 => omega / (2.0 * speed * value),
        b'Q' => 0.0,
        b'L' => value * omega / speed,
        b'W' => value * frequency / (8.685_889_6 * speed),
        _ => unreachable!("validated legacy attenuation unit"),
    };
    nepers += match volume {
        VolumeLoss::None => 0.0,
        VolumeLoss::Thorp => {
            let f2 = (frequency / 1000.0).powi(2);
            // The divisor literal in the reference Thorp path is single precision.
            (3.3e-3 + 0.11 * f2 / (1.0 + f2) + 44.0 * f2 / (4100.0 + f2) + 3e-4 * f2)
                / f64::from(8685.8896_f32)
        }
        VolumeLoss::FrancoisGarrison(parameters) => {
            francois_garrison(frequency / 1000.0, *parameters) / f64::from(8685.8896_f32)
        }
        VolumeLoss::Biological(layers) => layers
            .iter()
            .filter(|p| depth >= p[0] && depth <= p[1])
            .map(|p| {
                p[4] / ((1.0 - p[2].powi(2) / frequency.powi(2)).powi(2) + 1.0 / p[3].powi(2))
                    / f64::from(8685.8896_f32)
            })
            .sum(),
    };
    if nepers == 0.0 {
        0.0
    } else {
        nepers * (8.685_889_6 * speed) / frequency
    }
}

fn francois_garrison(f: f64, [t, salinity, ph, depth]: [f64; 4]) -> f64 {
    // Unsuffixed Fortran real literals are rounded to f32 before promotion.
    let speed = 1412.0
        + f64::from(3.21_f32) * t
        + f64::from(1.19_f32) * salinity
        + f64::from(0.0167_f32) * depth;
    let a1 = f64::from(8.86_f32) / speed * 10_f64.powf(f64::from(0.78_f32) * ph - 5.0);
    let f1 =
        f64::from(2.8_f32) * (salinity / 35.0).sqrt() * 10_f64.powf(4.0 - 1245.0 / (t + 273.0));
    let a2 = f64::from(21.44_f32) * salinity / speed * (1.0 + f64::from(0.025_f32) * t);
    let p2 = 1.0 - 1.37e-4 * depth + 6.2e-9 * depth.powi(2);
    let f2 = f64::from(8.17_f32) * 10_f64.powf(8.0 - 1990.0 / (t + 273.0))
        / (1.0 + f64::from(0.0018_f32) * (salinity - 35.0));
    let p3 = 1.0 - 3.83e-5 * depth + 4.9e-10 * depth.powi(2);
    let a3 = if t < 20.0 {
        4.937e-4 - 2.59e-5 * t + 9.11e-7 * t.powi(2) - 1.5e-8 * t.powi(3)
    } else {
        3.964e-4 - 1.146e-5 * t + 1.45e-7 * t.powi(2) - 6.5e-10 * t.powi(3)
    };
    a1 * (f1 * f * f) / (f1 * f1 + f * f)
        + a2 * p2 * (f2 * f * f) / (f2 * f2 + f * f)
        + a3 * p3 * f * f
}
