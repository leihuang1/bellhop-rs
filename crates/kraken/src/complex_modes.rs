// Adapted from Acoustics Toolbox v2023.5 Kraken/krakenc.f90,
// Kraken/InverseIterationMod.f90 and misc/RootFinderSecantMod.f90.
// Copyright (C) 2009 Michael B. Porter. GPL-3.0-or-later; see LICENSE.
//! Lossless N-profile fluid water column over a lossless fluid half-space,
//! vacuum surface, base mesh only, modes only.
use crate::profile::Profile;
use crate::solver::error;
use crate::{
    BottomBoundary, Case, DiagnosticReport, Interpolation, ModeSet, NormalMode, SurfaceBoundary,
};
use num_complex::Complex64;
use std::f64::consts::PI;

const MAX_ROOT_WORK: usize = 30_000_000;
const SECANT_RELATIVE_TOLERANCE: f64 = 1e-14;

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]
pub(super) fn solve(case: &Case) -> Result<ModeSet, DiagnosticReport> {
    if case.interpolation != Interpolation::N2Linear
        || case.surface_boundary != SurfaceBoundary::Vacuum
        || case.bottom_boundary != BottomBoundary::FluidHalfSpace
        || case.bottom_attenuation_db_per_wavelength != 0.0
        || case.max_range_m != 0.0
        || case.c_high_mps <= case.bottom_sound_speed_mps
        || (case.sound_speed_profile.len() > 7
            && case
                .sound_speed_profile
                .iter()
                .any(|p| p.sound_speed_mps != case.sound_speed_profile[0].sound_speed_mps))
    {
        return Err(error(
            "KR0302",
            "KRAKENC currently requires a constant or at most seven-point lossless N-profile with a vacuum surface and lossless fluid bottom, leaky spectral interval, and no mesh extrapolation",
            "mode_solver",
        ));
    }
    let omega = 2.0 * PI * case.frequency_hz;
    let needed = (case.water_depth_m * case.frequency_hz * 20.0
        / case.sound_speed_profile.last().unwrap().sound_speed_mps)
        .floor()
        .max(10.0);
    let n = if case.mesh_points == 0 {
        needed as usize
    } else {
        case.mesh_points
    };
    if !needed.is_finite() || n < needed as usize / 2 || n > crate::MAX_MESH_POINTS {
        return Err(error("KR0302", "unsupported KRAKENC mesh", "mesh_points"));
    }
    let h = case.water_depth_m / n as f64;
    let profile = Profile::new(case)?;
    let mut min_speed = f64::INFINITY;
    let b: Vec<_> = (0..=n)
        .map(|i| {
            let c = profile.speed((i as f64 * h).min(case.water_depth_m));
            min_speed = min_speed.min(c);
            -2.0 + h * h * (omega / c).powi(2)
        })
        .collect();
    let bottom_k2 = (omega / case.bottom_sound_speed_mps).powi(2);
    let low_k2 = (omega / case.c_high_mps).powi(2);
    let high_k2 = (omega / case.c_low_mps).powi(2);
    let water_k2 = (omega / min_speed).powi(2);
    if !low_k2.is_finite() || !water_k2.is_finite() || b.iter().any(|x| !x.is_finite()) {
        return Err(error(
            "KR0302",
            "mesh coefficients exceed numeric range",
            "mesh_points",
        ));
    }
    // The minimum water speed bounds the first root; deflation keeps
    // subsequent secants distinct.
    let guesses =
        (case.water_depth_m * (water_k2 - low_k2).max(0.0).sqrt() / PI).ceil() as usize + 1;
    // ponytail: worst-case 1000 secant steps per guess; a counted budget can
    // replace this ceiling when broader KRAKENC cases need it.
    if guesses > crate::MAX_MODE_LIMIT
        || guesses == 0
        || guesses
            .checked_mul(n)
            .and_then(|w| w.checked_mul(1000))
            .is_none_or(|w| w > MAX_ROOT_WORK)
    {
        return Err(error(
            "KR0302",
            "complex root work limit exceeded",
            "mesh_points",
        ));
    }
    let mut roots = Vec::new();
    let mut reached_lower_limit = false;
    for index in 1..=guesses {
        let vertical = (index as f64 - 0.5) * PI / case.water_depth_m;
        let guess = Complex64::new(water_k2 - vertical * vertical, 0.0);
        let root = secant(
            guess,
            &roots,
            &b,
            h,
            case.water_density_g_cm3,
            case.bottom_density_g_cm3,
            bottom_k2,
        )?;
        if root.re <= low_k2 {
            reached_lower_limit = true;
            break;
        }
        // Even excluded roots must be deflated so the next secant finds a new mode.
        if roots.iter().any(|&previous: &Complex64| {
            (root - previous).norm()
                < root.norm().max(previous.norm()) * b.len() as f64 * SECANT_RELATIVE_TOLERANCE
        }) {
            return Err(error(
                "KR0303",
                "complex root search repeated a mode",
                "phase_speed_limits",
            ));
        }
        roots.push(root);
    }
    if !reached_lower_limit {
        return Err(error(
            "KR0303",
            "complex root search cannot prove the full spectral interval",
            "phase_speed_limits",
        ));
    }
    let selected: Vec<_> = roots.into_iter().filter(|x| x.re <= high_k2).collect();
    if selected.is_empty() {
        return Err(error(
            "KR0301",
            "no complex modes inside spectral limits",
            "phase_speed_limits",
        ));
    }
    if selected
        .len()
        .checked_mul(case.mode_sample_depths_m.len())
        .is_none_or(|w| w > 5_000_000)
    {
        return Err(error(
            "KR0302",
            "mode shape sample limit exceeded",
            "mode_sample_depths_m",
        ));
    }
    let modes = selected
        .iter()
        .map(|&root| mode(case, &b, h, omega, bottom_k2, root))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ModeSet {
        frequency_hz: case.frequency_hz,
        sampled_depths_m: case.mode_sample_depths_m.clone(),
        modes,
    })
}

// The upstream PekerisRoot branch is not the principal complex square root.
fn pekeris_root(z: Complex64) -> Complex64 {
    if z.re >= 0.0 {
        z.sqrt()
    } else {
        Complex64::new(0.0, 1.0) * (-z).sqrt()
    }
}

fn dispersion(
    x: Complex64,
    roots: &[Complex64],
    b: &[f64],
    h: f64,
    water_rho: f64,
    bottom_rho: f64,
    bottom_k2: f64,
) -> (Complex64, i32) {
    let gamma = pekeris_root(x - bottom_k2);
    let mut prev = Complex64::new(-2.0 * bottom_rho, 0.0);
    let mut current = (b[b.len() - 1] - h * h * x) * bottom_rho - 2.0 * h * gamma * water_rho;
    // Keep the decimal scaling exponent: secant compares evaluations at
    // different scales after shooting and deflation (as in RootFinderSecant).
    let mut power = 0;
    for &coefficient in b[..b.len() - 1].iter().rev() {
        let next = (h * h * x - coefficient) * current - prev;
        prev = current;
        current = next;
        if current.re.abs() > 1e50 {
            prev *= 1e-50;
            current *= 1e-50;
            power += 50;
        }
    }
    // Vacuum top: Delta = -g = p1 (previous), not the next recurrence value.
    let mut value = prev;
    for &root in roots {
        value /= x - root;
        if value.re.abs() > 1e50 {
            value *= 1e-50;
            power += 50;
        }
        if value.re.abs() < 1e-50 && value.norm() > 0.0 {
            value *= 1e50;
            power -= 50;
        }
    }
    (value, power)
}

#[allow(clippy::cast_precision_loss)]
fn secant(
    mut x: Complex64,
    roots: &[Complex64],
    b: &[f64],
    h: f64,
    water_rho: f64,
    bottom_rho: f64,
    bottom_k2: f64,
) -> Result<Complex64, DiagnosticReport> {
    let tolerance = x.norm() * b.len() as f64 * SECANT_RELATIVE_TOLERANCE;
    let mut previous = x + 100.0 * tolerance;
    let (mut f_previous, mut previous_power) =
        dispersion(previous, roots, b, h, water_rho, bottom_rho, bottom_k2);
    for _ in 0..1000 {
        let (f, power) = dispersion(x, roots, b, h, water_rho, bottom_rho, bottom_k2);
        let numerator = f * (x - previous);
        let denominator = f - f_previous * 10_f64.powi(previous_power - power);
        let shift = if numerator.norm() >= (denominator * x).norm() {
            Complex64::new(0.1 * tolerance, 0.0)
        } else {
            numerator / denominator
        };
        let next = x - shift;
        if !next.re.is_finite() || !next.im.is_finite() {
            break;
        }
        if (next - x).norm() + (next - previous).norm() < tolerance {
            return Ok(next);
        }
        previous = x;
        f_previous = f;
        previous_power = power;
        x = next;
    }
    Err(error(
        "KR0303",
        "complex secant did not converge",
        "phase_speed_limits",
    ))
}

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::many_single_char_names
)]
fn mode(
    case: &Case,
    b: &[f64],
    h: f64,
    omega: f64,
    bottom_k2: f64,
    x: Complex64,
) -> Result<NormalMode, DiagnosticReport> {
    let n = b.len();
    let shift = h * h * x;
    // ponytail: shooting suffices for the 1 km water column; deeper
    // evanescent columns need scaled shooting or inverse iteration.
    let mut phi = vec![Complex64::new(0.0, 0.0); n];
    phi[1] = Complex64::new(1.0, 0.0);
    for i in 2..n {
        phi[i] = (shift - b[i - 1]) * phi[i - 1] - phi[i - 2];
        if !phi[i].re.is_finite() || !phi[i].im.is_finite() {
            return Err(error(
                "KR0303",
                "complex eigenfunction shooting overflow",
                "modes",
            ));
        }
    }
    let mut sq_norm = Complex64::new(0.0, 0.0);
    let mut slow = Complex64::new(0.0, 0.0);
    for (i, &value) in phi.iter().enumerate() {
        let weight = if i == 0 || i + 1 == n { 0.5 } else { 1.0 };
        let mass = weight * h / case.water_density_g_cm3 * value * value;
        sq_norm += mass;
        slow += mass * (b[i] + 2.0) / (omega * omega * h * h);
    }
    let gamma = (x - bottom_k2).sqrt();
    slow += phi[n - 1].powi(2)
        / (2.0 * gamma * case.bottom_density_g_cm3 * case.bottom_sound_speed_mps.powi(2));
    let x1 = x * 0.999_999_9;
    let x2 = x * 1.000_000_1;
    let derivative = (pekeris_root(x2 - bottom_k2) - pekeris_root(x1 - bottom_k2))
        / (case.bottom_density_g_cm3 * (x2 - x1));
    let norm = sq_norm + derivative * phi[n - 1].powi(2);
    let turning = (1..n)
        .find(|&i| (b[i] - shift).re + 2.0 > 0.0)
        .unwrap_or(n - 2);
    let mut scale = Complex64::new(1.0, 0.0) / norm.sqrt();
    if (scale * phi[turning]).re < 0.0 {
        scale = -scale;
    }
    for value in &mut phi {
        *value *= scale;
    }
    let group_speed = (Complex64::new(1.0, 0.0) / (scale * scale * slow * omega / x.sqrt())).re;
    let grid: Vec<_> = (0..n).map(|i| (i as f64 * h) as f32).collect();
    let eigenfunction: Vec<Complex64> = case
        .mode_sample_depths_m
        .iter()
        .map(|&depth| {
            let depth = depth as f32;
            let upper = grid.partition_point(|&z| z < depth).clamp(1, n - 1);
            let weight = (depth - grid[upper - 1]) / (grid[upper] - grid[upper - 1]);
            let first = phi[upper - 1];
            let difference = phi[upper] - first;
            Complex64::new(
                f64::from(first.re as f32 + weight * difference.re as f32),
                f64::from(first.im as f32 + weight * difference.im as f32),
            )
        })
        .collect();
    let k = x.sqrt();
    if !k.re.is_finite()
        || !k.im.is_finite()
        || k.im > 0.0
        || !group_speed.is_finite()
        || !norm.re.is_finite()
        || !norm.im.is_finite()
        || eigenfunction
            .iter()
            .any(|v: &Complex64| !v.re.is_finite() || !v.im.is_finite())
    {
        return Err(error(
            "KR0303",
            "invalid complex mode normalization",
            "modes",
        ));
    }
    Ok(NormalMode {
        horizontal_wavenumber_rad_per_m: k,
        phase_speed_mps: omega / k.re,
        group_speed_mps: group_speed,
        attenuation_nepers_per_m: -k.im,
        eigenfunction,
    })
}
