// Adapted from Acoustics Toolbox v2023.5 Kraken/krakenc.f90,
// Kraken/InverseIterationMod.f90 and misc/RootFinderSecantMod.f90.
// Copyright (C) 2009 Michael B. Porter. GPL-3.0-or-later; see LICENSE.
//! Lossless N/C-profile fluid water column over a fluid half-space,
//! vacuum surface and optional Richardson mesh extrapolation.
use crate::profile::Profile;
use crate::solver::error;
use crate::{
    BottomBoundary, Case, DiagnosticReport, Interpolation, ModeSet, NormalMode, SurfaceBoundary,
};
use num_complex::Complex64;
use std::f64::consts::PI;

const MAX_ROOT_WORK: usize = 300_000_000;
const SECANT_RELATIVE_TOLERANCE: f64 = 1e-14;

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]
pub(super) fn solve(case: &Case) -> Result<ModeSet, DiagnosticReport> {
    if !matches!(
        case.interpolation,
        Interpolation::N2Linear | Interpolation::CLinear
    ) || case.surface_boundary != SurfaceBoundary::Vacuum
        || case.bottom_boundary != BottomBoundary::FluidHalfSpace
    {
        return Err(error(
            "KR0302",
            "KRAKENC currently requires lossless N/C-profile water with a vacuum surface and fluid bottom with optional W attenuation",
            "mode_solver",
        ));
    }
    let omega = 2.0 * PI * case.frequency_hz;
    let profile = Profile::new(case)?;
    let bottom_c = Complex64::new(
        case.bottom_sound_speed_mps,
        case.bottom_attenuation_db_per_wavelength * case.bottom_sound_speed_mps
            / (8.685_889_6 * 2.0 * PI),
    );
    let bottom_k2 = (Complex64::new(omega, 0.0) / bottom_c).powi(2);
    let low_k2 = (omega / case.c_high_mps).powi(2);
    let high_k2 = (omega / case.c_low_mps).powi(2);
    if !low_k2.is_finite() || !bottom_k2.re.is_finite() || !bottom_k2.im.is_finite() {
        return Err(error(
            "KR0302",
            "mesh coefficients exceed numeric range",
            "mesh_points",
        ));
    }

    let mut table: Vec<Vec<Complex64>> = Vec::new();
    let mut modes = Vec::new();
    let mut work = 0;
    for set in 0..5 {
        let multiplier = 1_usize << set;
        let n = case.mesh_points_at(multiplier)?;
        let h = case.water_depth_m / n as f64;
        let mut min_speed = f64::INFINITY;
        let b: Vec<_> = (0..=n)
            .map(|i| {
                let c = profile.speed((i as f64 * h).min(case.water_depth_m));
                min_speed = min_speed.min(c);
                -2.0 + h * h * (omega / c).powi(2)
            })
            .collect();
        let water_k2 = (omega / min_speed).powi(2);
        if !water_k2.is_finite() || b.iter().any(|x| !x.is_finite()) {
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
        if guesses > crate::MAX_MODE_LIMIT || guesses == 0 {
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
            // ponytail: wide spectral intervals use the previous-root seed from
            // Fortran; revisit the 10x cutoff if a medium-width spectrum fails.
            let guess = if case.c_high_mps > 10.0 * case.c_low_mps {
                roots.last().map_or_else(
                    || Complex64::new(water_k2 - vertical * vertical, 0.0),
                    |&previous| previous * 1.000_01,
                )
            } else {
                Complex64::new(water_k2 - vertical * vertical, 0.0)
            };
            let root = secant(
                guess,
                &roots,
                &b,
                h,
                case.water_density_g_cm3,
                case.bottom_density_g_cm3,
                bottom_k2,
                &mut work,
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
        // Upstream computes shapes/group speeds only on the first mesh.
        if set == 0 {
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
            modes = selected
                .iter()
                .map(|&root| mode(case, &b, h, omega, bottom_c, bottom_k2, root))
                .collect::<Result<Vec<_>, _>>()?;
        } else if selected.len() != modes.len() {
            return Err(error(
                "KR0303",
                "mode count changed during KRAKENC mesh refinement; move spectral limits away from roots",
                "phase_speed_limits",
            ));
        }

        let key = 2 * modes.len() / 3;
        let previous = table.first().map(|row| row[key]);
        table.push(selected);
        // Richardson-extrapolate complex k² with mesh ratios 1, 2, 4, 8, 16.
        for j in (0..set).rev() {
            let denominator = (multiplier as f64 / (1_usize << j) as f64).powi(2) - 1.0;
            let (earlier, later) = table.split_at_mut(j + 1);
            for (value, &next) in earlier[j].iter_mut().zip(&later[0]) {
                *value = next - (*value - next) / denominator;
            }
        }
        let delta = previous.map_or(1e10, |x| (table[0][key] - x).norm());
        if delta * case.max_range_m < 1.0 {
            for (mode, &x) in modes.iter_mut().zip(&table[0]) {
                let k = x.sqrt();
                mode.horizontal_wavenumber_rad_per_m = k;
                mode.phase_speed_mps = omega / k.re;
                mode.attenuation_nepers_per_m = -k.im;
            }
            if modes.iter().any(|mode| {
                !mode.phase_speed_mps.is_finite()
                    || !mode.group_speed_mps.is_finite()
                    || !mode.attenuation_nepers_per_m.is_finite()
                    || !mode.horizontal_wavenumber_rad_per_m.re.is_finite()
                    || !mode.horizontal_wavenumber_rad_per_m.im.is_finite()
                    || mode
                        .eigenfunction
                        .iter()
                        .any(|v| !v.re.is_finite() || !v.im.is_finite())
            }) {
                return Err(error("KR0303", "invalid complex mode result", "modes"));
            }
            return Ok(ModeSet {
                frequency_hz: case.frequency_hz,
                sampled_depths_m: case.mode_sample_depths_m.clone(),
                modes,
            });
        }
    }
    Err(error(
        "KR0303",
        "KRAKENC eigenvalue extrapolation did not converge within five meshes",
        "max_range_m",
    ))
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
    bottom_k2: Complex64,
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

#[allow(clippy::cast_precision_loss, clippy::too_many_arguments)]
fn secant(
    mut x: Complex64,
    roots: &[Complex64],
    b: &[f64],
    h: f64,
    water_rho: f64,
    bottom_rho: f64,
    bottom_k2: Complex64,
    work: &mut usize,
) -> Result<Complex64, DiagnosticReport> {
    let mut evaluate = |x| {
        *work += b.len() + roots.len();
        if *work > MAX_ROOT_WORK {
            return Err(error(
                "KR0302",
                "complex root work limit exceeded",
                "mesh_points",
            ));
        }
        Ok(dispersion(x, roots, b, h, water_rho, bottom_rho, bottom_k2))
    };
    let tolerance = x.norm() * b.len() as f64 * SECANT_RELATIVE_TOLERANCE;
    let mut previous = x + 100.0 * tolerance;
    let (mut f_previous, mut previous_power) = evaluate(previous)?;
    for _ in 0..1000 {
        let (f, power) = evaluate(x)?;
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

// Complex tridiagonal inverse iteration from upstream InverseIterationMod.f90.
#[allow(
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::many_single_char_names
)]
fn inverse_iteration(d: &[Complex64], e: &[f64]) -> Result<Vec<Complex64>, DiagnosticReport> {
    let n = d.len();
    let eps3 = 100.0
        * f64::EPSILON
        * (d.iter().map(|v| v.re.abs() + v.im.abs()).sum::<f64>()
            + e[1..n].iter().map(|v| v.abs()).sum::<f64>());
    let eps4 = n as f64 * eps3;
    let mut a = vec![Complex64::new(0.0, 0.0); n];
    let mut b = a.clone();
    let mut c = a.clone();
    let mut multipliers = a.clone();
    let mut swapped = vec![false; n];
    let mut u = d[0];
    let mut v = Complex64::new(e[1], 0.0);
    for i in 1..n {
        if e[i].abs() >= u.norm() {
            let ratio = u / e[i];
            multipliers[i] = ratio;
            swapped[i] = true;
            a[i - 1] = e[i].into();
            b[i - 1] = d[i];
            c[i - 1] = e[i + 1].into();
            u = v - ratio * d[i];
            v = -ratio * e[i + 1];
        } else {
            let ratio = e[i] / u;
            multipliers[i] = ratio;
            a[i - 1] = u;
            b[i - 1] = v;
            u = d[i] - ratio * v;
            v = e[i + 1].into();
        }
    }
    a[n - 1] = if u == Complex64::new(0.0, 0.0) {
        Complex64::new(eps3, 0.0)
    } else {
        u
    };
    c[n - 2] = Complex64::new(0.0, 0.0);
    let mut phi = vec![Complex64::new(eps4 / (n as f64).sqrt(), 0.0); n];
    for _ in 0..3 {
        let mut next = Complex64::new(0.0, 0.0);
        let mut next2 = next;
        for i in (0..n).rev() {
            phi[i] = (phi[i] - b[i] * next - c[i] * next2) / a[i];
            next2 = next;
            next = phi[i];
        }
        let norm: f64 = phi.iter().map(|v| v.re.abs() + v.im.abs()).sum();
        if !norm.is_finite() || norm == 0.0 {
            break;
        }
        if norm >= 1.0 {
            return Ok(phi);
        }
        for value in &mut phi {
            *value *= eps4 / norm;
        }
        for i in 1..n {
            if swapped[i] {
                phi.swap(i, i - 1);
            }
            let before = phi[i - 1];
            phi[i] -= multipliers[i] * before;
        }
    }
    Err(error("KR0303", "complex inverse iteration failed", "modes"))
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
    bottom_c: Complex64,
    bottom_k2: Complex64,
    x: Complex64,
) -> Result<NormalMode, DiagnosticReport> {
    let n = b.len();
    let shift = h * h * x;
    let h_rho = h * case.water_density_g_cm3;
    let mut d: Vec<_> = b.iter().map(|&v| (v - shift) / h_rho).collect();
    let mut e = vec![1.0 / h_rho; n + 1];
    d[0] = Complex64::new(1.0, 0.0);
    e[1] = 0.0;
    d[n - 1] = (b[n - 1] - shift) / (2.0 * h_rho)
        - pekeris_root(x - bottom_k2) / case.bottom_density_g_cm3;
    let mut phi = inverse_iteration(&d, &e)?;
    let mut sq_norm = Complex64::new(0.0, 0.0);
    let mut slow = Complex64::new(0.0, 0.0);
    for (i, &value) in phi.iter().enumerate() {
        let weight = if i == 0 || i + 1 == n { 0.5 } else { 1.0 };
        let mass = weight * h / case.water_density_g_cm3 * value * value;
        sq_norm += mass;
        slow += mass * (b[i] + 2.0) / (omega * omega * h * h);
    }
    let gamma = (x - bottom_k2).sqrt();
    slow += phi[n - 1].powi(2) / (2.0 * gamma * case.bottom_density_g_cm3 * bottom_c.powi(2));
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
        // Pinned Fortran roots can carry positive imaginary roundoff near 1e-18.
        || k.im > k.re * 1e-12
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
