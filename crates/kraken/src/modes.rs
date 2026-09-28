// Adapted from Acoustics Toolbox v2023.5 Kraken/kraken.f90 and
// Kraken/InverseIterationMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE. Distributed without warranty.
//!
//! Real, single-fluid finite-difference path from KRAKEN v2023.5.
//! Sturm counts isolate modes; inverse iteration samples the first mesh;
//! Richardson extrapolation refines eigenvalues only, as in the reference.
use crate::solver::error;
use crate::{
    Case, DiagnosticReport, Interpolation, MAX_MESH_POINTS, MAX_MODE_LIMIT, ModeSet, NormalMode,
};
use num_complex::Complex64;
use std::f64::consts::PI;

const MAX_SHAPES: usize = 5_000_000;
// ponytail: fixed numerical work budget; parallel root searches only if real cases require it.
const MAX_WORK: usize = 300_000_000;
const ROOT_STEPS: usize = 64;

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]
pub(super) fn solve(case: &Case) -> Result<ModeSet, DiagnosticReport> {
    let omega = 2.0 * PI * case.frequency_hz;
    let bottom_k2 = (omega / case.bottom_sound_speed_mps).powi(2);
    // AttenMod::CRCI converts dB/wavelength to a positive imaginary sound speed.
    let bottom_c = Complex64::new(
        case.bottom_sound_speed_mps,
        case.bottom_attenuation_db_per_wavelength * case.bottom_sound_speed_mps
            / (8.685_889_6 * 2.0 * PI),
    );
    let bottom_complex_k2 = (Complex64::new(omega, 0.0) / bottom_c).powi(2);
    let last_speed = case.sound_speed_profile.last().unwrap().sound_speed_mps;
    let needed = (case.water_depth_m / (last_speed / case.frequency_hz / 20.0))
        .floor()
        .max(10.0);
    if !needed.is_finite()
        || !omega.is_finite()
        || !bottom_k2.is_finite()
        || !bottom_complex_k2.re.is_finite()
        || !bottom_complex_k2.im.is_finite()
        || needed > 2.0 * MAX_MESH_POINTS as f64
    {
        return Err(error(
            "KR0302",
            "frequency/profile requires an unsupported mesh",
            "mesh_points",
        ));
    }
    let base = if case.mesh_points == 0 {
        needed as usize
    } else {
        case.mesh_points
    };
    if base < needed as usize / 2 || base > MAX_MESH_POINTS {
        return Err(error(
            "KR0302",
            "mesh is too coarse or exceeds the mesh limit",
            "mesh_points",
        ));
    }
    let mut table: Vec<Vec<f64>> = Vec::new();
    let mut modes = Vec::new();
    let mut work = 0_usize;
    for set in 0..5 {
        let multiplier = 1 << set;
        let n = base
            .checked_mul(multiplier)
            .filter(|n| *n <= MAX_MESH_POINTS)
            .ok_or_else(|| {
                error(
                    "KR0302",
                    "refined mesh exceeds the mesh limit",
                    "mesh_points",
                )
            })?;
        let mesh = Mesh::new(case, n, omega, bottom_k2, bottom_complex_k2)?;
        let roots = mesh.roots(&mut work)?;
        if set == 0 {
            modes = roots
                .iter()
                .map(|&x| mesh.mode(x))
                .collect::<Result<Vec<_>, _>>()?;
        } else if roots.len() != modes.len() {
            return Err(error(
                "KR0303",
                "mode count changed during mesh refinement; move spectral limits away from roots",
                "phase_speed_limits",
            ));
        }
        let key = 2 * modes.len() / 3;
        let previous = table.first().map(|row| row[key]);
        table.push(roots);
        for j in (0..set).rev() {
            let denominator = (multiplier as f64 / f64::from(1 << j)).powi(2) - 1.0;
            let (earlier, later) = table.split_at_mut(j + 1);
            for (value, &next) in earlier[j].iter_mut().zip(&later[0]) {
                *value = next - (*value - next) / denominator;
            }
        }
        let delta = previous.map_or(1e10, |x| (table[0][key] - x).abs());
        if delta * case.max_range_m < 1.0 {
            for (mode, &x) in modes.iter_mut().zip(&table[0]) {
                // KRAKEN combines extrapolated real k² with the first-mesh loss perturbation.
                let loss_k2 = mode.horizontal_wavenumber_rad_per_m.powi(2).im;
                let k = Complex64::new(x, loss_k2).sqrt();
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
                return Err(error("KR0302", "non-finite mode result", "modes"));
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
        "eigenvalue extrapolation did not converge within five meshes",
        "max_range_m",
    ))
}

struct Mesh<'a> {
    case: &'a Case,
    h: f64,
    omega: f64,
    bottom_k2: f64,
    bottom_complex_k2: Complex64,
    b1: Vec<f64>,
    min_speed: f64,
}

impl<'a> Mesh<'a> {
    #[allow(clippy::cast_precision_loss)]
    fn new(
        case: &'a Case,
        n: usize,
        omega: f64,
        bottom_k2: f64,
        bottom_complex_k2: Complex64,
    ) -> Result<Self, DiagnosticReport> {
        let h = case.water_depth_m / n as f64;
        let mut min_speed = f64::INFINITY;
        let b1: Vec<_> = (0..=n)
            .map(|i| {
                let speed = sound_speed(case, (i as f64 * h).min(case.water_depth_m));
                min_speed = min_speed.min(speed);
                -2.0 + h * h * (omega * omega / (speed * speed))
            })
            .collect();
        if h <= 0.0 || !h.is_finite() || b1.iter().any(|x| !x.is_finite()) {
            return Err(error(
                "KR0302",
                "mesh coefficients exceed numeric range",
                "mesh_points",
            ));
        }
        Ok(Self {
            case,
            h,
            omega,
            bottom_k2,
            bottom_complex_k2,
            b1,
            min_speed,
        })
    }

    fn bottom_gamma(&self, x: f64) -> Complex64 {
        (Complex64::new(x, 0.0) - self.bottom_complex_k2).sqrt()
    }

    fn bottom_diagonal(&self, x: f64) -> f64 {
        (self.b1.last().unwrap() - self.h * self.h * x) * 0.5
            - self.h * self.case.water_density_g_cm3 / self.case.bottom_density_g_cm3
                * self.bottom_gamma(x).re
    }

    // Inertia of the symmetric tridiagonal acoustic operator A(x): number of roots above x.
    fn count(&self, x: f64) -> usize {
        let shift = self.h * self.h * x;
        let mut pivot = 0.0;
        let mut count = 0;
        for i in 1..self.b1.len() {
            let d = if i + 1 == self.b1.len() {
                self.bottom_diagonal(x)
            } else {
                self.b1[i] - shift
            };
            pivot = if i == 1 { d } else { d - 1.0 / pivot };
            if pivot.abs() < 1e-30 {
                pivot = -1e-30;
            }
            if pivot > 0.0 {
                count += 1;
            }
        }
        count
    }

    #[allow(clippy::float_cmp)]
    fn roots(&self, work: &mut usize) -> Result<Vec<f64>, DiagnosticReport> {
        // Preserve Solve1's lower search guard, including its cutoff exclusion.
        let low = 1.00001 * (self.omega / self.case.c_high_mps).powi(2);
        let high = (self.omega / self.case.c_low_mps.max(self.min_speed)).powi(2);
        if low >= high || !high.is_finite() || low <= self.bottom_k2 {
            return Err(error(
                "KR0301",
                "phase-speed limits contain no trapped modes",
                "phase_speed_limits",
            ));
        }
        let above = self.count(high);
        let count = self.count(low).saturating_sub(above);
        if count == 0 {
            return Err(error(
                "KR0301",
                "no modes inside the spectral limits",
                "phase_speed_limits",
            ));
        }
        if count > MAX_MODE_LIMIT {
            return Err(error(
                "KR0302",
                "mode count exceeds the limit",
                "mode_count",
            ));
        }
        if count
            .checked_mul(self.case.mode_sample_depths_m.len())
            .is_none_or(|n| n > MAX_SHAPES)
        {
            return Err(error(
                "KR0302",
                "mode shape sample limit exceeded",
                "mode_sample_depths_m",
            ));
        }
        *work = count
            .checked_mul(self.b1.len())
            .and_then(|n| n.checked_mul(ROOT_STEPS + 5))
            .and_then(|n| n.checked_add(*work))
            .filter(|n| *n <= MAX_WORK)
            .ok_or_else(|| error("KR0302", "modal mesh-work limit exceeded", "mesh_points"))?;
        let mut roots = Vec::with_capacity(count);
        for mode in 0..count {
            let mut left = low;
            let mut right = high;
            for _ in 0..ROOT_STEPS {
                let middle = left.midpoint(right);
                if middle == left || middle == right {
                    break;
                }
                if self.count(middle) > above + mode {
                    left = middle;
                } else {
                    right = middle;
                }
            }
            roots.push(left.midpoint(right));
        }
        Ok(roots)
    }

    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::many_single_char_names
    )]
    fn mode(&self, x: f64) -> Result<NormalMode, DiagnosticReport> {
        let n = self.b1.len();
        let h_rho = self.h * self.case.water_density_g_cm3;
        let shift = self.h * self.h * x;
        let mut d: Vec<_> = self.b1.iter().map(|b| (b - shift) / h_rho).collect();
        let mut e = vec![1.0 / h_rho; n + 1];
        d[0] = 1.0;
        e[1] = 0.0; // vacuum surface
        d[n - 1] = self.bottom_diagonal(x) / h_rho;
        let mut phi = inverse_iteration(&d, &e)?;
        let mut norm = 0.0;
        let mut slow = 0.0;
        for (i, &value) in phi.iter().enumerate() {
            let weight = if i == 0 || i + 1 == n { 0.5 } else { 1.0 };
            let mass = weight * self.h * value * value / self.case.water_density_g_cm3;
            norm += mass;
            slow += mass * (self.b1[i] + 2.0) / (self.omega * self.omega * self.h * self.h);
        }
        let gamma = (x - self.bottom_complex_k2.re).sqrt();
        let x1 = 0.999_999_9 * x;
        let x2 = 1.000_000_1 * x;
        let derivative = (self.bottom_gamma(x2).re - self.bottom_gamma(x1).re)
            / (self.case.bottom_density_g_cm3 * (x2 - x1));
        norm += derivative * phi[n - 1].powi(2);
        slow += phi[n - 1].powi(2)
            / (2.0
                * gamma
                * self.case.bottom_density_g_cm3
                * self.case.bottom_sound_speed_mps.powi(2));
        if norm <= 0.0 || !norm.is_finite() || !slow.is_finite() {
            return Err(error("KR0303", "invalid mode normalization", "modes"));
        }
        let turning = (1..n)
            .find(|&i| self.b1[i] - shift + 2.0 > 0.0)
            .unwrap_or(n - 2);
        let scale = if phi[turning] < 0.0 {
            -norm.sqrt().recip()
        } else {
            norm.sqrt().recip()
        };
        for value in &mut phi {
            *value *= scale;
        }
        // Vector.f90 stores a single-precision mesh and interpolates into complex32 .mod samples.
        let grid_depths: Vec<_> = (0..n).map(|i| (i as f64 * self.h) as f32).collect();
        let eigenfunction = self
            .case
            .mode_sample_depths_m
            .iter()
            .map(|&depth| {
                let depth = depth as f32;
                let upper = grid_depths
                    .partition_point(|&sample| sample < depth)
                    .clamp(1, n - 1);
                let weight = (depth - grid_depths[upper - 1])
                    / (grid_depths[upper] - grid_depths[upper - 1]);
                let value = phi[upper - 1] as f32 + weight * (phi[upper] - phi[upper - 1]) as f32;
                Complex64::new(f64::from(value), 0.0)
            })
            .collect();
        // BCImpedance returns the real admittance for mode finding and the
        // complex admittance for first-order attenuation (Normalize in kraken.f90).
        let loss_k2 =
            -self.bottom_gamma(x).im * phi[n - 1].powi(2) / self.case.bottom_density_g_cm3;
        let k = Complex64::new(x, loss_k2).sqrt();
        Ok(NormalMode {
            horizontal_wavenumber_rad_per_m: k,
            phase_speed_mps: self.omega / k.re,
            group_speed_mps: x.sqrt() * norm / (self.omega * slow),
            attenuation_nepers_per_m: -k.im,
            eigenfunction,
        })
    }
}

fn sound_speed(case: &Case, depth: f64) -> f64 {
    let points = &case.sound_speed_profile;
    let upper = points
        .partition_point(|p| p.depth_m < depth)
        .clamp(1, points.len() - 1);
    let a = points[upper - 1];
    let b = points[upper];
    let weight = (depth - a.depth_m) / (b.depth_m - a.depth_m);
    match case.interpolation {
        Interpolation::N2Linear => ((1.0 - weight) / a.sound_speed_mps.powi(2)
            + weight / b.sound_speed_mps.powi(2))
        .sqrt()
        .recip(),
        Interpolation::CLinear => (1.0 - weight) * a.sound_speed_mps + weight * b.sound_speed_mps,
    }
}

// Specialized real tridiagonal inverse iteration, with row interchanges, from InverseIterationMod.f90.
#[allow(clippy::cast_precision_loss, clippy::many_single_char_names)]
fn inverse_iteration(d: &[f64], e: &[f64]) -> Result<Vec<f64>, DiagnosticReport> {
    let n = d.len();
    let eps3 = 100.0
        * f64::EPSILON
        * (d.iter().map(|v| v.abs()).sum::<f64>() + e[1..n].iter().map(|v| v.abs()).sum::<f64>());
    let eps4 = n as f64 * eps3;
    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut multipliers = vec![0.0; n];
    let mut swapped = vec![false; n];
    let mut u = d[0];
    let mut v = e[1];
    for i in 1..n {
        if e[i].abs() >= u.abs() {
            let ratio = u / e[i];
            multipliers[i] = ratio;
            swapped[i] = true;
            a[i - 1] = e[i];
            b[i - 1] = d[i];
            c[i - 1] = e[i + 1];
            u = v - ratio * d[i];
            v = -ratio * e[i + 1];
        } else {
            let ratio = e[i] / u;
            multipliers[i] = ratio;
            a[i - 1] = u;
            b[i - 1] = v;
            u = d[i] - ratio * v;
            v = e[i + 1];
        }
    }
    a[n - 1] = if u == 0.0 { eps3 } else { u };
    c[n - 2] = 0.0;
    let mut phi = vec![eps4 / (n as f64).sqrt(); n];
    for _ in 0..3 {
        let mut next = 0.0;
        let mut next2 = 0.0;
        for i in (0..n).rev() {
            phi[i] = (phi[i] - b[i] * next - c[i] * next2) / a[i];
            next2 = next;
            next = phi[i];
        }
        let norm: f64 = phi.iter().map(|x| x.abs()).sum();
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
            phi[i] -= multipliers[i] * phi[i - 1];
        }
    }
    Err(error(
        "KR0303",
        "inverse iteration failed to converge",
        "modes",
    ))
}

#[cfg(test)]
mod tests {
    use super::{solve, sound_speed};
    use crate::{Case, Interpolation, legacy::load_case, pekeris};
    use std::path::Path;

    #[test]
    fn interpolation_and_mesh_work_limit() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/Pekeris");
        let mut input = load_case(root.with_extension("env"), root.with_extension("flp"))
            .unwrap()
            .into_definition();
        input.sound_speed_profile[1].sound_speed_mps = 1600.0;
        let n2 = Case::from_definition(input.clone()).unwrap();
        let expected = (0.5 / 1500.0_f64.powi(2) + 0.5 / 1600.0_f64.powi(2))
            .sqrt()
            .recip();
        assert!((sound_speed(&n2, 50.0) - expected).abs() < 1e-12);
        input.interpolation = Interpolation::CLinear;
        let linear = Case::from_definition(input.clone()).unwrap();
        assert!((sound_speed(&linear, 50.0) - 1550.0).abs() < 1e-12);
        input.frequency_hz = 1000.0;
        input.mesh_points = 100_000;
        let report = solve(&Case::from_definition(input).unwrap()).unwrap_err();
        assert!(report.diagnostics()[0].message.contains("work limit"));
    }

    #[test]
    fn finite_difference_modes_match_independent_pekeris_oracle() {
        for name in ["Pekeris", "PekerisFiltered", "PekerisDense"] {
            let root = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures")
                .join(name);
            let case = load_case(root.with_extension("env"), root.with_extension("flp")).unwrap();
            let expected = pekeris::solve(&case);
            let actual = solve(&case).unwrap();
            assert_eq!(actual.modes.len(), expected.modes.len());
            for (actual, expected) in actual.modes.iter().zip(&expected.modes) {
                assert!(
                    (actual.horizontal_wavenumber_rad_per_m
                        - expected.horizontal_wavenumber_rad_per_m)
                        .norm()
                        < 5e-10
                );
                let sign = if actual
                    .eigenfunction
                    .iter()
                    .zip(&expected.eigenfunction)
                    .map(|(a, b)| a.re * b.re)
                    .sum::<f64>()
                    < 0.0
                {
                    -1.0
                } else {
                    1.0
                };
                for (a, b) in actual.eigenfunction.iter().zip(&expected.eigenfunction) {
                    assert!((*a - sign * *b).norm() < 1e-6);
                }
            }
        }
    }
}
