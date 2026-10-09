// Adapted from Acoustics Toolbox v2023.5 misc/pchipMod.f90, splinec.f90, and munk.f90,
// Copyright (C) 2009 Michael B. Porter. GPL-3.0-or-later; see LICENSE.
// Constant-density fluid profiles. Complex speeds are converted before interpolation.
use crate::solver::{elastic::complex_speed, error};
use crate::{CaseDefinition, DiagnosticReport, Interpolation};
use num_complex::Complex64;
use std::f64::consts::PI;

pub(super) struct Profile<'a> {
    case: &'a CaseDefinition,
    layer: crate::solver::layers::Layer<'a>,
    minimum_speed: f64,
    // Local polynomial coefficients, c0 + t * (c1 + t * (c2 + t * c3)).
    cubic: Vec<[f64; 4]>,
    imaginary_speeds: Vec<f64>,
    imaginary_cubic: Vec<[f64; 4]>,
}

impl<'a> Profile<'a> {
    pub(super) fn new(case: &'a CaseDefinition) -> Result<Self, DiagnosticReport> {
        Self::new_layer(case, crate::solver::layers::iter(case).next().unwrap())
    }

    pub(super) fn new_layer(
        case: &'a CaseDefinition,
        layer: crate::solver::layers::Layer<'a>,
    ) -> Result<Self, DiagnosticReport> {
        Self::build(case, layer, false)
    }

    pub(crate) fn new_elastic_layer(
        case: &'a CaseDefinition,
        layer: crate::solver::layers::Layer<'a>,
    ) -> Result<Self, DiagnosticReport> {
        Self::build(case, layer, true)
    }

    fn build(
        case: &'a CaseDefinition,
        layer: crate::solver::layers::Layer<'a>,
        elastic: bool,
    ) -> Result<Self, DiagnosticReport> {
        let points = layer.points;
        let mut minimum_speed = if case.interpolation == Interpolation::AnalyticMunk {
            1500.0
        } else {
            points
                .iter()
                .map(|point| point.sound_speed_mps)
                .fold(f64::INFINITY, f64::min)
        };
        let mut cubic = Vec::new();
        if matches!(
            case.interpolation,
            Interpolation::Pchip | Interpolation::Spline
        ) {
            let h: Vec<_> = points
                .windows(2)
                .map(|p| p[1].depth_m - p[0].depth_m)
                .collect();
            let speeds: Vec<_> = points.iter().map(|p| p.sound_speed_mps).collect();
            cubic = cubic_coefficients(case.interpolation, &h, &speeds)?;
            for (&coefficients, &step) in cubic.iter().zip(&h) {
                let segment_minimum = segment_minimum(coefficients, step);
                if !segment_minimum.is_finite() {
                    return Err(invalid_profile());
                }
                minimum_speed = minimum_speed.min(segment_minimum);
            }
            if minimum_speed <= 0.0 {
                return Err(error(
                    "KR0302",
                    "interpolated sound speed must stay positive between profile points",
                    "sound_speed_profile",
                ));
            }
        }
        let imaginary_speeds: Vec<_> = if layer.loss.iter().any(|&a| a != 0.0) {
            points
                .iter()
                .zip(layer.loss)
                .map(|(p, &a)| {
                    if elastic {
                        // CRCI converts wavelength loss to neper/m, then multiplies by c²/omega.
                        let neper = a * case.frequency_hz / (8.685_889_6 * p.sound_speed_mps);
                        neper * p.sound_speed_mps.powi(2) / (2.0 * PI * case.frequency_hz)
                    } else {
                        complex_speed(p.sound_speed_mps, a).im
                    }
                })
                .collect()
        } else {
            Vec::new()
        };
        if imaginary_speeds.iter().any(|c| !c.is_finite()) {
            return Err(invalid_profile());
        }
        let mut imaginary_cubic = Vec::new();
        if !imaginary_speeds.is_empty() && !cubic.is_empty() {
            let h: Vec<_> = points
                .windows(2)
                .map(|p| p[1].depth_m - p[0].depth_m)
                .collect();
            imaginary_cubic = cubic_coefficients(case.interpolation, &h, &imaginary_speeds)?;
            for ((&real, &imaginary), &step) in cubic.iter().zip(&imaginary_cubic).zip(&h) {
                let difference = std::array::from_fn(|i| real[i] - imaginary[i]);
                let minimum = segment_minimum(imaginary, step);
                let margin = segment_minimum(difference, step);
                if !minimum.is_finite() || !margin.is_finite() || minimum < 0.0 || margin < 0.0 {
                    return Err(error(
                        "KR0302",
                        "interpolated complex speed requires 0 <= Im(c) <= Re(c)",
                        "water_attenuation_db_per_wavelength",
                    ));
                }
            }
        }
        Ok(Self {
            case,
            layer,
            minimum_speed,
            cubic,
            imaginary_speeds,
            imaginary_cubic,
        })
    }

    pub(super) fn minimum_speed(&self) -> f64 {
        self.minimum_speed
    }

    pub(crate) fn minimum_scaled_difference(&self, other: &Self, scale: f64) -> f64 {
        if self.cubic.is_empty() {
            self.layer
                .points
                .iter()
                .zip(other.layer.points)
                .map(|(a, b)| a.sound_speed_mps - scale * b.sound_speed_mps)
                .fold(f64::INFINITY, f64::min)
        } else {
            self.cubic
                .iter()
                .zip(&other.cubic)
                .zip(self.layer.points.windows(2))
                .map(|((a, b), p)| {
                    segment_minimum(
                        std::array::from_fn(|i| a[i] - scale * b[i]),
                        p[1].depth_m - p[0].depth_m,
                    )
                })
                .fold(f64::INFINITY, f64::min)
        }
    }

    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    pub(super) fn mesh_speed(&self, index: usize, intervals: usize) -> f64 {
        let analytic = self.case.interpolation == Interpolation::AnalyticMunk;
        // munk.f90 computes 5000.0 / N in f32 before assigning it to f64 h.
        let step = if analytic {
            f64::from(5000.0_f32 / intervals as f32)
        } else {
            (self.layer.bottom - self.layer.top) / intervals as f64
        };
        let depth = self.layer.top + index as f64 * step;
        self.speed(if analytic {
            depth
        } else {
            depth.min(self.layer.bottom)
        })
    }

    #[allow(clippy::cast_precision_loss)]
    pub(super) fn mesh_complex_speed(&self, index: usize, intervals: usize) -> Complex64 {
        if self.imaginary_speeds.is_empty() {
            return Complex64::new(self.mesh_speed(index, intervals), 0.0);
        }
        self.complex_speed(
            (self.layer.top
                + index as f64 * ((self.layer.bottom - self.layer.top) / intervals as f64))
                .min(self.layer.bottom),
            false,
        )
    }

    #[allow(clippy::cast_precision_loss)]
    pub(crate) fn mesh_elastic_speed(&self, index: usize, intervals: usize) -> Complex64 {
        let depth = if index == intervals {
            self.layer.bottom
        } else {
            self.layer.top
                + index as f64 * ((self.layer.bottom - self.layer.top) / intervals as f64)
        };
        self.complex_speed(depth, true)
    }

    fn complex_speed(&self, depth: f64, elastic: bool) -> Complex64 {
        let points = self.layer.points;
        let upper = points
            .partition_point(|p| p.depth_m < depth)
            .clamp(1, points.len() - 1);
        let a = Complex64::new(
            points[upper - 1].sound_speed_mps,
            self.imaginary_speeds.get(upper - 1).copied().unwrap_or(0.0),
        );
        let b = Complex64::new(
            points[upper].sound_speed_mps,
            self.imaginary_speeds.get(upper).copied().unwrap_or(0.0),
        );
        let t = depth - points[upper - 1].depth_m;
        let weight = t / (points[upper].depth_m - points[upper - 1].depth_m);
        match self.case.interpolation {
            Interpolation::N2Linear => {
                // n2Linear takes each reciprocal before weighting. Dividing
                // the weights instead changes branch-sensitive deflated roots.
                let one = Complex64::new(1.0, 0.0);
                let top = one / a.powi(2);
                let bottom = one / b.powi(2);
                let z = (1.0 - weight) * top + weight * bottom;
                let root = if elastic
                    || (self.case.mode_solver == crate::ModeSolver::Kraken
                        && crate::solver::elastic::has_layers(self.case))
                {
                    principal_root(z)
                } else {
                    z.sqrt()
                };
                one / root
            }
            Interpolation::CLinear => (1.0 - weight) * a + weight * b,
            Interpolation::Pchip | Interpolation::Spline => {
                let [c0, c1, c2, c3] = self
                    .imaginary_cubic
                    .get(upper - 1)
                    .copied()
                    .unwrap_or([0.0; 4]);
                Complex64::new(self.speed(depth), c0 + t * (c1 + t * (c2 + t * c3)))
            }
            Interpolation::AnalyticMunk => unreachable!(),
        }
    }

    pub(super) fn speed(&self, depth: f64) -> f64 {
        if self.case.interpolation == Interpolation::AnalyticMunk {
            // Acoustics Toolbox v2023.5 misc/munk.f90, medium 1.
            let x = 2.0 * (depth - 1300.0) / 1300.0;
            // munk.f90 declares eps as f64 but initializes it from an f32 literal.
            return 1500.0 * (1.0 + f64::from(0.00737_f32) * (x - 1.0 + (-x).exp()));
        }
        let points = self.layer.points;
        let upper = points
            .partition_point(|p| p.depth_m < depth)
            .clamp(1, points.len() - 1);
        let a = points[upper - 1];
        let b = points[upper];
        let weight = (depth - a.depth_m) / (b.depth_m - a.depth_m);
        match self.case.interpolation {
            Interpolation::N2Linear
                if self.case.bottom_boundary.is_tabulated()
                    || self.case.surface_boundary.is_tabulated()
                    || (self.case.mode_solver == crate::ModeSolver::Kraken
                        && crate::solver::elastic::has_layers(self.case)) =>
            {
                // Preserve n2Linear's complex divisions for table-root seeding.
                let one = num_complex::Complex64::new(1.0, 0.0);
                let top = one / num_complex::Complex64::new(a.sound_speed_mps, 0.0).powi(2);
                let bottom = one / num_complex::Complex64::new(b.sound_speed_mps, 0.0).powi(2);
                (one / ((1.0 - weight) * top + weight * bottom).sqrt()).re
            }
            Interpolation::N2Linear => ((1.0 - weight) / a.sound_speed_mps.powi(2)
                + weight / b.sound_speed_mps.powi(2))
            .sqrt()
            .recip(),
            Interpolation::CLinear => {
                (1.0 - weight) * a.sound_speed_mps + weight * b.sound_speed_mps
            }
            Interpolation::Pchip | Interpolation::Spline => {
                let [c0, c1, c2, c3] = self.cubic[upper - 1];
                let t = depth - a.depth_m;
                c0 + t * (c1 + t * (c2 + t * c3))
            }
            Interpolation::AnalyticMunk => unreachable!(),
        }
    }
}

// The pinned Fortran runtime uses an algebraic principal square root, not the
// polar sin/cos formula in num-complex. Preserve it on real elastic shooting paths.
pub(crate) fn principal_root(z: Complex64) -> Complex64 {
    if z.im == 0.0 || !z.re.is_finite() || !z.im.is_finite() {
        return z.sqrt();
    }
    let mut norm = z.norm();
    let squared = norm * norm;
    if squared.is_normal() {
        // Correct platform hypot rounding using an FMA residual; N² shooting is
        // sensitive to a one-ulp square-root difference at the same complex input.
        let imaginary_squared = z.im * z.im;
        let residual = z.re.mul_add(z.re, -squared)
            + imaginary_squared
            + z.im.mul_add(z.im, -imaginary_squared)
            - norm.mul_add(norm, -squared);
        norm = residual.mul_add(0.5 / norm, norm);
    }
    let part = (0.5 * norm + 0.5 * z.re.abs()).sqrt();
    if z.re >= 0.0 {
        Complex64::new(part, 0.5 * z.im / part)
    } else {
        Complex64::new(0.5 * z.im.abs() / part, part.copysign(z.im))
    }
}

pub(crate) fn cubic_coefficients(
    kind: Interpolation,
    h: &[f64],
    values: &[f64],
) -> Result<Vec<[f64; 4]>, DiagnosticReport> {
    let delta: Vec<_> = values
        .windows(2)
        .zip(h)
        .map(|(p, &step)| (p[1] - p[0]) / step)
        .collect();
    if delta.iter().any(|v| !v.is_finite()) {
        return Err(invalid_profile());
    }
    let mut slopes = if kind == Interpolation::Pchip {
        pchip_slopes(h, &delta)
    } else {
        spline_slopes(h, &delta)
    };
    if kind == Interpolation::Pchip && h.len() > 1 {
        for (index, slope) in slopes.iter_mut().enumerate().skip(1).take(h.len() - 1) {
            *slope = project(delta[index - 1], delta[index], *slope);
        }
    }
    let cubic: Vec<_> = values
        .windows(2)
        .zip(h)
        .zip(slopes.windows(2))
        .map(|((p, &step), slope)| {
            let difference = p[1] - p[0];
            [
                p[0],
                slope[0],
                (3.0 * difference - step * (2.0 * slope[0] + slope[1])) / step.powi(2),
                (step * (slope[0] + slope[1]) - 2.0 * difference) / step.powi(3),
            ]
        })
        .collect();
    if cubic.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid_profile());
    }
    Ok(cubic)
}

// Cubic extrema can lie between mesh nodes; use them for validation and trapping.
pub(crate) fn segment_minimum([c0, c1, c2, c3]: [f64; 4], step: f64) -> f64 {
    let evaluate = |t: f64| c0 + t * (c1 + t * (c2 + t * c3));
    let end = evaluate(step);
    if !end.is_finite() {
        return f64::NAN;
    }
    let mut minimum = c0.min(end);
    let roots = if c3 == 0.0 {
        if c2 == 0.0 {
            return minimum;
        }
        [-c1 / (2.0 * c2), 0.0]
    } else {
        let discriminant = c2 * c2 - 3.0 * c3 * c1;
        if !discriminant.is_finite() {
            return f64::NAN;
        }
        if discriminant < 0.0 {
            return minimum;
        }
        let root = discriminant.sqrt();
        let q = -c2 - root.copysign(c2);
        if q == 0.0 {
            [0.0, 0.0]
        } else {
            [q / (3.0 * c3), c1 / q]
        }
    };
    for t in roots {
        if t > 0.0 && t < step {
            let value = evaluate(t);
            if !value.is_finite() {
                return f64::NAN;
            }
            minimum = minimum.min(value);
        }
    }
    minimum
}

fn invalid_profile() -> DiagnosticReport {
    error(
        "KR0302",
        "profile interpolation exceeds numeric range",
        "sound_speed_profile",
    )
}

// Thomas elimination; rows hold [lower, diagonal, upper, rhs].
fn tridiagonal(rows: &mut [[f64; 4]]) -> Vec<f64> {
    for i in 1..rows.len() {
        let factor = rows[i][0] / rows[i - 1][1];
        rows[i][1] -= factor * rows[i - 1][2];
        rows[i][3] -= factor * rows[i - 1][3];
    }
    let mut answer = vec![0.0; rows.len()];
    for i in (0..rows.len()).rev() {
        answer[i] =
            (rows[i][3] - rows[i][2] * answer.get(i + 1).copied().unwrap_or(0.0)) / rows[i][1];
    }
    answer
}

fn pchip_slopes(h: &[f64], delta: &[f64]) -> Vec<f64> {
    if h.len() == 1 {
        return vec![delta[0]; 2];
    }
    let end = |a: usize, b: usize| {
        let raw = ((2.0 * h[a] + h[b]) * delta[a] - h[a] * delta[b]) / (h[a] + h[b]);
        if delta[a] * raw <= 0.0 {
            0.0
        } else if delta[a] * delta[b] <= 0.0 && raw.abs() > (3.0 * delta[a]).abs() {
            3.0 * delta[a]
        } else {
            raw
        }
    };
    let first = end(0, 1);
    let last = end(h.len() - 1, h.len() - 2);
    clamped_slopes(h, delta, first, last)
}

fn project(left: f64, right: f64, slope: f64) -> f64 {
    if left * right > 0.0 {
        if left > 0.0 {
            slope.clamp(0.0, 3.0 * left.min(right))
        } else {
            slope.clamp(3.0 * left.max(right), 0.0)
        }
    } else {
        0.0
    }
}

fn clamped_slopes(h: &[f64], delta: &[f64], first: f64, last: f64) -> Vec<f64> {
    let n = h.len() + 1;
    let mut rows = Vec::with_capacity(n - 2);
    for j in 1..n - 1 {
        let (left, right) = (h[j - 1], h[j]);
        let mut rhs = 3.0 * (right * delta[j - 1] + left * delta[j]);
        if j == 1 {
            rhs -= right * first;
        }
        if j == n - 2 {
            rhs -= left * last;
        }
        rows.push([right, 2.0 * (left + right), left, rhs]);
    }
    let mut slopes = Vec::with_capacity(n);
    slopes.push(first);
    slopes.extend(tridiagonal(&mut rows));
    slopes.push(last);
    slopes
}

fn spline_slopes(h: &[f64], delta: &[f64]) -> Vec<f64> {
    let n = h.len() + 1;
    if n == 2 {
        return vec![delta[0]; 2];
    }
    if n == 3 {
        let curvature = (delta[1] - delta[0]) / (h[0] + h[1]);
        return vec![
            delta[0] - h[0] * curvature,
            delta[0] + h[0] * curvature,
            delta[1] + h[1] * curvature,
        ];
    }
    // Not-a-knot: eliminate the endpoint curvatures from the interior system.
    let mut rows = Vec::with_capacity(n - 2);
    for j in 1..n - 1 {
        let (left, right) = (h[j - 1], h[j]);
        let mut row = [
            left,
            2.0 * (left + right),
            right,
            6.0 * (delta[j] - delta[j - 1]),
        ];
        if j == 1 {
            row[1] += left * (1.0 + left / right);
            row[2] -= left * left / right;
        }
        if j == n - 2 {
            row[0] -= right * right / left;
            row[1] += right * (1.0 + right / left);
        }
        rows.push(row);
    }
    let inner = tridiagonal(&mut rows);
    let mut m = Vec::with_capacity(n);
    m.push((1.0 + h[0] / h[1]) * inner[0] - h[0] / h[1] * inner[1]);
    m.extend(inner);
    m.push((1.0 + h[n - 2] / h[n - 3]) * m[n - 2] - h[n - 2] / h[n - 3] * m[n - 3]);
    let mut slopes: Vec<_> = (0..n - 1)
        .map(|i| delta[i] - h[i] * (2.0 * m[i] + m[i + 1]) / 6.0)
        .collect();
    slopes.push(delta[n - 2] + h[n - 2] * (m[n - 2] + 2.0 * m[n - 1]) / 6.0);
    slopes
}

#[cfg(test)]
mod tests {
    use super::{Profile, principal_root, spline_slopes};

    #[test]
    fn elastic_square_root_keeps_pinned_n2_rounding() {
        let root = principal_root(num_complex::Complex64::new(
            1.103_288_297_121_612_5e-7,
            -1.274_145_729_606_203_5e-9,
        ));
        assert_eq!(root.re.to_bits(), 3.321_633_758_344_346e-4_f64.to_bits());
        assert_eq!(root.im.to_bits(), (-1.917_950_355_612_498e-6_f64).to_bits());
    }
    use crate::{Case, Interpolation, legacy::load_case};
    use std::path::Path;

    #[test]
    fn analytic_mesh_sampling_preserves_the_reference_f32_step() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/MunkAnalytic");
        let case = load_case(root.with_extension("env"), root.with_extension("flp")).unwrap();
        let profile = Profile::new(&case).unwrap();
        let rounded_depth = 5000.0 * f64::from(5000.0_f32 / 5000.0_f32);
        assert_eq!(
            profile.mesh_speed(5000, 5000).to_bits(),
            profile.speed(rounded_depth).to_bits()
        );
        let rounded_depth = 3333.0 * f64::from(5000.0_f32 / 3333.0_f32);
        assert_eq!(
            profile.mesh_speed(3333, 3333).to_bits(),
            profile.speed(rounded_depth).to_bits()
        );
        assert_ne!(
            profile.mesh_speed(3333, 3333).to_bits(),
            profile.speed(5000.0).to_bits()
        );
    }

    #[test]
    fn short_spline_and_pchip_profiles() {
        assert_eq!(spline_slopes(&[1.0], &[3.0]), [3.0; 2]);
        assert_eq!(spline_slopes(&[1.0, 1.0], &[1.0, 3.0]), [0.0, 2.0, 4.0]);
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/Pekeris");
        let mut input = load_case(root.with_extension("env"), root.with_extension("flp"))
            .unwrap()
            .into_definition();
        for kind in [Interpolation::Pchip, Interpolation::Spline] {
            input.interpolation = kind;
            let case = Case::from_definition(input.clone()).unwrap();
            assert!((Profile::new(&case).unwrap().speed(50.0) - 1500.0).abs() < 1e-12);
        }
        input.sound_speed_profile.insert(
            1,
            crate::SoundSpeedPoint {
                depth_m: 40.0,
                sound_speed_mps: 1600.0,
            },
        );
        for kind in [Interpolation::Pchip, Interpolation::Spline] {
            input.interpolation = kind;
            let case = Case::from_definition(input.clone()).unwrap();
            let profile = Profile::new(&case).unwrap();
            for point in &input.sound_speed_profile {
                assert!((profile.speed(point.depth_m) - point.sound_speed_mps).abs() < 1e-10);
            }
            assert!((1500.0..1600.0).contains(&profile.speed(20.0)));
        }
        input.interpolation = Interpolation::Spline;
        input.sound_speed_profile = [0.0, 23.0, 61.0, 100.0]
            .map(|depth_m| crate::SoundSpeedPoint {
                depth_m,
                sound_speed_mps: 1500.0 + 2.0 * depth_m - 0.03 * depth_m.powi(2)
                    + 0.0001 * depth_m.powi(3),
            })
            .to_vec();
        let case = Case::from_definition(input.clone()).unwrap();
        let profile = Profile::new(&case).unwrap();
        let depth = 40.0_f64;
        let exact = 1500.0 + 2.0 * depth - 0.03 * depth.powi(2) + 0.0001 * depth.powi(3);
        assert!((profile.speed(depth) - exact).abs() < 1e-10);

        // A three-point spline can go negative between positive knots on a
        // narrow interval: reject it even when no mesh point samples the dip.
        input.sound_speed_profile = [(0.0, 1500.0), (1.0, 1500.0), (100.0, 1e8)]
            .map(|(depth_m, sound_speed_mps)| crate::SoundSpeedPoint {
                depth_m,
                sound_speed_mps,
            })
            .to_vec();
        assert_eq!(
            Case::from_definition(input).unwrap_err().diagnostics()[0].field,
            "sound_speed_profile"
        );
    }
}
