// Adapted from Acoustics Toolbox v2023.5 misc/pchipMod.f90 and splinec.f90,
// Copyright (C) 2009 Michael B. Porter. GPL-3.0-or-later; see LICENSE.
// Real, constant-density water-column interpolation only.
use crate::solver::error;
use crate::{CaseDefinition, DiagnosticReport, Interpolation};

pub(super) struct Profile<'a> {
    case: &'a CaseDefinition,
    minimum_speed: f64,
    // Local polynomial coefficients, c0 + t * (c1 + t * (c2 + t * c3)).
    cubic: Vec<[f64; 4]>,
}

impl<'a> Profile<'a> {
    pub(super) fn new(case: &'a CaseDefinition) -> Result<Self, DiagnosticReport> {
        let points = &case.sound_speed_profile;
        let mut minimum_speed = points
            .iter()
            .map(|point| point.sound_speed_mps)
            .fold(f64::INFINITY, f64::min);
        let mut cubic = Vec::new();
        if matches!(
            case.interpolation,
            Interpolation::Pchip | Interpolation::Spline
        ) {
            let h: Vec<_> = points
                .windows(2)
                .map(|p| p[1].depth_m - p[0].depth_m)
                .collect();
            let delta: Vec<_> = points
                .windows(2)
                .zip(&h)
                .map(|(p, &step)| (p[1].sound_speed_mps - p[0].sound_speed_mps) / step)
                .collect();
            if delta.iter().any(|v| !v.is_finite()) {
                return Err(invalid_profile());
            }
            let mut slopes = if case.interpolation == Interpolation::Pchip {
                pchip_slopes(&h, &delta)
            } else {
                spline_slopes(&h, &delta)
            };
            if case.interpolation == Interpolation::Pchip && h.len() > 1 {
                for (index, slope) in slopes.iter_mut().enumerate().skip(1).take(h.len() - 1) {
                    *slope = project(delta[index - 1], delta[index], *slope);
                }
            }
            cubic = points
                .windows(2)
                .zip(&h)
                .zip(slopes.windows(2))
                .map(|((pair, &step), slope)| {
                    let difference = pair[1].sound_speed_mps - pair[0].sound_speed_mps;
                    [
                        pair[0].sound_speed_mps,
                        slope[0],
                        (3.0 * difference - step * (2.0 * slope[0] + slope[1])) / step.powi(2),
                        (step * (slope[0] + slope[1]) - 2.0 * difference) / step.powi(3),
                    ]
                })
                .collect();
            if cubic.iter().flatten().any(|v| !v.is_finite()) {
                return Err(invalid_profile());
            }
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
        Ok(Self {
            case,
            minimum_speed,
            cubic,
        })
    }

    pub(super) fn minimum_speed(&self) -> f64 {
        self.minimum_speed
    }

    pub(super) fn speed(&self, depth: f64) -> f64 {
        let points = &self.case.sound_speed_profile;
        let upper = points
            .partition_point(|p| p.depth_m < depth)
            .clamp(1, points.len() - 1);
        let a = points[upper - 1];
        let b = points[upper];
        let weight = (depth - a.depth_m) / (b.depth_m - a.depth_m);
        match self.case.interpolation {
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
        }
    }
}

// Cubic extrema can lie between mesh nodes; use them for validation and trapping.
fn segment_minimum([c0, c1, c2, c3]: [f64; 4], step: f64) -> f64 {
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
    use super::{Profile, spline_slopes};
    use crate::{Case, Interpolation, legacy::load_case};
    use std::path::Path;

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
