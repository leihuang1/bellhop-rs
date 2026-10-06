// Adapted from Acoustics Toolbox v2023.5 misc/RefCoef.f90, PolyMod.f90
// and Kraken/BCImpedancecMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE.
//! Pinned tabulated bottom semantics.
use crate::{BottomBoundary, CaseDefinition, Interpolation, ModeSolver, SurfaceBoundary};
use num_complex::Complex64;

pub(crate) fn validate(case: &CaseDefinition) -> Result<(), &'static str> {
    if !case.bottom_boundary.is_tabulated() {
        return Ok(());
    }
    if case.mode_solver != ModeSolver::Krakenc
        || !matches!(
            case.interpolation,
            Interpolation::N2Linear | Interpolation::CLinear
        )
        || !case.additional_fluid_layers.is_empty()
        || case.surface_boundary != SurfaceBoundary::Vacuum
        || case.max_range_m != 0.0
        || case.mesh_reference_frequency_hz.is_some()
        || case
            .water_attenuation_db_per_wavelength
            .iter()
            .any(|&a| a != 0.0)
    {
        return Err(
            "tabulated bottoms require single-frequency KRAKENC lossless N/C water, vacuum top and RMax=0 (no B option)",
        );
    }
    validate_table(&case.bottom_boundary, case.frequency_hz)
}

pub(crate) fn validate_surface(case: &CaseDefinition) -> Result<(), &'static str> {
    if matches!(case.surface_boundary, SurfaceBoundary::Impedance { .. }) {
        return Err("top P/IRC is not supported");
    }
    if !case.surface_boundary.is_tabulated() {
        return Ok(());
    }
    if case.mode_solver != ModeSolver::Krakenc
        || !case.additional_fluid_layers.is_empty()
        || !matches!(
            case.interpolation,
            Interpolation::N2Linear | Interpolation::CLinear
        )
        || case.max_range_m != 0.0
        || case.mesh_reference_frequency_hz.is_some()
        || case
            .water_attenuation_db_per_wavelength
            .iter()
            .any(|&a| a != 0.0)
    {
        return Err(
            "TRC requires single-frequency KRAKENC lossless N/C water and RMax=0 (no B option)",
        );
    }
    if case
        .sound_speed_profile
        .last()
        .is_some_and(|p| case.c_low_mps < p.sound_speed_mps)
    {
        return Err(
            "TRC requires cLow at or above the reference inside-water speed (last SSP node); evanescent top-table roots are not supported",
        );
    }
    validate_table(&case.surface_boundary, case.frequency_hz)
}

#[allow(clippy::float_cmp)]
pub(crate) fn validate_table(
    boundary: &BottomBoundary,
    solve_frequency_hz: f64,
) -> Result<(), &'static str> {
    match boundary {
        BottomBoundary::Reflection(points) => {
            if !(2..=100_000).contains(&points.len())
                || points.iter().any(|p| {
                    !p.angle_degrees.is_finite()
                        || !(0.0..=90.0).contains(&p.angle_degrees)
                        || !p.magnitude.is_finite()
                        || p.magnitude < 0.0
                        || !p.phase_radians.is_finite()
                })
                || points
                    .windows(2)
                    .any(|p| p[0].angle_degrees >= p[1].angle_degrees)
            {
                return Err(
                    "BRC requires 2..=100000 ordered finite angles in 0..90, nonnegative magnitudes and finite unwrapped phases",
                );
            }
        }
        BottomBoundary::Impedance {
            frequency_hz,
            points,
        } => {
            if *frequency_hz != solve_frequency_hz {
                return Err("IRC frequency must equal the solve frequency");
            }
            if !(2..=100_000).contains(&points.len())
                || points.iter().any(|p| {
                    !p.wavenumber_squared.is_finite()
                        || p.wavenumber_squared < 0.0
                        || !p.f.re.is_finite()
                        || !p.f.im.is_finite()
                        || !p.g.re.is_finite()
                        || !p.g.im.is_finite()
                        || p.power.abs_diff(0) > 1000
                })
                || points
                    .windows(2)
                    .any(|p| p[0].wavenumber_squared >= p[1].wavenumber_squared)
            {
                return Err(
                    "IRC requires 2..=100000 ordered finite nonnegative squared wavenumbers, finite f/g and bounded powers",
                );
            }
            // Only neighbouring (up to three) points are rescaled together.
            if points
                .windows(2)
                .any(|p| p[0].power.abs_diff(p[1].power) > 100)
                || points.windows(3).any(|p| {
                    p.iter().map(|q| q.power).max().unwrap()
                        - p.iter().map(|q| q.power).min().unwrap()
                        > 100
                })
            {
                return Err("IRC neighbouring scale differences exceed 100");
            }
        }
        _ => {}
    }
    Ok(())
}

#[allow(clippy::cast_possible_truncation, clippy::many_single_char_names)]
pub(crate) fn impedance(
    boundary: &BottomBoundary,
    x: Complex64,
    water_k2: f64,
) -> (Complex64, Complex64, i32) {
    match boundary {
        BottomBoundary::Reflection(points) => {
            let kx = x.sqrt();
            let kz = (water_k2 - x).sqrt();
            let angle = kz.re.atan2(kx.re).to_degrees();
            // RefCoef uses REAL(theta) for bracket selection, double theta for weights.
            let query = f64::from(angle as f32);
            let r = if query < points[0].angle_degrees
                || query > points.last().unwrap().angle_degrees
            {
                Complex64::new(0.0, 0.0)
            } else {
                let i = points
                    .partition_point(|p| p.angle_degrees <= query)
                    .clamp(1, points.len() - 1)
                    - 1;
                let a = &points[i];
                let b = &points[i + 1];
                let weight = (angle - a.angle_degrees) / (b.angle_degrees - a.angle_degrees);
                Complex64::from_polar(
                    (1.0 - weight) * a.magnitude + weight * b.magnitude,
                    (1.0 - weight) * a.phase_radians + weight * b.phase_radians,
                )
            };
            (
                Complex64::new(1.0, 0.0),
                (1.0 + r) / (Complex64::i() * kz * (1.0 - r)),
                0,
            )
        }
        BottomBoundary::Impedance { points, .. } => {
            let endpoint = if x.re < points[0].wavenumber_squared {
                Some(&points[0])
            } else if x.re > points.last().unwrap().wavenumber_squared {
                points.last()
            } else {
                None
            };
            if let Some(p) = endpoint {
                return (p.f, p.g, p.power);
            }
            let left = points
                .partition_point(|p| p.wavenumber_squared <= x.re)
                .clamp(1, points.len() - 1)
                - 1;
            let slice = &points[left..(left + 3).min(points.len())];
            let mut f = [Complex64::new(0.0, 0.0); 3];
            let mut g = f;
            let mut h = f;
            for (i, p) in slice.iter().enumerate() {
                h[i] = p.wavenumber_squared - x;
                let scale = 10_f64.powi(p.power - slice[0].power);
                f[i] = p.f * scale;
                g[i] = p.g * scale;
            }
            for order in 1..slice.len() {
                for j in 0..slice.len() - order {
                    for values in [&mut f, &mut g] {
                        let a = h[j + order];
                        let b = h[j];
                        let left = values[j];
                        let right = values[j + 1];
                        // Preserve PolyZ's GNU Fortran 12.2 -ffast-math grouping.
                        let numerator = polynomial_numerator(a, left, b, right);
                        values[j] = numerator / (a - b);
                    }
                }
            }
            (f[0], g[0], slice[0].power)
        }
        _ => unreachable!(),
    }
}

fn polynomial_numerator(
    a: Complex64,
    left: Complex64,
    b: Complex64,
    right: Complex64,
) -> Complex64 {
    Complex64::new(
        (a.re * left.re + b.im * right.im) - (a.im * left.im + b.re * right.re),
        (a.re * left.im + a.im * left.re) - (b.re * right.im + b.im * right.re),
    )
}

#[cfg(test)]
mod tests {
    use super::{impedance, polynomial_numerator};
    use crate::{BottomBoundary, ImpedancePoint, ReflectionPoint};
    use num_complex::Complex64;

    #[test]
    fn irc_scales_selects_by_real_part_and_evaluates_complex_polynomials() {
        let points = [(1.0, 0), (2.0, 10), (3.0, 20)]
            .map(|(x, power)| ImpedancePoint {
                wavenumber_squared: x,
                f: Complex64::new(x * x, 0.0) * 10_f64.powi(-power),
                g: Complex64::new(1.0, x) * 10_f64.powi(-power),
                power,
            })
            .to_vec();
        let boundary = BottomBoundary::Impedance {
            frequency_hz: 50.0,
            points,
        };
        let x = Complex64::new(1.5, 0.2);
        let (f, g, power) = impedance(&boundary, x, 0.0);
        assert_eq!(power, 0);
        assert!((f - x * x).norm() < 1e-12);
        assert!((g - (1.0 + Complex64::i() * x)).norm() < 1e-12);
        let x = Complex64::new(2.5, 0.2);
        let (f, g, power) = impedance(&boundary, x, 0.0);
        assert_eq!(power, 10); // final interval uses only two points
        assert!((f * 1e10 - (5.0 * x - 6.0)).norm() < 1e-12);
        assert!((g * 1e10 - (1.0 + Complex64::i() * x)).norm() < 1e-12);
        let (f, g, power) = impedance(&boundary, Complex64::new(0.5, 2.0), 0.0);
        assert_eq!(
            (f, g, power),
            (Complex64::new(1.0, 0.0), Complex64::new(1.0, 1.0), 0)
        );
        let (f, g, power) = impedance(&boundary, Complex64::new(4.0, -2.0), 0.0);
        assert_eq!(power, 20);
        assert!((f * 1e20 - 9.0).norm() < 1e-12);
        assert!((g * 1e20 - Complex64::new(1.0, 3.0)).norm() < 1e-12);
    }

    #[test]
    fn brc_interpolates_unwrapped_phase_and_zeroes_outside_the_domain() {
        let mut points = vec![
            ReflectionPoint {
                angle_degrees: 10.0,
                magnitude: 0.5,
                phase_radians: 170_f64.to_radians(),
            },
            ReflectionPoint {
                angle_degrees: 30.0,
                magnitude: 1.0,
                phase_radians: 190_f64.to_radians(),
            },
        ];
        let x = Complex64::new(20_f64.to_radians().cos().powi(2), 0.0);
        let (_, g, _) = impedance(&BottomBoundary::Reflection(points.clone()), x, 1.0);
        let ratio = Complex64::i() * (1.0 - x).sqrt() * g;
        let reflection = (ratio - 1.0) / (ratio + 1.0);
        assert!((reflection - Complex64::new(-0.75, 0.0)).norm() < 1e-12);
        // The double angle is inside, but its f32 bracket query is outside.
        points[0].angle_degrees = 10.000_000_05;
        let x = Complex64::new(10.000_000_1_f64.to_radians().cos().powi(2), 0.0);
        let (_, g, _) = impedance(&BottomBoundary::Reflection(points), x, 1.0);
        assert!((g - 1.0 / (Complex64::i() * (1.0 - x).sqrt())).norm() < 1e-12);
    }

    #[test]
    fn polyz_preserves_pinned_real_component_grouping() {
        let numerator = polynomial_numerator(
            Complex64::new(1.0, 1.0),
            Complex64::new(1.0, 0.3),
            Complex64::new(0.7, 0.1),
            Complex64::new(1.0, 1.0),
        );
        // (1 + 0.1) - (0.3 + 0.7), not (1 - 0.3) - (0.7 - 0.1).
        assert_eq!(
            numerator.re.to_bits(),
            0.100_000_000_000_000_09_f64.to_bits()
        );
    }
}
