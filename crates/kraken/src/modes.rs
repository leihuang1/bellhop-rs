// Adapted from Acoustics Toolbox v2023.5 Kraken/kraken.f90 and
// Kraken/InverseIterationMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE. Distributed without warranty.
//!
//! Real, layered-fluid finite-difference path from KRAKEN v2023.5.
//! Sturm counts isolate modes; inverse iteration samples the first mesh;
//! Richardson extrapolation refines eigenvalues only, as in the reference.
use crate::profile::Profile;
use crate::solver::error;
use crate::{
    BottomBoundary, Case, DiagnosticReport, MAX_MODE_LIMIT, ModeSet, NormalMode, SurfaceBoundary,
};
use num_complex::Complex64;
use std::f64::consts::PI;

const MAX_SHAPES: usize = 5_000_000;
// ponytail: 2.5b conservative work covers BroadBand/MunkK at 500 Hz (2.28b);
// use measured root work or faster isolation if larger spectra are required.
const MAX_WORK: usize = 2_500_000_000;
const ROOT_STEPS: usize = 64;

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]
pub(super) fn solve(case: &Case) -> Result<ModeSet, DiagnosticReport> {
    let omega = 2.0 * PI * case.frequency_hz;
    let (bottom_k2, bottom_complex_k2) = if case.bottom_boundary.is_half_space() {
        // AttenMod::CRCI converts dB/wavelength to a positive imaginary sound speed.
        let bottom_c = Complex64::new(
            case.bottom_sound_speed_mps,
            case.bottom_attenuation_db_per_wavelength * case.bottom_sound_speed_mps
                / (8.685_889_6 * 2.0 * PI),
        );
        (
            (omega / case.bottom_sound_speed_mps).powi(2),
            (Complex64::new(omega, 0.0) / bottom_c).powi(2),
        )
    } else {
        (0.0, Complex64::new(0.0, 0.0))
    };
    let profiles = crate::layers::iter(case)
        .map(|layer| Profile::new_layer(case, layer))
        .collect::<Result<Vec<_>, _>>()?;
    if !omega.is_finite()
        || !bottom_k2.is_finite()
        || !bottom_complex_k2.re.is_finite()
        || !bottom_complex_k2.im.is_finite()
    {
        return Err(error(
            "KR0302",
            "frequency/profile requires an unsupported mesh",
            "mesh_points",
        ));
    }
    let mut table: Vec<Vec<f64>> = Vec::new();
    let mut modes = Vec::new();
    let mut work = 0_usize;
    for set in 0..5 {
        let multiplier = 1 << set;
        let layers = crate::layers::mesh_layers(case, multiplier)?;
        let mesh = Mesh::new(case, &profiles, layers, omega, bottom_k2, bottom_complex_k2)?;
        let roots = mesh.roots(&mut work)?;
        if set == 0 {
            modes = roots
                .iter()
                .map(|&x| mesh.mode(x))
                .collect::<Result<Vec<_>, _>>()?;
        } else if roots.len() != modes.len() {
            return Err(error(
                "KR0303",
                format!(
                    "mode count changed during mesh refinement ({} -> {} at multiplier {multiplier}); move spectral limits away from roots",
                    modes.len(),
                    roots.len()
                ),
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
    surface_k2: f64,
    surface_complex_k2: Complex64,
    b1: Vec<f64>,
    b1c: Vec<f64>,
    layers: Vec<crate::layers::MeshLayer>,
    min_speed: f64,
}

impl<'a> Mesh<'a> {
    #[allow(clippy::cast_precision_loss)]
    fn new(
        case: &'a Case,
        profiles: &[Profile<'_>],
        layers: Vec<crate::layers::MeshLayer>,
        omega: f64,
        bottom_k2: f64,
        bottom_complex_k2: Complex64,
    ) -> Result<Self, DiagnosticReport> {
        let h = layers[0].h;
        let mut min_speed = f64::INFINITY;
        let mut invalid_speed = false;
        let mut b1c = Vec::new();
        let mut b1 = Vec::new();
        for ((layer, profile), material) in
            layers.iter().zip(profiles).zip(crate::layers::iter(case))
        {
            let water_loss = material.loss.iter().any(|&a| a != 0.0);
            for i in 0..=layer.intervals {
                let speed = profile.mesh_complex_speed(i, layer.intervals);
                invalid_speed |= !speed.re.is_finite() || speed.re <= 0.0 || !speed.im.is_finite();
                min_speed = min_speed.min(speed.re);
                if water_loss {
                    let k2 = Complex64::new(omega * omega, 0.0) / speed.powi(2);
                    b1c.push(k2.im);
                    b1.push(-2.0 + layer.h * layer.h * k2.re);
                } else {
                    b1c.push(0.0);
                    b1.push(-2.0 + layer.h * layer.h * (omega * omega / (speed.re * speed.re)));
                }
            }
        }
        if invalid_speed {
            return Err(error(
                "KR0302",
                "interpolated sound speed is non-finite or non-positive",
                "sound_speed_profile",
            ));
        }
        if h <= 0.0 || !h.is_finite() || b1.iter().chain(&b1c).any(|x| !x.is_finite()) {
            return Err(error(
                "KR0302",
                "mesh coefficients exceed numeric range",
                "mesh_points",
            ));
        }
        let (surface_k2, surface_complex_k2) = if case.surface_boundary.is_half_space() {
            let c = Complex64::new(
                case.surface_sound_speed_mps,
                case.surface_attenuation_db_per_wavelength * case.surface_sound_speed_mps
                    / (8.685_889_6 * 2.0 * PI),
            );
            (
                (omega / case.surface_sound_speed_mps).powi(2),
                (Complex64::new(omega, 0.0) / c).powi(2),
            )
        } else {
            (0.0, Complex64::new(0.0, 0.0))
        };
        if !surface_k2.is_finite()
            || !surface_complex_k2.re.is_finite()
            || !surface_complex_k2.im.is_finite()
        {
            return Err(error(
                "KR0302",
                "surface impedance exceeds numeric range",
                "surface_boundary",
            ));
        }
        Ok(Self {
            case,
            h,
            omega,
            bottom_k2,
            bottom_complex_k2,
            surface_k2,
            surface_complex_k2,
            b1,
            b1c,
            layers,
            min_speed,
        })
    }

    fn bottom_gamma(&self, x: f64) -> Complex64 {
        (Complex64::new(x, 0.0) - self.bottom_complex_k2).sqrt()
    }

    fn bottom_admittance(&self, x: f64) -> f64 {
        if matches!(
            self.case.bottom_boundary,
            BottomBoundary::ElasticHalfSpace { .. }
        ) {
            let (f, g) = crate::elastic::half_space(
                &self.case.bottom_boundary,
                x.into(),
                self.omega,
                self.case.bottom_sound_speed_mps,
                self.case.bottom_density_g_cm3,
                self.case.bottom_attenuation_db_per_wavelength,
                true,
            );
            (f / g).re
        } else if self.case.bottom_boundary.is_half_space() {
            self.bottom_gamma(x).re / self.case.bottom_density_g_cm3
        } else {
            0.0
        }
    }

    fn surface_admittance(&self, x: f64) -> f64 {
        if self.case.surface_boundary == SurfaceBoundary::FluidHalfSpace {
            self.surface_gamma(x).re / self.case.surface_density_g_cm3
        } else {
            0.0
        }
    }

    fn bottom_diagonal(&self, x: f64) -> f64 {
        let diagonal = (self.b1.last().unwrap() - self.h * self.h * x) * 0.5;
        if self.case.bottom_boundary == BottomBoundary::FluidHalfSpace {
            diagonal
                - self.h * self.case.water_density_g_cm3 / self.case.bottom_density_g_cm3
                    * self.bottom_gamma(x).re
        } else {
            diagonal - self.h * self.case.water_density_g_cm3 * self.bottom_admittance(x)
        }
    }

    fn surface_gamma(&self, x: f64) -> Complex64 {
        (Complex64::new(x, 0.0) - self.surface_complex_k2).sqrt()
    }

    fn surface_diagonal(&self, x: f64) -> f64 {
        let diagonal = (self.b1[0] - self.h * self.h * x) * 0.5;
        if self.case.surface_boundary == SurfaceBoundary::FluidHalfSpace {
            diagonal
                - self.h * self.case.water_density_g_cm3 / self.case.surface_density_g_cm3
                    * self.surface_gamma(x).re
        } else {
            diagonal - self.h * self.case.water_density_g_cm3 * self.surface_admittance(x)
        }
    }

    fn layer_diagonal(&self, x: f64, medium: usize, i: usize) -> f64 {
        let layer = &self.layers[medium];
        let diagonal = (self.b1[layer.coefficient_start + i] - layer.h * layer.h * x)
            / (layer.h * layer.density);
        if medium == 0 && i == 0 {
            diagonal * 0.5 - self.surface_admittance(x)
        } else if i == layer.intervals {
            if let Some(next) = self.layers.get(medium + 1) {
                diagonal.midpoint(
                    (self.b1[next.coefficient_start] - next.h * next.h * x)
                        / (next.h * next.density),
                )
            } else {
                diagonal * 0.5 - self.bottom_admittance(x)
            }
        } else {
            diagonal
        }
    }

    // Inertia of the symmetric tridiagonal acoustic operator A(x): number of roots above x.
    fn count(&self, x: f64) -> usize {
        if crate::elastic::has_half_space(self.case) {
            return self.elastic_count(x);
        }
        if self.layers.len() > 1 {
            let mut pivot: f64 = 0.0;
            let mut count = 0;
            let mut first = true;
            for (medium, layer) in self.layers.iter().enumerate() {
                let start = usize::from(
                    medium > 0 || self.case.surface_boundary == SurfaceBoundary::Vacuum,
                );
                let end = layer.intervals + 1
                    - usize::from(
                        medium + 1 == self.layers.len()
                            && self.case.bottom_boundary == BottomBoundary::Vacuum,
                    );
                for i in start..end {
                    let diagonal = self.layer_diagonal(x, medium, i);
                    let e = 1.0 / (layer.h * layer.density);
                    pivot = if first {
                        diagonal
                    } else {
                        diagonal - e * e / pivot
                    };
                    first = false;
                    if pivot.abs() < 1e-30 {
                        pivot = -1e-30;
                    }
                    if pivot > 0.0 {
                        count += 1;
                    }
                }
            }
            return count;
        }
        let shift = self.h * self.h * x;
        let first = usize::from(self.case.surface_boundary == SurfaceBoundary::Vacuum);
        let end = self.b1.len() - usize::from(self.case.bottom_boundary == BottomBoundary::Vacuum);
        let mut pivot = 0.0;
        let mut count = 0;
        for i in first..end {
            let d = if i == 0 {
                self.surface_diagonal(x)
            } else if i + 1 == self.b1.len() {
                self.bottom_diagonal(x)
            } else {
                self.b1[i] - shift
            };
            pivot = if i == first { d } else { d - 1.0 / pivot };
            if pivot.abs() < 1e-30 {
                pivot = -1e-30;
            }
            if pivot > 0.0 {
                count += 1;
            }
        }
        count
    }

    // AcousticLayers/FUNCT includes the elastic impedance poles in ModeCount;
    // acoustic-matrix inertia alone misses the Scholte/interface branch.
    fn elastic_count(&self, x: f64) -> usize {
        let mut count = 0;
        let impedance = |boundary: &BottomBoundary, cp, density: f64, loss, k2: Complex64| {
            if matches!(boundary, BottomBoundary::ElasticHalfSpace { .. }) {
                crate::elastic::half_space(boundary, x.into(), self.omega, cp, density, loss, true)
            } else if boundary.is_half_space() {
                (
                    (Complex64::new(x, 0.0) - k2).sqrt().re.into(),
                    density.into(),
                )
            } else if *boundary == BottomBoundary::Vacuum {
                (1.0.into(), 0.0.into())
            } else {
                (0.0.into(), 1.0.into())
            }
        };
        let (f, g) = impedance(
            &self.case.bottom_boundary,
            self.case.bottom_sound_speed_mps,
            self.case.bottom_density_g_cm3,
            self.case.bottom_attenuation_db_per_wavelength,
            self.bottom_complex_k2,
        );
        if matches!(
            self.case.bottom_boundary,
            BottomBoundary::ElasticHalfSpace { .. }
        ) && g.re > 0.0
        {
            count += 1;
        }
        let (mut f, mut g) = (f.re, g.re);
        for layer in self.layers.iter().rev() {
            let shift = layer.h * layer.h * x;
            let mut p0 = 0.0;
            let mut p1 = -2.0 * g;
            let mut p2 = (self.b1[layer.coefficient_start + layer.intervals] - shift) * g
                - 2.0 * layer.h * f * layer.density;
            for &coefficient in self.b1
                [layer.coefficient_start..layer.coefficient_start + layer.intervals]
                .iter()
                .rev()
            {
                p0 = p1;
                p1 = p2;
                p2 = (shift - coefficient) * p1 - p0;
                if p0 * p1 <= 0.0 {
                    count += 1;
                }
                if p2.abs() > 1e50 {
                    p0 *= 1e-50;
                    p1 *= 1e-50;
                    p2 *= 1e-50;
                }
            }
            f = -(p2 - p0) / (2.0 * layer.h) / layer.density;
            g = -p1;
        }
        let (ft, gt) = impedance(
            &self.case.surface_boundary,
            self.case.surface_sound_speed_mps,
            self.case.surface_density_g_cm3,
            self.case.surface_attenuation_db_per_wavelength,
            self.surface_complex_k2,
        );
        if matches!(
            self.case.surface_boundary,
            BottomBoundary::ElasticHalfSpace { .. }
        ) && gt.re > 0.0
        {
            count += 1;
        }
        let delta = -f * gt.re - g * ft.re;
        if g * delta > 0.0 {
            count += 1;
        }
        count
    }

    #[allow(clippy::float_cmp)]
    fn roots(&self, work: &mut usize) -> Result<Vec<f64>, DiagnosticReport> {
        // Preserve Solve1's lower search guard, including its cutoff exclusion.
        let elastic = crate::elastic::has_half_space(self.case);
        let c_high = if elastic {
            crate::elastic::maximum_speed(self.case)
        } else {
            self.case.c_high_mps
        };
        let minimum = if elastic {
            crate::elastic::minimum_speed(self.case, self.min_speed)
        } else {
            self.min_speed
        };
        let low = 1.00001 * (self.omega / c_high).powi(2);
        let mut high = (self.omega / self.case.c_low_mps.max(minimum)).powi(2);
        if self.case.surface_boundary == SurfaceBoundary::Rigid
            && self.case.bottom_boundary == BottomBoundary::Rigid
            && self.case.c_low_mps <= self.min_speed
            && crate::layers::iter(self.case)
                .flat_map(|layer| layer.points)
                .all(|point| point.sound_speed_mps == self.min_speed)
        {
            // The constant-profile plane eigenvalue can round just above omega²/c²
            // in the finite-difference diagonal; keep the inclusive endpoint.
            high = high.max((self.b1[0] + 2.0) / (self.h * self.h)).next_up();
        }
        if low >= high
            || !high.is_finite()
            || (self.case.bottom_boundary == BottomBoundary::FluidHalfSpace
                && low <= self.bottom_k2)
            || (self.case.surface_boundary == SurfaceBoundary::FluidHalfSpace
                && low <= self.surface_k2)
        {
            return Err(error(
                "KR0301",
                "phase-speed limits contain no supported modes",
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
        clippy::many_single_char_names,
        clippy::too_many_lines
    )]
    fn mode(&self, x: f64) -> Result<NormalMode, DiagnosticReport> {
        let last = self.layers.last().unwrap();
        let n = last.node_start + last.intervals + 1;
        let h_rho = self.h * self.case.water_density_g_cm3;
        let shift = self.h * self.h * x;
        let mut d = Vec::with_capacity(n);
        let mut e = vec![0.0; n + 1];
        for (medium, layer) in self.layers.iter().enumerate() {
            for i in usize::from(medium > 0)..=layer.intervals {
                d.push(if self.layers.len() == 1 {
                    (self.b1[i] - shift) / h_rho
                } else {
                    self.layer_diagonal(x, medium, i)
                });
                if i > 0 {
                    e[layer.node_start + i] = 1.0 / (layer.h * layer.density);
                }
            }
        }
        e[n] = 1.0 / (last.h * last.density);
        if self.case.surface_boundary == SurfaceBoundary::Vacuum {
            d[0] = 1.0;
            e[1] = 0.0;
        } else if self.layers.len() == 1 {
            d[0] = self.surface_diagonal(x) / h_rho;
        }
        if self.case.bottom_boundary == BottomBoundary::Vacuum {
            d[n - 1] = 1.0;
            e[n - 1] = 0.0;
        } else if self.layers.len() == 1 {
            d[n - 1] = self.bottom_diagonal(x) / h_rho;
        }
        let mut phi = inverse_iteration(&d, &e)?;
        let mut norm = 0.0;
        let mut slow = 0.0;
        let mut volume_loss = 0.0;
        for layer in &self.layers {
            for i in 0..=layer.intervals {
                let value = phi[layer.node_start + i];
                let coefficient = layer.coefficient_start + i;
                let weight = if i == 0 || i == layer.intervals {
                    0.5
                } else {
                    1.0
                };
                let mass = weight * layer.h * value * value / layer.density;
                norm += mass;
                volume_loss += mass * self.b1c[coefficient];
                slow += mass * (self.b1[coefficient] + 2.0)
                    / (self.omega * self.omega * layer.h * layer.h);
            }
        }
        if self.case.surface_boundary == SurfaceBoundary::FluidHalfSpace {
            let gamma = (x - self.surface_complex_k2.re).sqrt();
            let x1 = 0.999_999_9 * x;
            let x2 = 1.000_000_1 * x;
            let derivative = (self.surface_gamma(x2).re - self.surface_gamma(x1).re)
                / (self.case.surface_density_g_cm3 * (x2 - x1));
            norm += derivative * phi[0].powi(2);
            slow += phi[0].powi(2)
                / (2.0
                    * gamma
                    * self.case.surface_density_g_cm3
                    * self.case.surface_sound_speed_mps.powi(2));
        }
        if self.case.bottom_boundary.is_half_space() {
            let gamma = (x - self.bottom_complex_k2.re).sqrt();
            let x1 = 0.999_999_9 * x;
            let x2 = 1.000_000_1 * x;
            let derivative = if matches!(
                self.case.bottom_boundary,
                BottomBoundary::ElasticHalfSpace { .. }
            ) {
                (self.bottom_admittance(x2) - self.bottom_admittance(x1)) / (x2 - x1)
            } else {
                (self.bottom_gamma(x2).re - self.bottom_gamma(x1).re)
                    / (self.case.bottom_density_g_cm3 * (x2 - x1))
            };
            norm += derivative * phi[n - 1].powi(2);
            slow += phi[n - 1].powi(2)
                / (2.0
                    * gamma
                    * self.case.bottom_density_g_cm3
                    * self.case.bottom_sound_speed_mps.powi(2));
        }
        if norm <= 0.0 || !norm.is_finite() || !slow.is_finite() {
            return Err(error("KR0303", "invalid mode normalization", "modes"));
        }
        let turning = self
            .layers
            .iter()
            .flat_map(|layer| (1..=layer.intervals).map(move |i| (layer, i)))
            .find(|(layer, i)| {
                self.b1[layer.coefficient_start + i] - layer.h * layer.h * x + 2.0 > 0.0
            })
            .map_or(n - 2, |(layer, i)| layer.node_start + i);
        let scale = if phi[turning] < 0.0 {
            -norm.sqrt().recip()
        } else {
            norm.sqrt().recip()
        };
        for value in &mut phi {
            *value *= scale;
        }
        // Vector.f90 stores a single-precision mesh and interpolates into complex32 .mod samples.
        let grid_depths = crate::layers::grid(&self.layers);
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
        // Pinned top-half-space perturbation takes a default-kind CMPLX/SQRT.
        let top_loss = if self.case.surface_boundary == SurfaceBoundary::FluidHalfSpace {
            let z = Complex64::new(x, 0.0) - self.surface_complex_k2;
            -f64::from(
                num_complex::Complex32::new(z.re as f32, z.im as f32)
                    .sqrt()
                    .im,
            ) * phi[0].powi(2)
                / self.case.surface_density_g_cm3
        } else {
            0.0
        };
        let loss_k2 = volume_loss * scale * scale
            + top_loss
            + if self.case.bottom_boundary == BottomBoundary::FluidHalfSpace {
                -self.bottom_gamma(x).im * phi[n - 1].powi(2) / self.case.bottom_density_g_cm3
            } else {
                0.0
            };
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
    use super::solve;
    use crate::{Case, Interpolation, legacy::load_case, pekeris, profile::Profile};
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
        assert!((Profile::new(&n2).unwrap().speed(50.0) - expected).abs() < 1e-12);
        input.interpolation = Interpolation::CLinear;
        let linear = Case::from_definition(input.clone()).unwrap();
        assert!((Profile::new(&linear).unwrap().speed(50.0) - 1550.0).abs() < 1e-12);
        input.frequency_hz = 1000.0;
        input.mesh_points = 1_000_000;
        let report = solve(&Case::from_definition(input).unwrap()).unwrap_err();
        assert!(report.diagnostics()[0].message.contains("work limit"));
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
    fn analytic_auto_mesh_matches_pinned_fortran_rounding() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/MunkAnalytic");
        let mut input = load_case(root.with_extension("env"), root.with_extension("flp"))
            .unwrap()
            .into_definition();
        input.mesh_points = 0;
        input.max_range_m = 200_000.0;
        let modes = solve(&Case::from_definition(input).unwrap()).unwrap();
        assert_eq!(modes.modes.len(), 102);
        // Unmodified upstream MunkAnalytic.env: ninth .mod wavenumber, after extrapolation.
        assert_eq!(
            modes.modes[8].horizontal_wavenumber_rad_per_m.re as f32,
            f32::from_bits(0x3e55_599f)
        );
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
