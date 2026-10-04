// Adapted from Acoustics Toolbox v2023.5 Kraken/BCImpedanceMod.f90 and
// BCImpedancecMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE.
//! Elastic half-space impedance and finite homogeneous solid-cap transfer.
use crate::{Boundary, CaseDefinition, DiagnosticReport, Interpolation, ModeSolver, error};
use num_complex::Complex64;
use std::f64::consts::PI;

#[allow(clippy::too_many_lines)]
pub(crate) fn validate(case: &CaseDefinition, diagnostics: &mut DiagnosticReport) {
    if case.mode_solver == ModeSolver::Kraken
        && has_layers(case)
        && matches!(case.surface_boundary, Boundary::ElasticHalfSpace { .. })
    {
        diagnostics.push(error(
            "surface_boundary",
            "KRAKEN elastic top half-space with finite solids is not validated",
        ));
    }
    for (boundary, cp, field) in [
        (
            &case.surface_boundary,
            case.surface_sound_speed_mps,
            "surface_boundary",
        ),
        (
            &case.bottom_boundary,
            case.bottom_sound_speed_mps,
            "bottom_boundary",
        ),
    ] {
        if let Boundary::ElasticHalfSpace {
            shear_sound_speed_mps: cs,
            shear_attenuation_db_per_wavelength: loss,
        } = boundary
        {
            // Positive bulk modulus: cp² > 4/3 cs², not merely cp > cs.
            if !cs.is_finite() || *cs <= 0.0 || cp <= (4.0_f64 / 3.0).sqrt() * cs {
                diagnostics.push(error(
                    field,
                    "elastic half-space requires positive shear speed and cp² > 4/3 cs²",
                ));
            }
            if !(0.0..=8.685_889_6 * 2.0 * PI).contains(loss) {
                diagnostics.push(error(
                    field,
                    "shear loss requires finite 0 <= Im(cs) <= Re(cs)",
                ));
            }
            if case.interpolation == Interpolation::AnalyticMunk
                || case.surface_boundary.is_tabulated()
                || case.bottom_boundary.is_tabulated()
            {
                diagnostics.push(error(
                    field,
                    "elastic half-spaces require N/C/P/S and smooth V/R/A boundaries",
                ));
            }
        }
    }
    if has_layers(case) {
        if case.interpolation == Interpolation::AnalyticMunk
            || case.surface_boundary.is_tabulated()
            || case.bottom_boundary.is_tabulated()
        {
            diagnostics.push(error(
                "top_elastic_layers",
                "finite elastic layers require N/C/P/S and smooth V/R/A boundaries",
            ));
        }
        for (name, layers, start, boundary) in [
            (
                "top_elastic_layers",
                &case.top_elastic_layers,
                0.0,
                &case.surface_boundary,
            ),
            (
                "bottom_elastic_layers",
                &case.bottom_elastic_layers,
                case.fluid_bottom_depth_m(),
                &case.bottom_boundary,
            ),
        ] {
            if !layers.is_empty() && *boundary == Boundary::FluidHalfSpace {
                diagnostics.push(error(name, "finite elastic layers require vacuum, rigid or elastic outer boundary; pinned acoustic half-space compound state is undefined"));
            }
            let mut top = start;
            for (index, layer) in layers.iter().enumerate() {
                let field = format!("{name}[{index}]");
                if !layer.bottom_depth_m.is_finite() || layer.bottom_depth_m <= top {
                    diagnostics.push(error(
                        format!("{field}.bottom_depth_m"),
                        "elastic interfaces must be finite and strictly increasing",
                    ));
                }
                let cp = layer.compressional_sound_speed_mps;
                let cs = layer.shear_sound_speed_mps;
                if !cp.is_finite()
                    || !cs.is_finite()
                    || cs <= 0.0
                    || cp <= (4.0_f64 / 3.0).sqrt() * cs
                    || !layer.density_g_cm3.is_finite()
                    || layer.density_g_cm3 <= 0.0
                {
                    diagnostics.push(error(
                        format!("{field}.material"),
                        "require finite positive density and shear speed, cp² > 4/3 cs²",
                    ));
                }
                for loss in [
                    layer.compressional_attenuation_db_per_wavelength,
                    layer.shear_attenuation_db_per_wavelength,
                ] {
                    if !(0.0..=8.685_889_6 * 2.0 * PI).contains(&loss) {
                        diagnostics.push(error(
                            format!("{field}.attenuation"),
                            "elastic loss requires finite 0 <= Im(c) <= Re(c)",
                        ));
                    }
                }
                if layer.mesh_points != 0
                    && !(10..=crate::MAX_MESH_POINTS).contains(&layer.mesh_points)
                {
                    diagnostics.push(error(
                        format!("{field}.mesh_points"),
                        "mesh points must be 0 or in 10..=1000000",
                    ));
                }
                top = layer.bottom_depth_m;
            }
        }
        if case.fluid_top_depth_m() >= case.water_depth_m {
            diagnostics.push(error(
                "water_depth_m",
                "first fluid bottom must lie below the elastic cap",
            ));
        }
    }
}

pub(crate) fn has_layers(case: &CaseDefinition) -> bool {
    !case.top_elastic_layers.is_empty() || !case.bottom_elastic_layers.is_empty()
}

pub(crate) fn has_half_space(case: &CaseDefinition) -> bool {
    [&case.surface_boundary, &case.bottom_boundary]
        .iter()
        .any(|boundary| matches!(boundary, Boundary::ElasticHalfSpace { .. }))
}

pub(crate) fn minimum_speed(case: &CaseDefinition, water: f64) -> f64 {
    let mut minimum = water;
    for (boundary, cp) in [
        (&case.surface_boundary, case.surface_sound_speed_mps),
        (&case.bottom_boundary, case.bottom_sound_speed_mps),
    ] {
        match boundary {
            Boundary::ElasticHalfSpace {
                shear_sound_speed_mps,
                ..
            } => {
                minimum = minimum.min(*shear_sound_speed_mps);
            }
            Boundary::FluidHalfSpace => minimum = minimum.min(cp),
            _ => {}
        }
    }
    for layer in case
        .top_elastic_layers
        .iter()
        .chain(&case.bottom_elastic_layers)
    {
        minimum = minimum.min(layer.shear_sound_speed_mps);
    }
    // Initialize reduces cMin for Scholte waves; 0.85 is a default-kind literal.
    f64::from(0.85_f32) * minimum
}

pub(crate) fn maximum_speed(case: &CaseDefinition) -> f64 {
    [&case.surface_boundary, &case.bottom_boundary].iter().fold(
        case.c_high_mps,
        |maximum, boundary| match boundary {
            Boundary::ElasticHalfSpace {
                shear_sound_speed_mps,
                ..
            } => maximum.min(*shear_sound_speed_mps),
            _ => maximum,
        },
    )
}

pub(crate) fn complex_speed(speed: f64, loss: f64) -> Complex64 {
    Complex64::new(speed, loss * speed / (8.685_889_6 * 2.0 * PI))
}

pub(crate) fn half_space(
    boundary: &Boundary,
    x: Complex64,
    omega: f64,
    cp: f64,
    density: f64,
    loss: f64,
    real: bool,
) -> (Complex64, Complex64) {
    let Boundary::ElasticHalfSpace {
        shear_sound_speed_mps: cs,
        shear_attenuation_db_per_wavelength: shear_loss,
    } = boundary
    else {
        unreachable!()
    };
    // BCImpedanceMod ignores half-space elastic attenuation even with ComplexFlag.
    if real {
        let y = real_half_space(x.re, omega, cp, *cs, density);
        return (Complex64::from(omega * omega * y[3]), Complex64::from(y[1]));
    }
    let cp = complex_speed(cp, loss);
    let cs = complex_speed(*cs, *shear_loss);
    let shear = x - omega * omega / cs.powi(2);
    let pressure = x - omega * omega / cp.powi(2);
    let gs = crate::complex_modes::pekeris_root(shear);
    let gp = crate::complex_modes::pekeris_root(pressure);
    let mu = density * cs.powi(2);
    let f = omega * omega * gp * (x - shear);
    let g = ((shear + x).powi(2) - 4.0 * gs * gp * x) * mu;
    (f, g)
}

// BCImpedanceMod's real branch, including pinned GNU grouping: x - (x - q)
// simplifies to q, rather than introducing another subtraction near a pole.
fn real_half_space(x: f64, omega: f64, cp: f64, cs: f64, density: f64) -> [f64; 5] {
    let shear_k2 = omega * omega / (cs * cs);
    let s = x - shear_k2;
    let p = x - omega * omega / (cp * cp);
    let gs = if s < 0.0 { 0.0 } else { s.sqrt() };
    let gp = if p < 0.0 { 0.0 } else { p.sqrt() };
    let product = gs * gp;
    let mu = density * (cs * cs);
    [
        (product - x) / mu,
        mu * ((x + s).powi(2) - (4.0 * x) * product),
        (shear_k2 - 2.0 * x) + 2.0 * product,
        gp * shear_k2,
        gs * -shear_k2,
    ]
}

/// Constant-material compound-matrix mesh, outside the acoustic pressure unknowns.
pub(crate) struct SolidMesh {
    intervals: usize,
    h: f64,
    coefficients: Vec<([Complex64; 4], f64)>,
    real: bool,
}

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn intervals(
    case: &crate::Case,
    layer: &crate::ElasticLayer,
    top: f64,
    multiplier: usize,
) -> Result<usize, DiagnosticReport> {
    let reference = case
        .mesh_reference_frequency_hz
        .unwrap_or(case.frequency_hz);
    let needed = ((layer.bottom_depth_m - top) / (layer.shear_sound_speed_mps / reference / 20.0))
        .floor()
        .max(10.0);
    let base = if layer.mesh_points == 0 {
        needed
    } else {
        layer.mesh_points as f64
    };
    let scaled = if case.mesh_reference_frequency_hz.is_some() {
        (base * multiplier as f64 * case.frequency_hz / reference).floor()
    } else {
        base * multiplier as f64
    };
    if !needed.is_finite()
        || base < (needed as usize / 2) as f64
        || !scaled.is_finite()
        || !(10.0..=crate::MAX_MESH_POINTS as f64).contains(&scaled)
    {
        return Err(crate::solver::error(
            "KR0302",
            "elastic mesh is too coarse or exceeds the mesh limit",
            "elastic_layers.mesh_points",
        ));
    }
    Ok(scaled as usize)
}

pub(crate) fn mesh_points(
    case: &crate::Case,
    multiplier: usize,
) -> Result<usize, DiagnosticReport> {
    let mut total = 0;
    for (layers, mut top) in [
        (&case.top_elastic_layers, 0.0),
        (&case.bottom_elastic_layers, case.fluid_bottom_depth_m()),
    ] {
        for layer in layers {
            total += intervals(case, layer, top, multiplier)?;
            top = layer.bottom_depth_m;
            if total > crate::MAX_MESH_POINTS {
                return Err(crate::solver::error(
                    "KR0302",
                    "total elastic mesh exceeds the mesh limit",
                    "mesh_points",
                ));
            }
        }
    }
    Ok(total)
}

fn sampled_speed(
    case: &crate::Case,
    layer: &crate::ElasticLayer,
    top: f64,
    n: usize,
    speed: f64,
    attenuation: f64,
) -> Result<Vec<Complex64>, crate::DiagnosticReport> {
    let points = [top, layer.bottom_depth_m].map(|depth_m| crate::SoundSpeedPoint {
        depth_m,
        sound_speed_mps: speed,
    });
    let loss = [attenuation; 2];
    let profile = crate::profile::Profile::new_layer(
        case,
        crate::layers::Layer {
            top,
            bottom: layer.bottom_depth_m,
            density: layer.density_g_cm3,
            points: &points,
            loss: &loss,
            mesh_points: n,
        },
    )?;
    Ok((0..=n).map(|i| profile.mesh_complex_speed(i, n)).collect())
}

impl SolidMesh {
    #[allow(clippy::cast_precision_loss, clippy::too_many_lines)] // Keep pinned coefficient grouping together.
    pub(crate) fn build(
        case: &crate::Case,
        multiplier: usize,
        top_side: bool,
        real: bool,
    ) -> Result<Vec<Self>, DiagnosticReport> {
        let layers = if top_side {
            &case.top_elastic_layers
        } else {
            &case.bottom_elastic_layers
        };
        let mut top = if top_side {
            0.0
        } else {
            case.fluid_bottom_depth_m()
        };
        let mut meshes = Vec::new();
        let omega = 2.0 * PI * case.frequency_hz;
        for layer in layers {
            let n = intervals(case, layer, top, multiplier)?;
            let h = (layer.bottom_depth_m - top) / n as f64;
            let cp = complex_speed(
                layer.compressional_sound_speed_mps,
                layer.compressional_attenuation_db_per_wavelength,
            )
            .powi(2);
            let cs = complex_speed(
                layer.shear_sound_speed_mps,
                layer.shear_attenuation_db_per_wavelength,
            )
            .powi(2);
            let two_h = 2.0 * h;
            let rho = layer.density_g_cm3;
            let coefficients = if real {
                // Homogeneous N/C material still has node-specific SSP rounding.
                // Initialize uses Re(c²), not the half-space's Re(c)².
                let cp = sampled_speed(
                    case,
                    layer,
                    top,
                    n,
                    layer.compressional_sound_speed_mps,
                    layer.compressional_attenuation_db_per_wavelength,
                )?;
                let cs = sampled_speed(
                    case,
                    layer,
                    top,
                    n,
                    layer.shear_sound_speed_mps,
                    layer.shear_attenuation_db_per_wavelength,
                )?;
                (0..=n)
                    .map(|i| {
                        let weight = if i == n {
                            1.0
                        } else {
                            (top + i as f64 * h - top) / (layer.bottom_depth_m - top)
                        };
                        let density = if matches!(
                            case.interpolation,
                            crate::Interpolation::N2Linear | crate::Interpolation::CLinear
                        ) {
                            (1.0 - weight) * rho + weight * rho
                        } else {
                            rho
                        };
                        let pressure_squared = cp[i].powi(2).re;
                        let cs_real2 = cs[i].re * cs[i].re;
                        let cs_imag2 = cs[i].im * cs[i].im;
                        let shear_squared = cs_real2 - cs_imag2;
                        // Preserve Initialize's pinned reassociation of Re(cP²-cS²).
                        let difference = pressure_squared - cs_real2 + cs_imag2;
                        let b = [
                            two_h / (density * shear_squared),
                            two_h / (density * pressure_squared),
                            (shear_squared * difference) * (4.0 * density) * two_h
                                / pressure_squared,
                            two_h * (pressure_squared - 2.0 * shear_squared) / pressure_squared,
                        ]
                        .map(Complex64::from);
                        (b, (omega * omega * density) * two_h)
                    })
                    .collect::<Vec<_>>()
            } else {
                vec![(
                    [
                        two_h / (rho * cs),
                        two_h / (rho * cp),
                        4.0 * two_h * rho * cs * (cp - cs) / cp,
                        two_h * (cp - 2.0 * cs) / cp,
                    ],
                    two_h * omega * omega * rho,
                )]
            };
            if coefficients.iter().any(|(b, rho)| {
                !rho.is_finite() || b.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
            }) {
                return Err(crate::solver::error(
                    "KR0302",
                    "non-finite elastic mesh coefficients",
                    "elastic_layers",
                ));
            }
            meshes.push(Self {
                intervals: n,
                h,
                coefficients,
                real,
            });
            top = layer.bottom_depth_m;
        }
        Ok(meshes)
    }

    pub(crate) fn points(&self) -> usize {
        self.intervals + 1
    }

    pub(crate) fn spacing(&self) -> f64 {
        self.h
    }

    fn step(&self, x: Complex64, y: [Complex64; 5], node: usize, euler: bool) -> [Complex64; 5] {
        let (b, rho) = self.coefficients[node.min(self.coefficients.len() - 1)];
        if self.real {
            // GNU's real ElasticUP/DN group Euler and midpoint steps differently.
            // Keep KRAKENC's separate complex arithmetic unchanged.
            let [b1, b2, b3, b4] = b.map(|b| b.re);
            let y = y.map(|y| y.re);
            let two_x = 2.0 * x.re;
            let four_h_x = self.h * (4.0 * x.re);
            let xb3 = x.re * b3 - rho;
            let fourth = if euler {
                xb3 * y[0] + (b2 * y[1] - b4 * (two_x * y[2]))
            } else {
                b2 * y[1] + (xb3 * y[0] - two_x * (b4 * y[2]))
            };
            let fifth = if euler {
                (rho * y[0] - four_h_x * y[2]) - b1 * y[1]
            } else {
                (rho * y[0] - b1 * y[1]) - four_h_x * y[2]
            };
            return [
                b1 * y[3] - b2 * y[4],
                -(rho * y[3] + xb3 * y[4]),
                2.0 * self.h * y[3] + b4 * y[4],
                fourth,
                fifth,
            ]
            .map(Complex64::from);
        }
        let [b1, b2, b3, b4] = b;
        let xb3 = x * b3 - rho;
        [
            b1 * y[3] - b2 * y[4],
            -rho * y[3] - xb3 * y[4],
            2.0 * self.h * y[3] + b4 * y[4],
            xb3 * y[0] + b2 * y[1] - 2.0 * x * b4 * y[2],
            rho * y[0] - b1 * y[1] - 4.0 * self.h * x * y[2],
        ]
    }

    fn transfer(
        &self,
        x: Complex64,
        mut y: [Complex64; 5],
        top: bool,
        power: &mut i32,
    ) -> [Complex64; 5] {
        let sign = if top { 1.0 } else { -1.0 };
        let delta = self.step(x, y, if top { 0 } else { self.intervals }, true);
        let mut z = std::array::from_fn(|i| y[i] + sign * 0.5 * delta[i]);
        let mut previous = y;
        for step in 0..self.intervals {
            previous = y;
            y = z;
            let node = if top {
                step + 1
            } else {
                self.intervals - step - 1
            };
            let delta = self.step(x, y, node, false);
            z = std::array::from_fn(|i| previous[i] + sign * delta[i]);
            if step + 1 != self.intervals {
                let scale = if z[1].re.abs() < 1e-50 {
                    *power -= 50;
                    1e50
                } else if z[1].re.abs() > 1e50 {
                    *power += 50;
                    1e-50
                } else {
                    1.0
                };
                for i in 0..5 {
                    y[i] *= scale;
                    z[i] *= scale;
                }
            }
        }
        std::array::from_fn(|i| {
            if self.real {
                ((previous[i] + z[i]) + 2.0 * y[i]) * 0.25
            } else {
                (previous[i] + 2.0 * y[i] + z[i]) / 4.0
            }
        })
    }
}

/// Returns the effective Robin data after propagating through a finite solid cap.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cap_impedance(
    boundary: &Boundary,
    x: Complex64,
    omega: f64,
    cp: f64,
    density: f64,
    loss: f64,
    meshes: &[SolidMesh],
    top: bool,
    real: bool,
) -> (Complex64, Complex64, i32) {
    let mut y = match boundary {
        Boundary::Vacuum => [1.0.into(), 0.0.into(), 0.0.into(), 0.0.into(), 0.0.into()],
        Boundary::Rigid => [0.0.into(), 1.0.into(), 0.0.into(), 0.0.into(), 0.0.into()],
        Boundary::ElasticHalfSpace {
            shear_sound_speed_mps: cs,
            shear_attenuation_db_per_wavelength: shear_loss,
        } => {
            if real {
                real_half_space(x.re, omega, cp, *cs, density).map(Complex64::from)
            } else {
                let cp = complex_speed(cp, loss);
                let cs = complex_speed(*cs, *shear_loss);
                let s = x - omega * omega / cs.powi(2);
                let p = x - omega * omega / cp.powi(2);
                let gs = crate::complex_modes::pekeris_root(s);
                let gp = crate::complex_modes::pekeris_root(p);
                let mu = density * cs.powi(2);
                [
                    (gs * gp - x) / mu,
                    ((s + x).powi(2) - 4.0 * gs * gp * x) * mu,
                    2.0 * gs * gp - s - x,
                    gp * (x - s),
                    gs * (s - x),
                ]
            }
        }
        _ => unreachable!("cap outer boundary was validated"),
    };
    let mut power = 0;
    if top {
        for mesh in meshes {
            y = mesh.transfer(x, y, true, &mut power);
        }
    } else {
        for mesh in meshes.iter().rev() {
            y = mesh.transfer(x, y, false, &mut power);
        }
    }
    (omega * omega * y[3], y[1], power)
}
