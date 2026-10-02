// Adapted from Acoustics Toolbox v2023.5 Kraken/BCImpedanceMod.f90 and
// BCImpedancecMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE.
//! Elastic half-space impedance. Finite solid layers are not implemented here yet.
use crate::{Boundary, CaseDefinition, DiagnosticReport, Interpolation, ModeSolver, error};
use num_complex::Complex64;
use std::f64::consts::PI;

#[allow(clippy::too_many_lines)]
pub(crate) fn validate(case: &CaseDefinition, diagnostics: &mut DiagnosticReport) {
    if case.mode_solver == ModeSolver::Kraken
        && matches!(case.surface_boundary, Boundary::ElasticHalfSpace { .. })
    {
        diagnostics.push(error(
            "surface_boundary",
            "elastic top requires KRAKENC; KRAKEN top-elastic mode counting is not validated",
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
        if case.mode_solver == ModeSolver::Kraken && !case.additional_fluid_layers.is_empty() {
            diagnostics.push(error("additional_fluid_layers", "KRAKEN finite elasticity requires one contiguous fluid layer; multi-fluid elastic secant parity is not validated (use KRAKENC)"));
        }
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
    let cp = complex_speed(cp, if real { 0.0 } else { loss });
    let cs = complex_speed(*cs, if real { 0.0 } else { *shear_loss });
    let shear = x - omega * omega / cs.powi(2);
    let pressure = x - omega * omega / cp.powi(2);
    let gs = if real {
        Complex64::new(shear.sqrt().re, 0.0)
    } else {
        crate::complex_modes::pekeris_root(shear)
    };
    let gp = if real {
        Complex64::new(pressure.sqrt().re, 0.0)
    } else {
        crate::complex_modes::pekeris_root(pressure)
    };
    let mu = density * cs.powi(2);
    let f = omega * omega * gp * (x - shear);
    let g = ((shear + x).powi(2) - 4.0 * gs * gp * x) * mu;
    (f, g)
}

/// Constant-material compound-matrix mesh, outside the acoustic pressure unknowns.
pub(crate) struct SolidMesh {
    intervals: usize,
    h: f64,
    b: [Complex64; 4],
    rho: f64,
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

impl SolidMesh {
    #[allow(clippy::cast_precision_loss)]
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
            // Initialize uses Re(c²) for finite KRAKEN solids, unlike Re(c)² in its half-space.
            // The stiffness changes with loss, but no elastic absorption perturbation is added.
            let cp = if real { Complex64::new(cp.re, 0.0) } else { cp };
            let cs = if real { Complex64::new(cs.re, 0.0) } else { cs };
            let two_h = 2.0 * h;
            let rho = layer.density_g_cm3;
            let b = [
                two_h / (rho * cs),
                two_h / (rho * cp),
                4.0 * two_h * rho * cs * (cp - cs) / cp,
                two_h * (cp - 2.0 * cs) / cp,
            ];
            if b.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
                || !(two_h * omega * omega * rho).is_finite()
            {
                return Err(crate::solver::error(
                    "KR0302",
                    "non-finite elastic mesh coefficients",
                    "elastic_layers",
                ));
            }
            meshes.push(Self {
                intervals: n,
                h,
                b,
                rho: two_h * omega * omega * rho,
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

    fn step(&self, x: Complex64, y: [Complex64; 5]) -> [Complex64; 5] {
        let [b1, b2, b3, b4] = self.b;
        let xb3 = x * b3 - self.rho;
        [
            b1 * y[3] - b2 * y[4],
            -self.rho * y[3] - xb3 * y[4],
            2.0 * self.h * y[3] + b4 * y[4],
            xb3 * y[0] + b2 * y[1] - 2.0 * x * b4 * y[2],
            self.rho * y[0] - b1 * y[1] - 4.0 * self.h * x * y[2],
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
        let delta = self.step(x, y);
        let mut z = std::array::from_fn(|i| y[i] + sign * 0.5 * delta[i]);
        let mut previous = y;
        for step in 0..self.intervals {
            previous = y;
            y = z;
            let delta = self.step(x, y);
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
        std::array::from_fn(|i| (previous[i] + 2.0 * y[i] + z[i]) / 4.0)
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
            let cp = complex_speed(cp, if real { 0.0 } else { loss });
            let cs = complex_speed(*cs, if real { 0.0 } else { *shear_loss });
            let s = x - omega * omega / cs.powi(2);
            let p = x - omega * omega / cp.powi(2);
            let gs = if real {
                Complex64::new(s.sqrt().re, 0.0)
            } else {
                crate::complex_modes::pekeris_root(s)
            };
            let gp = if real {
                Complex64::new(p.sqrt().re, 0.0)
            } else {
                crate::complex_modes::pekeris_root(p)
            };
            let mu = density * cs.powi(2);
            [
                (gs * gp - x) / mu,
                ((s + x).powi(2) - 4.0 * gs * gp * x) * mu,
                2.0 * gs * gp - s - x,
                gp * (x - s),
                gs * (s - x),
            ]
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
