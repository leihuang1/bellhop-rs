// Adapted from Acoustics Toolbox v2023.5 Kraken/BCImpedanceMod.f90 and
// BCImpedancecMod.f90, Copyright (C) 2009 Michael B. Porter.
// GPL-3.0-or-later; see LICENSE.
//! Elastic half-space impedance. Finite solid layers are not implemented here yet.
use crate::{Boundary, CaseDefinition, DiagnosticReport, Interpolation, ModeSolver, error};
use num_complex::Complex64;
use std::f64::consts::PI;

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
