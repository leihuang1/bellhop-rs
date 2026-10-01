//! Finite fluid-layer geometry, shared by both backends. Interfaces have two
//! material samples but one pressure unknown; no averaging of density or mesh step.
use crate::{
    Case, CaseDefinition, DiagnosticReport, MAX_MESH_POINTS, MAX_VECTOR_LENGTH, SoundSpeedPoint,
    error, profile::Profile, solver,
};

pub(crate) const MAX_LAYERS: usize = 500;

#[derive(Clone, Copy)]
pub(crate) struct Layer<'a> {
    pub top: f64,
    pub bottom: f64,
    pub density: f64,
    pub points: &'a [SoundSpeedPoint],
    pub loss: &'a [f64],
    pub mesh_points: usize,
}

impl CaseDefinition {
    /// Bottom of the entire finite fluid stack, excluding the half-space.
    #[must_use]
    pub fn total_depth_m(&self) -> f64 {
        self.additional_fluid_layers
            .last()
            .map_or(self.water_depth_m, |layer| layer.bottom_depth_m)
    }
}

pub(crate) fn iter(case: &CaseDefinition) -> impl Iterator<Item = Layer<'_>> {
    std::iter::once(Layer {
        top: 0.0,
        bottom: case.water_depth_m,
        density: case.water_density_g_cm3,
        points: &case.sound_speed_profile,
        loss: &case.water_attenuation_db_per_wavelength,
        mesh_points: case.mesh_points,
    })
    .chain(
        case.additional_fluid_layers
            .iter()
            .enumerate()
            .map(|(index, layer)| Layer {
                top: if index == 0 {
                    case.water_depth_m
                } else {
                    case.additional_fluid_layers[index - 1].bottom_depth_m
                },
                bottom: layer.bottom_depth_m,
                density: layer.density_g_cm3,
                points: &layer.sound_speed_profile,
                loss: &layer.attenuation_db_per_wavelength,
                mesh_points: layer.mesh_points,
            }),
    )
}

#[allow(clippy::float_cmp)]
pub(crate) fn validate(case: &CaseDefinition, diagnostics: &mut DiagnosticReport) {
    if case.additional_fluid_layers.len() >= MAX_LAYERS {
        diagnostics.push(error(
            "additional_fluid_layers",
            "at most 500 finite fluid layers",
        ));
    }
    let profile_values: usize = iter(case).map(|layer| layer.points.len()).sum();
    let loss_values: usize = iter(case).map(|layer| layer.loss.len()).sum();
    if profile_values > MAX_VECTOR_LENGTH || loss_values > MAX_VECTOR_LENGTH {
        diagnostics.push(error(
            "additional_fluid_layers",
            "total fluid profile storage exceeds the limit",
        ));
    }
    for (index, layer) in iter(case).enumerate().skip(1) {
        let field = format!("additional_fluid_layers[{}]", index - 1);
        if !layer.top.is_finite() || !layer.bottom.is_finite() || layer.bottom <= layer.top {
            diagnostics.push(error(
                format!("{field}.bottom_depth_m"),
                "layer interfaces must be finite and strictly increasing",
            ));
        }
        if !layer.density.is_finite() || layer.density <= 0.0 {
            diagnostics.push(error(
                format!("{field}.density_g_cm3"),
                "density must be finite and positive",
            ));
        }
        let invalid = !(2..=MAX_VECTOR_LENGTH).contains(&layer.points.len())
            || layer.points.iter().any(|p| {
                !p.depth_m.is_finite() || !p.sound_speed_mps.is_finite() || p.sound_speed_mps <= 0.0
            })
            || layer.points.first().is_none_or(|p| p.depth_m != layer.top)
            || layer
                .points
                .last()
                .is_none_or(|p| p.depth_m != layer.bottom)
            || layer
                .points
                .windows(2)
                .any(|p| p[1].depth_m <= p[0].depth_m);
        let bad_loss = !layer.loss.is_empty()
            && (layer.loss.len() != layer.points.len()
                || layer
                    .loss
                    .iter()
                    .any(|&a| !(0.0..=8.685_889_6 * 2.0 * std::f64::consts::PI).contains(&a)));
        if invalid {
            diagnostics.push(error(format!("{field}.sound_speed_profile"), "require increasing absolute depths from the previous interface to this layer bottom and positive sound speeds"));
        }
        if bad_loss {
            diagnostics.push(error(
                format!("{field}.attenuation_db_per_wavelength"),
                "require one finite nonnegative loss per node with Im(c) <= Re(c)",
            ));
        }
        if layer.mesh_points != 0 && !(10..=MAX_MESH_POINTS).contains(&layer.mesh_points) {
            diagnostics.push(error(
                format!("{field}.mesh_points"),
                "mesh points must be 0 or in 10..=1000000",
            ));
        }
        if !invalid
            && !bad_loss
            && let Err(report) = Profile::new_layer(case, layer)
        {
            let d = &report.diagnostics()[0];
            let suffix = if d.field == "water_attenuation_db_per_wavelength" {
                "attenuation_db_per_wavelength"
            } else {
                "sound_speed_profile"
            };
            diagnostics.push(error(format!("{field}.{suffix}"), &d.message));
        }
    }
}

pub(crate) fn minimum_speed(case: &CaseDefinition) -> Result<f64, DiagnosticReport> {
    let mut minimum = Profile::new(case)?.minimum_speed();
    for layer in iter(case).skip(1) {
        // Malformed additional layers have already been diagnosed before interpolation.
        if layer.points.len() >= 2
            && (layer.loss.is_empty() || layer.loss.len() == layer.points.len())
        {
            minimum = minimum.min(Profile::new_layer(case, layer)?.minimum_speed());
        }
    }
    Ok(minimum)
}

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
pub(crate) fn mesh_intervals(
    case: &Case,
    multiplier: usize,
) -> Result<Vec<usize>, DiagnosticReport> {
    let reference = case
        .mesh_reference_frequency_hz
        .unwrap_or(case.frequency_hz);
    let mut intervals = Vec::new();
    let mut total = 0;
    for (index, layer) in iter(case).enumerate() {
        let field = if index == 0 {
            "mesh_points".into()
        } else {
            format!("additional_fluid_layers[{}].mesh_points", index - 1)
        };
        let last_speed = layer.points.last().map_or(1500.0, |p| p.sound_speed_mps);
        let needed = ((layer.bottom - layer.top) / (last_speed / reference / 20.0))
            .floor()
            .max(10.0);
        let base = if layer.mesh_points == 0 {
            needed
        } else {
            layer.mesh_points as f64
        };
        if !needed.is_finite()
            || base < (needed as usize / 2) as f64
            || base > MAX_MESH_POINTS as f64
        {
            return Err(solver::error(
                "KR0302",
                "mesh is too coarse or exceeds the mesh limit",
                &field,
            ));
        }
        let scaled = if case.mesh_reference_frequency_hz.is_some() {
            (base * multiplier as f64 * case.frequency_hz / reference).floor()
        } else {
            base * multiplier as f64
        };
        if !scaled.is_finite() || !(10.0..=MAX_MESH_POINTS as f64).contains(&scaled) {
            return Err(solver::error(
                "KR0302",
                "scaled or refined mesh exceeds the mesh limits",
                &field,
            ));
        }
        total += scaled as usize;
        if total > MAX_MESH_POINTS {
            return Err(solver::error(
                "KR0302",
                "total fluid mesh exceeds the mesh limit",
                "mesh_points",
            ));
        }
        intervals.push(scaled as usize);
    }
    Ok(intervals)
}

pub(crate) struct MeshLayer {
    pub coefficient_start: usize,
    pub node_start: usize,
    pub intervals: usize,
    pub h: f64,
    pub density: f64,
}

#[allow(clippy::cast_precision_loss)]
pub(crate) fn mesh_layers(
    case: &Case,
    multiplier: usize,
) -> Result<Vec<MeshLayer>, DiagnosticReport> {
    let mut coefficient_start = 0;
    let mut node_start = 0;
    let mut layers = Vec::new();
    for (layer, intervals) in iter(case).zip(mesh_intervals(case, multiplier)?) {
        layers.push(MeshLayer {
            coefficient_start,
            node_start,
            intervals,
            h: (layer.bottom - layer.top) / intervals as f64,
            density: layer.density,
        });
        coefficient_start += intervals + 1;
        node_start += intervals;
    }
    Ok(layers)
}

#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub(crate) fn grid(layers: &[MeshLayer]) -> Vec<f32> {
    let mut grid = vec![0.0];
    for layer in layers {
        let top = grid[layer.node_start];
        grid.extend((1..=layer.intervals).map(|i| top + (i as f64 * layer.h) as f32));
    }
    grid
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_layer_scales_before_truncation_on_each_mesh() {
        let cases = crate::legacy::load_frequency_cases_from_sources(
            include_str!("../tests/fixtures/LayeredFluidPower.env"),
            include_str!("../tests/fixtures/LayeredFluidPower.flp"),
            std::path::Path::new("layers.env"),
            std::path::Path::new("layers.flp"),
            crate::ModeSolver::Krakenc,
        )
        .unwrap();
        for (index, expected) in [
            (0, [vec![151, 124], vec![303, 249], vec![606, 498]]),
            (1, [vec![101, 83], vec![202, 166], vec![404, 332]]),
            (2, [vec![126, 103], vec![252, 207], vec![505, 415]]),
            (3, [vec![101, 83], vec![202, 166], vec![404, 332]]),
        ] {
            for (multiplier, intervals) in [1, 2, 4].into_iter().zip(expected) {
                assert_eq!(
                    mesh_intervals(&cases[index], multiplier).unwrap(),
                    intervals
                );
            }
        }
    }
}
