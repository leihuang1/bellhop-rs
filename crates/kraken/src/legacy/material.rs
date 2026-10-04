//! Parsed legacy materials own raw absorption and its power law until a solve
//! frequency is selected. Canonical Case fields never hold legacy loss units.
use crate::attenuation::{VolumeLoss, db_per_wavelength};
use crate::{Boundary, CaseDefinition, ElasticLayer, FluidLayer, SoundSpeedPoint};

type PowerLaw = (f64, f64, f64);

pub(super) struct Fluid {
    pub bottom_depth_m: f64,
    pub density_g_cm3: f64,
    pub mesh_points: usize,
    pub points: Vec<SoundSpeedPoint>,
    pub attenuation: Vec<f64>,
    pub power_law: PowerLaw,
}

pub(super) struct Elastic {
    pub bottom_depth_m: f64,
    pub compressional_speed: f64,
    pub shear_speed: f64,
    pub density_g_cm3: f64,
    pub compressional_attenuation: f64,
    pub shear_attenuation: f64,
    pub mesh_points: usize,
    pub power_law: PowerLaw,
}

pub(super) struct HalfSpace {
    pub values: [f64; 6],
    pub power_law: PowerLaw,
}

pub(super) struct Materials {
    pub water: Fluid,
    pub additional: Vec<Fluid>,
    pub top: Vec<Elastic>,
    pub bottom: Vec<Elastic>,
    pub surface_half_space: HalfSpace,
    pub bottom_half_space: HalfSpace,
    pub unit: u8,
    pub volume: VolumeLoss,
}

struct Loss<'a> {
    unit: u8,
    volume: &'a VolumeLoss,
    frequency: f64,
}

impl Loss<'_> {
    fn convert(&self, raw: f64, depth: f64, speed: f64, power: PowerLaw) -> f64 {
        db_per_wavelength(
            self.unit,
            self.volume,
            raw,
            depth,
            speed,
            self.frequency,
            power,
        )
    }
}

impl Fluid {
    fn canonical(&self, loss: &Loss<'_>) -> FluidLayer {
        let mut attenuation: Vec<_> = self
            .points
            .iter()
            .zip(&self.attenuation)
            .map(|(p, &a)| loss.convert(a, p.depth_m, p.sound_speed_mps, self.power_law))
            .collect();
        if attenuation.iter().all(|&a| a == 0.0) {
            attenuation.clear();
        }
        FluidLayer {
            bottom_depth_m: self.bottom_depth_m,
            density_g_cm3: self.density_g_cm3,
            mesh_points: self.mesh_points,
            sound_speed_profile: self.points.clone(),
            attenuation_db_per_wavelength: attenuation,
        }
    }
}

impl Elastic {
    fn canonical(&self, top: f64, loss: &Loss<'_>) -> ElasticLayer {
        ElasticLayer {
            bottom_depth_m: self.bottom_depth_m,
            compressional_sound_speed_mps: self.compressional_speed,
            shear_sound_speed_mps: self.shear_speed,
            density_g_cm3: self.density_g_cm3,
            compressional_attenuation_db_per_wavelength: loss.convert(
                self.compressional_attenuation,
                top,
                self.compressional_speed,
                self.power_law,
            ),
            shear_attenuation_db_per_wavelength: loss.convert(
                self.shear_attenuation,
                top,
                self.shear_speed,
                self.power_law,
            ),
            mesh_points: self.mesh_points,
        }
    }
}

impl HalfSpace {
    fn apply(&self, boundary: &mut Boundary, attenuation: &mut f64, loss: &Loss<'_>) {
        if boundary.is_half_space() {
            // UpdateHSLoss excludes depth-local biological layers using HUGE.
            *attenuation = loss.convert(self.values[4], f64::MAX, self.values[1], self.power_law);
        }
        if let Boundary::ElasticHalfSpace {
            shear_attenuation_db_per_wavelength,
            ..
        } = boundary
        {
            *shear_attenuation_db_per_wavelength =
                loss.convert(self.values[5], f64::MAX, self.values[2], self.power_law);
        }
    }
}

impl Materials {
    pub fn fluid_bottom(&self) -> f64 {
        self.additional
            .last()
            .map_or(self.water.bottom_depth_m, |layer| layer.bottom_depth_m)
    }

    pub fn fluid_top(&self) -> f64 {
        self.top.last().map_or(0.0, |layer| layer.bottom_depth_m)
    }

    pub fn input_values(&self) -> usize {
        7 * (self.top.len() + self.bottom.len())
            + self.water.points.len()
            + self.water.attenuation.len()
            + self
                .additional
                .iter()
                .map(|layer| 1 + layer.points.len() + layer.attenuation.len())
                .sum::<usize>()
    }

    pub fn case_definition(&self, template: &CaseDefinition, frequency: f64) -> CaseDefinition {
        let mut input = template.clone();
        input.frequency_hz = frequency;
        let loss = Loss {
            unit: self.unit,
            volume: &self.volume,
            frequency,
        };
        let water = self.water.canonical(&loss);
        input.water_depth_m = water.bottom_depth_m;
        input.water_density_g_cm3 = water.density_g_cm3;
        input.mesh_points = water.mesh_points;
        input.sound_speed_profile = water.sound_speed_profile;
        input.water_attenuation_db_per_wavelength = water.attenuation_db_per_wavelength;
        input.additional_fluid_layers = self
            .additional
            .iter()
            .map(|layer| layer.canonical(&loss))
            .collect();
        for (raw, layers, mut top) in [
            (&self.top, &mut input.top_elastic_layers, 0.0),
            (
                &self.bottom,
                &mut input.bottom_elastic_layers,
                self.fluid_bottom(),
            ),
        ] {
            layers.clear();
            for layer in raw {
                layers.push(layer.canonical(top, &loss));
                top = layer.bottom_depth_m;
            }
        }
        input.surface_sound_speed_mps = self.surface_half_space.values[1];
        input.surface_density_g_cm3 = self.surface_half_space.values[3];
        input.bottom_sound_speed_mps = self.bottom_half_space.values[1];
        input.bottom_density_g_cm3 = self.bottom_half_space.values[3];
        self.surface_half_space.apply(
            &mut input.surface_boundary,
            &mut input.surface_attenuation_db_per_wavelength,
            &loss,
        );
        self.bottom_half_space.apply(
            &mut input.bottom_boundary,
            &mut input.bottom_attenuation_db_per_wavelength,
            &loss,
        );
        input
    }
}
