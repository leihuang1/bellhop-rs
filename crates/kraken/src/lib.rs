#![forbid(unsafe_code)]

use std::error::Error;
use std::fmt;
use std::ops::Deref;
use std::path::{Path, PathBuf};

use num_complex::Complex64;

mod attenuation;
mod complex_modes;
pub mod legacy;
mod modes;
#[cfg(test)]
mod pekeris;
mod profile;
mod reflection;
mod solver;

const MAX_FIELD_SAMPLES: usize = 1_000_000;
const MAX_VECTOR_LENGTH: usize = 100_000;
const MAX_MESH_POINTS: usize = 1_000_000;
const MAX_MODE_LIMIT: usize = 20_000;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Diagnostic {
    pub code: &'static str,
    pub message: String,
    pub field: String,
    pub path: PathBuf,
    pub line: usize,
    pub column: usize,
}

impl Diagnostic {
    pub(crate) fn new(
        code: &'static str,
        message: impl Into<String>,
        field: impl Into<String>,
        path: impl Into<PathBuf>,
        line: usize,
        column: usize,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            field: field.into(),
            path: path.into(),
            line,
            column,
        }
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}:{}:{}: error[{}]: {} ({})",
            self.path.display(),
            self.line,
            self.column,
            self.code,
            self.message,
            self.field
        )
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct DiagnosticReport {
    diagnostics: Vec<Diagnostic>,
}

impl DiagnosticReport {
    #[must_use]
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    fn push(&mut self, diagnostic: Diagnostic) {
        self.diagnostics.push(diagnostic);
    }

    pub(crate) fn one(diagnostic: Diagnostic) -> Self {
        Self {
            diagnostics: vec![diagnostic],
        }
    }
}

impl fmt::Display for DiagnosticReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, diagnostic) in self.diagnostics.iter().enumerate() {
            if index > 0 {
                writeln!(formatter)?;
            }
            write!(formatter, "{diagnostic}")?;
        }
        Ok(())
    }
}

impl Error for DiagnosticReport {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Interpolation {
    N2Linear,
    CLinear,
    Pchip,
    Spline,
    /// The fixed 5000 m Munk profile from the pinned Fortran `misc/munk.f90`.
    AnalyticMunk,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ModeSolver {
    Kraken,
    Krakenc,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SourceGeometry {
    Line,
    Point,
}

/// Smooth fluid boundary; material values live in the corresponding case fields.
#[derive(Clone, Debug, PartialEq)]
pub enum Boundary {
    Vacuum,
    FluidHalfSpace,
    Rigid,
    /// KRAKENC F (top TRC or bottom BRC): magnitude and unwrapped phase versus grazing angle.
    Reflection(Vec<ReflectionPoint>),
    /// KRAKENC bottom P, scaled impedance functions versus squared wavenumber.
    Impedance {
        frequency_hz: f64,
        points: Vec<ImpedancePoint>,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReflectionPoint {
    pub angle_degrees: f64,
    pub magnitude: f64,
    pub phase_radians: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpedancePoint {
    pub wavenumber_squared: f64,
    pub f: Complex64,
    pub g: Complex64,
    pub power: i32,
}

/// Top V/R/A/F boundary. Top P is explicitly unsupported.
pub type SurfaceBoundary = Boundary;
/// Bottom V/R/A/F/P boundary.
pub type BottomBoundary = Boundary;

impl Boundary {
    pub(crate) fn is_tabulated(&self) -> bool {
        matches!(self, Self::Reflection(_) | Self::Impedance { .. })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SoundSpeedPoint {
    pub depth_m: f64,
    pub sound_speed_mps: f64,
}

/// Unvalidated input for a single fluid layer with smooth boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct CaseDefinition {
    pub title: String,
    pub mode_solver: ModeSolver,
    pub frequency_hz: f64,
    /// Nominal frequency for scaling the mesh; `None` uses `frequency_hz`.
    /// Legacy broadband inputs retain freq0 here, including for automatic meshes.
    pub mesh_reference_frequency_hz: Option<f64>,
    pub water_depth_m: f64,
    pub interpolation: Interpolation,
    pub surface_boundary: SurfaceBoundary,
    /// Zero unless the surface is a fluid half-space.
    pub surface_sound_speed_mps: f64,
    /// Zero unless the surface is a fluid half-space.
    pub surface_density_g_cm3: f64,
    /// Surface half-space loss at the solve frequency, in dB/wavelength.
    pub surface_attenuation_db_per_wavelength: f64,
    /// Empty for `AnalyticMunk`; otherwise depths from the surface to the interface.
    pub sound_speed_profile: Vec<SoundSpeedPoint>,
    pub water_density_g_cm3: f64,
    /// Empty for lossless water, otherwise one dB/wavelength value per SSP node.
    /// Values are evaluated at `frequency_hz` before complex SSP interpolation.
    pub water_attenuation_db_per_wavelength: Vec<f64>,
    pub bottom_boundary: BottomBoundary,
    /// Zero for a rigid or tabulated bottom (no half-space material).
    pub bottom_sound_speed_mps: f64,
    /// Zero for a rigid or tabulated bottom (no half-space material).
    pub bottom_density_g_cm3: f64,
    /// Bottom half-space attenuation in dB per wavelength (0 for lossless, rigid or table).
    pub bottom_attenuation_db_per_wavelength: f64,
    pub source_geometry: SourceGeometry,
    /// Mesh intervals at the reference frequency, or 0 for 20 per last-SSP wavelength.
    pub mesh_points: usize,
    pub c_low_mps: f64,
    pub c_high_mps: f64,
    /// Eigenvalue extrapolation convergence control; 0 selects the base mesh only.
    pub max_range_m: f64,
    /// Depths at which the legacy mode file samples each eigenfunction.
    pub mode_sample_depths_m: Vec<f64>,
    pub mode_limit: usize,
    pub source_depths_m: Vec<f64>,
    pub receiver_depths_m: Vec<f64>,
    pub receiver_ranges_m: Vec<f64>,
    pub receiver_offsets_m: Vec<f64>,
}

/// A validated case. Its definition is read-only after construction.
///
/// ```compile_fail
/// fn change_frequency(mut case: kraken::Case) {
///     case.frequency_hz = 100.0;
/// }
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct Case(CaseDefinition);

impl Case {
    /// Validate a candidate fluid-profile case.
    ///
    /// # Errors
    ///
    /// Returns every invalid field in the case definition.
    #[allow(clippy::too_many_lines, clippy::float_cmp)]
    pub fn from_definition(definition: CaseDefinition) -> Result<Self, DiagnosticReport> {
        let mut diagnostics = DiagnosticReport::default();

        for (field, value) in [
            ("frequency_hz", definition.frequency_hz),
            ("water_depth_m", definition.water_depth_m),
            ("water_density_g_cm3", definition.water_density_g_cm3),
            (
                "surface_sound_speed_mps",
                definition.surface_sound_speed_mps,
            ),
            ("surface_density_g_cm3", definition.surface_density_g_cm3),
            (
                "surface_attenuation_db_per_wavelength",
                definition.surface_attenuation_db_per_wavelength,
            ),
            ("bottom_sound_speed_mps", definition.bottom_sound_speed_mps),
            ("bottom_density_g_cm3", definition.bottom_density_g_cm3),
            (
                "bottom_attenuation_db_per_wavelength",
                definition.bottom_attenuation_db_per_wavelength,
            ),
            ("c_low_mps", definition.c_low_mps),
            ("c_high_mps", definition.c_high_mps),
            ("max_range_m", definition.max_range_m),
        ] {
            if !value.is_finite() {
                diagnostics.push(error(field, "value must be finite"));
            }
        }
        if definition.frequency_hz <= 0.0 {
            diagnostics.push(error("frequency_hz", "frequency must be positive"));
        }
        if definition
            .mesh_reference_frequency_hz
            .is_some_and(|frequency| !frequency.is_finite() || frequency <= 0.0)
        {
            diagnostics.push(error(
                "mesh_reference_frequency_hz",
                "mesh reference frequency must be finite and positive",
            ));
        }
        if definition.water_depth_m <= 0.0 {
            diagnostics.push(error("water_depth_m", "water depth must be positive"));
        }
        let profile = &definition.sound_speed_profile;
        let analytic = definition.interpolation == Interpolation::AnalyticMunk;
        let profile_invalid = if analytic {
            !profile.is_empty()
        } else {
            !(2..=MAX_VECTOR_LENGTH).contains(&profile.len())
                || profile.iter().any(|point| {
                    !point.depth_m.is_finite()
                        || !point.sound_speed_mps.is_finite()
                        || point.sound_speed_mps <= 0.0
                })
                || profile.first().is_none_or(|point| point.depth_m != 0.0)
                || profile
                    .last()
                    .is_none_or(|point| point.depth_m != definition.water_depth_m)
                || profile
                    .windows(2)
                    .any(|pair| pair[1].depth_m <= pair[0].depth_m)
        };
        if profile_invalid {
            diagnostics.push(error(
                "sound_speed_profile",
                if analytic {
                    "analytic Munk profile has no point records"
                } else {
                    "require finite increasing depths from 0 to water depth and positive sound speeds"
                },
            ));
        }
        let water_loss = &definition.water_attenuation_db_per_wavelength;
        let loss_invalid = !water_loss.is_empty()
            && (analytic
                || water_loss.len() != profile.len()
                || water_loss
                    .iter()
                    .any(|&a| !(0.0..=8.685_889_6 * 2.0 * std::f64::consts::PI).contains(&a)));
        if loss_invalid {
            diagnostics.push(error(
                "water_attenuation_db_per_wavelength",
                "require one finite nonnegative loss per SSP node, Im(c) <= Re(c); analytic Munk remains lossless",
            ));
        }
        if analytic {
            if definition.water_depth_m != 5000.0 {
                diagnostics.push(error(
                    "water_depth_m",
                    "analytic Munk profile requires 5000 m",
                ));
            }
            if definition.water_density_g_cm3 != 1.0 {
                diagnostics.push(error(
                    "water_density_g_cm3",
                    "analytic Munk profile requires density 1",
                ));
            }
        }
        if !profile_invalid && !loss_invalid {
            match profile::Profile::new(&definition) {
                Ok(profile) => {
                    if definition.mode_solver == ModeSolver::Kraken
                        && definition.surface_boundary == SurfaceBoundary::FluidHalfSpace
                        && definition.surface_sound_speed_mps > 0.0
                        && definition.surface_sound_speed_mps <= profile.minimum_speed()
                    {
                        diagnostics.push(error("surface_sound_speed_mps", "trapped modes require a surface half-space faster than the minimum water sound speed"));
                    }
                    if definition.mode_solver == ModeSolver::Kraken
                        && definition.bottom_boundary == BottomBoundary::FluidHalfSpace
                        && definition.bottom_sound_speed_mps > 0.0
                        && definition.bottom_sound_speed_mps <= profile.minimum_speed()
                    {
                        diagnostics.push(error(
                            "bottom_sound_speed_mps",
                            "trapped modes require a bottom faster than the minimum water sound speed",
                        ));
                    }
                }
                Err(report) => diagnostics.push(error(
                    &report.diagnostics()[0].field,
                    &report.diagnostics()[0].message,
                )),
            }
        }
        if let Err(message) = reflection::validate(&definition) {
            diagnostics.push(error("bottom_boundary", message));
        }
        if let Err(message) = reflection::validate_surface(&definition) {
            diagnostics.push(error("surface_boundary", message));
        }
        if definition.surface_boundary == SurfaceBoundary::FluidHalfSpace {
            if definition.surface_sound_speed_mps <= 0.0 || definition.surface_density_g_cm3 <= 0.0
            {
                diagnostics.push(error(
                    "surface_boundary",
                    "surface half-space requires positive sound speed and density",
                ));
            }
            if !(0.0..=8.685_889_6 * 2.0 * std::f64::consts::PI)
                .contains(&definition.surface_attenuation_db_per_wavelength)
            {
                diagnostics.push(error(
                    "surface_attenuation_db_per_wavelength",
                    "surface loss requires 0 <= Im(c) <= Re(c)",
                ));
            }
            if definition.mode_solver == ModeSolver::Kraken
                && definition.c_high_mps > definition.surface_sound_speed_mps
            {
                diagnostics.push(error(
                    "phase_speed_limits",
                    "KRAKEN does not support leaky modes above the surface half-space speed",
                ));
            }
        } else {
            for (field, value) in [
                (
                    "surface_sound_speed_mps",
                    definition.surface_sound_speed_mps,
                ),
                ("surface_density_g_cm3", definition.surface_density_g_cm3),
                (
                    "surface_attenuation_db_per_wavelength",
                    definition.surface_attenuation_db_per_wavelength,
                ),
            ] {
                if value != 0.0 {
                    diagnostics.push(error(
                        field,
                        "non-half-space surface has no half-space material",
                    ));
                }
            }
        }
        if definition.bottom_boundary != BottomBoundary::FluidHalfSpace {
            for (field, value) in [
                ("bottom_sound_speed_mps", definition.bottom_sound_speed_mps),
                ("bottom_density_g_cm3", definition.bottom_density_g_cm3),
                (
                    "bottom_attenuation_db_per_wavelength",
                    definition.bottom_attenuation_db_per_wavelength,
                ),
            ] {
                if value != 0.0 {
                    diagnostics.push(error(
                        field,
                        "non-half-space bottom has no half-space material",
                    ));
                }
            }
        } else if definition.bottom_attenuation_db_per_wavelength < 0.0
            || definition.bottom_attenuation_db_per_wavelength
                > 8.685_889_6 * 2.0 * std::f64::consts::PI
        {
            diagnostics.push(error(
                "bottom_attenuation_db_per_wavelength",
                "bottom attenuation must yield a complex speed with imaginary part no larger than real part",
            ));
        }
        if definition.mode_solver == ModeSolver::Kraken
            && definition.bottom_boundary == BottomBoundary::FluidHalfSpace
            && definition.c_high_mps > definition.bottom_sound_speed_mps
        {
            diagnostics.push(error(
                "phase_speed_limits",
                "leaky modes above the bottom sound speed are not supported yet",
            ));
        }
        if definition.water_density_g_cm3 <= 0.0 {
            diagnostics.push(error("water_density_g_cm3", "density must be positive"));
        }
        if definition.bottom_boundary == BottomBoundary::FluidHalfSpace {
            if definition.bottom_sound_speed_mps <= 0.0 {
                diagnostics.push(error(
                    "bottom_sound_speed_mps",
                    "fluid bottom sound speed must be positive",
                ));
            }
            if definition.bottom_density_g_cm3 <= 0.0 {
                diagnostics.push(error("bottom_density_g_cm3", "density must be positive"));
            }
        }
        if definition.mesh_points != 0 && !(10..=MAX_MESH_POINTS).contains(&definition.mesh_points)
        {
            diagnostics.push(error(
                "mesh_points",
                format!("mesh points must be 0 (automatic) or in 10..={MAX_MESH_POINTS}"),
            ));
        }
        if definition.c_low_mps <= 0.0 || definition.c_high_mps <= definition.c_low_mps {
            diagnostics.push(error("phase_speed_limits", "require 0 < c_low < c_high"));
        }
        if definition.max_range_m < 0.0 {
            diagnostics.push(error("max_range_m", "maximum range must be non-negative"));
        }
        if definition.mode_limit == 0 || definition.mode_limit > MAX_MODE_LIMIT {
            diagnostics.push(error(
                "mode_limit",
                format!("mode limit must be in 1..={MAX_MODE_LIMIT}"),
            ));
        }
        if definition.mode_sample_depths_m.is_empty()
            || definition.mode_sample_depths_m.len() > MAX_VECTOR_LENGTH
            || definition.mode_sample_depths_m.iter().any(|depth| {
                !depth.is_finite() || *depth < 0.0 || *depth > definition.water_depth_m
            })
            || definition
                .mode_sample_depths_m
                .windows(2)
                .any(|pair| pair[1] <= pair[0])
        {
            diagnostics.push(error(
                "mode_sample_depths_m",
                "mode depths must be finite, increasing, and inside the water column",
            ));
        }
        for (field, values) in [
            ("source_depths_m", &definition.source_depths_m),
            ("receiver_depths_m", &definition.receiver_depths_m),
            ("receiver_ranges_m", &definition.receiver_ranges_m),
        ] {
            if values.is_empty() || values.len() > MAX_VECTOR_LENGTH {
                diagnostics.push(error(field, "vector length is outside the supported range"));
            }
            if values.iter().any(|value| !value.is_finite()) {
                diagnostics.push(error(field, "values must be finite"));
            }
            if field != "receiver_ranges_m"
                && values
                    .iter()
                    .any(|depth| *depth < 0.0 || *depth > definition.water_depth_m)
            {
                diagnostics.push(error(field, "depths must lie in water"));
            }
        }
        if definition
            .receiver_ranges_m
            .iter()
            .any(|range| *range < 0.0)
            || definition
                .receiver_ranges_m
                .windows(2)
                .any(|pair| pair[1] <= pair[0])
        {
            diagnostics.push(error(
                "receiver_ranges_m",
                "receiver ranges must be non-negative and strictly increasing",
            ));
        }
        if definition.receiver_offsets_m.len() != definition.receiver_depths_m.len()
            || definition
                .receiver_offsets_m
                .iter()
                .any(|offset| !offset.is_finite())
        {
            diagnostics.push(error(
                "receiver_offsets_m",
                "receiver offsets must be finite and match receiver-depth count",
            ));
        }
        if definition.source_geometry == SourceGeometry::Point
            && definition.receiver_ranges_m.first().is_some_and(|&range| {
                definition
                    .receiver_offsets_m
                    .iter()
                    .any(|&offset| !(range + offset).is_finite() || range + offset < 0.0)
            })
        {
            diagnostics.push(error(
                "receiver_offsets_m",
                "point-source effective ranges must be finite and non-negative",
            ));
        }
        if let (Some(first), Some(last)) = (
            definition.mode_sample_depths_m.first(),
            definition.mode_sample_depths_m.last(),
        ) && definition
            .source_depths_m
            .iter()
            .chain(&definition.receiver_depths_m)
            .any(|depth| {
                let extension = if definition.mode_sample_depths_m.len() >= 2 {
                    1500.0 / definition.frequency_hz
                } else {
                    0.0
                };
                *depth < first - extension || *depth > last + extension
            })
        {
            diagnostics.push(error(
                "mode_sample_depths_m",
                "FIELD depths must lie within one reference wavelength of the mode-sample interval",
            ));
        }
        let field_samples = definition
            .source_depths_m
            .len()
            .checked_mul(definition.receiver_depths_m.len())
            .and_then(|count| count.checked_mul(definition.receiver_ranges_m.len()));
        if field_samples.is_none_or(|count| count > MAX_FIELD_SAMPLES) {
            diagnostics.push(error(
                "field_grid",
                format!("field grid may not exceed {MAX_FIELD_SAMPLES} samples"),
            ));
        }

        if diagnostics.diagnostics.is_empty() {
            Ok(Self(definition))
        } else {
            Err(diagnostics)
        }
    }

    #[must_use]
    pub fn into_definition(self) -> CaseDefinition {
        self.0
    }

    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub(crate) fn mesh_points_at(&self, multiplier: usize) -> Result<usize, DiagnosticReport> {
        let reference = self
            .mesh_reference_frequency_hz
            .unwrap_or(self.frequency_hz);
        // Analytic INIT leaves the reference's last-read sound speed at 1500.
        let last_speed = self
            .sound_speed_profile
            .last()
            .map_or(1500.0, |p| p.sound_speed_mps);
        let needed = (self.water_depth_m / (last_speed / reference / 20.0))
            .floor()
            .max(10.0);
        let base = if self.mesh_points == 0 {
            needed
        } else {
            self.mesh_points as f64
        };
        if !needed.is_finite()
            || base < (needed as usize / 2) as f64
            || base > MAX_MESH_POINTS as f64
        {
            return Err(solver::error(
                "KR0302",
                "mesh is too coarse or exceeds the mesh limit",
                "mesh_points",
            ));
        }
        // INT(NG * NV * freq / freq0), not INT(NG * freq / freq0) * NV.
        let scaled = if self.mesh_reference_frequency_hz.is_some() {
            (base * multiplier as f64 * self.frequency_hz / reference).floor()
        } else {
            base * multiplier as f64
        };
        if !scaled.is_finite() || !(10.0..=MAX_MESH_POINTS as f64).contains(&scaled) {
            return Err(solver::error(
                "KR0302",
                "scaled or refined mesh exceeds the mesh limits",
                "mesh_points",
            ));
        }
        Ok(scaled as usize)
    }
}

fn error(field: impl Into<String>, message: impl Into<String>) -> Diagnostic {
    Diagnostic::new("KR0201", message, field, Path::new("<case>"), 1, 1)
}

impl Deref for Case {
    type Target = CaseDefinition;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct NormalMode {
    pub horizontal_wavenumber_rad_per_m: Complex64,
    pub phase_speed_mps: f64,
    pub group_speed_mps: f64,
    pub attenuation_nepers_per_m: f64,
    /// Complex pressure eigenfunction values, ordered like `ModeSet::sampled_depths_m`.
    pub eigenfunction: Vec<Complex64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModeSet {
    pub frequency_hz: f64,
    pub sampled_depths_m: Vec<f64>,
    pub modes: Vec<NormalMode>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PressureField {
    pub source_depths_m: Vec<f64>,
    pub receiver_depths_m: Vec<f64>,
    pub receiver_ranges_m: Vec<f64>,
    pub receiver_offsets_m: Vec<f64>,
    /// Row-major values indexed as `[source_depth][receiver_depth][receiver_range]`.
    /// FIELD uses reference single-precision modal products and accumulation;
    /// final values are promoted to `Complex64` for the result model.
    pub pressure: Vec<Complex64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResult {
    pub modes: ModeSet,
    pub field: PressureField,
}

/// Compute normal modes and the coherent range-independent line- or point-source field.
///
/// # Errors
///
/// Returns a diagnostic if the phase-speed range contains no supported modes or
/// the field exceeds the bounded modal-work limit.
pub fn solve(case: &Case) -> Result<SimulationResult, DiagnosticReport> {
    solver::solve(case)
}

/// Compute complex normal modes without synthesizing a FIELD.
///
/// # Errors
///
/// Returns a diagnostic when the fluid complex-mode subset cannot be solved.
pub fn solve_complex_modes(case: &Case) -> Result<ModeSet, DiagnosticReport> {
    if case.mode_solver != ModeSolver::Krakenc {
        return Err(DiagnosticReport::one(error(
            "mode_solver",
            "requires KRAKENC input",
        )));
    }
    complex_modes::solve(case)
}
