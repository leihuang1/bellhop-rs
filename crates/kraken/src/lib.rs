#![forbid(unsafe_code)]

use std::error::Error;
use std::fmt;
use std::ops::Deref;
use std::path::{Path, PathBuf};

use num_complex::Complex64;

pub mod legacy;
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

/// Unvalidated input for the initial homogeneous-fluid Pekeris slice.
#[derive(Clone, Debug, PartialEq)]
pub struct CaseDefinition {
    pub title: String,
    pub frequency_hz: f64,
    pub water_depth_m: f64,
    pub water_sound_speed_mps: f64,
    pub water_density_g_cm3: f64,
    pub bottom_sound_speed_mps: f64,
    pub bottom_density_g_cm3: f64,
    pub mesh_points: usize,
    pub c_low_mps: f64,
    pub c_high_mps: f64,
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
#[derive(Clone, Debug, PartialEq)]
pub struct Case(CaseDefinition);

impl Case {
    /// Validate a candidate Pekeris case.
    ///
    /// # Errors
    ///
    /// Returns every invalid field in the case definition.
    #[allow(clippy::too_many_lines)]
    pub fn from_definition(definition: CaseDefinition) -> Result<Self, DiagnosticReport> {
        let mut diagnostics = DiagnosticReport::default();

        for (field, value) in [
            ("frequency_hz", definition.frequency_hz),
            ("water_depth_m", definition.water_depth_m),
            ("water_sound_speed_mps", definition.water_sound_speed_mps),
            ("water_density_g_cm3", definition.water_density_g_cm3),
            ("bottom_sound_speed_mps", definition.bottom_sound_speed_mps),
            ("bottom_density_g_cm3", definition.bottom_density_g_cm3),
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
        if definition.water_depth_m <= 0.0 {
            diagnostics.push(error("water_depth_m", "water depth must be positive"));
        }
        if definition.water_sound_speed_mps <= 0.0
            || definition.bottom_sound_speed_mps <= definition.water_sound_speed_mps
        {
            diagnostics.push(error(
                "bottom_sound_speed_mps",
                "a trapped Pekeris waveguide requires a faster fluid bottom",
            ));
        }
        if definition.water_density_g_cm3 <= 0.0 || definition.bottom_density_g_cm3 <= 0.0 {
            diagnostics.push(error("density", "densities must be positive"));
        }
        if !(10..=MAX_MESH_POINTS).contains(&definition.mesh_points) {
            diagnostics.push(error(
                "mesh_points",
                format!("mesh points must be in 10..={MAX_MESH_POINTS}"),
            ));
        }
        if definition.c_low_mps <= 0.0 || definition.c_high_mps <= definition.c_low_mps {
            diagnostics.push(error("phase_speed_limits", "require 0 < c_low < c_high"));
        }
        if definition.max_range_m <= 0.0 {
            diagnostics.push(error("max_range_m", "maximum range must be positive"));
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
        }
        if definition
            .source_depths_m
            .iter()
            .chain(&definition.receiver_depths_m)
            .any(|depth| *depth < 0.0 || *depth > definition.water_depth_m)
        {
            diagnostics.push(error(
                "field_depths_m",
                "source and receiver depths must lie in water",
            ));
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
        if let (Some(first), Some(last)) = (
            definition.mode_sample_depths_m.first(),
            definition.mode_sample_depths_m.last(),
        ) && definition
            .source_depths_m
            .iter()
            .chain(&definition.receiver_depths_m)
            .any(|depth| depth < first || depth > last)
        {
            diagnostics.push(error(
                "mode_sample_depths_m",
                "FIELD depths must be covered by the legacy mode-sample depths",
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
    pub pressure: Vec<Complex64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResult {
    pub modes: ModeSet,
    pub field: PressureField,
}

/// Compute normal modes and the coherent range-independent line-source field.
///
/// # Errors
///
/// Returns a diagnostic if the phase-speed range contains no trapped modes or
/// the field exceeds the bounded modal-work limit.
pub fn solve(case: &Case) -> Result<SimulationResult, DiagnosticReport> {
    solver::solve(case)
}
