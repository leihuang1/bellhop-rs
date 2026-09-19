mod auxiliary;
mod env;
mod records;

use std::fs;
use std::path::Path;

use crate::case::validate_environment;
use crate::diagnostic::{Diagnostic, DiagnosticReport, LoadOutcome, SourceLocation};
use crate::model::{BoundaryCondition, EnvironmentCase};
use crate::{Case, CaseDefinition};

use auxiliary::{
    BoundarySide, parse_boundary_shape, parse_internal_reflection_coefficients,
    parse_range_dependent_sound_speed, parse_reflection_coefficients, parse_source_beam_pattern,
};

/// Loads and validates a two-dimensional legacy BELLHOP `.env` file.
///
/// This lower-level function parses only the primary environment. Use
/// [`load_case`] to resolve and validate its auxiliary files.
///
/// # Errors
///
/// Returns structured diagnostics when the file cannot be read, its legacy
/// records are malformed, it requests a three-dimensional option, or semantic
/// validation fails.
pub fn load_env(path: &Path) -> Result<LoadOutcome<EnvironmentCase>, DiagnosticReport> {
    env::parse(&read_environment(path)?, path)
}

/// Loads a complete two-dimensional legacy BELLHOP case.
///
/// Required auxiliary files are resolved beside the `.env` file using the
/// same stem. Independent auxiliary-file failures are collected into one
/// diagnostic report.
///
/// # Errors
///
/// Returns structured diagnostics when the environment or any required
/// `.ssp`, `.ati`, `.bty`, `.brc`, `.trc`, `.irc`, or `.sbp` input is missing,
/// malformed, inconsistent, or unsupported.
#[allow(clippy::too_many_lines)]
pub fn load_case(path: &Path) -> Result<LoadOutcome<Case>, DiagnosticReport> {
    let source = read_environment(path)?;
    let env::ParsedEnvironment {
        value: environment,
        mut diagnostics,
        mut locations,
    } = env::parse_unvalidated(&source, path)?;
    let mut preflight = DiagnosticReport::default();
    validate_environment(&environment, &locations, &mut preflight);
    if preflight.has_errors() {
        diagnostics.extend(preflight.diagnostics().iter().cloned());
        return Err(diagnostics);
    }

    let range_dependent_sound_speed = if environment
        .top_options
        .interpolation
        .needs_range_dependent_file()
    {
        let depths: Vec<f64> = environment
            .sound_speed
            .points
            .iter()
            .map(|point| point.depth_m)
            .collect();
        collect_auxiliary(
            read_auxiliary(path, "ssp", "range_dependent_sound_speed").and_then(
                |(source, path)| parse_range_dependent_sound_speed(&source, &path, &depths),
            ),
            &mut diagnostics,
        )
    } else {
        None
    };

    let altimetry = if environment.top_boundary.has_shape_file {
        collect_auxiliary(
            read_auxiliary(path, "ati", "altimetry").and_then(|(source, path)| {
                parse_boundary_shape(
                    &source,
                    &path,
                    BoundarySide::Top,
                    environment.sound_speed.top_depth_m,
                )
            }),
            &mut diagnostics,
        )
    } else {
        None
    };

    let bathymetry = if environment.bottom_boundary.has_shape_file {
        collect_auxiliary(
            read_auxiliary(path, "bty", "bathymetry").and_then(|(source, path)| {
                parse_boundary_shape(
                    &source,
                    &path,
                    BoundarySide::Bottom,
                    environment.sound_speed.bottom_depth_m,
                )
            }),
            &mut diagnostics,
        )
    } else {
        None
    };

    let bottom_reflection = if matches!(
        environment.bottom_boundary.condition,
        BoundaryCondition::ReflectionCoefficientFile
    ) {
        collect_auxiliary(
            read_auxiliary(path, "brc", "bottom_reflection").and_then(|(source, path)| {
                parse_reflection_coefficients(&source, &path, "bottom_reflection")
            }),
            &mut diagnostics,
        )
    } else {
        None
    };

    let top_reflection = if matches!(
        environment.top_boundary.condition,
        BoundaryCondition::ReflectionCoefficientFile
    ) {
        collect_auxiliary(
            read_auxiliary(path, "trc", "top_reflection").and_then(|(source, path)| {
                parse_reflection_coefficients(&source, &path, "top_reflection")
            }),
            &mut diagnostics,
        )
    } else {
        None
    };

    let internal_reflection = if matches!(
        environment.bottom_boundary.condition,
        BoundaryCondition::PrecalculatedReflectionCoefficient
    ) {
        let table = collect_auxiliary(
            read_auxiliary(path, "irc", "internal_reflection")
                .and_then(|(source, path)| parse_internal_reflection_coefficients(&source, &path)),
            &mut diagnostics,
        );
        if table.is_some() {
            locations.insert(
                "internal_reflection.frequency",
                SourceLocation::file(path.with_extension("irc")),
            );
        }
        table
    } else {
        None
    };

    let source_beam_pattern = if environment.run.has_source_beam_pattern {
        collect_auxiliary(
            read_auxiliary(path, "sbp", "source_beam_pattern")
                .and_then(|(source, path)| parse_source_beam_pattern(&source, &path)),
            &mut diagnostics,
        )
    } else {
        None
    };

    if diagnostics.has_errors() {
        diagnostics.extend(preflight.diagnostics().iter().cloned());
        return Err(diagnostics);
    }

    Case::from_definition_with(
        CaseDefinition {
            environment,
            range_dependent_sound_speed,
            altimetry,
            bathymetry,
            bottom_reflection,
            top_reflection,
            internal_reflection,
            source_beam_pattern,
        },
        &locations,
        diagnostics,
    )
}

fn read_environment(path: &Path) -> Result<String, DiagnosticReport> {
    if path.extension().and_then(|extension| extension.to_str()) != Some("env") {
        return Err(DiagnosticReport::from_diagnostic(Diagnostic::error(
            "BH0002",
            "input path must name a .env file",
            "input",
            SourceLocation::file(path),
        )));
    }
    fs::read_to_string(path)
        .map_err(|error| DiagnosticReport::from_diagnostic(Diagnostic::io(path, &error)))
}

fn read_auxiliary(
    environment_path: &Path,
    extension: &str,
    field: &'static str,
) -> Result<(String, std::path::PathBuf), Diagnostic> {
    let path = environment_path.with_extension(extension);
    fs::read_to_string(&path)
        .map(|source| (source, path.clone()))
        .map_err(|error| {
            Diagnostic::error(
                "BH0001",
                format!("unable to read required .{extension} input: {error}"),
                field,
                SourceLocation::file(path),
            )
        })
}

fn collect_auxiliary<T>(
    result: Result<T, Diagnostic>,
    diagnostics: &mut DiagnosticReport,
) -> Option<T> {
    match result {
        Ok(value) => Some(value),
        Err(diagnostic) => {
            diagnostics.push(diagnostic);
            None
        }
    }
}
