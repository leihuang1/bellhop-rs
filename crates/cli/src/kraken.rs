#![forbid(unsafe_code)]

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Subcommand, ValueEnum};
use output::kraken::{
    DEFAULT_MAX_OUTPUT_BYTES, RunError, run_json_directory, run_legacy_directory,
};

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Run legacy .env/.flp or self-contained .json input and safely update a result directory.
    Run {
        case: PathBuf,
        /// Legacy FIELD geometry; defaults to the same-stem .flp. Invalid for JSON.
        #[arg(long)]
        flp: Option<PathBuf>,
        /// Legacy engine (default kraken); for JSON, assert the document's engine.
        #[arg(long, value_enum)]
        solver: Option<Solver>,
        /// Result directory, default `<case-stem>` in the current directory. Only owned files are replaced.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Numerical output format; HTTP is unaffected.
        #[arg(long, value_enum, default_value = "legacy")]
        format: crate::Format,
        /// Cumulative file/payload quota, checked before datasets and after flushes.
        #[arg(long, default_value_t = DEFAULT_MAX_OUTPUT_BYTES, value_parser = clap::value_parser!(u64).range(1..))]
        max_output_bytes: u64,
    },
    /// Export validated legacy inputs and all resources to self-contained JSON on stdout.
    Export {
        case: PathBuf,
        #[arg(long)]
        flp: Option<PathBuf>,
        #[arg(long, value_enum)]
        solver: Option<Solver>,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Solver {
    Kraken,
    Krakenc,
}

impl From<Solver> for kraken::ModeSolver {
    fn from(value: Solver) -> Self {
        match value {
            Solver::Kraken => Self::Kraken,
            Solver::Krakenc => Self::Krakenc,
        }
    }
}

pub fn execute(command: Command) -> ExitCode {
    let result = match command {
        Command::Run {
            case,
            flp,
            solver,
            output,
            format,
            max_output_bytes,
        } => {
            let output = output.unwrap_or_else(|| {
                PathBuf::from(
                    case.file_stem()
                        .unwrap_or_else(|| std::ffi::OsStr::new("kraken")),
                )
            });
            let result = if is_json(&case) {
                reject_json_flp(flp.as_deref()).and_then(|()| {
                    run_json_directory(
                        &case,
                        &output,
                        solver.map(Into::into),
                        format.into(),
                        max_output_bytes,
                    )
                })
            } else {
                run_legacy_directory(
                    &case,
                    &flp.unwrap_or_else(|| case.with_extension("flp")),
                    &output,
                    solver.unwrap_or(Solver::Kraken).into(),
                    format.into(),
                    max_output_bytes,
                )
            };
            result.map(|summary| {
                println!(
                    "wrote {} ({} frequencies, {} modes, {} pressures)",
                    output.display(),
                    summary.frequency_count,
                    summary.mode_count,
                    summary.pressure_count
                );
            })
        }
        Command::Export { case, flp, solver } => export(&case, flp, solver),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(match error {
                RunError::Input(_) => 2,
                RunError::Simulation { .. } => 3,
                RunError::Output(_) => 4,
            })
        }
    }
}

fn is_json(path: &Path) -> bool {
    path.extension()
        .is_some_and(|extension| extension == "json")
}

fn reject_json_flp(flp: Option<&Path>) -> Result<(), RunError> {
    if flp.is_some() {
        Err(RunError::Input(
            "error[KR0202]: --flp is invalid for self-contained JSON (flp)".into(),
        ))
    } else {
        Ok(())
    }
}

fn export(path: &Path, flp: Option<PathBuf>, solver: Option<Solver>) -> Result<(), RunError> {
    let cases = if is_json(path) {
        reject_json_flp(flp.as_deref())?;
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .and_then(|file| {
                file.take(kraken::legacy::MAX_INPUT_BYTES + 1)
                    .read_to_end(&mut bytes)
            })
            .map_err(|error| {
                RunError::Input(format!("{}: error[KR0001]: {error}", path.display()))
            })?;
        let cases = kraken::json::load_case_document_named(&bytes, path)
            .map_err(|report| RunError::Input(report.to_string()))?;
        if solver.is_some_and(|solver| {
            kraken::ModeSolver::from(solver) != cases[0].profiles()[0].mode_solver
        }) {
            return Err(RunError::Input(
                "error[KR0202]: --solver must match the JSON document (mode_solver)".into(),
            ));
        }
        cases
    } else {
        kraken::legacy::load_field_cases(
            path,
            flp.unwrap_or_else(|| path.with_extension("flp")),
            solver.unwrap_or(Solver::Kraken).into(),
        )
        .map_err(|report| RunError::Input(report.to_string()))?
    };
    let document = kraken::json::export_case_document(&cases)
        .map_err(|report| RunError::Input(report.to_string()))?;
    // Serialize completely before touching stdout: invalid/oversized exports emit no partial document.
    let mut bytes =
        serde_json::to_vec(&document).map_err(|error| RunError::Input(error.to_string()))?;
    bytes.push(b'\n');
    if bytes.len() as u64 > kraken::legacy::MAX_INPUT_BYTES {
        return Err(RunError::Input(
            "error[KR0101]: exported JSON exceeds 1 MiB (json)".into(),
        ));
    }
    std::io::stdout()
        .lock()
        .write_all(&bytes)
        .map_err(|error| RunError::Output(error.to_string()))
}
