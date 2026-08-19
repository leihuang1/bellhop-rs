#![forbid(unsafe_code)]

use std::path::PathBuf;
use std::process::ExitCode;

use bellhop::Severity;
use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(
    name = "bellhop",
    version = concat!(env!("CARGO_PKG_VERSION"), " (Rust implementation)"),
    about = "Modern two-dimensional BELLHOP"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Parse and validate a legacy BELLHOP case without running it.
    Validate {
        /// Path to the primary legacy environment file.
        case: PathBuf,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Validate { case } => validate(&case),
    }
}

fn validate(path: &std::path::Path) -> ExitCode {
    match bellhop::legacy::load_case(path) {
        Ok(outcome) => {
            for diagnostic in outcome.warnings {
                eprintln!("{diagnostic}");
            }
            println!(
                "valid 2D BELLHOP case: {} ({:?}, {} launch angles)",
                path.display(),
                outcome.value.environment.run.kind,
                outcome.value.environment.trace.launch_angles_degrees.len()
            );
            ExitCode::SUCCESS
        }
        Err(report) => {
            for diagnostic in report.diagnostics() {
                eprintln!("{diagnostic}");
            }
            if report
                .diagnostics()
                .iter()
                .any(|diagnostic| diagnostic.severity == Severity::Error)
            {
                ExitCode::from(2)
            } else {
                ExitCode::SUCCESS
            }
        }
    }
}
