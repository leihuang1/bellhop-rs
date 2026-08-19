#![forbid(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bellhop::Severity;
use bellhop::solver::{SimulationLimits, run as run_simulation};
use clap::{Parser, Subcommand};

mod output;

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
    /// Run a case and write a versioned HDF5 result.
    Run {
        /// Path to the primary legacy environment file.
        case: PathBuf,
        /// Result path; defaults to `<case-stem>.h5` in the current directory.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Replace an existing result.
        #[arg(long)]
        overwrite: bool,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Validate { case } => validate(&case),
        Command::Run {
            case,
            output,
            overwrite,
        } => run(&case, output.as_deref(), overwrite),
    }
}

fn validate(path: &Path) -> ExitCode {
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

fn run(path: &Path, requested_output: Option<&Path>, overwrite: bool) -> ExitCode {
    let output_path = requested_output.map_or_else(|| default_output_path(path), Path::to_path_buf);
    if output_path.exists() && !overwrite {
        eprintln!(
            "error[BH0401]: output already exists: {}; pass --overwrite to replace it",
            output_path.display()
        );
        return ExitCode::from(4);
    }

    let case = match bellhop::legacy::load_case(path) {
        Ok(outcome) => {
            for diagnostic in outcome.warnings {
                eprintln!("{diagnostic}");
            }
            outcome.value
        }
        Err(report) => {
            render_report(&report);
            return ExitCode::from(2);
        }
    };
    let result = match run_simulation(&case, SimulationLimits::default()) {
        Ok(result) => result,
        Err(report) => {
            render_report(&report);
            return ExitCode::from(3);
        }
    };

    let temporary_path = temporary_output_path(&output_path);
    if temporary_path.exists() {
        eprintln!(
            "error[BH0402]: temporary output already exists: {}",
            temporary_path.display()
        );
        return ExitCode::from(4);
    }
    if let Err(error) = output::write_hdf5(&temporary_path, path, &result) {
        let _ = fs::remove_file(&temporary_path);
        eprintln!("error[BH0402]: {error}");
        return ExitCode::from(4);
    }
    if let Err(error) = fs::rename(&temporary_path, &output_path) {
        let _ = fs::remove_file(&temporary_path);
        eprintln!(
            "error[BH0402]: unable to atomically install {}: {error}",
            output_path.display()
        );
        return ExitCode::from(4);
    }

    let ray_count: usize = result.sources.iter().map(|source| source.rays.len()).sum();
    let eigenray_count: usize = result
        .eigenray_sources
        .iter()
        .flat_map(|source| &source.receivers)
        .map(|receiver| receiver.eigenrays.len())
        .sum();
    let arrival_count: usize = result
        .arrival_sources
        .iter()
        .flat_map(|source| &source.receivers)
        .map(|receiver| receiver.arrivals.len())
        .sum();
    let field_sample_count: usize = result
        .field_sources
        .iter()
        .map(|source| source.samples.len())
        .sum();
    let point_count: usize = result
        .sources
        .iter()
        .flat_map(|source| &source.rays)
        .chain(
            result
                .eigenray_sources
                .iter()
                .flat_map(|source| &source.receivers)
                .flat_map(|receiver| &receiver.eigenrays),
        )
        .map(|ray| ray.points.len())
        .sum();
    println!(
        "wrote {} ({ray_count} rays, {eigenray_count} eigenrays, {arrival_count} arrivals, {field_sample_count} field samples, {point_count} trajectory points)",
        output_path.display()
    );
    ExitCode::SUCCESS
}

fn default_output_path(input: &Path) -> PathBuf {
    let stem = input
        .file_stem()
        .map_or_else(|| "bellhop".into(), std::ffi::OsStr::to_os_string);
    PathBuf::from(stem).with_extension("h5")
}

fn temporary_output_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".tmp");
    PathBuf::from(name)
}

fn render_report(report: &bellhop::DiagnosticReport) {
    for diagnostic in report.diagnostics() {
        eprintln!("{diagnostic}");
    }
}
