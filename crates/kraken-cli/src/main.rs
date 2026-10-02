#![forbid(unsafe_code)]

use std::path::PathBuf;
use std::process::ExitCode;

use bellhop_hdf5::kraken::{DEFAULT_MAX_OUTPUT_BYTES, RunError, run_legacy};
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Parser)]
#[command(
    name = "kraken",
    version,
    about = "Supported 2D KRAKEN/KRAKENC modes and FIELD"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Run a supported legacy .env/.flp pair and atomically write KRAKEN HDF5 v1.
    Run {
        /// Legacy .env input; JSON inputs are not supported.
        case: PathBuf,
        /// FIELD geometry; defaults to the same-stem .flp.
        #[arg(long)]
        flp: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = Solver::Kraken)]
        solver: Solver,
        /// Defaults to <case-stem>.h5 in the current directory.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Atomically replace an existing result, never an input file.
        #[arg(long)]
        overwrite: bool,
        /// Cumulative file/payload quota, checked before datasets and after flushes.
        #[arg(long, default_value_t = DEFAULT_MAX_OUTPUT_BYTES, value_parser = clap::value_parser!(u64).range(1..))]
        max_output_bytes: u64,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Solver {
    Kraken,
    Krakenc,
}

fn main() -> ExitCode {
    let Command::Run {
        case,
        flp,
        solver,
        output,
        overwrite,
        max_output_bytes,
    } = Cli::parse().command;
    let flp = flp.unwrap_or_else(|| case.with_extension("flp"));
    let output = output.unwrap_or_else(|| {
        PathBuf::from(
            case.file_stem()
                .unwrap_or_else(|| std::ffi::OsStr::new("kraken")),
        )
        .with_extension("h5")
    });
    let solver = match solver {
        Solver::Kraken => kraken::ModeSolver::Kraken,
        Solver::Krakenc => kraken::ModeSolver::Krakenc,
    };
    match run_legacy(&case, &flp, &output, solver, overwrite, max_output_bytes) {
        Ok(summary) => {
            println!(
                "wrote {} ({} frequencies, {} modes, {} pressures)",
                output.display(),
                summary.frequency_count,
                summary.mode_count,
                summary.pressure_count
            );
            ExitCode::SUCCESS
        }
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
