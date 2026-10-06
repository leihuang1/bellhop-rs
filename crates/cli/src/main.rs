#![forbid(unsafe_code)]

mod bellhop;
mod kraken;

use std::process::ExitCode;

use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(
    name = "pelagic",
    version,
    about = "Two-dimensional underwater acoustics"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// BELLHOP ray, arrival and pressure-field calculations.
    #[command(subcommand)]
    Bellhop(bellhop::Command),
    /// KRAKEN/KRAKENC normal modes and FIELD calculations.
    #[command(subcommand)]
    Kraken(kraken::Command),
}

fn main() -> ExitCode {
    match Cli::parse().command {
        Command::Bellhop(command) => bellhop::execute(command),
        Command::Kraken(command) => kraken::execute(command),
    }
}
