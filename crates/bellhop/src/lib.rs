#![forbid(unsafe_code)]

pub mod diagnostic;
pub mod legacy;
pub mod model;

pub use diagnostic::{Diagnostic, DiagnosticReport, LoadOutcome, Severity};
pub use model::{Case, EnvironmentCase};
