#![forbid(unsafe_code)]

pub mod case;
pub mod diagnostic;
pub mod json;
pub mod legacy;
pub mod model;
pub mod solver;

pub use case::{Case, CaseDefinition};
pub use diagnostic::{Diagnostic, DiagnosticReport, LoadOutcome, Severity};
pub use model::EnvironmentCase;
