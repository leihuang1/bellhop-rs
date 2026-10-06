use std::error::Error;
use std::fmt;
use std::path::PathBuf;

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
    pub(crate) diagnostics: Vec<Diagnostic>,
}

impl DiagnosticReport {
    #[must_use]
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    pub(crate) fn push(&mut self, diagnostic: Diagnostic) {
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
