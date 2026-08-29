use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use bellhop::legacy::load_case;
use bellhop::solver::{SimulationLimits, run};

static NEXT_DIRECTORY: AtomicUsize = AtomicUsize::new(0);

/// Minimal single-ray environment whose run type can be swapped.
const BASIC: &str = "'Rejection case'\n\
                     100.0\n\
                     1\n\
                     'CRW'\n\
                     2 0.0 100.0\n\
                     0.0 1500.0 /\n\
                     100.0 1500.0 /\n\
                     'V' 0.0\n\
                     1\n\
                     50.0 /\n\
                     1\n\
                     50.0 /\n\
                     1\n\
                     1.0 /\n\
                     'R'\n\
                     1\n\
                     30.0 /\n\
                     1.0 101.0 1.0\n";

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new() -> Self {
        let sequence = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "bellhop-rejection-test-{}-{sequence}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn write_case(&self, name: &str, environment: &str) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, environment).unwrap();
        path
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn error_codes(report: &bellhop::diagnostic::DiagnosticReport) -> Vec<&str> {
    report
        .diagnostics()
        .iter()
        .map(|diagnostic| diagnostic.code)
        .collect()
}

#[test]
fn ray_centered_geometric_gaussian_influence_is_rejected() {
    // Reference BELLHOP v2023.5 has no `b` influence branch; the option must
    // be rejected instead of silently falling back to another beam family.
    let directory = TemporaryDirectory::new();
    let path = directory.write_case("gaussian-b.env", &BASIC.replace("'R'\n", "'Ab'\n"));

    let case = load_case(&path).unwrap().value;
    let report = run(&case, SimulationLimits::default()).unwrap_err();
    assert_eq!(error_codes(&report), ["BH0301"]);
    assert!(report.to_string().contains("not yet supported"));
}

#[test]
fn reflection_table_generation_boundary_is_rejected() {
    // BELLHOP v2023.5's `W` option advertises IRC-generation but has no
    // generator; the solver rejects it up front.
    let directory = TemporaryDirectory::new();
    let path = directory.write_case("write-w.env", &BASIC.replace("'V' 0.0\n", "'W' 0.0\n"));
    let case = load_case(&path).unwrap().value;
    let report = run(&case, SimulationLimits::default()).unwrap_err();
    assert_eq!(error_codes(&report), ["BH0301"]);
    assert!(report.to_string().contains("table generation"));
}

#[test]
fn top_internal_reflection_table_is_rejected() {
    // The .irc format defines a bottom impedance only; a top `P` condition
    // cannot be satisfied by any auxiliary file. The top condition is the
    // second character of the SSPOPT record.
    let directory = TemporaryDirectory::new();
    let path = directory.write_case("top-p.env", &BASIC.replace("'CRW'", "'CPW'"));

    let report = load_case(&path).unwrap_err();
    assert_eq!(error_codes(&report), ["BH0202"]);
    assert!(report.to_string().contains("bottom impedance"));
}

#[test]
fn three_dimensional_runs_are_rejected() {
    let directory = TemporaryDirectory::new();
    let path = directory.write_case("three-d.env", &BASIC.replace("'R'\n", "'R    3'\n"));

    let report = load_case(&path).unwrap_err();
    assert_eq!(error_codes(&report), ["BH0202"]);
    assert!(report.to_string().contains("three-dimensional"));
}
