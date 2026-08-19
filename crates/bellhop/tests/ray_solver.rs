use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use bellhop::legacy::load_case;
use bellhop::solver::{RayTermination, SimulationLimits, run};

static NEXT_DIRECTORY: AtomicUsize = AtomicUsize::new(0);

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new() -> Self {
        let sequence = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "bellhop-ray-test-{}-{sequence}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn join(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn traces_official_munk_ray_fan_deterministically() {
    let case = load_case(&fixture("MunkB_ray.env")).unwrap().value;
    let first = run(&case, SimulationLimits::default()).unwrap();
    let second = run(&case, SimulationLimits::default()).unwrap();

    assert_eq!(first, second);
    assert_eq!(first.sources.len(), 2);
    assert_eq!(first.sources[0].rays.len(), 41);
    assert!(
        first
            .sources
            .iter()
            .flat_map(|source| &source.rays)
            .all(|ray| !ray.points.is_empty())
    );
    assert!(
        first
            .sources
            .iter()
            .flat_map(|source| &source.rays)
            .all(|ray| ray.termination != RayTermination::StepLimit)
    );
}

#[test]
fn rejects_ray_fans_above_the_configured_limit() {
    let case = load_case(&fixture("MunkB_ray.env")).unwrap().value;
    let report = run(
        &case,
        SimulationLimits {
            max_rays: 1,
            ..SimulationLimits::default()
        },
    )
    .unwrap_err();
    assert_eq!(report.diagnostics()[0].code, "BH0303");
}

#[test]
fn applies_tabulated_and_material_boundary_reflections() {
    for (name, bottom_options, bottom_data) in [
        ("table", "'F' 0.0", ""),
        ("elastic", "'A' 0.0", "100.0 1800.0 600.0 2.0 0.1 0.2 /\n"),
        ("grain", "'G' 0.0", "100.0 1.5\n"),
    ] {
        let directory = TemporaryDirectory::new();
        let environment = directory.join(&format!("{name}.env"));
        fs::write(
            &environment,
            format!(
                "'Boundary reflection'\n\
                 100.0\n\
                 1\n\
                 'CRW'\n\
                 2 0.0 100.0\n\
                 0.0 1500.0 /\n\
                 100.0 1500.0 /\n\
                 {bottom_options}\n\
                 {bottom_data}\
                 1\n\
                 50.0 /\n\
                 2\n\
                 0.0 100.0 /\n\
                 2\n\
                 0.0 1.0 /\n\
                 'R'\n\
                 1\n\
                 30.0 /\n\
                 1.0 101.0 1.0\n"
            ),
        )
        .unwrap();
        if name == "table" {
            fs::write(
                directory.join("table.brc"),
                "2\n0.0 0.001 0.0\n90.0 0.001 0.0\n",
            )
            .unwrap();
        }

        let case = load_case(&environment).unwrap().value;
        let result = run(&case, SimulationLimits::default()).unwrap();
        let ray = &result.sources[0].rays[0];
        assert!(ray.bottom_bounces >= 1, "{name}");
        if name == "table" {
            assert_eq!(ray.termination, RayTermination::LostEnergy);
        }
    }
}
