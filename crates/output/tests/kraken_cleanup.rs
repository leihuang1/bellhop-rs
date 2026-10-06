use std::fs;
use std::path::Path;

use kraken::ModeSolver;
use output::kraken::{DEFAULT_MAX_OUTPUT_BYTES, RunError, run_legacy};

#[test]
fn failed_runs_release_scratch_before_retry_in_the_same_process() {
    let directory = std::env::temp_dir().join(format!("kraken-cleanup-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let fixture =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../kraken/tests/fixtures/PekerisBroadband");
    let source = fs::read_to_string(fixture.with_extension("env")).unwrap();
    let env = directory.join("case.env");
    let flp = directory.join("case.flp");
    let output = directory.join("result.h5");
    let scratch = directory.join("result.h5.tmp");
    fs::copy(fixture.with_extension("flp"), &flp).unwrap();
    let run = |quota| run_legacy(&env, &flp, &output, ModeSolver::Kraken, true, quota);
    for simulation_failure in [false, true] {
        fs::write(&output, b"previous complete result").unwrap();
        fs::write(
            &env,
            if simulation_failure {
                source.replace("75.0 50.0 62.5 /", "75.0 5.0 62.5 /")
            } else {
                source.clone()
            },
        )
        .unwrap();
        let error = run(if simulation_failure {
            DEFAULT_MAX_OUTPUT_BYTES
        } else {
            4096
        })
        .unwrap_err();
        if simulation_failure {
            assert!(matches!(
                error,
                RunError::Simulation {
                    frequency_index: 1,
                    ..
                }
            ));
        } else {
            assert!(
                matches!(error, RunError::Output(ref message) if message.contains("HDF5 file exceeds"))
            );
        }
        assert!(!scratch.exists(), "failed run left its scratch file");
        assert_eq!(fs::read(&output).unwrap(), b"previous complete result");
        // Do not rely on subprocess exit to release HDF5 resources before retry.
        fs::write(&env, &source).unwrap();
        assert_eq!(run(DEFAULT_MAX_OUTPUT_BYTES).unwrap().frequency_count, 3);
        assert!(!scratch.exists());
    }
    fs::remove_dir_all(directory).unwrap();
}
