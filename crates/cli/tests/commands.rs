use std::process::Command;

#[test]
fn unified_cli_routes_solver_commands_without_old_root_aliases_or_new_formats() {
    for args in [
        vec!["--help"],
        vec!["bellhop", "--help"],
        vec!["bellhop", "validate", "--help"],
        vec!["bellhop", "export", "--help"],
        vec!["bellhop", "run", "--help"],
        vec!["kraken", "--help"],
        vec!["kraken", "export", "--help"],
        vec!["kraken", "run", "--help"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_pelagic"))
            .args(&args)
            .output()
            .unwrap();
        assert!(output.status.success(), "{args:?}");
        assert!(String::from_utf8_lossy(&output.stdout).contains("pelagic"));
    }
    for args in [
        vec!["run", "case.env"],
        vec!["bellhop", "run", "case.env", "--format", "legacy"],
        vec!["kraken", "run", "case.env", "--format", "both"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_pelagic"))
            .args(&args)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2), "{args:?}");
    }
}
