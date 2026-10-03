use std::fs;
use std::path::Path;

use kraken::{ModeSolver, input, legacy};

#[test]
fn retained_snapshots_match_parsing_and_ignore_later_file_changes() {
    let root = std::env::temp_dir().join(format!("kraken-snapshots-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let env = root.join("profiles.env");
    let flp = root.join("geometry.flp");
    let mut source = String::new();
    for name in [
        "PekerisComplexBlank",
        "TabRefBrcC",
        "TabRefIrcC",
        "FluidTrcC",
    ] {
        source.push_str(&fs::read_to_string(fixtures.join(name).with_extension("env")).unwrap());
    }
    // Deliberate CRLF snapshot: preserve transport bytes, not just equivalent definitions.
    fs::write(&env, source.replace('\n', "\r\n")).unwrap();
    let field = fs::read_to_string(fixtures.join("TabRefBrcC.flp")).unwrap();
    fs::write(
        &flp,
        field
            .replace("'R OC'", "'RA*C'")
            .replace("1\n0.0 /", "4\n0.0 0.4 0.8 1.2 /"),
    )
    .unwrap();
    for (name, ext) in [
        ("TabRefBrcC", "brc"),
        ("TabRefIrcC", "irc"),
        ("FluidTrcC", "trc"),
    ] {
        fs::copy(
            fixtures.join(name).with_extension(ext),
            env.with_extension(ext),
        )
        .unwrap();
    }
    fs::copy(fixtures.join("FieldPattern.sbp"), flp.with_extension("sbp")).unwrap();
    // Neither ENV's SBP nor FLP's table stem is consumed.
    fs::write(env.with_extension("sbp"), b"poisoned unused source pattern").unwrap();
    fs::write(flp.with_extension("trc"), b"poisoned unused table").unwrap();
    let loaded = input::load_legacy(&env, &flp, ModeSolver::Krakenc).unwrap();
    assert_eq!(
        loaded
            .snapshots()
            .iter()
            .map(input::InputSnapshot::role)
            .collect::<Vec<_>>(),
        ["env", "flp", "trc", "brc", "irc", "sbp"]
    );
    assert_eq!(
        legacy::load_field_cases(&env, &flp, ModeSolver::Krakenc).unwrap(),
        loaded.cases()
    );
    for snapshot in loaded.snapshots() {
        assert_eq!(
            snapshot.source().as_bytes(),
            fs::read(snapshot.path()).unwrap()
        );
        fs::write(snapshot.path(), b"changed after loading").unwrap();
    }
    assert!(loaded.snapshots()[0].source().contains("\r\n"));
    let expected = legacy::load_field_cases_with_resources(
        loaded.snapshots()[0].source(),
        loaded.snapshots()[1].source(),
        &env,
        &flp,
        ModeSolver::Krakenc,
        [2, 3, 4].map(|index| Some(loaded.snapshots()[index].source())),
        Some(loaded.snapshots()[5].source()),
    )
    .unwrap();
    assert_eq!(expected, loaded.cases());
    fs::remove_dir_all(root).unwrap();
}
