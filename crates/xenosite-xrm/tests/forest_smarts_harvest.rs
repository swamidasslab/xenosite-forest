//! Forest SMARTS harvest is offline data for validation — namer never imports forest.

use std::fs;
use std::path::PathBuf;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct ForestSmartsRow {
    smirks: String,
    reactant_smarts: String,
    forest_tag: Option<String>,
    rule_class: Option<String>,
}

fn harvest_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/candidates/forest-smarts.jsonl")
}

#[test]
fn forest_smarts_harvest_is_checked_in_and_well_formed() {
    let text = fs::read_to_string(harvest_path()).expect(
        "run: python3 crates/xenosite-xrm/tools/harvest_forest_smarts.py",
    );
    let mut n = 0usize;
    let mut with_tags = 0usize;
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: ForestSmartsRow = serde_json::from_str(line).expect(line);
        assert!(row.smirks.contains(">>"), "expected SMIRKS: {}", row.smirks);
        assert!(
            !row.reactant_smarts.is_empty(),
            "empty reactant SMARTS for {}",
            row.smirks
        );
        if let Some(tag) = &row.forest_tag {
            assert!(
                tag.starts_with("forest.rule:") || tag.starts_with("forest.pattern:"),
                "opaque forest tag required: {tag}"
            );
            with_tags += 1;
        }
        assert!(row.rule_class.is_some(), "missing rule_class on {}", row.smirks);
        n += 1;
    }
    assert!(n >= 50, "expected dozens of Forest SMIRKS, got {n}");
    assert!(with_tags >= 40, "expected most rows to carry forest tags, got {with_tags}");
}

#[test]
fn namer_still_has_no_forest_code_dependency() {
    // Boundary: harvest may read Forest source text; the crate must not link forest.
    let cargo = fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml"),
    )
    .unwrap();
    assert!(!cargo.contains("xenosite-forest"));
}
