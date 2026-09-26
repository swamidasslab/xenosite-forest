//! Forest SMARTS harvest validates XRM coverage — namer never imports forest.

use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;

use serde::Deserialize;
use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

#[derive(Debug, Deserialize)]
struct ForestSmartsRow {
    smirks: String,
    reactant_smarts: String,
    forest_tag: Option<String>,
    rule_class: Option<String>,
}

#[derive(Debug, Deserialize)]
struct CoverageReport {
    n_smirks: usize,
    n_rule_classes: usize,
    classes_with_sssom: usize,
    classes_assignment_taggable: usize,
    classes: Vec<CoverageClass>,
}

#[derive(Debug, Deserialize)]
struct CoverageClass {
    rule_class: String,
    sssom_rule_mapped: bool,
    #[allow(dead_code)]
    assignment_taggable: bool,
    #[serde(default)]
    patterns_sssom_mapped: Vec<String>,
}

fn data(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

#[test]
fn forest_smarts_harvest_is_checked_in_and_well_formed() {
    let text = fs::read_to_string(data("data/candidates/forest-smarts.jsonl")).expect(
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
        assert!(
            row.rule_class.is_some(),
            "missing rule_class on {}",
            row.smirks
        );
        n += 1;
    }
    assert!(n >= 50, "expected dozens of Forest SMIRKS, got {n}");
    assert!(
        with_tags >= 40,
        "expected most rows to carry forest tags, got {with_tags}"
    );
}

#[test]
fn forest_coverage_report_meets_policy() {
    let path = data("data/candidates/forest-coverage.json");
    let text = fs::read_to_string(&path).unwrap_or_else(|_| {
        panic!(
            "missing {}; run harvest_forest_smarts.py or validate_forest_coverage.py",
            path.display()
        )
    });
    let report: CoverageReport = serde_json::from_str(&text).expect("coverage json");
    assert!(report.n_smirks >= 50, "smirks={}", report.n_smirks);
    assert!(
        report.n_rule_classes >= 15,
        "rule classes={}",
        report.n_rule_classes
    );
    // Most Forest rule classes should already have SSSOM and/or tag pathways.
    let sssom_frac = report.classes_with_sssom as f64 / report.n_rule_classes as f64;
    let tag_frac =
        report.classes_assignment_taggable as f64 / report.n_rule_classes as f64;
    assert!(
        sssom_frac >= 0.5,
        "SSSOM class coverage too low: {}/{} ({sssom_frac:.2})",
        report.classes_with_sssom,
        report.n_rule_classes
    );
    assert!(
        tag_frac >= 0.35,
        "assignment-tag class coverage too low: {}/{} ({tag_frac:.2})",
        report.classes_assignment_taggable,
        report.n_rule_classes
    );
}

#[test]
fn opaque_forest_tags_from_harvest_name_without_importing_forest() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    let text = fs::read_to_string(data("data/candidates/forest-smarts.jsonl")).unwrap();
    let report: CoverageReport = serde_json::from_str(
        &fs::read_to_string(data("data/candidates/forest-coverage.json")).unwrap(),
    )
    .unwrap();

    // Prefer tags that SSSOM already maps (stable naming path).
    let mut tags: BTreeSet<String> = BTreeSet::new();
    for c in &report.classes {
        if c.sssom_rule_mapped {
            tags.insert(format!("forest.rule:{}", c.rule_class));
        }
        for p in &c.patterns_sssom_mapped {
            tags.insert(p.clone());
        }
    }
    // Also include a few raw harvest tags that assignments listen for.
    for line in text.lines().take(200) {
        if line.is_empty() {
            continue;
        }
        let row: ForestSmartsRow = serde_json::from_str(line).unwrap();
        if let Some(tag) = row.forest_tag {
            if tag.contains("Hydroxylation")
                || tag.contains("Epoxidation")
                || tag.contains("Glucuronidation")
                || tag.contains("NDealkylation")
                || tag.contains("Glutathionation")
            {
                tags.insert(tag);
            }
        }
    }

    assert!(tags.len() >= 10, "expected mapped tags to probe, got {}", tags.len());
    let mut named = 0usize;
    for tag in &tags {
        let terms = namer
            .name_smiles("C", "C", &[tag.as_str()])
            .unwrap_or_else(|e| panic!("tag {tag} must not error: {e}"));
        // SSSOM-mapped Forest-map concepts or chemist emits should surface.
        if !terms.is_empty() {
            named += 1;
        }
    }
    assert!(
        named as f64 / tags.len() as f64 >= 0.6,
        "too few Forest tags produced names: {named}/{}",
        tags.len()
    );
}

#[test]
fn namer_still_has_no_forest_code_dependency() {
    // Boundary: harvest may read Forest source text; the crate must not link forest.
    let cargo = fs::read_to_string(data("Cargo.toml")).unwrap();
    let deps = cargo
        .split("[dependencies]")
        .nth(1)
        .unwrap_or("")
        .split('[')
        .next()
        .unwrap_or("");
    assert!(
        !deps.contains("xenosite-forest") && !deps.contains("xenosite_forest"),
        "dependencies must not include forest crates:\n{deps}"
    );
}
