//! Tagging-layer competency fixtures exported from data/competency/cq.yml.

use std::fs;
use std::path::PathBuf;

use serde::Deserialize;
use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

#[derive(Debug, Deserialize)]
struct Fixture {
    id: String,
    reactant_smiles: String,
    product_smiles: String,
    #[serde(default)]
    tags: Vec<String>,
    #[serde(default)]
    expected_labels: Vec<String>,
    #[serde(default)]
    expected_excludes_labels: Vec<String>,
    #[serde(default)]
    expect_site_map: Vec<u16>,
    #[serde(default)]
    expect_site_label_substring: Option<String>,
}

fn fixtures() -> Vec<Fixture> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data/competency/fixtures/reactions.jsonl");
    fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("missing {}; run tools/run_competency_tests.py", path.display()))
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| serde_json::from_str(l).expect(l))
        .collect()
}

#[test]
fn tagging_competency_fixtures() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    let cases = fixtures();
    assert!(
        cases.len() >= 10,
        "expected tagging fixtures, got {}",
        cases.len()
    );
    for fx in cases {
        let tag_refs: Vec<&str> = fx.tags.iter().map(String::as_str).collect();
        let terms = namer
            .name_smiles(&fx.reactant_smiles, &fx.product_smiles, &tag_refs)
            .unwrap_or_else(|e| panic!("{}: {e}", fx.id));
        let labs: Vec<&str> = terms.iter().map(|t| t.pref_label.as_str()).collect();
        for exp in &fx.expected_labels {
            assert!(
                labs.contains(&exp.as_str()),
                "{} missing label {exp:?} in {labs:?}",
                fx.id
            );
        }
        for bad in &fx.expected_excludes_labels {
            assert!(
                !labs.contains(&bad.as_str()),
                "{} unexpectedly has {bad:?} in {labs:?}",
                fx.id
            );
        }
        if !fx.expect_site_map.is_empty() {
            let hit = terms.iter().find(|t| {
                fx.expected_labels.iter().any(|l| l == &t.pref_label)
                    && t.site.map_nums == fx.expect_site_map
            });
            assert!(
                hit.is_some(),
                "{} expected site maps {:?} on one of {:?}; terms={:?}",
                fx.id,
                fx.expect_site_map,
                fx.expected_labels,
                terms
                    .iter()
                    .map(|t| (&t.pref_label, &t.site.map_nums))
                    .collect::<Vec<_>>()
            );
        }
        if let Some(sub) = &fx.expect_site_label_substring {
            let bundles = namer
                .annotate_smiles(&fx.reactant_smiles, &fx.product_smiles, &tag_refs)
                .unwrap_or_else(|e| panic!("{} annotate: {e}", fx.id));
            let ok = bundles.iter().any(|b| {
                b.site_label
                    .as_deref()
                    .map(|s| s.contains(sub.as_str()))
                    .unwrap_or(false)
                    || terms.iter().any(|t| t.site.display_suffix().contains(sub))
            });
            assert!(
                ok,
                "{} expected site_label/display containing {sub:?}; bundles={:?} terms={:?}",
                fx.id,
                bundles
                    .iter()
                    .map(|b| &b.site_label)
                    .collect::<Vec<_>>(),
                terms
                    .iter()
                    .map(|t| t.site.display_suffix())
                    .collect::<Vec<_>>()
            );
        }
    }
}
