//! Emit actual namer labels for curated gold reactions (competency scoring).

use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

fn gold_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/competency/gold/curated_reactions.tsv")
}

fn out_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/competency/gold/actual_tags.jsonl")
}

#[test]
fn emit_gold_actual_tags() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    let text = fs::read_to_string(gold_path()).expect("curated_reactions.tsv");
    let mut lines = text.lines();
    let header = lines.next().expect("header");
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .unwrap_or_else(|| panic!("missing column {name}"))
    };
    let i_id = idx("id");
    let i_r = idx("reactant_smiles");
    let i_p = idx("product_smiles");
    let i_tags = idx("tags");

    let mut out = File::create(out_path()).unwrap();
    let mut n = 0usize;
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        let id = parts[i_id];
        let reactant = parts[i_r];
        let product = parts[i_p];
        let tags_field = parts.get(i_tags).copied().unwrap_or("");
        let tags: Vec<&str> = tags_field
            .split('|')
            .map(str::trim)
            .filter(|t| !t.is_empty())
            .collect();
        let terms = namer
            .name_smiles(reactant, product, &tags)
            .unwrap_or_else(|e| panic!("{id}: {e}"));
        let labels: Vec<String> = terms.iter().map(|t| t.pref_label.clone()).collect();
        let site_maps: Vec<Vec<u16>> = terms
            .iter()
            .filter(|t| !t.site.map_nums.is_empty())
            .map(|t| t.site.map_nums.clone())
            .collect();
        let external_mappings: Vec<String> = terms
            .iter()
            .flat_map(|t| t.inter_matches.iter().map(|m| m.target.0.clone()))
            .collect();
        let row = serde_json::json!({
            "id": id,
            "labels": labels,
            "site_maps": site_maps,
            "external_mappings": external_mappings,
        });
        writeln!(out, "{row}").unwrap();
        n += 1;
    }
    assert!(n >= 8, "expected curated gold rows, got {n}");
}
