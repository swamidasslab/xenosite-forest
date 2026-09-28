//! Compare products: EpoxideHydration vs Epox→Open vs OH×2.
use std::collections::BTreeSet;

use xenosite_forest::rules::{epoxidation, epoxide_opening, hydroxylation};
use xenosite_forest::{
    FindPathConfig, HeapScoreMode, MatchScoreSpec, PathCounters, RuleSet, accept_all_rules,
    accept_all_sites, canon_of, epoxide_hydration, find_path_with, parse_mol, phase_one,
};

fn emit_products(set: &RuleSet, smi: &str) -> BTreeSet<String> {
    let mol = parse_mol(smi).unwrap();
    set.metabolize(&mol, accept_all_rules, accept_all_sites, true)
        .filter_map(|r| r.ok())
        .flat_map(|e| e.products.into_iter())
        .filter_map(|p| canon_of(&p).ok())
        .collect()
}

fn two_hop_oh(smi: &str) -> BTreeSet<String> {
    let oh = hydroxylation();
    let mol = parse_mol(smi).unwrap();
    let mut out = BTreeSet::new();
    let first: Vec<_> = oh
        .metabolize(&mol, accept_all_rules, accept_all_sites, true)
        .filter_map(|r| r.ok())
        .collect();
    for e in first {
        for p in &e.products {
            let mid = parse_mol(p).unwrap();
            for e2 in oh
                .metabolize(&mid, accept_all_rules, accept_all_sites, true)
                .filter_map(|r| r.ok())
            {
                for p2 in e2.products {
                    if let Ok(c) = canon_of(&p2) {
                        out.insert(c);
                    }
                }
            }
        }
    }
    out
}

fn epox_then_open(smi: &str) -> BTreeSet<String> {
    let epox = epoxidation();
    let open = epoxide_opening();
    let mol = parse_mol(smi).unwrap();
    let mut out = BTreeSet::new();
    for e in epox
        .metabolize(&mol, accept_all_rules, accept_all_sites, true)
        .filter_map(|r| r.ok())
    {
        for p in &e.products {
            let mid = parse_mol(p).unwrap();
            for e2 in open
                .metabolize(&mid, accept_all_rules, accept_all_sites, true)
                .filter_map(|r| r.ok())
            {
                for p2 in e2.products {
                    if let Ok(c) = canon_of(&p2) {
                        out.insert(c);
                    }
                }
            }
        }
    }
    out
}

fn section(title: &str, smi: &str) {
    println!("## {title}  ({smi})");
    let hyd = emit_products(&epoxide_hydration(), smi);
    let eo = epox_then_open(smi);
    let oh2 = two_hop_oh(smi);
    println!("  EpoxideHydration (1 hop): {hyd:?}");
    println!("  Epoxidation→EpoxideOpening: {eo:?}");
    println!("  Hydroxylation×2:           {oh2:?}");
    let hyd_eo: BTreeSet<_> = hyd.intersection(&eo).cloned().collect();
    let hyd_oh: BTreeSet<_> = hyd.intersection(&oh2).cloned().collect();
    let eo_oh: BTreeSet<_> = eo.intersection(&oh2).cloned().collect();
    println!("  ∩ hyd∩epox→open: {hyd_eo:?}");
    println!("  ∩ hyd∩OH×2:      {hyd_oh:?}");
    println!("  ∩ epox→open∩OH×2:{eo_oh:?}");
    println!();
}

fn path_probe(label: &str, reactant: &str, target: &str) {
    let mut counters = PathCounters::default();
    let rules = phase_one();
    let hits = find_path_with(
        reactant,
        target,
        &rules,
        &mut counters,
        FindPathConfig {
            max_paths: 5,
            max_nodes: 400,
            use_atom_diff: true,
            lazy_closer: true,
            heap_score: HeapScoreMode::Match(MatchScoreSpec::product_both()),
            ..FindPathConfig::default()
        },
        |_| true,
    );
    match hits {
        Ok(iter) => {
            let found = iter.collect_all().unwrap();
            println!(
                "## path {label}  hits={} bill={}",
                found.len(),
                counters.billed()
            );
            for (i, h) in found.iter().enumerate() {
                let steps: Vec<_> = h
                    .steps
                    .iter()
                    .map(|s| format!("{}:{}", s.leaf_rule().unwrap_or("?"), s.pattern_name))
                    .collect();
                println!("  plan{}: {}", i + 1, steps.join(" → "));
            }
        }
        Err(e) => println!("## path {label}  err={e:?}"),
    }
    println!();
}

fn main() {
    section("ethene", "C=C");
    section("propene", "CC=C");
    section("benzene", "c1ccccc1");
    section("phenol", "Oc1ccccc1");
    path_probe("ethene→glycol", "C=C", "OCCO");
    path_probe("benzene→catechol", "c1ccccc1", "Oc1ccccc1O");
    path_probe("benzene→dihydrodiol", "c1ccccc1", "OC1C=CC=CC1O");
}
