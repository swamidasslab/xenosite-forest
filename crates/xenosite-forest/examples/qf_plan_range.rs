//! Enumerate QuinoneFormation emission plans (canonical step rule sequences).
use std::collections::{BTreeMap, BTreeSet};

use xenosite_forest::rules::quinone_formation;
use xenosite_forest::{accept_all_rules, accept_all_sites, canon_of, parse_mol};

fn plan_key(plan: &[xenosite_forest::Step]) -> String {
    plan.iter()
        .map(|s| {
            let sites: Vec<_> = s.site.iter().map(|a| format!("{a:?}")).collect();
            format!("{}[{}]", s.rule, sites.join(","))
        })
        .collect::<Vec<_>>()
        .join(" → ")
}

fn plan_rules(plan: &[xenosite_forest::Step]) -> String {
    plan.iter()
        .map(|s| s.rule.as_str())
        .collect::<Vec<_>>()
        .join(" → ")
}

fn probe(smi: &str) {
    let set = quinone_formation();
    let mol = parse_mol(smi).unwrap();
    let emissions: Vec<_> = set
        .metabolize(&mol, accept_all_rules, accept_all_sites, true)
        .filter_map(|r| r.ok())
        .collect();
    println!("## {smi}  emissions={}", emissions.len());
    let mut by_shape: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for e in &emissions {
        let shape = plan_rules(&e.plan);
        let detail = format!(
            "pattern={} site={:?} plan={}",
            e.pattern_name,
            e.site_atoms,
            plan_key(&e.plan)
        );
        by_shape.entry(shape).or_default().insert(detail);
        let products: Vec<_> = e.products.iter().filter_map(|p| canon_of(p).ok()).collect();
        println!(
            "  [{}] {}  products={products:?}",
            e.pattern_name,
            plan_rules(&e.plan)
        );
    }
    println!("  shapes:");
    for (shape, details) in &by_shape {
        println!("    {shape}  (n={})", details.len());
    }
    println!();
}

fn main() {
    for smi in [
        "c1ccccc1",       // benzene: add_carbonyl_o ×2
        "Oc1ccccc1",      // phenol: phenol end + add O / etc
        "Oc1ccc(O)cc1",   // hydroquinone: DH only (no prep)
        "Oc1ccc(N)cc1",   // aminophenol
        "Cc1ccccc1",      // toluene: methide
        "Clc1ccccc1",     // chlorobenzene: replace_halogen
        "COc1ccccc1",     // anisole: dealkylate
        "Oc1ccc(Cl)cc1",  // para-chlorophenol
        "Oc1ccc(OC)cc1",  // MeOPhOH
        "c1ccc2ccccc2c1", // naphthalene
    ] {
        probe(smi);
    }
}
