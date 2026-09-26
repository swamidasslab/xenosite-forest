//! Depth-1 PhaseOne products: dump `rule\tpattern\tcsmi` (one line per unique
//! child under chematic `stable_csmi_key`, else display CSMI).
//!
//! Pair with `tests/forest/probe_d1_diff.py` — both sides rekey through
//! chematic so RDKit vs Chematic display spelling cannot fake a mismatch.

use std::collections::BTreeMap;
use std::env;

use xenosite_forest::mol::stable_csmi_key_of;
use xenosite_forest::{
    EnumConfig, ForestMol, ProductGraphConfig, enumerate_metabolites, phase_one, product_layer,
};

fn key_of_smiles(s: &str) -> String {
    stable_csmi_key_of(s).unwrap_or_else(|| s.to_string())
}

fn key_of_mol(mol: &ForestMol) -> String {
    mol.stable_csmi_key()
        .map(|k| k.as_ref().to_string())
        .unwrap_or_else(|| mol.csmi().as_ref().to_string())
}

fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    let mode = args.first().map(String::as_str).unwrap_or("dump");
    match mode {
        "dump" => {
            let name = args.get(1).map(String::as_str).unwrap_or("mol");
            let smi = args.get(2).expect("dump NAME SMILES");
            dump(name, smi);
        }
        "rekey" => {
            // stdin: rule\tpattern\tsmiles  →  rule\tpattern\tstable_key
            use std::io::{self, BufRead};
            for line in io::stdin().lock().lines() {
                let line = line.unwrap();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let parts: Vec<_> = line.splitn(3, '\t').collect();
                if parts.len() != 3 {
                    eprintln!("bad line: {line}");
                    continue;
                }
                println!("{}\t{}\t{}", parts[0], parts[1], key_of_smiles(parts[2]));
            }
        }
        "compare" => {
            // args: compare NAME py.tsv rs.tsv
            let name = args.get(1).map(String::as_str).unwrap_or("mol");
            let py_path = args.get(2).expect("compare NAME py.tsv rs.tsv");
            let rs_path = args.get(3).expect("compare NAME py.tsv rs.tsv");
            compare(name, py_path, rs_path);
        }
        other => {
            eprintln!("unknown mode {other}; use dump|rekey|compare");
            std::process::exit(2);
        }
    }
}

fn dump(name: &str, smi: &str) {
    let set = phase_one();
    let layer = ProductGraphConfig {
        target: None,
        max_nodes: usize::MAX,
        max_depth: usize::MAX,
    };
    let mol = ForestMol::parse(smi).unwrap();
    // Prefer product_layer (one hop, rule on hop) over enumerate.
    let children = product_layer(&mol, &set, &layer).unwrap();
    // key -> (rule, pattern) first-seen; also track multi-rule
    let mut by_key: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for ch in &children {
        let k = key_of_mol(&ch.child);
        by_key
            .entry(k)
            .or_default()
            .entry(ch.hop.rule.clone())
            .or_insert_with(|| ch.hop.pattern_name.clone());
    }
    eprintln!(
        "# rust {name}: layer_hops={} unique_keys={}",
        children.len(),
        by_key.len()
    );
    // Also confirm bfs d1 count
    let n_bfs = enumerate_metabolites(smi, &set, EnumConfig::bfs(1))
        .unwrap()
        .map(|h| h.unwrap())
        .count();
    eprintln!("# rust {name}: bfs_d1={n_bfs}");
    for (k, rules) in &by_key {
        // one row per (key, rule)
        for (rule, pat) in rules {
            println!("{rule}\t{pat}\t{k}");
        }
    }
}

fn load_tsv(path: &str) -> BTreeMap<String, BTreeMap<String, String>> {
    // key -> rule -> pattern
    let mut out: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for line in std::fs::read_to_string(path).unwrap().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = line.splitn(3, '\t').collect();
        if parts.len() != 3 {
            continue;
        }
        out.entry(parts[2].to_string())
            .or_default()
            .entry(parts[0].to_string())
            .or_insert_with(|| parts[1].to_string());
    }
    out
}

fn compare(name: &str, py_path: &str, rs_path: &str) {
    let py = load_tsv(py_path);
    let rs = load_tsv(rs_path);
    let py_keys: std::collections::BTreeSet<_> = py.keys().cloned().collect();
    let rs_keys: std::collections::BTreeSet<_> = rs.keys().cloned().collect();
    let both: Vec<_> = py_keys.intersection(&rs_keys).cloned().collect();
    let only_py: Vec<_> = py_keys.difference(&rs_keys).cloned().collect();
    let only_rs: Vec<_> = rs_keys.difference(&py_keys).cloned().collect();

    println!("=== {name} depth-1 (chematic-stable keys) ===");
    println!(
        "keys: py={} rs={} both={} only_py={} only_rs={}",
        py_keys.len(),
        rs_keys.len(),
        both.len(),
        only_py.len(),
        only_rs.len()
    );

    // Per-rule tallies: for each side, count keys that that rule touches,
    // split by only_py / only_rs / both.
    let mut rule_stats: BTreeMap<String, [usize; 6]> = BTreeMap::new();
    // [py_only, py_both, py_total, rs_only, rs_both, rs_total]
    let bump = |map: &mut BTreeMap<String, [usize; 6]>, rule: &str, idx: usize| {
        map.entry(rule.to_string()).or_insert([0; 6])[idx] += 1;
    };

    for k in &only_py {
        for rule in py.get(k).unwrap().keys() {
            bump(&mut rule_stats, rule, 0);
            bump(&mut rule_stats, rule, 2);
        }
    }
    for k in &both {
        for rule in py.get(k).unwrap().keys() {
            bump(&mut rule_stats, rule, 1);
            bump(&mut rule_stats, rule, 2);
        }
        for rule in rs.get(k).unwrap().keys() {
            bump(&mut rule_stats, rule, 4);
            bump(&mut rule_stats, rule, 5);
        }
    }
    for k in &only_rs {
        for rule in rs.get(k).unwrap().keys() {
            bump(&mut rule_stats, rule, 3);
            bump(&mut rule_stats, rule, 5);
        }
    }

    println!(
        "{:<28} {:>7} {:>7} {:>7}  {:>7} {:>7} {:>7}",
        "rule", "pyOnly", "pyBoth", "pyTot", "rsOnly", "rsBoth", "rsTot"
    );
    for (rule, s) in &rule_stats {
        println!(
            "{:<28} {:>7} {:>7} {:>7}  {:>7} {:>7} {:>7}",
            rule, s[0], s[1], s[2], s[3], s[4], s[5]
        );
    }

    // Sample only_py / only_rs with rules
    let show = |label: &str, keys: &[String], side: &BTreeMap<String, BTreeMap<String, String>>| {
        println!("\n--- {label} (up to 12) ---");
        for k in keys.iter().take(12) {
            let rules: Vec<_> = side[k]
                .iter()
                .map(|(r, p)| format!("{r}/{p}"))
                .collect();
            println!("  {}  {}", rules.join(","), k);
        }
        if keys.len() > 12 {
            println!("  … {} more", keys.len() - 12);
        }
    };
    show("only_py", &only_py, &py);
    show("only_rs", &only_rs, &rs);

    // Rule-pair mismatches on shared keys (same product, different rule labels)
    let mut rule_mismatch = 0usize;
    for k in &both {
        let pr: std::collections::BTreeSet<_> = py[k].keys().cloned().collect();
        let rr: std::collections::BTreeSet<_> = rs[k].keys().cloned().collect();
        if pr != rr {
            rule_mismatch += 1;
            if rule_mismatch <= 8 {
                println!(
                    "shared-key rule mismatch: py={pr:?} rs={rr:?} key={k}"
                );
            }
        }
    }
    println!("\nshared keys with different rule-sets: {rule_mismatch}");
}
