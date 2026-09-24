//! Dump multipath find_path plans (SMILES + step plans).
//!
//! cargo run -p xenosite-forest --example multipath_plans --release

use xenosite_forest::{FindPathConfig, PathCounters, PathOutcome, find_path_with, phase_one};

fn fmt_plan(hit: &PathOutcome) -> String {
    let steps: Vec<String> = hit
        .plan
        .iter()
        .map(|s| {
            let site: Vec<_> = s
                .site
                .iter()
                .map(|a| match a {
                    xenosite_forest::PlanAtom::Index(i) => i.to_string(),
                    other => format!("{other:?}"),
                })
                .collect();
            let orbit = if s.orbit.is_empty() {
                String::new()
            } else {
                format!(" orbit={:?}", s.orbit)
            };
            format!("{}@[{}]{}", s.rule, site.join(","), orbit)
        })
        .collect();
    let precedes: Vec<String> = hit
        .plan
        .precedes()
        .iter()
        .map(|(a, b)| format!("{a}≺{b}"))
        .collect();
    let maybe: Vec<_> = hit.maybe().sides();
    format!(
        "  steps: {}\n  precedes: {}\n  maybe: {:?}\n  n_lin: {}",
        steps.join(" → "),
        if precedes.is_empty() {
            "(none)".into()
        } else {
            precedes.join(", ")
        },
        maybe,
        hit.plan.n_linearizations()
    )
}

fn dump(name: &str, start: &str, target: &str, max_paths: usize, max_nodes: usize) {
    println!("=== {name} ===");
    println!("reactant: {start}");
    println!("target:   {target}");
    let mut counters = PathCounters::default();
    let hits = find_path_with(
        start,
        target,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths,
            max_nodes,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap();
    println!(
        "hits: {}  (billed={}, nodes={})",
        hits.len(),
        counters.billed(),
        counters.nodes
    );
    for (i, hit) in hits.iter().enumerate() {
        println!("\n-- plan {i} → {}", hit.smiles);
        println!("{}", fmt_plan(hit));
        if i > 0 {
            let ov = hits[0].plan.linearization_overlap(&hit.plan);
            println!("  overlap with plan0: {ov}");
        }
    }
    println!();
}

fn main() {
    dump(
        "N,N-dimethylbenzylamine → benzaldehyde",
        "CN(C)Cc1ccccc1",
        "O=Cc1ccccc1",
        4,
        200,
    );
    dump(
        "t-butyl benzoate → benzoic acid (hydrolysis)",
        "c1ccccc1C(=O)OC(C)(C)C",
        "O=C(O)c1ccccc1",
        3,
        80,
    );
    dump(
        "anisole → phenol",
        "COc1ccccc1",
        "Oc1ccccc1",
        4,
        100,
    );
    dump(
        "hydroquinone diacetate → hydroquinone",
        "CC(=O)Oc1ccc(OC(C)=O)cc1",
        "Oc1ccc(O)cc1",
        4,
        150,
    );
    dump(
        "dimethoxy-PEA → catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
        4,
        300,
    );
}
