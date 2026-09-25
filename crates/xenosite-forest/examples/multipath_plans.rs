//! Dump multipath find_path plans on larger / hard PhaseOne cases.
//! Tiny door mols (anisole, TBA-benzyl) are not useful here.
//!
//! ```text
//! cargo run -p xenosite-forest --example multipath_plans --release
//! cargo run -p xenosite-forest --example multipath_plans --release -- --hard
//! ```

use std::env;
use std::time::Instant;

use xenosite_forest::{FindPathConfig, PathCounters, PathOutcome, find_path_with, phase_one};

const MAX_NODES: usize = 800;
const MAX_PATHS: usize = 4;

/// Mid→large scaffolds (same set as find_path_bench `--larger`).
const LARGER: &[(&str, &str, &str)] = &[
    (
        "tBu-bis-ND→dialdehyde",
        "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
        "O=Cc1ccc(C=O)cc1",
    ),
    (
        "macrocycle-ND→aminoK",
        "C1CCCCCCNC2CCCC(CC2)NCCCC1",
        "NC1CCCC(=O)CC1",
    ),
    (
        "tribenzyl→PhCHO",
        "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
        "O=Cc1ccccc1",
    ),
    (
        "MeO-diphenyl→catechol",
        "COc1ccc(Cc2ccc(OC)cc2)cc1",
        "Oc1ccc(Cc2ccc(O)cc2)cc1",
    ),
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    (
        "eugenol→allyl-quinone",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
];

/// ≥3–8 PhaseOne hops (same set as find_path_bench `--hard`).
const HARD: &[(&str, &str, &str)] = &[
    (
        "trimethoxy-PEA→catechol",
        "COc1cc(OC)c(OC)c(CCN)c1",
        "NCCc1cc(O)c(O)c(O)c1",
    ),
    (
        "eugenol-MeO→allylQ",
        "COc1cc(CC=C)cc(OC)c1O",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "bisMeO-naph→1,2NQ",
        "COc1ccc2c(OC)cccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
    ),
    (
        "tetraMeO-biphenyl→tetraOH",
        "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC",
        "Oc1ccc(-c2ccc(O)c(O)c2)cc1O",
    ),
    (
        "veratrole-allyl→allylQ",
        "COc1ccc(CC=C)c(OC)c1OC",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "tetraMeO-naph→polyOH-NQ",
        "COc1cc(OC)c2c(OC)cc(OC)cc2c1",
        "O=C1C=C(O)C(=O)c2c(O)cc(O)cc12",
    ),
];

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
            let orbit = if s.orbit.len() <= 1 {
                String::new()
            } else {
                format!(" |Ω|={}", s.orbit.len())
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
        "  {} step(s): {}\n  precedes: {}\n  maybe ({}): {:?}\n  n_lin={}",
        hit.plan.len(),
        steps.join(" → "),
        if precedes.is_empty() {
            "(none)".into()
        } else {
            precedes.join(", ")
        },
        maybe.len(),
        maybe,
        hit.plan.n_linearizations()
    )
}

fn dump(name: &str, start: &str, target: &str) {
    println!("=== {name} ===");
    println!("reactant: {start}");
    println!("target:   {target}");
    let mut counters = PathCounters::default();
    let t0 = Instant::now();
    let hits = find_path_with(
        start,
        target,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: MAX_PATHS,
            max_nodes: MAX_NODES,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    let n_lin_sum: usize = hits.iter().map(|h| h.plan.n_linearizations()).sum();
    let mut skeleton_pairs = 0usize;
    let mut same_lin_pairs = 0usize;
    for i in 0..hits.len() {
        for j in (i + 1)..hits.len() {
            if hits[i].plan.same_linearizations(&hits[j].plan) {
                same_lin_pairs += 1;
            }
            if hits[i].plan.same_rule_maybe_skeleton(&hits[j].plan) {
                skeleton_pairs += 1;
            }
        }
    }
    println!(
        "hits={}  n_lin_sum={n_lin_sum}  same_lin_pairs={same_lin_pairs}  skeleton_pairs={skeleton_pairs}  nodes={}  mol_edits={}  drop_dup={}  signal_contained={}  plan_drops={}  wall={ms:.0} ms",
        hits.len(),
        counters.nodes,
        counters.mol_edits,
        counters.dropped_duplicate_plan,
        counters.signal_contained_plan,
        counters.plan_drops()
    );
    if hits.is_empty() {
        println!("  (no hit under budget)\n");
        return;
    }
    for (i, hit) in hits.iter().enumerate() {
        println!("\n-- plan {i} → {}", hit.smiles);
        println!("{}", fmt_plan(hit));
        if i > 0 {
            let ov = hits[0].plan.linearization_overlap(&hit.plan);
            let same = hits[0].plan.same_linearizations(&hit.plan);
            let sk = hits[0].plan.same_rule_maybe_skeleton(&hit.plan);
            println!("  vs plan0: overlap={ov} same_linearizations={same} skeleton={sk}");
        }
    }
    println!();
}

fn main() {
    let hard = env::args().any(|a| a == "--hard");
    let cases = if hard { HARD } else { LARGER };
    println!(
        "multipath plans  max_paths={MAX_PATHS}  max_nodes={MAX_NODES}  set={}\n",
        if hard { "HARD" } else { "LARGER" }
    );
    for (name, start, target) in cases {
        dump(name, start, target);
    }
}
