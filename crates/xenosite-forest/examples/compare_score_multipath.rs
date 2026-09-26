//! Emit multiple plans under close-both / close×improve-both / improve-both.
//!
//! ```text
//! cargo run -p xenosite-forest --example compare_score_multipath --release
//! cargo run -p xenosite-forest --example compare_score_multipath --release -- --hard
//! ```

use std::env;
use std::time::Instant;

use xenosite_forest::{
    FindPathConfig, HeapScoreMode, MatchCombine, MatchMetric, MatchScoreSpec, PathCounters,
    PathOutcome, find_path_with, phase_one,
};

const MAX_NODES: usize = 800;
const MAX_PATHS: usize = 5;

const MID: &[(&str, &str, &str)] = &[
    (
        "eugenol→allyl-quinone",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    (
        "MeOPhOH→hydroxyquinone",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        r"C(#C/C=C/C=O)C(C)(C)C",
    ),
    (
        "2-MeO-naph→1,2-NQ",
        "COc1ccc2ccccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
    ),
];

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

fn scores() -> [HeapScoreMode; 3] {
    [
        HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Close,
            metric: MatchMetric::Both,
        }),
        HeapScoreMode::Match(MatchScoreSpec::product_both()),
        HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Improve,
            metric: MatchMetric::Both,
        }),
    ]
}

fn step_line(hit: &PathOutcome) -> String {
    hit.steps
        .iter()
        .map(|s| {
            format!(
                "{}:{}",
                s.leaf_rule().unwrap_or("?"),
                s.pattern_name
            )
        })
        .collect::<Vec<_>>()
        .join(" → ")
}

fn run_case(name: &str, reactant: &str, target: &str, score: HeapScoreMode) {
    let mut counters = PathCounters::default();
    let t0 = Instant::now();
    let hits = find_path_with(
        reactant,
        target,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: MAX_PATHS,
            max_nodes: MAX_NODES,
            use_atom_diff: true,
            lazy_closer: true,
            heap_score: score,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    let hops: Vec<usize> = hits.iter().map(|h| h.steps.len()).collect();
    println!(
        "  [{:<20}] hits={} hops={hops:?} nodes={} edits={} bill={} wall={ms:.0}ms",
        score.label(),
        hits.len(),
        counters.nodes,
        counters.mol_edits,
        counters.billed(),
    );
    for (i, hit) in hits.iter().enumerate() {
        println!(
            "    plan{} ({}h n_lin={}): {}",
            i + 1,
            hit.steps.len(),
            hit.plan.n_linearizations(),
            step_line(hit)
        );
    }
    if hits.is_empty() {
        println!("    (no hit)");
    }
    let _ = name;
}

fn main() {
    let hard = env::args().any(|a| a == "--hard");
    let cases = if hard { HARD } else { MID };
    let title = if hard { "HARD" } else { "MID" };
    println!(
        "compare_score_multipath — {title}  max_paths={MAX_PATHS} max_nodes={MAX_NODES}\n\
         scores: close-both | close×improve-both | improve-both\n\
         close = nearer child; improve = hop cut; close×improve = both factors; both = atom+formula\n"
    );
    for &(name, reactant, target) in cases {
        println!("## {name}");
        for mode in scores() {
            run_case(name, reactant, target, mode);
        }
        println!();
    }
}
