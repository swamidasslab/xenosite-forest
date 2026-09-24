//! Apples-to-apples `find_path` wall times vs Python live PhaseOne.
//!
//! Same mid-size / multi-edit cases as `tests/forest/bench_find_path_h2h.py`.
//!
//! ```text
//! cargo run -p xenosite-forest --example find_path_bench --release
//! ```
//!
//! Pair with:
//! ```text
//! uv run python tests/forest/bench_find_path_rust_h2h.py
//! ```

use std::time::Instant;

use xenosite_forest::{FindPathConfig, PathCounters, canon_of, find_path_with, phase_one};

const MAX_NODES: usize = 800;
const MAX_PATHS: usize = 1;
const REPEATS: u32 = 5;

/// (name, reactant, target). Targets are chematic-reachable CSMI spellings
/// (TBA keeps E stereo; non-stereo spelling is a different canonaut string).
const CASES: &[(&str, &str, &str)] = &[
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
        "triPh-butyl→OH",
        "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
        "Oc1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
    ),
    (
        "MeO-diphenyl→catechol",
        "COc1ccc(Cc2ccc(OC)cc2)cc1",
        "Oc1ccc(Cc2ccc(O)cc2)cc1",
    ),
];

#[derive(Clone, Debug)]
struct Row {
    hit: bool,
    seconds: f64,
    steps: usize,
    nodes: usize,
    mol_edits: usize,
    billed: usize,
}

fn run_one(reactant: &str, target: &str, use_atom_diff: bool) -> Row {
    let want = canon_of(target).unwrap_or_else(|e| panic!("{target}: {e}"));
    let set = phase_one();
    let config = FindPathConfig {
        max_paths: MAX_PATHS,
        max_nodes: MAX_NODES,
        use_atom_diff,
    };

    // Warmup
    {
        let mut c = PathCounters::default();
        let _ = find_path_with(reactant, target, &set, &mut c, config, |_| true);
    }

    let mut best: Option<Row> = None;
    for _ in 0..REPEATS {
        let mut counters = PathCounters::default();
        let t0 = Instant::now();
        let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
            .unwrap_or_else(|e| panic!("{reactant} → {target}: {e}"));
        let seconds = t0.elapsed().as_secs_f64();
        let hit = hits.first().is_some_and(|h| h.smiles == want);
        let steps = hits.first().map(|h| h.steps.len()).unwrap_or(0);
        let row = Row {
            hit,
            seconds,
            steps,
            nodes: counters.nodes,
            mol_edits: counters.mol_edits,
            billed: counters.billed(),
        };
        best = Some(match best {
            None => row,
            Some(prev) if row.seconds < prev.seconds => row,
            Some(prev) => prev,
        });
    }
    best.expect("repeats")
}

fn print_table(title: &str, cases: &[(&str, &str, &str)], use_atom_diff: bool) {
    println!("\n=== {title} (atom_diff={use_atom_diff}) ===");
    println!(
        "{:<28} {:>4} {:>9} {:>5} {:>6} {:>7} {:>6}",
        "case", "hit", "seconds", "steps", "nodes", "edits", "bill"
    );
    let mut total = 0.0;
    for &(name, reactant, target) in cases {
        let row = run_one(reactant, target, use_atom_diff);
        total += row.seconds;
        println!(
            "{:<28} {:>4} {:>9.3} {:>5} {:>6} {:>7} {:>6}",
            name,
            if row.hit { "ok" } else { "MISS" },
            row.seconds,
            row.steps,
            row.nodes,
            row.mol_edits,
            row.billed
        );
    }
    println!("{:<28} {:>4} {:>9.3}", "TOTAL", "", total);
}

fn main() {
    let larger = std::env::args().any(|a| a == "--larger");
    println!(
        "Rust find_path PhaseOne  max_nodes={MAX_NODES}  max_paths={MAX_PATHS}  best-of-{REPEATS}"
    );
    println!("(release; chematic door; provisional atom_diff optional)");

    let cases = if larger { LARGER } else { CASES };
    let title = if larger {
        "larger HA≈17–26"
    } else {
        "mid-size / multi-edit"
    };

    // Default live Python uses filters (atom_diff). Report both modes.
    print_table(title, cases, false);
    print_table(title, cases, true);
}
