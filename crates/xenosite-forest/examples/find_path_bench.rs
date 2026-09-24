//! Apples-to-apples `find_path` wall times vs Python live PhaseOne.
//!
//! Same mid-size / multi-edit cases as `tests/forest/bench_find_path_h2h.py`.
//! Default runs **filter-only** (atom_diff on). Unfiltered is opt-in and
//! budgeted — hard nofilter can take minutes.
//!
//! ```text
//! cargo run -p xenosite-forest --example find_path_bench --release
//! cargo run -p xenosite-forest --example find_path_bench --release -- --larger
//! cargo run -p xenosite-forest --example find_path_bench --release -- --hard
//! cargo run -p xenosite-forest --example find_path_bench --release -- --hard --nofilter
//! ```
//!
//! Flags: `--filter-only` (default), `--nofilter`, `--eager`, `--budget-secs N`
//! (skip remaining rows once wall exceeds N; default 30 for filter, 60 with
//! `--nofilter`), `--paths N` (emit up to N plans; default 1).
//!
//! Pair with:
//! ```text
//! uv run python tests/forest/bench_find_path_rust_h2h.py
//! uv run python tests/forest/bench_find_path_rust_h2h.py --larger
//! uv run python tests/forest/bench_find_path_rust_h2h.py --hard
//! ```

use std::time::{Duration, Instant};

use xenosite_forest::{FindPathConfig, PathCounters, canon_of, find_path_with, phase_one};

const MAX_NODES: usize = 800;
const DEFAULT_MAX_PATHS: usize = 1;
const REPEATS: u32 = 5;
/// Best-of repeats when running the expensive unfiltered table.
const NOFILTER_REPEATS: u32 = 1;

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

/// Larger scaffolds with ≥3–8 PhaseOne hops (dealk + OH + DH / QF).
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

#[derive(Clone, Debug)]
struct Row {
    hit: bool,
    hits: usize,
    seconds: f64,
    steps: usize,
    nodes: usize,
    mol_edits: usize,
    billed: usize,
}

fn run_one(reactant: &str, target: &str, config: FindPathConfig, repeats: u32) -> Row {
    let want = canon_of(target).unwrap_or_else(|e| panic!("{target}: {e}"));
    let set = phase_one();

    // Warmup
    {
        let mut c = PathCounters::default();
        let _ = find_path_with(reactant, target, &set, &mut c, config, |_| true).unwrap().collect_all().unwrap();
    }

    let mut best: Option<Row> = None;
    for _ in 0..repeats {
        let mut counters = PathCounters::default();
        let t0 = Instant::now();
        let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
            .unwrap_or_else(|e| panic!("{reactant} → {target}: {e}")).collect_all().unwrap();
        let seconds = t0.elapsed().as_secs_f64();
        let hit = hits.first().is_some_and(|h| h.smiles == want);
        let steps = hits.first().map(|h| h.steps.len()).unwrap_or(0);
        let row = Row {
            hit,
            hits: hits.len(),
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

fn print_table(
    title: &str,
    cases: &[(&str, &str, &str)],
    config: FindPathConfig,
    repeats: u32,
    budget: Duration,
) {
    println!(
        "\n=== {title} (atom_diff={}, lazy_closer={}, max_paths={}) ===",
        config.use_atom_diff, config.lazy_closer, config.max_paths
    );
    println!(
        "{:<32} {:>4} {:>5} {:>9} {:>5} {:>6} {:>7} {:>6}",
        "case", "hit", "hits", "seconds", "steps", "nodes", "edits", "bill"
    );
    let mut total = 0.0;
    let suite_t0 = Instant::now();
    for &(name, reactant, target) in cases {
        if suite_t0.elapsed() >= budget {
            println!("{name:<32} SKIP  (budget {:.0}s)", budget.as_secs_f64());
            continue;
        }
        let row = run_one(reactant, target, config, repeats);
        total += row.seconds;
        println!(
            "{:<32} {:>4} {:>5} {:>9.3} {:>5} {:>6} {:>7} {:>6}",
            name,
            if row.hit { "ok" } else { "MISS" },
            row.hits,
            row.seconds,
            row.steps,
            row.nodes,
            row.mol_edits,
            row.billed
        );
    }
    println!("{:<32} {:>4} {:>5} {:>9.3}", "TOTAL", "", "", total);
}

fn parse_budget(args: &[String], default_secs: u64) -> Duration {
    args.windows(2)
        .find(|w| w[0] == "--budget-secs")
        .and_then(|w| w[1].parse().ok())
        .map(Duration::from_secs)
        .unwrap_or_else(|| Duration::from_secs(default_secs))
}

fn parse_paths(args: &[String]) -> usize {
    args.windows(2)
        .find(|w| w[0] == "--paths")
        .and_then(|w| w[1].parse().ok())
        .unwrap_or(DEFAULT_MAX_PATHS)
        .max(1)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let larger = args.iter().any(|a| a == "--larger");
    let hard = args.iter().any(|a| a == "--hard");
    // Default: filter-only. Unfiltered (slow) is --nofilter.
    let nofilter = args.iter().any(|a| a == "--nofilter");
    let filter_only = !nofilter || args.iter().any(|a| a == "--filter-only");
    let budget = parse_budget(&args, if nofilter { 60 } else { 30 });
    let max_paths = parse_paths(&args);

    println!(
        "Rust find_path PhaseOne  max_nodes={MAX_NODES}  max_paths={max_paths}  best-of-{REPEATS}"
    );
    println!(
        "(release; tagged ForestMol; filter-only={filter_only}; budget={}s)",
        budget.as_secs()
    );

    let (cases, title) = if hard {
        (HARD, "hard HA≈14–20 · multi-step (≥3–8 hops)")
    } else if larger {
        (LARGER, "larger HA≈17–26")
    } else {
        (CASES, "mid-size / multi-edit")
    };

    let base = FindPathConfig {
        max_paths,
        max_nodes: MAX_NODES,
        ..FindPathConfig::default()
    };
    if nofilter && !args.iter().any(|a| a == "--filter-only") {
        print_table(
            title,
            cases,
            FindPathConfig {
                use_atom_diff: false,
                lazy_closer: false,
                ..base
            },
            NOFILTER_REPEATS,
            budget,
        );
    }
    print_table(
        title,
        cases,
        FindPathConfig {
            use_atom_diff: true,
            lazy_closer: true,
            ..base
        },
        REPEATS,
        budget,
    );
    if hard || args.iter().any(|a| a == "--eager") {
        print_table(
            title,
            cases,
            FindPathConfig {
                use_atom_diff: true,
                lazy_closer: false,
                ..base
            },
            REPEATS,
            budget,
        );
    }
}
