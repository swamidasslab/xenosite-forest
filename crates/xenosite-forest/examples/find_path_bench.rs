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
//! `--nofilter`), `--paths N` (emit up to N plans; default 1),
//! `--score LABEL` (`soft` or `combine-metric` e.g. `close×improve-both`),
//! `--matrix` (run all 9 match variants + soft; summary ranked by miss/bill/time),
//! `--mcs-extend` (opt-in: atom_diff uses MCS+placeable grow; default bare MCS).
//!
//! Pair with:
//! ```text
//! uv run python tests/forest/bench_find_path_rust_h2h.py
//! uv run python tests/forest/bench_find_path_rust_h2h.py --larger
//! uv run python tests/forest/bench_find_path_rust_h2h.py --hard
//! ```

use std::time::{Duration, Instant};

use xenosite_forest::{
    FindPathConfig, HeapScoreMode, MatchCombine, MatchMetric, MatchScoreSpec, PathCounters,
    canon_of, find_path_with, phase_one, set_use_mcs_extend, use_mcs_extend,
};

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
        let _ = find_path_with(reactant, target, &set, &mut c, config, |_| true)
            .unwrap()
            .collect_all()
            .unwrap();
    }

    let mut best: Option<Row> = None;
    for _ in 0..repeats {
        let mut counters = PathCounters::default();
        let t0 = Instant::now();
        let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
            .unwrap_or_else(|e| panic!("{reactant} → {target}: {e}"))
            .collect_all()
            .unwrap();
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

fn parse_score_label(label: &str) -> HeapScoreMode {
    match label {
        "soft" | "soft-stack" | "legacy" => HeapScoreMode::SoftStack,
        "match" | "default" | "log-neg-pc" => HeapScoreMode::match_log_neg_pc(),
        "add" | "add-both" => HeapScoreMode::match_add(),
        "product"
        | "product-both"
        | "match-product"
        | "close×improve"
        | "close×improve-both"
        | "close-x-improve"
        | "close-x-improve-both" => HeapScoreMode::Match(MatchScoreSpec::product_both()),
        "close" | "dist" | "close-both" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Close,
            metric: MatchMetric::Both,
        }),
        "improve" | "imp" | "improve-both" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Improve,
            metric: MatchMetric::Both,
        }),
        "atom" | "add-atom" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Add,
            metric: MatchMetric::Atom,
        }),
        "formula" | "add-formula" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Add,
            metric: MatchMetric::Formula,
        }),
        "product-atom" | "close×improve-atom" | "close-x-improve-atom" => {
            HeapScoreMode::Match(MatchScoreSpec {
                combine: MatchCombine::Product,
                metric: MatchMetric::Atom,
            })
        }
        "product-formula" | "close×improve-formula" | "close-x-improve-formula" => {
            HeapScoreMode::Match(MatchScoreSpec {
                combine: MatchCombine::Product,
                metric: MatchMetric::Formula,
            })
        }
        "close-atom" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Close,
            metric: MatchMetric::Atom,
        }),
        "close-formula" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Close,
            metric: MatchMetric::Formula,
        }),
        "improve-atom" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Improve,
            metric: MatchMetric::Atom,
        }),
        "improve-formula" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::Improve,
            metric: MatchMetric::Formula,
        }),
        "lin-neg-c" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::LinNegC,
            metric: MatchMetric::Both,
        }),
        "lin-neg-pc" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::LinNegPC,
            metric: MatchMetric::Both,
        }),
        "lin-p-neg2c" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::LinPNeg2C,
            metric: MatchMetric::Both,
        }),
        "lin-neg-p-neg2c" => HeapScoreMode::Match(MatchScoreSpec {
            combine: MatchCombine::LinNegPNeg2C,
            metric: MatchMetric::Both,
        }),
        other => panic!(
            "unknown --score {other} (soft|add-both|product-both|log-neg-pc|close×improve-both|lin-*|…)"
        ),
    }
}

fn parse_score(args: &[String]) -> HeapScoreMode {
    args.windows(2)
        .find(|w| w[0] == "--score")
        .map(|w| parse_score_label(&w[1]))
        .unwrap_or_else(HeapScoreMode::match_log_neg_pc)
}

#[derive(Clone, Debug)]
struct SuiteSummary {
    label: String,
    total_secs: f64,
    total_bill: usize,
    misses: usize,
    rows: Vec<(String, Row)>,
}

fn run_suite(
    cases: &[(&str, &str, &str)],
    config: FindPathConfig,
    repeats: u32,
    budget: Duration,
) -> SuiteSummary {
    let label = config.heap_score.label().to_string();
    let mut total_secs = 0.0;
    let mut total_bill = 0;
    let mut misses = 0;
    let mut rows = Vec::new();
    let suite_t0 = Instant::now();
    for &(name, reactant, target) in cases {
        if suite_t0.elapsed() >= budget {
            misses += 1;
            rows.push((
                name.to_string(),
                Row {
                    hit: false,
                    hits: 0,
                    seconds: 0.0,
                    steps: 0,
                    nodes: 0,
                    mol_edits: 0,
                    billed: 0,
                },
            ));
            continue;
        }
        let row = run_one(reactant, target, config, repeats);
        total_secs += row.seconds;
        total_bill += row.billed;
        if !row.hit {
            misses += 1;
        }
        rows.push((name.to_string(), row));
    }
    SuiteSummary {
        label,
        total_secs,
        total_bill,
        misses,
        rows,
    }
}

fn print_suite_detail(title: &str, summary: &SuiteSummary, config: &FindPathConfig) {
    println!(
        "\n=== {title} (atom_diff={}, lazy_closer={}, max_paths={}, score={}) ===",
        config.use_atom_diff, config.lazy_closer, config.max_paths, summary.label
    );
    println!(
        "{:<32} {:>4} {:>5} {:>9} {:>5} {:>6} {:>7} {:>6}",
        "case", "hit", "hits", "seconds", "steps", "nodes", "edits", "bill"
    );
    for (name, row) in &summary.rows {
        if row.hits == 0 && row.seconds == 0.0 && !row.hit {
            println!("{name:<32} SKIP");
            continue;
        }
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
    println!(
        "{:<32} {:>4} {:>5} {:>9.3} miss={} bill={}",
        "TOTAL", "", "", summary.total_secs, summary.misses, summary.total_bill
    );
}

fn print_matrix_summary(summaries: &[SuiteSummary]) {
    println!("\n=== score matrix summary ===");
    println!(
        "{:<18} {:>9} {:>8} {:>5}",
        "score", "seconds", "bill", "miss"
    );
    let mut ranked: Vec<&SuiteSummary> = summaries.iter().collect();
    ranked.sort_by(|a, b| {
        a.misses
            .cmp(&b.misses)
            .then_with(|| a.total_bill.cmp(&b.total_bill))
            .then_with(|| {
                a.total_secs
                    .partial_cmp(&b.total_secs)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    for s in ranked {
        println!(
            "{:<18} {:>9.3} {:>8} {:>5}",
            s.label, s.total_secs, s.total_bill, s.misses
        );
    }
}

fn print_table(
    title: &str,
    cases: &[(&str, &str, &str)],
    config: FindPathConfig,
    repeats: u32,
    budget: Duration,
) {
    let summary = run_suite(cases, config, repeats, budget);
    print_suite_detail(title, &summary, &config);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let larger = args.iter().any(|a| a == "--larger");
    let hard = args.iter().any(|a| a == "--hard");
    let matrix = args.iter().any(|a| a == "--matrix");
    // Default: filter-only. Unfiltered (slow) is --nofilter.
    let nofilter = args.iter().any(|a| a == "--nofilter");
    let filter_only = !nofilter || args.iter().any(|a| a == "--filter-only");
    let budget = parse_budget(&args, if nofilter { 60 } else { 30 });
    let max_paths = parse_paths(&args);
    let heap_score = parse_score(&args);
    let mcs_extend_on = args.iter().any(|a| a == "--mcs-extend");
    set_use_mcs_extend(mcs_extend_on);

    println!(
        "Rust find_path PhaseOne  max_nodes={MAX_NODES}  max_paths={max_paths}  best-of-{REPEATS}"
    );

    let (cases, title) = if hard {
        (HARD, "hard HA≈14–20 · multi-step (≥3–8 hops)")
    } else if larger {
        (LARGER, "larger HA≈17–26")
    } else {
        (CASES, "mid-size / multi-edit")
    };

    if matrix {
        println!(
            "(release; filter-only={filter_only}; mcs_extend={}; score matrix 3×3 + soft; budget={}s)",
            use_mcs_extend(),
            budget.as_secs()
        );
        let mut modes: Vec<HeapScoreMode> = MatchScoreSpec::matrix()
            .into_iter()
            .map(HeapScoreMode::Match)
            .collect();
        modes.push(HeapScoreMode::SoftStack);
        let mut summaries = Vec::new();
        for mode in modes {
            let config = FindPathConfig {
                max_paths,
                max_nodes: MAX_NODES,
                use_atom_diff: true,
                lazy_closer: true,
                heap_score: mode,
                ..FindPathConfig::default()
            };
            let summary = run_suite(cases, config, REPEATS, budget);
            print_suite_detail(title, &summary, &config);
            summaries.push(summary);
        }
        print_matrix_summary(&summaries);
        return;
    }

    println!(
        "(release; tagged ForestMol; filter-only={filter_only}; heap_score={}; mcs_extend={}; budget={}s)",
        heap_score.label(),
        use_mcs_extend(),
        budget.as_secs()
    );

    let base = FindPathConfig {
        max_paths,
        max_nodes: MAX_NODES,
        heap_score,
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
