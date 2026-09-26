//! Ablate find_path knobs one-at-a-time vs baseline.
//!
//! Baseline = defaults (eager closer, product-both, novel-site on, skeleton
//! twins on, atom_diff on, full PhaseOne). Each variant flips one lever.
//!
//! ```text
//! cargo run -p xenosite-forest --example find_path_ablate --release
//! cargo run -p xenosite-forest --example find_path_ablate --release -- --hard
//! cargo run -p xenosite-forest --example find_path_ablate --release -- --hard --paths 5
//! ```
//!
//! Flags: `--hard`, `--paths N` (default 1), `--budget-secs N` (default 60),
//! `--repeats N` (default 3), `--nofilter` (include unfiltered variant).

use std::time::{Duration, Instant};

use xenosite_forest::{
    FindPathConfig, HeapScoreMode, MatchCombine, MatchMetric, MatchScoreSpec, PathCounters,
    RuleSet, canon_of, find_path_with, leaf_rule, phase_one,
};

const MAX_NODES: usize = 800;

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

/// PhaseOne leaf order minus EpoxideHydration.
const PHASE_ONE_NO_EPOX_HYD: &[&str] = &[
    "Hydroxylation",
    "Epoxidation",
    "SulfurOxidation",
    "NitrogenOxidation",
    "Dehydrogenation",
    "QuinoneFormation",
    "Dephosphorylation",
    "EpoxideOpening",
    "Hydrolysis",
    "Dehydration",
    "Hydrogenation",
    "NitrogenReduction",
    "OxygenReduction",
    "ReductiveDehalogenation",
    "SulfurReduction",
    "Dealkylation",
    "OxidativeDehalogenation",
];

#[derive(Clone, Copy)]
enum RulesKind {
    PhaseOne,
    NoEpoxHyd,
}

#[derive(Clone, Copy)]
struct Variant {
    label: &'static str,
    rules: RulesKind,
    config: FindPathConfig,
}

#[derive(Clone, Debug)]
struct Row {
    hit: bool,
    hits: usize,
    seconds: f64,
    billed: usize,
    drop_dup: usize,
    drop_exact: usize,
    drop_skel: usize,
}

#[derive(Clone, Debug)]
struct Suite {
    label: String,
    misses: usize,
    total_bill: usize,
    total_secs: f64,
    total_drop_dup: usize,
    total_drop_exact: usize,
    total_drop_skel: usize,
    rows: Vec<(String, Row)>,
}

fn phase_one_without_epox_hyd() -> RuleSet {
    let members: Vec<RuleSet> = PHASE_ONE_NO_EPOX_HYD
        .iter()
        .map(|n| leaf_rule(n).unwrap_or_else(|| panic!("missing leaf {n}")))
        .collect();
    RuleSet::compose(Some("PhaseOneNoEpoxHyd".into()), members)
}

fn run_one(
    reactant: &str,
    target: &str,
    rules: &RuleSet,
    config: FindPathConfig,
    repeats: u32,
) -> Row {
    let want = canon_of(target).unwrap_or_else(|e| panic!("{target}: {e}"));
    {
        let mut c = PathCounters::default();
        let _ = find_path_with(reactant, target, rules, &mut c, config, |_| true)
            .unwrap()
            .collect_all()
            .unwrap();
    }
    let mut best: Option<Row> = None;
    for _ in 0..repeats {
        let mut counters = PathCounters::default();
        let t0 = Instant::now();
        let hits = find_path_with(reactant, target, rules, &mut counters, config, |_| true)
            .unwrap_or_else(|e| panic!("{reactant} → {target}: {e}"))
            .collect_all()
            .unwrap();
        let seconds = t0.elapsed().as_secs_f64();
        let hit = hits.first().is_some_and(|h| h.smiles == want);
        let row = Row {
            hit,
            hits: hits.len(),
            seconds,
            billed: counters.billed(),
            drop_dup: counters.dropped_duplicate_plan,
            drop_exact: counters.dropped_exact_plan,
            drop_skel: counters.dropped_skeleton_twin,
        };
        best = Some(match best {
            None => row,
            Some(prev) if row.seconds < prev.seconds => row,
            Some(prev) => prev,
        });
    }
    best.expect("repeats")
}

fn run_suite(
    cases: &[(&str, &str, &str)],
    variant: Variant,
    rules_cache: &(RuleSet, RuleSet),
    repeats: u32,
    budget: Duration,
) -> Suite {
    let rules = match variant.rules {
        RulesKind::PhaseOne => &rules_cache.0,
        RulesKind::NoEpoxHyd => &rules_cache.1,
    };
    let mut misses = 0;
    let mut total_bill = 0;
    let mut total_secs = 0.0;
    let mut total_drop_dup = 0;
    let mut total_drop_exact = 0;
    let mut total_drop_skel = 0;
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
                    billed: 0,
                    drop_dup: 0,
                    drop_exact: 0,
                    drop_skel: 0,
                },
            ));
            continue;
        }
        let row = run_one(reactant, target, rules, variant.config, repeats);
        total_secs += row.seconds;
        total_bill += row.billed;
        total_drop_dup += row.drop_dup;
        total_drop_exact += row.drop_exact;
        total_drop_skel += row.drop_skel;
        if !row.hit {
            misses += 1;
        }
        rows.push((name.to_string(), row));
    }
    Suite {
        label: variant.label.to_string(),
        misses,
        total_bill,
        total_secs,
        total_drop_dup,
        total_drop_exact,
        total_drop_skel,
        rows,
    }
}

fn print_suite(suite: &Suite, detail: bool) {
    if detail {
        println!("\n=== {} ===", suite.label);
        println!(
            "{:<32} {:>4} {:>5} {:>9} {:>6} {:>8} {:>8} {:>8}",
            "case", "hit", "hits", "seconds", "bill", "drop_dup", "exact", "skel"
        );
        for (name, row) in &suite.rows {
            if row.hits == 0 && row.seconds == 0.0 && !row.hit {
                println!("{name:<32} SKIP");
                continue;
            }
            println!(
                "{:<32} {:>4} {:>5} {:>9.3} {:>6} {:>8} {:>8} {:>8}",
                name,
                if row.hit { "ok" } else { "MISS" },
                row.hits,
                row.seconds,
                row.billed,
                row.drop_dup,
                row.drop_exact,
                row.drop_skel
            );
        }
    }
}

fn print_summary(suites: &[Suite], baseline: &Suite) {
    println!("\n=== ablation summary (Δbill / Δsecs vs baseline) ===");
    println!(
        "{:<18} {:>5} {:>8} {:>9} {:>8} {:>9} {:>8} {:>8} {:>8}",
        "variant", "miss", "bill", "Δbill", "seconds", "Δsecs", "drop_dup", "exact", "skel"
    );
    for s in suites {
        let db = s.total_bill as i64 - baseline.total_bill as i64;
        let ds = s.total_secs - baseline.total_secs;
        println!(
            "{:<18} {:>5} {:>8} {:>+8} {:>9.3} {:>+9.3} {:>8} {:>8} {:>8}",
            s.label,
            s.misses,
            s.total_bill,
            db,
            s.total_secs,
            ds,
            s.total_drop_dup,
            s.total_drop_exact,
            s.total_drop_skel
        );
    }
}

fn parse_u64(args: &[String], flag: &str, default: u64) -> u64 {
    args.windows(2)
        .find(|w| w[0] == flag)
        .and_then(|w| w[1].parse().ok())
        .unwrap_or(default)
}

fn parse_usize(args: &[String], flag: &str, default: usize) -> usize {
    args.windows(2)
        .find(|w| w[0] == flag)
        .and_then(|w| w[1].parse().ok())
        .unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let hard = args.iter().any(|a| a == "--hard");
    let nofilter = args.iter().any(|a| a == "--nofilter");
    let detail = args.iter().any(|a| a == "--detail");
    let max_paths = parse_usize(&args, "--paths", 1).max(1);
    let repeats = parse_u64(&args, "--repeats", 3) as u32;
    let budget = Duration::from_secs(parse_u64(&args, "--budget-secs", 60));

    let (cases, title) = if hard {
        (HARD, "hard")
    } else {
        (MID, "mid")
    };

    let base_cfg = FindPathConfig {
        max_paths,
        max_nodes: MAX_NODES,
        ..FindPathConfig::default()
    };

    let mut variants = vec![
        Variant {
            label: "baseline",
            rules: RulesKind::PhaseOne,
            config: base_cfg,
        },
        Variant {
            label: "lazy",
            rules: RulesKind::PhaseOne,
            config: FindPathConfig {
                lazy_closer: true,
                ..base_cfg
            },
        },
        Variant {
            label: "soft",
            rules: RulesKind::PhaseOne,
            config: FindPathConfig {
                heap_score: HeapScoreMode::SoftStack,
                ..base_cfg
            },
        },
        Variant {
            label: "close-atom",
            rules: RulesKind::PhaseOne,
            config: FindPathConfig {
                heap_score: HeapScoreMode::Match(MatchScoreSpec {
                    combine: MatchCombine::Close,
                    metric: MatchMetric::Atom,
                }),
                ..base_cfg
            },
        },
        Variant {
            label: "no-skeleton",
            rules: RulesKind::PhaseOne,
            config: FindPathConfig {
                drop_skeleton_twins: false,
                ..base_cfg
            },
        },
        Variant {
            label: "no-epox-hyd",
            rules: RulesKind::NoEpoxHyd,
            config: base_cfg,
        },
    ];
    if nofilter {
        variants.insert(
            1,
            Variant {
                label: "no-filter",
                rules: RulesKind::PhaseOne,
                config: FindPathConfig {
                    use_atom_diff: false,
                    lazy_closer: false,
                    ..base_cfg
                },
            },
        );
    }

    println!(
        "find_path_ablate  {title}  paths={max_paths}  nodes={MAX_NODES}  \
         best-of-{repeats}  budget={}s",
        budget.as_secs()
    );

    let rules_cache = (phase_one(), phase_one_without_epox_hyd());
    let mut suites = Vec::new();
    for v in variants {
        let suite = run_suite(cases, v, &rules_cache, repeats, budget);
        print_suite(&suite, detail);
        suites.push(suite);
    }
    let baseline = suites[0].clone();
    print_summary(&suites, &baseline);
}
