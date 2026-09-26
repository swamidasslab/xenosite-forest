//! Benchmark Rust `bfs` / `dfs` metabolite enumeration.
//!
//! Pair with `uv run python tests/forest/bench_enumerate.py`.
//!
//! ```text
//! cargo run -p xenosite-forest --example enumerate_bench --release
//! cargo run -p xenosite-forest --example enumerate_bench --release -- --phase-one
//! ```

use std::env;
use std::time::Instant;

use xenosite_forest::{
    EnumConfig, EnumOrder, RuleSet, enumerate_metabolites, hydroxylation, phase_one,
};

const REPEATS: usize = 5;

struct Case {
    name: &'static str,
    smiles: &'static str,
    depth: usize,
    order: EnumOrder,
}

const OH_CASES: &[Case] = &[
    Case {
        name: "OH ethane d1 bfs",
        smiles: "CC",
        depth: 1,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "OH ethane d1 dfs",
        smiles: "CC",
        depth: 1,
        order: EnumOrder::Dfs,
    },
    Case {
        name: "OH ethane d2 bfs",
        smiles: "CC",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "OH ethane d2 dfs",
        smiles: "CC",
        depth: 2,
        order: EnumOrder::Dfs,
    },
    Case {
        name: "OH benzene d2 bfs",
        smiles: "c1ccccc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "OH butylbenzene d2 bfs",
        smiles: "c1ccc(CCCC)cc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "OH butylbenzene d2 dfs",
        smiles: "c1ccc(CCCC)cc1",
        depth: 2,
        order: EnumOrder::Dfs,
    },
    Case {
        name: "OH toluene d2 bfs",
        smiles: "Cc1ccccc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
];

const PHASE_ONE_CASES: &[Case] = &[
    Case {
        name: "P1 ethane d2 bfs",
        smiles: "CC",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 anisole d1 bfs",
        smiles: "COc1ccccc1",
        depth: 1,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 anisole d2 bfs",
        smiles: "COc1ccccc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 anisole d2 dfs",
        smiles: "COc1ccccc1",
        depth: 2,
        order: EnumOrder::Dfs,
    },
    Case {
        name: "P1 eugenol d1 bfs",
        smiles: "COc1ccc(CC=C)cc1O",
        depth: 1,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 eugenol d2 bfs",
        smiles: "COc1ccc(CC=C)cc1O",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 veratrole d2 bfs",
        smiles: "COc1ccc(OC)cc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
    Case {
        name: "P1 phenacetin d2 bfs",
        smiles: "CCOc1ccc(NC(C)=O)cc1",
        depth: 2,
        order: EnumOrder::Bfs,
    },
];

fn run_one(set: &RuleSet, case: &Case) -> (usize, f64) {
    let config = EnumConfig {
        order: case.order,
        max_depth: case.depth,
        max_nodes: 0,
    };
    // Warmup
    let _ = enumerate_metabolites(case.smiles, set, config.clone())
        .unwrap()
        .count();
    let mut best = f64::INFINITY;
    let mut n = 0usize;
    for _ in 0..REPEATS {
        let t0 = Instant::now();
        n = enumerate_metabolites(case.smiles, set, config.clone())
            .unwrap()
            .map(|h| h.unwrap())
            .count();
        let secs = t0.elapsed().as_secs_f64();
        if secs < best {
            best = secs;
        }
    }
    (n, best)
}

fn bench_suite(title: &str, set: &RuleSet, cases: &[Case]) {
    println!("\n=== {title} (best-of-{REPEATS}) ===");
    println!(
        "{:<28} {:>6} {:>12} {:>8}",
        "case", "n", "seconds", "µs/hit"
    );
    let mut total_n = 0usize;
    let mut total_t = 0.0;
    for case in cases {
        let (n, secs) = run_one(set, case);
        let us_per = if n > 0 {
            secs * 1e6 / n as f64
        } else {
            0.0
        };
        println!(
            "{:<28} {:>6} {:>12.6} {:>8.1}",
            case.name, n, secs, us_per
        );
        total_n += n;
        total_t += secs;
    }
    println!(
        "{:<28} {:>6} {:>12.6}",
        "TOTAL", total_n, total_t
    );
}

fn main() {
    let phase_one_only = env::args().any(|a| a == "--phase-one");
    let oh_only = env::args().any(|a| a == "--oh");
    println!("Rust enumerate bfs/dfs  (unlimited nodes)");

    if !phase_one_only {
        bench_suite("Hydroxylation-only", &hydroxylation(), OH_CASES);
    }
    if !oh_only {
        bench_suite("PhaseOne", &phase_one(), PHASE_ONE_CASES);
    }
}
