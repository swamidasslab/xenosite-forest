//! Scan MetXBioDB Phase I pairs for hard `find_path` / `find_path_partial` cases.
//!
//! Input TSV from RDKit InChI→SMILES (see `artifacts/metx_phase1_pairs.tsv`):
//! ```text
//! uv run --extra rdkit python …  # regenerates TSV
//! cargo run -p xenosite-forest --example metx_hard_cases --release -- \
//!   artifacts/metx_phase1_pairs.tsv 200 1.5
//! ```
//! Args: `[tsv] [max_nodes=200] [timeout_secs=1.5] [limit=0]` (`limit` 0 = all).

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::{Duration, Instant};

use xenosite_forest::{
    FindPathConfig, PathCounters, atom_diff, canon_of, find_path_partial, parse_mol, phase_one,
    residual_from_diff,
};

#[derive(Clone)]
struct Pair {
    biot_id: String,
    substrate: String,
    product: String,
    reaction_type: String,
    reactant_smi: String,
    product_smi: String,
}

#[derive(Clone)]
struct RowOut {
    pair: Pair,
    hit: bool,
    secs: f64,
    nodes: usize,
    edits: usize,
    bill: usize,
    timed_out: bool,
    residual_cost: Option<usize>,
    residual_cats: Vec<String>,
    closest_smi: Option<String>,
    root_cost: usize,
    root_n_extra: usize,
    parse_ok: bool,
    skipped_equal: bool,
}

fn load_tsv(path: &Path) -> Vec<Pair> {
    let f = File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().expect("empty tsv").unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .unwrap_or_else(|| panic!("missing column {name} in {header}"))
    };
    let i_id = idx("biot_id");
    let i_sub = idx("substrate");
    let i_prod = idx("product");
    let i_rt = idx("reaction_type");
    let i_rs = idx("reactant_smi");
    let i_ps = idx("product_smi");
    let mut out = Vec::new();
    for line in lines {
        let line = line.unwrap();
        if line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < cols.len() {
            continue;
        }
        out.push(Pair {
            biot_id: c[i_id].to_string(),
            substrate: c[i_sub].to_string(),
            product: c[i_prod].to_string(),
            reaction_type: c[i_rt].to_string(),
            reactant_smi: c[i_rs].to_string(),
            product_smi: c[i_ps].to_string(),
        });
    }
    out
}

fn run_one(pair: &Pair, config: FindPathConfig) -> RowOut {
    let mut row = RowOut {
        pair: pair.clone(),
        hit: false,
        secs: 0.0,
        nodes: 0,
        edits: 0,
        bill: 0,
        timed_out: false,
        residual_cost: None,
        residual_cats: Vec::new(),
        closest_smi: None,
        root_cost: 0,
        root_n_extra: 0,
        parse_ok: false,
        skipped_equal: false,
    };
    let Ok(rc) = canon_of(&pair.reactant_smi) else {
        return row;
    };
    let Ok(pc) = canon_of(&pair.product_smi) else {
        return row;
    };
    row.parse_ok = true;
    if rc == pc {
        row.skipped_equal = true;
        return row;
    }
    let Ok(rmol) = parse_mol(&pair.reactant_smi) else {
        row.parse_ok = false;
        return row;
    };
    let Ok(pmol) = parse_mol(&pair.product_smi) else {
        row.parse_ok = false;
        return row;
    };
    let diff0 = atom_diff(&rmol, &pmol);
    row.root_cost = diff0.cost();
    row.root_n_extra = diff0.n_extra;
    let set = phase_one();
    let mut counters = PathCounters::default();
    let t0 = Instant::now();
    let out = match find_path_partial(
        pair.reactant_smi.as_str(),
        pair.product_smi.as_str(),
        &set,
        &mut counters,
        config,
        None,
        |_| true,
    ) {
        Ok(o) => o,
        Err(_) => {
            row.secs = t0.elapsed().as_secs_f64();
            row.nodes = counters.nodes;
            row.edits = counters.mol_edits;
            row.bill = counters.billed();
            row.timed_out = counters.timed_out;
            return row;
        }
    };
    row.secs = t0.elapsed().as_secs_f64();
    row.nodes = counters.nodes;
    row.edits = counters.mol_edits;
    row.bill = counters.billed();
    row.timed_out = counters.timed_out;
    row.hit = !out.exact.is_empty();
    if let Some(p) = out.partials.first() {
        row.residual_cost = Some(p.residual.cost);
        row.residual_cats = p.residual.categories.clone();
        row.closest_smi = Some(p.smiles.clone());
    } else if !row.hit {
        // No partial tracked — still report root residual categories.
        let r = residual_from_diff(&diff0, Some(&rmol), Some(&pmol));
        row.residual_cost = Some(r.cost);
        row.residual_cats = r.categories;
    }
    row
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let tsv = args
        .get(1)
        .map(Path::new)
        .unwrap_or(Path::new("artifacts/metx_phase1_pairs.tsv"));
    let max_nodes: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(200);
    let timeout_secs: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.5);
    let limit: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);

    let mut pairs = load_tsv(tsv);
    if limit > 0 {
        pairs.truncate(limit);
    }
    let config = FindPathConfig {
        max_paths: 1,
        max_nodes,
        use_atom_diff: true,
        lazy_closer: true,
        timeout: Some(Duration::from_secs_f64(timeout_secs)),
        ..FindPathConfig::default()
    };

    println!(
        "MetX Phase I hard-case scan  n={} max_nodes={max_nodes} timeout={timeout_secs}s  tsv={}",
        pairs.len(),
        tsv.display()
    );

    let mut results = Vec::new();
    let wall0 = Instant::now();
    for (i, pair) in pairs.iter().enumerate() {
        let row = run_one(pair, config);
        if (i + 1) % 50 == 0 || i + 1 == pairs.len() {
            let hits = results.iter().filter(|r: &&RowOut| r.hit).count()
                + usize::from(row.hit);
            let misses = results
                .iter()
                .filter(|r| r.parse_ok && !r.skipped_equal && !r.hit)
                .count()
                + usize::from(row.parse_ok && !row.skipped_equal && !row.hit);
            eprintln!(
                "[{}/{}] elapsed={:.1}s hits={hits} misses={misses} last_bill={} last_hit={}",
                i + 1,
                pairs.len(),
                wall0.elapsed().as_secs_f64(),
                row.bill,
                row.hit
            );
            let _ = std::io::stderr().flush();
        }
        results.push(row);
    }

    let parse_fail = results.iter().filter(|r| !r.parse_ok).count();
    let equal = results.iter().filter(|r| r.skipped_equal).count();
    let runnable: Vec<&RowOut> = results
        .iter()
        .filter(|r| r.parse_ok && !r.skipped_equal)
        .collect();
    let hits = runnable.iter().filter(|r| r.hit).count();
    let misses: Vec<&RowOut> = runnable.iter().copied().filter(|r| !r.hit).collect();
    let timed = runnable.iter().filter(|r| r.timed_out).count();

    println!(
        "\nsummary: runnable={} hits={} misses={} timed_out={} parse_fail={} equal_skip={} wall={:.1}s",
        runnable.len(),
        hits,
        misses.len(),
        timed,
        parse_fail,
        equal,
        wall0.elapsed().as_secs_f64()
    );

    // Hard misses: high bill, or high residual, or timeout.
    let mut hard = misses.clone();
    hard.sort_by(|a, b| {
        b.bill
            .cmp(&a.bill)
            .then(b.residual_cost.unwrap_or(0).cmp(&a.residual_cost.unwrap_or(0)))
            .then(b.secs.partial_cmp(&a.secs).unwrap())
    });

    println!("\n=== top 40 hard misses (by bill, then residual cost) ===");
    println!(
        "{:<14} {:>5} {:>6} {:>5} {:>5} {:>6} {:<4} {:>4} {:>3} {}",
        "biot_id", "bill", "secs", "nodes", "edits", "r_cost", "t/o", "root", "xtr", "name / reaction"
    );
    for r in hard.iter().take(40) {
        let rcost = r.residual_cost.unwrap_or(0);
        let name = if r.pair.substrate.is_empty() {
            r.pair.biot_id.clone()
        } else {
            format!("{} → {}", r.pair.substrate, r.pair.product)
        };
        println!(
            "{:<14} {:>5} {:>6.2} {:>5} {:>5} {:>6} {:<4} {:>4} {:>3} {} | {}",
            r.pair.biot_id,
            r.bill,
            r.secs,
            r.nodes,
            r.edits,
            rcost,
            if r.timed_out { "Y" } else { "" },
            r.root_cost,
            r.root_n_extra,
            name,
            r.pair.reaction_type
        );
        if !r.residual_cats.is_empty() {
            println!("               cats={:?} closest={:?}", r.residual_cats, r.closest_smi);
        }
    }

    // Also: high-bill hits (slow successes) — useful for multipath cost.
    let mut slow_hits: Vec<&RowOut> = runnable.iter().copied().filter(|r| r.hit).collect();
    slow_hits.sort_by(|a, b| b.bill.cmp(&a.bill).then(b.secs.partial_cmp(&a.secs).unwrap()));
    println!("\n=== top 15 expensive hits (by bill) ===");
    for r in slow_hits.iter().take(15) {
        println!(
            "{:<14} bill={:<5} secs={:.2} nodes={}  {} → {} | {}",
            r.pair.biot_id,
            r.bill,
            r.secs,
            r.nodes,
            r.pair.substrate,
            r.pair.product,
            r.pair.reaction_type
        );
    }

    // Dump full miss table for later.
    let dump = Path::new("artifacts/metx_hard_misses.tsv");
    let mut f = File::create(dump).unwrap();
    writeln!(
        f,
        "biot_id\tbill\tsecs\tnodes\tedits\ttimed_out\troot_cost\troot_n_extra\tresidual_cost\tcategories\tclosest\treactant_smi\tproduct_smi\tsubstrate\tproduct\treaction_type"
    )
    .unwrap();
    for r in &hard {
        writeln!(
            f,
            "{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.pair.biot_id,
            r.bill,
            r.secs,
            r.nodes,
            r.edits,
            r.timed_out,
            r.root_cost,
            r.root_n_extra,
            r.residual_cost.unwrap_or(0),
            r.residual_cats.join(","),
            r.closest_smi.as_deref().unwrap_or(""),
            r.pair.reactant_smi,
            r.pair.product_smi,
            r.pair.substrate.replace('\t', " "),
            r.pair.product.replace('\t', " "),
            r.pair.reaction_type.replace('\t', " "),
        )
        .unwrap();
    }
    println!("\nwrote {} ({} misses)", dump.display(), hard.len());
}
