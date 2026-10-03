//! Scan MetXBioDB Phase I pairs for hard `find_path` / `find_path_partial` cases.
//!
//! Input TSV from RDKit InChI→SMILES (see `artifacts/metx_phase1_pairs.tsv`):
//! ```text
//! cargo run -p xenosite-forest --example metx_hard_cases --release -- \
//!   artifacts/metx_phase1_pairs.tsv 200 1.5 0 1
//! ```
//! Args: `[tsv] [max_nodes=200] [timeout_secs=1.5] [limit=0] [normalize_tautomer=1]`
//! (`limit` 0 = all; `normalize_tautomer` defaults **on** for this scanner —
//! pass `0` to search given forms as-is).
//!
//! Reactant and product are stereo-stripped at the call site before search /
//! CSMI compare (enantiomers and E/Z match). No library canon change.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use chematic::chem::remove_stereo;
use rayon::prelude::*;
use xenosite_forest::{
    FindPathConfig, PathCounters, as_forest_mol, atom_diff, find_path_partial, phase_one,
    residual_from_diff,
};

/// Stereo-free CSMI for MetX matching (call-site only).
fn nostereo_smi(s: &str) -> Option<String> {
    let fm = as_forest_mol(s).ok()?;
    let stripped = remove_stereo(fm.mol());
    Some(as_forest_mol(stripped).ok()?.csmi().as_ref().to_string())
}

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
    /// Exact empty but a partial was flushed (closest residual).
    has_partial: bool,
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
        has_partial: false,
    };
    let Some(rc) = nostereo_smi(&pair.reactant_smi) else {
        return row;
    };
    let Some(pc) = nostereo_smi(&pair.product_smi) else {
        return row;
    };
    row.parse_ok = true;
    if rc == pc {
        row.skipped_equal = true;
        return row;
    }
    let Ok(rmol) = as_forest_mol(rc.as_str()) else {
        row.parse_ok = false;
        return row;
    };
    let Ok(pmol) = as_forest_mol(pc.as_str()) else {
        row.parse_ok = false;
        return row;
    };
    let diff0 = atom_diff(rmol.mol(), pmol.mol());
    row.root_cost = diff0.cost();
    row.root_n_extra = diff0.n_extra;
    let set = phase_one();
    let mut counters = PathCounters::default();
    let t0 = Instant::now();
    // Search on stereo-stripped CSMI spellings.
    let out = match find_path_partial(
        rc.as_str(),
        pc.as_str(),
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
    let (exact, partials): (Vec<_>, Vec<_>) = out.into_iter().partition(|o| o.is_exact());
    row.hit = !exact.is_empty();
    row.has_partial = !partials.is_empty();
    if let Some(p) = partials.first() {
        row.residual_cost = Some(p.residual.cost);
        row.residual_cats = p.residual.categories.clone();
        row.closest_smi = Some(p.smiles.clone());
    } else if !row.hit {
        // No partial tracked — still report root residual categories.
        let r = residual_from_diff(&diff0, Some(rmol.mol()), Some(pmol.mol()));
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
    let normalize_tautomer: bool = args
        .get(5)
        .and_then(|s| s.parse::<u8>().ok())
        .map(|v| v != 0)
        .unwrap_or(true);

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
        normalize_tautomer,
        ..FindPathConfig::default()
    };

    println!(
        "MetX Phase I hard-case scan  n={} max_nodes={max_nodes} timeout={timeout_secs}s normalize_tautomer={normalize_tautomer} nostereo=1 rayon  tsv={}",
        pairs.len(),
        tsv.display()
    );

    let wall0 = Instant::now();
    let done = AtomicUsize::new(0);
    let hit_count = AtomicUsize::new(0);
    let miss_count = AtomicUsize::new(0);
    let n = pairs.len();
    let progress = Mutex::new(());

    let results: Vec<RowOut> = pairs
        .par_iter()
        .map(|pair| {
            let row = run_one(pair, config);
            let i = done.fetch_add(1, Ordering::Relaxed) + 1;
            if row.parse_ok && !row.skipped_equal {
                if row.hit {
                    hit_count.fetch_add(1, Ordering::Relaxed);
                } else {
                    miss_count.fetch_add(1, Ordering::Relaxed);
                }
            }
            if i.is_multiple_of(50) || i == n {
                let _g = progress.lock().unwrap();
                eprintln!(
                    "[{i}/{n}] elapsed={:.1}s hits={} misses={} last_bill={} last_hit={}",
                    wall0.elapsed().as_secs_f64(),
                    hit_count.load(Ordering::Relaxed),
                    miss_count.load(Ordering::Relaxed),
                    row.bill,
                    row.hit
                );
                let _ = std::io::stderr().flush();
            }
            row
        })
        .collect();

    let parse_fail = results.iter().filter(|r| !r.parse_ok).count();
    let equal = results.iter().filter(|r| r.skipped_equal).count();
    let runnable: Vec<&RowOut> = results
        .iter()
        .filter(|r| r.parse_ok && !r.skipped_equal)
        .collect();
    let hits = runnable.iter().filter(|r| r.hit).count();
    let misses: Vec<&RowOut> = runnable.iter().copied().filter(|r| !r.hit).collect();
    let timed = runnable.iter().filter(|r| r.timed_out).count();
    let partial_only = runnable.iter().filter(|r| !r.hit && r.has_partial).count();
    let partial_zero = runnable
        .iter()
        .filter(|r| !r.hit && r.has_partial && r.residual_cost == Some(0))
        .count();

    println!(
        "\nsummary: runnable={} hits={} misses={} partial_flush={} partial_cost0={} timed_out={} parse_fail={} equal_skip={} normalize_tautomer={normalize_tautomer} wall={:.1}s",
        runnable.len(),
        hits,
        misses.len(),
        partial_only,
        partial_zero,
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
            .then(
                b.residual_cost
                    .unwrap_or(0)
                    .cmp(&a.residual_cost.unwrap_or(0)),
            )
            .then(b.secs.partial_cmp(&a.secs).unwrap())
    });

    println!("\n=== top 40 hard misses (by bill, then residual cost) ===");
    println!(
        "{:<14} {:>5} {:>6} {:>5} {:>5} {:>6} {:<4} {:>4} {:>3} name / reaction",
        "biot_id", "bill", "secs", "nodes", "edits", "r_cost", "t/o", "root", "xtr"
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
            println!(
                "               cats={:?} closest={:?}",
                r.residual_cats, r.closest_smi
            );
        }
    }

    // Also: high-bill hits (slow successes) — useful for multipath cost.
    let mut slow_hits: Vec<&RowOut> = runnable.iter().copied().filter(|r| r.hit).collect();
    slow_hits.sort_by(|a, b| {
        b.bill
            .cmp(&a.bill)
            .then(b.secs.partial_cmp(&a.secs).unwrap())
    });
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

    // Dump full miss table for later (nostereo call-site).
    let dump_name = if normalize_tautomer {
        "artifacts/metx_hard_misses_nostereo_tautnorm.tsv"
    } else {
        "artifacts/metx_hard_misses_nostereo.tsv"
    };
    let dump = Path::new(dump_name);
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

    // Thrash: high bill chasing dearomatize / saturate while still far from target.
    // Heuristic: miss + bill>=200 + (dearomatize|needs_oxygen) + residual_cost>=3.
    let thrash: Vec<&&RowOut> = hard
        .iter()
        .filter(|r| {
            r.bill >= 200
                && r.residual_cost.unwrap_or(0) >= 3
                && r.residual_cats.iter().any(|c| {
                    c == "dearomatize" || c == "needs_oxygen" || c == "extra_target_heavies"
                })
        })
        .collect();
    let thrash_name = if normalize_tautomer {
        "artifacts/metx_hard_thrash_nostereo_tautnorm.tsv"
    } else {
        "artifacts/metx_hard_thrash_nostereo.tsv"
    };
    let thrash_path = Path::new(thrash_name);
    let mut tf = File::create(thrash_path).unwrap();
    writeln!(
        tf,
        "biot_id\tbill\tsecs\troot_cost\tresidual_cost\tcategories\treaction_type\treactant_smi\tproduct_smi\tclosest"
    )
    .unwrap();
    for r in &thrash {
        writeln!(
            tf,
            "{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.pair.biot_id,
            r.bill,
            r.secs,
            r.root_cost,
            r.residual_cost.unwrap_or(0),
            r.residual_cats.join(","),
            r.pair.reaction_type.replace('\t', " "),
            r.pair.reactant_smi,
            r.pair.product_smi,
            r.closest_smi.as_deref().unwrap_or(""),
        )
        .unwrap();
    }
    println!(
        "wrote {} ({} thrash misses: bill>=200, rcost>=3, dearomatize/O/extra)",
        thrash_path.display(),
        thrash.len()
    );
    println!("\n=== top 20 thrash misses ===");
    for r in thrash.iter().take(20) {
        println!(
            "{:<14} bill={:<5} root={} rcost={} cats={:?} | {}",
            r.pair.biot_id,
            r.bill,
            r.root_cost,
            r.residual_cost.unwrap_or(0),
            r.residual_cats,
            r.pair.reaction_type
        );
    }
}
