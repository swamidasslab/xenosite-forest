//! Phase + clone cost profile for hard `find_path` cases.
//!
//! Matches production lazy closer: HA gate at enqueue, full cost on pop.
//!
//! ```text
//! cargo run -p xenosite-forest --example find_path_profile --release
//! ```

use std::collections::{BinaryHeap, HashSet};
use std::time::{Duration, Instant};

use xenosite_forest::atom_diff::{
    atom_diff, candidate_could_help_on, candidate_order_key, pair_could_help,
};
use xenosite_forest::mol::{Molecule, canon_of, canon_smiles, parse_mol};
use xenosite_forest::pattern::Emission;
use xenosite_forest::rules::phase_one;
use xenosite_forest::ruleset::RuleSet;
use xenosite_forest::{FindPathConfig, PathCounters, find_path_with};

const HARD: &[(&str, &str, &str)] = &[
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
    (
        "trimethoxy-PEA→catechol",
        "COc1cc(OC)c(OC)c(CCN)c1",
        "NCCc1cc(O)c(O)c(O)c1",
    ),
];

#[derive(Default, Clone)]
struct Timers {
    parse: Duration,
    atom_diff: Duration,
    discover: Duration,
    filter: Duration,
    materialize: Duration,
    /// Pop-time cost verify rejects (paid MCS, no expand).
    reject: Duration,
    enqueue: Duration,
    nodes: usize,
    rejected: usize,
    mol_edits: usize,
    cand_kept: usize,
    pair_kept: usize,
}

impl Timers {
    fn total_accounted(&self) -> Duration {
        self.parse
            + self.atom_diff
            + self.discover
            + self.filter
            + self.materialize
            + self.reject
            + self.enqueue
    }
}

#[derive(Clone, Eq, PartialEq)]
struct HeapItem {
    target_hit: bool,
    seq: usize,
    smiles: String,
    heavy: usize,
    steps_len: usize,
    parent_cost: Option<usize>,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .target_hit
            .cmp(&self.target_hit)
            .then_with(|| other.seq.cmp(&self.seq))
    }
}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

fn heavy_atoms(s: &str) -> usize {
    parse_mol(s).map(|m| m.atom_count()).unwrap_or(0)
}

fn expand_timed(
    ruleset: &RuleSet,
    mol: &Molecule,
    target: &Molecule,
    diff: &xenosite_forest::AtomDiff,
    t: &mut Timers,
) -> Vec<Emission> {
    let t0 = Instant::now();
    let mut candidates = ruleset.candidates(mol).unwrap();
    let all_pairs = ruleset.pair_candidates(mol).unwrap();
    t.discover += t0.elapsed();

    let t0 = Instant::now();
    candidates.retain(|c| candidate_could_help_on(c, diff, Some(mol), Some(target)));
    candidates.sort_by_key(|c| candidate_order_key(c, diff));
    let mut kept_pairs = all_pairs;
    kept_pairs.retain(|p| pair_could_help(p, diff, mol, target));
    t.filter += t0.elapsed();
    t.cand_kept += candidates.len();
    t.pair_kept += kept_pairs.len();

    let t0 = Instant::now();
    let mut out = Vec::new();
    for c in candidates {
        t.mol_edits += 1;
        if let Some(em) = c.emit(mol).unwrap() {
            out.push(em);
        }
    }
    for p in kept_pairs {
        if let Some(em) = p.emit(mol).unwrap() {
            t.mol_edits += 1;
            out.push(Emission {
                site: em.site,
                pattern_name: em.pattern_name,
                rule_path: vec![None],
                products: em.products,
                plan: Vec::new(),
            });
        }
    }
    t.materialize += t0.elapsed();
    out
}

fn profile_search(reactant: &str, target: &str) -> (Timers, Duration, bool, usize) {
    let wall0 = Instant::now();
    let mut t = Timers::default();
    let set = phase_one();
    let start = parse_mol(reactant).unwrap();
    let start_csmi = canon_smiles(&start);
    let target_csmi = canon_of(target).unwrap();
    let target_mol = parse_mol(&target_csmi).unwrap();
    let target_ha = target_mol.atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi == target_csmi,
        seq,
        smiles: start_csmi.clone(),
        heavy: start.atom_count(),
        steps_len: 0,
        parent_cost: None,
    });
    seq += 1;
    let mut seen = HashSet::new();
    seen.insert(start_csmi);
    let mut found_steps = 0usize;
    let mut hit = false;
    const MAX_NODES: usize = 800;

    while let Some(item) = heap.pop() {
        if t.nodes >= MAX_NODES || hit {
            break;
        }
        if item.smiles == target_csmi {
            t.nodes += 1;
            hit = true;
            found_steps = item.steps_len;
            break;
        }

        let t0 = Instant::now();
        let mol = parse_mol(&item.smiles).unwrap();
        t.parse += t0.elapsed();

        let t0 = Instant::now();
        let diff = atom_diff(&mol, &target_mol);
        let cost = diff.cost();
        t.atom_diff += t0.elapsed();

        if let Some(pc) = item.parent_cost {
            if cost >= pc {
                let t0 = Instant::now();
                t.rejected += 1;
                t.reject += t0.elapsed();
                continue;
            }
        }
        t.nodes += 1;

        let emissions = expand_timed(&set, &mol, &target_mol, &diff, &mut t);

        for emission in emissions {
            let t0 = Instant::now();
            let products = &emission.products;
            let kept = {
                let mut best: Option<(String, usize, i32)> = None;
                for p in products {
                    let ha = heavy_atoms(p);
                    let dist = (ha as i32 - target_ha as i32).abs();
                    let is_hit = p == &target_csmi;
                    if is_hit {
                        best = Some((p.clone(), ha, -1));
                        break;
                    }
                    match &best {
                        None => best = Some((p.clone(), ha, dist)),
                        Some((_, _, d)) if dist < *d => best = Some((p.clone(), ha, dist)),
                        _ => {}
                    }
                }
                best
            };
            let Some((kept, child_ha, _)) = kept else {
                t.enqueue += t0.elapsed();
                continue;
            };
            let target_hit = kept == target_csmi;
            // Lazy: no HA gate — cost verified on pop.
            if seen.contains(&kept) && !target_hit {
                t.enqueue += t0.elapsed();
                continue;
            }
            seen.insert(kept.clone());
            heap.push(HeapItem {
                target_hit,
                seq,
                smiles: kept,
                heavy: child_ha,
                steps_len: item.steps_len + 1,
                parent_cost: Some(cost),
            });
            seq += 1;
            t.enqueue += t0.elapsed();
            if target_hit {
                hit = true;
                found_steps = item.steps_len + 1;
                break;
            }
        }
    }
    (t, wall0.elapsed(), hit, found_steps)
}

fn pct(part: Duration, whole: Duration) -> f64 {
    if whole.as_secs_f64() == 0.0 {
        0.0
    } else {
        100.0 * part.as_secs_f64() / whole.as_secs_f64()
    }
}

fn microbench_clones() {
    println!("\n=== microbench: clone / parse / diff unit cost ===");
    let mol = parse_mol("COc1cc(OC)c2c(OC)cc(OC)cc2c1").unwrap();
    let target = parse_mol("O=C1C=C(O)C(=O)c2c(O)cc(O)cc12").unwrap();
    let n = 2_000usize;
    let t0 = Instant::now();
    for _ in 0..n {
        let _ = mol.clone();
    }
    let mol_clone = t0.elapsed();
    let t0 = Instant::now();
    for _ in 0..100 {
        let _ = atom_diff(&mol, &target);
    }
    let diff = t0.elapsed();
    println!(
        "Molecule.clone  {n}×  {:.3} µs/op",
        mol_clone.as_secs_f64() * 1e6 / n as f64
    );
    println!(
        "atom_diff       100×  {:.3} ms/op",
        diff.as_secs_f64() * 1e3 / 100.0
    );
}

fn main() {
    println!("Rust find_path phase profile (release, atom_diff + lazy closer)\n");

    {
        let mut c = PathCounters::default();
        let _ = find_path_with(
            HARD[0].1,
            HARD[0].2,
            &phase_one(),
            &mut c,
            FindPathConfig::default(),
            |_| true,
        );
    }

    println!(
        "{:<28} {:>7} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>5} {:>4}",
        "case", "wall", "parse%", "diff%", "disc%", "filt%", "mat%", "enq%", "hit", "rej"
    );

    for &(name, r, tgt) in HARD {
        let mut best: Option<(Timers, Duration, bool, usize)> = None;
        for _ in 0..3 {
            let row = profile_search(r, tgt);
            best = Some(match best {
                None => row,
                Some(prev) if row.1 < prev.1 => row,
                Some(prev) => prev,
            });
        }
        let (tm, wall, hit, steps) = best.unwrap();
        let w = wall;
        println!(
            "{:<28} {:>6.3}s {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>3}/{steps} {:>4}",
            name,
            w.as_secs_f64(),
            pct(tm.parse, w),
            pct(tm.atom_diff, w),
            pct(tm.discover, w),
            pct(tm.filter, w),
            pct(tm.materialize, w),
            pct(tm.enqueue, w),
            if hit { "ok" } else { "NO" },
            tm.rejected,
        );
        println!(
            "  nodes={} edits={} cand_kept={} pair_kept={} accounted={:.1}%",
            tm.nodes,
            tm.mol_edits,
            tm.cand_kept,
            tm.pair_kept,
            pct(tm.total_accounted(), w),
        );
    }

    microbench_clones();

    println!("\n=== production wall lazy vs eager (best of 3) ===");
    let set = phase_one();
    for &(name, r, tgt) in HARD {
        for (label, lazy) in [("lazy", true), ("eager", false)] {
            let config = FindPathConfig {
                max_paths: 1,
                max_nodes: 800,
                use_atom_diff: true,
                lazy_closer: lazy,
            };
            let mut best = f64::MAX;
            let mut bill = 0usize;
            for _ in 0..3 {
                let mut c = PathCounters::default();
                let t0 = Instant::now();
                let _ = find_path_with(r, tgt, &set, &mut c, config, |_| true).unwrap();
                let sec = t0.elapsed().as_secs_f64();
                if sec < best {
                    best = sec;
                    bill = c.billed();
                }
            }
            println!("{name:<28} {label:<5} {best:.3}s  bill={bill}");
        }
    }
}
