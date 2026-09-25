//! Phase profile for the **tagged `ForestMol`** `find_path` door.
//!
//! Same walk as production: structure/`csmi` cache, `Atom.tag` through
//! `adopt_product`, lazy cost closer on pop, eager tag-lift child MCS.
//! Does **not** re-parse walk CSMI strings.
//!
//! ```text
//! cargo run -p xenosite-forest --example find_path_profile --release
//! ```

use std::collections::{BinaryHeap, HashSet};
use std::rc::Rc;
use std::time::{Duration, Instant};

use xenosite_forest::atom_diff::{
    atom_diff, atom_diff_for_child, candidate_could_help_on, candidate_order_key, pair_could_help,
};
use xenosite_forest::forest_mol::ForestMol;
use xenosite_forest::mol::{Molecule, canon_of, parse_mol};
use xenosite_forest::rules::phase_one;
use xenosite_forest::ruleset::RuleSet;
use xenosite_forest::{Candidate, FindPathConfig, HeapScoreMode, PathCounters, find_path_with};

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
    csmi: Duration,
    atom_diff: Duration,
    lift: Duration,
    discover: Duration,
    filter: Duration,
    materialize: Duration,
    reject: Duration,
    enqueue: Duration,
    nodes: usize,
    rejected: usize,
    mol_edits: usize,
    cand_kept: usize,
    pair_kept: usize,
    csmi_calls: usize,
}

impl Timers {
    fn total_accounted(&self) -> Duration {
        self.csmi
            + self.atom_diff
            + self.lift
            + self.discover
            + self.filter
            + self.materialize
            + self.reject
            + self.enqueue
    }
}

struct Walk {
    mol: ForestMol,
    steps_len: usize,
    parent_cost: Option<usize>,
}

struct HeapItem {
    target_hit: bool,
    seq: usize,
    walk: Walk,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.target_hit == other.target_hit && self.seq == other.seq
    }
}
impl Eq for HeapItem {}
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

struct ForestEmission {
    products: Vec<ForestMol>,
}

fn keep_fragment(
    products: &[ForestMol],
    target_csmi: &str,
    target_ha: usize,
) -> Option<(ForestMol, String)> {
    let mut best: Option<(usize, Rc<str>, usize)> = None;
    for (i, mol) in products.iter().enumerate() {
        let csmi = mol.csmi();
        let cost = if csmi.as_ref() == target_csmi {
            0
        } else {
            1 + mol.heavy_atom_count().abs_diff(target_ha)
        };
        match &best {
            None => best = Some((i, Rc::clone(&csmi), cost)),
            Some((_, _, bc)) if cost < *bc => best = Some((i, Rc::clone(&csmi), cost)),
            Some((_, kept, bc)) if cost == *bc && csmi.as_ref() < kept.as_ref() => {
                best = Some((i, Rc::clone(&csmi), cost));
            }
            _ => {}
        }
    }
    let (i, csmi, _) = best?;
    Some((products[i].clone(), csmi.as_ref().to_string()))
}

fn expand_timed(
    ruleset: &RuleSet,
    parent: &ForestMol,
    target: &Molecule,
    diff: &xenosite_forest::AtomDiff,
    t: &mut Timers,
) -> Vec<ForestEmission> {
    let mol = parent.mol();
    let t0 = Instant::now();
    let mut candidates = ruleset
        .candidates(mol)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    // Nested sets: full pair_candidates (same as production walk).
    let mut pairs = ruleset
        .pair_candidates(mol)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    t.discover += t0.elapsed();

    let t0 = Instant::now();
    candidates.retain(|c| candidate_could_help_on(c, diff, Some(mol), Some(target)));
    candidates.sort_by_key(|c| candidate_order_key(c, diff));
    pairs.retain(|p| pair_could_help(p, diff, mol, target));
    t.filter += t0.elapsed();
    t.cand_kept += candidates.len();
    t.pair_kept += pairs.len();

    let t0 = Instant::now();
    let mut out = Vec::new();
    for c in candidates {
        t.mol_edits += 1;
        let pieces = c.materialize_mols(mol).unwrap();
        if pieces.is_empty() {
            continue;
        }
        let products: Vec<_> = pieces
            .into_iter()
            .map(|piece| parent.adopt_product(piece))
            .collect();
        out.push(ForestEmission { products });
    }
    for pair in pairs {
        let pieces = pair.materialize_mols(mol).unwrap();
        if pieces.is_empty() {
            continue;
        }
        t.mol_edits += 1;
        let products: Vec<_> = pieces
            .into_iter()
            .map(|piece| parent.adopt_product(piece))
            .collect();
        let _ = ruleset.canonical_plan(
            mol,
            &pair.plan_site_atoms(),
            Some(&[&pair.left.effect, &pair.right.effect]),
        );
        out.push(ForestEmission { products });
    }
    t.materialize += t0.elapsed();
    out
}

/// Instrumented copy of production tagged-walk search (lazy or eager closer).
fn profile_search(reactant: &str, target: &str, lazy: bool) -> (Timers, Duration, bool, usize) {
    let wall0 = Instant::now();
    let mut t = Timers::default();
    let set = phase_one();

    let start = ForestMol::parse(reactant).unwrap();
    let t0 = Instant::now();
    let start_csmi = start.csmi();
    t.csmi += t0.elapsed();
    t.csmi_calls += 1;

    let target_csmi = canon_of(target).unwrap();
    let target_mol = parse_mol(&target_csmi).unwrap();
    let target_ha = ForestMol::parse(&target_csmi).unwrap().heavy_atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi.as_ref() == target_csmi.as_str(),
        seq,
        walk: Walk {
            mol: start,
            steps_len: 0,
            parent_cost: None,
        },
    });
    seq += 1;
    let mut seen = HashSet::new();
    seen.insert(start_csmi.as_ref().to_string());
    let mut hit = false;
    let mut found_steps = 0usize;
    const MAX_NODES: usize = 800;

    while let Some(item) = heap.pop() {
        if t.nodes >= MAX_NODES || hit {
            break;
        }
        let walk = item.walk;

        let t0 = Instant::now();
        let here = walk.mol.csmi();
        t.csmi += t0.elapsed();
        t.csmi_calls += 1;

        if here.as_ref() == target_csmi.as_str() {
            t.nodes += 1;
            hit = true;
            found_steps = walk.steps_len;
            break;
        }

        let t0 = Instant::now();
        let diff = atom_diff(walk.mol.mol(), &target_mol);
        let cost = diff.cost();
        t.atom_diff += t0.elapsed();

        if lazy {
            if let Some(pc) = walk.parent_cost {
                if cost >= pc {
                    t.rejected += 1;
                    t.reject += Duration::ZERO; // MCS already in atom_diff
                    continue;
                }
            }
        }
        t.nodes += 1;

        let emissions = expand_timed(&set, &walk.mol, &target_mol, &diff, &mut t);

        for emission in emissions {
            let t0 = Instant::now();
            let Some((kept, kept_csmi)) =
                keep_fragment(&emission.products, &target_csmi, target_ha)
            else {
                t.enqueue += t0.elapsed();
                continue;
            };
            // keep_fragment called csmi on products
            t.csmi_calls += emission.products.len();

            let target_hit = kept_csmi == target_csmi;
            let allow = if lazy {
                true
            } else if let Some(pc) = walk.parent_cost {
                let t1 = Instant::now();
                let child_cost = atom_diff_for_child(&walk.mol, &diff, &kept, &target_mol).cost();
                t.lift += t1.elapsed();
                target_hit || child_cost < pc
            } else {
                true
            };
            if !allow {
                t.enqueue += t0.elapsed();
                continue;
            }
            if seen.contains(&kept_csmi) && !target_hit {
                t.enqueue += t0.elapsed();
                continue;
            }
            seen.insert(kept_csmi);
            heap.push(HeapItem {
                target_hit,
                seq,
                walk: Walk {
                    mol: kept,
                    steps_len: walk.steps_len + 1,
                    parent_cost: Some(cost),
                },
            });
            seq += 1;
            t.enqueue += t0.elapsed();
            if target_hit {
                hit = true;
                found_steps = walk.steps_len + 1;
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

fn microbench() {
    println!("\n=== microbench (tagged ForestMol door) ===");
    let parent = ForestMol::parse("COc1cc(OC)c2c(OC)cc(OC)cc2c1").unwrap();
    let target = parse_mol("O=C1C=C(O)C(=O)c2c(O)cc(O)cc12").unwrap();
    let set = phase_one();

    let n = 2_000usize;
    let t0 = Instant::now();
    for _ in 0..n {
        let _ = parent.mol().clone();
    }
    let mol_clone = t0.elapsed();

    let t0 = Instant::now();
    let first = parent.csmi();
    let first_t = t0.elapsed();
    let t0 = Instant::now();
    for _ in 0..n {
        let _ = parent.csmi();
    }
    let cached = t0.elapsed();
    assert!(Rc::ptr_eq(&first, &parent.csmi()));

    let t0 = Instant::now();
    for _ in 0..100 {
        let _ = atom_diff(parent.mol(), &target);
    }
    let diff = t0.elapsed();

    // Adopt cost: one real hydroxylation-style product if any candidate exists.
    let cands = set
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let adopt = if let Some(c) = cands.first() {
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        if let Some(piece) = pieces.into_iter().next() {
            let t0 = Instant::now();
            for _ in 0..50 {
                let _ = parent.adopt_product(piece.clone());
            }
            Some(t0.elapsed())
        } else {
            None
        }
    } else {
        None
    };

    println!(
        "Molecule.clone           {n}×  {:.3} µs/op",
        mol_clone.as_secs_f64() * 1e6 / n as f64
    );
    println!(
        "ForestMol.csmi (cold)          {:.3} µs",
        first_t.as_secs_f64() * 1e6
    );
    println!(
        "ForestMol.csmi (cached)  {n}×  {:.3} ns/op",
        cached.as_secs_f64() * 1e9 / n as f64
    );
    println!(
        "atom_diff                100×  {:.3} ms/op",
        diff.as_secs_f64() * 1e3 / 100.0
    );
    if let Some(a) = adopt {
        println!(
            "adopt_product            50×  {:.3} µs/op",
            a.as_secs_f64() * 1e6 / 50.0
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let budget_secs: u64 = args
        .windows(2)
        .find(|w| w[0] == "--budget-secs")
        .and_then(|w| w[1].parse().ok())
        .unwrap_or(45);
    let budget = Duration::from_secs(budget_secs);
    let profile_t0 = Instant::now();

    println!(
        "Rust find_path phase profile (tagged ForestMol, atom_diff + lazy closer; budget {budget_secs}s)\n"
    );

    // Warm production API once.
    {
        let mut c = PathCounters::default();
        let _ = find_path_with(
            HARD[0].1,
            HARD[0].2,
            &phase_one(),
            &mut c,
            FindPathConfig::default(),
            |_: &Candidate| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
    }

    println!(
        "{:<28} {:>7} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>5} {:>4}",
        "case", "wall", "csmi%", "diff%", "lift%", "disc%", "filt%", "mat%", "enq%", "hit", "rej"
    );

    for &(name, r, tgt) in HARD {
        if profile_t0.elapsed() >= budget {
            println!("{name:<28} SKIP  (budget)");
            continue;
        }
        let mut best: Option<(Timers, Duration, bool, usize)> = None;
        for _ in 0..3 {
            let row = profile_search(r, tgt, true);
            best = Some(match best {
                None => row,
                Some(prev) if row.1 < prev.1 => row,
                Some(prev) => prev,
            });
        }
        let (tm, wall, hit, steps) = best.unwrap();
        let w = wall;
        println!(
            "{:<28} {:>6.3}s {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>5.1} {:>3}/{steps} {:>4}",
            name,
            w.as_secs_f64(),
            pct(tm.csmi, w),
            pct(tm.atom_diff, w),
            pct(tm.lift, w),
            pct(tm.discover, w),
            pct(tm.filter, w),
            pct(tm.materialize, w),
            pct(tm.enqueue, w),
            if hit { "ok" } else { "NO" },
            tm.rejected,
        );
        println!(
            "  nodes={} edits={} cand_kept={} pair_kept={} csmi_calls={} accounted={:.1}%",
            tm.nodes,
            tm.mol_edits,
            tm.cand_kept,
            tm.pair_kept,
            tm.csmi_calls,
            pct(tm.total_accounted(), w),
        );
    }

    if profile_t0.elapsed() < budget {
        microbench();
    }

    println!("\n=== production wall lazy vs eager (best of 3) ===");
    let set = phase_one();
    for &(name, r, tgt) in HARD {
        if profile_t0.elapsed() >= budget {
            println!("{name:<28} SKIP  (budget)");
            continue;
        }
        for (label, lazy) in [("lazy", true), ("eager", false)] {
            let config = FindPathConfig {
                max_paths: 1,
                max_nodes: 800,
                use_atom_diff: true,
                lazy_closer: lazy,
                heap_score: HeapScoreMode::MatchProduct,
            };
            let mut best = f64::MAX;
            let mut bill = 0usize;
            for _ in 0..3 {
                let mut c = PathCounters::default();
                let t0 = Instant::now();
                let _ = find_path_with(r, tgt, &set, &mut c, config, |_: &Candidate| true)
                    .unwrap()
                    .collect_all()
                    .unwrap();
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
