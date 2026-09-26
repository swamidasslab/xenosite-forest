//! Shell-cost gate ablation vs live `candidate_could_help` / `pair_could_help`.
//!
//! **Never** scores `|current − target|` alone (δ = 0). Site-shell cost is always
//! `|projected − target|` with `projected = current + editδ` (δ = product − reactant).
//!
//! Chains:
//! - **A** — legacy `at_sites` Σ|δ| baseline; optional jump to projected residual.
//! - **B** — projected residual (keep when residual drops vs no-edit); try +leave
//!   and +|Δaromatic| (aromatic bit mismatch 0/1 between projected and target).
//!
//! ```text
//! cargo run -p xenosite-forest --example shell_cost_bench --release
//! ```

use std::collections::BTreeSet;
use std::time::Instant;

use xenosite_forest::{
    FindPathConfig, ForestMol, PathCounters, SiteShellCostOpts, aligned_shells, atom_diff,
    candidate_could_help_on, edit_shells, find_path_with, molecule_shells, pair_could_help,
    parse_mol, phase_one, site_atoms_with_leave, site_shell_cost_opts,
};

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
        "COc1ccc(CC=C)cc1OC",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CostKind {
    /// Legacy: Σ|δ| on kept site atoms only (no cleaved count).
    AtSitesLegacy,
    /// Σ|δ| + cleaved + added counts.
    AtSitesCleaved,
    /// Σ|projected − target|; keep iff residual drops vs no-edit residual.
    /// No-edit residual is only the progress baseline — never a gate score.
    ProjResidualDrop,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Mode {
    leave: bool,
    /// Include |Δaromatic| ∈ {0,1} per atom when comparing projected vs target.
    aromatic_delta: bool,
    kind: CostKind,
}

impl Mode {
    fn at_sites_baseline() -> Self {
        Self {
            leave: false,
            aromatic_delta: false,
            kind: CostKind::AtSitesLegacy,
        }
    }

    fn proj_baseline() -> Self {
        Self {
            leave: false,
            aromatic_delta: false,
            kind: CostKind::ProjResidualDrop,
        }
    }

    fn label(self) -> String {
        let mut parts = Vec::new();
        parts.push(match self.kind {
            CostKind::AtSitesLegacy => "at_sites Σ|δ| (legacy)",
            CostKind::AtSitesCleaved => "at_sites Σ|δ|+cleaved",
            CostKind::ProjResidualDrop => "Σ|proj−tgt| keep-if-drop",
        });
        if self.aromatic_delta {
            parts.push("+|Δaromatic|");
        }
        if self.leave {
            parts.push("+leave");
        }
        parts.join(" ")
    }

    fn opts(self) -> SiteShellCostOpts {
        SiteShellCostOpts {
            dearomatic: self.aromatic_delta,
        }
    }
}

#[derive(Clone, Debug)]
struct SuiteStats {
    label: String,
    mid: GateStats,
    hard: GateStats,
}

#[derive(Clone, Debug, Default)]
struct GateStats {
    tp: usize,
    fp: usize,
    tn: usize,
    fn_: usize,
    top1: usize,
    known: usize,
}

impl GateStats {
    fn prec(&self) -> f64 {
        let d = self.tp + self.fp;
        if d == 0 {
            0.0
        } else {
            self.tp as f64 / d as f64
        }
    }
    fn rec(&self) -> f64 {
        let d = self.tp + self.fn_;
        if d == 0 {
            0.0
        } else {
            self.tp as f64 / d as f64
        }
    }
    fn agree(&self) -> f64 {
        let n = self.tp + self.fp + self.tn + self.fn_;
        if n == 0 {
            0.0
        } else {
            (self.tp + self.tn) as f64 / n as f64
        }
    }
}

#[derive(Clone, Debug)]
struct FnHit {
    suite: &'static str,
    case: String,
    kind: &'static str,
    name: String,
    before: f64,
    after: f64,
    atoms: Vec<usize>,
}

fn site_atoms_cand(c: &xenosite_forest::Candidate) -> Vec<usize> {
    let mut atoms: Vec<usize> = c
        .pattern
        .site_map
        .iter()
        .filter_map(|m| c.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.push(c.site);
    }
    for &i in &c.orbit {
        atoms.push(i);
    }
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn expand_atoms(
    mol: &xenosite_forest::Molecule,
    atoms: &[usize],
    leave: bool,
    leave_count: Option<usize>,
    cleavage_bonds: &BTreeSet<(usize, usize)>,
) -> Vec<usize> {
    if leave {
        site_atoms_with_leave(mol, atoms, leave_count, cleavage_bonds)
    } else {
        let mut a = atoms.to_vec();
        a.sort_unstable();
        a.dedup();
        a
    }
}

fn residual(
    mode: Mode,
    delta: Option<&xenosite_forest::AlignedShells>,
    cur: &xenosite_forest::MoleculeShells,
    tgt: &xenosite_forest::MoleculeShells,
    map: &std::collections::BTreeMap<usize, usize>,
    atoms: &[usize],
) -> f64 {
    site_shell_cost_opts(cur, delta, tgt, map, atoms, mode.opts())
}

fn residual_after_edit(
    mode: Mode,
    parent: &ForestMol,
    child: &ForestMol,
    cur: &xenosite_forest::MoleculeShells,
    tgt: &xenosite_forest::MoleculeShells,
    map: &std::collections::BTreeMap<usize, usize>,
    atoms: &[usize],
) -> f64 {
    let edit = edit_shells(parent, child);
    residual(mode, Some(&edit), cur, tgt, map, atoms)
}

fn site_cost_legacy(mode: Mode, align: &xenosite_forest::AlignedShells, atoms: &[usize]) -> f64 {
    match mode.kind {
        CostKind::AtSitesLegacy => align
            .at_sites(atoms)
            .atoms
            .values()
            .map(|e| e.abs_delta())
            .sum::<usize>() as f64,
        CostKind::AtSitesCleaved => align.at_sites(atoms).cost() as f64,
        CostKind::ProjResidualDrop => 0.0, // unused
    }
}

fn eval_suite(
    suite: &'static str,
    cases: &[(&str, &str, &str)],
    mode: Mode,
    fns: &mut Vec<FnHit>,
    collect_fn: bool,
) -> GateStats {
    let mut stats = GateStats::default();
    let set = phase_one();
    let proj = matches!(mode.kind, CostKind::ProjResidualDrop);
    for &(name, reactant, target) in cases {
        let ra = parse_mol(reactant).unwrap();
        let parent = ForestMol::parse(reactant).unwrap();
        let rb = parse_mol(target).unwrap();
        let align = aligned_shells(&ra, &rb);
        let cur = molecule_shells(&ra);
        let tgt = molecule_shells(&rb);
        let ad = atom_diff(&ra, &rb);
        let map = ad.mapping.clone();

        let mut scored: Vec<(f64, bool, Vec<usize>)> = Vec::new();

        for c in set.candidates(&ra).collect::<Result<Vec<_>, _>>().unwrap() {
            let leave_n = c.pattern.effect.leave_count.map(|n| n as usize);
            let atoms = expand_atoms(
                &ra,
                &site_atoms_cand(&c),
                mode.leave && c.pattern.effect.cleaves,
                leave_n,
                &ad.cleavage_bonds,
            );
            let gate = candidate_could_help_on(&c, &ad, Some(&ra), Some(&rb));
            let (score, rank_cost) = if proj {
                let before = residual(mode, None, &cur, &tgt, &map, &atoms);
                let Ok(pieces) = c.materialize_mols(&ra) else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                let child = parent.adopt_product(pieces[0].clone());
                let after = residual_after_edit(mode, &parent, &child, &cur, &tgt, &map, &atoms);
                let keep = before > after + 1e-12;
                if collect_fn && gate && !keep {
                    fns.push(FnHit {
                        suite,
                        case: name.into(),
                        kind: "cand",
                        name: format!("{}:{}", c.leaf_rule().unwrap_or("?"), c.pattern.name),
                        before,
                        after,
                        atoms: atoms.clone(),
                    });
                }
                (if keep { 1.0 } else { 0.0 }, before - after)
            } else {
                let cost = site_cost_legacy(mode, &align, &atoms);
                (cost, cost)
            };
            scored.push((rank_cost, gate, atoms));
            let shell = score > 1e-12;
            match (gate, shell) {
                (true, true) => stats.tp += 1,
                (false, false) => stats.tn += 1,
                (false, true) => stats.fp += 1,
                (true, false) => stats.fn_ += 1,
            }
        }

        for p in set
            .pair_candidates(&ra)
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
        {
            let mut atoms = p.plan_site_atoms();
            atoms.sort_unstable();
            atoms.dedup();
            let gate = pair_could_help(&p, &ad, &ra, &rb);
            let (score, rank_cost) = if proj {
                let before = residual(mode, None, &cur, &tgt, &map, &atoms);
                let Ok(pieces) = p.materialize_mols(&ra) else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                let child = parent.adopt_product(pieces[0].clone());
                let after = residual_after_edit(mode, &parent, &child, &cur, &tgt, &map, &atoms);
                let keep = before > after + 1e-12;
                if collect_fn && gate && !keep {
                    fns.push(FnHit {
                        suite,
                        case: name.into(),
                        kind: "pair",
                        name: p.pattern_name.clone(),
                        before,
                        after,
                        atoms: atoms.clone(),
                    });
                }
                (if keep { 1.0 } else { 0.0 }, before - after)
            } else {
                let cost = site_cost_legacy(mode, &align, &atoms);
                (cost, cost)
            };
            scored.push((rank_cost, gate, atoms));
            let shell = score > 1e-12;
            match (gate, shell) {
                (true, true) => stats.tp += 1,
                (false, false) => stats.tn += 1,
                (false, true) => stats.fp += 1,
                (true, false) => stats.fn_ += 1,
            }
        }

        // hop0 rank among gate-kept. Higher rank_cost = more progress (proj) or cost.
        let mut counters = PathCounters::default();
        let config = FindPathConfig {
            max_nodes: 800,
            max_paths: 1,
            use_atom_diff: true,
            ..FindPathConfig::default()
        };
        let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
            .unwrap()
            .collect_all()
            .unwrap();
        if let Some(hit) = hits.first() {
            if let Some(step) = hit.steps.first() {
                let mut atoms = if step.site_orbit.is_empty() {
                    vec![step.site]
                } else {
                    step.site_orbit.clone()
                };
                atoms = expand_atoms(
                    &ra,
                    &atoms,
                    mode.leave && (!step.sides.is_empty() || ad.has_cleavage()),
                    None,
                    &ad.cleavage_bonds,
                );
                let hop_rank_val = if proj {
                    scored
                        .iter()
                        .find(|(_, g, a)| *g && a == &atoms)
                        .map(|(c, _, _)| *c)
                        .unwrap_or(0.0)
                } else {
                    site_cost_legacy(mode, &align, &atoms)
                };
                let mut gate_costs: Vec<f64> = scored
                    .iter()
                    .filter(|(_, g, _)| *g)
                    .map(|(c, _, _)| *c)
                    .collect();
                gate_costs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                if !gate_costs.is_empty() {
                    stats.known += 1;
                    let rank = gate_costs
                        .iter()
                        .filter(|&&c| c > hop_rank_val + 1e-9)
                        .count()
                        + 1;
                    if rank == 1 {
                        stats.top1 += 1;
                    }
                }
            }
        }
    }
    stats
}

fn run_mode(mode: Mode, collect_fn: bool) -> (SuiteStats, Vec<FnHit>) {
    let mut fns = Vec::new();
    let mid = eval_suite("mid", MID, mode, &mut fns, collect_fn);
    let hard = eval_suite("hard", HARD, mode, &mut fns, collect_fn);
    (
        SuiteStats {
            label: mode.label(),
            mid,
            hard,
        },
        fns,
    )
}

fn score_key(s: &SuiteStats) -> (f64, f64, f64) {
    let agree = s.mid.agree() + s.hard.agree();
    let prec = s.mid.prec() + s.hard.prec();
    let top = {
        let k = s.mid.known + s.hard.known;
        if k == 0 {
            0.0
        } else {
            (s.mid.top1 + s.hard.top1) as f64 / k as f64
        }
    };
    (agree, prec, top)
}

fn print_row(tag: &str, s: &SuiteStats) {
    println!(
        "{:<28}  mid {:.2}/{:.2}/{:.2} top{}/{}   hard {:.2}/{:.2}/{:.2} top{}/{}   [{}]",
        tag,
        s.mid.prec(),
        s.mid.rec(),
        s.mid.agree(),
        s.mid.top1,
        s.mid.known,
        s.hard.prec(),
        s.hard.rec(),
        s.hard.agree(),
        s.hard.top1,
        s.hard.known,
        s.label
    );
}

fn candidates_from(base: Mode) -> Vec<(String, Mode)> {
    let mut out = Vec::new();
    if !base.leave {
        out.push((
            "+leave".into(),
            Mode {
                leave: true,
                ..base
            },
        ));
    }
    match base.kind {
        CostKind::AtSitesLegacy => {
            out.push((
                "+cleaved count".into(),
                Mode {
                    kind: CostKind::AtSitesCleaved,
                    ..base
                },
            ));
            out.push((
                "→ Σ|proj−tgt| keep-if-drop".into(),
                Mode {
                    kind: CostKind::ProjResidualDrop,
                    aromatic_delta: false,
                    leave: base.leave,
                },
            ));
        }
        CostKind::AtSitesCleaved => {
            out.push((
                "→ Σ|proj−tgt| keep-if-drop".into(),
                Mode {
                    kind: CostKind::ProjResidualDrop,
                    aromatic_delta: false,
                    leave: base.leave,
                },
            ));
        }
        CostKind::ProjResidualDrop => {
            if !base.aromatic_delta {
                out.push((
                    "+|Δaromatic|".into(),
                    Mode {
                        aromatic_delta: true,
                        ..base
                    },
                ));
            }
        }
    }
    out
}

fn run_greedy_chain(
    start: Mode,
    start_tag: &str,
    history: &mut Vec<(String, SuiteStats)>,
) -> (Mode, SuiteStats) {
    let mut current = start;
    let (start_stats, _) = run_mode(current, false);
    print_row(start_tag, &start_stats);
    history.push((start_tag.into(), start_stats.clone()));
    let mut best = start_stats;
    let mut best_mode = current;

    println!("\n### Single adds from {start_tag}");
    for (tag, mode) in candidates_from(start) {
        if mode == start {
            continue;
        }
        let (s, _) = run_mode(mode, false);
        print_row(&tag, &s);
        history.push((format!("{start_tag} {tag}"), s));
    }

    println!("\n### Greedy keep-if-helps ({start_tag})");
    for _ in 0..6 {
        let mut round_best: Option<(String, Mode, SuiteStats)> = None;
        for (tag, mode) in candidates_from(current) {
            if mode == current {
                continue;
            }
            let (s, _) = run_mode(mode, false);
            print_row(&format!("  try {tag}"), &s);
            if score_key(&s) > score_key(&best) {
                match &round_best {
                    None => round_best = Some((tag, mode, s)),
                    Some((_, _, prev)) if score_key(&s) > score_key(prev) => {
                        round_best = Some((tag, mode, s));
                    }
                    _ => {}
                }
            }
        }
        if let Some((tag, mode, s)) = round_best {
            println!("  KEEP {tag}  → {}", s.label);
            history.push((format!("KEEP {tag}"), s.clone()));
            best = s;
            best_mode = mode;
            current = mode;
        } else {
            println!("  (no improving add)");
            break;
        }
    }
    (best_mode, best)
}

fn print_ledger(history: &[(String, SuiteStats)]) {
    println!("\n## Ledger (all measured)");
    println!(
        "{:<42} {:>14} {:>5}  {:>14} {:>5}  config",
        "step", "mid P/R/A", "top", "hard P/R/A", "top"
    );
    for (tag, s) in history {
        println!(
            "{:<42} {:>4.2}/{:.2}/{:.2}  {:>2}/{:<2}  {:>4.2}/{:.2}/{:.2}  {:>2}/{:<2}  {}",
            tag,
            s.mid.prec(),
            s.mid.rec(),
            s.mid.agree(),
            s.mid.top1,
            s.mid.known,
            s.hard.prec(),
            s.hard.rec(),
            s.hard.agree(),
            s.hard.top1,
            s.hard.known,
            s.label
        );
    }
}

fn print_fn_report(fns: &[FnHit]) {
    println!(
        "\n## False negatives (live gate=true, but |proj−tgt| did not drop)\n\
         These cut recall: the edit is allowed by atom_diff, yet site shells after\n\
         the applied edit are not closer to the target than before.\n"
    );
    if fns.is_empty() {
        println!("  (none)");
        return;
    }
    println!(
        "{:<5} {:<28} {:<6} {:<36} {:>8} {:>8} atoms",
        "suite", "case", "kind", "rule:pattern", "before", "after"
    );
    for h in fns {
        println!(
            "{:<5} {:<28} {:<6} {:<36} {:>8.3} {:>8.3} {:?}",
            h.suite, h.case, h.kind, h.name, h.before, h.after, h.atoms
        );
    }
    println!("\n  count: {}", fns.len());
}

fn main() {
    println!(
        "shell_cost_bench — greedy ablation vs live gate\n\
         Never uses |current−target| as a score (that was a bug).\n\
         Proj cost = |projected−target|, projected = current + editδ.\n\
         |Δaromatic| = optional 0/1 per-atom aromatic mismatch (was labeled “dear”).\n\
         columns: P/R/Agree  hop0-top1/known\n"
    );
    let t0 = Instant::now();
    let mut history: Vec<(String, SuiteStats)> = Vec::new();

    println!("## Chain A — from legacy at_sites baseline");
    let (best_a_mode, best_a) =
        run_greedy_chain(Mode::at_sites_baseline(), "A BASELINE", &mut history);

    println!("\n## Chain B — from Σ|proj−tgt| keep-if-drop (shells only)");
    let (best_b_mode, best_b) = run_greedy_chain(Mode::proj_baseline(), "B PROJ", &mut history);

    print_ledger(&history);

    // FN dump for the proj baseline (explains recall < 1).
    let (proj_stats, fns) = run_mode(Mode::proj_baseline(), true);
    println!(
        "\n## Recall check — B PROJ  mid R={:.2} (FN={})  hard R={:.2} (FN={})",
        proj_stats.mid.rec(),
        proj_stats.mid.fn_,
        proj_stats.hard.rec(),
        proj_stats.hard.fn_
    );
    print_fn_report(&fns);

    let (best_mode, best, which) = if score_key(&best_b) > score_key(&best_a) {
        (best_b_mode, best_b, "B")
    } else {
        (best_a_mode, best_a, "A")
    };
    println!(
        "\nBEST overall [{which}]: {}  (agree_sum={:.3} prec_sum={:.3})",
        best_mode.label(),
        best.mid.agree() + best.hard.agree(),
        best.mid.prec() + best.hard.prec()
    );
    println!("wall {:.2}s", t0.elapsed().as_secs_f64());
}
