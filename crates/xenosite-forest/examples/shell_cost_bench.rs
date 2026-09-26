//! Incremental shell-cost ablation from baseline.
//!
//! Baseline = `at_sites` raw Σ|δ| (aromatic + n0/n1/n2, missing=0), no leave
//! expand. Each step tries unused features one at a time; keeps the add that
//! best improves (agree, then precision, then hop0 top-1) on mid+hard vs the
//! live gate. Recall stays ~1.0 across trials.
//!
//! ```text
//! cargo run -p xenosite-forest --example shell_cost_bench --release
//! ```

use std::collections::BTreeSet;
use std::time::Instant;

use xenosite_forest::{
    FindPathConfig, PathCounters, SiteShellCostOpts, aligned_shells, atom_diff,
    candidate_could_help_on, find_path_with, molecule_shells, pair_could_help, parse_mol,
    phase_one, site_atoms_with_leave, site_shell_cost_opts,
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
    /// Legacy baseline: Σ|δ| on kept site atoms only (no cleaved count).
    AtSitesLegacy,
    /// Current `at_sites().cost()` = Σ|δ| + cleaved + added.
    AtSitesCleaved,
    /// Multiset Σ |current−target| via [`AtomNeighborhood::l1`] (raw site shell cost).
    MultisetRaw,
    /// Multiset Σ norm|current−target|; shells only.
    MultisetNormShells,
    /// Multiset Σ norm|current−target|; shells + |Δaromatic|.
    MultisetNormDear,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Mode {
    leave: bool,
    kind: CostKind,
}

impl Mode {
    fn baseline() -> Self {
        Self {
            leave: false,
            kind: CostKind::AtSitesLegacy,
        }
    }

    fn label(self) -> String {
        let mut parts = Vec::new();
        parts.push(match self.kind {
            CostKind::AtSitesLegacy => "at_sites Σ|δ| (legacy)",
            CostKind::AtSitesCleaved => "at_sites Σ|δ|+cleaved",
            CostKind::MultisetRaw => "site_shell Σ|cur−tgt|",
            CostKind::MultisetNormShells => "site_shell norm shells",
            CostKind::MultisetNormDear => "site_shell norm+dear",
        });
        if self.leave {
            parts.push("+leave");
        }
        parts.join(" ")
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

fn neighborhood_bag_raw_l1(
    a: &[xenosite_forest::AtomNeighborhood],
    b: &[xenosite_forest::AtomNeighborhood],
) -> f64 {
    let mut unused = b.to_vec();
    let mut cost = 0usize;
    for env in a {
        if let Some((i, _)) = unused.iter().enumerate().min_by_key(|(_, x)| env.l1(x)) {
            let other = unused.swap_remove(i);
            cost += env.l1(&other);
        } else {
            cost += env.abs_delta();
        }
    }
    for other in &unused {
        cost += other.abs_delta();
    }
    cost as f64
}

fn multiset_raw_cost(
    cur: &xenosite_forest::MoleculeShells,
    tgt: &xenosite_forest::MoleculeShells,
    map: &std::collections::BTreeMap<usize, usize>,
    atoms: &[usize],
) -> f64 {
    let mut projected = Vec::new();
    let mut target_envs = Vec::new();
    for &r in atoms {
        let Some(c) = cur.atoms.get(&r) else {
            continue;
        };
        if let Some(&t) = map.get(&r) {
            if let Some(te) = tgt.atoms.get(&t) {
                target_envs.push(te.clone());
            }
        }
        projected.push(c.clone());
    }
    neighborhood_bag_raw_l1(&projected, &target_envs)
}

fn site_cost(
    mode: Mode,
    align: &xenosite_forest::AlignedShells,
    cur: &xenosite_forest::MoleculeShells,
    tgt: &xenosite_forest::MoleculeShells,
    map: &std::collections::BTreeMap<usize, usize>,
    atoms: &[usize],
) -> f64 {
    match mode.kind {
        CostKind::AtSitesLegacy => align
            .at_sites(atoms)
            .atoms
            .values()
            .map(|e| e.abs_delta())
            .sum::<usize>() as f64,
        CostKind::AtSitesCleaved => align.at_sites(atoms).cost() as f64,
        CostKind::MultisetRaw => multiset_raw_cost(cur, tgt, map, atoms),
        CostKind::MultisetNormShells => site_shell_cost_opts(
            cur,
            None,
            tgt,
            map,
            atoms,
            SiteShellCostOpts { dearomatic: false },
        ),
        CostKind::MultisetNormDear => site_shell_cost_opts(
            cur,
            None,
            tgt,
            map,
            atoms,
            SiteShellCostOpts { dearomatic: true },
        ),
    }
}

fn eval_suite(cases: &[(&str, &str, &str)], mode: Mode) -> GateStats {
    let mut stats = GateStats::default();
    let set = phase_one();
    for &(_name, reactant, target) in cases {
        let ra = parse_mol(reactant).unwrap();
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
            let cost = site_cost(mode, &align, &cur, &tgt, &map, &atoms);
            scored.push((cost, gate, atoms));
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
            let cost = site_cost(mode, &align, &cur, &tgt, &map, &atoms);
            scored.push((cost, gate, atoms));
        }

        for &(cost, gate, _) in &scored {
            let shell = cost > 1e-12;
            match (gate, shell) {
                (true, true) => stats.tp += 1,
                (false, false) => stats.tn += 1,
                (false, true) => stats.fp += 1,
                (true, false) => stats.fn_ += 1,
            }
        }

        // hop0 rank among gate-kept (needs a path).
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
                let cost = site_cost(mode, &align, &cur, &tgt, &map, &atoms);
                let mut gate_costs: Vec<f64> = scored
                    .iter()
                    .filter(|(_, g, _)| *g)
                    .map(|(c, _, _)| *c)
                    .collect();
                gate_costs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                if !gate_costs.is_empty() {
                    stats.known += 1;
                    let rank = gate_costs.iter().filter(|&&c| c > cost + 1e-9).count() + 1;
                    if rank == 1 {
                        stats.top1 += 1;
                    }
                }
            }
        }
    }
    stats
}

fn run_mode(mode: Mode) -> SuiteStats {
    SuiteStats {
        label: mode.label(),
        mid: eval_suite(MID, mode),
        hard: eval_suite(HARD, mode),
    }
}

fn score_key(s: &SuiteStats) -> (f64, f64, f64) {
    // Higher better: combined agree, then precision, then hop0 top1 rate.
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
                kind: base.kind,
            },
        ));
    }
    match base.kind {
        CostKind::AtSitesLegacy => {
            out.push((
                "+cleaved count".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::AtSitesCleaved,
                },
            ));
            out.push((
                "→ site_shell Σ|cur−tgt|".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetRaw,
                },
            ));
            out.push((
                "→ site_shell norm shells".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
            out.push((
                "→ site_shell norm+dear".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormDear,
                },
            ));
        }
        CostKind::AtSitesCleaved => {
            out.push((
                "→ site_shell Σ|cur−tgt|".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetRaw,
                },
            ));
            out.push((
                "→ site_shell norm shells".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
            out.push((
                "→ site_shell norm+dear".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormDear,
                },
            ));
        }
        CostKind::MultisetRaw => {
            out.push((
                "+normalize shells".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
            out.push((
                "+normalize+dearomatic".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormDear,
                },
            ));
        }
        CostKind::MultisetNormShells => {
            out.push((
                "+dearomatic".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormDear,
                },
            ));
        }
        CostKind::MultisetNormDear => {
            out.push((
                "−dearomatic (shells only)".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
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
    let start_stats = run_mode(current);
    print_row(start_tag, &start_stats);
    history.push((start_tag.into(), start_stats.clone()));
    let mut best = start_stats;
    let mut best_mode = current;

    println!("\n### Single adds from {start_tag}");
    for (tag, mode) in candidates_from(start) {
        if mode == start {
            continue;
        }
        let s = run_mode(mode);
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
            let s = run_mode(mode);
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
        "{:<36}  {:>22}  {:>22}  config",
        "step", "mid P/R/A top", "hard P/R/A top"
    );
    for (tag, s) in history {
        println!(
            "{:<36}  {:>5.2}/{:.2}/{:.2} {:>2}/{:<2}  {:>5.2}/{:.2}/{:.2} {:>2}/{:<2}  {}",
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

fn main() {
    println!(
        "shell_cost_bench — greedy ablation vs live gate\n\
         Chain A: legacy at_sites Σ|δ|\n\
         Chain B: site_shell_cost Σ|current−target| (multiset raw)\n\
         columns: P/R/Agree  hop0-top1/known\n"
    );
    let t0 = Instant::now();
    let mut history: Vec<(String, SuiteStats)> = Vec::new();

    println!("## Chain A — from legacy at_sites baseline");
    let (best_a_mode, best_a) = run_greedy_chain(Mode::baseline(), "A BASELINE", &mut history);

    println!("\n## Chain B — from site_shell_cost Σ|current−target|");
    let site_start = Mode {
        leave: false,
        kind: CostKind::MultisetRaw,
    };
    let (best_b_mode, best_b) = run_greedy_chain(site_start, "B SITE_SHELL", &mut history);

    print_ledger(&history);

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
