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
    /// Matches pre-SiteShellBag `shell_cost_bench` numbers.
    AtSitesLegacy,
    /// Current `at_sites().cost()` = Σ|δ| + cleaved + added.
    AtSitesCleaved,
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
            CostKind::MultisetNormShells => "multiset norm shells",
            CostKind::MultisetNormDear => "multiset norm+dear",
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
                "→ multiset norm shells".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
            out.push((
                "→ multiset norm+dear".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormDear,
                },
            ));
        }
        CostKind::AtSitesCleaved => {
            out.push((
                "→ multiset norm shells".into(),
                Mode {
                    leave: base.leave,
                    kind: CostKind::MultisetNormShells,
                },
            ));
            out.push((
                "→ multiset norm+dear".into(),
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

fn main() {
    println!(
        "shell_cost_bench — greedy ablation from baseline vs live gate\n\
         columns: suite  P/R/Agree  hop0-top1/known\n"
    );
    let t0 = Instant::now();

    let mut history: Vec<(String, SuiteStats)> = Vec::new();
    let mut current = Mode::baseline();
    let baseline = run_mode(current);
    print_row("BASELINE", &baseline);
    history.push(("BASELINE".into(), baseline.clone()));
    let mut best = baseline;
    let mut best_mode = current;

    // Also record every single-add from baseline for the ledger.
    println!("\n## Single adds from baseline");
    for (tag, mode) in candidates_from(Mode::baseline()) {
        let s = run_mode(mode);
        print_row(&tag, &s);
        history.push((tag, s));
    }

    println!("\n## Greedy keep-if-helps");
    loop {
        let mut improved = false;
        let mut round_best: Option<(String, Mode, SuiteStats)> = None;
        for (tag, mode) in candidates_from(current) {
            // Skip modes already equal to current.
            if mode == current {
                continue;
            }
            let s = run_mode(mode);
            print_row(&format!("  try {tag}"), &s);
            let better = score_key(&s) > score_key(&best);
            if better {
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
            improved = true;
        }
        if !improved {
            println!("  (no improving add)");
            break;
        }
        // Stop if no further unused feature directions.
        if candidates_from(current)
            .into_iter()
            .all(|(_, m)| m == current || history.iter().any(|(_, h)| h.label == m.label()))
        {
            // Still allow one more round of tries; break only when none help.
        }
        // Cap rounds.
        if history.len() > 12 {
            break;
        }
    }

    println!("\n## Ledger (all measured)");
    println!(
        "{:<28}  {:>22}  {:>22}  config",
        "step", "mid P/R/A top", "hard P/R/A top"
    );
    for (tag, s) in &history {
        println!(
            "{:<28}  {:>5.2}/{:.2}/{:.2} {:>2}/{:<2}  {:>5.2}/{:.2}/{:.2} {:>2}/{:<2}  {}",
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
    println!(
        "\nBEST: {}  (agree_sum={:.3} prec_sum={:.3})",
        best_mode.label(),
        best.mid.agree() + best.hard.agree(),
        best.mid.prec() + best.hard.prec()
    );
    println!("wall {:.2}s", t0.elapsed().as_secs_f64());
}
