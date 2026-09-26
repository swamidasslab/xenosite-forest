//! Shell-cost ablation: residual keep-if-drop is the decision under test.
//!
//! The live atom_diff gate is a **comparison**, not gold. Agreement with gate is
//! secondary; `path_miss` is when `find_path` hop0 did not drop residual.
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
    parse_mol, phase_one, site_atoms_with_leave, site_shell_cost, site_shell_cost_leave,
    site_shell_cost_opts,
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
    /// H bags one shell closer; drop n2 (on-demand view).
    h_closer_no_n2: bool,
    kind: CostKind,
}

impl Mode {
    fn at_sites_baseline() -> Self {
        Self {
            leave: false,
            aromatic_delta: false,
            h_closer_no_n2: false,
            kind: CostKind::AtSitesLegacy,
        }
    }

    fn proj_baseline() -> Self {
        Self {
            leave: false,
            aromatic_delta: false,
            h_closer_no_n2: false,
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
        if self.h_closer_no_n2 {
            parts.push("H→closer/no-n2");
        }
        if self.aromatic_delta {
            parts.push("+|Δaromatic|");
        }
        if self.leave {
            parts.push("+leave");
        }
        parts.join(" ")
    }

    fn opts_for_effect(self, dearomatizes: bool) -> SiteShellCostOpts {
        SiteShellCostOpts {
            dearomatic: self.aromatic_delta || dearomatizes,
            h_closer_no_n2: self.h_closer_no_n2,
        }
    }
}

#[derive(Clone, Debug)]
struct SuiteStats {
    label: String,
    mid: ShellStats,
    hard: ShellStats,
}

#[derive(Clone, Debug, Default)]
struct ShellStats {
    /// Residual drops (keep-if-drop decision under test).
    drop: usize,
    /// Residual does not drop.
    nodrop: usize,
    /// Live gate keeps.
    gate_keep: usize,
    /// Residual and gate agree (both keep or both refuse).
    agree: usize,
    /// find_path hop0 among residual-drop candidates is rank 1.
    top1: usize,
    /// Cases where find_path returned a hop0 we could score.
    known: usize,
    /// find_path hop0 where residual did **not** drop (true miss for shell gate).
    path_miss: usize,
}

impl ShellStats {
    fn drop_rate(&self) -> f64 {
        let n = self.drop + self.nodrop;
        if n == 0 {
            0.0
        } else {
            self.drop as f64 / n as f64
        }
    }
    fn agree_rate(&self) -> f64 {
        let n = self.drop + self.nodrop;
        if n == 0 {
            0.0
        } else {
            self.agree as f64 / n as f64
        }
    }
}

#[derive(Clone, Debug)]
struct DisHit {
    suite: &'static str,
    case: String,
    kind: &'static str,
    name: String,
    before: f64,
    after: f64,
    atoms: Vec<usize>,
}

fn site_atoms_cand(c: &xenosite_forest::Candidate) -> Vec<usize> {
    // PatternInfo.site_map only — same atoms as candidate_could_help / emit.
    // Unique-edit orbit is the collapsed class of primary-map hits, not the site.
    let mut atoms: Vec<usize> = c
        .pattern
        .site_map
        .iter()
        .filter_map(|m| c.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.push(c.site);
    }
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn expand_atoms(
    mol: &xenosite_forest::Molecule,
    atoms: &[usize],
    mapped: &[usize],
    leave: bool,
    leave_count: Option<usize>,
    cleavage_bonds: &BTreeSet<(usize, usize)>,
) -> Vec<usize> {
    if leave {
        site_atoms_with_leave(mol, atoms, leave_count, cleavage_bonds, mapped)
    } else {
        let mut a = atoms.to_vec();
        a.sort_unstable();
        a.dedup();
        a
    }
}

struct ResidArgs<'a> {
    atoms: &'a [usize],
    leave: &'a [usize],
    dearomatizes: bool,
}

fn residual(
    mode: Mode,
    delta: Option<&xenosite_forest::AlignedShells>,
    cur: &xenosite_forest::MoleculeShells,
    tgt: &xenosite_forest::MoleculeShells,
    map: &std::collections::BTreeMap<usize, usize>,
    args: ResidArgs<'_>,
) -> f64 {
    let opts = mode.opts_for_effect(args.dearomatizes);
    if args.leave.is_empty() {
        site_shell_cost_opts(cur, delta, tgt, map, args.atoms, opts)
    } else {
        site_shell_cost_leave(cur, delta, tgt, map, args.atoms, args.leave, opts)
    }
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

/// Leave debt for residual: sized leave = expand − site; open leave = bond ends
/// the applied edit actually drops (do not treat the kept end as leave).
fn leave_debt(
    site: &[usize],
    expanded: &[usize],
    leave_count: Option<usize>,
    edit: Option<&xenosite_forest::AlignedShells>,
) -> (Vec<usize>, Vec<usize>) {
    match leave_count {
        Some(_) => {
            let leave: Vec<usize> = expanded
                .iter()
                .copied()
                .filter(|a| !site.contains(a))
                .collect();
            (expanded.to_vec(), leave)
        }
        None => {
            let Some(edit) = edit else {
                // No edit yet: bond ends beyond site are provisional leave.
                let leave: Vec<usize> = expanded
                    .iter()
                    .copied()
                    .filter(|a| !site.contains(a))
                    .collect();
                return (expanded.to_vec(), leave);
            };
            let leave: Vec<usize> = expanded
                .iter()
                .copied()
                .filter(|a| !edit.alignment.contains_key(a))
                .collect();
            let mut atoms = site.to_vec();
            atoms.extend(&leave);
            atoms.sort_unstable();
            atoms.dedup();
            (atoms, leave)
        }
    }
}

fn record_decision(
    stats: &mut ShellStats,
    diss: &mut Vec<DisHit>,
    collect_dis: bool,
    drop: bool,
    gate: bool,
    meta: (&'static str, &str, &'static str, String, f64, f64, Vec<usize>),
) {
    let (suite, case, kind, name, before, after, atoms) = meta;
    if drop {
        stats.drop += 1;
    } else {
        stats.nodrop += 1;
    }
    if gate {
        stats.gate_keep += 1;
    }
    if drop == gate {
        stats.agree += 1;
    }
    // Disagreement dump: residual refused while gate kept (not a "false negative").
    if collect_dis && !drop && gate {
        diss.push(DisHit {
            suite,
            case: case.into(),
            kind,
            name,
            before,
            after,
            atoms,
        });
    }
}

fn eval_suite(
    suite: &'static str,
    cases: &[(&str, &str, &str)],
    mode: Mode,
    diss: &mut Vec<DisHit>,
    collect_dis: bool,
) -> ShellStats {
    let mut stats = ShellStats::default();
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

        // (rank_cost, dropped, atoms) — rank among residual-drop, not gate.
        let mut scored: Vec<(f64, bool, Vec<usize>)> = Vec::new();

        for c in set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
        {
            let leave_n = c.pattern.effect.leave_count.map(|n| n as usize);
            let site = site_atoms_cand(&c);
            let mapped: Vec<usize> = c.mapped.values().copied().collect();
            // Cleaving effects always include leave. Ablation `mode.leave` too.
            let expanded = expand_atoms(
                parent.mol(),
                &site,
                &mapped,
                c.pattern.effect.cleaves || mode.leave,
                leave_n,
                &ad.cleavage_bonds,
            );
            let gate = candidate_could_help_on(&c, &ad, Some(parent.mol()), Some(&rb));
            let dear = c.pattern.effect.dearomatizes;
            let (rank_cost, dropped, atoms) = if proj {
                let Ok(pieces) = c.materialize_mols(parent.mol()) else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                let child = parent.adopt_product(pieces[0].clone());
                let edit = edit_shells(&parent, &child);
                let (atoms, leave_only) = if c.pattern.effect.cleaves || mode.leave {
                    leave_debt(&site, &expanded, leave_n, Some(&edit))
                } else {
                    (expanded.clone(), Vec::new())
                };
                let before = residual(
                    mode,
                    None,
                    &cur,
                    &tgt,
                    &map,
                    ResidArgs {
                        atoms: &atoms,
                        leave: &leave_only,
                        dearomatizes: dear,
                    },
                );
                let after = residual(
                    mode,
                    Some(&edit),
                    &cur,
                    &tgt,
                    &map,
                    ResidArgs {
                        atoms: &atoms,
                        leave: &leave_only,
                        dearomatizes: dear,
                    },
                );
                let drop = before > after + 1e-12;
                record_decision(
                    &mut stats,
                    diss,
                    collect_dis,
                    drop,
                    gate,
                    (
                        suite,
                        name,
                        "cand",
                        format!("{}:{}", c.leaf_rule().unwrap_or("?"), c.pattern.name),
                        before,
                        after,
                        atoms.clone(),
                    ),
                );
                (before - after, drop, atoms)
            } else {
                let cost = site_cost_legacy(mode, &align, &expanded);
                let drop = cost > 1e-12;
                record_decision(
                    &mut stats,
                    diss,
                    collect_dis,
                    drop,
                    gate,
                    (
                        suite,
                        name,
                        "cand",
                        format!("{}:{}", c.leaf_rule().unwrap_or("?"), c.pattern.name),
                        cost,
                        0.0,
                        expanded.clone(),
                    ),
                );
                (cost, drop, expanded)
            };
            scored.push((rank_cost, dropped, atoms));
        }

        for p in set
            .pair_candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
        {
            let mut atoms = p.plan_site_atoms();
            atoms.sort_unstable();
            atoms.dedup();
            let gate = pair_could_help(&p, &ad, parent.mol(), &rb);
            let dear = p.effect.dearomatizes;
            let (rank_cost, dropped, atoms) = if proj {
                let Ok(pieces) = p.materialize_mols(parent.mol()) else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                let child = parent.adopt_product(pieces[0].clone());
                let edit = edit_shells(&parent, &child);
                let before = residual(
                    mode,
                    None,
                    &cur,
                    &tgt,
                    &map,
                    ResidArgs {
                        atoms: &atoms,
                        leave: &[],
                        dearomatizes: dear,
                    },
                );
                let after = residual(
                    mode,
                    Some(&edit),
                    &cur,
                    &tgt,
                    &map,
                    ResidArgs {
                        atoms: &atoms,
                        leave: &[],
                        dearomatizes: dear,
                    },
                );
                let drop = before > after + 1e-12;
                record_decision(
                    &mut stats,
                    diss,
                    collect_dis,
                    drop,
                    gate,
                    (
                        suite,
                        name,
                        "pair",
                        p.pattern_name.clone(),
                        before,
                        after,
                        atoms.clone(),
                    ),
                );
                (before - after, drop, atoms)
            } else {
                let cost = site_cost_legacy(mode, &align, &atoms);
                let drop = cost > 1e-12;
                record_decision(
                    &mut stats,
                    diss,
                    collect_dis,
                    drop,
                    gate,
                    (
                        suite,
                        name,
                        "pair",
                        p.pattern_name.clone(),
                        cost,
                        0.0,
                        atoms.clone(),
                    ),
                );
                (cost, drop, atoms)
            };
            scored.push((rank_cost, dropped, atoms));
        }

        // hop0: path_miss if residual did not drop; top1 among residual-drop.
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
                    &atoms,
                    mode.leave && (!step.sides.is_empty() || ad.has_cleavage()),
                    None,
                    &ad.cleavage_bonds,
                );
                let hop = scored.iter().find(|(_, _, a)| a == &atoms);
                let (hop_rank_val, hop_dropped) = if let Some((c, d, _)) = hop {
                    (*c, *d)
                } else if proj {
                    (0.0, false)
                } else {
                    (site_cost_legacy(mode, &align, &atoms), true)
                };
                if proj && !hop_dropped {
                    stats.path_miss += 1;
                }
                let mut drop_costs: Vec<f64> = scored
                    .iter()
                    .filter(|(_, d, _)| *d)
                    .map(|(c, _, _)| *c)
                    .collect();
                drop_costs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                if hop_dropped && !drop_costs.is_empty() {
                    stats.known += 1;
                    let rank = drop_costs
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

fn run_mode(mode: Mode, collect_dis: bool) -> (SuiteStats, Vec<DisHit>) {
    let mut diss = Vec::new();
    let mid = eval_suite("mid", MID, mode, &mut diss, collect_dis);
    let hard = eval_suite("hard", HARD, mode, &mut diss, collect_dis);
    (
        SuiteStats {
            label: mode.label(),
            mid,
            hard,
        },
        diss,
    )
}

/// Prefer residual keep-if-drop modes, then fewer path misses, then gate agree, then hop0 top1.
fn score_key(s: &SuiteStats) -> (bool, i64, f64, f64) {
    let proj = s.label.contains("keep-if-drop");
    let miss = -((s.mid.path_miss + s.hard.path_miss) as i64);
    let agree = s.mid.agree_rate() + s.hard.agree_rate();
    let top = {
        let k = s.mid.known + s.hard.known;
        if k == 0 {
            0.0
        } else {
            (s.mid.top1 + s.hard.top1) as f64 / k as f64
        }
    };
    (proj, miss, agree, top)
}

fn print_row(tag: &str, s: &SuiteStats) {
    println!(
        "{:<28}  mid drop{:.2} agree{:.2} miss{} top{}/{}   hard drop{:.2} agree{:.2} miss{} top{}/{}   [{}]",
        tag,
        s.mid.drop_rate(),
        s.mid.agree_rate(),
        s.mid.path_miss,
        s.mid.top1,
        s.mid.known,
        s.hard.drop_rate(),
        s.hard.agree_rate(),
        s.hard.path_miss,
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
                    h_closer_no_n2: false,
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
                    h_closer_no_n2: false,
                    leave: base.leave,
                },
            ));
        }
        CostKind::ProjResidualDrop => {
            if !base.h_closer_no_n2 {
                out.push((
                    "H→closer/no-n2".into(),
                    Mode {
                        h_closer_no_n2: true,
                        ..base
                    },
                ));
            }
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
        "{:<42} {:>16} {:>5} {:>4}  {:>16} {:>5} {:>4}  config",
        "step", "mid drop/agree", "top", "miss", "hard drop/agree", "top", "miss"
    );
    for (tag, s) in history {
        println!(
            "{:<42} {:>5.2}/{:<5.2}  {:>2}/{:<2} {:>4}  {:>5.2}/{:<5.2}  {:>2}/{:<2} {:>4}  {}",
            tag,
            s.mid.drop_rate(),
            s.mid.agree_rate(),
            s.mid.top1,
            s.mid.known,
            s.mid.path_miss,
            s.hard.drop_rate(),
            s.hard.agree_rate(),
            s.hard.top1,
            s.hard.known,
            s.hard.path_miss,
            s.label
        );
    }
}

fn print_dis_report(diss: &[DisHit], mid: &ShellStats, hard: &ShellStats) {
    println!(
        "\n## Gate-kept but residual nodrop (disagreement — gate is not gold)\n\
         path_miss counts find_path hop0 where residual did not drop.\n"
    );
    if diss.is_empty() {
        println!("  (none)");
    } else {
        println!(
            "{:<5} {:<28} {:<6} {:<36} {:>8} {:>8} atoms",
            "suite", "case", "kind", "rule:pattern", "before", "after"
        );
        for h in diss {
            println!(
                "{:<5} {:<28} {:<6} {:<36} {:>8.3} {:>8.3} {:?}",
                h.suite, h.case, h.kind, h.name, h.before, h.after, h.atoms
            );
        }
    }
    println!(
        "\n  disagreements={}  mid nodrop={} path_miss={}  hard nodrop={} path_miss={}",
        diss.len(),
        mid.nodrop,
        mid.path_miss,
        hard.nodrop,
        hard.path_miss
    );
}

fn all_heavy(shells: &xenosite_forest::MoleculeShells) -> Vec<usize> {
    let mut v: Vec<usize> = shells.atoms.keys().copied().collect();
    v.sort_unstable();
    v
}

fn fmt_shell_path(costs: &[f64]) -> String {
    costs
        .iter()
        .map(|c| format!("{c:.1}"))
        .collect::<Vec<_>>()
        .join("→")
}

fn fmt_atom_path(costs: &[usize]) -> String {
    costs
        .iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join("→")
}

/// Live `find_path` cost traces on mid/hard — atom_diff field cost and full-mol
/// shell residual (|current−target| at each node; display only, not a gate).
fn print_find_path_costs(label: &str, cases: &[(&str, &str, &str)]) {
    let set = phase_one();
    println!("\n## find_path costs — {label}");
    println!(
        "{:<28} {:>3} {:>4}  {:<28}  {:<18} mono",
        "case", "hit", "hops", "shell_path", "atom_path"
    );
    for &(name, reactant, target) in cases {
        let rb = parse_mol(target).unwrap();
        let tgt = molecule_shells(&rb);
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
        let Some(hit) = hits.first() else {
            println!("{name:<28} MISS {:>4}  —  —  —", 0);
            continue;
        };
        let mut shell_costs = Vec::new();
        let mut atom_costs = Vec::new();
        let mut shell_mono = true;
        let mut atom_mono = true;
        let mut cur = reactant.to_string();
        let mut prev_s: Option<f64> = None;
        let mut prev_a: Option<usize> = None;
        for step in &hit.steps {
            let mol = parse_mol(&cur).unwrap();
            let cur_s = molecule_shells(&mol);
            let ad = atom_diff(&mol, &rb);
            let sh = site_shell_cost(&cur_s, None, &tgt, &ad.mapping, &all_heavy(&cur_s));
            let at = ad.cost();
            if let Some(ps) = prev_s {
                if sh > ps + 1e-9 {
                    shell_mono = false;
                }
            }
            if let Some(pa) = prev_a {
                if at > pa {
                    atom_mono = false;
                }
            }
            prev_s = Some(sh);
            prev_a = Some(at);
            shell_costs.push(sh);
            atom_costs.push(at);
            cur = step.product.clone();
        }
        shell_costs.push(0.0);
        atom_costs.push(0);
        let mono = format!(
            "s{} a{}",
            if shell_mono { "↓" } else { "↑" },
            if atom_mono { "↓" } else { "↑" }
        );
        println!(
            "{:<28} {:>3} {:>4}  {:<28}  {:<18} {mono}",
            name,
            "ok",
            hit.steps.len(),
            fmt_shell_path(&shell_costs),
            fmt_atom_path(&atom_costs),
        );
    }
}

fn main() {
    println!(
        "shell_cost_bench — residual keep-if-drop (gate is comparison, not gold)\n\
         Never uses |current−target| as a score (that was a bug).\n\
         Proj cost = |projected−target|, projected = current + editδ.\n\
         columns: drop_rate agree_rate path_miss hop0-top1/known\n"
    );
    let t0 = Instant::now();
    let mut history: Vec<(String, SuiteStats)> = Vec::new();

    println!("## Chain A — from legacy at_sites baseline");
    let (best_a_mode, best_a) =
        run_greedy_chain(Mode::at_sites_baseline(), "A BASELINE", &mut history);

    println!("\n## Chain B — from Σ|proj−tgt| keep-if-drop (shells only)");
    let (best_b_mode, best_b) = run_greedy_chain(Mode::proj_baseline(), "B PROJ", &mut history);

    print_ledger(&history);

    let (proj_stats, diss) = run_mode(Mode::proj_baseline(), true);
    println!(
        "\n## Path check — B PROJ  mid miss={} drop={:.2} agree={:.2}  hard miss={} drop={:.2} agree={:.2}",
        proj_stats.mid.path_miss,
        proj_stats.mid.drop_rate(),
        proj_stats.mid.agree_rate(),
        proj_stats.hard.path_miss,
        proj_stats.hard.drop_rate(),
        proj_stats.hard.agree_rate(),
    );
    print_dis_report(&diss, &proj_stats.mid, &proj_stats.hard);

    // Head-to-head: full n0/n1/n2 vs H-closer/no-n2 (PatternInfo unchanged).
    println!("\n## Shell view: full n0/n1/n2 vs H→closer/no-n2");
    let views = [
        (
            "full n0/n1/n2",
            Mode::proj_baseline(),
        ),
        (
            "full +leave",
            Mode {
                leave: true,
                ..Mode::proj_baseline()
            },
        ),
        (
            "H→closer/no-n2",
            Mode {
                h_closer_no_n2: true,
                ..Mode::proj_baseline()
            },
        ),
        (
            "H→closer/no-n2 +leave",
            Mode {
                h_closer_no_n2: true,
                leave: true,
                ..Mode::proj_baseline()
            },
        ),
    ];
    for (tag, mode) in views {
        let (s, _) = run_mode(mode, false);
        print_row(tag, &s);
        history.push((format!("view {tag}"), s));
    }

    // Per-case live find_path cost traces (dropped when gate-as-gold was removed).
    print_find_path_costs("MID", MID);
    print_find_path_costs("HARD", HARD);

    let (best_mode, best, which) = if score_key(&best_b) > score_key(&best_a) {
        (best_b_mode, best_b, "B")
    } else {
        (best_a_mode, best_a, "A")
    };
    println!(
        "\nBEST overall [{which}]: {}  (miss_sum={} agree_sum={:.3})",
        best_mode.label(),
        best.mid.path_miss + best.hard.path_miss,
        best.mid.agree_rate() + best.hard.agree_rate()
    );
    println!("wall {:.2}s", t0.elapsed().as_secs_f64());
}
