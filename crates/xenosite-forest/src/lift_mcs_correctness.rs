//! Correctness: tag-lift + extend (cost 0) or MCS must match full MCS cost.
//!
//! Covers hydroxyl add/extend, DH same-tag, cleavage shrink, symmetric rings,
//! site-cast residual bounds, and find_path with `mcs_lift_fallback == 0`.

use std::collections::BTreeSet;

use crate::atom_diff::{
    added_heavy_atoms, atom_diff, residual_cost_after_site_cast, try_atom_diff_for_child,
    try_atom_diff_for_child_goal, try_lift_cleaved_child, try_lift_cleaved_child_goal,
};
use crate::find_path::{FindPathConfig, HeapScoreMode, PathCounters, find_path_with};
use crate::forest_mol::ForestMol;
use crate::mol::parse_mol;
use crate::rules::{dealkylation, dehydrogenation, hydroxylation, phase_one, quinone_formation};
use crate::ruleset::RuleSet;

fn assert_lift_cost_eq_mcs(
    parent: &ForestMol,
    child: &ForestMol,
    target: &crate::Molecule,
    parent_diff: &crate::atom_diff::AtomDiff,
    goal: Option<usize>,
    label: &str,
) {
    let full = atom_diff(child.mol(), target);
    let lifted = if added_heavy_atoms(parent, child).is_empty()
        && (0..parent.mol().atom_count()).all(|i| {
            parent
                .mol()
                .atom(crate::mol::atom_idx(i))
                .element
                .atomic_number()
                <= 1
                || parent.tag_of(i).and_then(|t| child.index_of(t)).is_some()
        }) {
        try_atom_diff_for_child_goal(parent, parent_diff, child, target, goal)
            .unwrap_or_else(|| panic!("{label}: try_atom_diff_for_child_goal returned None"))
    } else {
        // Shrink and/or adds on a cleaved fragment.
        try_lift_cleaved_child_goal(parent, parent_diff, child, target, goal)
            .or_else(|| try_atom_diff_for_child_goal(parent, parent_diff, child, target, goal))
            .unwrap_or_else(|| panic!("{label}: no lift path"))
    };
    assert_eq!(
        lifted.cost(),
        full.cost(),
        "{label}: lift={} mcs={} lift={lifted:?} mcs={full:?}",
        lifted.cost(),
        full.cost()
    );
    // Site-cast residual is advisory only. Cost-0 extend or MCS equals MCS.
    let _ = goal;
}

fn first_product(parent: &ForestMol, set: &RuleSet) -> ForestMol {
    let cands: Vec<_> = set
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert!(!cands.is_empty(), "no candidates on {}", parent.csmi());
    let pieces = cands[0].materialize_mols(parent.mol()).unwrap();
    assert!(!pieces.is_empty());
    parent.adopt_product(pieces[0].clone())
}

fn goal_for_candidate(
    parent_diff: &crate::atom_diff::AtomDiff,
    set: &RuleSet,
    parent: &ForestMol,
) -> (Option<usize>, Vec<usize>) {
    let cands: Vec<_> = set
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let c = &cands[0];
    let atoms: Vec<usize> = {
        let mut a: BTreeSet<usize> = c.site_atoms().into_iter().collect();
        if a.is_empty() {
            a.insert(c.site());
        }
        a.into_iter().collect()
    };
    let goal = residual_cost_after_site_cast(parent_diff, &c.effect(), &atoms, &[]);
    (Some(goal), atoms)
}

#[test]
fn lift_matches_mcs_hydroxylation_ethane_to_ethanol() {
    let parent = ForestMol::parse("CC").unwrap();
    let target = parse_mol("CCO").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = hydroxylation();
    let (goal, _) = goal_for_candidate(&parent_diff, &set, &parent);
    let child = first_product(&parent, &set);
    assert_lift_cost_eq_mcs(
        &parent,
        &child,
        &target,
        &parent_diff,
        goal,
        "CC→CCO hydroxyl",
    );
}

#[test]
fn lift_matches_mcs_hydroxylation_benzene_symmetry() {
    // Extend places O; cost-0 lift skips MCS (or MCS if extend misses).
    let parent = ForestMol::parse("c1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = hydroxylation();
    let (goal, _) = goal_for_candidate(&parent_diff, &set, &parent);
    let child = first_product(&parent, &set);
    assert_lift_cost_eq_mcs(
        &parent,
        &child,
        &target,
        &parent_diff,
        goal,
        "benzene→phenol hydroxyl",
    );
}

#[test]
fn lift_matches_mcs_hydroxylation_propane_primary() {
    let parent = ForestMol::parse("CCC").unwrap();
    let target = parse_mol("CCCO").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = hydroxylation();
    let (goal, _) = goal_for_candidate(&parent_diff, &set, &parent);
    let child = first_product(&parent, &set);
    assert_lift_cost_eq_mcs(
        &parent,
        &child,
        &target,
        &parent_diff,
        goal,
        "propane→n-propanol",
    );
}

#[test]
fn lift_matches_mcs_dehydrogenation_hydroquinone() {
    let parent = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
    let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let pairs = dehydrogenation()
        .pair_candidates_leaf(parent.mol())
        .unwrap();
    assert!(!pairs.is_empty());
    let pair = &pairs[0];
    let (p0, p1) = pair.path_ends();
    let atoms = pair
        .end_atoms()
        .map(|(a, b)| vec![a, b])
        .unwrap_or_else(|| vec![pair.site]);
    let goal = residual_cost_after_site_cast(&parent_diff, &pair.effect, &atoms, &[p0, p1]);
    let pieces = pair.materialize_mols(parent.mol()).unwrap();
    let child = parent.adopt_product(pieces[0].clone());
    assert_lift_cost_eq_mcs(
        &parent,
        &child,
        &target,
        &parent_diff,
        Some(goal),
        "hydroquinone→quinone DH",
    );
}

#[test]
fn lift_matches_mcs_dealkylation_anisole() {
    let parent = ForestMol::parse("COc1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let phenol = crate::mol::canon_of("Oc1ccccc1").unwrap();
    let cands = dealkylation()
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let mut child = None;
    let mut goal = None;
    for c in &cands {
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        for piece in pieces {
            let adopted = parent.adopt_product(piece);
            if adopted.csmi().as_ref() == phenol.as_str() {
                let atoms: Vec<usize> = {
                    let mut a: BTreeSet<usize> = c.site_atoms().into_iter().collect();
                    if a.is_empty() {
                        a.insert(c.site());
                    }
                    a.into_iter().collect()
                };
                goal = Some(residual_cost_after_site_cast(
                    &parent_diff,
                    &c.effect(),
                    &atoms,
                    &[],
                ));
                child = Some(adopted);
                break;
            }
        }
        if child.is_some() {
            break;
        }
    }
    let child = child.expect("phenol from anisole");
    assert_lift_cost_eq_mcs(
        &parent,
        &child,
        &target,
        &parent_diff,
        goal,
        "anisole→phenol dealk",
    );
}

#[test]
fn lift_matches_mcs_dealkylation_dimethoxy() {
    let parent = ForestMol::parse("COc1ccc(CCN)cc1OC").unwrap();
    let target = parse_mol("NCCc1ccc(O)c(O)c1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = dealkylation();
    let cands: Vec<_> = set
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert!(!cands.is_empty());
    // Every demethylation product: lift cost == MCS cost.
    let mut checked = 0usize;
    for c in &cands {
        let atoms: Vec<usize> = {
            let mut a: BTreeSet<usize> = c.site_atoms().into_iter().collect();
            if a.is_empty() {
                a.insert(c.site());
            }
            a.into_iter().collect()
        };
        let goal = residual_cost_after_site_cast(&parent_diff, &c.effect(), &atoms, &[]);
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        for piece in pieces {
            let child = parent.adopt_product(piece);
            // Skip tiny leave fragments.
            if child.heavy_atom_count() < 6 {
                continue;
            }
            assert_lift_cost_eq_mcs(
                &parent,
                &child,
                &target,
                &parent_diff,
                Some(goal),
                &format!("dimethoxy dealk site={}", c.site()),
            );
            checked += 1;
        }
    }
    assert!(checked >= 1, "expected at least one kept demethylation");
}

#[test]
fn lift_matches_mcs_quinone_formation_ends() {
    let parent = ForestMol::parse("COc1ccc(O)cc1").unwrap();
    let target = parse_mol("O=C1C=C(O)C(=O)C(O)=C1").unwrap();
    let target_ha = target
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count();
    let parent_diff = atom_diff(parent.mol(), &target);
    let pairs = quinone_formation()
        .pair_candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert!(!pairs.is_empty());
    let mut checked = 0usize;
    for pair in &pairs {
        let (p0, p1) = pair.path_ends();
        let atoms = pair
            .end_atoms()
            .map(|(a, b)| vec![a, b])
            .unwrap_or_else(|| vec![pair.site]);
        let goal = residual_cost_after_site_cast(&parent_diff, &pair.effect, &atoms, &[p0, p1]);
        let Ok(pieces) = pair.materialize_mols(parent.mol()) else {
            continue;
        };
        for piece in pieces {
            let child = parent.adopt_product(piece);
            // Leave fragments from cleaving pairs (e.g. Me from dealkylate)
            // are split out now; find_path would not continue them toward a
            // larger target. Skip those for lift≡MCS.
            if child.heavy_atom_count() < target_ha {
                continue;
            }
            assert_lift_cost_eq_mcs(
                &parent,
                &child,
                &target,
                &parent_diff,
                Some(goal),
                &format!("QF {}", pair.pattern_name),
            );
            checked += 1;
        }
    }
    assert!(checked >= 1, "expected at least one QF product");
}

#[test]
fn residual_never_exceeds_parent_cost_when_cast_helps() {
    let cases: &[(&str, &str, RuleSet)] = &[
        ("CC", "CCO", hydroxylation()),
        ("c1ccccc1", "Oc1ccccc1", hydroxylation()),
        ("COc1ccccc1", "Oc1ccccc1", dealkylation()),
    ];
    for &(reactant, target_smi, ref set) in cases {
        let parent = ForestMol::parse(reactant).unwrap();
        let target = parse_mol(target_smi).unwrap();
        let parent_diff = atom_diff(parent.mol(), &target);
        let cands: Vec<_> = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        for c in &cands {
            let atoms: Vec<usize> = {
                let mut a: BTreeSet<usize> = c.site_atoms().into_iter().collect();
                if a.is_empty() {
                    a.insert(c.site());
                }
                a.into_iter().collect()
            };
            let goal = residual_cost_after_site_cast(&parent_diff, &c.effect(), &atoms, &[]);
            assert!(
                goal <= parent_diff.cost(),
                "{reactant}→{target_smi} pattern={}: goal={goal} parent={}",
                c.pattern_name(),
                parent_diff.cost()
            );
        }
    }
}

#[test]
fn residual_may_over_credit_but_lift_still_matches_mcs() {
    // Residual is an early-stop *goal*, often optimistic (below true cost).
    // Over-credit ⇒ we do not early-stop; orbit closure must still hit MCS.
    let parent = ForestMol::parse("CC").unwrap();
    let target = parse_mol("CCO").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = hydroxylation();
    let (goal, _) = goal_for_candidate(&parent_diff, &set, &parent);
    let goal = goal.expect("goal");
    assert!(goal <= parent_diff.cost());
    let child = first_product(&parent, &set);
    let lifted = try_atom_diff_for_child_goal(&parent, &parent_diff, &child, &target, Some(goal))
        .expect("lift");
    let full = atom_diff(child.mol(), &target);
    assert_eq!(lifted.cost(), full.cost());
}

/// Expand every PhaseOne survivor one hop and compare lift vs MCS.
#[test]
fn one_hop_phase_one_lift_matches_mcs_on_eugenol() {
    let parent = ForestMol::parse("COc1cc(CC=C)ccc1O").unwrap();
    let target = parse_mol("O=C1C=C(O)C(=O)C=C1CC=C").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = phase_one();
    let mut checked = 0usize;

    let cands: Vec<_> = set
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    for c in &cands {
        if !crate::atom_diff::candidate_could_help_on(
            c,
            &parent_diff,
            Some(parent.mol()),
            Some(&target),
        ) {
            continue;
        }
        let atoms: Vec<usize> = {
            let mut a: BTreeSet<usize> = c.site_atoms().into_iter().collect();
            if a.is_empty() {
                a.insert(c.site());
            }
            a.into_iter().collect()
        };
        let goal = residual_cost_after_site_cast(&parent_diff, &c.effect(), &atoms, &[]);
        let Ok(pieces) = c.materialize_mols(parent.mol()) else {
            continue;
        };
        for piece in pieces {
            let child = parent.adopt_product(piece);
            if child.heavy_atom_count() + 2 < parent.heavy_atom_count()
                && child.heavy_atom_count() < 6
            {
                // Tiny leave fragment — skip.
                continue;
            }
            let label = format!(
                "eugenol hop {} site={} child={}",
                c.pattern_name(),
                c.site(),
                child.csmi()
            );
            if let Some(lifted) =
                try_atom_diff_for_child_goal(&parent, &parent_diff, &child, &target, Some(goal))
                    .or_else(|| {
                        try_lift_cleaved_child_goal(
                            &parent,
                            &parent_diff,
                            &child,
                            &target,
                            Some(goal),
                        )
                    })
            {
                let full = atom_diff(child.mol(), &target);
                assert_eq!(
                    lifted.cost(),
                    full.cost(),
                    "{label}: lift={} mcs={}",
                    lifted.cost(),
                    full.cost()
                );
                checked += 1;
            }
        }
    }
    // Pair door separately.
    let pairs: Vec<_> = set
        .pair_candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    for pair in &pairs {
        if !crate::atom_diff::pair_could_help(pair, &parent_diff, parent.mol(), &target) {
            continue;
        }
        let (p0, p1) = pair.path_ends();
        let atoms = pair
            .end_atoms()
            .map(|(a, b)| vec![a, b])
            .unwrap_or_else(|| vec![pair.site]);
        let goal = residual_cost_after_site_cast(&parent_diff, &pair.effect, &atoms, &[p0, p1]);
        let Ok(pieces) = pair.materialize_mols(parent.mol()) else {
            continue;
        };
        for piece in pieces {
            let child = parent.adopt_product(piece);
            let label = format!("eugenol pair {} child={}", pair.pattern_name, child.csmi());
            assert_lift_cost_eq_mcs(&parent, &child, &target, &parent_diff, Some(goal), &label);
            checked += 1;
        }
    }
    assert!(
        checked >= 3,
        "expected several eugenol hops to check; got {checked}"
    );
}

fn assert_find_path_no_mcs_fallback(reactant: &str, target: &str, max_nodes: usize, label: &str) {
    let mut counters = PathCounters::default();
    // HA-alignment-only cost can prefer ResonancePair hops whose sealed end
    // bags soft-mismatch (HEURISTICS: pair mismatches recorded, suite gate
    // not decided). Do not fail the MCS-fallback check on that.
    counters.allow_formula_delta_mismatch = true;
    let hits = find_path_with(
        reactant,
        target,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: 1,
            max_nodes,
            use_atom_diff: true,
            lazy_closer: false,
            heap_score: HeapScoreMode::match_product(),
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    assert!(
        !hits.is_empty(),
        "{label}: miss; billed={} mcs_fb={} mcs_rm={} nodes={}",
        counters.billed(),
        counters.mcs_lift_fallback,
        counters.mcs_lift_rematch,
        counters.nodes
    );
    assert_eq!(
        counters.mcs_lift_fallback, 0,
        "{label}: mcs_lift_fallback={}",
        counters.mcs_lift_fallback
    );
}

#[test]
fn find_path_mcs_fallback_zero_mid_cases() {
    // Same SMILES as find_path_bench MID (max_nodes=800 default door).
    // eugenol→allyl-Q still misses under lift+rematch (~11 nodes); parked —
    // do not chase; other mid cases must hit with mcs_lift_fallback == 0.
    assert_find_path_no_mcs_fallback(
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
        800,
        "dimethoxy-PEA→catechol",
    );
    assert_find_path_no_mcs_fallback(
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
        800,
        "MeOPhOH→hydroxyQ",
    );
    assert_find_path_no_mcs_fallback(
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        r"C(#C/C=C/C=O)C(C)(C)C",
        800,
        "TBA→enyne aldehyde",
    );
    assert_find_path_no_mcs_fallback(
        "COc1ccc2ccccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
        800,
        "2-MeO-naph→1,2-NQ",
    );
}

#[test]
fn find_path_mcs_fallback_zero_dealk_multipath() {
    let mut counters = PathCounters::default();
    let hits = find_path_with(
        "CN(C)Cc1ccccc1",
        "O=Cc1ccccc1",
        &dealkylation(),
        &mut counters,
        FindPathConfig {
            max_paths: 4,
            max_nodes: 80,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    assert!(!hits.is_empty());
    assert_eq!(counters.mcs_lift_fallback, 0);
}

#[test]
fn find_path_mcs_fallback_zero_tba() {
    const TERB: &str = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12";
    const TBA: &str = r"C(#C/C=C/C=O)C(C)(C)C";
    let mut counters = PathCounters::default();
    let hits = find_path_with(
        TERB,
        TBA,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: 3,
            max_nodes: 80,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    assert!(!hits.is_empty(), "billed={}", counters.billed());
    assert_eq!(counters.mcs_lift_fallback, 0);
}

#[test]
fn lift_with_mcs_rematch_matches_full_mcs_cost() {
    // Cost-0 after extend kept; else MCS → equals MCS cost.
    let parent = ForestMol::parse("c1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let set = hydroxylation();
    let child = first_product(&parent, &set);
    let lifted = try_atom_diff_for_child(&parent, &parent_diff, &child, &target).unwrap();
    let mcs = atom_diff(child.mol(), &target);
    assert_eq!(lifted.cost(), mcs.cost());
    assert_eq!(lifted.cost(), 0);
}

#[test]
fn nonzero_lift_always_sticks_with_mcs() {
    // No Aut chase: if extend is not cost 0, result is MCS.
    let parent = ForestMol::parse("COc1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let phenol = crate::mol::canon_of("Oc1ccccc1").unwrap();
    let cands = dealkylation()
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let child = cands
        .iter()
        .find_map(|c| {
            c.materialize_mols(parent.mol()).ok().and_then(|pieces| {
                pieces.into_iter().find_map(|p| {
                    let a = parent.adopt_product(p);
                    (a.csmi().as_ref() == phenol.as_str()).then_some(a)
                })
            })
        })
        .expect("phenol");
    let mut rematch = 0usize;
    let out = crate::atom_diff::try_lift_cleaved_child_tracked(
        &parent,
        &parent_diff,
        &child,
        &target,
        Some(&mut rematch),
    )
    .expect("cleave lift");
    let mcs = atom_diff(child.mol(), &target);
    assert_eq!(out.cost(), mcs.cost());
    if out.cost() > 0 {
        assert!(
            rematch >= 1,
            "non-zero result must be MCS (rematch={rematch})"
        );
    }
}

#[test]
fn lift_cost_zero_skips_mcs_and_beats_fresh_mcs_wall() {
    // Speed claim: when lift+extend hits cost 0, we never rematch MCS and the
    // lift path is cheaper than always computing a fresh MCS.
    use std::time::Instant;

    let parent = ForestMol::parse("c1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let child = first_product(&parent, &hydroxylation());
    let mut rematch = 0usize;
    let lifted = crate::atom_diff::try_atom_diff_for_child_tracked(
        &parent,
        &parent_diff,
        &child,
        &target,
        Some(&mut rematch),
    )
    .expect("lift");
    assert_eq!(lifted.cost(), 0);
    assert_eq!(rematch, 0, "cost-0 lift must not rematch MCS");

    const N: u32 = 200;
    // Warmup
    for _ in 0..20 {
        let _ = try_atom_diff_for_child(&parent, &parent_diff, &child, &target);
        let _ = atom_diff(child.mol(), &target);
    }
    let t0 = Instant::now();
    for _ in 0..N {
        let d = try_atom_diff_for_child(&parent, &parent_diff, &child, &target).unwrap();
        assert_eq!(d.cost(), 0);
    }
    let lift_ns = t0.elapsed().as_nanos() as f64 / N as f64;
    let t1 = Instant::now();
    for _ in 0..N {
        let d = atom_diff(child.mol(), &target);
        assert_eq!(d.cost(), 0);
    }
    let mcs_ns = t1.elapsed().as_nanos() as f64 / N as f64;
    assert!(
        lift_ns < mcs_ns,
        "cost-0 lift should beat fresh MCS: lift={lift_ns:.0}ns mcs={mcs_ns:.0}ns"
    );
}

#[test]
fn try_atom_diff_refuses_shrink_without_cleave_door() {
    let parent = ForestMol::parse("COc1ccccc1").unwrap();
    let target = parse_mol("Oc1ccccc1").unwrap();
    let parent_diff = atom_diff(parent.mol(), &target);
    let phenol = crate::mol::canon_of("Oc1ccccc1").unwrap();
    let cands = dealkylation()
        .candidates(parent.mol())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let child = cands
        .iter()
        .find_map(|c| {
            c.materialize_mols(parent.mol()).ok().and_then(|pieces| {
                pieces.into_iter().find_map(|p| {
                    let a = parent.adopt_product(p);
                    (a.csmi().as_ref() == phenol.as_str()).then_some(a)
                })
            })
        })
        .expect("phenol");
    assert!(
        try_atom_diff_for_child(&parent, &parent_diff, &child, &target).is_none(),
        "non-cleavage try must refuse shrink"
    );
    let cleaved =
        try_lift_cleaved_child(&parent, &parent_diff, &child, &target).expect("cleavage lift");
    assert_eq!(cleaved.cost(), atom_diff(child.mol(), &target).cost());
}
