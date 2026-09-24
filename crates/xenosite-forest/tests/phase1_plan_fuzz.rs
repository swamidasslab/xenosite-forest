//! Port of Python phase-I / find_path plan fuzz (`proptest`).
//!
//! Sources:
//! - `tests/forest/test_phase1_steps_fuzz.py`
//! - `tests/forest/test_find_path_phase1_plan_fuzz.py`
//!
//! Rust has no `canonical_emitted_sites` knob yet; strategies still draw a
//! bool where Python did so the property surface stays parallel.
//!
//! Case count: `PROPTEST_CASES`, or `XENOSITE_FUZZ_EXAMPLES` when set (Python
//! parity). Defaults match the Hypothesis suites (4).

use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;
use xenosite_forest::rules::{
    dealkylation, dehydrogenation, epoxidation, hydroxylation, n_dealkylation, phase_one,
    quinone_formation,
};
use xenosite_forest::{
    FindPathConfig, PathCounters, PathOutcome, accept_all_rules, accept_all_sites, canon_of,
    find_path_with, parse_mol,
};

fn fuzz_config(default_cases: u32) -> ProptestConfig {
    let cases = std::env::var("XENOSITE_FUZZ_EXAMPLES")
        .ok()
        .and_then(|raw| raw.parse().ok())
        .unwrap_or(default_cases);
    ProptestConfig::with_cases(cases)
}

fn corpus_smiles() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("c1ccccc1"),
        Just("Oc1ccccc1"),
        Just("Oc1ccc(O)cc1"),
        Just("Nc1ccc(O)cc1"),
        Just("C=C"),
        Just("CCN"),
    ]
}

fn find_path_corpus() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("c1ccccc1"),
        Just("Oc1ccccc1"),
        Just("Oc1ccc(O)cc1"),
        Just("Cc1ccccc1"),
        Just("Nc1ccc(O)cc1"),
        Just("COc1ccccc1"),
        Just("CN(C)Cc1ccccc1"),
        Just("CCO"),
        Just("CCN"),
        Just("C=C"),
        Just("COc1ccc(O)cc1"),
        Just("COc1ccc(CC=C)cc1O"),
        Just("COc1ccc2ccccc2c1"),
        Just("N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1"),
        Just("c1ccccc1C(=O)OC(C)(C)C"),
        Just("CC(=O)Oc1ccc(OC(C)=O)cc1"),
    ]
}

/// Known multipath-friendly reactant→target pairs (bench / Maybe cases).
fn multipath_pairs() -> impl Strategy<Value = (&'static str, &'static str)> {
    prop_oneof![
        Just(("CN(C)Cc1ccccc1", "O=Cc1ccccc1")),
        Just(("COc1ccccc1", "Oc1ccccc1")),
        Just(("c1ccccc1C(=O)OC(C)(C)C", "O=C(O)c1ccccc1")),
        Just(("CC(=O)Oc1ccc(OC(C)=O)cc1", "Oc1ccc(O)cc1")),
        Just(("COc1ccc(O)cc1", "O=C1C=C(O)C(=O)C(O)=C1")),
        Just(("Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1")),
        Just(("c1ccccc1", "O=C1C=CC(=O)C=C1")),
        Just(("Oc1ccccc1", "O=C1C=CC(=O)C=C1")),
        Just(("COc1ccc(CC=C)cc1O", "O=C1C=CC(=O)C(CC=C)=C1")),
        Just(("COc1ccc(CCN)cc1OC", "NCCc1ccc(O)c(O)c1")),
        Just(("N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1", "O=Cc1ccccc1")),
        Just(("COc1ccc2ccccc2c1", "O=C1C(=O)c2ccccc2C=C1")),
        Just((
            "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
            r"C(#C/C=C/C=O)C(C)(C)C",
        )),
    ]
}

fn heavy_atom_count(mol: &xenosite_forest::Molecule) -> usize {
    mol.atoms()
        .filter(|(_, atom)| atom.element.atomic_number() > 1)
        .count()
}

fn assert_plan_elementary(outcome: &PathOutcome, target: &str) {
    assert!(!outcome.smiles.contains('.'));
    assert_eq!(&outcome.smiles, target);
    assert!(!outcome.plan.is_empty());
    let names: Vec<_> = outcome.plan.iter().map(|s| s.rule.as_str()).collect();
    assert!(
        !names.contains(&"QuinoneFormation"),
        "opaque QF step in plan {names:?}"
    );
}

/// HEURISTICS: multipath hits must not share total orders (no Deps reordering).
fn assert_no_linearization_overlap(hits: &[PathOutcome]) {
    for (i, a) in hits.iter().enumerate() {
        for (j, b) in hits.iter().enumerate().skip(i + 1) {
            assert!(
                !a.plan.same_linearizations(&b.plan),
                "hit{i} same_linearizations hit{j}: {:?} vs {:?}",
                a.plan.iter().map(|s| &s.rule).collect::<Vec<_>>(),
                b.plan.iter().map(|s| &s.rule).collect::<Vec<_>>()
            );
            let ov = a.plan.linearization_overlap(&b.plan);
            assert_eq!(
                ov,
                0,
                "hit{i} linearization_overlap={ov} with hit{j}: {:?} vs {:?}",
                a.plan.iter().map(|s| &s.rule).collect::<Vec<_>>(),
                b.plan.iter().map(|s| &s.rule).collect::<Vec<_>>()
            );
        }
    }
}

proptest! {
    #![proptest_config(fuzz_config(4))]

    /// Python `test_benzene_quinone_plan_ends_in_dehydrogenation`.
    #[test]
    fn benzene_quinone_plan_ends_in_dehydrogenation(
        _canonical_emitted_sites in any::<bool>(),
    ) {
        let mol = parse_mol("c1ccccc1").unwrap();
        let rule = quinone_formation();
        let mut plans = Vec::new();
        for emission in rule
            .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap()
        {
            if !emission.plan.is_empty() {
                plans.push(emission.plan);
            }
            if plans.len() >= 4 {
                break;
            }
        }
        prop_assert!(!plans.is_empty());
        for steps in &plans {
            prop_assert_eq!(steps.last().map(|s| s.rule.as_str()), Some("Dehydrogenation"));
            prop_assert!(
                steps[..steps.len() - 1]
                    .iter()
                    .any(|s| s.rule == "Hydroxylation"),
                "expected OH prep before DH, got {:?}",
                steps.iter().map(|s| &s.rule).collect::<Vec<_>>()
            );
        }
    }

    /// Python `test_epoxidation_plan_is_one_step`.
    #[test]
    fn epoxidation_plan_is_one_step(_canonical_emitted_sites in any::<bool>()) {
        let mol = parse_mol("C=C").unwrap();
        let rule = epoxidation();
        let emissions = rule
            .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap();
        prop_assume!(!emissions.is_empty());
        let emission = &emissions[0];
        let names: Vec<_> = emission.plan.iter().map(|s| s.rule.as_str()).collect();
        prop_assert_eq!(names, vec!["Epoxidation"]);
        prop_assert!(emission.products.iter().all(|p| !p.contains('.')));
    }

    /// Python `test_ndealkylation_plan_is_one_step`.
    #[test]
    fn ndealkylation_plan_is_one_step(_canonical_emitted_sites in any::<bool>()) {
        let mol = parse_mol("CCN").unwrap();
        let rule = n_dealkylation();
        let emissions = rule
            .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap();
        prop_assume!(!emissions.is_empty());
        let emission = &emissions[0];
        let names: Vec<_> = emission.plan.iter().map(|s| s.rule.as_str()).collect();
        prop_assert_eq!(names, vec!["NDealkylation"]);
        prop_assert!(emission.products.iter().all(|p| !p.contains('.')));
    }

    /// Python `test_fuzz_quinone_plans_end_in_dehydrogenation`.
    #[test]
    fn fuzz_quinone_plans_end_in_dehydrogenation(
        smiles in corpus_smiles(),
        _canonical_emitted_sites in any::<bool>(),
    ) {
        let mol = parse_mol(smiles).unwrap();
        let rule = quinone_formation();
        let mut seen = 0usize;
        for emission in rule
            .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap()
        {
            if emission.plan.is_empty() {
                continue;
            }
            seen += 1;
            prop_assert_eq!(
                emission.plan.last().map(|s| s.rule.as_str()),
                Some("Dehydrogenation")
            );
            prop_assert!(emission.products.iter().all(|p| !p.contains('.')));
            if seen >= 4 {
                break;
            }
        }
        prop_assume!(seen > 0);
    }

    /// Python `test_benzene_quinone_plan_is_elementary`.
    #[test]
    fn benzene_quinone_plan_is_elementary(_canonical_emitted_sites in any::<bool>()) {
        let target = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "c1ccccc1",
            &target,
            &phase_one(),
            &mut counters,
            FindPathConfig {
                max_nodes: 400,
                max_paths: 1,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap().collect_all().unwrap();
        prop_assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_plan_elementary(&hits[0], &target);
        let names: Vec<_> = hits[0].plan.iter().map(|s| s.rule.as_str()).collect();
        prop_assert_eq!(
            names.iter().filter(|&&n| n == "Hydroxylation").count(),
            2,
            "plan={:?}",
            names
        );
        prop_assert_eq!(
            names.iter().filter(|&&n| n == "Dehydrogenation").count(),
            1,
            "plan={:?}",
            names
        );
    }
}

proptest! {
    #![proptest_config(fuzz_config(4))]

    /// Python `test_fuzz_created_target_hit_or_honest_miss`.
    #[test]
    fn fuzz_created_target_hit_or_honest_miss(
        start in find_path_corpus(),
        _canonical_emitted_sites in any::<bool>(),
        index in any::<prop::sample::Index>(),
    ) {
        let mol = parse_mol(start).unwrap();
        let start_csmi = canon_of(start).unwrap();
        let mut seen = std::collections::HashSet::new();
        seen.insert(start_csmi);
        let create_rules = [quinone_formation(), hydroxylation(), dehydrogenation()];
        let mut pool = Vec::new();
        let start_ha = heavy_atom_count(&mol);
        'outer: for rule in &create_rules {
            for emission in rule
                .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap()
            {
                for product in &emission.products {
                    if product.is_empty() || product.contains('.') || !seen.insert(product.clone())
                    {
                        continue;
                    }
                    let product_mol = parse_mol(product).unwrap();
                    let min_ha = 4.max(start_ha / 3);
                    if heavy_atom_count(&product_mol) < min_ha {
                        continue;
                    }
                    pool.push(product.clone());
                    if pool.len() >= 8 {
                        break 'outer;
                    }
                }
            }
        }
        prop_assume!(!pool.is_empty());
        let target = &pool[index.index(pool.len())];
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            start,
            target,
            &phase_one(),
            &mut counters,
            FindPathConfig {
                max_nodes: 150,
                max_paths: 1,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap().collect_all().unwrap();
        if let Some(hit) = hits.first() {
            assert_plan_elementary(hit, target);
        }
    }
}

proptest! {
    #![proptest_config(fuzz_config(12))]

    /// Multipath hits for fixed reactant→target pairs share no linearizations.
    #[test]
    fn fuzz_multipath_pairs_no_linearization_overlap(
        (start, target) in multipath_pairs(),
    ) {
        let want = canon_of(target).unwrap();
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            start,
            &want,
            &phase_one(),
            &mut counters,
            FindPathConfig {
                max_nodes: 250,
                max_paths: 4,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap().collect_all().unwrap();
        prop_assume!(!hits.is_empty());
        prop_assert!(hits.iter().all(|h| h.smiles == want));
        if hits.len() >= 2 {
            assert_no_linearization_overlap(&hits);
        }
    }

    /// Create a Phase-I product, find_path with max_paths>1, assert no overlap.
    #[test]
    fn fuzz_created_target_multipath_no_linearization_overlap(
        start in find_path_corpus(),
        index in any::<prop::sample::Index>(),
    ) {
        let mol = parse_mol(start).unwrap();
        let start_csmi = canon_of(start).unwrap();
        let mut seen = std::collections::HashSet::new();
        seen.insert(start_csmi);
        let create_rules = [
            quinone_formation(),
            hydroxylation(),
            dehydrogenation(),
            dealkylation(),
        ];
        let mut pool = Vec::new();
        let start_ha = heavy_atom_count(&mol);
        'outer: for rule in &create_rules {
            for emission in rule
                .metabolize(&mol, accept_all_rules, accept_all_sites, true).collect::<Result<Vec<_>, _>>().unwrap()
            {
                for product in &emission.products {
                    if product.is_empty() || product.contains('.') || !seen.insert(product.clone())
                    {
                        continue;
                    }
                    let product_mol = parse_mol(product).unwrap();
                    let min_ha = 4.max(start_ha / 3);
                    if heavy_atom_count(&product_mol) < min_ha {
                        continue;
                    }
                    pool.push(product.clone());
                    if pool.len() >= 10 {
                        break 'outer;
                    }
                }
            }
        }
        prop_assume!(!pool.is_empty());
        let target = &pool[index.index(pool.len())];
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            start,
            target,
            &phase_one(),
            &mut counters,
            FindPathConfig {
                max_nodes: 200,
                max_paths: 4,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap().collect_all().unwrap();
        if hits.len() < 2 {
            return Ok(());
        }
        prop_assert!(hits.iter().all(|h| h.smiles == *target));
        assert_no_linearization_overlap(&hits);
    }
}
