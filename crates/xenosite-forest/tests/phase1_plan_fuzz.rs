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
    dehydrogenation, epoxidation, hydroxylation, n_dealkylation, phase_one, quinone_formation,
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
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap()
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
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
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
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
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
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap()
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
        .unwrap();
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
                .metabolize(&mol, accept_all_rules, accept_all_sites, true)
                .unwrap()
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
        .unwrap();
        if let Some(hit) = hits.first() {
            assert_plan_elementary(hit, target);
        }
    }
}
