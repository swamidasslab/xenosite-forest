//! Fuzz property: apply a leaf → product must appear in `find_path_ms1`.
//!
//! Pattern (same shape as [`phase1_plan_fuzz`]):
//! 1. Draw a reactant SMILES and a catalog leaf name.
//! 2. Materialize every distinct first-fragment product of one apply.
//! 3. For each product, run `find_path_ms1` at that product's [M+H]⁺ with
//!    `ApplyN { arms: [leaf], count: 1 }`.
//! 4. Assert the product CSMI is among the hits, and **every emitted hit's
//!    final product** (hit smiles + last path-step product) lies within tol
//!    of the target m/z. Intermediate steps need not match the target mass;
//!    they must be monotonically closer.
//!
//! Case count: `XENOSITE_FUZZ_EXAMPLES` when set, else 8.

use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;
use xenosite_forest::mass::{Ms1Adduct, mz_abs_error, mz_of_mol, mz_within};
use xenosite_forest::{
    ApplyN, ForestMol, Ms1Config, PathCounters, PathOutcome, canon_smiles, find_path_ms1,
    leaf_rule,
};

fn fuzz_config(default_cases: u32) -> ProptestConfig {
    let cases = std::env::var("XENOSITE_FUZZ_EXAMPLES")
        .ok()
        .and_then(|raw| raw.parse().ok())
        .unwrap_or(default_cases);
    ProptestConfig::with_cases(cases)
}

/// Substrates that Phase I leaves often touch.
fn reactant_corpus() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("CC"),
        Just("CCC"),
        Just("CCO"),
        Just("C=C"),
        Just("C#C"),
        Just("C/C=C/C"),
        Just("C1OC1"),
        Just("CC1OC1C"),
        Just("c1ccccc1"),
        Just("Cc1ccccc1"),
        Just("Oc1ccccc1"),
        Just("CCCCc1ccccc1"),
        Just("CCS"),
        Just("CSC"),
        Just("CN"),
        Just("CCN"),
        Just("CN(C)c1ccccc1"),
        Just("COc1ccccc1"),
        Just("COc1ccc(OC)cc1"),
        Just("CC(O)C"),
        Just("Clc1ccccc1"),
        Just("Brc1ccccc1"),
        Just("CC(=O)OC"),
        Just("CC(=O)NCC"),
        Just("CS(=O)C"),
        Just("O=Nc1ccccc1"),
        Just("c1ccc2c(c1)C(=O)c1ccccc1C2=O"),
        Just("c1ccc2c(c1)OCO2"),
    ]
}

/// Catalog leaves with reliable one-hop mass change on the corpus.
fn leaf_corpus() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("Hydroxylation"),
        Just("Epoxidation"),
        Just("EpoxideHydration"),
        Just("EpoxideOpening"),
        Just("SulfurOxidation"),
        Just("NitrogenOxidation"),
        Just("Dehydrogenation"),
        Just("Dealkylation"),
        Just("NDealkylation"),
        Just("Hydrogenation"),
        Just("OxidativeDehalogenation"),
        Just("ReductiveDehalogenation"),
        Just("Hydrolysis"),
        Just("Dehydration"),
        Just("SulfurReduction"),
        Just("NitrogenReduction"),
        Just("OxygenReduction"),
        Just("BenzodioxoleReduction"),
    ]
}

fn products_from_apply(reactant: &str, leaf: &str) -> Vec<(String, f64)> {
    let Some(set) = leaf_rule(leaf) else {
        return Vec::new();
    };
    let Ok(parent) = ForestMol::parse(reactant) else {
        return Vec::new();
    };
    let mol = parent.mol();
    let mut out = Vec::new();
    for cand in set.candidates(mol) {
        let Ok(cand) = cand else {
            continue;
        };
        let Ok(pieces) = cand.materialize_mols(mol) else {
            continue;
        };
        if pieces.is_empty() {
            continue;
        }
        let child = parent.adopt_product(pieces[0].clone());
        let Some(mz) = mz_of_mol(child.mol(), Ms1Adduct::MPlusH) else {
            continue;
        };
        out.push((child.csmi().as_ref().to_string(), mz));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out.dedup_by(|a, b| a.0 == b.0);
    out
}

/// Every emitted hit's product(s) must satisfy the target m/z.
///
/// Replays each plan linearization: sole products must all lie within tol;
/// cleavage leave fragments are exempt (hit species must appear at target mz).
fn assert_hits_mz_ok(hits: &[PathOutcome], reactant: &str, mz: f64, tol_da: f64) {
    let start = ForestMol::parse(reactant).expect("reactant");
    let start_mz = mz_of_mol(start.mol(), Ms1Adduct::MPlusH).expect("start mz");
    let reactant_mol = start.mol();
    for h in hits {
        let mut prev_err = mz_abs_error(start_mz, mz);
        for (i, step) in h.steps.iter().enumerate() {
            let step_mol = ForestMol::parse(&step.product).expect("step product");
            let step_mz = mz_of_mol(step_mol.mol(), Ms1Adduct::MPlusH).expect("step mz");
            let err = mz_abs_error(step_mz, mz);
            assert!(
                err <= prev_err + 1e-6,
                "step {i} mz err grew: {prev_err} → {err} ({})",
                step.product
            );
            prev_err = err;
            let last = i + 1 == h.steps.len();
            if last {
                assert!(
                    mz_within(step_mz, mz, tol_da),
                    "last step {} mz={step_mz} outside tol of {mz}",
                    step.product
                );
                assert_eq!(step.product, h.smiles, "last step must equal hit smiles");
            }
        }
        let hit_mol = ForestMol::parse(&h.smiles).expect("hit");
        let hit_mz = mz_of_mol(hit_mol.mol(), Ms1Adduct::MPlusH).expect("hit mz");
        assert!(
            mz_within(hit_mz, mz, tol_da),
            "hit {} mz={hit_mz} outside tol of {mz}",
            h.smiles
        );

        let mut any_replay = false;
        for (li, lin) in h.plan.linearizations().into_iter().enumerate() {
            let products = lin.apply(reactant_mol).unwrap_or_default();
            if products.is_empty() {
                continue;
            }
            any_replay = true;
            if products.len() == 1 {
                let pmz = mz_of_mol(&products[0], Ms1Adduct::MPlusH).expect("plan product mz");
                assert!(
                    mz_within(pmz, mz, tol_da),
                    "plan lin{li} product mz={pmz} outside tol of {mz}"
                );
            } else {
                // Mass-identical isomers all at target, or cleavage leave
                // fragments off-target (then hit species must appear at mz).
                let all_at_target = products.iter().all(|p| {
                    mz_of_mol(p, Ms1Adduct::MPlusH)
                        .is_some_and(|pmz| mz_within(pmz, mz, tol_da))
                });
                if !all_at_target {
                    let hit_ok = products.iter().any(|p| {
                        canon_smiles(p) == h.smiles
                            && mz_of_mol(p, Ms1Adduct::MPlusH)
                                .is_some_and(|pmz| mz_within(pmz, mz, tol_da))
                    });
                    assert!(
                        hit_ok,
                        "plan lin{li} missing hit {} at target mz; got {:?}",
                        h.smiles,
                        products.iter().map(canon_smiles).collect::<Vec<_>>()
                    );
                }
            }
        }
        assert!(
            any_replay || h.plan.steps().is_empty(),
            "emitted plan for {} must replay from reactant; plan={:?}",
            h.smiles,
            h.plan
        );
    }
}

proptest! {
    #![proptest_config(fuzz_config(8))]

    /// Apply leaf once; every materialized product CSMI must appear in
    /// `find_path_ms1` hits at that product's [M+H]⁺, and every emitted hit's
    /// final product must satisfy the m/z target.
    #[test]
    fn fuzz_apply_product_recoverable_by_ms1(
        reactant in reactant_corpus(),
        leaf in leaf_corpus(),
    ) {
        let products = products_from_apply(reactant, leaf);
        prop_assume!(!products.is_empty());

        let set = leaf_rule(leaf).expect(leaf);
        let pools = [ApplyN::new([leaf], 1)];
        for (csmi, mz) in &products {
            let mut counters = PathCounters::default();
            let hits = find_path_ms1(
                reactant,
                &set,
                &pools,
                &mut counters,
                Ms1Config {
                    mz: *mz,
                    tol_da: 0.001,
                    adduct: Ms1Adduct::MPlusH,
                    max_paths: 16,
                    max_nodes: 800,
                },
            )
            .expect("find_path_ms1");
            prop_assert!(
                hits.iter().any(|h| h.smiles == *csmi),
                "{reactant} + {leaf} → {csmi} (mz={mz}) missing from {:?}; billed={}",
                hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
                counters.billed()
            );
            prop_assert!(!hits.is_empty());
            assert_hits_mz_ok(&hits, reactant, *mz, 0.001);
        }
    }
}

proptest! {
    #![proptest_config(fuzz_config(4))]

    /// Two successive hydroxylations: chain-apply, then recover under ApplyN
    /// count=2. Harder than one-hop; OR pool is still a single leaf name.
    #[test]
    fn fuzz_two_hydroxylations_recoverable_by_ms1(
        reactant in prop_oneof![
            Just("CC"),
            Just("CCC"),
            Just("c1ccccc1"),
            Just("Cc1ccccc1"),
            Just("CCCCc1ccccc1"),
            Just("CCO"),
        ],
    ) {
        let set = leaf_rule("Hydroxylation").expect("Hydroxylation");
        let first_products = products_from_apply(reactant, "Hydroxylation");
        prop_assume!(!first_products.is_empty());

        // Chain: first apply → second apply on that child.
        let (first_csmi, _) = &first_products[0];
        let mid = ForestMol::parse(first_csmi).expect("mid");
        let mut second = None;
        for cand in set.candidates(mid.mol()) {
            let Ok(cand) = cand else { continue };
            let Ok(pieces) = cand.materialize_mols(mid.mol()) else {
                continue;
            };
            if pieces.is_empty() {
                continue;
            }
            let child = mid.adopt_product(pieces[0].clone());
            let mz = mz_of_mol(child.mol(), Ms1Adduct::MPlusH).expect("mz");
            second = Some((child.csmi().as_ref().to_string(), mz));
            break;
        }
        prop_assume!(second.is_some());
        let (csmi, mz) = second.unwrap();

        let pools = [ApplyN::new(["Hydroxylation"], 2)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            reactant,
            &set,
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 16,
                max_nodes: 2000,
            },
        )
        .expect("find_path_ms1");
        prop_assert!(
            hits.iter().any(|h| h.smiles == csmi),
            "{reactant} OH×2 → {csmi} missing from {:?}; billed={}",
            hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
            counters.billed()
        );
        assert_hits_mz_ok(&hits, reactant, mz, 0.001);
    }
}
