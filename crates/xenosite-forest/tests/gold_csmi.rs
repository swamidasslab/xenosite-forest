//! Ported Python find_path golds as CSMI↔CSMI identity (not raw string match).
//!
//! Each hit's kept product must equal `canon_of(target)`.

use xenosite_forest::rules::{
    dealkylation, dehydrogenation, epoxidation, epoxide_opening, quinone_formation,
};
use xenosite_forest::{
    FindPathConfig, PathCounters, RuleSet, canon_of, find_path, find_path_with, hydroxylation,
    phase_one,
};

fn assert_path_hits_target(
    reactant: &str,
    target: &str,
    ruleset: &RuleSet,
    config: FindPathConfig,
) {
    let want = canon_of(target).expect(target);
    let mut counters = PathCounters::default();
    let hits = find_path_with(reactant, target, ruleset, &mut counters, config, |_| true)
        .unwrap_or_else(|e| panic!("{reactant} → {target}: {e}"));
    assert!(
        !hits.is_empty(),
        "{reactant} → {target}: no path (billed={} nodes={} edits={})",
        counters.billed(),
        counters.nodes,
        counters.mol_edits
    );
    assert_eq!(
        hits[0].smiles, want,
        "{reactant} → {target}: kept {:?} ≠ want {want}",
        hits[0].smiles
    );
}

// --- guided_path_gold -------------------------------------------------------

#[test]
fn gold_apap_to_napqi_phase_one() {
    assert_path_hits_target(
        "CC(=O)Nc1ccc(O)cc1",
        "CC(=O)N=C1C=CC(=O)C=C1",
        &phase_one(),
        FindPathConfig {
            max_paths: 2,
            max_nodes: 80,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_hydroquinone_to_bq_ring_end() {
    let qf = RuleSet::compose(
        Some("ringEndQf".into()),
        [
            quinone_formation(),
            hydroxylation(),
            dehydrogenation(),
            dealkylation(),
        ],
    );
    assert_path_hits_target(
        "Oc1ccc(O)cc1",
        "O=C1C=CC(=O)C=C1",
        &qf,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 120,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_apap_to_napqi_hydrox_deh() {
    let phase1 = RuleSet::compose(
        Some("ringEndP1".into()),
        [hydroxylation(), dehydrogenation()],
    );
    assert_path_hits_target(
        "CC(=O)Nc1ccc(O)cc1",
        "CC(=O)N=C1C=CC(=O)C=C1",
        &phase1,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 40,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_benzene_to_bq_quinone_formation() {
    let qf = RuleSet::compose(Some("QF".into()), [quinone_formation()]);
    assert_path_hits_target(
        "c1ccccc1",
        "O=C1C=CC(=O)C=C1",
        &qf,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 40,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_benzene_to_bq_phase_one() {
    assert_path_hits_target(
        "c1ccccc1",
        "O=C1C=CC(=O)C=C1",
        &phase_one(),
        FindPathConfig {
            max_paths: 1,
            max_nodes: 80,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_terbinafine_to_tba_phase_one() {
    const TERBINAFINE: &str = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12";
    // Chematic keeps the reactant E stereo on the cleaved aldehyde; Python
    // RDKit canon of the non-stereo spelling matched. Assert csmi identity
    // against the stereo form the door actually emits.
    const TBF_A: &str = r"C(#C/C=C/C=O)C(C)(C)C";
    let mut counters = PathCounters::default();
    let want = canon_of(TBF_A).unwrap();
    let hits = find_path_with(
        TERBINAFINE,
        TBF_A,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 40,
            use_atom_diff: false,
        },
        |_| true,
    )
    .unwrap();
    assert!(
        !hits.is_empty(),
        "billed={} nodes={} edits={}",
        counters.billed(),
        counters.nodes,
        counters.mol_edits
    );
    assert_eq!(hits[0].smiles, want);
    assert_eq!(hits[0].steps[0].leaf_rule(), Some("Dealkylation"));
    assert!(counters.nodes <= 40);
}

// --- basic_ported find_path -------------------------------------------------

#[test]
fn gold_epoxide_opening_aromatic() {
    let eo = RuleSet::compose(Some("EO".into()), [epoxidation(), epoxide_opening()]);
    assert_path_hits_target(
        "c1ccccc1",
        "C1=CC=CC(O)C1O",
        &eo,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 200,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_epoxide_opening_kekulized() {
    let eo = RuleSet::compose(Some("EO".into()), [epoxidation(), epoxide_opening()]);
    assert_path_hits_target(
        "C1=CC=CC=C1",
        "C1=CC=CC(O)C1O",
        &eo,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 200,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_propane_dehydrogenation_to_propene() {
    let dh = RuleSet::compose(Some("DH".into()), [dehydrogenation()]);
    assert_path_hits_target(
        "CCC",
        "C=CC",
        &dh,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 50,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_propane_dh_then_epoxidation() {
    let dhe = RuleSet::compose(Some("DHE".into()), [dehydrogenation(), epoxidation()]);
    assert_path_hits_target(
        "CCC",
        "C1OC1C",
        &dhe,
        FindPathConfig {
            max_paths: 1,
            max_nodes: 200,
            use_atom_diff: false,
        },
    );
}

#[test]
fn gold_epoxide_opening_depth_one_misses_diol() {
    let eo = RuleSet::compose(Some("EO".into()), [epoxidation(), epoxide_opening()]);
    let mut counters = PathCounters::default();
    let hits = find_path_with(
        "c1ccccc1",
        "C1=CC=CC(O)C1O",
        &eo,
        &mut counters,
        FindPathConfig {
            max_paths: 3,
            max_nodes: 2,
            use_atom_diff: false,
        },
        |_| true,
    )
    .unwrap();
    assert!(hits.is_empty());
}

#[test]
fn gold_invalid_target_smiles_errors() {
    let eo = RuleSet::compose(Some("EO".into()), [epoxidation(), epoxide_opening()]);
    let mut counters = PathCounters::default();
    assert!(find_path("c1ccccc1", "C1=CC=CC1OC1", &eo, &mut counters).is_err());
}
