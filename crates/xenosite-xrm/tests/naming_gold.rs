//! Golden reactant→product naming cases.

use xenosite_xrm::{MappedReaction, Namer, DEFAULT_MANIFEST};

fn namer() -> Namer {
    Namer::from_manifest(DEFAULT_MANIFEST).unwrap()
}

fn labels(terms: &[xenosite_xrm::Term]) -> Vec<&str> {
    terms.iter().map(|t| terms_label(t)).collect()
}

fn terms_label(t: &xenosite_xrm::Term) -> &str {
    t.pref_label.as_str()
}

#[test]
fn aliphatic_hydroxylation_ethane_to_ethanol() {
    let terms = namer().name_smiles("CC", "CCO", &[]).unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"aliphatic hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"phase I"), "{labs:?}");
    assert!(labs.contains(&"stable oxygenation"), "{labs:?}");
    assert!(!labs.contains(&"aromatic hydroxylation"), "{labs:?}");
    assert!(terms[0].specificity.depth >= terms.last().unwrap().specificity.depth);
}

#[test]
fn aromatic_hydroxylation_benzene_to_phenol() {
    let terms = namer().name_smiles("c1ccccc1", "Oc1ccccc1", &[]).unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"aromatic hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"phase I"), "{labs:?}");
}

#[test]
fn epoxidation_ethene() {
    let terms = namer().name_smiles("C=C", "C1CO1", &[]).unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"epoxidation"), "{labs:?}");
    assert!(labs.contains(&"phase I"), "{labs:?}");
    assert!(
        !labs.contains(&"hydroxylation"),
        "epoxide must not be labeled hydroxylation: {labs:?}"
    );
}

#[test]
fn alcohol_dehydrogenation() {
    let terms = namer().name_smiles("CCO", "CC=O", &[]).unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"dehydrogenation"), "{labs:?}");
    assert!(labs.contains(&"phase I"), "{labs:?}");
}

#[test]
fn glucuronidation_via_opaque_tag() {
    // Structural conjugation SMARTS are not required; opaque tags are sufficient.
    let terms = namer()
        .name_smiles("CCO", "CCO", &["forest.rule:Glucuronidation"])
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"glucuronidation"), "{labs:?}");
    assert!(labs.contains(&"phase II"), "{labs:?}");
}

#[test]
fn gsh_michael_via_pattern_tag() {
    let terms = namer()
        .name(
            &MappedReaction::new("C=CC=O", "C=CC=O")
                .with_tags(["forest.pattern:Glutathionation/michael"]),
        )
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"Michael glutathionation"), "{labs:?}");
    assert!(labs.contains(&"glutathionation"), "{labs:?}");
    assert!(labs.contains(&"phase II"), "{labs:?}");
}

#[test]
fn terms_carry_ontology_and_inter_matches() {
    let terms = namer().name_smiles("CC", "CCO", &[]).unwrap();
    let oh = terms
        .iter()
        .find(|t| t.pref_label == "hydroxylation")
        .expect("hydroxylation term");
    assert_eq!(oh.ontology.id.as_str(), "xrm:scheme");
    assert!(oh.ontology.label.contains("Xenobiotic"));
    assert!(
        oh.inter_matches
            .iter()
            .any(|l| l.target.as_str() == "mop:0000673"),
        "{:?}",
        oh.inter_matches
    );
    assert!(
        oh.inter_matches
            .iter()
            .any(|l| l.target.as_str() == "forest.rule:Hydroxylation"),
        "forest mapping should appear as inter-ontology relatedMatch"
    );
    assert!(!oh.path_labels.is_empty());
    assert!(oh.broader.iter().any(|b| b.predicate == "broader"));
}

#[test]
fn no_enzyme_strings_in_primary_labels() {
    let namer = namer();
    for c in namer.thesaurus.concepts.values() {
        let l = c.pref_label.to_ascii_lowercase();
        assert!(!l.contains("cyp"), "{}", c.pref_label);
        assert!(!l.contains("cytochrome"), "{}", c.pref_label);
        assert!(!l.starts_with("ec "), "{}", c.pref_label);
    }
}

#[test]
fn tag_hydroxylation_without_structural_match_still_names() {
    // Same formula pair that is NOT a hydroxylation chemically, but caller
    // supplied an opaque forest tag — namer trusts tags as data.
    let terms = namer()
        .name_smiles("C", "C", &["forest.rule:Hydroxylation"])
        .unwrap();
    assert!(labels(&terms).contains(&"hydroxylation"));
}
