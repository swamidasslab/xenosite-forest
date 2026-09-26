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
    assert!(labs.contains(&"oxidation"), "{labs:?}");
    // Cross-cutting spines auto-tag alongside the chemist type.
    assert!(labs.contains(&"net oxidation"), "{labs:?}");
    assert!(labs.contains(&"oxygen gain"), "{labs:?}");
    assert!(labs.contains(&"stable oxygen addition"), "{labs:?}");
    assert!(labs.contains(&"single-metabolite transformation"), "{labs:?}");
    assert!(labs.contains(&"aliphatic site"), "{labs:?}");
    // Forest-map class parents are a parallel hierarchy, not on the chemist path.
    assert!(!labs.contains(&"stable oxygenation"), "{labs:?}");
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
fn parallel_hierarchies_exist() {
    let namer = namer();
    let forest_map = namer.thesaurus.get("xrm:9000000").expect("forest map");
    assert_eq!(forest_map.pref_label, "Metabolic Forest map");
    assert_eq!(
        forest_map.broader.first().map(|b| b.as_str()),
        Some("xrm:0000000")
    );

    let aromatic = namer.thesaurus.get("xrm:8000000").expect("aromatic impact");
    assert_eq!(
        aromatic.pref_label,
        "aromatic and conjugated-system impact"
    );

    // Forest-map class parents use full names and hang under the Forest map, not phase I.
    for (id, label) in [
        ("xrm:0000010", "stable oxygenation"),
        ("xrm:0000011", "unstable oxygenation"),
    ] {
        let c = namer.thesaurus.get(id).expect(id);
        assert_eq!(c.pref_label, label);
        assert!(
            c.broader.iter().any(|b| b.as_str() == "xrm:9000000"),
            "{id} broader={:?}",
            c.broader
        );
        assert!(
            !c.broader.iter().any(|b| b.as_str() == "xrm:0000001"),
            "{id} must not be under chemist phase I"
        );
        for alt in &c.alt_labels {
            assert!(
                !matches!(alt.as_str(), "SO" | "UO" | "DH" | "HD" | "RD"),
                "{id} must not use abbreviations as altLabel: {alt}"
            );
        }
    }

    // Chemist hydroxylation is under oxidation, related to Forest Hydroxylation rule.
    let oh = namer.thesaurus.get("xrm:0000100").unwrap();
    assert!(oh.broader.iter().any(|b| b.as_str() == "xrm:0000021"));
    let rule = namer.thesaurus.get("xrm:9100100").unwrap();
    assert_eq!(rule.pref_label, "Hydroxylation rule");
    assert!(rule.broader.iter().any(|b| b.as_str() == "xrm:0000010"));

    let forest = &namer.mappings.by_subject;
    assert!(forest["xrm:0000010"]
        .iter()
        .any(|m| m.object_id.as_str() == "forest.ruleset:SO"));
    assert!(forest["xrm:9100100"]
        .iter()
        .any(|m| m.object_id.as_str() == "forest.rule:Hydroxylation"));
}

#[test]
fn forest_tag_emits_forest_map_rule() {
    let terms = namer()
        .name_smiles("C", "C", &["forest.rule:Hydroxylation"])
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"Hydroxylation rule"), "{labs:?}");
    assert!(labs.contains(&"stable oxygenation"), "{labs:?}");
    assert!(labs.contains(&"Metabolic Forest map"), "{labs:?}");
}

#[test]
fn dearomatization_tag_reaches_aromatic_impact_hierarchy() {
    let terms = namer()
        .name_smiles("c1ccccc1", "C1=CC=CC=C1", &["chem:dearomatization"])
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"dearomatization"), "{labs:?}");
    assert!(labs.contains(&"aromaticity loss"), "{labs:?}");
    assert!(
        labs.contains(&"aromatic and conjugated-system impact"),
        "{labs:?}"
    );
}

#[test]
fn ambiguity_tags_stack_with_positive_terms() {
    let terms = namer()
        .name_smiles(
            "c1ccccc1",
            "Oc1ccccc1",
            &[
                "chem:aromatic-hydroxylation",
                "chem:regio-ambiguity",
                "chem:competing-type",
            ],
        )
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"aromatic hydroxylation"), "{labs:?}");
    assert!(labs.contains(&"regiochemical ambiguity"), "{labs:?}");
    assert!(labs.contains(&"competing-type ambiguity"), "{labs:?}");
    assert!(labs.contains(&"site-of-metabolism ambiguity"), "{labs:?}");
    assert!(
        labs.contains(&"ambiguity and underspecification"),
        "{labs:?}"
    );
}

#[test]
fn parallel_spine_roots_present() {
    let namer = namer();
    for (id, label) in [
        ("xrm:6000000", "ambiguity and underspecification"),
        ("xrm:7000000", "redox polarity"),
        ("xrm:7100000", "site atom class"),
        ("xrm:7200000", "bond-edit topology"),
        ("xrm:7300000", "metabolite cardinality"),
        ("xrm:7400000", "oxygenation outcome"),
        ("xrm:7500000", "electrophile role"),
        ("xrm:7600000", "ring fate"),
        ("xrm:7700000", "formula-delta class"),
        ("xrm:7800000", "site aromaticity"),
        ("xrm:7900000", "pathway-step role"),
        ("xrm:8000000", "aromatic and conjugated-system impact"),
        ("xrm:9000000", "Metabolic Forest map"),
    ] {
        let c = namer.thesaurus.get(id).unwrap_or_else(|| panic!("{id}"));
        assert_eq!(c.pref_label, label);
        assert!(
            c.broader.iter().any(|b| b.as_str() == "xrm:0000000"),
            "{id} should hang under root"
        );
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

#[test]
fn site_localized_tags_disambiguate_multi_change() {
    // Two methyls; only map 1 is asserted as N/A hydroxylation via localized tag.
    // Structural SMARTS also fire per mapped carbon when present.
    let terms = namer()
        .name_smiles(
            "[CH3:1]c1ccc([CH3:2])cc1",
            "O[CH2:1]c1ccc([CH3:2])cc1",
            &["chem:benzylic-hydroxylation@1"],
        )
        .unwrap();
    let benzylic: Vec<_> = terms
        .iter()
        .filter(|t| t.pref_label == "benzylic hydroxylation")
        .collect();
    assert_eq!(benzylic.len(), 1, "{:?}", labels(&terms));
    assert_eq!(benzylic[0].site.map_nums, vec![1]);

    // Localized chem tag alone pins site 2 (no structural SMARTS required).
    let terms = namer()
        .name_smiles(
            "[CH3:1]c1ccc([CH3:2])cc1",
            "[CH3:1]c1ccc([CH3:2])cc1",
            &["chem:aromatic-hydroxylation@2"],
        )
        .unwrap();
    let aryl: Vec<_> = terms
        .iter()
        .filter(|t| t.pref_label == "aromatic hydroxylation")
        .collect();
    assert_eq!(aryl.len(), 1, "{:?}", labels(&terms));
    assert_eq!(aryl[0].site.map_nums, vec![2]);
}
