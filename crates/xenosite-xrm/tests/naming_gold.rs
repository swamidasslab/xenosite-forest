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
    assert!(labs.contains(&"metabolism phase"), "{labs:?}");
    assert!(labs.contains(&"oxidation"), "{labs:?}");
    // Rainbow phase I family + medchem liability (ChatGPT spines).
    assert!(labs.contains(&"stable oxygenation"), "{labs:?}");
    assert!(labs.contains(&"phase I reaction family"), "{labs:?}");
    assert!(labs.contains(&"metabolic soft spot"), "{labs:?}");
    // Nested structural-delta facets still auto-tag.
    assert!(labs.contains(&"net oxidation"), "{labs:?}");
    assert!(labs.contains(&"oxygen gain"), "{labs:?}");
    assert!(labs.contains(&"stable oxygen addition"), "{labs:?}");
    assert!(labs.contains(&"single-metabolite transformation"), "{labs:?}");
    assert!(labs.contains(&"aliphatic site"), "{labs:?}");
    // Forest ruleset alias is not on the pure structural path.
    assert!(!labs.contains(&"stable oxygenation ruleset"), "{labs:?}");
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
    // Enzyme names are allowed only under biological context (orthogonal spine).
    let bio = "xrm:2100000";
    for c in namer.thesaurus.concepts.values() {
        let under_bio = c.broader.iter().any(|b| b.as_str() == bio)
            || c.id.as_str() == bio
            || namer
                .thesaurus
                .ancestors(c.id.as_str())
                .iter()
                .any(|a| a.as_str() == bio);
        if under_bio {
            continue;
        }
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

    // Rainbow chemist classes hang under phase I reaction family, not Forest map.
    for (id, label) in [
        ("xrm:0000010", "stable oxygenation"),
        ("xrm:0000011", "unstable oxygenation"),
    ] {
        let c = namer.thesaurus.get(id).expect(id);
        assert_eq!(c.pref_label, label);
        assert!(
            c.broader.iter().any(|b| b.as_str() == "xrm:1200000"),
            "{id} broader={:?}",
            c.broader
        );
        assert!(
            !c.broader.iter().any(|b| b.as_str() == "xrm:9000000"),
            "{id} must not be under Forest map"
        );
        for alt in &c.alt_labels {
            assert!(
                !matches!(alt.as_str(), "SO" | "UO" | "DH" | "HD" | "RD"),
                "{id} must not use abbreviations as altLabel: {alt}"
            );
        }
    }

    // Forest ruleset aliases use unabbreviated *ruleset* labels.
    for (id, label) in [
        ("xrm:9000010", "stable oxygenation ruleset"),
        ("xrm:9000011", "unstable oxygenation ruleset"),
    ] {
        let c = namer.thesaurus.get(id).expect(id);
        assert_eq!(c.pref_label, label);
        assert!(c.broader.iter().any(|b| b.as_str() == "xrm:9000000"));
    }

    // Chemist hydroxylation is under carbon oxidation + Rainbow stable oxygenation.
    let oh = namer.thesaurus.get("xrm:0000100").unwrap();
    assert!(oh.broader.iter().any(|b| b.as_str() == "xrm:0000021"));
    assert!(oh.broader.iter().any(|b| b.as_str() == "xrm:0000010"));
    let rule = namer.thesaurus.get("xrm:9100100").unwrap();
    assert_eq!(rule.pref_label, "Hydroxylation rule");
    assert!(rule.broader.iter().any(|b| b.as_str() == "xrm:9000010"));

    let forest = &namer.mappings.by_subject;
    assert!(forest["xrm:9000010"]
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
    assert!(labs.contains(&"stable oxygenation ruleset"), "{labs:?}");
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
        ("xrm:1000000", "metabolism phase"),
        ("xrm:1100000", "chemical transformation"),
        ("xrm:1200000", "phase I reaction family"),
        ("xrm:1300000", "phase II conjugation family"),
        ("xrm:1400000", "medchem liability"),
        ("xrm:1500000", "reactive metabolite family"),
        ("xrm:1600000", "site type"),
        ("xrm:1700000", "structural delta"),
        ("xrm:1800000", "product status"),
        ("xrm:1900000", "rule provenance"),
        ("xrm:2000000", "evidence"),
        ("xrm:2100000", "biological context"),
        ("xrm:2200000", "leaving group"),
        ("xrm:6000000", "ambiguity and underspecification"),
        ("xrm:9000000", "Metabolic Forest map"),
    ] {
        let c = namer.thesaurus.get(id).unwrap_or_else(|| panic!("{id}"));
        assert_eq!(c.pref_label, label);
        assert!(
            c.broader.iter().any(|b| b.as_str() == "xrm:0000000"),
            "{id} should hang under root"
        );
    }
    // Legacy facets nest under the ChatGPT spines (not peer roots).
    let redox = namer.thesaurus.get("xrm:7000000").unwrap();
    assert!(redox.broader.iter().any(|b| b.as_str() == "xrm:1700000"));
    let aromatic = namer.thesaurus.get("xrm:8000000").unwrap();
    assert!(aromatic.broader.iter().any(|b| b.as_str() == "xrm:1700000"));
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

#[test]
fn leaving_group_tags_with_dealkylation() {
    let terms = namer()
        .name_smiles("C", "C", &["chem:N-demethylation"])
        .unwrap();
    let labs = labels(&terms);
    assert!(labs.contains(&"N-demethylation"), "{labs:?}");
    assert!(labs.contains(&"methyl leaving group"), "{labs:?}");
    assert!(labs.contains(&"leaving group"), "{labs:?}");
}

#[test]
fn annotation_bundle_facets_without_combinatorial_concepts() {
    let n = namer();
    let bundles = n.annotate_smiles("CC", "CCO", &[]).unwrap();
    assert!(!bundles.is_empty(), "{bundles:?}");
    // At least one bundle carries phase I + a transformation + structural delta.
    let hit = bundles.iter().find(|b| {
        b.phase.as_deref() == Some("xrm:0000001")
            && b.transformation.is_some()
            && !b.structural_delta.is_empty()
    });
    assert!(hit.is_some(), "{bundles:?}");
    let b = hit.unwrap();
    assert!(
        b.xmet_properties.iter().any(|p| p == "xmet:hasPhase"),
        "{:?}",
        b.xmet_properties
    );
    assert!(
        b.xmet_properties
            .iter()
            .any(|p| p == "xmet:hasTransformation"),
        "{:?}",
        b.xmet_properties
    );
    // site_label is a display string, not a minted concept id.
    let label = b.site_label.as_deref().unwrap_or("");
    assert!(!label.is_empty(), "{b:?}");
    assert!(!label.starts_with("xrm:"), "{label}");
    assert!(
        !n.thesaurus
            .concepts
            .values()
            .any(|c| c.pref_label == "benzylic_phase_1_hydroxylation_clearance_liability")
    );
}

#[test]
fn annotation_bundle_site_label_uses_map_template() {
    let bundles = namer()
        .annotate_smiles(
            "[CH3:1]c1ccc([CH3:2])cc1",
            "O[CH2:1]c1ccc([CH3:2])cc1",
            &["chem:benzylic-hydroxylation@1"],
        )
        .unwrap();
    let localized = bundles
        .iter()
        .find(|b| b.site.map_nums == vec![1])
        .expect("site @1 bundle");
    let label = localized.site_label.as_deref().unwrap();
    assert!(label.contains('@') || label.contains("benzylic"), "{label}");
    assert!(
        localized
            .xmet_properties
            .iter()
            .any(|p| p == "xmet:localizesToSite"),
        "{:?}",
        localized.xmet_properties
    );
}
