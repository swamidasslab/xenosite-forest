//! Thesaurus / mapping / assignment schema tests.

use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

#[test]
fn default_manifest_loads_and_validates() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).expect("load default manifest");
    assert!(namer.thesaurus.concepts.contains_key("xrm:0000001"));
    assert!(namer.thesaurus.concepts.contains_key("xrm:0000002"));
    namer.thesaurus.validate().unwrap();
}

#[test]
fn phase_i_and_ii_present_with_mesh_crosswalks() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    let p1 = namer.thesaurus.get("xrm:0000001").unwrap();
    assert_eq!(p1.pref_label, "phase I");
    let p2 = namer.thesaurus.get("xrm:0000002").unwrap();
    assert_eq!(p2.pref_label, "phase II");

    let mesh = &namer.mappings.by_subject;
    assert!(mesh["xrm:0000001"]
        .iter()
        .any(|m| m.object_id.as_str() == "mesh:D050216"));
    assert!(mesh["xrm:0000002"]
        .iter()
        .any(|m| m.object_id.as_str() == "mesh:D050217"));
}

#[test]
fn hydroxylation_has_mop_exact_match() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    assert!(namer.mappings.by_subject["xrm:0000100"]
        .iter()
        .any(|m| m.object_id.as_str() == "mop:0000673"));
}

#[test]
fn forest_mappings_are_opaque_curies_only() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    for rows in namer.mappings.by_subject.values() {
        for m in rows {
            if let Some(p) = m.object_id.prefix() {
                if p.starts_with("forest") {
                    assert!(
                        m.object_id.as_str().contains(':'),
                        "forest CURIE must stay opaque: {}",
                        m.object_id
                    );
                }
            }
        }
    }
}

#[test]
fn assignment_emit_ids_exist_in_thesaurus() {
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).unwrap();
    for rule in &namer.assignments.rules {
        for id in &rule.emit {
            assert!(
                namer.thesaurus.get(id).is_some(),
                "assignment {} emits unknown {}",
                rule.id,
                id
            );
        }
    }
}
