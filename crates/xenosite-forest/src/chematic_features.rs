//! Chematic upstream feature probes (tag carry through chem rebuilds).
//!
//! Chematic-chem rebuilds (`MoleculeBuilder::new` without
//! `copy_atom_tags_from`) drop caller-tag sidecars on
//! [`normalize_zwitterion`](chematic::chem::normalize_zwitterion),
//! [`remove_hydrogens`](chematic::chem::remove_hydrogens), and
//! [`canonical_tautomer`](chematic::chem::canonical_tautomer). Forest
//! [`crate::normalize`] snapshots survivor tags and restamps once after the
//! pick until that is fixed upstream.
//!
//! These tests assert the **desired** chematic behavior and are
//! `#[ignore]` (= xfail). When they pass without ignore:
//!
//! 1. Drop the begin/end restamp in [`crate::normalize::chematic_tautomer_pick`].
//! 2. Remove the `#[ignore]` attributes here (or delete this module).
//!
//! Run: `cargo test -p xenosite-forest chematic_features -- --ignored`

use chematic::chem::{canonical_tautomer, normalize_zwitterion, remove_hydrogens};
use chematic::core::Element;
use chematic::smiles::parse;

use crate::chematic_tags::{get_label, set_label};
use crate::labels::Tag;
use crate::mol::{Molecule, atom_idx};

fn stamp_unique_tags(mol: &mut Molecule) {
    for i in 0..mol.atom_count() {
        set_label(mol, atom_idx(i), Some(Tag::new(100 + i as u16).unwrap()));
    }
}

fn tags_in_order(mol: &Molecule) -> Vec<Option<u16>> {
    (0..mol.atom_count())
        .map(|i| get_label(mol, atom_idx(i)).map(Tag::get))
        .collect()
}

fn survivor_tags(mol: &Molecule) -> Vec<Option<u16>> {
    (0..mol.atom_count())
        .filter(|&i| {
            let a = mol.atom(atom_idx(i));
            !(a.element == Element::H && a.isotope.is_none())
        })
        .map(|i| get_label(mol, atom_idx(i)).map(Tag::get))
        .collect()
}

fn any_caller_tag(mol: &Molecule) -> bool {
    (0..mol.atom_count()).any(|i| get_label(mol, atom_idx(i)).is_some())
}

/// Desired: chematic stages keep caller tags. Ignored until upstream does.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "chematic-chem rebuilds drop atom_tags (no copy_atom_tags_from); remove ignore + Forest restamp when fixed"]
    fn remove_hydrogens_preserves_survivor_tags() {
        let mut mol = parse("[H]OC([H])=C([H])[H]").unwrap();
        stamp_unique_tags(&mut mol);
        let before = survivor_tags(&mol);
        assert!(before.iter().any(|t| t.is_some()), "stamped");
        let out = remove_hydrogens(&mol);
        assert_eq!(tags_in_order(&out), before);
    }

    #[test]
    #[ignore = "chematic-chem rebuilds drop atom_tags (no copy_atom_tags_from); remove ignore + Forest restamp when fixed"]
    fn normalize_zwitterion_preserves_tags_on_rebuild() {
        let mut mol = parse("[NH3+]CC(=O)[O-]").unwrap();
        stamp_unique_tags(&mut mol);
        let before = tags_in_order(&mol);
        let out = normalize_zwitterion(&mol);
        assert_ne!(
            chematic::smiles::canonical_smiles(&mol),
            chematic::smiles::canonical_smiles(&out),
            "fixture must rebuild (not a no-op clone)"
        );
        assert_eq!(out.atom_count(), mol.atom_count());
        assert_eq!(tags_in_order(&out), before);
    }

    #[test]
    #[ignore = "chematic-chem rebuilds drop atom_tags (no copy_atom_tags_from); remove ignore + Forest restamp when fixed"]
    fn canonical_tautomer_preserves_tags_on_rebuild() {
        let mut mol = parse("OC=C").unwrap();
        stamp_unique_tags(&mut mol);
        let before = tags_in_order(&mol);
        let out = canonical_tautomer(&mol);
        assert_ne!(
            chematic::smiles::canonical_smiles(&mol),
            chematic::smiles::canonical_smiles(&out),
            "fixture must rebuild (enol → keto)"
        );
        assert_eq!(out.atom_count(), mol.atom_count());
        assert_eq!(tags_in_order(&out), before);
    }

    #[test]
    #[ignore = "chematic-chem rebuilds drop atom_tags (no copy_atom_tags_from); remove ignore + Forest restamp when fixed"]
    fn full_pick_stream_preserves_survivor_tags_without_forest_restamp() {
        let mut mol = parse("[H]OC([H])=C([H])[H]").unwrap();
        stamp_unique_tags(&mut mol);
        let before = survivor_tags(&mol);
        let mol = normalize_zwitterion(&mol);
        let mol = remove_hydrogens(&mol);
        let out = canonical_tautomer(&mol);
        assert!(
            any_caller_tag(&out),
            "upstream pick must keep at least one caller tag"
        );
        assert_eq!(tags_in_order(&out), before);
    }
}
