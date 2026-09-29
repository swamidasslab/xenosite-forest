//! Sync forest [`Tag`] with chematic caller tags (`1..=u16::MAX`).
//!
//! Chematic: `set_tag(..., None)` or `Some(0)` clears; `atom_tag` returns only
//! non-zero labels. Our [`Tag`] is always non-zero, so clear is `Option::None`.
//!
//! Some chematic rebuilds omit `copy_atom_tags_from`. Prefer [`preserving`]
//! wrappers for index-stable edits rather than remapping after the fact.
//! See also [`crate::chematic_features`] for chematic-chem pick stages.

use chematic::core::{AtomIdx, Element, Molecule};

use crate::labels::Tag;
use crate::mol::{aromatize as aromatize_raw, atom_idx};

/// Push a sidecar label onto chematic (`None` clears; never writes `0`).
pub fn set_label(mol: &mut Molecule, idx: AtomIdx, tag: Option<Tag>) {
    mol.set_tag(idx, tag.map(Tag::get));
}

/// Read a chematic tag (`None` if missing or cleared).
pub fn get_label(mol: &Molecule, idx: AtomIdx) -> Option<Tag> {
    mol.atom_tag(idx).map(Tag::from_nonzero)
}

fn copy_all_labels(from: &Molecule, to: &mut Molecule) {
    let n = from.atom_count().min(to.atom_count());
    for i in 0..n {
        set_label(to, atom_idx(i), get_label(from, atom_idx(i)));
    }
}

/// Index-stable chematic mutations that preserve caller tags.
///
/// Valid only when atom count and order are unchanged (or new atoms are
/// appended). Do not use after SMIRKS reorder — use isotope/map tracing there.
pub mod preserving {
    use super::*;

    pub fn with_atom_aromatic(mol: &Molecule, idx: AtomIdx, aromatic: bool) -> Molecule {
        let mut out = mol.with_atom_aromatic(idx, aromatic);
        copy_all_labels(mol, &mut out);
        out
    }

    pub fn with_atom_charge(mol: &Molecule, idx: AtomIdx, charge: i8) -> Molecule {
        let mut out = mol.with_atom_charge(idx, charge);
        copy_all_labels(mol, &mut out);
        out
    }

    pub fn with_atom_element(mol: &Molecule, idx: AtomIdx, element: Element) -> Molecule {
        let mut out = mol.with_atom_element(idx, element);
        copy_all_labels(mol, &mut out);
        out
    }

    pub fn with_atom_explicit_h(mol: &Molecule, idx: AtomIdx, h: u8) -> Molecule {
        let mut out = crate::kekule::with_atom_explicit_h(mol, idx, h);
        copy_all_labels(mol, &mut out);
        out
    }

    /// RDKit-parity aromatize, then restore tags by index (perception rebuild).
    pub fn aromatize(mol: &Molecule) -> Molecule {
        let mut out = aromatize_raw(mol);
        copy_all_labels(mol, &mut out);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{atom_idx, parse_mol};
    use chematic::core::{Atom, Element};

    fn stamp(mol: &mut Molecule) -> Vec<Option<Tag>> {
        for i in 0..mol.atom_count() {
            set_label(mol, atom_idx(i), Some(Tag::new((i as u16) + 1).unwrap()));
        }
        tags_of(mol)
    }

    fn tags_of(mol: &Molecule) -> Vec<Option<Tag>> {
        (0..mol.atom_count())
            .map(|i| get_label(mol, atom_idx(i)))
            .collect()
    }

    fn zs(mol: &Molecule) -> Vec<u8> {
        (0..mol.atom_count())
            .map(|i| mol.atom(atom_idx(i)).element.atomic_number())
            .collect()
    }

    /// Index-copy is only valid when heavy-atom order is unchanged.
    fn assert_same_atom_order(before: &Molecule, after: &Molecule) {
        assert_eq!(
            before.atom_count(),
            after.atom_count(),
            "atom count changed — cannot copy tags by index"
        );
        assert_eq!(
            zs(before),
            zs(after),
            "element order changed — index tag copy would mis-assign"
        );
    }

    // --- Correct behavior (green when the op already keeps tags / order) ---

    #[test]
    fn with_atom_aromatic_preserves_tags_and_order() {
        let mut mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let before = stamp(&mut mol);
        let out = mol.with_atom_aromatic(atom_idx(0), false);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    fn with_atom_charge_preserves_tags_and_order() {
        let mut mol = parse_mol("CCN").unwrap();
        let before = stamp(&mut mol);
        let out = mol.with_atom_charge(atom_idx(2), 1);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
        assert_eq!(out.atom(atom_idx(2)).charge, 1);
    }

    #[test]
    fn with_atom_element_preserves_tags_and_atom_indexes() {
        let mut mol = parse_mol("CCl").unwrap();
        let before = stamp(&mut mol);
        let out = mol.with_atom_element(atom_idx(1), Element::O);
        assert_eq!(mol.atom_count(), out.atom_count());
        assert_eq!(tags_of(&out), before);
        assert_eq!(out.atom(atom_idx(1)).element.atomic_number(), 8);
        assert_eq!(zs(&mol)[0], zs(&out)[0]);
    }

    #[test]
    fn with_atom_added_preserves_existing_tags() {
        let mut mol = parse_mol("CC").unwrap();
        let before = stamp(&mut mol);
        let (next, _) = mol.with_atom_added(Atom::organic(Element::O));
        assert_eq!(&tags_of(&next)[..before.len()], before.as_slice());
    }

    #[test]
    fn fragments_preserves_tags() {
        let mut mol = parse_mol("CC.OO").unwrap();
        let before = stamp(&mut mol);
        let frags = mol.fragments();
        assert_eq!(frags.len(), 2);
        let kept: usize = frags.iter().map(|f| tags_of(f).iter().filter(|t| t.is_some()).count()).sum();
        assert_eq!(kept, before.iter().filter(|t| t.is_some()).count());
    }

    #[test]
    fn preserving_with_atom_aromatic_keeps_tags_and_order() {
        let mut mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::with_atom_aromatic(&mol, atom_idx(1), false);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
        assert!(!out.atom(atom_idx(1)).aromatic);
    }

    #[test]
    fn preserving_with_atom_charge_keeps_tags_and_order() {
        let mut mol = parse_mol("CCN").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::with_atom_charge(&mol, atom_idx(2), 1);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    fn preserving_with_atom_element_keeps_tags() {
        let mut mol = parse_mol("CCl").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::with_atom_element(&mol, atom_idx(1), Element::O);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    #[ignore = "xfail: chematic — MoleculeBuilder without copy_atom_tags_from drops caller tags; keep to track the pitfall (Forest kekule::with_atom_explicit_h copies)"]
    fn molecule_builder_without_copy_atom_tags_drops_tags() {
        use chematic::core::MoleculeBuilder;
        let mut mol = parse_mol("CCO").unwrap();
        let before = stamp(&mut mol);
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in mol.atoms() {
            builder.add_atom(atom.clone());
        }
        for (_, bond) in mol.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        // Deliberately omit copy_atom_tags_from — the chematic-chem pitfall.
        let out = builder.build();
        assert_same_atom_order(&mol, &out);
        assert_eq!(
            tags_of(&out),
            before,
            "desired: builder carries tags without an explicit copy"
        );
    }

    #[test]
    fn with_atom_explicit_h_preserves_tags_and_order() {
        let mut mol = parse_mol("CCO").unwrap();
        let before = stamp(&mut mol);
        let out = crate::kekule::with_atom_explicit_h(&mol, atom_idx(0), 2);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    fn preserving_with_atom_explicit_h_keeps_tags_and_order() {
        let mut mol = parse_mol("CCO").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::with_atom_explicit_h(&mol, atom_idx(0), 2);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    fn preserving_aromatize_keeps_tags_and_order() {
        let mut mol = parse_mol("C1=CC=CC=C1").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::aromatize(&mol);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    #[test]
    fn preserving_aromatize_on_already_aromatic_keeps_tags() {
        let mut mol = parse_mol("c1ccccc1").unwrap();
        let before = stamp(&mut mol);
        let out = preserving::aromatize(&mol);
        assert_same_atom_order(&mol, &out);
        assert_eq!(tags_of(&out), before);
    }

    // --- Desired behavior; currently broken → xfail until fixed ---

    #[test]
    fn raw_aromatize_preserves_tags_and_order() {
        let mol = parse_mol("C1=CC=CC=C1").unwrap();
        // parse_mol already aromatized; force kekulé then raw aromatize.
        let mut kek = mol.clone();
        for i in 0..kek.atom_count() {
            if kek.atom(atom_idx(i)).aromatic {
                kek = kek.with_atom_aromatic(atom_idx(i), false);
            }
        }
        let before = stamp(&mut kek);
        let out = aromatize_raw(&kek);
        assert_same_atom_order(&kek, &out);
        assert_eq!(
            tags_of(&out),
            before,
            "raw aromatize must keep caller tags without Forest restamp"
        );
    }

    #[test]
    #[ignore = "xfail: chematic — ResonancePair materialize products lose caller tags (rebuild path); keep ignored to track until chematic/path carries tags"]
    fn pair_materialize_preserves_parent_tags_on_hydroquinone() {
        use crate::ForestMol;
        use crate::rules::dehydrogenation;
        let parent = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let parent_tags: Vec<_> = (0..parent.mol().atom_count())
            .map(|i| parent.tag_of(i))
            .collect();
        assert!(parent_tags.iter().all(|t| t.is_some()));
        let pairs = dehydrogenation()
            .candidates(&parent)
            .filter(|c| matches!(c, Ok(s) if s.is_pair()))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!pairs.is_empty());
        let pieces = pairs[0].materialize_mols().unwrap();
        assert!(!pieces.is_empty());
        let piece_tags = tags_of(&pieces[0]);
        assert!(
            piece_tags.iter().any(|t| t.is_some()),
            "pair materialize must keep caller tags: {piece_tags:?}"
        );
        for t in piece_tags.iter().flatten() {
            assert!(
                parent_tags.iter().flatten().any(|p| p == t),
                "product tag {t:?} not from parent {parent_tags:?}"
            );
        }
        let child = parent.from_edit_product(pieces[0].clone());
        for t in parent_tags.iter().flatten() {
            assert!(
                child.index_of(*t).is_some(),
                "parent tag {t:?} missing on adopted child"
            );
        }
    }
}
