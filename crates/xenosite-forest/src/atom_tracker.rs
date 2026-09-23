//! Atom tracker POC on chematic `Atom.tag` + SMILES visit order.
//!
//! Stamps a unique integer on each atom. Apply copies it (patch); born atoms
//! are `None` until [`AtomTracker::adopt_born`] mints ids. Write/parse uses
//! [`write_with_order`] so tags ride the visit permutation onto the new mol
//! without isotopes or a sidecar.
//!
//! This is the production-shaped seam once Forest stops carrying
//! [`crate::labels`] beside the graph. `ForestMol`'s sidecar stays for now.

use chematic::core::Molecule;
use chematic::smiles::{parse, write_with_order};

use crate::labels::Tag;
use crate::mol::{ForestError, atom_idx};

/// Mints and follows chematic atom tags through rewrite.
#[derive(Debug, Clone)]
pub struct AtomTracker {
    next: u32,
}

impl AtomTracker {
    /// Stamp every atom with a fresh tag. Existing tags are overwritten.
    pub fn stamp(mol: &mut Molecule) -> Self {
        let n = mol.atom_count();
        for i in 0..n {
            mol.set_tag(atom_idx(i), Some(i as u32));
        }
        Self { next: n as u32 }
    }

    pub fn next_tag(&self) -> u32 {
        self.next
    }

    pub fn tag_of(mol: &Molecule, idx: usize) -> Option<Tag> {
        mol.atom(atom_idx(idx)).tag.map(Tag)
    }

    pub fn index_of(mol: &Molecule, tag: Tag) -> Option<usize> {
        (0..mol.atom_count()).find(|&i| mol.atom(atom_idx(i)).tag == Some(tag.0))
    }

    /// Mint tags for atoms that still have `None` (born in apply / fragments).
    pub fn adopt_born(&mut self, mol: &mut Molecule) {
        for i in 0..mol.atom_count() {
            let idx = atom_idx(i);
            if mol.atom(idx).tag.is_none() {
                mol.set_tag(idx, Some(self.next));
                self.next += 1;
            }
        }
    }

    /// Parent index → product index for atoms whose tags survived.
    ///
    /// `out[src] = Some(dst)` when parent's tag at `src` appears on the
    /// product; `None` if that atom left this piece. Born product atoms are
    /// not represented here — find them as product indexes with no parent.
    pub fn src_to_new(parent: &Molecule, product: &Molecule) -> Vec<Option<usize>> {
        (0..parent.atom_count())
            .map(|src| Self::tag_of(parent, src).and_then(|tag| Self::index_of(product, tag)))
            .collect()
    }

    /// Non-canonical write, then parse, with tags remapped by visit order.
    ///
    /// `order[k]` is the old mol index of the k-th atom in the SMILES string;
    /// after parse, new index `k` holds that atom's tag.
    pub fn write_parse(mol: &Molecule) -> Result<(String, Molecule), ForestError> {
        let (smi, order) = write_with_order(mol);
        let mut fresh = parse(&smi).map_err(|err| ForestError::Parse(err.to_string()))?;
        if fresh.atom_count() != order.len() {
            return Err(ForestError::Parse(format!(
                "write/parse atom count mismatch: order {} vs parse {}",
                order.len(),
                fresh.atom_count()
            )));
        }
        for (new_i, &old_idx) in order.iter().enumerate() {
            let tag = mol.atom(old_idx).tag;
            fresh.set_tag(atom_idx(new_i), tag);
        }
        Ok((smi, fresh))
    }
}

/// True when every tagged parent atom that still exists has the same element.
pub fn tags_agree_elements(parent: &Molecule, product: &Molecule) -> bool {
    for src in 0..parent.atom_count() {
        let Some(tag) = AtomTracker::tag_of(parent, src) else {
            continue;
        };
        let Some(dst) = AtomTracker::index_of(product, tag) else {
            continue;
        };
        if parent.atom(atom_idx(src)).element.atomic_number()
            != product.atom(atom_idx(dst)).element.atomic_number()
        {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use chematic::core::{Atom, BondOrder, Element, MoleculeBuilder};
    use chematic::rxn::{apply_reaction_match, find_reaction_matches};
    use chematic::smiles::{parse, write_with_order};

    use super::*;
    use crate::mol::atom_usize;

    #[test]
    fn stamp_sets_unique_tags() {
        let mut mol = parse("CCO").unwrap();
        let tracker = AtomTracker::stamp(&mut mol);
        assert_eq!(tracker.next_tag(), 3);
        assert_eq!(AtomTracker::tag_of(&mol, 0), Some(Tag(0)));
        assert_eq!(AtomTracker::tag_of(&mol, 2), Some(Tag(2)));
        assert_eq!(AtomTracker::index_of(&mol, Tag(1)), Some(1));
    }

    #[test]
    fn follow_apply_keeps_carbon_tags_and_mints_born_oxygen() {
        let mut parent = parse("CCO").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let c0 = AtomTracker::tag_of(&parent, 0).unwrap();
        let c1 = AtomTracker::tag_of(&parent, 1).unwrap();
        let o = AtomTracker::tag_of(&parent, 2).unwrap();

        let matches = find_reaction_matches("[C:1]>>[C:1]O", &[&parent]).unwrap();
        let products = apply_reaction_match("[C:1]>>[C:1]O", &[&parent], &matches[0], true)
            .unwrap()
            .expect("valence");
        let mut product = products[0].clone();

        let map = AtomTracker::src_to_new(&parent, &product);
        // All three reactant atoms still present somewhere.
        assert!(map.iter().all(|d| d.is_some()));
        assert_eq!(AtomTracker::tag_of(&product, map[0].unwrap()), Some(c0));
        assert_eq!(AtomTracker::tag_of(&product, map[1].unwrap()), Some(c1));
        assert_eq!(AtomTracker::tag_of(&product, map[2].unwrap()), Some(o));

        let born: Vec<_> = (0..product.atom_count())
            .filter(|&i| AtomTracker::tag_of(&product, i).is_none())
            .collect();
        assert_eq!(born.len(), 1);
        assert_eq!(product.atom(atom_idx(born[0])).element.atomic_number(), 8);

        tracker.adopt_born(&mut product);
        assert!(AtomTracker::tag_of(&product, born[0]).is_some());
        assert_eq!(tracker.next_tag(), 4);
        assert!(tags_agree_elements(&parent, &product));
    }

    #[test]
    fn write_parse_remaps_tags_by_visit_order() {
        // Scrambled ethanol: O=0, methyl=1, center=2 → visit [0,2,1].
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let me = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::C));
        b.add_bond(c, o, BondOrder::Single).unwrap();
        b.add_bond(c, me, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let _tracker = AtomTracker::stamp(&mut mol);

        let (_smi, order) = write_with_order(&mol);
        let idxs: Vec<_> = order.iter().map(|a| a.0).collect();
        assert_ne!(idxs, vec![0, 1, 2]);

        let tag_at_old: Vec<_> = (0..mol.atom_count())
            .map(|i| AtomTracker::tag_of(&mol, i).unwrap())
            .collect();

        let (_smi, fresh) = AtomTracker::write_parse(&mol).unwrap();
        assert_eq!(fresh.atom_count(), mol.atom_count());
        for (new_i, &old) in order.iter().enumerate() {
            assert_eq!(
                AtomTracker::tag_of(&fresh, new_i).unwrap(),
                tag_at_old[atom_usize(old)]
            );
        }
        // Lookup by tag still finds the same element after reparse.
        for tag in [Tag(0), Tag(1), Tag(2)] {
            let old_i = AtomTracker::index_of(&mol, tag).unwrap();
            let new_i = AtomTracker::index_of(&fresh, tag).unwrap();
            assert_eq!(
                mol.atom(atom_idx(old_i)).element.atomic_number(),
                fresh.atom(atom_idx(new_i)).element.atomic_number()
            );
        }
    }

    #[test]
    fn apply_then_write_parse_preserves_identity() {
        let mut parent = parse("CC").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let matches = find_reaction_matches("[C:1]>>[C:1]O", &[&parent]).unwrap();
        let products = apply_reaction_match("[C:1]>>[C:1]O", &[&parent], &matches[0], true)
            .unwrap()
            .expect("valence");
        let mut product = products[0].clone();
        tracker.adopt_born(&mut product);

        let before: Vec<_> = (0..product.atom_count())
            .map(|i| {
                (
                    AtomTracker::tag_of(&product, i).unwrap(),
                    product.atom(atom_idx(i)).element.atomic_number(),
                )
            })
            .collect();

        let (_smi, fresh) = AtomTracker::write_parse(&product).unwrap();
        for (tag, z) in before {
            let i = AtomTracker::index_of(&fresh, tag).unwrap();
            assert_eq!(fresh.atom(atom_idx(i)).element.atomic_number(), z);
        }
    }
}
