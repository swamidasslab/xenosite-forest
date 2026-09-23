//! Atom tracker POC on chematic `Atom.tag` + SMILES visit order.
//!
//! Stamps a unique integer on each atom. Apply copies it (patch); born atoms
//! are `None` until [`AtomTracker::adopt_born`] mints ids. Write/parse uses
//! visit-order helpers so tags ride the permutation onto the new mol without
//! isotopes or a sidecar.
//!
//! This is the production-shaped seam once Forest stops carrying
//! [`crate::labels`] beside the graph. `ForestMol`'s sidecar stays for now.

use chematic::core::{Molecule, MoleculeBuilder};
use chematic::smiles::{canonical_smiles_with_order, parse, write_with_order};

use crate::labels::Tag;
use crate::mol::{ForestError, atom_idx};

/// Mints and follows chematic atom tags through rewrite.
#[derive(Debug, Clone)]
pub struct AtomTracker {
    next: u32,
}

impl AtomTracker {
    /// Stamp every atom with a fresh tag. Existing tags are overwritten.
    ///
    /// The counter restarts at `atom_count()` — prior high ids from a shared
    /// tracker are not preserved. Callers that chain applies should keep the
    /// same [`AtomTracker`] and use [`Self::adopt_born`], not re-stamp.
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

    /// Tag at `idx`, or `None` if out of range or untagged.
    pub fn tag_of(mol: &Molecule, idx: usize) -> Option<Tag> {
        mol.atom_opt(atom_idx(idx))?.tag.map(Tag)
    }

    /// First index holding `tag`, if any. Duplicate tags are ambiguous.
    pub fn index_of(mol: &Molecule, tag: Tag) -> Option<usize> {
        (0..mol.atom_count()).find(|&i| mol.atom(atom_idx(i)).tag == Some(tag.0))
    }

    /// Snapshot of every atom's tag (parallel to index order).
    pub fn all_tags(mol: &Molecule) -> Vec<Option<Tag>> {
        (0..mol.atom_count())
            .map(|i| Self::tag_of(mol, i))
            .collect()
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
    pub fn write_parse(mol: &Molecule) -> Result<(String, Molecule), ForestError> {
        Self::remap_write_parse(mol, write_with_order)
    }

    /// Canonical write, then parse, with tags remapped by canonical DFS visit.
    pub fn canonical_write_parse(mol: &Molecule) -> Result<(String, Molecule), ForestError> {
        Self::remap_write_parse(mol, canonical_smiles_with_order)
    }

    fn remap_write_parse(
        mol: &Molecule,
        write: fn(&Molecule) -> (String, Vec<chematic::core::AtomIdx>),
    ) -> Result<(String, Molecule), ForestError> {
        if mol.atom_count() == 0 {
            // `write`/`canonical_smiles` return ""; `parse("")` is EmptyInput.
            return Ok((String::new(), MoleculeBuilder::new().build()));
        }
        let (smi, order) = write(mol);
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
    use std::collections::{BTreeMap, HashSet};

    use chematic::core::{Atom, BondOrder, Element, MoleculeBuilder};
    use chematic::rxn::{apply_reaction_match, find_reaction_matches};
    use chematic::smiles::{canonical_smiles_with_order, parse, write, write_with_order};

    use super::*;
    use crate::mol::atom_usize;
    use crate::smirks::apply_smirks_at;

    fn hydroxylate_first(mol: &Molecule) -> Molecule {
        let matches = find_reaction_matches("[C:1]>>[C:1]O", &[mol]).unwrap();
        assert!(!matches.is_empty(), "expected a carbon match");
        apply_reaction_match("[C:1]>>[C:1]O", &[mol], &matches[0], true)
            .unwrap()
            .expect("valence")
            .into_iter()
            .next()
            .unwrap()
    }

    fn assert_unique_tags(mol: &Molecule) {
        let mut seen = HashSet::new();
        for tag in AtomTracker::all_tags(mol).into_iter().flatten() {
            assert!(seen.insert(tag), "duplicate tag {tag:?}");
        }
    }

    fn assert_identity_roundtrip(mol: &Molecule) {
        assert_unique_tags(mol);
        let before: BTreeMap<Tag, u8> = (0..mol.atom_count())
            .filter_map(|i| {
                AtomTracker::tag_of(mol, i)
                    .map(|t| (t, mol.atom(atom_idx(i)).element.atomic_number()))
            })
            .collect();
        let (_smi, fresh) = AtomTracker::write_parse(mol).unwrap();
        assert_eq!(fresh.atom_count(), mol.atom_count());
        for (tag, z) in &before {
            let i = AtomTracker::index_of(&fresh, *tag).expect("tag survived write/parse");
            assert_eq!(fresh.atom(atom_idx(i)).element.atomic_number(), *z);
        }
        let (_csmi, cfresh) = AtomTracker::canonical_write_parse(mol).unwrap();
        for (tag, z) in &before {
            let i =
                AtomTracker::index_of(&cfresh, *tag).expect("tag survived canonical write/parse");
            assert_eq!(cfresh.atom(atom_idx(i)).element.atomic_number(), *z);
        }
    }

    // --- API surface / boundaries ------------------------------------------------

    #[test]
    fn stamp_empty_molecule() {
        // parse("") is EmptyInput; builder yields a true empty graph.
        let mut mol = MoleculeBuilder::new().build();
        let tracker = AtomTracker::stamp(&mut mol);
        assert_eq!(mol.atom_count(), 0);
        assert_eq!(tracker.next_tag(), 0);
        assert!(AtomTracker::all_tags(&mol).is_empty());
        let (smi, fresh) = AtomTracker::write_parse(&mol).unwrap();
        assert_eq!(smi, "");
        assert_eq!(fresh.atom_count(), 0);
        let (csmi, cfresh) = AtomTracker::canonical_write_parse(&mol).unwrap();
        assert_eq!(csmi, "");
        assert_eq!(cfresh.atom_count(), 0);
    }

    #[test]
    fn stamp_single_atom_and_lookups() {
        let mut mol = parse("C").unwrap();
        let tracker = AtomTracker::stamp(&mut mol);
        assert_eq!(tracker.next_tag(), 1);
        assert_eq!(AtomTracker::tag_of(&mol, 0), Some(Tag(0)));
        assert_eq!(AtomTracker::tag_of(&mol, 1), None); // OOB
        assert_eq!(AtomTracker::index_of(&mol, Tag(0)), Some(0));
        assert_eq!(AtomTracker::index_of(&mol, Tag(99)), None);
        assert_eq!(AtomTracker::all_tags(&mol), vec![Some(Tag(0))]);
    }

    #[test]
    fn stamp_overwrites_and_restarts_counter() {
        let mut mol = parse("CCO").unwrap();
        let mut tracker = AtomTracker::stamp(&mut mol);
        mol.set_tag(atom_idx(0), Some(900));
        tracker.adopt_born(&mut mol); // no-op; all tagged
        assert_eq!(tracker.next_tag(), 3);
        // Re-stamp is destructive: ids restart at 0..n.
        let tracker2 = AtomTracker::stamp(&mut mol);
        assert_eq!(
            AtomTracker::all_tags(&mol),
            vec![Some(Tag(0)), Some(Tag(1)), Some(Tag(2))]
        );
        assert_eq!(tracker2.next_tag(), 3);
        assert_ne!(AtomTracker::tag_of(&mol, 0), Some(Tag(900)));
    }

    #[test]
    fn adopt_born_is_idempotent_and_skips_tagged() {
        let mut mol = parse("CC").unwrap();
        let mut tracker = AtomTracker::stamp(&mut mol);
        let before = tracker.next_tag();
        tracker.adopt_born(&mut mol);
        assert_eq!(tracker.next_tag(), before);
        mol.set_tag(atom_idx(1), None);
        tracker.adopt_born(&mut mol);
        assert_eq!(AtomTracker::tag_of(&mol, 1), Some(Tag(before)));
        assert_eq!(tracker.next_tag(), before + 1);
        tracker.adopt_born(&mut mol);
        assert_eq!(tracker.next_tag(), before + 1);
    }

    #[test]
    fn duplicate_tags_make_index_of_ambiguous() {
        let mut mol = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        mol.set_tag(atom_idx(2), Some(0)); // collide with atom 0
        assert_eq!(AtomTracker::index_of(&mol, Tag(0)), Some(0)); // first only
        let tagged: Vec<_> = AtomTracker::all_tags(&mol)
            .into_iter()
            .flatten()
            .filter(|t| *t == Tag(0))
            .collect();
        assert_eq!(tagged.len(), 2);
    }

    #[test]
    fn src_to_new_identity_and_missing_parent_tags() {
        let mut mol = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        let map = AtomTracker::src_to_new(&mol, &mol);
        assert_eq!(map, vec![Some(0), Some(1), Some(2)]);

        mol.set_tag(atom_idx(1), None);
        let map = AtomTracker::src_to_new(&mol, &mol);
        assert_eq!(map[1], None);
        assert_eq!(map[0], Some(0));
    }

    #[test]
    fn write_parse_preserves_none_tags() {
        let mut mol = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        mol.set_tag(atom_idx(1), None);
        let (_smi, fresh) = AtomTracker::write_parse(&mol).unwrap();
        let none_count = AtomTracker::all_tags(&fresh)
            .into_iter()
            .filter(|t| t.is_none())
            .count();
        assert_eq!(none_count, 1);
        assert!(AtomTracker::index_of(&fresh, Tag(0)).is_some());
        assert!(AtomTracker::index_of(&fresh, Tag(2)).is_some());
    }

    #[test]
    fn atom_map_cleared_on_apply_tags_kept() {
        let mut parent = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut parent);
        let mut b = MoleculeBuilder::new();
        for i in 0..parent.atom_count() {
            let mut atom = parent.atom(atom_idx(i)).clone();
            atom.atom_map = Some(10 + i as u16);
            b.add_atom(atom);
        }
        for (_, bond) in parent.bonds() {
            b.add_bond(bond.atom1, bond.atom2, bond.order).unwrap();
        }
        let parent = b.build();
        assert!(parent.atoms().any(|(_, a)| a.atom_map.is_some()));
        assert!(parent.atoms().all(|(_, a)| a.tag.is_some()));

        let product = hydroxylate_first(&parent);
        assert!(product.atoms().all(|(_, a)| a.atom_map.is_none()));
        let surviving = AtomTracker::src_to_new(&parent, &product)
            .into_iter()
            .flatten()
            .count();
        assert_eq!(surviving, parent.atom_count());
    }

    // --- Apply / fragments / cleavage -------------------------------------------

    #[test]
    fn follow_apply_keeps_carbon_tags_and_mints_born_oxygen() {
        let mut parent = parse("CCO").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let c0 = AtomTracker::tag_of(&parent, 0).unwrap();
        let c1 = AtomTracker::tag_of(&parent, 1).unwrap();
        let o = AtomTracker::tag_of(&parent, 2).unwrap();

        let mut product = hydroxylate_first(&parent);
        let map = AtomTracker::src_to_new(&parent, &product);
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
    fn cleavage_leaves_some_src_to_new_none_per_fragment() {
        let mut parent = parse("CN").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let hits = crate::smarts::smarts_matches(&parent, "[#6H3:1][#7:2]").unwrap();
        let products = apply_smirks_at("[C:1][N:2]>>[N:2].[C:1](=O)O", &parent, &hits[0]).unwrap();
        assert_eq!(products.len(), 2);

        let mut covered = HashSet::new();
        for mut frag in products {
            let map = AtomTracker::src_to_new(&parent, &frag);
            let survivors: Vec<_> = map
                .iter()
                .copied()
                .enumerate()
                .filter_map(|(s, d)| d.map(|d| (s, d)))
                .collect();
            assert!(!survivors.is_empty());
            assert!(survivors.len() < parent.atom_count());
            for (src, _) in &survivors {
                covered.insert(*src);
            }
            tracker.adopt_born(&mut frag);
            assert_unique_tags(&frag);
            assert!(tags_agree_elements(&parent, &frag));
        }
        assert_eq!(covered.len(), parent.atom_count());
    }

    #[test]
    fn fragments_preserve_tags_without_apply() {
        let mut mol = parse("CCO.N").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        let frags = mol.fragments();
        assert_eq!(frags.len(), 2);
        let mut tags = HashSet::new();
        for frag in &frags {
            for t in AtomTracker::all_tags(frag).into_iter().flatten() {
                assert!(tags.insert(t), "tag leaked across fragments: {t:?}");
            }
        }
        assert_eq!(tags.len(), mol.atom_count());
    }

    #[test]
    fn multi_step_apply_chain_keeps_shared_tracker() {
        let mut mol = parse("CCC").unwrap();
        let mut tracker = AtomTracker::stamp(&mut mol);
        for step in 0..3 {
            let parent_tags: HashSet<_> =
                AtomTracker::all_tags(&mol).into_iter().flatten().collect();
            let mut product = hydroxylate_first(&mol);
            let map = AtomTracker::src_to_new(&mol, &product);
            assert_eq!(map.iter().flatten().count(), mol.atom_count());
            tracker.adopt_born(&mut product);
            let product_tags: HashSet<_> = AtomTracker::all_tags(&product)
                .into_iter()
                .flatten()
                .collect();
            assert!(parent_tags.is_subset(&product_tags));
            assert_eq!(product_tags.len(), mol.atom_count() + 1);
            assert_eq!(tracker.next_tag() as usize, 3 + step + 1);
            assert_unique_tags(&product);
            mol = product;
        }
        assert_identity_roundtrip(&mol);
    }

    // --- Write order / stereo / aromatic / disconnected -------------------------

    #[test]
    fn write_parse_remaps_tags_by_visit_order() {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let me = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::C));
        b.add_bond(c, o, BondOrder::Single).unwrap();
        b.add_bond(c, me, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let _tracker = AtomTracker::stamp(&mut mol);

        let (smi, order) = write_with_order(&mol);
        assert_eq!(smi, write(&mol));
        let idxs: Vec<_> = order.iter().map(|a| a.0).collect();
        assert_ne!(idxs, vec![0, 1, 2]);

        let tag_at_old: Vec<_> = (0..mol.atom_count())
            .map(|i| AtomTracker::tag_of(&mol, i).unwrap())
            .collect();
        let (_smi, fresh) = AtomTracker::write_parse(&mol).unwrap();
        for (new_i, &old) in order.iter().enumerate() {
            assert_eq!(
                AtomTracker::tag_of(&fresh, new_i).unwrap(),
                tag_at_old[atom_usize(old)]
            );
        }
    }

    #[test]
    fn canonical_write_parse_differs_from_index_order_on_ethanol() {
        let mut mol = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        let (csmi, order) = canonical_smiles_with_order(&mol);
        let (_s, fresh) = AtomTracker::canonical_write_parse(&mol).unwrap();
        assert!(!csmi.is_empty());
        let idxs: Vec<_> = order.iter().map(|a| a.0).collect();
        if idxs != (0..mol.atom_count() as u32).collect::<Vec<_>>() {
            for (new_i, &old) in order.iter().enumerate() {
                assert_eq!(
                    AtomTracker::tag_of(&fresh, new_i),
                    AtomTracker::tag_of(&mol, atom_usize(old))
                );
            }
        }
        assert_identity_roundtrip(&mol);
    }

    #[test]
    fn aromatic_ring_and_stereo_roundtrip() {
        for smi in [
            "c1ccccc1",
            "c1ccc(O)cc1",
            "C[C@H](O)C",
            "C/C=C/C",
            "CC(=O)Oc1ccccc1C(=O)O", // aspirin
        ] {
            let mut mol = parse(smi).unwrap();
            let _ = AtomTracker::stamp(&mut mol);
            assert_identity_roundtrip(&mol);
        }
    }

    #[test]
    fn disconnected_components_roundtrip() {
        let mut mol = parse("CCO.N.O=C=O").unwrap();
        let _ = AtomTracker::stamp(&mut mol);
        assert_eq!(mol.fragments().len(), 3);
        assert_identity_roundtrip(&mol);
    }

    // --- Larger molecules / stress ----------------------------------------------

    #[test]
    fn large_druglike_and_polyaromatic_roundtrips() {
        let c40: String = std::iter::repeat_n('C', 40).collect();
        let cases = [
            ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
            ("ibuprofen", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
            ("butylbenzene", "c1ccc(CCCC)cc1"),
            ("triph_butyl", "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3"),
            ("macrocycle_nd", "C1CCCCCCNC2CCCC(CC2)NCCCC1"),
            ("c40", c40.as_str()),
        ];
        for (name, smi) in cases {
            let mut mol = parse(smi).unwrap_or_else(|e| panic!("{name}: parse {e}"));
            assert!(
                mol.atom_count() >= 8 || name == "c40",
                "{name} too small for stress"
            );
            let tracker = AtomTracker::stamp(&mut mol);
            assert_eq!(tracker.next_tag() as usize, mol.atom_count());
            assert_unique_tags(&mol);
            assert_identity_roundtrip(&mol);
        }
    }

    #[test]
    fn large_mol_apply_then_write_parse() {
        let mut parent = parse("c1ccc(CCCC)cc1").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let n = parent.atom_count();
        let mut product = hydroxylate_first(&parent);
        let map = AtomTracker::src_to_new(&parent, &product);
        assert_eq!(map.iter().flatten().count(), n);
        tracker.adopt_born(&mut product);
        assert_eq!(product.atom_count(), n + 1);
        assert_unique_tags(&product);
        assert!(tags_agree_elements(&parent, &product));
        assert_identity_roundtrip(&product);
    }

    #[test]
    fn scrambled_builder_chain_then_apply() {
        let mut b = MoleculeBuilder::new();
        let leaves: Vec<_> = (0..3).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let center = b.add_atom(Atom::new(Element::C));
        for leaf in leaves {
            b.add_bond(center, leaf, BondOrder::Single).unwrap();
        }
        let mut mol = b.build();
        let mut tracker = AtomTracker::stamp(&mut mol);
        let (_s, order) = write_with_order(&mol);
        let mut product = hydroxylate_first(&mol);
        tracker.adopt_born(&mut product);
        assert_identity_roundtrip(&product);
        assert!(!order.is_empty());
    }

    #[test]
    fn tags_agree_elements_false_when_mutated() {
        let mut parent = parse("CCO").unwrap();
        let _ = AtomTracker::stamp(&mut parent);
        let mut product = parent.clone();
        product.set_element(atom_idx(0), Element::N);
        assert!(!tags_agree_elements(&parent, &product));
    }

    #[test]
    fn apply_then_write_parse_preserves_identity() {
        let mut parent = parse("CC").unwrap();
        let mut tracker = AtomTracker::stamp(&mut parent);
        let mut product = hydroxylate_first(&parent);
        tracker.adopt_born(&mut product);
        assert_identity_roundtrip(&product);
    }
}
