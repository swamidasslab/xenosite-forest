//! Smoke that the vendored chematic patch works: `Atom.tag` and visit order.

#[cfg(test)]
mod tests {
    use chematic::core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic::rxn::{apply_reaction_match, find_reaction_matches};
    use chematic::smiles::{canonical_smiles_with_order, parse, write, write_with_order};

    #[test]
    fn write_with_order_matches_write_and_can_differ_from_index_order() {
        let mol = parse("CCO").unwrap();
        let (smi, order) = write_with_order(&mol);
        assert_eq!(smi, write(&mol));
        assert_eq!(order.len(), mol.atom_count());

        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let me = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::C));
        b.add_bond(c, o, BondOrder::Single).unwrap();
        b.add_bond(c, me, BondOrder::Single).unwrap();
        let scrambled = b.build();
        let (_s, ord) = write_with_order(&scrambled);
        let idxs: Vec<_> = ord.iter().map(|a| a.0).collect();
        assert_ne!(idxs, vec![0, 1, 2]);
    }

    #[test]
    fn canonical_smiles_with_order_returns_visit_permutation() {
        let mol = parse("CCO").unwrap();
        let (csmi, order) = canonical_smiles_with_order(&mol);
        assert!(!csmi.is_empty());
        assert_eq!(order.len(), mol.atom_count());
        let mut seen = order.clone();
        seen.sort_by_key(|a| a.0);
        assert_eq!(
            seen,
            (0..mol.atom_count() as u32)
                .map(AtomIdx)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn atom_tag_survives_apply_atom_map_does_not() {
        let mut mol = parse("CCO").unwrap();
        for i in 0..mol.atom_count() {
            mol.set_tag(AtomIdx(i as u32), Some(200 + i as u32));
            // Leave any existing maps alone; apply must clear them.
        }
        let matches = find_reaction_matches("[C:1]>>[C:1]O", &[&mol]).expect("match");
        assert!(!matches.is_empty());
        let products =
            apply_reaction_match("[C:1]>>[C:1]O", &[&mol], &matches[0], true).expect("apply");
        let product = &products.expect("valence")[0];
        let tags: Vec<_> = product.atoms().map(|(_, a)| a.tag).collect();
        assert!(
            tags.contains(&Some(200)) || tags.contains(&Some(201)) || tags.contains(&Some(202))
        );
        assert!(product.atoms().all(|(_, a)| a.atom_map.is_none()));
        // Born oxygen has no tag.
        assert!(
            product
                .atoms()
                .any(|(_, a)| a.tag.is_none() && a.element.atomic_number() == 8)
        );
    }
}
