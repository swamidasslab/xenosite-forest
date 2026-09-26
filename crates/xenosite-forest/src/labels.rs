//! Atom identity beside the chematic graph.
//!
//! Chematic on crates.io has no userdata; `atom_map` is a SMIRKS match key
//! (apply clears it; canonical write emits `:n`). This repo vendors chematic
//! with `Atom.tag` and SMILES visit-order helpers — see
//! [`crate::chematic_vendor`] and `patches/README.md`.
//!
//! [`ForestMol`] still keeps a tag sidecar for derisk continuity. Tags are
//! remapped by a correspondence from each rewrite when one is supplied.
//! Tests can also recover maps with public `set_isotope` plus write/parse.

#[cfg(test)]
use std::collections::HashSet;

#[cfg(test)]
use chematic::core::{BondIdx, Molecule};

#[cfg(test)]
use crate::mol::{atom_idx, atom_usize};

/// Stable atom id in one `ForestMol` copy tree. Not a SMIRKS map number.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Tag(pub u32);

/// Stamp a fresh tag on every atom, in current index order.
pub fn stamp(n: usize) -> (Vec<Option<Tag>>, u32) {
    let labels: Vec<Option<Tag>> = (0..n).map(|i| Some(Tag(i as u32))).collect();
    (labels, n as u32)
}

/// Carry tags through an apply: `src_to_new[src] = Some(dst)` or `None` if
/// that atom left this piece. Unmapped product atoms get new tags.
pub fn remap_apply(
    parent: &[Option<Tag>],
    src_to_new: &[Option<usize>],
    product_n: usize,
    mut next: u32,
) -> (Vec<Option<Tag>>, u32) {
    let mut labels = vec![None; product_n];
    for (src, dst) in src_to_new.iter().copied().enumerate() {
        let Some(dst) = dst else {
            continue;
        };
        if dst < product_n {
            labels[dst] = parent.get(src).copied().flatten();
        }
    }
    for slot in &mut labels {
        if slot.is_none() {
            *slot = Some(Tag(next));
            next += 1;
        }
    }
    (labels, next)
}

/// `old_at_new[new_idx] = old_idx` after a permutation that keeps atom count.
pub fn remap_permute(parent: &[Option<Tag>], old_at_new: &[usize]) -> Vec<Option<Tag>> {
    old_at_new
        .iter()
        .map(|&old| parent.get(old).copied().flatten())
        .collect()
}

/// Index-stable edit (clone or append). New atoms at the end get new tags.
pub fn remap_index_stable(
    parent: &[Option<Tag>],
    product_n: usize,
    mut next: u32,
) -> (Vec<Option<Tag>>, u32) {
    let mut labels: Vec<Option<Tag>> = parent.iter().copied().take(product_n).collect();
    labels.resize(product_n, None);
    for slot in &mut labels {
        if slot.is_none() {
            *slot = Some(Tag(next));
            next += 1;
        }
    }
    (labels, next)
}

/// Non-canonical SMILES visit order: a **test clone** of chematic `write` DFS.
///
/// Not a public Forest seam. If this disagrees with isotope write/parse, the
/// writer changed; do not "fix" it by reaching into chematic-smiles.
#[cfg(test)]
fn write_visit_order(mol: &Molecule) -> Vec<usize> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }
    let ring_bonds = ring_closure_bonds(mol);
    let mut written = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for start in 0..n {
        if written[start] {
            continue;
        }
        write_chain(mol, &ring_bonds, &mut written, &mut order, start, None);
    }
    order
}

#[cfg(test)]
fn ring_closure_bonds(mol: &Molecule) -> HashSet<BondIdx> {
    let n = mol.atom_count();
    let mut ring_bonds = HashSet::new();
    let mut visited = vec![false; n];
    let mut in_stack = vec![false; n];
    for start in 0..n {
        if visited[start] {
            continue;
        }
        mark_rings(
            mol,
            &mut ring_bonds,
            &mut visited,
            &mut in_stack,
            start,
            None,
        );
    }
    ring_bonds
}

#[cfg(test)]
fn mark_rings(
    mol: &Molecule,
    ring_bonds: &mut HashSet<BondIdx>,
    visited: &mut [bool],
    in_stack: &mut [bool],
    atom: usize,
    from_bond: Option<BondIdx>,
) {
    visited[atom] = true;
    in_stack[atom] = true;
    for (neighbor, bidx) in mol.neighbors(atom_idx(atom)) {
        if Some(bidx) == from_bond || ring_bonds.contains(&bidx) {
            continue;
        }
        let other = atom_usize(neighbor);
        if !visited[other] {
            mark_rings(mol, ring_bonds, visited, in_stack, other, Some(bidx));
        } else if in_stack[other] {
            ring_bonds.insert(bidx);
        }
    }
    in_stack[atom] = false;
}

#[cfg(test)]
fn write_chain(
    mol: &Molecule,
    ring_bonds: &HashSet<BondIdx>,
    written: &mut [bool],
    order: &mut Vec<usize>,
    atom: usize,
    from_atom: Option<usize>,
) {
    written[atom] = true;
    order.push(atom);
    let children: Vec<usize> = mol
        .neighbors(atom_idx(atom))
        .filter_map(|(nb, bidx)| {
            let other = atom_usize(nb);
            (Some(other) != from_atom && !written[other] && !ring_bonds.contains(&bidx))
                .then_some(other)
        })
        .collect();
    for child in children {
        write_chain(mol, ring_bonds, written, order, child, Some(atom));
    }
}

/// Recover `src_to_new` when each parent atom has a unique isotope.
///
/// Test probe only. Production tags are the sidecar; chematic apply does
/// not return `src_to_new`, but it does clone `Atom` (so isotope survives).
#[cfg(test)]
fn src_to_new_from_isotopes(parent: &Molecule, product: &Molecule) -> Vec<Option<usize>> {
    let mut src_to_new = vec![None; parent.atom_count()];
    for (src, atom) in parent.atoms() {
        let Some(iso) = atom.isotope else {
            continue;
        };
        let dst = product
            .atoms()
            .find_map(|(idx, other)| (other.isotope == Some(iso)).then_some(atom_usize(idx)));
        src_to_new[atom_usize(src)] = dst;
    }
    src_to_new
}

/// `old_at_new[new_idx] = old_idx` using unique isotopes on both mols.
#[cfg(test)]
fn permute_from_isotopes(old: &Molecule, new: &Molecule) -> Vec<usize> {
    let n = new.atom_count();
    (0..n)
        .map(|new_i| {
            let iso = new.atom(atom_idx(new_i)).isotope;
            old.atoms()
                .find_map(|(idx, atom)| (atom.isotope == iso).then_some(atom_usize(idx)))
                .expect("isotope identity lost across rewrite")
        })
        .collect()
}

#[cfg(test)]
fn stamp_unique_isotopes(mol: &mut Molecule) {
    for i in 0..mol.atom_count() {
        mol.set_isotope(atom_idx(i), Some(200 + i as u16));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forest_mol::ForestMol;
    use crate::mol::{canon_of, parse_mol};
    use crate::smarts::smarts_matches;
    use crate::smirks::apply_smirks_at;
    use chematic::core::{Atom, BondOrder, Element, MoleculeBuilder};
    use chematic::smiles::{canonical_atom_order, canonical_smiles, write};

    fn identity(n: usize) -> Vec<usize> {
        (0..n).collect()
    }

    /// Ethanol with stored indexes `O=0`, methyl `=1`, center `=2` (bonds 0–2–1).
    /// Plain `write` DFS from 0 visits `[0, 2, 1]`, not `0,1,2`.
    fn ethanol_scrambled() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let oxygen = builder.add_atom(Atom::organic(Element::O));
        let methyl = builder.add_atom(Atom::organic(Element::C));
        let center = builder.add_atom(Atom::organic(Element::C));
        builder.add_bond(oxygen, center, BondOrder::Single).unwrap();
        builder.add_bond(center, methyl, BondOrder::Single).unwrap();
        builder.build()
    }

    fn tag_of_z(mol: &ForestMol, z: u8) -> Tag {
        (0..mol.mol().atom_count())
            .find(|&i| mol.mol().atom(atom_idx(i)).element.atomic_number() == z)
            .and_then(|i| mol.tag_of(i))
            .expect("element present")
    }

    #[test]
    fn parse_stamps_a_tag_per_atom() {
        let mol = ForestMol::parse("CCO").unwrap();
        assert_eq!(mol.tag_of(0), Some(Tag(0)));
        assert_eq!(mol.tag_of(1), Some(Tag(1)));
        assert_eq!(mol.tag_of(2), Some(Tag(2)));
        assert_eq!(mol.index_of(Tag(1)), Some(1));
        assert_eq!(mol.tag_of(3), None);
    }

    #[test]
    fn plain_write_visit_is_not_index_order() {
        let mol = ethanol_scrambled();
        let visit = write_visit_order(&mol);
        assert_eq!(visit, vec![0, 2, 1]);
        assert_ne!(visit, identity(mol.atom_count()));
    }

    #[test]
    fn canonical_atom_order_is_not_plain_write_visit() {
        let mol = parse_mol("CC(C)C").unwrap();
        let visit = write_visit_order(&mol);
        let ranks = canonical_atom_order(&mol);
        assert_ne!(
            ranks, visit,
            "canonical_atom_order is Morgan rank sort, not SMILES DFS"
        );
        assert_eq!(
            visit[0], 0,
            "plain write starts at the first unwritten index"
        );
        assert_ne!(ranks[0], visit[0]);
    }

    #[test]
    fn canonical_smiles_dfs_is_not_canonical_atom_order() {
        let mut mol = parse_mol("c1ccccc1O").unwrap();
        stamp_unique_isotopes(&mut mol);
        let smiles = canonical_smiles(&mol);
        let parsed = parse_mol(&smiles).unwrap();
        let dfs = permute_from_isotopes(&mol, &parsed);
        let ranks = canonical_atom_order(&mol);
        assert_ne!(
            dfs, ranks,
            "canonical SMILES visit is a DFS, not the rank-sorted atom list"
        );
        assert_eq!(parsed.atom_count(), mol.atom_count());
    }

    #[test]
    fn plain_write_then_parse_remaps_tags_by_visit_order() {
        let parent = ForestMol::new(ethanol_scrambled());
        let mut probe = parent.mol().clone();
        stamp_unique_isotopes(&mut probe);
        let smiles = write(&probe);
        let parsed = parse_mol(&smiles).unwrap();
        let visit = permute_from_isotopes(&probe, &parsed);
        assert_eq!(visit, write_visit_order(parent.mol()));
        assert_eq!(visit, vec![0, 2, 1]);

        let child = parent.after_permute(&parsed, &visit);
        for (new_i, &old_i) in visit.iter().enumerate() {
            assert_eq!(child.tag_of(new_i), parent.tag_of(old_i));
            assert_eq!(
                child.mol().atom(atom_idx(new_i)).element.atomic_number(),
                parent.mol().atom(atom_idx(old_i)).element.atomic_number()
            );
        }
        let oxygen = tag_of_z(&parent, 8);
        let new_o = child.index_of(oxygen).unwrap();
        assert_eq!(child.mol().atom(atom_idx(new_o)).element.atomic_number(), 8);
        assert_eq!(new_o, 0, "plain write started at oxygen (index 0)");
    }

    #[test]
    fn canonical_smiles_then_parse_remaps_tags_by_dfs_not_rank_order() {
        let parent = ForestMol::parse("c1ccccc1O").unwrap();
        let mut probe = parent.mol().clone();
        stamp_unique_isotopes(&mut probe);
        let smiles = canonical_smiles(&probe);
        let parsed = parse_mol(&smiles).unwrap();
        let dfs = permute_from_isotopes(&probe, &parsed);
        assert_ne!(dfs, canonical_atom_order(parent.mol()));

        let child = parent.after_permute(&parsed, &dfs);
        let oxygen = tag_of_z(&parent, 8);
        let new_o = child.index_of(oxygen).unwrap();
        assert_eq!(child.mol().atom(atom_idx(new_o)).element.atomic_number(), 8);
        assert_eq!(child.tag_of(new_o), Some(oxygen));
    }

    #[test]
    fn apply_clears_atom_maps_so_they_cannot_be_tags() {
        let mol = parse_mol("[CH3:1][CH2:2][OH:3]").unwrap();
        assert!(mol.atoms().any(|(_, atom)| atom.atom_map.is_some()));
        let hits = smarts_matches(&mol, "[#6h3:1]").unwrap();
        let products = apply_smirks_at("[C:1]>>[C:1]O", &mol, &hits[0]).unwrap();
        assert!(!products.is_empty());
        for product in &products {
            assert!(
                product.atoms().all(|(_, atom)| atom.atom_map.is_none()),
                "chematic build_product clears atom_map"
            );
        }
    }

    #[test]
    fn apply_then_index_of_tag_finds_the_same_atom() {
        let parent = ForestMol::parse("CCO").unwrap();
        let hits = parent.smarts_matches("[#6h3:1]").unwrap();
        let mapped = &hits[0];

        let mut probe = parent.mol().clone();
        stamp_unique_isotopes(&mut probe);
        let probed = apply_smirks_at("[C:1]>>[C:1]O", &probe, mapped).unwrap();
        assert_eq!(probed.len(), 1);
        let src_to_new = src_to_new_from_isotopes(&probe, &probed[0]);
        let carried = src_to_new.iter().filter(|d| d.is_some()).count();
        assert_eq!(carried, parent.mol().atom_count());

        let child = parent.from_apply(probed[0].clone(), &src_to_new);
        for (src, dst) in src_to_new.iter().copied().enumerate() {
            let tag = parent.tag_of(src).unwrap();
            let dst = dst.unwrap();
            assert_eq!(child.index_of(tag), Some(dst));
            assert_eq!(
                child.mol().atom(atom_idx(dst)).isotope,
                probe.atom(atom_idx(src)).isotope
            );
        }
        let born = (0..child.mol().atom_count())
            .filter(|&i| parent.index_of(child.tag_of(i).unwrap()).is_none())
            .collect::<Vec<_>>();
        assert_eq!(born.len(), 1);
        assert_eq!(
            child.mol().atom(atom_idx(born[0])).element.atomic_number(),
            8
        );

        let clean = apply_smirks_at("[C:1]>>[C:1]O", parent.mol(), mapped).unwrap();
        assert_eq!(
            canon_of(&chematic::smiles::canonical_smiles(&clean[0])).unwrap(),
            canon_of("OCCO").unwrap()
        );
    }

    #[test]
    fn index_stable_product_keeps_tags_on_old_atoms() {
        let parent = ForestMol::parse("CCO").unwrap();
        let tag = parent.tag_of(0).unwrap();
        let child = parent.edit_copy();
        assert_eq!(child.tag_of(0), Some(tag));
        assert_eq!(child.index_of(tag), Some(0));
        // `edit_copy` / `product` uses a fresh tag_gen (born atoms continue
        // the counter value, not the Rc). `copy_mol` shares the gen.
        assert!(!child.shares_tag_gen(&parent));
        assert!(parent.copy_mol().shares_tag_gen(&parent));
    }
}
