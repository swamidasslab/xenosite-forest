//! Product gate: chematic valence plus forest's two-double nitrogen drop.
//!
//! Also: generic closed-shell H fill from bond orders + charge (emit path),
//! not rule-named chemistry.

use std::collections::BTreeSet;

use chematic::core::BondOrder;
use chematic::perception::validate_valence;

use crate::kekule::{bond_order_sums, with_atom_explicit_h};
use crate::mol::{Molecule, atom_idx, atom_usize};

/// True when a nitrogen has two double bonds (`C=[N+]=C` is not an iminium).
pub fn nitrogen_two_doubles(mol: &Molecule) -> bool {
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 7 {
            continue;
        }
        let mut doubles = 0;
        for (_nbr, bond_idx) in mol.neighbors(idx) {
            if mol.bond(bond_idx).order == BondOrder::Double {
                doubles += 1;
                if doubles >= 2 {
                    return true;
                }
            }
        }
    }
    false
}

/// Forest `_sanitize_piece`: valence ok and not the two-double nitrogen.
pub fn accept_product(mol: &Molecule) -> bool {
    if !validate_valence(mol).is_empty() {
        return false;
    }
    !nitrogen_two_doubles(mol)
}

/// Set explicit H so atom valence is complete (closed shell).
///
/// Generic emit-path bookkeeping after graph edits (forced doubles, cleavage,
/// saturate). Organic targets: C 4, N 3, O 2, adjusted by formal charge.
/// Not rule-named chemistry — Kekulé matching stays charge/H-agnostic.
pub fn fill_closed_shell_h(mol: &mut Molecule, atom: usize) {
    let idx = atom_idx(atom);
    let a = mol.atom(idx);
    let z = a.element.atomic_number();
    let charge = a.charge as i16;
    let mut bond_sum = 0.0_f32;
    for (_nbr, bidx) in mol.neighbors(idx) {
        bond_sum += match mol.bond(bidx).order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::Quadruple => 4.0,
            _ => 1.0,
        };
    }
    let target = match z {
        6 => 4 + charge,
        7 => 3 + charge,
        8 => 2 + charge,
        _ => return,
    };
    let need = target - bond_sum.round() as i16;
    if need >= 0 {
        *mol = with_atom_explicit_h(mol, idx, need as u8);
    }
}

/// Atoms whose neighbor set changed between `parent` and `edited` (same index
/// layout), plus atoms born on `edited`. Perception — not edit-token names.
pub fn skeleton_changed_between(parent: &Molecule, edited: &Molecule) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    let n = parent.atom_count().min(edited.atom_count());
    for i in 0..n {
        let pn: BTreeSet<_> = parent
            .neighbors(atom_idx(i))
            .map(|(nbr, _)| atom_usize(nbr))
            .collect();
        let en: BTreeSet<_> = edited
            .neighbors(atom_idx(i))
            .map(|(nbr, _)| atom_usize(nbr))
            .collect();
        if pn != en {
            out.insert(i);
        }
    }
    for i in n..edited.atom_count() {
        out.insert(i);
    }
    out
}

/// Atoms that need emit-path closed-shell settle and charge-baseline reset:
/// neighbor-set change **or** incident bond-order sum change (forced
/// single→double, new leaf, cleavage). Perception — not edit tokens.
pub fn edited_valence_atoms(parent: &Molecule, edited: &Molecule) -> BTreeSet<usize> {
    let mut out = skeleton_changed_between(parent, edited);
    let before = bond_order_sums(parent);
    let after = bond_order_sums(edited);
    for (&a, &v) in &after {
        if before.get(&a).copied() != Some(v) {
            out.insert(a);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn carbamazepine_style_two_double_nitrogen_is_refused() {
        let mol = parse_mol("C=[N+]=C").unwrap();
        assert!(nitrogen_two_doubles(&mol));
        assert!(!accept_product(&mol));
    }

    #[test]
    fn iminium_is_kept() {
        let mol = parse_mol("C[N+](C)=C").unwrap();
        assert!(!nitrogen_two_doubles(&mol));
        assert!(accept_product(&mol));
    }

    #[test]
    fn ethanol_is_kept() {
        let mol = parse_mol("CCO").unwrap();
        assert!(accept_product(&mol));
    }

    #[test]
    fn fill_closed_shell_h_makes_methane_from_isolated_carbon() {
        let mut mol = parse_mol("C").unwrap();
        // Force zero explicit H then refill.
        fill_closed_shell_h(&mut mol, 0);
        assert!(accept_product(&mol));
    }

    #[test]
    fn edited_valence_atoms_sees_bond_order_change() {
        let parent = parse_mol("CCO").unwrap();
        let mut edited = parent.clone();
        // Raise C–O to double (phenol/enol-style end edit).
        let (bi, _) = edited
            .bond_between(atom_idx(1), atom_idx(2))
            .expect("C-O");
        edited.set_bond_order(bi, BondOrder::Double);
        let settle = edited_valence_atoms(&parent, &edited);
        assert!(settle.contains(&1) && settle.contains(&2));
        fill_closed_shell_h(&mut edited, 2);
        assert_eq!(edited.atom(atom_idx(2)).charge, 0);
    }
}
