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

/// Forest `_sanitize_piece`: valence ok, not the two-double nitrogen, and
/// closed-shell (chematic does not model radicals — refuse underfilled C/N/O).
pub fn accept_product(mol: &Molecule) -> bool {
    if !validate_valence(mol).is_empty() {
        return false;
    }
    if nitrogen_two_doubles(mol) {
        return false;
    }
    if oxygen_oxonium(mol) {
        return false;
    }
    closed_shell(mol)
}

/// Protonated carbonyl / phenol oxonium (`=[OH+]`): closed-shell prefer
/// neutral (HEURISTICS C10 — not a sanitize rescue, a refuse).
fn oxygen_oxonium(mol: &Molecule) -> bool {
    for (_idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 8 && atom.charge > 0 {
            return true;
        }
    }
    false
}

/// Cumulative `C=C=O` (ketene carbon: C with double bonds to C and O).
pub fn has_ketene(mol: &Molecule) -> bool {
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        let mut to_c = false;
        let mut to_o = false;
        for (nbr, bidx) in mol.neighbors(idx) {
            if mol.bond(bidx).order != BondOrder::Double {
                continue;
            }
            match mol.atom(nbr).element.atomic_number() {
                6 => to_c = true,
                8 => to_o = true,
                _ => {}
            }
        }
        if to_c && to_o {
            return true;
        }
    }
    false
}

/// True when any atom is aromatic.
pub fn has_aromatic_atom(mol: &Molecule) -> bool {
    mol.atoms().any(|(_, a)| a.aromatic)
}

/// O-leave cleavage must not fully dearomatize into a ketene (benzoic acid →
/// `O=C=C1CCCCC1`) or drop a full aromatic sextet (nitroarene β-elim junk).
/// Ring-open dealk ketenes that keep an aromatic piece pass. Leave fragments
/// (water) are not judged — only heavy products.
pub fn accept_o_leave_product(parent: &Molecule, product: &Molecule) -> bool {
    if !accept_product(product) {
        return false;
    }
    let heavy = product
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count();
    if heavy <= 1 {
        return true;
    }
    let parent_arom = parent.atoms().filter(|(_, a)| a.aromatic).count();
    let product_arom = product.atoms().filter(|(_, a)| a.aromatic).count();
    if parent_arom >= 6 && product_arom + 6 <= parent_arom {
        return false;
    }
    if has_aromatic_atom(parent) && !has_aromatic_atom(product) && has_ketene(product) {
        return false;
    }
    true
}

/// Organic C/N/O atoms have enough bonds+H for a closed shell (no radicals).
fn closed_shell(mol: &Molecule) -> bool {
    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        let target = match z {
            6 => 4 + atom.charge as i16,
            7 => 3 + atom.charge as i16,
            8 => 2 + atom.charge as i16,
            _ => continue,
        };
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
        let h = atom.hydrogen_count.unwrap_or_else(|| {
            // Implicit H: chematic inference when unset.
            mol.implicit_hydrogen_count(idx)
        }) as i16;
        let used = bond_sum.round() as i16 + h;
        if used < target {
            return false;
        }
    }
    true
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
    fn radical_carbon_is_refused() {
        // Formyl radical — chematic valence may pass; closed-shell gate refuses.
        let mol = parse_mol("[C]=O").unwrap();
        assert!(!accept_product(&mol), "open-shell [C]=O must be refused");
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
