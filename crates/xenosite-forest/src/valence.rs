//! Product gate: chematic valence plus forest's two-double nitrogen drop.

use chematic::core::BondOrder;
use chematic::perception::validate_valence;

use crate::mol::Molecule;

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
}
