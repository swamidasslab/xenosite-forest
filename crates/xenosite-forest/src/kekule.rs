//! Kekulé parents for ResonanceRule: match aromatic, react on a bond-order parent.

use std::collections::BTreeMap;

use chematic::core::{BondOrder, apply_kekule, kekulize};
use chematic::perception::find_sssr;
use chematic::smarts::{BondPrimitive, BondQuery, parse_smarts};

use crate::mol::{ForestError, Molecule, atom_idx, atom_usize};

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn bond_order_map(mol: &Molecule) -> BTreeMap<(usize, usize), BondOrder> {
    mol.bonds()
        .map(|(_, bond)| {
            (
                bond_key(atom_usize(bond.atom1), atom_usize(bond.atom2)),
                bond.order,
            )
        })
        .collect()
}

fn implied_order(query: &BondQuery) -> Option<f32> {
    match query {
        BondQuery::Primitive(BondPrimitive::Double) => Some(2.0),
        BondQuery::Primitive(BondPrimitive::Single) => Some(1.0),
        BondQuery::Primitive(BondPrimitive::Aromatic) => Some(1.5),
        BondQuery::Or(left, right) => match (implied_order(left), implied_order(right)) {
            (Some(2.0), Some(1.5)) | (Some(1.5), Some(2.0)) => Some(2.0),
            (Some(1.0), Some(1.5)) | (Some(1.5), Some(1.0)) => Some(1.0),
            (a, b) => a.or(b),
        },
        BondQuery::Any | BondQuery::Primitive(BondPrimitive::Any) => Some(1.0),
        BondQuery::And(left, right) => implied_order(left).or_else(|| implied_order(right)),
        BondQuery::Not(_) => None,
        BondQuery::Primitive(_) => None,
    }
}

/// Bond order the reactant SMARTS asks for between maps 1 and 2.
pub fn smirks_mapped_bond_order(smirks: &str) -> Option<f32> {
    let reactant = smirks.split(">>").next()?;
    let query = parse_smarts(reactant).ok()?;
    let mut idxs = BTreeMap::new();
    for (i, atom) in query.atoms.iter().enumerate() {
        if let Some(mapno) = atom.atom_map {
            idxs.insert(mapno, i);
        }
    }
    let left = *idxs.get(&1)?;
    let right = *idxs.get(&2)?;
    let bond = query.bonds.iter().find(|bond| {
        (bond.atom1 == left && bond.atom2 == right) || (bond.atom1 == right && bond.atom2 == left)
    })?;
    implied_order(&bond.query)
}

fn clear_aromatic(mol: &Molecule) -> Molecule {
    let mut out = mol.clone();
    for (idx, atom) in mol.atoms() {
        if atom.aromatic {
            out = out.with_atom_aromatic(idx, false);
        }
    }
    out
}

fn apply_orders(mol: &Molecule, orders: &BTreeMap<(usize, usize), BondOrder>) -> Molecule {
    let mut out = mol.clone();
    for (bond_idx, bond) in mol.bonds() {
        let key = bond_key(atom_usize(bond.atom1), atom_usize(bond.atom2));
        if let Some(order) = orders.get(&key) {
            if bond.order != *order {
                out.set_bond_order(bond_idx, *order);
            }
        }
    }
    clear_aromatic(&out)
}

fn swap_even_ring(
    orders: &BTreeMap<(usize, usize), BondOrder>,
    ring: &[usize],
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    if ring.len() < 4 || ring.len() % 2 != 0 {
        return None;
    }
    let mut next = orders.clone();
    let n = ring.len();
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        let key = bond_key(a, b);
        let order = next.get(&key)?;
        let flipped = match order {
            BondOrder::Single => BondOrder::Double,
            BondOrder::Double => BondOrder::Single,
            _ => return None,
        };
        next.insert(key, flipped);
    }
    Some(next)
}

/// Kekulé forms with stable atom indexes (default matching plus even-ring swaps).
pub fn kekule_forms(mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
    if !mol.atoms().any(|(_, atom)| atom.aromatic) {
        return Ok(vec![mol.clone()]);
    }
    let assignment = kekulize(mol).map_err(|err| ForestError::Kekule(err.to_string()))?;
    // apply_kekule keeps aromatic atom flags; ResonanceRule parents are
    // integer-order graphs, so aliphatic SMIRKS ([S:1], [C:1]=[C:2]) can match.
    let primary = clear_aromatic(&apply_kekule(mol, &assignment));
    let mut forms = vec![primary.clone()];
    let base_orders = bond_order_map(&primary);
    let rings = find_sssr(mol);
    for ring in rings.rings() {
        let atoms: Vec<usize> = ring.iter().copied().map(atom_usize).collect();
        if let Some(swapped) = swap_even_ring(&base_orders, &atoms) {
            if swapped != base_orders {
                forms.push(apply_orders(mol, &swapped));
            }
        }
    }
    Ok(forms)
}

/// ResonanceRule parent: Kekulé form whose maps 1–2 have the SMARTS-implied order.
pub fn reactant_parent(
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
    smirks: &str,
) -> Result<Molecule, ForestError> {
    let Some(&left) = mapped.get(&1) else {
        return Ok(mol.clone());
    };
    let Some(&right) = mapped.get(&2) else {
        return Ok(mol.clone());
    };
    let Some((_, bond)) = mol.bond_between(atom_idx(left), atom_idx(right)) else {
        return Ok(mol.clone());
    };
    let aromatic = bond.order == BondOrder::Aromatic
        || mol.atom(atom_idx(left)).aromatic
        || mol.atom(atom_idx(right)).aromatic;
    if !aromatic {
        return Ok(mol.clone());
    }
    let want = smirks_mapped_bond_order(smirks).unwrap_or(2.0);
    let want_order = if want >= 1.5 {
        BondOrder::Double
    } else {
        BondOrder::Single
    };
    for form in kekule_forms(mol)? {
        if let Some((_, bond)) = form.bond_between(atom_idx(left), atom_idx(right)) {
            if bond.order == want_order {
                return Ok(form);
            }
        }
    }
    Ok(mol.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, canon_smiles, parse_mol};
    use crate::smarts::smarts_matches;
    use crate::smirks::apply_smirks_at;

    #[test]
    fn benzene_has_two_kekule_forms() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let forms = kekule_forms(&mol).unwrap();
        assert!(forms.len() >= 2);
        let orders: Vec<_> = forms.iter().map(bond_order_map).collect();
        assert_ne!(orders[0], orders[1]);
        for form in &forms {
            let doubles = form
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Double)
                .count();
            assert_eq!(doubles, 3);
        }
    }

    #[test]
    fn thiophene_s_oxidation_picks_single_s_c_parent() {
        let mol = parse_mol("c1ccsc1").unwrap();
        let smarts = "[#6:2]1=,:[#6:3][#6:4]=,:[#6:5][#16;v2,v4:1]1";
        let apply = "[S:1]>>[S+:1][O-]";
        let hits = smarts_matches(&mol, smarts).unwrap();
        assert!(!hits.is_empty());
        let parent = reactant_parent(&mol, &hits[0], &format!("{smarts}>>[S:1]")).unwrap();
        let s = hits[0][&1];
        let c = hits[0][&2];
        let (_, bond) = parent.bond_between(atom_idx(s), atom_idx(c)).unwrap();
        assert_eq!(bond.order, BondOrder::Single);
        let products = apply_smirks_at(apply, &parent, &hits[0]).unwrap();
        assert!(!products.is_empty());
        let want = canon_of("[O-][s+]1cccc1").unwrap();
        assert!(
            products
                .iter()
                .any(|p| canon_of(&canon_smiles(p)).unwrap() == want),
            "want {want}, got {:?}",
            products.iter().map(canon_smiles).collect::<Vec<_>>()
        );
    }

    #[test]
    fn epoxidation_picks_double_parent_on_benzene() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let smirks = "[#6:1]=[#6,#7:2]>>[*:1]1-[*:2][O]1";
        let hits = smarts_matches(&mol, "[#6:1]=,:[#6:2]").unwrap();
        assert_eq!(hits.len(), 6);
        let parent = reactant_parent(&mol, &hits[0], smirks).unwrap();
        let a = hits[0][&1];
        let b = hits[0][&2];
        let (_, bond) = parent.bond_between(atom_idx(a), atom_idx(b)).unwrap();
        assert_eq!(bond.order, BondOrder::Double);
    }
}
