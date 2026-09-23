//! ResonancePair path edit: hydroquinone → benzoquinone.

use std::collections::{HashSet, VecDeque};

use chematic::core::BondOrder;

use crate::kekule::kekule_forms;
use crate::mol::{ForestError, Molecule, atom_idx, atom_usize, canon_smiles};
use crate::smarts::smarts_matches;
use crate::valence::accept_product;

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn order_value(order: BondOrder) -> Option<i32> {
    match order {
        BondOrder::Single => Some(1),
        BondOrder::Double => Some(2),
        _ => None,
    }
}

fn current_orders(mol: &Molecule) -> std::collections::HashMap<(usize, usize), i32> {
    mol.bonds()
        .filter_map(|(_, bond)| {
            order_value(bond.order).map(|order| {
                (
                    bond_key(atom_usize(bond.atom1), atom_usize(bond.atom2)),
                    order,
                )
            })
        })
        .collect()
}

fn system_neighbors(
    mol: &Molecule,
    system: &HashSet<usize>,
) -> std::collections::HashMap<usize, Vec<usize>> {
    let mut neighbors = std::collections::HashMap::new();
    for &i in system {
        let mut nbrs = Vec::new();
        for (nbr, _) in mol.neighbors(atom_idx(i)) {
            let j = atom_usize(nbr);
            if system.contains(&j) {
                nbrs.push(j);
            }
        }
        neighbors.insert(i, nbrs);
    }
    neighbors
}

fn alternating_from(
    bond_map: &std::collections::HashMap<(usize, usize), i32>,
    start: usize,
    end: usize,
    neighbors: &std::collections::HashMap<usize, Vec<usize>>,
    first: i32,
) -> Vec<Vec<usize>> {
    let mut queue = VecDeque::from([(start, first, vec![start])]);
    let mut seen = HashSet::from([(start, first)]);
    let mut found = Vec::new();
    while let Some((node, want, path)) = queue.pop_front() {
        let next_want = if want == 2 { 1 } else { 2 };
        let Some(nbrs) = neighbors.get(&node) else {
            continue;
        };
        for &nbr in nbrs {
            if path.contains(&nbr) {
                continue;
            }
            let Some(&order) = bond_map.get(&bond_key(node, nbr)) else {
                continue;
            };
            if order != want {
                continue;
            }
            let mut nxt = path.clone();
            nxt.push(nbr);
            if nbr == end {
                found.push(nxt);
                continue;
            }
            if seen.insert((nbr, next_want)) {
                queue.push_back((nbr, next_want, nxt));
            }
        }
    }
    found
}

fn flip_path(mol: &mut Molecule, path: &[usize]) -> bool {
    let mut flipped = false;
    for window in path.windows(2) {
        let a = atom_idx(window[0]);
        let b = atom_idx(window[1]);
        let Some((bond_idx, bond)) = mol.bond_between(a, b) else {
            return false;
        };
        match bond.order {
            BondOrder::Double => {
                mol.set_bond_order(bond_idx, BondOrder::Single);
                flipped = true;
            }
            BondOrder::Single => {
                mol.set_bond_order(bond_idx, BondOrder::Double);
                flipped = true;
            }
            _ => {}
        }
    }
    flipped
}

/// Path dehydrogenation of para-hydroquinone. Derisks ResonancePairRule.
pub fn dehydrogenate_hydroquinone(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    let hits = smarts_matches(mol, "[#6:1][#8H1:2]")?;
    if hits.len() < 2 {
        return Ok(Vec::new());
    }
    let carbons: Vec<usize> = hits.iter().map(|h| h[&1]).collect();
    let oxygens: Vec<usize> = hits.iter().map(|h| h[&2]).collect();
    let mut products = Vec::new();
    let mut seen = std::collections::BTreeSet::new();

    for form in kekule_forms(mol)? {
        let aromatic: HashSet<usize> = form
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(atom_usize(idx)))
            .chain(form.atoms().map(|(idx, _)| atom_usize(idx)))
            .collect();
        let neighbors = system_neighbors(&form, &aromatic);
        let bond_map = current_orders(&form);
        for i in 0..carbons.len() {
            for j in (i + 1)..carbons.len() {
                let start = carbons[i];
                let end = carbons[j];
                let mut paths = alternating_from(&bond_map, start, end, &neighbors, 2);
                paths.extend(alternating_from(&bond_map, end, start, &neighbors, 2));
                for path in paths {
                    // Odd bond count ⇔ even atom count. Even walks are not a pair flip.
                    if path.len() % 2 == 1 {
                        continue;
                    }
                    let mut rw = form.clone();
                    for atom in rw.atoms().map(|(idx, _)| idx).collect::<Vec<_>>() {
                        if rw.atom(atom).aromatic {
                            rw = rw.with_atom_aromatic(atom, false);
                        }
                    }
                    if !flip_path(&mut rw, &path) {
                        continue;
                    }
                    for &oxygen_i in &[i, j] {
                        let carbon = atom_idx(carbons[oxygen_i]);
                        let oxygen = atom_idx(oxygens[oxygen_i]);
                        if let Some((bond_idx, _)) = rw.bond_between(carbon, oxygen) {
                            rw.set_bond_order(bond_idx, BondOrder::Double);
                        }
                    }
                    if !accept_product(&rw) {
                        continue;
                    }
                    let smiles = canon_smiles(&rw);
                    if seen.insert(smiles.clone()) {
                        products.push(smiles);
                    }
                }
            }
        }
    }
    Ok(products)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};

    #[test]
    fn hydroquinone_yields_benzoquinone() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let products = dehydrogenate_hydroquinone(&mol).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            products.iter().any(|s| canon_of(s).unwrap() == want),
            "want {want}, got {products:?}"
        );
    }
}
