//! Hydroxylation: unique-edit then add OH (graph edit, SMIRKS dialect aside).

use chematic::core::{Atom, BondOrder, Element};

use crate::mol::{ForestError, Molecule, atom_idx, canon_smiles};
use crate::unique_edit::unique_atom_sites;
use crate::valence::accept_product;

const H: &str = "[#6h1:1]";
const H2: &str = "[#6h2,#6h3:1]";

fn add_hydroxyl(mol: &Molecule, carbon: usize) -> Result<Molecule, ForestError> {
    let (mut product, oxygen) = mol.with_atom_added(Atom::organic(Element::O));
    product
        .add_bond(atom_idx(carbon), oxygen, BondOrder::Single)
        .map_err(|err| ForestError::Smirks(err.to_string()))?;
    Ok(product)
}

/// Unique hydroxylation products as canonical SMILES.
pub fn hydroxylate(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    let mut products = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for reactant in [H, H2] {
        for mapped in unique_atom_sites(mol, reactant)? {
            let Some(&carbon) = mapped.get(&1) else {
                continue;
            };
            let product = add_hydroxyl(mol, carbon)?;
            if !accept_product(&product) {
                continue;
            }
            let smiles = canon_smiles(&product);
            if seen.insert(smiles.clone()) {
                products.push(smiles);
            }
        }
    }
    Ok(products)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use std::collections::BTreeSet;

    fn canon_set(smiles: impl IntoIterator<Item = impl AsRef<str>>) -> BTreeSet<String> {
        smiles
            .into_iter()
            .map(|s| canon_of(s.as_ref()).unwrap())
            .collect()
    }

    #[test]
    fn ethane_yields_ethanol_once() {
        let mol = parse_mol("CC").unwrap();
        assert_eq!(canon_set(hydroxylate(&mol).unwrap()), canon_set(["CCO"]));
    }

    #[test]
    fn benzene_yields_phenol_once() {
        let mol = parse_mol("c1ccccc1").unwrap();
        assert_eq!(
            canon_set(hydroxylate(&mol).unwrap()),
            canon_set(["Oc1ccccc1"])
        );
    }

    #[test]
    fn propane_yields_primary_and_secondary_alcohols() {
        let mol = parse_mol("CCC").unwrap();
        assert_eq!(
            canon_set(hydroxylate(&mol).unwrap()),
            canon_set(["CCCO", "CC(C)O"])
        );
    }
}
