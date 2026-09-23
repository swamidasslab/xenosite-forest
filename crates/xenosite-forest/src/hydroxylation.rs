//! Hydroxylation: the first real rule, unique-edit then one-match SMIRKS.

use crate::ForestError;
use crate::mol::{Molecule, canon_smiles};
use crate::smirks::apply_smirks_at;
use crate::unique_edit::unique_atom_sites;

const H: &str = "[#6h1:1]>>[*:1]O";
const H2: &str = "[#6h2,#6h3:1]>>[*:1]O";

/// Unique hydroxylation products as canonical SMILES.
pub fn hydroxylate(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    let mut products = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for (smirks, reactant) in [(H, "[#6h1:1]"), (H2, "[#6h2,#6h3:1]")] {
        for mapped in unique_atom_sites(mol, reactant)? {
            for product in apply_smirks_at(smirks, mol, &mapped)? {
                let smiles = canon_smiles(&product);
                if seen.insert(smiles.clone()) {
                    products.push(smiles);
                }
            }
        }
    }
    Ok(products)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use std::collections::BTreeSet;

    #[test]
    fn ethane_yields_ethanol_once() {
        let mol = parse_mol("CC").unwrap();
        assert_eq!(hydroxylate(&mol).unwrap(), vec!["CCO".to_string()]);
    }

    #[test]
    fn benzene_yields_phenol_once() {
        let mol = parse_mol("c1ccccc1").unwrap();
        assert_eq!(hydroxylate(&mol).unwrap(), vec!["Oc1ccccc1".to_string()]);
    }

    #[test]
    fn propane_yields_primary_and_secondary_alcohols() {
        let mol = parse_mol("CCC").unwrap();
        let got: BTreeSet<String> = hydroxylate(&mol).unwrap().into_iter().collect();
        assert_eq!(
            got,
            BTreeSet::from(["CCCO".to_string(), "CC(C)O".to_string()])
        );
    }
}
