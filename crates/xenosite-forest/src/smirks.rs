//! Apply one SMIRKS match (chematic's replacement for isotope-pinned RunReactants).

use std::collections::BTreeMap;

use chematic::rxn::{apply_reaction_match, find_reaction_matches};

use crate::mol::{ForestError, Molecule, atom_usize};
use crate::valence::accept_product;

/// Run `smirks` on the match whose atom maps equal `mapped`.
///
/// Products are split into connected fragments and dropped when the forest
/// valence gate refuses them.
pub fn apply_smirks_at(
    smirks: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<Vec<Molecule>, ForestError> {
    let matches = find_reaction_matches(smirks, &[mol])
        .map_err(|err| ForestError::Smirks(err.to_string()))?;
    for reaction_match in matches {
        let positions = reaction_match
            .atom_map_positions(smirks)
            .map_err(|err| ForestError::Smirks(err.to_string()))?;
        let same = mapped.iter().all(|(&mapno, &atom)| {
            positions
                .get(&mapno)
                .is_some_and(|(_slot, idx)| atom_usize(*idx) == atom)
        });
        if !same {
            continue;
        }
        let products = apply_reaction_match(smirks, &[mol], &reaction_match, true)
            .map_err(|err| ForestError::Smirks(err.to_string()))?;
        let Some(products) = products else {
            return Ok(Vec::new());
        };
        let mut pieces = Vec::new();
        for product in products {
            for frag in product.fragments() {
                if accept_product(&frag) {
                    pieces.push(frag);
                }
            }
        }
        return Ok(pieces);
    }
    Ok(Vec::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_smiles, parse_mol};
    use crate::smarts::smarts_matches;

    #[test]
    fn hydroxylation_smirks_on_named_atom() {
        let mol = parse_mol("CC").unwrap();
        let hits = smarts_matches(&mol, "[#6h3:1]").unwrap();
        assert_eq!(hits.len(), 2);
        let products = apply_smirks_at("[#6h3:1]>>[*:1]O", &mol, &hits[0]).unwrap();
        assert_eq!(products.len(), 1);
        assert_eq!(canon_smiles(&products[0]), "CCO");
    }

    #[test]
    fn cleavage_smirks_yields_two_fragments() {
        let mol = parse_mol("CN").unwrap();
        let hits = smarts_matches(&mol, "[#6H3:1][#7:2]").unwrap();
        assert_eq!(hits.len(), 1);
        let products =
            apply_smirks_at("[#6H3:1][#7:2]>>([*:2].[*:1](=O)O)", &mol, &hits[0]).unwrap();
        let smiles: Vec<String> = products.iter().map(canon_smiles).collect();
        assert_eq!(smiles.len(), 2);
        assert!(smiles.iter().any(|s| s == "N" || s == "N"));
        assert!(smiles.iter().any(|s| s.contains('O')));
    }
}
