//! Hydroxylation: unique-edit then add OH (graph edit, SMIRKS dialect aside).
//!
//! Patterns live on [`crate::rules::hydroxylation`].

use crate::ForestMol;
use crate::mol::ForestError;
use crate::ruleset::{accept_all_rules, accept_all_sites};

pub use crate::rules::hydroxylation;

/// Unique hydroxylation products as canonical SMILES (explicit CSMI downgrade).
pub fn hydroxylate(mol: &ForestMol) -> Result<Vec<String>, ForestError> {
    Ok(hydroxylation()
        .metabolize(mol, accept_all_rules, accept_all_sites, true)
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flat_map(|emission| emission.product_csmis())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::canon_of;
    use std::collections::BTreeSet;

    fn canon_set(smiles: impl IntoIterator<Item = impl AsRef<str>>) -> BTreeSet<String> {
        smiles
            .into_iter()
            .map(|s| canon_of(s.as_ref()).unwrap())
            .collect()
    }

    #[test]
    fn ethane_yields_ethanol_once() {
        let mol = ForestMol::parse("CC").unwrap();
        assert_eq!(canon_set(hydroxylate(&mol).unwrap()), canon_set(["CCO"]));
    }

    #[test]
    fn benzene_yields_phenol_once() {
        let mol = ForestMol::parse("c1ccccc1").unwrap();
        assert_eq!(
            canon_set(hydroxylate(&mol).unwrap()),
            canon_set(["Oc1ccccc1"])
        );
    }

    #[test]
    fn propane_yields_primary_and_secondary_alcohols() {
        let mol = ForestMol::parse("CCC").unwrap();
        assert_eq!(
            canon_set(hydroxylate(&mol).unwrap()),
            canon_set(["CCCO", "CC(C)O"])
        );
    }
}
