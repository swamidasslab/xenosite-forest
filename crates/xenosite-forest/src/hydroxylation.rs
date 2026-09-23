//! Hydroxylation: unique-edit then add OH (graph edit, SMIRKS dialect aside).
//!
//! The patterns are [`PatternInfo`] data. [`hydroxylate`] runs them through
//! [`crate::ruleset::RuleSet::metabolize`].

use crate::mol::{ForestError, Molecule};
use crate::pattern::PatternInfo;
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};

const H: &str = "[#6h1:1]";
const H2: &str = "[#6h2,#6h3:1]";

/// Hydroxylation as a [`RuleSet`]: two `PatternInfo` rows, `h` then `h2`.
pub fn hydroxylation() -> RuleSet {
    RuleSet::new(
        Some("Hydroxylation".into()),
        [
            PatternInfo::hydroxyl("h", H),
            PatternInfo::hydroxyl("h2", H2),
        ],
    )
}

/// Unique hydroxylation products as canonical SMILES.
pub fn hydroxylate(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    Ok(hydroxylation()
        .metabolize(mol, accept_all_rules, accept_all_sites, true)?
        .into_iter()
        .flat_map(|emission| emission.products)
        .collect())
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
