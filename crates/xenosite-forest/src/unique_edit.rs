//! One-site unique-edit: collapse SMARTS hits that share map ranks.

use std::collections::{BTreeMap, BTreeSet};

use crate::ForestError;
use crate::mol::{Molecule, ranks};
use crate::smarts::smarts_matches;

pub fn map_rank_key(ranks: &[usize], mapped: &BTreeMap<u16, usize>) -> Vec<(u16, usize)> {
    mapped
        .iter()
        .map(|(&mapno, &idx)| (mapno, ranks[idx]))
        .collect()
}

/// Unique atom sites for a SMARTS, keyed by `(mapno, topological rank)`.
pub fn unique_atom_sites(
    mol: &Molecule,
    smarts: &str,
) -> Result<Vec<BTreeMap<u16, usize>>, ForestError> {
    let ranks = ranks(mol);
    let mut seen = BTreeSet::new();
    let mut unique = Vec::new();
    for mapped in smarts_matches(mol, smarts)? {
        let key = map_rank_key(&ranks, &mapped);
        if seen.insert(key) {
            unique.push(mapped);
        }
    }
    Ok(unique)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn benzene_hydroxylation_is_one_site() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let sites = unique_atom_sites(&mol, "[#6h1:1]").unwrap();
        assert_eq!(sites.len(), 1);
    }

    #[test]
    fn propane_has_two_unique_carbons() {
        let mol = parse_mol("CCC").unwrap();
        let sites = unique_atom_sites(&mol, "[#6h2,#6h3:1]").unwrap();
        assert_eq!(sites.len(), 2);
    }

    #[test]
    fn ethane_has_one_unique_carbon() {
        let mol = parse_mol("CC").unwrap();
        let sites = unique_atom_sites(&mol, "[#6h3:1]").unwrap();
        assert_eq!(sites.len(), 1);
    }
}
