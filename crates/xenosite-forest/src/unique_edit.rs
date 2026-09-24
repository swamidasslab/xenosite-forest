//! One-site unique-edit: collapse SMARTS hits that share map ranks.

use std::collections::{BTreeMap, BTreeSet};

use crate::ForestError;
use crate::mol::{Molecule, ranks};
use crate::pattern::SiteKind;
use crate::smarts::smarts_matches;

pub fn map_rank_key(ranks: &[usize], mapped: &BTreeMap<u16, usize>) -> Vec<(u16, usize)> {
    mapped
        .iter()
        .map(|(&mapno, &idx)| (mapno, ranks[idx]))
        .collect()
}

/// Unordered bond ends as sorted topological ranks (`site_kind="bond"`).
pub fn bond_rank_key(ranks: &[usize], site_atoms: &[usize]) -> Vec<usize> {
    let mut key: Vec<usize> = site_atoms.iter().map(|&i| ranks[i]).collect();
    key.sort_unstable();
    key
}

fn site_atoms(mapped: &BTreeMap<u16, usize>, site_map: &[u16]) -> Vec<usize> {
    site_map
        .iter()
        .filter_map(|m| mapped.get(m).copied())
        .collect()
}

/// Dedup key for a SMARTS hit. Bond / atom_pair use undirected ranks.
pub fn site_rank_key(
    ranks: &[usize],
    mapped: &BTreeMap<u16, usize>,
    site_kind: SiteKind,
    site_map: &[u16],
) -> Vec<(u16, usize)> {
    match site_kind {
        SiteKind::Bond | SiteKind::AtomPair => {
            let atoms = site_atoms(mapped, site_map);
            bond_rank_key(ranks, &atoms)
                .into_iter()
                .enumerate()
                .map(|(i, r)| (i as u16, r))
                .collect()
        }
        SiteKind::Atom | SiteKind::DirectedBond => map_rank_key(ranks, mapped),
    }
}

/// Unique atom sites for a SMARTS, keyed by `(mapno, topological rank)`.
pub fn unique_atom_sites(
    mol: &Molecule,
    smarts: &str,
) -> Result<Vec<BTreeMap<u16, usize>>, ForestError> {
    unique_sites(mol, smarts, SiteKind::Atom, &[1])
}

/// Unique sites for a pattern's SMARTS / site kind.
pub fn unique_sites(
    mol: &Molecule,
    smarts: &str,
    site_kind: SiteKind,
    site_map: &[u16],
) -> Result<Vec<BTreeMap<u16, usize>>, ForestError> {
    let ranks = ranks(mol);
    let mut seen = BTreeSet::new();
    let mut unique = Vec::new();
    for mapped in smarts_matches(mol, smarts)? {
        let key = site_rank_key(&ranks, &mapped, site_kind, site_map);
        if seen.insert(key) {
            unique.push(mapped);
        }
    }
    Ok(unique)
}

/// Unique sites across Kekulé forms, keyed on the aromatic parent's ranks.
///
/// Shared `seen` across forms so a bond that is double in more than one
/// writing is one site (Python Epoxidation).
pub type FormSiteHit = (BTreeMap<u16, usize>, usize);

pub fn unique_sites_on_forms(
    context: &Molecule,
    forms: &[Molecule],
    smarts: &str,
    site_kind: SiteKind,
    site_map: &[u16],
) -> Result<Vec<FormSiteHit>, ForestError> {
    let ranks = ranks(context);
    let mut seen = BTreeSet::new();
    let mut unique = Vec::new();
    for (form_i, form) in forms.iter().enumerate() {
        for mapped in smarts_matches(form, smarts)? {
            let key = site_rank_key(&ranks, &mapped, site_kind, site_map);
            if seen.insert(key) {
                unique.push((mapped, form_i));
            }
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

    #[test]
    fn bond_rank_key_is_undirected() {
        let ranks = vec![0, 1, 2];
        assert_eq!(
            bond_rank_key(&ranks, &[0, 1]),
            bond_rank_key(&ranks, &[1, 0])
        );
    }
}
