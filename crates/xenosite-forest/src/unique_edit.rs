//! One-site unique-edit: collapse SMARTS hits that share map ranks.
//!
//! Still emits one representative per class. The collapsed orbit of primary
//! site atoms is passed down on [`UniqueSite::orbit`] / [`crate::pattern::SiteInfo`]
//! so multiplicity is visible and plans can test site-class equivalence.

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

/// One unique-edit survivor plus the primary-map orbit that collapsed into it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UniqueSite {
    pub mapped: BTreeMap<u16, usize>,
    /// Primary-map atom indexes that share this unique-edit key (sorted).
    /// Length is the topological multiplicity; the representative is
    /// [`Self::mapped`]'s primary atom (also in this list).
    pub orbit: Vec<usize>,
}

impl UniqueSite {
    pub fn primary(&self, site_map: &[u16]) -> Option<usize> {
        let map = site_map.first().copied().unwrap_or(1);
        self.mapped.get(&map).copied()
    }
}

/// True when `a` and `b` name the same unique-edit class given passed orbits.
///
/// Plan equivalence under automorphism reads this: two steps with the same
/// rule are site-equivalent when each primary lies in the other's orbit (or
/// they are equal). Orbits of length 1 still work (`[a]` contains only `a`).
pub fn same_site_orbit(a: usize, a_orbit: &[usize], b: usize, b_orbit: &[usize]) -> bool {
    if a == b {
        return true;
    }
    a_orbit.contains(&b) || b_orbit.contains(&a)
}

fn primary_atom(mapped: &BTreeMap<u16, usize>, site_map: &[u16]) -> Option<usize> {
    let map = site_map.first().copied().unwrap_or(1);
    mapped.get(&map).copied()
}

type RankKey = Vec<(u16, usize)>;
type OrbitGroup = (BTreeMap<u16, usize>, BTreeSet<usize>);
type FormOrbitGroup = (BTreeMap<u16, usize>, BTreeSet<usize>, usize);

fn finish_groups(
    order: Vec<RankKey>,
    mut groups: BTreeMap<RankKey, OrbitGroup>,
) -> Vec<UniqueSite> {
    order
        .into_iter()
        .filter_map(|key| {
            let (mapped, orbit) = groups.remove(&key)?;
            Some(UniqueSite {
                mapped,
                orbit: orbit.into_iter().collect(),
            })
        })
        .collect()
}

fn record_hit(
    order: &mut Vec<RankKey>,
    groups: &mut BTreeMap<RankKey, OrbitGroup>,
    key: RankKey,
    mapped: BTreeMap<u16, usize>,
    primary: Option<usize>,
) {
    if let Some((_, orbit)) = groups.get_mut(&key) {
        if let Some(p) = primary {
            orbit.insert(p);
        }
        return;
    }
    let mut orbit = BTreeSet::new();
    if let Some(p) = primary {
        orbit.insert(p);
    }
    order.push(key.clone());
    groups.insert(key, (mapped, orbit));
}

/// Unique atom sites for a SMARTS, keyed by `(mapno, topological rank)`.
pub fn unique_atom_sites(
    mol: &Molecule,
    smarts: &str,
) -> Result<Vec<BTreeMap<u16, usize>>, ForestError> {
    Ok(unique_sites(mol, smarts, SiteKind::Atom, &[1])?
        .into_iter()
        .map(|u| u.mapped)
        .collect())
}

/// Unique atom sites with collapsed orbits (benzene → one site, orbit len 6).
pub fn unique_atom_sites_with_orbits(
    mol: &Molecule,
    smarts: &str,
) -> Result<Vec<UniqueSite>, ForestError> {
    unique_sites(mol, smarts, SiteKind::Atom, &[1])
}

/// Unique sites for a pattern's SMARTS / site kind.
pub fn unique_sites(
    mol: &Molecule,
    smarts: &str,
    site_kind: SiteKind,
    site_map: &[u16],
) -> Result<Vec<UniqueSite>, ForestError> {
    let ranks = ranks(mol);
    let mut order = Vec::new();
    let mut groups: BTreeMap<RankKey, OrbitGroup> = BTreeMap::new();
    for mapped in smarts_matches(mol, smarts)? {
        let key = site_rank_key(&ranks, &mapped, site_kind, site_map);
        let primary = primary_atom(&mapped, site_map);
        record_hit(&mut order, &mut groups, key, mapped, primary);
    }
    Ok(finish_groups(order, groups))
}

/// Unique sites across Kekulé forms, keyed on the aromatic parent's ranks.
///
/// Shared `seen` across forms so a bond that is double in more than one
/// writing is one site (Python Epoxidation).
pub type FormSiteHit = (UniqueSite, usize);

pub fn unique_sites_on_forms(
    context: &Molecule,
    forms: &[Molecule],
    smarts: &str,
    site_kind: SiteKind,
    site_map: &[u16],
) -> Result<Vec<FormSiteHit>, ForestError> {
    let ranks = ranks(context);
    let mut order = Vec::new();
    let mut groups: BTreeMap<RankKey, FormOrbitGroup> = BTreeMap::new();
    for (form_i, form) in forms.iter().enumerate() {
        for mapped in smarts_matches(form, smarts)? {
            let key = site_rank_key(&ranks, &mapped, site_kind, site_map);
            let primary = primary_atom(&mapped, site_map);
            if let Some((_, orbit, _)) = groups.get_mut(&key) {
                if let Some(p) = primary {
                    orbit.insert(p);
                }
                continue;
            }
            let mut orbit = BTreeSet::new();
            if let Some(p) = primary {
                orbit.insert(p);
            }
            order.push(key.clone());
            groups.insert(key, (mapped, orbit, form_i));
        }
    }
    Ok(order
        .into_iter()
        .filter_map(|key| {
            let (mapped, orbit, form_i) = groups.remove(&key)?;
            Some((
                UniqueSite {
                    mapped,
                    orbit: orbit.into_iter().collect(),
                },
                form_i,
            ))
        })
        .collect())
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
    fn benzene_orbit_has_six_carbons() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let sites = unique_atom_sites_with_orbits(&mol, "[#6h1:1]").unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].orbit.len(), 6);
        let primary = sites[0].primary(&[1]).unwrap();
        assert!(sites[0].orbit.contains(&primary));
    }

    #[test]
    fn ethane_orbit_has_two_carbons() {
        let mol = parse_mol("CC").unwrap();
        let sites = unique_atom_sites_with_orbits(&mol, "[#6h3:1]").unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].orbit.len(), 2);
    }

    #[test]
    fn propane_has_two_unique_carbons() {
        let mol = parse_mol("CCC").unwrap();
        let sites = unique_atom_sites(&mol, "[#6h2,#6h3:1]").unwrap();
        assert_eq!(sites.len(), 2);
    }

    #[test]
    fn propane_methyls_share_an_orbit_of_two() {
        let mol = parse_mol("CCC").unwrap();
        let sites = unique_atom_sites_with_orbits(&mol, "[#6h2,#6h3:1]").unwrap();
        assert_eq!(sites.len(), 2);
        let mut sizes: Vec<_> = sites.iter().map(|s| s.orbit.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1, 2]);
    }

    #[test]
    fn ethane_has_one_unique_carbon() {
        let mol = parse_mol("CC").unwrap();
        let sites = unique_atom_sites(&mol, "[#6h3:1]").unwrap();
        assert_eq!(sites.len(), 1);
    }

    #[test]
    fn same_site_orbit_reads_passed_orbits() {
        assert!(same_site_orbit(0, &[0, 1, 2], 2, &[2]));
        assert!(same_site_orbit(2, &[2], 0, &[0, 1, 2]));
        assert!(!same_site_orbit(0, &[0, 1], 3, &[3, 4]));
        assert!(same_site_orbit(5, &[5], 5, &[5]));
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
