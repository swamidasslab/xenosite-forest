//! ResonancePair path edits driven by endpoint [`PatternInfo`] data.
//!
//! Match each [`Edit::PairEndpoint`] SMARTS, join two anchors by an odd
//! alternating path on a Kekulé form, apply the named end edits, flip the
//! path. Methide is an effect field — two methide ends are allowed.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use chematic::core::{Atom, BondOrder, Element};
use chematic::perception::find_sssr;

use crate::kekule::{conjugated_component, kekule_forms};
use crate::mol::{ForestError, Molecule, aromatize, atom_idx, atom_usize, canon_smiles};
use crate::pattern::{Edit, PatternInfo};
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

fn current_orders(mol: &Molecule) -> HashMap<(usize, usize), i32> {
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

fn system_neighbors(mol: &Molecule, system: &HashSet<usize>) -> HashMap<usize, Vec<usize>> {
    let mut neighbors = HashMap::new();
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
    bond_map: &HashMap<(usize, usize), i32>,
    start: usize,
    end: usize,
    neighbors: &HashMap<usize, Vec<usize>>,
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

fn clear_aromatic(mol: &mut Molecule) {
    for atom in mol.atoms().map(|(idx, _)| idx).collect::<Vec<_>>() {
        if mol.atom(atom).aromatic {
            *mol = mol.with_atom_aromatic(atom, false);
        }
    }
}

/// Capability OR, resolved against whether the conjugated system is aromatic.
fn merge_dearomatizes(left: &PatternInfo, right: &PatternInfo, system_aromatic: bool) -> bool {
    (left.effect.dearomatizes || right.effect.dearomatizes) && system_aromatic
}

/// True when the kekulized system is aromatic again after sanitize.
///
/// Other aromatic systems may stay. The kekulized one did not if any of its
/// aromatic atoms is no longer aromatic, or a non-aromatic double or triple
/// bond still touches it (carbonyl, exocyclic methide). HEURISTICS approved.
fn system_stayed_aromatic(parent: &Molecule, product: &Molecule, system: &HashSet<usize>) -> bool {
    let aromatic_idxs: Vec<usize> = system
        .iter()
        .copied()
        .filter(|&i| parent.atom(atom_idx(i)).aromatic)
        .collect();
    if aromatic_idxs.len() < 2 {
        return false;
    }
    if aromatic_idxs
        .iter()
        .any(|&i| !product.atom(atom_idx(i)).aromatic)
    {
        return false;
    }
    let aromatic_set: HashSet<usize> = aromatic_idxs.into_iter().collect();
    for (_, bond) in product.bonds() {
        let aromatic_bond = bond.order == BondOrder::Aromatic
            || (product.atom(bond.atom1).aromatic && product.atom(bond.atom2).aromatic);
        if aromatic_bond {
            continue;
        }
        match bond.order {
            BondOrder::Double | BondOrder::Triple => {}
            _ => continue,
        }
        let left = atom_usize(bond.atom1);
        let right = atom_usize(bond.atom2);
        if aromatic_set.contains(&left) || aromatic_set.contains(&right) {
            return false;
        }
    }
    true
}

fn ring_sets(mol: &Molecule) -> HashMap<usize, BTreeSet<usize>> {
    let mut out: HashMap<usize, BTreeSet<usize>> = HashMap::new();
    for (ring_i, ring) in find_sssr(mol).rings().iter().enumerate() {
        for &atom in ring {
            out.entry(atom_usize(atom)).or_default().insert(ring_i);
        }
    }
    out
}

fn same_rings(rings: &HashMap<usize, BTreeSet<usize>>, a: usize, b: usize) -> bool {
    rings.get(&a) == rings.get(&b)
}

fn site_atom(mapped: &BTreeMap<u16, usize>, pattern: &PatternInfo) -> Option<usize> {
    mapped.get(&pattern.primary_map()).copied()
}

fn edit_end(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let Edit::PairEndpoint(edit) = &pattern.edit else {
        return false;
    };
    match edit.as_str() {
        "single_to_double" => edit_single_to_double(mol, mapped, pattern, rings),
        "add_carbonyl_o" => edit_add_carbonyl_o(mol, mapped),
        "replace_halogen" => edit_replace_halogen(mol, mapped),
        "iminium" => {
            if !edit_single_to_double(mol, mapped, pattern, rings) {
                return false;
            }
            let Some(&n) = mapped.get(&2) else {
                return false;
            };
            *mol = mol.with_atom_charge(atom_idx(n), 1);
            true
        }
        "dealkylate" => edit_dealkylate(mol, mapped, pattern, rings),
        "keep" => mapped.contains_key(&1),
        _ => false,
    }
}

fn edit_single_to_double(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let (Some(&a), Some(&b)) = (mapped.get(&1), mapped.get(&2)) else {
        return false;
    };
    if pattern.skip_same_rings && same_rings(rings, a, b) {
        return false;
    }
    let Some((bond_idx, _)) = mol.bond_between(atom_idx(a), atom_idx(b)) else {
        return false;
    };
    mol.set_bond_order(bond_idx, BondOrder::Double);
    true
}

fn edit_add_carbonyl_o(mol: &mut Molecule, mapped: &BTreeMap<u16, usize>) -> bool {
    let Some(&carbon) = mapped.get(&1) else {
        return false;
    };
    let (mut next, oxygen) = mol.with_atom_added(Atom::organic(Element::O));
    if next
        .add_bond(atom_idx(carbon), oxygen, BondOrder::Double)
        .is_err()
    {
        return false;
    }
    *mol = next;
    true
}

fn edit_replace_halogen(mol: &mut Molecule, mapped: &BTreeMap<u16, usize>) -> bool {
    let (Some(&carbon), Some(&halogen)) = (mapped.get(&1), mapped.get(&2)) else {
        return false;
    };
    *mol = mol.with_atom_element(atom_idx(halogen), Element::O);
    *mol = mol.with_atom_charge(atom_idx(halogen), 0);
    let Some((bond_idx, _)) = mol.bond_between(atom_idx(carbon), atom_idx(halogen)) else {
        return false;
    };
    mol.set_bond_order(bond_idx, BondOrder::Double);
    true
}

fn edit_dealkylate(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let (Some(&hetero), Some(&alkyl)) = (mapped.get(&2), mapped.get(&3)) else {
        return false;
    };
    if !edit_single_to_double(mol, mapped, pattern, rings) {
        return false;
    }
    let Some((bond_idx, _)) = mol.bond_between(atom_idx(hetero), atom_idx(alkyl)) else {
        return false;
    };
    let _ = bond_idx;
    // Chematic: remove the hetero–alkyl bond by rebuilding without it is heavy;
    // use with_bond_removed if available.
    if let Some((bi, _)) = mol.bond_between(atom_idx(hetero), atom_idx(alkyl)) {
        *mol = mol.with_bond_removed(bi);
    }
    true
}

/// One pair emission before RuleSet packaging.
#[derive(Clone, Debug)]
pub struct PairEmission {
    pub site: usize,
    pub pattern_name: String,
    pub products: Vec<String>,
}

/// Discovered pair site before path flip / product CSMI.
///
/// Carries merged [`crate::pattern::Effect`] so a search can filter without
/// running the edit. [`PairCandidate::materialize`] finds alternating paths
/// and applies end edits.
#[derive(Clone, Debug)]
pub struct PairCandidate {
    pub site: usize,
    pub pattern_name: String,
    pub left: PatternInfo,
    pub right: PatternInfo,
    /// Merged end effects (dearomatizes resolved against system aromaticity).
    pub effect: crate::pattern::Effect,
    map1: BTreeMap<u16, usize>,
    map2: BTreeMap<u16, usize>,
    start: usize,
    end: usize,
    system: HashSet<usize>,
}

impl PairCandidate {
    /// Discovery site atoms for each end (Python `end_atoms`).
    pub fn end_atoms(&self) -> Option<(usize, usize)> {
        let a = site_atom(&self.map1, &self.left)?;
        let b = site_atom(&self.map2, &self.right)?;
        Some((a, b))
    }

    /// Conjugated-system anchors for the alternating path (Python `path_ends`).
    pub fn path_ends(&self) -> (usize, usize) {
        (self.start, self.end)
    }

    pub fn materialize_mols(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        let forms = kekule_forms(mol)?;
        let rings = ring_sets(mol);
        let neighbors = system_neighbors(mol, &self.system);
        let mut products = Vec::new();
        let mut local_csmi = BTreeSet::new();
        for form in &forms {
            let bond_map = current_orders(form);
            let mut paths = alternating_from(&bond_map, self.start, self.end, &neighbors, 2);
            paths.extend(alternating_from(
                &bond_map, self.end, self.start, &neighbors, 2,
            ));
            for path in paths {
                if path.len() % 2 == 1 {
                    continue;
                }
                let mut rw = form.clone();
                clear_aromatic(&mut rw);
                if !edit_end(&mut rw, &self.map1, &self.left, &rings) {
                    continue;
                }
                if !edit_end(&mut rw, &self.map2, &self.right, &rings) {
                    continue;
                }
                if !flip_path(&mut rw, &path) {
                    continue;
                }
                if !accept_product(&rw) {
                    continue;
                }
                let checked = aromatize(&rw);
                if self.effect.dearomatizes && system_stayed_aromatic(mol, &checked, &self.system) {
                    continue;
                }
                let smiles = canon_smiles(&checked);
                if local_csmi.insert(smiles) {
                    products.push(checked);
                }
            }
        }
        Ok(products)
    }

    pub fn materialize(&self, mol: &Molecule) -> Result<Vec<String>, ForestError> {
        Ok(self
            .materialize_mols(mol)?
            .iter()
            .map(canon_smiles)
            .collect())
    }

    /// Elementary plan site atoms (pair ends), falling back to discovery site.
    pub fn plan_site_atoms(&self) -> Vec<usize> {
        match self.end_atoms() {
            Some((a, b)) => vec![a, b],
            None => vec![self.site],
        }
    }

    pub fn emit(&self, mol: &Molecule) -> Result<Option<PairEmission>, ForestError> {
        let products = self.materialize(mol)?;
        if products.is_empty() {
            return Ok(None);
        }
        Ok(Some(PairEmission {
            site: self.site,
            pattern_name: self.pattern_name.clone(),
            products,
        }))
    }
}

type EndpointHit<'a> = (BTreeMap<u16, usize>, &'a PatternInfo);

fn merge_effect_fields(
    left: &PatternInfo,
    right: &PatternInfo,
    system_aromatic: bool,
) -> crate::pattern::Effect {
    let adds = match (&left.effect.adds, &right.effect.adds) {
        (None, None) => None,
        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
        (Some(a), Some(b)) => Some(format!("{a}{b}")),
    };
    let removes = match (&left.effect.removes, &right.effect.removes) {
        (None, None) => None,
        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
        (Some(a), Some(b)) => Some(format!("{a}{b}")),
    };
    crate::pattern::Effect {
        adds,
        removes,
        cleaves: left.effect.cleaves || right.effect.cleaves,
        leave_count: left.effect.leave_count.or(right.effect.leave_count),
        methide: left.effect.methide || right.effect.methide,
        dearomatizes: merge_dearomatizes(left, right, system_aromatic),
        partner: left
            .effect
            .partner
            .clone()
            .or_else(|| right.effect.partner.clone()),
    }
}

/// Discover pair sites without applying path flips.
pub fn pair_candidates(
    mol: &Molecule,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairCandidate>, ForestError> {
    let endpoints: Vec<&PatternInfo> = endpoints
        .iter()
        .filter(|p| matches!(p.edit, Edit::PairEndpoint(_)))
        .collect();
    if endpoints.is_empty() {
        return Ok(Vec::new());
    }

    let mut hits: HashMap<usize, Vec<EndpointHit<'_>>> = HashMap::new();
    for pattern in &endpoints {
        for mapped in smarts_matches(mol, &pattern.smarts)? {
            let Some(&anchor) = mapped.get(&1) else {
                continue;
            };
            hits.entry(anchor).or_default().push((mapped, *pattern));
        }
    }
    if hits.len() < 2 {
        return Ok(Vec::new());
    }

    let mut systems: Vec<HashSet<usize>> = Vec::new();
    let mut covered = HashSet::new();
    for &anchor in hits.keys() {
        if !covered.insert(anchor) {
            continue;
        }
        let (atoms, _) = conjugated_component(mol, anchor);
        let set: HashSet<usize> = atoms.iter().copied().collect();
        covered.extend(&set);
        if hits.keys().filter(|a| set.contains(a)).count() >= 2 {
            systems.push(set);
        }
    }
    if systems.is_empty() {
        let all: HashSet<usize> = mol.atoms().map(|(i, _)| atom_usize(i)).collect();
        if hits.keys().filter(|a| all.contains(a)).count() >= 2 {
            systems.push(all);
        }
    }

    let mut out = Vec::new();
    let mut seen_sig: BTreeSet<(usize, usize, String, String, usize)> = BTreeSet::new();
    let gens = crate::orbits::atom_bond_generators(mol);
    let n_atoms = mol.atom_count();

    for system in &systems {
        let anchors: Vec<usize> = hits
            .keys()
            .copied()
            .filter(|a| system.contains(a))
            .collect();
        if anchors.len() < 2 {
            continue;
        }
        let system_aromatic = system.iter().any(|&i| mol.atom(atom_idx(i)).aromatic);
        let mut sorted_anchors = anchors;
        sorted_anchors.sort_unstable();
        for i in 0..sorted_anchors.len() {
            for j in (i + 1)..sorted_anchors.len() {
                let start = sorted_anchors[i];
                let end = sorted_anchors[j];
                let Some(hits_a) = hits.get(&start) else {
                    continue;
                };
                let Some(hits_b) = hits.get(&end) else {
                    continue;
                };
                for (map1, info1) in hits_a {
                    for (map2, info2) in hits_b {
                        let (Some(site_a), Some(site_b)) =
                            (site_atom(map1, info1), site_atom(map2, info2))
                        else {
                            continue;
                        };
                        if site_a == site_b {
                            continue;
                        }
                        let (n1, n2) = if info1.name <= info2.name {
                            (info1.name.clone(), info2.name.clone())
                        } else {
                            (info2.name.clone(), info1.name.clone())
                        };
                        let sa = site_a.min(site_b);
                        let sb = site_a.max(site_b);
                        let orbit =
                            crate::orbits::atom_pair_orbit_id_with_gens(&gens, n_atoms, sa, sb);
                        if !seen_sig.insert((sa, sb, n1.clone(), n2.clone(), orbit)) {
                            continue;
                        }
                        out.push(PairCandidate {
                            site: sa,
                            pattern_name: format!("{n1}+{n2}"),
                            left: (*info1).clone(),
                            right: (*info2).clone(),
                            effect: merge_effect_fields(info1, info2, system_aromatic),
                            map1: map1.clone(),
                            map2: map2.clone(),
                            start,
                            end,
                            system: system.clone(),
                        });
                    }
                }
            }
        }
    }
    Ok(out)
}

/// Run ResonancePair metabolize for the given endpoint patterns.
pub fn pair_metabolize(
    mol: &Molecule,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairEmission>, ForestError> {
    let mut emissions = Vec::new();
    let mut seen_csmi: BTreeSet<BTreeSet<String>> = BTreeSet::new();
    for candidate in pair_candidates(mol, endpoints)? {
        let Some(emission) = candidate.emit(mol)? else {
            continue;
        };
        let key: BTreeSet<String> = emission.products.iter().cloned().collect();
        if !seen_csmi.insert(key) {
            continue;
        }
        emissions.push(emission);
    }
    Ok(emissions)
}

/// Path dehydrogenation of para-hydroquinone (door smoke for pair metabolize).
pub fn dehydrogenate_hydroquinone(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    let phenol_end = PatternInfo {
        name: "phenol_end".into(),
        smarts: "[#6:1]-[#8H:2]".into(),
        site_kind: crate::pattern::SiteKind::AtomPair,
        site_map: vec![2],
        edit: Edit::PairEndpoint("single_to_double".into()),
        effect: crate::pattern::Effect {
            removes: Some("H".into()),
            dearomatizes: true,
            ..Default::default()
        },
        skip_same_rings: false,
    };
    Ok(pair_metabolize(mol, &[phenol_end])?
        .into_iter()
        .flat_map(|e| e.products)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use crate::rules::{dehydrogenation, quinone_formation};

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

    #[test]
    fn dehydrogenation_rule_pair_door_on_hydroquinone() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        assert!(!endpoints.is_empty());
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| { e.products.iter().any(|p| canon_of(p).unwrap() == want) }),
            "{emissions:?}"
        );
    }

    #[test]
    fn quinone_formation_add_carbonyl_on_benzene() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let endpoints: Vec<_> = quinone_formation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions.iter().any(|e| {
                e.pattern_name.contains("add_carbonyl_o")
                    && e.products.iter().any(|p| canon_of(p).unwrap() == want)
            }),
            "expected para-quinone from two add_carbonyl_o; got {emissions:?}"
        );
    }

    #[test]
    fn pair_candidates_defer_path_flip() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints).unwrap();
        assert!(!cands.is_empty());
        assert!(
            cands.iter().any(|c| c.effect.dearomatizes),
            "aromatic hydroquinone pair should resolve dearomatizes"
        );
        // Refuse before materialize.
        let kept: Vec<_> = cands
            .into_iter()
            .filter(|c| !c.effect.dearomatizes)
            .collect();
        assert!(kept.is_empty());
        // Accept and materialize.
        let cands = pair_candidates(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(cands.iter().any(|c| {
            c.materialize(&mol)
                .unwrap()
                .iter()
                .any(|p| canon_of(p).unwrap() == want)
        }));
    }
}
