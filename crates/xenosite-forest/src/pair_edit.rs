//! ResonancePair path edits driven by endpoint [`PatternInfo`] data.
//!
//! Match each [`Edit::PairEndpoint`] SMARTS, join two anchors by an odd
//! alternating path on a Kekulé form, apply the named end edits, flip the
//! path. Methide is an effect field — two methide ends are allowed.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use chematic::core::{Atom, BondOrder, Element};
use chematic::perception::find_sssr;

use crate::kekule::{conjugated_component, kekule_forms};
use crate::mol::{ForestError, Molecule, atom_idx, atom_usize, canon_smiles};
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

fn system_neighbors(
    mol: &Molecule,
    system: &HashSet<usize>,
) -> HashMap<usize, Vec<usize>> {
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

type EndpointHit<'a> = (BTreeMap<u16, usize>, &'a PatternInfo);

/// Run ResonancePair metabolize for the given endpoint patterns.
pub fn pair_metabolize(
    mol: &Molecule,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairEmission>, ForestError> {
    let endpoints: Vec<&PatternInfo> = endpoints
        .iter()
        .filter(|p| matches!(p.edit, Edit::PairEndpoint(_)))
        .collect();
    if endpoints.is_empty() {
        return Ok(Vec::new());
    }

    // hits[map1_atom] = list of (mapped, pattern)
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

    // Conjugated systems covering at least two anchors.
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
        // Fall back: whole-molecule atom set for small systems.
        let all: HashSet<usize> = mol.atoms().map(|(i, _)| atom_usize(i)).collect();
        if hits.keys().filter(|a| all.contains(a)).count() >= 2 {
            systems.push(all);
        }
    }

    let forms = kekule_forms(mol)?;
    let rings = ring_sets(mol);
    let mut emissions = Vec::new();
    let mut seen_sig: BTreeSet<(usize, usize, String, String)> = BTreeSet::new();
    let mut seen_csmi: BTreeSet<BTreeSet<String>> = BTreeSet::new();

    for system in &systems {
        let anchors: Vec<usize> = hits.keys().copied().filter(|a| system.contains(a)).collect();
        if anchors.len() < 2 {
            continue;
        }
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
                        if !seen_sig.insert((sa, sb, n1.clone(), n2.clone())) {
                            continue;
                        }
                        let neighbors = system_neighbors(mol, system);
                        let mut products = Vec::new();
                        let mut local_csmi = BTreeSet::new();
                        for form in &forms {
                            let bond_map = current_orders(form);
                            let mut paths = alternating_from(&bond_map, start, end, &neighbors, 2);
                            paths.extend(alternating_from(
                                &bond_map, end, start, &neighbors, 2,
                            ));
                            for path in paths {
                                // Odd bond count ⇔ even atom count.
                                if path.len() % 2 == 1 {
                                    continue;
                                }
                                let mut rw = form.clone();
                                clear_aromatic(&mut rw);
                                if !edit_end(&mut rw, map1, info1, &rings) {
                                    continue;
                                }
                                if !edit_end(&mut rw, map2, info2, &rings) {
                                    continue;
                                }
                                if !flip_path(&mut rw, &path) {
                                    continue;
                                }
                                if !accept_product(&rw) {
                                    continue;
                                }
                                let smiles = canon_smiles(&rw);
                                if local_csmi.insert(smiles.clone()) {
                                    products.push(smiles);
                                }
                            }
                        }
                        if products.is_empty() {
                            continue;
                        }
                        let key: BTreeSet<String> = products.iter().cloned().collect();
                        if !seen_csmi.insert(key) {
                            continue;
                        }
                        emissions.push(PairEmission {
                            site: sa,
                            pattern_name: format!("{n1}+{n2}"),
                            products,
                        });
                    }
                }
            }
        }
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
            .cloned()
            .filter(|p| matches!(p.edit, Edit::PairEndpoint(_)))
            .collect();
        assert!(!endpoints.is_empty());
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions.iter().any(|e| {
                e.products
                    .iter()
                    .any(|p| canon_of(p).unwrap() == want)
            }),
            "{emissions:?}"
        );
    }

    #[test]
    fn quinone_formation_add_carbonyl_on_benzene() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let endpoints: Vec<_> = quinone_formation()
            .patterns()
            .into_iter()
            .cloned()
            .filter(|p| matches!(p.edit, Edit::PairEndpoint(_)))
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
}
