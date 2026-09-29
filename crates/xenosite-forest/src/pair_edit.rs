//! ResonancePair path edits driven by endpoint [`PatternInfo`] data.
//!
//! Match each [`Edit::PairEndpoint`] SMARTS, join two anchors by an odd
//! alternating path on a Kekulé form, apply the named end edits, flip the
//! path. Methide is an effect field — two methide ends are allowed.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::rc::Rc;

use chematic::core::{Atom, BondOrder, Element};
use chematic::perception::find_sssr;

use crate::kekule::{conjugated_component_ext, kekule_forms};
use crate::mol::{ForestError, Molecule, aromatize, atom_idx, atom_usize, canon_smiles};
use crate::pattern::{Edit, PatternInfo};
use crate::smarts::smarts_matches;
use crate::valence::accept_product;
use crate::ForestMol;

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

/// Legacy Tautomerization gate: no atom may carry two double bonds (drops allenes).
fn at_most_one_double_per_atom(mol: &Molecule) -> bool {
    for (idx, _) in mol.atoms() {
        let mut doubles = 0;
        for (_, bond_idx) in mol.neighbors(idx) {
            if mol.bond(bond_idx).order == BondOrder::Double {
                doubles += 1;
                if doubles > 1 {
                    return false;
                }
            }
        }
    }
    true
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

fn is_pair_endpoint(pattern: &PatternInfo, name: &str) -> bool {
    matches!(&pattern.edit, Edit::PairEndpoint(e) if e == name)
}

fn is_tautomer_extend(pattern: &PatternInfo) -> bool {
    is_pair_endpoint(pattern, "tautomer_extend")
}

fn is_tautomer_far(pattern: &PatternInfo) -> bool {
    is_pair_endpoint(pattern, "tautomer_far")
}

/// H-donor atom (map 2) for a tautomer_extend endpoint, if mapped.
fn tautomer_h_donor(mapped: &BTreeMap<u16, usize>, pattern: &PatternInfo) -> Option<usize> {
    if !is_tautomer_extend(pattern) {
        return None;
    }
    mapped.get(&2).copied()
}

/// Prepend or append `h_donor` onto an alternating path that ends at `anchor`.
fn extend_path_with_h(path: &[usize], anchor: usize, h_donor: usize) -> Option<Vec<usize>> {
    if path.is_empty() || path.contains(&h_donor) {
        return None;
    }
    if path[0] == anchor {
        let mut full = Vec::with_capacity(path.len() + 1);
        full.push(h_donor);
        full.extend_from_slice(path);
        Some(full)
    } else if path[path.len() - 1] == anchor {
        let mut full = path.to_vec();
        full.push(h_donor);
        Some(full)
    } else {
        None
    }
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
        // Path flip carries the tautomer; end edits only validate maps.
        "tautomer_extend" => mapped.contains_key(&1) && mapped.contains_key(&2),
        "tautomer_far" | "keep" => mapped.contains_key(&1),
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

/// Soft rename — pair SOMs are [`crate::candidate::DeferredSite`].
pub type PairCandidate = crate::candidate::DeferredSite;

/// Apply path flip for a resonance-pair match (used by [`DeferredSite::materialize_mols`]).
///
/// Tautomer pairs (`tautomer_extend` on exactly one end) reuse the same
/// alternating-path discovery and [`flip_path`], after extending the path by
/// the H-donor (map 2). Odd atom counts are allowed only for that extension.
pub(crate) fn materialize_pair_mols(
    mol: &Molecule,
    left: &PatternInfo,
    right: &PatternInfo,
    effect: &crate::pattern::Effect,
    map1: &BTreeMap<u16, usize>,
    map2: &BTreeMap<u16, usize>,
    start: usize,
    end: usize,
    system: &HashSet<usize>,
) -> Result<Vec<Molecule>, ForestError> {
    let left_h = tautomer_h_donor(map1, left);
    let right_h = tautomer_h_donor(map2, right);
    // PatternInfo edits decide the door: exactly one tautomer_extend + one tautomer_far.
    let tautomer_anchor_h = match (left_h, right_h, is_tautomer_far(left), is_tautomer_far(right)) {
        (Some(h), None, false, true) => Some((map1.get(&1).copied(), h)),
        (None, Some(h), true, false) => Some((map2.get(&1).copied(), h)),
        _ => None,
    };
    let is_tautomer = tautomer_anchor_h.is_some();
    // extend×extend, far×far, or extend paired with a non-far end: not a tautomer site.
    if (is_tautomer_extend(left) || is_tautomer_extend(right) || is_tautomer_far(left) || is_tautomer_far(right))
        && !is_tautomer
    {
        return Ok(Vec::new());
    }

    let forms = kekule_forms(mol)?;
    let rings = ring_sets(mol);
    let neighbors = system_neighbors(mol, system);
    let mut products = Vec::new();
    let mut local_csmi = BTreeSet::new();
    for form in &forms {
        let bond_map = current_orders(form);
        let mut paths = alternating_from(&bond_map, start, end, &neighbors, 2);
        paths.extend(alternating_from(&bond_map, end, start, &neighbors, 2));
        for path in paths {
            // Ordinary ResonancePair: even atom count (odd bonds). Tautomer
            // PatternInfo edits allow either parity — extension by the H-donor
            // makes the flipped path; long-range enol↔ketone needs odd conjugated paths.
            if !is_tautomer && path.len() % 2 == 1 {
                continue;
            }
            let flip_atoms = if let Some((Some(anchor), h_donor)) = tautomer_anchor_h {
                match extend_path_with_h(&path, anchor, h_donor) {
                    Some(full) => full,
                    None => continue,
                }
            } else if is_tautomer {
                continue;
            } else {
                path.clone()
            };
            let mut rw = form.clone();
            clear_aromatic(&mut rw);
            if !edit_end(&mut rw, map1, left, &rings) {
                continue;
            }
            if !edit_end(&mut rw, map2, right, &rings) {
                continue;
            }
            if !flip_path(&mut rw, &flip_atoms) {
                continue;
            }
            if is_tautomer {
                // Path flip relocates H via bond orders; bracket H (e.g. [nH])
                // must not stick on the donor or the product fails valence
                // (imine→amine on tacrine). Clear path atoms so implicit H
                // recomputes — same as RDKit SanitizeMol after swap_bonds.
                for &i in &flip_atoms {
                    rw.set_hydrogen_count(atom_idx(i), None);
                }
            }
            if is_tautomer && !at_most_one_double_per_atom(&rw) {
                continue;
            }
            if !accept_product(&rw) {
                continue;
            }
            let checked = aromatize(&rw);
            if effect.dearomatizes && system_stayed_aromatic(mol, &checked, system) {
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
    let left_delta = left.effect.resolved_delta_formula();
    let right_delta = right.effect.resolved_delta_formula();
    crate::pattern::Effect {
        adds,
        removes,
        delta_formula: crate::pattern::merge_delta_formula(&left_delta, &right_delta),
        leave_formula: {
            let mut leave = left.effect.leave_formula.clone();
            for (el, n) in &right.effect.leave_formula {
                *leave.entry(el.clone()).or_insert(0) += n;
            }
            leave.retain(|_, n| *n != 0);
            leave
        },
        cleaves: left.effect.cleaves || right.effect.cleaves,
        leave_count: left.effect.leave_count.or(right.effect.leave_count),
        methide: left.effect.methide || right.effect.methide,
        dearomatizes: merge_dearomatizes(left, right, system_aromatic),
        partner: left
            .effect
            .partner
            .clone()
            .or_else(|| right.effect.partner.clone()),
        when: None,
    }
}

/// One ResonancePair endpoint hit before composition into a pair site.
///
/// Pair rules discover these; [`compose_pair_sites`] joins compatible ends into
/// ordinary [`DeferredSite`] candidates.
#[derive(Clone)]
pub(crate) struct EndCandidate {
    mol: Rc<ForestMol>,
    pub pattern: PatternInfo,
    pub mapped: BTreeMap<u16, usize>,
    /// Conjugated-system anchor (SMARTS map 1).
    pub anchor: usize,
    /// Primary site atom for this endpoint pattern.
    pub site: usize,
}

impl std::fmt::Debug for EndCandidate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EndCandidate")
            .field("pattern", &self.pattern.name)
            .field("anchor", &self.anchor)
            .field("site", &self.site)
            .finish()
    }
}

impl EndCandidate {
    #[allow(dead_code)]
    pub fn mol(&self) -> &Molecule {
        self.mol.mol()
    }

    pub fn forest_rc(&self) -> &Rc<ForestMol> {
        &self.mol
    }
}

/// Match each [`Edit::PairEndpoint`] pattern into endpoint hits (no pairing yet).
pub(crate) fn end_candidates(
    mol: Rc<ForestMol>,
    endpoints: &[PatternInfo],
) -> Result<Vec<EndCandidate>, ForestError> {
    let chemistry = mol.mol();
    let mut out = Vec::new();
    for pattern in endpoints.iter().filter(|p| matches!(p.edit, Edit::PairEndpoint(_))) {
        for mapped in smarts_matches(chemistry, &pattern.smarts)? {
            let Some(&anchor) = mapped.get(&1) else {
                continue;
            };
            let Some(site) = site_atom(&mapped, pattern) else {
                continue;
            };
            out.push(EndCandidate {
                mol: Rc::clone(&mol),
                pattern: pattern.clone(),
                mapped,
                anchor,
                site,
            });
        }
    }
    Ok(out)
}

/// Compose endpoint hits into ResonancePair [`DeferredSite`] candidates.
///
/// Joins two ends that share a conjugated system (or the whole mol fallback),
/// dedupes by site-pair orbit + pattern names. Standard composition used by
/// [`crate::stream::Candidates`].
pub(crate) fn compose_pair_sites(
    ends: &[EndCandidate],
) -> Result<Vec<crate::candidate::DeferredSite>, ForestError> {
    use crate::candidate::DeferredSite;

    if ends.len() < 2 {
        return Ok(Vec::new());
    }
    let forest = Rc::clone(ends[0].forest_rc());
    let mol = forest.mol();

    let mut by_anchor: HashMap<usize, Vec<&EndCandidate>> = HashMap::new();
    for end in ends {
        by_anchor.entry(end.anchor).or_default().push(end);
    }
    if by_anchor.len() < 2 {
        return Ok(Vec::new());
    }

    let chain_conjugate = ends.iter().any(|e| e.pattern.chain_conjugate);

    let mut systems: Vec<HashSet<usize>> = Vec::new();
    let mut covered = HashSet::new();
    for &anchor in by_anchor.keys() {
        if !covered.insert(anchor) {
            continue;
        }
        let (atoms, _) = conjugated_component_ext(mol, anchor, chain_conjugate);
        let set: HashSet<usize> = atoms.iter().copied().collect();
        covered.extend(&set);
        if by_anchor.keys().filter(|a| set.contains(a)).count() >= 2 {
            systems.push(set);
        }
    }
    if systems.is_empty() {
        let all: HashSet<usize> = mol.atoms().map(|(i, _)| atom_usize(i)).collect();
        if by_anchor.keys().filter(|a| all.contains(a)).count() >= 2 {
            systems.push(all);
        }
    }

    let mut out = Vec::new();
    let mut seen_sig: BTreeSet<(usize, usize, String, String, usize)> = BTreeSet::new();
    let gens = crate::orbits::atom_bond_generators(mol);
    let n_atoms = mol.atom_count();

    for system in &systems {
        let anchors: Vec<usize> = by_anchor
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
                let Some(hits_a) = by_anchor.get(&start) else {
                    continue;
                };
                let Some(hits_b) = by_anchor.get(&end) else {
                    continue;
                };
                for left in hits_a {
                    for right in hits_b {
                        if left.site == right.site {
                            continue;
                        }
                        let (n1, n2) = if left.pattern.name <= right.pattern.name {
                            (left.pattern.name.clone(), right.pattern.name.clone())
                        } else {
                            (right.pattern.name.clone(), left.pattern.name.clone())
                        };
                        let sa = left.site.min(right.site);
                        let sb = left.site.max(right.site);
                        let orbit =
                            crate::orbits::atom_pair_orbit_id_with_gens(&gens, n_atoms, sa, sb);
                        if !seen_sig.insert((sa, sb, n1.clone(), n2.clone(), orbit)) {
                            continue;
                        }
                        out.push(DeferredSite::pair(
                            Rc::clone(&forest),
                            sa,
                            format!("{n1}+{n2}"),
                            left.pattern.clone(),
                            right.pattern.clone(),
                            merge_effect_fields(&left.pattern, &right.pattern, system_aromatic),
                            left.mapped.clone(),
                            right.mapped.clone(),
                            start,
                            end,
                            system.clone(),
                            Vec::new(),
                        ));
                    }
                }
            }
        }
    }
    Ok(out)
}

/// Discover composed pair sites: [`end_candidates`] then [`compose_pair_sites`].
pub(crate) fn compose_candidates_from_endpoints(
    mol: Rc<ForestMol>,
    endpoints: &[PatternInfo],
) -> Result<Vec<crate::candidate::DeferredSite>, ForestError> {
    let ends = end_candidates(mol, endpoints)?;
    compose_pair_sites(&ends)
}

/// Run ResonancePair metabolize for the given endpoint patterns.
pub fn pair_metabolize(
    mol: &ForestMol,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairEmission>, ForestError> {
    let mut emissions = Vec::new();
    let mut seen_csmi: BTreeSet<BTreeSet<String>> = BTreeSet::new();
    for site in compose_candidates_from_endpoints(Rc::new(mol.copy_mol()), endpoints)? {
        let Some(emission) = site.apply()? else {
            continue;
        };
        // Explicit CSMI downgrade for PairEmission smoke rows.
        let products = emission.product_csmis();
        let key: BTreeSet<String> = products.iter().cloned().collect();
        if !seen_csmi.insert(key) {
            continue;
        }
        emissions.push(PairEmission {
            site: emission.site,
            pattern_name: emission.pattern_name,
            products,
        });
    }
    Ok(emissions)
}

/// Path dehydrogenation of para-hydroquinone (door smoke for pair metabolize).
pub fn dehydrogenate_hydroquinone(mol: &ForestMol) -> Result<Vec<String>, ForestError> {
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
        }
        .sealed(),
        possibilities: Vec::new(),
        skip_same_rings: false,
        chain_conjugate: false,
        cleave_side_group: None,
        search_bias: 0,
    };
    Ok(pair_metabolize(mol, &[phenol_end])?
        .into_iter()
        .flat_map(|e| e.products)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::canon_of;
    use crate::rules::{dehydrogenation, quinone_formation};

    #[test]
    fn hydroquinone_yields_benzoquinone() {
        let mol = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let products = dehydrogenate_hydroquinone(&mol).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            products.iter().any(|s| canon_of(s).unwrap() == want),
            "want {want}, got {products:?}"
        );
    }

    #[test]
    fn dehydrogenation_rule_pair_door_on_hydroquinone() {
        let mol = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
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
        let mol = ForestMol::parse("c1ccccc1").unwrap();
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
    fn composed_pair_sites_defer_path_flip() {
        let mol = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = compose_candidates_from_endpoints(Rc::new(mol.copy_mol()), &endpoints).unwrap();
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
        let cands = compose_candidates_from_endpoints(Rc::new(mol.copy_mol()), &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(cands.iter().any(|c| {
            c.apply()
                .ok()
                .flatten()
                .map(|e| e.product_csmis())
                .into_iter()
                .flatten()
                .any(|p| canon_of(&p).unwrap() == want)
        }));
    }

    #[test]
    fn tautomerization_cyclohexanone_enol() {
        use crate::rules::tautomerization;
        let mol = ForestMol::parse("O=C1CCCCC1").unwrap();
        let endpoints: Vec<_> = tautomerization()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("OC1=CCCCC1").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| e.products.iter().any(|p| canon_of(p).unwrap() == want)),
            "want {want}; got {emissions:?}"
        );
    }

    #[test]
    fn tautomerization_enol_to_ketone() {
        use crate::rules::tautomerization;
        let mol = ForestMol::parse("OC1=CCCCC1").unwrap();
        let endpoints: Vec<_> = tautomerization()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1CCCCC1").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| e.products.iter().any(|p| canon_of(p).unwrap() == want)),
            "want {want}; got {emissions:?}"
        );
    }

    #[test]
    fn tautomerization_tacrine_amine_imine() {
        use crate::rules::tautomerization;
        let amine = "Nc1c2c(nc3ccccc13)CCCC2";
        let imine = "N=c1c2c([nH]c3ccccc13)CCCC2";
        let mol = ForestMol::parse(amine).unwrap();
        let endpoints: Vec<_> = tautomerization()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of(imine).unwrap();
        let products: Vec<_> = emissions
            .iter()
            .flat_map(|e| e.products.iter().cloned())
            .collect();
        assert!(
            products.iter().any(|p| canon_of(p).unwrap() == want),
            "tacrine amine→imine; got {products:?}"
        );
    }

    #[test]
    fn tautomerization_tacrine_imine_amine() {
        use crate::rules::tautomerization;
        let amine = "Nc1c2c(nc3ccccc13)CCCC2";
        let imine = "N=c1c2c([nH]c3ccccc13)CCCC2";
        let mol = ForestMol::parse(imine).unwrap();
        let endpoints: Vec<_> = tautomerization()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of(amine).unwrap();
        let products: Vec<_> = emissions
            .iter()
            .flat_map(|e| e.products.iter().cloned())
            .collect();
        assert!(
            products.iter().any(|p| canon_of(p).unwrap() == want),
            "tacrine imine→amine; got {products:?}"
        );
    }

    #[test]
    fn tautomerization_tacrine_imine_manual_path_flip() {
        // Legacy-winning path [0,1,2,3,4] on a D-S-D-S kekulé form must yield aromatic amine.
        let imine = ForestMol::parse("N=c1c2c([nH]c3ccccc13)CCCC2").unwrap();
        let want = canon_of("Nc1c2c(nc3ccccc13)CCCC2").unwrap();
        let forms = kekule_forms(imine.mol()).unwrap();
        let path = [0usize, 1, 2, 3, 4];
        let mut any = false;
        for form in &forms {
            let mut rw = form.clone();
            clear_aromatic(&mut rw);
            assert!(flip_path(&mut rw, &path), "flip must touch bonds");
            for &i in &path {
                rw.set_hydrogen_count(atom_idx(i), None);
            }
            if !at_most_one_double_per_atom(&rw) {
                continue;
            }
            if !accept_product(&rw) {
                continue;
            }
            let checked = aromatize(&rw);
            let smiles = canon_smiles(&checked);
            eprintln!("smiles={smiles}");
            if smiles == want {
                any = true;
            }
        }
        assert!(any, "manual legacy path flip should aromatize to amine");
    }

    #[test]
    fn tautomerization_long_range_polyene() {
        use crate::rules::tautomerization;
        let mol = ForestMol::parse("ClCC=CC=CC=CC=CO").unwrap();
        let endpoints: Vec<_> = tautomerization()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("ClCCC=CC=CC=CC=O").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| e.products.iter().any(|p| canon_of(p).unwrap() == want)),
            "want {want}; got {emissions:?}"
        );
    }


    #[test]
    fn tautomerization_find_path_tacrine_7oh() {
        use crate::find_path::{FindPathConfig, PathCounters, find_path_with};
        use crate::rules::default_ruleset;
        let r = "N=c1c2c([nH]c3ccccc13)CCCC2";
        let p = "Nc1c2c(nc3c(O)cccc13)CCCC2";
        let set = default_ruleset();
        let mut c = PathCounters::default();
        let cfg = FindPathConfig {
            max_paths: 1,
            max_nodes: 400,
            use_atom_diff: true,
            ..FindPathConfig::default()
        };
        let hits = find_path_with(r, p, &set, &mut c, cfg, |_| true)
            .unwrap()
            .collect_all()
            .unwrap();
        assert_eq!(
            c.formula_delta_mismatch, 0,
            "tautomer multi-product alternatives must not sum as cleavage; {:?}",
            c.formula_delta_mismatches
        );
        assert!(
            !hits.is_empty(),
            "tacrine→7-OH should hit under default_ruleset (with Tautomerization); bill={} nodes={}",
            c.billed(),
            c.nodes
        );
        let step_names: Vec<_> = hits[0]
            .steps
            .iter()
            .map(|s| s.pattern_name.as_str())
            .collect();
        assert!(
            step_names.len() <= 2,
            "expect tautomer↔OH (≤2 steps), not quinoid+DH detour; bill={} steps={step_names:?}",
            c.billed()
        );
        assert!(
            !step_names.iter().any(|n| *n == "amine"),
            "DH amine must not repair stuck [nH]; bill={} steps={step_names:?}",
            c.billed()
        );
        eprintln!(
            "tacrine→7-OH bill={} steps={step_names:?}",
            c.billed()
        );
    }
}
