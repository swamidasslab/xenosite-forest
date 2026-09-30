//! ResonancePair edits as π-matching constraints (no path flip on ordinary pairs).
//!
//! Match each [`Edit::PairEndpoint`] SMARTS, apply end edits, perceive forced
//! doubles / saturate sites, complete the residual Kekulé matching, then mark
//! surviving aromatic atoms via the cyclic 2-core test (HEURISTICS). Methide
//! is an effect field — two methide ends are allowed.
//!
//! Tautomerization endpoints still use alternating-path flip (extend × far only).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::rc::Rc;

use chematic::core::{Atom, BondOrder, Element};
use chematic::perception::find_sssr;

use crate::kekule::{
    KekuleConfig, KekuleConstraints, PiGraph, bond_order_sums, conjugated_component_ext,
    kekule_forms, move_charge_with_bonds,
};
use crate::forest::{formula_delta, molecule_formula};
use crate::chematic_tags::preserving;
use crate::mol::{ForestError, Molecule, atom_idx, atom_usize, canon_smiles};
use crate::pattern::{Edit, PatternInfo};
use crate::smarts::smarts_matches;
use crate::valence::{accept_pair_product, accept_product, fill_closed_shell_h};
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
            *mol = preserving::with_atom_aromatic(mol, atom, false);
        }
    }
}

/// Demote aromatic bonds touching `system` to single so matching can assign.
fn demote_aromatic_bonds(mol: &mut Molecule, system: &HashSet<usize>) {
    let idxs: Vec<_> = mol
        .bonds()
        .filter_map(|(idx, bond)| {
            if bond.order != BondOrder::Aromatic {
                return None;
            }
            let a = atom_usize(bond.atom1);
            let b = atom_usize(bond.atom2);
            if system.contains(&a) || system.contains(&b) {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    for idx in idxs {
        mol.set_bond_order(idx, BondOrder::Single);
    }
}

/// `keep` + adds H → saturate site atoms (leave must-match, +H later).
fn saturate_sites(pattern: &PatternInfo, mapped: &BTreeMap<u16, usize>) -> BTreeSet<usize> {
    let Edit::PairEndpoint(edit) = &pattern.edit else {
        return BTreeSet::new();
    };
    if edit.as_str() != "keep" {
        return BTreeSet::new();
    }
    let adds = pattern.effect.adds.as_deref().unwrap_or("");
    if !adds.contains('H') {
        return BTreeSet::new();
    }
    mapped.get(&1).copied().into_iter().collect()
}

/// Residual π graph: conjugated edges among system atoms minus saturate sites.
fn residual_pi_graph(
    parent: &Molecule,
    system: &HashSet<usize>,
    saturate: &BTreeSet<usize>,
) -> PiGraph {
    let keep: BTreeSet<usize> = system
        .iter()
        .copied()
        .filter(|a| !saturate.contains(a))
        .collect();
    if let Some(&seed) = keep.iter().next().or_else(|| system.iter().next()) {
        let conj = PiGraph::conjugated(parent, seed);
        let residual = conj.restrict_atoms(&keep);
        if !residual.bonds.is_empty() || residual.atoms.len() <= 1 {
            return residual;
        }
    }
    let mut bonds = BTreeSet::new();
    for (_, bond) in parent.bonds() {
        let a = atom_usize(bond.atom1);
        let b = atom_usize(bond.atom2);
        if keep.contains(&a) && keep.contains(&b) {
            match bond.order {
                BondOrder::Single
                | BondOrder::Double
                | BondOrder::Triple
                | BondOrder::Aromatic => {
                    bonds.insert(if a < b { (a, b) } else { (b, a) });
                }
                _ => {}
            }
        }
    }
    PiGraph::new(keep, bonds)
}

/// Non-site mapped atoms when [`Effect::exclusive_partner`] is set.
fn exclusive_partner_atoms(
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
) -> BTreeSet<usize> {
    if !pattern.effect.exclusive_partner {
        return BTreeSet::new();
    }
    let Some(site) = site_atom(mapped, pattern) else {
        return BTreeSet::new();
    };
    mapped
        .values()
        .copied()
        .filter(|&idx| idx != site)
        .collect()
}

/// True when an exclusive partner atom appears on the other end's map.
fn shared_exclusive_partner(
    map1: &BTreeMap<u16, usize>,
    info1: &PatternInfo,
    map2: &BTreeMap<u16, usize>,
    info2: &PatternInfo,
) -> bool {
    let exclusive1 = exclusive_partner_atoms(map1, info1);
    let exclusive2 = exclusive_partner_atoms(map2, info2);
    if exclusive1.is_empty() && exclusive2.is_empty() {
        return false;
    }
    let other1: BTreeSet<usize> = map2.values().copied().collect();
    let other2: BTreeSet<usize> = map1.values().copied().collect();
    !exclusive1.is_disjoint(&other1) || !exclusive2.is_disjoint(&other2)
}

/// Exocyclic hetero leaf double (quinone/amide C=O). Ring C=N rematch freely.
fn is_exocyclic_hetero_leaf(mol: &Molecule, a: usize, b: usize) -> bool {
    let z_a = mol.atom(atom_idx(a)).element.atomic_number();
    let z_b = mol.atom(atom_idx(b)).element.atomic_number();
    if z_a == 6 && z_b == 6 {
        return false;
    }
    let heavy_deg = |i: usize| mol.neighbors(atom_idx(i)).count();
    heavy_deg(a) == 1 || heavy_deg(b) == 1
}

fn charge_baseline_atoms(parent: &Molecule, edited: &Molecule) -> BTreeSet<usize> {
    crate::valence::edited_valence_atoms(parent, edited)
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

/// Re-apply +1 on iminium N (map 2) after bond-order / charge travel.
fn reapply_iminium_charges(
    mol: &mut Molecule,
    left: &PatternInfo,
    right: &PatternInfo,
    map1: &BTreeMap<u16, usize>,
    map2: &BTreeMap<u16, usize>,
) {
    for (pat, map) in [(left, map1), (right, map2)] {
        if matches!(&pat.edit, Edit::PairEndpoint(e) if e == "iminium") {
            if let Some(&n) = map.get(&2) {
                mol.set_charge(atom_idx(n), 1);
            }
        }
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
            *mol = preserving::with_atom_charge(mol, atom_idx(n), 1);
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
    *mol = preserving::with_atom_element(mol, atom_idx(halogen), Element::O);
    *mol = preserving::with_atom_charge(mol, atom_idx(halogen), 0);
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

/// C18 constraint materialize (no path flip). Used for ordinary ResonancePair
/// ends; tautomer pairs keep the path-flip door in [`materialize_pair_mols`].
fn materialize_pair_constraints(
    mol: &Molecule,
    left: &PatternInfo,
    right: &PatternInfo,
    effect: &crate::pattern::Effect,
    map1: &BTreeMap<u16, usize>,
    map2: &BTreeMap<u16, usize>,
    system: &HashSet<usize>,
) -> Result<Vec<Molecule>, ForestError> {
        let rings = ring_sets(mol);
        let mut rw = mol.clone();
        clear_aromatic(&mut rw);
        demote_aromatic_bonds(&mut rw, system);
        if !edit_end(&mut rw, map1, left, &rings) {
            return Ok(Vec::new());
        }
        if !edit_end(&mut rw, map2, right, &rings) {
            return Ok(Vec::new());
        }

        let saturate = saturate_sites(left, map1)
            .union(&saturate_sites(right, map2))
            .copied()
            .collect::<BTreeSet<_>>();
        // Emit-path closed-shell settle + charge baseline (perception of
        // neighbor / bond-sum change — not edit-token names). Kekulé matcher
        // stays generic; config is the behavior contract.
        let mut settle = charge_baseline_atoms(mol, &rw);
        for &a in &settle {
            fill_closed_shell_h(&mut rw, a);
        }
        let residual = residual_pi_graph(mol, system, &saturate);
        // Parent doubles + cumulated degree (N=C=O central C): needed before
        // saturate→residual demote so framework edges are not consumed as
        // styrene-style residual π (else both cumulated doubles fall and H on
        // PhNCO emits empty / [S-] junk instead of O=CNAr).
        let mut parent_double_deg: HashMap<usize, usize> = HashMap::new();
        let parent_doubles: BTreeSet<(usize, usize)> = mol
            .bonds()
            .filter_map(|(_, bond)| {
                (bond.order == BondOrder::Double).then(|| {
                    let a = atom_usize(bond.atom1);
                    let b = atom_usize(bond.atom2);
                    *parent_double_deg.entry(a).or_default() += 1;
                    *parent_double_deg.entry(b).or_default() += 1;
                    bond_key(a, b)
                })
            })
            .collect();
        let is_cumulated_framework = |a: usize, b: usize| -> bool {
            let edge = bond_key(a, b);
            parent_doubles.contains(&edge)
                && (parent_double_deg.get(&a).copied().unwrap_or(0) >= 2
                    || parent_double_deg.get(&b).copied().unwrap_or(0) >= 2)
        };
        // Saturate sites consume incident π bonds before residual rematch:
        // (1) shared saturate–saturate edge (aldehyde C=O, ethene, amide);
        // (2) saturate→residual edges (styrene vinyl C=C when only the CH2
        //     end is a path_end — else rematch keeps C=C and mints allenes).
        //     Skip cumulated framework edges here: shared-edge demote (1) is
        //     the only way to saturate one half of N=C=O / N=C=S / N=C=N.
        let mut shared_edge_saturated = false;
        let demote_pi = |rw: &mut Molecule, a: usize, b: usize| -> bool {
            let Some((bond_idx, bond)) = rw.bond_between(atom_idx(a), atom_idx(b)) else {
                return false;
            };
            if !matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
            ) {
                return false;
            }
            rw.set_bond_order(
                bond_idx,
                if bond.order == BondOrder::Triple {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            );
            true
        };
        if saturate.len() == 2 {
            let mut ends = saturate.iter().copied();
            let a = ends.next().unwrap();
            let b = ends.next().unwrap();
            if demote_pi(&mut rw, a, b) {
                shared_edge_saturated = true;
            }
        }
        let mut residual_pi_saturated = false;
        for &s in &saturate {
            let nbrs: Vec<usize> = rw
                .neighbors(atom_idx(s))
                .map(|(n, _)| atom_usize(n))
                .collect();
            for n in nbrs {
                if !residual.atoms.contains(&n) {
                    continue;
                }
                if is_cumulated_framework(s, n) {
                    continue;
                }
                if demote_pi(&mut rw, s, n) {
                    residual_pi_saturated = true;
                }
            }
        }
        if shared_edge_saturated || residual_pi_saturated {
            for &atom in &saturate {
                fill_closed_shell_h(&mut rw, atom);
            }
            settle.extend(saturate.iter().copied());
        }
        let all_forced = residual.perceive_forced_doubles(&rw);
        // Cumulated parent doubles (N=C=O: central C has two doubles) cannot
        // exclusive-seed under one-partner matching — drop those atoms from the
        // residual (fixed framework). Hetero parent doubles (other quinone C=O)
        // exclusive-seed; pure C=C rematch freely after blanking residual π.
        // Edit-new doubles (add_carbonyl O, phenol C=O) remain forced seeds.
        let mut edit_forced = BTreeSet::new();
        let mut framework_forced = BTreeSet::new();
        for &edge in &all_forced {
            let (a, b) = edge;
            if is_cumulated_framework(a, b) {
                // N=C=O: strip from residual (one-partner matching cannot
                // express cumulated demand).
                framework_forced.insert(edge);
            } else if !parent_doubles.contains(&edge) {
                // Edit-new leaves (add_carbonyl O, phenol C=O): exclusive-seed.
                edit_forced.insert(edge);
            } else {
                // Surviving parent double: exclusive-seed only exocyclic hetero
                // leaves (degree-1 O/S/N on the edge — quinone C=O, amide).
                // Ring-embedded C=N / N=N rematch freely like C=C (HEURISTICS).
                if is_exocyclic_hetero_leaf(&rw, a, b) {
                    edit_forced.insert(edge);
                }
            }
        }
        let residual_match = residual.after_forced_doubles(&framework_forced);
        // Blank surviving parent doubles that are not exclusive-seeded (ring
        // C=C / C=N / N=N) so atom_must_be_matched does not lock them. Skip
        // edit-forced / exocyclic leaves; forced constraints re-assert them.
        for &(a, b) in &residual_match.bonds {
            let edge = bond_key(a, b);
            if edit_forced.contains(&edge) {
                continue;
            }
            if !parent_doubles.contains(&edge) {
                continue;
            }
            let _ = demote_pi(&mut rw, a, b);
        }
        // Empty residual: one-edge path_end (shared π saturated) may emit.
        // Vacuous keep+keep, or two carbonyl carbons that only demoted
        // leaf C=O into an empty residual (glyoxal → glycol), refuse.
        if edit_forced.is_empty() && residual_match.bonds.is_empty() {
            if shared_edge_saturated {
                let parent_csmi = canon_smiles(mol);
                let product_csmi = canon_smiles(&rw);
                if product_csmi != parent_csmi {
                    let checked = preserving::aromatize(&rw);
                    let mut products = Vec::new();
                    let mut local_csmi = BTreeSet::new();
                    for frag in checked.fragments() {
                        if !accept_pair_product(&frag, left, right) {
                            continue;
                        }
                        let smiles = canon_smiles(&frag);
                        if smiles == parent_csmi {
                            continue;
                        }
                        if local_csmi.insert(smiles) {
                            products.push(frag);
                        }
                    }
                    return Ok(products);
                }
            }
            return Ok(Vec::new());
        }
        let config = KekuleConfig::for_constraints();
        let constraints = KekuleConstraints::new()
            .with_forced(edit_forced.clone())
            .with_saturate(saturate.clone());
        // Caller already perceived forced doubles — do not double-count from mol.
        let match_cfg = config.explicit_forced_only();
        let assignments = residual_match.all_assignments(&rw, &constraints, &match_cfg);
        if assignments.is_empty() {
            return Ok(Vec::new());
        }

        // Residual aromaticity: cyclic 2-core after dropping demand-consumed
        // atoms (forced doubles / saturate) — HEURISTICS, not sanitize flags.
        let arom = residual.after_forced_doubles(&all_forced);
        let aromatic_core = arom.aromatic_2core_atoms(&rw, &config);
        if effect.dearomatizes {
            let edited: HashSet<usize> = system
                .iter()
                .copied()
                .filter(|i| mol.atom(atom_idx(*i)).aromatic)
                .collect();
            if !edited.is_empty() && edited.iter().all(|i| aromatic_core.contains(i)) {
                // Capability claimed dearomatization; whole edited system still aromatic.
                return Ok(Vec::new());
            }
        }

        let mut products = Vec::new();
        let mut local_csmi = BTreeSet::new();
        let parent_csmi = canon_smiles(mol);
        // Charge/H follow parent π bond sums (aromatic = 1.5). Edited-valence
        // atoms use the post-edit baseline so leave fragments stay closed-shell
        // (CH4 not [CH5]; phenol O not [OH+]) — no sanitize rescue.
        let mut before = bond_order_sums(mol);
        let post_edit = bond_order_sums(&rw);
        for a in settle {
            if let Some(&v) = post_edit.get(&a) {
                before.insert(a, v);
            }
        }
        let was_aromatic: HashSet<usize> = mol
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(atom_usize(idx)))
            .collect();
        for assignment in assignments {
            let mut product = rw.clone();
            for (&(left, right), &order) in &assignment {
                if let Some((bond_idx, _)) =
                    product.bond_between(atom_idx(left), atom_idx(right))
                {
                    product.set_bond_order(bond_idx, order);
                }
            }
            for &atom in &residual.atoms {
                if product.atom(atom_idx(atom)).aromatic {
                    product = preserving::with_atom_aromatic(&product, atom_idx(atom), false);
                }
            }
            // Saturate sites + residual π atoms: H from final bond orders.
            // Early settle fill (pre-rematch, all aromatic→single) can leave H
            // on a hetero that rematch then doubles — pyridine para H minted
            // `C1C=CC[NH]=C1` instead of neutral `C1=CCN=CC1`.
            // Saturate sites: H from final bond orders (closed shell), not a
            // blind +1 on top of the pre-match settle — that minted [CH3] on
            // aromatic path_end hydrogenation (benzene → cyclohexadiene).
            for &atom in &saturate {
                fill_closed_shell_h(&mut product, atom);
            }
            // Residual heteros only: early settle fill (pre-rematch, aromatic→
            // single) can leave H on N/O that rematch then doubles — pyridine
            // para H minted `C1C=CC[NH]=C1` instead of neutral `C1=CCN=CC1`.
            // Carbons stay on move_charge H-travel (aromatic bond-sum delta).
            for &atom in &residual.atoms {
                let z = product.atom(atom_idx(atom)).element.atomic_number();
                if z == 7 || z == 8 {
                    fill_closed_shell_h(&mut product, atom);
                }
            }
            move_charge_with_bonds(&mut product, &before, &was_aromatic);
            reapply_iminium_charges(&mut product, left, right, map1, map2);
            // Stamp surviving aromatic 2-core; leave the rest localized.
            for &atom in &aromatic_core {
                product = preserving::with_atom_aromatic(&product, atom_idx(atom), true);
            }
            // Perception finish (RDKit-parity aromaticity). Dearomatize refuse
            // already used the 2-core gate above — not sanitize flags.
            let checked = preserving::aromatize(&product);
            for frag in checked.fragments() {
                if !accept_pair_product(&frag, left, right) {
                    continue;
                }
                let smiles = canon_smiles(&frag);
                // keep+keep / path rematch must change the molecule — identity
                // is vacuous (esters, pyrrole, ethene) and not a metabolite.
                if smiles == parent_csmi {
                    continue;
                }
                if local_csmi.insert(smiles) {
                    products.push(frag);
                }
            }
        }
        Ok(products)
}

/// Native ResonancePairRule: odd bond-count alternating paths, end edits, flip.
fn materialize_pair_path_flip(
    mol: &Molecule,
    left: &PatternInfo,
    right: &PatternInfo,
    effect: &crate::pattern::Effect,
    map1: &BTreeMap<u16, usize>,
    map2: &BTreeMap<u16, usize>,
    start: usize,
    end: usize,
    system: &HashSet<usize>,
    tautomer_anchor_h: Option<(Option<usize>, usize)>,
) -> Result<Vec<Molecule>, ForestError> {
    let is_tautomer = tautomer_anchor_h.is_some();
    let forms = kekule_forms(mol)?;
    let rings = ring_sets(mol);
    let neighbors = system_neighbors(mol, system);
    let mut products = Vec::new();
    let mut local_csmi = BTreeSet::new();
    let parent_csmi = canon_smiles(mol);
    for form in &forms {
        let bond_map = current_orders(form);
        let mut paths = alternating_from(&bond_map, start, end, &neighbors, 2);
        paths.extend(alternating_from(&bond_map, end, start, &neighbors, 2));
        for path in paths {
            // Ordinary ResonancePair: odd bond count (native gate).
            if !is_tautomer && (path.len() - 1) % 2 != 1 {
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
            reapply_iminium_charges(&mut rw, left, right, map1, map2);
            if is_tautomer {
                for &i in &flip_atoms {
                    rw.set_hydrogen_count(atom_idx(i), None);
                    fill_closed_shell_h(&mut rw, i);
                }
            }
            if is_tautomer && !at_most_one_double_per_atom(&rw) {
                continue;
            }
            if !accept_pair_product(&rw, left, right) {
                continue;
            }
            let checked = preserving::aromatize(&rw);
            if effect.dearomatizes && system_stayed_aromatic(mol, &checked, system) {
                continue;
            }
            // Tautomerization is formula-neutral (empty Effect bags). Path-flip
            // on aromatic CH can mint cyclohexadiene (+2H) — refuse those; keep
            // enol↔ketone / phenol↔quinone-methide that conserve the formula.
            let parent_f = if is_tautomer {
                Some(molecule_formula(mol))
            } else {
                None
            };
            for frag in checked.fragments() {
                if !accept_pair_product(&frag, left, right) {
                    continue;
                }
                if let Some(ref pf) = parent_f {
                    let delta = formula_delta(pf, &molecule_formula(&frag));
                    if delta.counts.values().any(|&n| n != 0) {
                        continue;
                    }
                }
                let smiles = canon_smiles(&frag);
                if smiles == parent_csmi {
                    continue;
                }
                if local_csmi.insert(smiles) {
                    products.push(frag);
                }
            }
        }
    }
    Ok(products)
}

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

    if is_tautomer {
        return materialize_pair_path_flip(
            mol,
            left,
            right,
            effect,
            map1,
            map2,
            start,
            end,
            system,
            tautomer_anchor_h,
        );
    }

    materialize_pair_constraints(mol, left, right, effect, map1, map2, system)
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
        breaks_ring: left.effect.breaks_ring || right.effect.breaks_ring,
        leave_count: left.effect.leave_count.or(right.effect.leave_count),
        methide: left.effect.methide || right.effect.methide,
        exclusive_partner: left.effect.exclusive_partner || right.effect.exclusive_partner,
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
/// Joins two ends that share a conjugated system; dedupes by site pair, pattern
/// names, and path-end anchors (C18 / Python `pair_site_signature`).
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
    // Native ResonancePairRule: no conjugated system → no endpoint pairs (no
    // whole-molecule fallback).

    let mut out = Vec::new();
    let mut seen_sig: BTreeSet<(usize, usize, String, String, usize, usize)> = BTreeSet::new();

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
                        if shared_exclusive_partner(
                            &left.mapped,
                            &left.pattern,
                            &right.mapped,
                            &right.pattern,
                        ) {
                            continue;
                        }
                        let left_res = left.pattern.resolve_for_match(mol, &left.mapped);
                        let right_res = right.pattern.resolve_for_match(mol, &right.mapped);
                        let (n1, n2) = if left_res.name <= right_res.name {
                            (left_res.name.clone(), right_res.name.clone())
                        } else {
                            (right_res.name.clone(), left_res.name.clone())
                        };
                        let sa = left.site.min(right.site);
                        let sb = left.site.max(right.site);
                        let pe_a = start.min(end);
                        let pe_b = start.max(end);
                        if !seen_sig.insert((sa, sb, n1.clone(), n2.clone(), pe_a, pe_b)) {
                            continue;
                        }
                        out.push(DeferredSite::pair(
                            Rc::clone(&forest),
                            sa,
                            format!("{n1}+{n2}"),
                            left_res.clone(),
                            right_res.clone(),
                            merge_effect_fields(&left_res, &right_res, system_aromatic),
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
    use crate::candidate::DeferredSite;
    use crate::forest_mol::ForestMol;
    use crate::mol::{atom_usize, canon_of, canon_smiles, parse_mol};
    use crate::rules::{dehydrogenation, quinone_formation};

    fn pair_candidates(fm: &ForestMol, endpoints: &[PatternInfo]) -> Vec<DeferredSite> {
        compose_candidates_from_endpoints(Rc::new(fm.copy_mol()), endpoints).unwrap()
    }

    fn materialize_csmis(c: &DeferredSite) -> Vec<String> {
        c.materialize_mols()
            .unwrap_or_default()
            .into_iter()
            .map(|m| canon_of(&canon_smiles(&m)).unwrap())
            .collect()
    }

    fn qf_pair_endpoints() -> Vec<PatternInfo> {
        quinone_formation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect()
    }

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
    fn quinone_formation_add_carbonyl_on_phnco_keeps_nco() {
        let mol = ForestMol::parse("O=C=Nc1ccccc1").unwrap();
        let endpoints = qf_pair_endpoints();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = [
            canon_of("O=C=NC1=CC(=O)C=CC1=O").unwrap(),
            canon_of("O=C=NC1=CC=CC(=O)C1=O").unwrap(),
            canon_of("O=C=NC1=CC(=O)C(=O)C=C1").unwrap(),
        ];
        for w in &want {
            assert!(
                emissions
                    .iter()
                    .any(|e| e.products.iter().any(|p| canon_of(p).unwrap() == *w)),
                "missing {w}; got {emissions:?}"
            );
        }
    }

    #[test]
    fn quinone_formation_pyridine_para_emits_pyridinedione() {
        let mol = ForestMol::parse("c1ccncc1").unwrap();
        let endpoints = qf_pair_endpoints();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let got: BTreeSet<_> = emissions
            .iter()
            .flat_map(|e| e.products.iter().map(|p| canon_of(p).unwrap()))
            .collect();
        let want = [
            canon_of("O=C1C=CC(=O)N=C1").unwrap(),
            canon_of("O=C1C=CC=NC1=O").unwrap(),
            canon_of("O=C1C=CN=CC1=O").unwrap(),
        ];
        assert!(
            want.iter().any(|w| got.contains(w)),
            "pyridine QF should emit a pyridinedione; got {got:?}"
        );
    }

    #[test]
    fn pair_candidates_defer_materialize() {
        let mol = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints);
        assert!(!cands.is_empty());
        // Phenol site_map is O (not aromatic): resolve clears dearomatizes before
        // merge — same as Python resolve_effect + merge_effects. Materialize still
        // emits the quinone from the π constraints.
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(cands.iter().any(|c| {
            materialize_csmis(c).iter().any(|p| p == &want)
        }));
    }

    #[test]
    fn quinone_dealkylate_splits_fragments_like_find_path() {
        use crate::product_graph::{ProductGraphConfig, product_layer};
        let set = quinone_formation();
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let layer = ProductGraphConfig {
            target: None,
            max_nodes: usize::MAX,
            max_depth: usize::MAX,
        };
        let children = product_layer(&parent, &set, &layer).unwrap();
        let want_q = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        let want_me = canon_of("C").unwrap();
        let csmi: Vec<String> = children
            .iter()
            .map(|c| c.child.csmi().as_ref().to_string())
            .collect();
        assert!(
            csmi.iter().any(|s| canon_of(s).unwrap() == want_q),
            "quinone fragment missing: {csmi:?}"
        );
        assert!(
            csmi.iter().any(|s| canon_of(s).unwrap() == want_me),
            "methyl fragment missing: {csmi:?}"
        );
        assert!(
            csmi.iter().all(|s| !s.contains('.')),
            "disconnected CSMI should be split before yield: {csmi:?}"
        );
        let with_both = children.iter().any(|c| {
            c.hop.rule == "QuinoneFormation"
                && c.hop.cleaves
                && c.hop.products.len() >= 2
                && c.hop.products.iter().any(|p| canon_of(p).unwrap() == want_q)
                && c.hop.products.iter().any(|p| canon_of(p).unwrap() == want_me)
        });
        assert!(
            with_both,
            "cleaving hop must name QuinoneFormation and list both fragments"
        );
    }

    #[test]
    fn quinone_formation_dealkylates_phenyl_ncx_ahead_of_python() {
        let endpoints = qf_pair_endpoints();
        for (smi, leave, imines) in [
            ("O=C=Nc1ccccc1", "C=O", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
            ("S=C=Nc1ccccc1", "C=S", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
            ("N=C=Nc1ccccc1", "C=N", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
        ] {
            let mol = ForestMol::parse(smi).unwrap();
            let emissions = pair_metabolize(&mol, &endpoints).unwrap();
            let products: BTreeSet<String> = emissions
                .iter()
                .flat_map(|e| e.products.iter().cloned())
                .map(|s| canon_of(&s).unwrap())
                .collect();
            let want_leave = canon_of(leave).unwrap();
            assert!(products.contains(&want_leave), "{smi}: missing leave {leave}");
            for imine in imines {
                let want = canon_of(imine).unwrap();
                assert!(products.contains(&want), "{smi}: missing {imine}");
            }
        }
    }

    #[test]
    fn quinone_formation_declares_exclusive_partner_on_hetero_arms() {
        let set = quinone_formation();
        let by_name: BTreeMap<_, _> = set
            .patterns()
            .into_iter()
            .map(|p| (p.name.as_str(), p))
            .collect();
        for name in [
            "single_to_double",
            "replace_halogen",
            "iminium",
            "dealkylate",
        ] {
            let p = by_name.get(name).unwrap_or_else(|| panic!("missing {name}"));
            assert!(p.effect.exclusive_partner, "{name} should set exclusive_partner");
        }
        let methide = by_name.get("methide_end").unwrap();
        assert!(
            !methide.effect.exclusive_partner,
            "methide alkyl partner stays off exclusive_partner"
        );
        assert_eq!(methide.effect.partner.as_deref(), Some("C"));
    }

    #[test]
    fn shared_exclusive_partner_detects_bridging_atom() {
        let std = qf_pair_endpoints()
            .into_iter()
            .find(|p| p.name == "single_to_double")
            .unwrap();
        assert!(std.effect.exclusive_partner);
        let mut map1 = BTreeMap::new();
        map1.insert(1, 0);
        map1.insert(2, 10);
        let mut map2 = BTreeMap::new();
        map2.insert(1, 5);
        map2.insert(2, 10);
        assert!(shared_exclusive_partner(&map1, &std, &map2, &std));
        map2.insert(2, 11);
        assert!(!shared_exclusive_partner(&map1, &std, &map2, &std));
    }

    #[test]
    fn bridging_n_pair_candidates_do_not_share_exclusive_partner() {
        for smiles in [
            "c1ccc(N(C)c2ccccc2)cc1",
            "c1ccc2c(c1)Nc1ccccc1C2",
            "c1ccc2c(c1)Nc1ccccc1O2",
        ] {
            let mol = ForestMol::parse(smiles).unwrap();
            let endpoints = qf_pair_endpoints();
            for c in pair_candidates(&mol, &endpoints) {
                let map2 = c.pair.as_ref().unwrap().map2.clone();
                assert!(
                    !shared_exclusive_partner(&c.mapped, &c.left, &map2, &c.right),
                    "{smiles}: survivor {} still shares exclusive partner",
                    c.pattern_name
                );
            }
        }
    }

    #[test]
    fn catechol_identical_o_partners_still_emit() {
        let mol = ForestMol::parse("Oc1ccccc1O").unwrap();
        let cands = pair_candidates(&mol, &qf_pair_endpoints());
        let saw_phenol_pair = cands.iter().any(|c| {
            c.left.effect.partner.as_deref() == Some("O")
                && c.right.effect.partner.as_deref() == Some("O")
                && c.left.effect.exclusive_partner
                && c.right.effect.exclusive_partner
        });
        let saw_distinct_o = cands.iter().any(|c| {
            c.left.name == "single_to_double"
                && c.right.name == "single_to_double"
                && c.mapped.get(&2) != c.pair.as_ref().unwrap().map2.get(&2)
                && c.mapped.get(&2).is_some()
                && c.pair.as_ref().unwrap().map2.get(&2).is_some()
        });
        assert!(
            saw_phenol_pair || saw_distinct_o,
            "expected ortho catechol pair with distinct O partners; got {:?}",
            cands
                .iter()
                .map(|c| (
                    c.pattern_name.as_str(),
                    c.mapped.clone(),
                    c.pair.as_ref().map(|p| p.map2.clone()),
                    c.left.effect.partner.clone(),
                    c.right.effect.partner.clone(),
                ))
                .collect::<Vec<_>>()
        );
        assert!(!cands.is_empty(), "catechol should still emit pair candidates");
    }

    #[test]
    fn constraint_products_are_closed_shell() {
        let endpoints = qf_pair_endpoints();
        for smi in ["Oc1ccc(O)cc1", "COc1ccccc1", "O=C=Nc1ccccc1"] {
            let mol = ForestMol::parse(smi).unwrap();
            for c in pair_candidates(&mol, &endpoints) {
                for p in materialize_csmis(&c) {
                    assert!(!p.contains("[C]"), "{smi}: radical in {p}");
                    assert!(!p.contains("[CH5]") && !p.contains("[CH6]"), "{smi}: bad methyl {p}");
                    assert!(!p.contains("[OH+]"), "{smi}: bad charge {p}");
                }
            }
        }
    }

    #[test]
    fn apap_amine_phenol_pair_emits_quinone_imine() {
        let mol = ForestMol::parse("CC(=O)Nc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints);
        let want = canon_of("CC(=O)N=C1C=CC(=O)C=C1").unwrap();
        let products: Vec<String> = cands
            .iter()
            .filter(|c| {
                (c.left.name.contains("amine") && c.right.name.contains("phenol"))
                    || (c.left.name.contains("phenol") && c.right.name.contains("amine"))
            })
            .flat_map(materialize_csmis)
            .collect();
        assert!(
            products.iter().any(|p| p == &want),
            "APAP DH pair should emit quinone-imine; got {products:?}"
        );
    }

    /// 1,4-naphthalenediol → 1,4-naphthoquinone: fused system partially
    /// collapses; product keeps one aromatic ring (`c2ccccc12`).
    #[test]
    fn naphthalene_diol_dh_emits_14_naphthoquinone() {
        let mol = ForestMol::parse("Oc1ccc(O)c2ccccc12").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints);
        let want = canon_of("O=C1C=CC(=O)c2ccccc12").unwrap();
        let products: Vec<String> = cands
            .iter()
            .filter(|c| c.left.name.contains("phenol") && c.right.name.contains("phenol"))
            .flat_map(materialize_csmis)
            .collect();
        assert!(
            products.iter().any(|p| p == &want),
            "1,4-naphthalenediol DH should emit naphthoquinone; got {products:?}"
        );
        let product = cands
            .iter()
            .filter(|c| c.left.name.contains("phenol") && c.right.name.contains("phenol"))
            .flat_map(|c| c.materialize_mols().unwrap_or_default())
            .find(|m| canon_of(&canon_smiles(m)).unwrap() == want)
            .expect("naphthoquinone mol");
        let aromatic: Vec<_> = product
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(atom_usize(idx)))
            .collect();
        assert_eq!(
            aromatic.len(),
            6,
            "fused ring remains aromatic after partial collapse: {aromatic:?}"
        );
    }

    #[test]
    fn hydrogenation_benzene_path_end_emits_cyclohexadiene() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("c1ccccc1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints);
        assert!(!cands.is_empty());
        let want_13 = canon_of("C1=CCC=CC1").unwrap();
        let want_14 = canon_of("C1=CCCC=C1").unwrap();
        let csmi: Vec<_> = cands.iter().flat_map(materialize_csmis).collect();
        assert!(csmi.iter().any(|p| p == &want_13) && csmi.iter().any(|p| p == &want_14));
    }

    #[test]
    fn hydrogenation_ethene_path_end_refuses_identity() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("C=C").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let products: Vec<String> = pair_candidates(&mol, &endpoints)
            .iter()
            .flat_map(materialize_csmis)
            .collect();
        let parent = canon_of("C=C").unwrap();
        assert!(products.iter().all(|p| p != &parent));
        assert!(products.iter().any(|p| p == &canon_of("CC").unwrap()));
    }

    #[test]
    fn hydrogenation_acetaldehyde_path_end_emits_ethanol() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("CC=O").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let want = canon_of("CCO").unwrap();
        let products: Vec<String> = pair_candidates(&mol, &endpoints)
            .iter()
            .flat_map(materialize_csmis)
            .collect();
        assert!(products.iter().any(|p| p == &want));
    }

    #[test]
    fn hydrogenation_apap_amide_path_end_emits_hemiaminal() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("CC(=O)Nc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let want = canon_of("CC(O)Nc1ccc(O)cc1").unwrap();
        let products: Vec<String> = pair_candidates(&mol, &endpoints)
            .iter()
            .flat_map(materialize_csmis)
            .collect();
        assert!(products.iter().any(|p| p == &want), "got {products:?}");
    }

    #[test]
    fn hydrogenation_styrene_vinyl_ring_path_emits_exocyclic() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("C=Cc1ccccc1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let want = canon_of("CC=C1C=CC=CC1").unwrap();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        assert!(emissions.iter().any(|e| {
            e.products.iter().any(|p| canon_of(p).unwrap() == want)
        }));
    }

    #[test]
    fn hydrogenation_phnco_saturates_one_cumulated_double() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("O=C=Nc1ccccc1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let got: BTreeSet<_> = pair_metabolize(&mol, &endpoints)
            .unwrap()
            .into_iter()
            .flat_map(|e| e.products)
            .map(|p| canon_of(&p).unwrap())
            .collect();
        let amide = canon_of("O=CNc1ccccc1").unwrap();
        let iminol = canon_of("OC=Nc1ccccc1").unwrap();
        assert!(got.contains(&amide) || got.contains(&iminol), "got {got:?}");
    }

    #[test]
    fn hydrogenation_pyridine_para_emits_neutral_dihydropyridine() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("c1ccncc1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let got: BTreeSet<_> = pair_metabolize(&mol, &endpoints)
            .unwrap()
            .into_iter()
            .flat_map(|e| e.products)
            .map(|p| canon_of(&p).unwrap())
            .collect();
        let want = canon_of("C1=CCN=CC1").unwrap();
        assert!(got.contains(&want), "got {got:?}");
        assert!(!got.iter().any(|p| p.contains("[NH+") || p.contains("[nH+]")));
    }

    #[test]
    fn hydrogenation_benzoquinone_para_o_emits_hydroquinone() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("O=C1C=CC(=O)C=C1").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let want_hq = canon_of("Oc1ccc(O)cc1").unwrap();
        let got: Vec<_> = pair_metabolize(&mol, &endpoints)
            .unwrap()
            .into_iter()
            .flat_map(|e| e.products)
            .map(|p| canon_of(&p).unwrap())
            .collect();
        assert!(got.iter().any(|p| p == &want_hq), "got {got:?}");
    }

    #[test]
    fn hydrogenation_glyoxal_refuses_glycol_from_carbon_pair() {
        use crate::rules::hydrogenation;
        let mol = ForestMol::parse("O=CC=O").unwrap();
        let endpoints: Vec<_> = hydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let glycol = canon_of("OCCO").unwrap();
        assert!(emissions.iter().all(|e| {
            e.products.iter().all(|p| canon_of(p).unwrap() != glycol)
        }));
        let want = canon_of("OC=CO").unwrap();
        assert!(emissions.iter().any(|e| {
            e.products.iter().any(|p| canon_of(p).unwrap() == want)
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
            let checked = preserving::aromatize(&rw);
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
    #[ignore = "xfail: formula_delta_mismatch on tautomer multi-product alts until Effect filters re-enabled; remove ignore when green"]
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
