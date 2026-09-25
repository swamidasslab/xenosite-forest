//! Reactant→target local diff for candidate filtering.
//!
//! Parity door for Python `find_path.atom_diff` / `_site_could_help` /
//! `_pattern_could_help`. Aligns via chematic MCS with `BondCompare::Any`
//! (Python `rdFMCS` CompareAny). Multi-placement views merge so a step that
//! helps any ring is not refused. Filters read [`crate::pattern::Effect`]
//! on deferred candidates — no filter closures required.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use chematic::core::BondOrder;
use chematic::perception::ring_atom_flags;
use chematic::smarts::{BondCompare, McsConfig, find_matches, find_mcs_with_config};

use crate::candidate::Candidate;
use crate::mol::{Molecule, atom_idx, atom_usize, ranks};
use crate::pair_edit::PairCandidate;
use crate::pattern::Effect;

fn hydrogens(mol: &Molecule, idx: usize) -> i32 {
    let atom = mol.atom(atom_idx(idx));
    let explicit = atom.hydrogen_count.unwrap_or(0) as i32;
    let implicit = mol.implicit_hydrogen_count(atom_idx(idx)) as i32;
    explicit + implicit
}

fn order_value(order: BondOrder) -> f32 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Aromatic => 1.5,
        _ => 1.0,
    }
}

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn heavy_atom_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count()
}

fn formula_oxygen(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() == 8)
        .count()
}

/// Local reactant-to-target differences (Python `AtomDiff`).
#[derive(Clone, Debug, Default)]
pub struct AtomDiff {
    pub mapping: BTreeMap<usize, usize>,
    /// Every placement (best first). Filters / cost use the merge.
    pub mappings: Vec<BTreeMap<usize, usize>>,
    pub cleaved: BTreeSet<usize>,
    pub cleavage_bonds: BTreeSet<(usize, usize)>,
    pub loses_aromaticity: BTreeSet<usize>,
    pub h_delta: BTreeMap<usize, i32>,
    pub bond_raises: BTreeSet<(usize, usize)>,
    pub bond_order_mismatches: usize,
    pub n_extra: usize,
    pub reactant_heavy: usize,
    pub target_heavy: usize,
    /// Per-view field costs; [`Self::cost`] is the minimum.
    view_costs: Vec<usize>,
}

impl AtomDiff {
    pub fn target_smaller(&self) -> bool {
        self.target_heavy < self.reactant_heavy
    }

    pub fn has_cleavage(&self) -> bool {
        !self.cleaved.is_empty() || !self.cleavage_bonds.is_empty()
    }

    pub fn h_loss(&self) -> bool {
        self.h_delta.values().any(|&d| d < 0)
    }

    pub fn h_gain(&self) -> bool {
        self.h_delta.values().any(|&d| d > 0)
    }

    pub fn cost(&self) -> usize {
        if !self.view_costs.is_empty() {
            return *self.view_costs.iter().min().unwrap_or(&0);
        }
        self.field_cost()
    }

    fn field_cost(&self) -> usize {
        // MCS map gaps only. H is not a special cost term: `formula_l1` counts
        // it like any element; `h_delta` stays on the struct for filters.
        3 * self.cleaved.len() + 3 * self.n_extra + 3 * self.cleavage_bonds.len()
    }

    /// True when ``atoms`` is the bond that separates kept from gone.
    pub fn site_is_cleavage(&self, atoms: &[usize]) -> bool {
        if atoms.is_empty() {
            return false;
        }
        let set: HashSet<usize> = atoms.iter().copied().collect();
        for &(a, b) in &self.cleavage_bonds {
            let bond: HashSet<usize> = [a, b].into_iter().collect();
            if bond.is_subset(&set) || set.is_subset(&bond) {
                return true;
            }
        }
        false
    }

    fn mappings_slice(&self) -> &[BTreeMap<usize, usize>] {
        if self.mappings.is_empty() {
            std::slice::from_ref(&self.mapping)
        } else {
            self.mappings.as_slice()
        }
    }
}

/// Bond order of an unmapped target O attached to the MCS image of `r`, if any.
fn unmapped_oxygen_order(
    target: &Molecule,
    mapping: &BTreeMap<usize, usize>,
    r: usize,
) -> Option<f32> {
    let &t = mapping.get(&r)?;
    let image: HashSet<usize> = mapping.values().copied().collect();
    for (nbr, bidx) in target.neighbors(atom_idx(t)) {
        let n = atom_usize(nbr);
        if target.atom(nbr).element.atomic_number() != 8 || image.contains(&n) {
            continue;
        }
        return Some(order_value(target.bond(bidx).order));
    }
    None
}

/// Reactant atom (or same-rank orbit mate) needs an O under this placement.
pub fn atom_needs_oxygen(
    reactant: &Molecule,
    target: &Molecule,
    mapping: &BTreeMap<usize, usize>,
    atom: usize,
) -> bool {
    let rank = ranks(reactant);
    let want = rank.get(atom).copied().unwrap_or(usize::MAX);
    for (i, &ri) in rank.iter().enumerate() {
        if ri == want && unmapped_oxygen_order(target, mapping, i).is_some() {
            return true;
        }
    }
    false
}

/// Attachment is carbonyl-like (bond order ≥ 1.5) under this placement.
pub fn atom_needs_carbonyl(
    reactant: &Molecule,
    target: &Molecule,
    mapping: &BTreeMap<usize, usize>,
    atom: usize,
) -> bool {
    let rank = ranks(reactant);
    let want = rank.get(atom).copied().unwrap_or(usize::MAX);
    for (i, &ri) in rank.iter().enumerate() {
        if ri == want {
            if let Some(order) = unmapped_oxygen_order(target, mapping, i) {
                if order >= 1.5 {
                    return true;
                }
            }
        }
    }
    false
}

/// Any MCS placement has an unmapped target O on the mapped core.
pub fn any_needs_oxygen(target: &Molecule, diff: &AtomDiff) -> bool {
    for mapping in diff.mappings_slice() {
        for &r in mapping.keys() {
            if unmapped_oxygen_order(target, mapping, r).is_some() {
                return true;
            }
        }
    }
    false
}

fn mapping_score(reactant: &Molecule, target: &Molecule, aligned: &BTreeMap<usize, usize>) -> i32 {
    let mut score = 0i32;
    for (&r, &t) in aligned {
        let ra = reactant.atom(atom_idx(r));
        let ta = target.atom(atom_idx(t));
        if ra.aromatic == ta.aromatic {
            score += 1;
        } else {
            score -= 2;
        }
    }
    for (&ri, &ti) in aligned {
        for (rj, bond) in reactant.neighbors(atom_idx(ri)) {
            let rj = atom_usize(rj);
            let Some(&tj) = aligned.get(&rj) else {
                continue;
            };
            if rj < ri {
                continue;
            }
            let Some((_, tbond)) = target.bond_between(atom_idx(ti), atom_idx(tj)) else {
                score -= 4;
                continue;
            };
            let delta = (order_value(tbond.order) - order_value(reactant.bond(bond).order)).abs();
            if delta < 0.2 {
                score += 1;
            } else {
                score -= 1;
            }
        }
    }
    score
}

/// All MCS placements: one alignment per distinct reactant atom set.
fn mappings(reactant: &Molecule, target: &Molecule) -> Vec<BTreeMap<usize, usize>> {
    let cfg = McsConfig {
        bond_compare: BondCompare::Any,
        // match_bonds=false so the QueryMolecule embeds on both aromatic and
        // kekulé writings (VF2 would otherwise refuse quinone targets).
        match_bonds: false,
        ..McsConfig::default()
    };
    let query = find_mcs_with_config(&[reactant, target], &cfg);
    if query.atom_count() == 0 {
        return Vec::new();
    }
    let r_hits = find_matches(&query, reactant);
    let t_hits = find_matches(&query, target);
    let mut by_size: HashMap<usize, Vec<BTreeMap<usize, usize>>> = HashMap::new();
    for emb in &t_hits {
        by_size
            .entry(emb.len())
            .or_default()
            .push(emb.iter().map(|(&q, &t)| (q, atom_usize(t))).collect());
    }
    if r_hits.is_empty() || by_size.is_empty() {
        return Vec::new();
    }

    let mut best: HashMap<BTreeSet<usize>, (i32, BTreeMap<usize, usize>)> = HashMap::new();
    for r_emb in &r_hits {
        let r_map: BTreeMap<usize, usize> =
            r_emb.iter().map(|(&q, &t)| (q, atom_usize(t))).collect();
        let Some(mates) = by_size.get(&r_map.len()) else {
            continue;
        };
        let key: BTreeSet<usize> = r_map.values().copied().collect();
        for t_map in mates {
            let mut aligned = BTreeMap::new();
            for (&q, &r_idx) in &r_map {
                if let Some(&t_idx) = t_map.get(&q) {
                    aligned.insert(r_idx, t_idx);
                }
            }
            if aligned.len() != r_map.len() {
                continue;
            }
            let score = mapping_score(reactant, target, &aligned);
            match best.get(&key) {
                Some((held, _)) if *held >= score => {}
                _ => {
                    best.insert(key.clone(), (score, aligned));
                }
            }
        }
    }
    let mut ranked: Vec<_> = best.into_values().collect();
    ranked.sort_by(|a, b| b.0.cmp(&a.0));
    ranked.into_iter().map(|(_, m)| m).collect()
}

fn expand_by_rank(mol: &Molecule, sites: BTreeSet<usize>) -> BTreeSet<usize> {
    let rank = ranks(mol);
    let mut out = BTreeSet::new();
    for s in sites {
        let r = rank.get(s).copied().unwrap_or(usize::MAX);
        for (i, &ri) in rank.iter().enumerate() {
            if ri == r {
                out.insert(i);
            }
        }
    }
    out
}

fn expand_h_delta(mol: &Molecule, h_delta: BTreeMap<usize, i32>) -> BTreeMap<usize, i32> {
    let rank = ranks(mol);
    let mut out = BTreeMap::new();
    for (&atom, &delta) in &h_delta {
        let r = rank.get(atom).copied().unwrap_or(usize::MAX);
        for (i, &ri) in rank.iter().enumerate() {
            if ri == r {
                // Keep the more negative (stronger H loss) when ranks collide.
                out.entry(i)
                    .and_modify(|d| {
                        if delta < *d {
                            *d = delta;
                        }
                    })
                    .or_insert(delta);
            }
        }
    }
    out
}

fn diff_for(reactant: &Molecule, target: &Molecule, mapping: &BTreeMap<usize, usize>) -> AtomDiff {
    let image: HashSet<usize> = mapping.values().copied().collect();
    let r_rings = ring_atom_flags(reactant);
    let t_rings = ring_atom_flags(target);

    let mut cleaved = BTreeSet::new();
    for (r_idx, atom) in reactant.atoms() {
        let r = atom_usize(r_idx);
        if atom.element.atomic_number() == 1 {
            continue;
        }
        if !mapping.contains_key(&r) {
            cleaved.insert(r);
        }
    }

    let mut cleavage_bonds = BTreeSet::new();
    let mut loses_aromaticity = BTreeSet::new();
    let mut h_delta = BTreeMap::new();
    let mut bond_raises = BTreeSet::new();
    let mut bond_order_mismatches = 0usize;

    for (&r_idx, &t_idx) in mapping {
        let ra = reactant.atom(atom_idx(r_idx));
        let ta = target.atom(atom_idx(t_idx));
        if ra.aromatic && !ta.aromatic {
            loses_aromaticity.insert(r_idx);
        }
        h_delta.insert(r_idx, hydrogens(target, t_idx) - hydrogens(reactant, r_idx));
    }

    for (_, bond) in reactant.bonds() {
        let i = atom_usize(bond.atom1);
        let j = atom_usize(bond.atom2);
        let i_mapped = mapping.contains_key(&i);
        let j_mapped = mapping.contains_key(&j);
        if !i_mapped || !j_mapped {
            if i_mapped != j_mapped {
                cleavage_bonds.insert(bond_key(i, j));
            }
            continue;
        }
        let ti = mapping[&i];
        let tj = mapping[&j];
        let Some((_, tbond)) = target.bond_between(atom_idx(ti), atom_idx(tj)) else {
            cleavage_bonds.insert(bond_key(i, j));
            continue;
        };
        let delta = order_value(tbond.order) - order_value(bond.order);
        if delta.abs() >= 0.2 {
            bond_order_mismatches += 1;
        }
        if delta >= 0.2 {
            bond_raises.insert(bond_key(i, j));
        }
    }

    // MCS may map exocyclic onto ring (PhCH2OH CH2 → quinone C); bridge is still the cut.
    for (&r_idx, &t_idx) in mapping {
        let r_in = r_rings.get(r_idx).copied().unwrap_or(false);
        let t_in = t_rings.get(t_idx).copied().unwrap_or(false);
        if r_in == t_in {
            continue;
        }
        for (nbr, _) in reactant.neighbors(atom_idx(r_idx)) {
            let n = atom_usize(nbr);
            let n_in = r_rings.get(n).copied().unwrap_or(false);
            if r_in == n_in {
                continue;
            }
            cleavage_bonds.insert(bond_key(r_idx, n));
        }
    }

    let n_extra = target
        .atoms()
        .filter(|(idx, a)| a.element.atomic_number() > 1 && !image.contains(&atom_usize(*idx)))
        .count();

    let loses_aromaticity = expand_by_rank(reactant, loses_aromaticity);
    let h_delta = expand_h_delta(reactant, h_delta);

    let mut diff = AtomDiff {
        mapping: mapping.clone(),
        mappings: vec![mapping.clone()],
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta,
        bond_raises,
        bond_order_mismatches,
        n_extra,
        reactant_heavy: heavy_atom_count(reactant),
        target_heavy: heavy_atom_count(target),
        view_costs: Vec::new(),
    };
    diff.view_costs = vec![diff.field_cost()];
    diff
}

fn merge_views(mut views: Vec<AtomDiff>) -> AtomDiff {
    let costs: Vec<usize> = views.iter().map(|v| v.field_cost()).collect();
    let mappings: Vec<_> = views.iter().map(|v| v.mapping.clone()).collect();
    let mut primary = views.remove(0);
    primary.view_costs = costs;
    primary.mappings = mappings;
    if views.is_empty() {
        return primary;
    }
    for view in &views {
        primary.cleaved.extend(&view.cleaved);
        primary.cleavage_bonds.extend(&view.cleavage_bonds);
        primary.bond_raises.extend(&view.bond_raises);
        primary.loses_aromaticity.extend(&view.loses_aromaticity);
        for (&atom, &delta) in &view.h_delta {
            primary
                .h_delta
                .entry(atom)
                .and_modify(|d| {
                    if delta < *d {
                        *d = delta;
                    }
                })
                .or_insert(delta);
        }
        primary.bond_order_mismatches = primary
            .bond_order_mismatches
            .max(view.bond_order_mismatches);
        primary.n_extra = primary.n_extra.max(view.n_extra);
    }
    primary
}

/// Pair reactant atoms with target atoms and record the local change.
pub fn atom_diff(reactant: &Molecule, target: &Molecule) -> AtomDiff {
    let maps = mappings(reactant, target);
    if maps.is_empty() {
        return AtomDiff {
            reactant_heavy: heavy_atom_count(reactant),
            target_heavy: heavy_atom_count(target),
            ..AtomDiff::default()
        };
    }
    let views: Vec<_> = maps.iter().map(|m| diff_for(reactant, target, m)).collect();
    merge_views(views)
}

/// Build [`AtomDiff`] from known reactant→target mappings (no MCS).
pub fn atom_diff_from_mappings(
    reactant: &Molecule,
    target: &Molecule,
    maps: Vec<BTreeMap<usize, usize>>,
) -> AtomDiff {
    if maps.is_empty() {
        return AtomDiff {
            reactant_heavy: heavy_atom_count(reactant),
            target_heavy: heavy_atom_count(target),
            ..AtomDiff::default()
        };
    }
    let views: Vec<_> = maps.iter().map(|m| diff_for(reactant, target, m)).collect();
    merge_views(views)
}

/// Lift parent MCS mappings onto a tagged child via surviving atom tags.
///
/// Removed atoms drop out of each map (shrink). Added atoms are left
/// unmapped — callers extend via [`extend_mapping_where_possible`] and
/// product automorphism generators ([`best_diff_from_lifted_maps`]).
pub fn lift_mappings(
    parent: &crate::forest_mol::ForestMol,
    child: &crate::forest_mol::ForestMol,
    parent_maps: &[BTreeMap<usize, usize>],
) -> Option<Vec<BTreeMap<usize, usize>>> {
    if !parent.shares_tag_gen(child) || parent_maps.is_empty() {
        return None;
    }
    let mut lifted = Vec::new();
    for parent_map in parent_maps {
        let mut child_map = BTreeMap::new();
        for (&parent_idx, &target_idx) in parent_map {
            let Some(tag) = parent.tag_of(parent_idx) else {
                continue;
            };
            let Some(child_idx) = child.index_of(tag) else {
                continue;
            };
            child_map.insert(child_idx, target_idx);
        }
        if !child_map.is_empty() {
            lifted.push(child_map);
        }
    }
    if lifted.is_empty() {
        None
    } else {
        Some(lifted)
    }
}

/// Unmapped heavy atoms on `mol` (candidates to grow the alignment).
fn unmapped_heavy_atoms(mol: &Molecule, mapping: &BTreeMap<usize, usize>) -> Vec<usize> {
    let mut out: Vec<usize> = (0..mol.atom_count())
        .filter(|&i| mol.atom(atom_idx(i)).element.atomic_number() > 1 && !mapping.contains_key(&i))
        .collect();
    out.sort_unstable();
    out
}

/// Grow a product→target map where the gap is placeable (diff `n_extra` /
/// unmapped heavies adjacent to the mapped core).
///
/// Same placement rules as [`extend_mapping_for_added`]: free same-element
/// target atoms bonded to every mapped neighbor image (with single-neighbor
/// rematch when MCS orientation is flipped).
pub fn extend_mapping_where_possible(
    child: &Molecule,
    target: &Molecule,
    mapping: &mut BTreeMap<usize, usize>,
) {
    let pending = unmapped_heavy_atoms(child, mapping);
    if pending.is_empty() {
        return;
    }
    // Nothing to place onto if the target image already covers every heavy.
    let image: HashSet<usize> = mapping.values().copied().collect();
    let target_free = (0..target.atom_count())
        .any(|t| target.atom(atom_idx(t)).element.atomic_number() > 1 && !image.contains(&t));
    if !target_free {
        return;
    }
    extend_mapping_for_added(child, target, mapping, &pending);
}

/// Best [`AtomDiff`] from lifted seeds: **extend**; keep only a **cost-0** lift.
///
/// 1. Extend each tag-lifted seed where the mapping gap is placeable.
/// 2. If any extended map has `field_cost == 0`, keep it (guaranteed — skip MCS).
/// 3. Otherwise do **not** Aut-chase or trust a non-zero lift — return a fresh
///    MCS (caller counts `mcs_lift_rematch`). No double work.
fn best_diff_from_lifted_maps(
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    seed_maps: Vec<BTreeMap<usize, usize>>,
) -> Option<(AtomDiff, bool)> {
    if seed_maps.is_empty() {
        return None;
    }
    let mol = child.mol();
    let mut zero_maps: Vec<BTreeMap<usize, usize>> = Vec::new();
    let mut any_extended = false;

    for mut seed in seed_maps {
        extend_mapping_where_possible(mol, target, &mut seed);
        if seed.is_empty() {
            continue;
        }
        any_extended = true;
        if diff_for(mol, target, &seed).field_cost() == 0 {
            zero_maps.push(seed);
        }
    }
    if !any_extended {
        return None;
    }
    if !zero_maps.is_empty() {
        return Some((atom_diff_from_mappings(mol, target, zero_maps), false));
    }
    // Not sure → MCS only (skip Aut / non-zero lift).
    Some((atom_diff(mol, target), true))
}

/// Lower bound sketch on child [`AtomDiff::cost`] after this effect at `site`.
///
/// Casts the effect onto the parent MCS-gap cost: clear `n_extra` and cleavage
/// bond + leave heavies the effect is declared to fix on
/// `site_atoms ∪ path_ends`. Same weights as [`AtomDiff::field_cost`].
/// H is not a cost term (`formula_l1` / filter `h_delta`). Useful for tests /
/// filters. Not used to decide lift vs MCS.
pub fn residual_cost_after_site_cast(
    parent: &AtomDiff,
    effect: &Effect,
    site_atoms: &[usize],
    path_ends: &[usize],
) -> usize {
    let scope: HashSet<usize> = site_atoms.iter().chain(path_ends.iter()).copied().collect();

    let mut cleaved = parent.cleaved.clone();
    let mut cleavage_bonds = parent.cleavage_bonds.clone();
    let mut n_extra = parent.n_extra;

    if effect_adds_oxygen(effect) {
        let o_delta = effect.delta_formula.get("O").copied().unwrap_or(0).max(0) as usize;
        let o_place = o_delta.max(1);
        n_extra = n_extra.saturating_sub(o_place);
    }
    if effect.cleaves {
        cleavage_bonds.retain(|&(a, b)| !(scope.contains(&a) && scope.contains(&b)));
        let leave_heavies = effect
            .leave_count
            .map(|n| n as usize)
            .or_else(|| {
                let n = effect
                    .leave_formula
                    .values()
                    .filter(|&&c| c > 0)
                    .map(|&c| c as usize)
                    .sum::<usize>();
                if n > 0 { Some(n) } else { None }
            })
            .unwrap_or(0);
        if leave_heavies > 0 {
            let mut drop = leave_heavies;
            cleaved.retain(|_| {
                if drop > 0 {
                    drop -= 1;
                    false
                } else {
                    true
                }
            });
        }
    }

    3 * cleaved.len() + 3 * n_extra + 3 * cleavage_bonds.len()
}

/// Heavy child atoms whose tags are not on `parent` (local additions).
pub fn added_heavy_atoms(
    parent: &crate::forest_mol::ForestMol,
    child: &crate::forest_mol::ForestMol,
) -> Vec<usize> {
    let mut out = Vec::new();
    for i in 0..child.mol().atom_count() {
        if child.mol().atom(atom_idx(i)).element.atomic_number() <= 1 {
            continue;
        }
        let Some(tag) = child.tag_of(i) else {
            out.push(i);
            continue;
        };
        if parent.index_of(tag).is_none() {
            out.push(i);
        }
    }
    out
}

/// Place locally added child atoms onto free target atoms of the same element.
///
/// Rules only add/remove atoms at the edit site. Prefer a free target atom
/// bonded to every mapped neighbor's image. If a single-neighbor add does not
/// fit the current MCS orientation, rematch that neighbor onto a same-element
/// target atom adjacent to the candidate (swap when needed).
pub fn extend_mapping_for_added(
    child: &Molecule,
    target: &Molecule,
    mapping: &mut BTreeMap<usize, usize>,
    added: &[usize],
) {
    let mut image: HashSet<usize> = mapping.values().copied().collect();
    let mut pending: Vec<usize> = added
        .iter()
        .copied()
        .filter(|i| !mapping.contains_key(i))
        .collect();
    pending.sort_unstable();
    let mut progress = true;
    while progress {
        progress = false;
        let mut still = Vec::new();
        for child_idx in pending {
            if mapping.contains_key(&child_idx) {
                continue;
            }
            match place_added_atom(child, target, mapping, &mut image, child_idx) {
                Some(t_cand) => {
                    mapping.insert(child_idx, t_cand);
                    image.insert(t_cand);
                    progress = true;
                }
                None => still.push(child_idx),
            }
        }
        pending = still;
    }
}

fn place_added_atom(
    child: &Molecule,
    target: &Molecule,
    mapping: &mut BTreeMap<usize, usize>,
    image: &mut HashSet<usize>,
    child_idx: usize,
) -> Option<usize> {
    let z = child.atom(atom_idx(child_idx)).element.atomic_number();
    let mut mapped_nbrs: Vec<usize> = child
        .neighbors(atom_idx(child_idx))
        .map(|(nbr, _)| atom_usize(nbr))
        .filter(|n| mapping.contains_key(n))
        .collect();
    mapped_nbrs.sort_unstable();
    mapped_nbrs.dedup();
    if mapped_nbrs.is_empty() {
        return None;
    }

    // Free target atom of matching Z adjacent to every mapped neighbor image.
    let mut direct = Vec::new();
    for t_cand in 0..target.atom_count() {
        if image.contains(&t_cand) {
            continue;
        }
        if target.atom(atom_idx(t_cand)).element.atomic_number() != z {
            continue;
        }
        let ok = mapped_nbrs.iter().all(|&n| {
            let t_n = mapping[&n];
            target
                .bond_between(atom_idx(t_cand), atom_idx(t_n))
                .is_some()
        });
        if ok {
            direct.push(t_cand);
        }
    }
    if let Some(&best) = direct.iter().min() {
        return Some(best);
    }

    // Single-neighbor add: MCS orientation may be flipped on a symmetric
    // parent. Rematch the neighbor onto a same-element target atom adjacent
    // to a free candidate for the new atom.
    if mapped_nbrs.len() != 1 {
        return None;
    }
    let n = mapped_nbrs[0];
    let z_n = child.atom(atom_idx(n)).element.atomic_number();
    let old_t_n = mapping[&n];
    let mut repairs = Vec::new();
    for t_cand in 0..target.atom_count() {
        if image.contains(&t_cand) {
            continue;
        }
        if target.atom(atom_idx(t_cand)).element.atomic_number() != z {
            continue;
        }
        for (t_n_idx, _) in target.neighbors(atom_idx(t_cand)) {
            let t_n = atom_usize(t_n_idx);
            if target.atom(atom_idx(t_n)).element.atomic_number() != z_n {
                continue;
            }
            if t_n == old_t_n {
                repairs.push((t_cand, None));
                continue;
            }
            if !image.contains(&t_n) {
                repairs.push((t_cand, Some((t_n, None))));
                continue;
            }
            let Some((&m, _)) = mapping.iter().find(|&(_, &img)| img == t_n) else {
                continue;
            };
            if child.atom(atom_idx(m)).element.atomic_number() != z_n {
                continue;
            }
            // Swap images of n and m so n sits next to t_cand.
            repairs.push((t_cand, Some((t_n, Some(m)))));
        }
    }
    repairs.sort_unstable_by_key(|(t, repair)| (*t, repair.is_some()));
    let (t_cand, repair) = repairs.into_iter().next()?;
    match repair {
        None => {}
        Some((t_n, None)) => {
            image.remove(&old_t_n);
            mapping.insert(n, t_n);
            image.insert(t_n);
        }
        Some((t_n, Some(m))) => {
            mapping.insert(n, t_n);
            mapping.insert(m, old_t_n);
        }
    }
    Some(t_cand)
}

/// Child [`AtomDiff`] via tag-lift + extend when the result is cost `0`.
///
/// Same heavy-tag set (no shrink). Only a cost-0 extended lift is kept;
/// otherwise MCS ([`try_atom_diff_for_child_tracked`] counts). Cleavage
/// shrinks: use [`try_lift_cleaved_child`].
pub fn try_atom_diff_for_child(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> Option<AtomDiff> {
    try_atom_diff_for_child_tracked(parent, parent_diff, child, target, None)
}

/// Like [`try_atom_diff_for_child`], optionally counting MCS rematches.
pub fn try_atom_diff_for_child_tracked(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    mut mcs_rematch: Option<&mut usize>,
) -> Option<AtomDiff> {
    // Any heavy parent atom missing on the child → shrink; closer cost from a
    // non-cleavage lift is unreliable. Cleavage expand uses
    // [`try_lift_cleaved_child`] instead.
    for i in 0..parent.mol().atom_count() {
        if parent.mol().atom(atom_idx(i)).element.atomic_number() <= 1 {
            continue;
        }
        let Some(tag) = parent.tag_of(i) else {
            continue;
        };
        child.index_of(tag)?;
    }
    let parent_maps = if parent_diff.mappings.is_empty() {
        vec![parent_diff.mapping.clone()]
    } else {
        parent_diff.mappings.clone()
    };
    let lifted = lift_mappings(parent, child, &parent_maps)?;
    let (diff, used_mcs) = best_diff_from_lifted_maps(child, target, lifted)?;
    if used_mcs {
        if let Some(c) = mcs_rematch.as_mut() {
            **c += 1;
        }
    }
    Some(diff)
}

/// Deprecated alias: `goal_cost` is ignored (cost-0 lift or MCS).
pub fn try_atom_diff_for_child_goal(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    _goal_cost: Option<usize>,
) -> Option<AtomDiff> {
    try_atom_diff_for_child(parent, parent_diff, child, target)
}

/// Deprecated alias for [`try_atom_diff_for_child_tracked`].
pub fn try_atom_diff_for_child_goal_tracked(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    _goal_cost: Option<usize>,
    mcs_rematch: Option<&mut usize>,
) -> Option<AtomDiff> {
    try_atom_diff_for_child_tracked(parent, parent_diff, child, target, mcs_rematch)
}

/// Diff when a child shares tags but carries none of the parent MCS image
/// (typical discarded cleavage side). High cost; no MCS rematch.
fn unmapped_child_diff(child: &Molecule, target: &Molecule) -> AtomDiff {
    let mut cleaved = BTreeSet::new();
    for (idx, atom) in child.atoms() {
        if atom.element.atomic_number() > 1 {
            cleaved.insert(atom_usize(idx));
        }
    }
    let mut diff = AtomDiff {
        cleaved,
        n_extra: heavy_atom_count(target),
        reactant_heavy: heavy_atom_count(child),
        target_heavy: heavy_atom_count(target),
        ..AtomDiff::default()
    };
    diff.view_costs = vec![diff.field_cost()];
    diff
}

/// Tag-lift after a cleavage shrink (and optional local adds on a fragment).
///
/// Removed atoms drop out of the parent MCS map via [`lift_mappings`];
/// [`extend_mapping_where_possible`] may hit cost `0` (kept); otherwise MCS.
///
/// When the child shares tags but inherits **no** mapped atoms (discarded
/// cleavage side vs an MCS that lives on the other fragment), returns a
/// high-cost unmapped shell (no MCS). `None` only when tags are not shared.
pub fn try_lift_cleaved_child(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> Option<AtomDiff> {
    try_lift_cleaved_child_tracked(parent, parent_diff, child, target, None)
}

/// Like [`try_lift_cleaved_child`], optionally counting MCS rematches.
pub fn try_lift_cleaved_child_tracked(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    mut mcs_rematch: Option<&mut usize>,
) -> Option<AtomDiff> {
    if !parent.shares_tag_gen(child) {
        return None;
    }
    let parent_maps = if parent_diff.mappings.is_empty() {
        vec![parent_diff.mapping.clone()]
    } else {
        parent_diff.mappings.clone()
    };
    match lift_mappings(parent, child, &parent_maps) {
        Some(lifted) => {
            let (diff, used_mcs) = best_diff_from_lifted_maps(child, target, lifted)?;
            if used_mcs {
                if let Some(c) = mcs_rematch.as_mut() {
                    **c += 1;
                }
            }
            Some(diff)
        }
        None => Some(unmapped_child_diff(child.mol(), target)),
    }
}

/// Deprecated alias: `goal_cost` is ignored.
pub fn try_lift_cleaved_child_goal(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    _goal_cost: Option<usize>,
) -> Option<AtomDiff> {
    try_lift_cleaved_child(parent, parent_diff, child, target)
}

/// Deprecated alias for [`try_lift_cleaved_child_tracked`].
pub fn try_lift_cleaved_child_goal_tracked(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    _goal_cost: Option<usize>,
    mcs_rematch: Option<&mut usize>,
) -> Option<AtomDiff> {
    try_lift_cleaved_child_tracked(parent, parent_diff, child, target, mcs_rematch)
}

/// Child diff after cleavage: prefer cost-0 tag-lift+extend; else MCS
/// (`mcs_fallback` when lift is impossible).
///
/// Callers apply the expand gate (`cost < parent` / target hit) on the result.
pub fn atom_diff_after_cleavage(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> AtomDiff {
    atom_diff_after_cleavage_tracked(parent, parent_diff, child, target, None, None)
}

/// Like [`atom_diff_after_cleavage`], with optional counters.
///
/// - `mcs_fallback`: lift impossible (no shared tags).
/// - `mcs_rematch`: extended lift was not cost 0 → used MCS.
pub fn atom_diff_after_cleavage_tracked(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
    mut mcs_fallback: Option<&mut usize>,
    mcs_rematch: Option<&mut usize>,
) -> AtomDiff {
    if let Some(lifted) =
        try_lift_cleaved_child_tracked(parent, parent_diff, child, target, mcs_rematch)
    {
        return lifted;
    }
    if let Some(c) = mcs_fallback.as_mut() {
        **c += 1;
    }
    atom_diff(child.mol(), target)
}

/// Child [`AtomDiff`] via tag-lift + generators when possible; else full MCS.
///
/// Prefer [`try_atom_diff_for_child`] + an explicit MCS-fallback counter in
/// search. This helper keeps the profile / unit-test door.
pub fn atom_diff_for_child(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> AtomDiff {
    try_atom_diff_for_child(parent, parent_diff, child, target)
        .unwrap_or_else(|| atom_diff(child.mol(), target))
}

pub fn effect_adds_oxygen(effect: &Effect) -> bool {
    effect.adds.as_deref().is_some_and(|a| a.contains('O'))
}

/// Dehydration-style: effect removes O / OH (inverse of hydration).
pub fn effect_removes_oxygen(effect: &Effect) -> bool {
    effect.removes.as_deref().is_some_and(|r| r.contains('O'))
}

fn effect_adds_h(effect: &Effect) -> bool {
    effect.adds.as_deref().is_some_and(|a| a.contains('H'))
}

fn effect_removes_h(effect: &Effect) -> bool {
    effect.removes.as_deref().is_some_and(|r| r.contains('H'))
}

/// Same as [`pattern_could_help`], with optional live mol / target for O gates.
pub fn pattern_could_help_on(
    effect: &Effect,
    diff: &AtomDiff,
    mol: Option<&Molecule>,
    target: Option<&Molecule>,
) -> bool {
    let can_cleave = effect.cleaves;
    if diff.target_smaller() && !can_cleave {
        return false;
    }
    if effect_adds_oxygen(effect) && !can_cleave {
        match target {
            Some(t) if !any_needs_oxygen(t, diff) => return false,
            None => {}
            Some(_) => {}
        }
    }
    if effect_adds_oxygen(effect) && !can_cleave && !effect.dearomatizes {
        if let Some(m) = mol {
            let _ = m;
        }
    }
    if can_cleave && !diff.has_cleavage() {
        return false;
    }
    if effect.dearomatizes && !effect_adds_h(effect) && diff.loses_aromaticity.is_empty() {
        return false;
    }
    let drops_h_only = !can_cleave && effect.adds.is_none() && effect_removes_h(effect);
    if drops_h_only && !diff.h_loss() && diff.loses_aromaticity.is_empty() {
        return false;
    }
    let adds_h_only = !can_cleave && effect_adds_h(effect) && effect.removes.is_none();
    if adds_h_only && !diff.h_gain() {
        return false;
    }
    true
}

/// Pattern-level gate (Python `_pattern_could_help`).
pub fn pattern_could_help(effect: &Effect, diff: &AtomDiff) -> bool {
    pattern_could_help_on(effect, diff, None, None)
}

/// Live-mol oxygen gate (Python formula check).
pub fn pattern_could_help_mol(
    effect: &Effect,
    diff: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> bool {
    if !pattern_could_help_on(effect, diff, Some(mol), Some(target)) {
        return false;
    }
    if effect_adds_oxygen(effect)
        && !effect.cleaves
        && !effect.dearomatizes
        && formula_oxygen(mol) >= formula_oxygen(target)
    {
        return false;
    }
    true
}

fn site_atoms(candidate: &Candidate) -> Vec<usize> {
    let mut atoms: Vec<usize> = candidate
        .pattern
        .site_map
        .iter()
        .filter_map(|m| candidate.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.push(candidate.site);
    }
    atoms
}

fn leaving_heavy_counts(mol: &Molecule, atoms: &[usize]) -> Option<(usize, usize)> {
    if atoms.len() != 2 {
        return None;
    }
    let left = atoms[0];
    let right = atoms[1];
    mol.bond_between(atom_idx(left), atom_idx(right))?;
    let side = |start: usize, blocked: usize| -> usize {
        let mut seen = HashSet::from([start]);
        let mut stack = vec![start];
        while let Some(idx) = stack.pop() {
            for (nbr, _) in mol.neighbors(atom_idx(idx)) {
                let n = atom_usize(nbr);
                if n == blocked || !seen.insert(n) {
                    continue;
                }
                stack.push(n);
            }
        }
        seen.into_iter()
            .filter(|&i| mol.atom(atom_idx(i)).element.atomic_number() > 1)
            .count()
    };
    Some((side(left, right), side(right, left)))
}

fn alkyl_bond_raises(mol: &Molecule, atom_idx_u: usize, diff: &AtomDiff) -> bool {
    for (nbr, _) in mol.neighbors(atom_idx(atom_idx_u)) {
        let n = atom_usize(nbr);
        let atom = mol.atom(atom_idx(n));
        if atom.element.atomic_number() != 6 || atom.aromatic {
            continue;
        }
        if diff.bond_raises.contains(&bond_key(atom_idx_u, n)) {
            return true;
        }
    }
    false
}

/// Dearomatizing removes-H edit (DH / QF path ends) — not cleavage or OH.
pub fn is_dehydrogenation_effect(effect: &Effect) -> bool {
    effect_removes_h(effect)
        && effect.dearomatizes
        && !effect.cleaves
        && !effect_adds_oxygen(effect)
}

fn heavy_neighbor_idxs(mol: &Molecule, atom: usize) -> BTreeSet<usize> {
    mol.neighbors(atom_idx(atom))
        .filter(|(n, _)| mol.atom(*n).element.atomic_number() > 1)
        .map(|(n, _)| atom_usize(n))
        .collect()
}

fn dh_site_neighbors_match(
    mol: &Molecule,
    target: &Molecule,
    atom: usize,
    mapping: &BTreeMap<usize, usize>,
) -> bool {
    let Some(&t_idx) = mapping.get(&atom) else {
        return false;
    };
    let mut imaged = BTreeSet::new();
    for n in heavy_neighbor_idxs(mol, atom) {
        let Some(&t_n) = mapping.get(&n) else {
            return false;
        };
        imaged.insert(t_n);
    }
    imaged == heavy_neighbor_idxs(target, t_idx)
}

/// True when some MCS view has matching heavy neighbors for this site atom.
pub fn dh_neighbors_match_any_view(
    mol: &Molecule,
    target: &Molecule,
    atom: usize,
    diff: &AtomDiff,
) -> bool {
    for mapping in &diff.mappings {
        if !mapping.contains_key(&atom) {
            continue;
        }
        if dh_site_neighbors_match(mol, target, atom, mapping) {
            return true;
        }
    }
    false
}

/// After DH: each **product** end's heavy neighbors match the target.
///
/// `product_ends` are already mapped from the reactant site onto the product
/// (via forest tags). Pre-application reactant connectivity is not enough.
///
/// Both ends must match under the **same** MCS mapping — not each end's best
/// independent top-group / view match (optimistic split vs pair topology).
pub fn dh_product_ends_match(
    product: &Molecule,
    product_ends: &[usize],
    target: &Molecule,
) -> bool {
    if product_ends.is_empty() {
        return true;
    }
    let diff = atom_diff(product, target);
    for mapping in &diff.mappings {
        if product_ends.iter().all(|&atom| {
            mapping.contains_key(&atom) && dh_site_neighbors_match(product, target, atom, mapping)
        }) {
            return true;
        }
    }
    false
}

/// Site-level gate (Python `_site_could_help`) for a deferred candidate.
pub fn candidate_could_help(candidate: &Candidate, diff: &AtomDiff) -> bool {
    candidate_could_help_on(candidate, diff, None, None)
}

fn scope_could_help(
    effect: &Effect,
    atoms: &[usize],
    path_ends: &[usize],
    diff: &AtomDiff,
    mol: Option<&Molecule>,
) -> bool {
    let mut scope: HashSet<usize> = atoms.iter().copied().collect();
    scope.extend(path_ends.iter().copied());
    // Carbonyl / imine reduction sites the heteroatom; H change is on the
    // partner heavy atom. Include double-bond neighbors so adds-H sees whether
    // applying helps (undo would not).
    if let Some(m) = mol {
        extend_h_edit_partners(m, &mut scope);
    }
    if effect.dearomatizes && !scope.iter().any(|a| diff.loses_aromaticity.contains(a)) {
        return false;
    }
    if effect_removes_h(effect) && !effect_adds_oxygen(effect) && !effect.cleaves {
        let loses_h = scope
            .iter()
            .any(|a| diff.h_delta.get(a).copied().unwrap_or(0) < 0);
        let loses_ar = scope.iter().any(|a| diff.loses_aromaticity.contains(a));
        if !loses_h && !loses_ar {
            return false;
        }
    }
    // Adding H helps only where the target needs more H. Then undoing (remove H)
    // would move away from the target — so it would not help.
    if effect_adds_h(effect) && !effect_adds_oxygen(effect) && !effect.cleaves {
        let gains = scope
            .iter()
            .any(|a| diff.h_delta.get(a).copied().unwrap_or(0) > 0);
        if !gains {
            return false;
        }
    }
    true
}

/// Double-bond partners of O/N in `scope` (carbonyl / imine reduction).
fn extend_h_edit_partners(mol: &Molecule, scope: &mut HashSet<usize>) {
    let seeds: Vec<usize> = scope.iter().copied().collect();
    for a in seeds {
        let z = mol.atom(atom_idx(a)).element.atomic_number();
        if z != 7 && z != 8 {
            continue;
        }
        for (nbr, bidx) in mol.neighbors(atom_idx(a)) {
            let order = mol.bond(bidx).order;
            if matches!(
                order,
                chematic::core::BondOrder::Double | chematic::core::BondOrder::Aromatic
            ) {
                scope.insert(atom_usize(nbr));
            }
        }
    }
}

/// Full site gate with live mol (leave_count / methide partner).
///
/// Non-cleavage: the whole site must help under **one** MCS placement (same
/// rule for atom / bond / pair). Merged top-rank unions across placements are
/// optimistic when the site has more than one atom.
pub fn candidate_could_help_on(
    candidate: &Candidate,
    diff: &AtomDiff,
    mol: Option<&Molecule>,
    target: Option<&Molecule>,
) -> bool {
    let effect = &candidate.pattern.effect;
    let ok = match (mol, target) {
        (Some(m), Some(t)) => pattern_could_help_mol(effect, diff, m, t),
        _ => pattern_could_help_on(effect, diff, mol, target),
    };
    if !ok {
        return false;
    }
    let atoms = site_atoms(candidate);
    if effect.cleaves {
        if !diff.site_is_cleavage(&atoms) {
            return false;
        }
        if let (Some(n), Some(m)) = (effect.leave_count, mol) {
            if let Some((a, b)) = leaving_heavy_counts(m, &atoms) {
                if a.min(b) != n as usize {
                    return false;
                }
            }
        }
        return true;
    }

    match (mol, target) {
        (Some(m), Some(t)) => {
            for mapping in diff.mappings_slice() {
                let view = diff_for(m, t, mapping);
                if candidate_could_help_on_view(candidate, &atoms, &view, Some(m), Some(t)) {
                    return true;
                }
            }
            false
        }
        _ => candidate_could_help_on_view(candidate, &atoms, diff, mol, target),
    }
}

fn candidate_could_help_on_view(
    candidate: &Candidate,
    atoms: &[usize],
    view: &AtomDiff,
    mol: Option<&Molecule>,
    target: Option<&Molecule>,
) -> bool {
    let effect = &candidate.pattern.effect;
    if effect_adds_oxygen(effect) && !effect.dearomatizes {
        let (Some(m), Some(t)) = (mol, target) else {
            return false;
        };
        let oxygen_sites: Vec<_> = atoms
            .iter()
            .copied()
            .filter(|&a| atom_needs_oxygen(m, t, &view.mapping, a))
            .collect();
        if oxygen_sites.is_empty() {
            return false;
        }
        if oxygen_sites.iter().all(|&a| {
            atom_needs_carbonyl(m, t, &view.mapping, a) && view.loses_aromaticity.contains(&a)
        }) {
            return false;
        }
    }

    // Single-site methide: alkyl partner needs an exocyclic C–C bond raise.
    // Pair ends use [`pair_could_help`] (per-end partner), not this merge.
    if effect.partner.as_deref() == Some("C") {
        if let Some(m) = mol {
            if !atoms.iter().any(|&a| alkyl_bond_raises(m, a, view)) {
                return false;
            }
        }
    }

    scope_could_help(effect, atoms, &[], view, mol)
}

/// Pair-site gate (Python `_site_could_help` when ``"ends"`` is on the info).
///
/// Same one-placement rule as [`candidate_could_help_on`]. Per-end oxygen /
/// ``partner == "C"`` read each end's effect (pair data), not the merged span.
pub fn pair_could_help(
    pair: &PairCandidate,
    diff: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> bool {
    let effect = &pair.effect;
    if !pattern_could_help_mol(effect, diff, mol, target) {
        return false;
    }
    let Some((end_a, end_b)) = pair.end_atoms() else {
        return false;
    };
    let atoms = [end_a, end_b];
    if effect.cleaves {
        if !diff.site_is_cleavage(&atoms) {
            return false;
        }
        if let Some(n) = effect.leave_count {
            if let Some((a, b)) = leaving_heavy_counts(mol, &atoms) {
                if a.min(b) != n as usize {
                    return false;
                }
            }
        }
        return true;
    }

    let (p0, p1) = pair.path_ends();
    let mappings = if diff.mappings.is_empty() {
        std::slice::from_ref(&diff.mapping)
    } else {
        diff.mappings.as_slice()
    };
    for mapping in mappings {
        let view = diff_for(mol, target, mapping);
        if !pair_ends_match_view(pair, end_a, end_b, &view, mol, target) {
            continue;
        }
        if scope_could_help(effect, &atoms, &[p0, p1], &view, Some(mol)) {
            return true;
        }
    }
    false
}

fn pair_ends_match_view(
    pair: &PairCandidate,
    end_a: usize,
    end_b: usize,
    view: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> bool {
    let ends = [(&pair.left.effect, end_a), (&pair.right.effect, end_b)];
    for (end, atom) in ends {
        if effect_adds_oxygen(end) && !atom_needs_oxygen(mol, target, &view.mapping, atom) {
            return false;
        }
        if end.partner.as_deref() == Some("C") && !alkyl_bond_raises(mol, atom, view) {
            return false;
        }
    }
    true
}

/// Sort key for expand: cleavage / dearom / oxygen first (Python `order_key`).
///
/// Fourth field is negated site H-progress (higher progress sorts earlier): adds-H
/// at atoms that need H, removes-H where the target loses H. Carbonyl / imine
/// reduction partners are included when `mol` is given so OxygenReduction's
/// O-only site still scores the carbon that gains H.
pub fn candidate_order_key(candidate: &Candidate, diff: &AtomDiff) -> (u8, u8, u8, String) {
    let (a, b, c, _progress, name) = candidate_order_key_on(candidate, diff, None, None);
    (a, b, c, name)
}

/// Same as [`candidate_order_key`], with live mol / target for H-progress and O preference.
pub fn candidate_order_key_on(
    candidate: &Candidate,
    diff: &AtomDiff,
    mol: Option<&Molecule>,
    target: Option<&Molecule>,
) -> (u8, u8, u8, i32, String) {
    let effect = &candidate.pattern.effect;
    let want_cleave = diff.target_smaller() || diff.has_cleavage();
    let want_dear = !diff.loses_aromaticity.is_empty();
    let want_oxy = target.is_some_and(|t| any_needs_oxygen(t, diff));
    let cleave = if effect.cleaves { 0 } else { 1 };
    let dear = if effect.dearomatizes { 0 } else { 1 };
    let oxy = if effect_adds_oxygen(effect) { 0 } else { 1 };
    let primary = if want_cleave { cleave } else { 0 };
    let secondary = if want_dear { dear } else { 0 };
    let tertiary = if want_oxy { oxy } else { 0 };
    let atoms = site_atoms(candidate);
    let progress = site_h_progress(effect, &atoms, &[], diff, mol);
    // Negate so ascending sort prefers higher progress (apply helps more).
    (
        primary,
        secondary,
        tertiary,
        -progress,
        candidate.pattern.name.clone(),
    )
}

/// How much applying this H-direction effect helps at `atoms` ∪ path ends
/// (plus carbonyl partners). Positive ⇒ applying moves toward the target;
/// undo would move away (and would not pass [`scope_could_help`]).
pub fn site_h_progress(
    effect: &Effect,
    atoms: &[usize],
    path_ends: &[usize],
    diff: &AtomDiff,
    mol: Option<&Molecule>,
) -> i32 {
    if effect.cleaves || effect_adds_oxygen(effect) {
        return 0;
    }
    let mut scope: HashSet<usize> = atoms.iter().copied().collect();
    scope.extend(path_ends.iter().copied());
    if let Some(m) = mol {
        extend_h_edit_partners(m, &mut scope);
    }
    let mut progress = 0i32;
    if effect_adds_h(effect) {
        for a in &scope {
            let d = diff.h_delta.get(a).copied().unwrap_or(0);
            if d > 0 {
                progress += d;
            }
        }
    } else if effect_removes_h(effect) {
        for a in &scope {
            let d = diff.h_delta.get(a).copied().unwrap_or(0);
            if d < 0 {
                progress += -d;
            }
        }
    }
    progress
}

/// Best [`site_h_progress`] under any one MCS placement (not merged h_delta).
pub fn site_h_progress_best_placement(
    effect: &Effect,
    atoms: &[usize],
    path_ends: &[usize],
    diff: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> i32 {
    let mappings = if diff.mappings.is_empty() {
        std::slice::from_ref(&diff.mapping)
    } else {
        diff.mappings.as_slice()
    };
    mappings
        .iter()
        .map(|mapping| {
            let view = diff_for(mol, target, mapping);
            site_h_progress(effect, atoms, path_ends, &view, Some(mol))
        })
        .max()
        .unwrap_or(0)
}

/// Pair emit: H-progress for both ends under the best consistent placement.
pub fn pair_site_h_progress(
    pair: &PairCandidate,
    diff: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> i32 {
    let atoms = pair
        .end_atoms()
        .map(|(a, b)| [a, b])
        .unwrap_or([pair.site, pair.site]);
    let (p0, p1) = pair.path_ends();
    site_h_progress_best_placement(&pair.effect, &atoms, &[p0, p1], diff, mol, target)
}

/// Keep predicate for [`crate::find_path::find_path_with`] from an [`AtomDiff`].
pub fn keep_against_diff(diff: &AtomDiff) -> impl Fn(&Candidate) -> bool + '_ {
    move |c: &Candidate| candidate_could_help(c, diff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use crate::rules::{dealkylation, hydroxylation};

    #[test]
    fn ethane_to_ethanol_needs_oxygen_on_carbon() {
        let reactant = parse_mol("CC").unwrap();
        let target = parse_mol("CCO").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(any_needs_oxygen(&target, &diff), "{diff:?}");
        assert!(!diff.has_cleavage());
        let set = hydroxylation();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(cands.iter().any(|c| candidate_could_help_on(
            c,
            &diff,
            Some(&reactant),
            Some(&target)
        )));
    }

    #[test]
    fn anisole_to_phenol_is_cleavage() {
        let reactant = parse_mol("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.has_cleavage() || diff.target_smaller(), "{diff:?}");
        let set = dealkylation();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(
            cands.iter().any(|c| c.pattern.effect.cleaves
                && candidate_could_help_on(c, &diff, Some(&reactant), Some(&target))),
            "dealkylation should survive filter; diff={diff:?}"
        );
    }

    #[test]
    fn meoph_oh_mcs_covers_ring() {
        let reactant = parse_mol("COc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=C(O)C(=O)C(O)=C1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(
            diff.mapping.len() >= 6,
            "BondCompare::Any MCS should cover the ring; got {}",
            diff.mapping.len()
        );
        assert!(
            !diff.loses_aromaticity.is_empty() || any_needs_oxygen(&target, &diff),
            "{diff:?}"
        );
    }

    #[test]
    fn meoph_oh_pair_end_filter_matches_python() {
        use crate::rules::{phase_one, quinone_formation};
        let reactant = parse_mol("COc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=C(O)C(=O)C(O)=C1").unwrap();
        let diff = atom_diff(&reactant, &target);
        let qf = quinone_formation();
        let pairs = qf
            .pair_candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let kept: Vec<_> = pairs
            .iter()
            .filter(|p| pair_could_help(p, &diff, &reactant, &target))
            .collect();
        // Python QuinoneFormation filter_sites keeps 3 pair sites here.
        assert_eq!(
            kept.len(),
            3,
            "pair_kept={} of {}; names={:?}",
            kept.len(),
            pairs.len(),
            kept.iter().map(|p| &p.pattern_name).collect::<Vec<_>>()
        );
        assert!(
            kept.iter()
                .all(|p| !p.left.effect.methide && !p.right.effect.methide),
            "methide ends must fail alkyl_bond_raises on MeOPhOH"
        );
        let _ = phase_one();
    }

    #[test]
    fn dimethoxy_keeps_dealkylation() {
        let reactant = parse_mol("COc1ccc(CCN)cc1OC").unwrap();
        let target = parse_mol("NCCc1ccc(O)c(O)c1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.target_smaller());
        assert!(diff.has_cleavage(), "{diff:?}");
        let set = dealkylation();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let kept: Vec<_> = cands
            .iter()
            .filter(|c| candidate_could_help_on(c, &diff, Some(&reactant), Some(&target)))
            .collect();
        assert!(
            !kept.is_empty(),
            "expected dealkylation survivors; cands={} diff={diff:?}",
            cands.len()
        );
    }

    #[test]
    fn hydroquinone_to_quinone_loses_aromaticity_and_h() {
        let reactant = parse_mol("Oc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(
            !diff.loses_aromaticity.is_empty() || diff.h_loss(),
            "{diff:?}"
        );
    }

    #[test]
    fn hydrogenation_refused_toward_quinone() {
        use crate::rules::hydrogenation;
        let reactant = parse_mol("c1ccccc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let diff = atom_diff(&reactant, &target);
        let set = hydrogenation();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        for c in &cands {
            if c.pattern.effect.adds.as_deref() == Some("HH") && !diff.h_gain() {
                assert!(
                    !pattern_could_help(&c.pattern.effect, &diff) || c.pattern.effect.cleaves,
                    "pattern {} should not help toward quinone",
                    c.pattern.name
                );
            }
        }
    }

    #[test]
    fn target_smaller_refuses_non_cleaving() {
        let reactant = parse_mol("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.target_smaller());
        let hydroxyl = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            ..Effect::default()
        };
        assert!(!pattern_could_help(&hydroxyl, &diff));
        let cleave = Effect {
            adds: Some("O".into()),
            cleaves: true,
            ..Effect::default()
        };
        assert!(pattern_could_help(&cleave, &diff));
    }

    #[test]
    fn oxygen_reduction_carbonyl_partners_allow_toward_alcohol() {
        // OR sites the heteroatom (map 1); H change is also on the partner
        // carbon. Partners must be in scope so adds-H / progress see the full
        // helpful delta (apply helps; undo toward the carbonyl would not).
        use crate::rules::oxygen_reduction;
        let reactant = parse_mol("CC=O").unwrap();
        let target = parse_mol("CCO").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.h_gain(), "{diff:?}");
        let set = oxygen_reduction();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let carbonyl: Vec<_> = cands
            .iter()
            .filter(|c| c.pattern.name == "carbonyl")
            .collect();
        assert!(!carbonyl.is_empty(), "expected carbonyl OR candidates");
        assert!(
            carbonyl.iter().any(|c| candidate_could_help_on(
                c,
                &diff,
                Some(&reactant),
                Some(&target)
            )),
            "OR carbonyl should help CC=O→CCO when partners expand scope"
        );
    }

    #[test]
    fn oxygen_reduction_refused_when_undoing_toward_carbonyl() {
        // Still a carbonyl match, but the target wants fewer H on that carbon
        // (oxidation to acid) — adding H would undo progress.
        use crate::rules::oxygen_reduction;
        let reactant = parse_mol("CC=O").unwrap();
        let target = parse_mol("CC(=O)O").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.h_loss() || !diff.h_gain(), "{diff:?}");
        let set = oxygen_reduction();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let carbonyl: Vec<_> = cands
            .iter()
            .filter(|c| c.pattern.name == "carbonyl")
            .collect();
        assert!(
            !carbonyl.is_empty(),
            "need a carbonyl that OR can still match while undoing"
        );
        for c in &carbonyl {
            assert!(
                !candidate_could_help_on(c, &diff, Some(&reactant), Some(&target)),
                "OR toward acid (H loss at carbonyl C) must not help"
            );
        }
    }

    #[test]
    fn order_key_h_progress_prefers_adds_h_toward_alcohol() {
        use crate::rules::oxygen_reduction;
        let reactant = parse_mol("CC=O").unwrap();
        let target = parse_mol("CCO").unwrap();
        let diff = atom_diff(&reactant, &target);
        let set = oxygen_reduction();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let c = cands
            .iter()
            .find(|c| c.pattern.name == "carbonyl")
            .expect("carbonyl");
        let (_a, _b, _c, prog_no_mol, _) = candidate_order_key_on(c, &diff, None, None);
        let (_a, _b, _c, prog_with, _) =
            candidate_order_key_on(c, &diff, Some(&reactant), Some(&target));
        // Negated progress: with partners, progress > 0 ⇒ key more negative.
        assert!(
            prog_with < prog_no_mol,
            "partners should raise H-progress (lower negated key); with={prog_with} without={prog_no_mol}"
        );
    }

    #[test]
    fn tag_lift_dh_matches_full_mcs_cost() {
        use crate::forest_mol::ForestMol;
        use crate::rules::dehydrogenation;

        let parent = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let parent_diff = atom_diff(parent.mol(), &target);
        let pairs = dehydrogenation()
            .pair_candidates_leaf(parent.mol())
            .unwrap();
        assert!(!pairs.is_empty());
        let pieces = pairs[0].materialize_mols(parent.mol()).unwrap();
        assert!(!pieces.is_empty());
        let child = parent.adopt_product(pieces[0].clone());
        assert!(
            added_heavy_atoms(&parent, &child).is_empty(),
            "DH does not add heavy atoms"
        );
        let lifted = try_atom_diff_for_child(&parent, &parent_diff, &child, &target)
            .expect("same-atom-count tag lift should not need MCS");
        let full = atom_diff(child.mol(), &target);
        assert_eq!(
            lifted.cost(),
            full.cost(),
            "lifted={lifted:?} full={full:?}"
        );
    }

    #[test]
    fn tag_lift_extends_added_oxygen() {
        use crate::forest_mol::ForestMol;
        use crate::hydroxylation::hydroxylation;

        let parent = ForestMol::parse("CC").unwrap();
        let target = parse_mol("CCO").unwrap();
        let parent_diff = atom_diff(parent.mol(), &target);
        let cands = hydroxylation()
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let pieces = cands[0].materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        assert_eq!(child.mol().atom_count(), 3);
        let added = added_heavy_atoms(&parent, &child);
        assert_eq!(added.len(), 1);
        assert_eq!(
            child.mol().atom(atom_idx(added[0])).element.atomic_number(),
            8
        );
        // Add: lift + extend where diff shows a free target O; no MCS.
        let lifted = try_atom_diff_for_child(&parent, &parent_diff, &child, &target)
            .expect("generator lift + extend should place added oxygen");
        assert_eq!(lifted.cost(), 0, "{lifted:?}");
        let via = atom_diff_for_child(&parent, &parent_diff, &child, &target);
        assert_eq!(via.cost(), 0, "{via:?}");
        assert!(child.shares_tag_gen(&parent));
    }

    #[test]
    fn site_cast_residual_is_lower_bound_for_hydroxylation() {
        use crate::forest_mol::ForestMol;
        use crate::hydroxylation::hydroxylation;

        let parent = ForestMol::parse("CC").unwrap();
        let target = parse_mol("CCO").unwrap();
        let parent_diff = atom_diff(parent.mol(), &target);
        assert!(parent_diff.cost() > 0, "{parent_diff:?}");
        let cands = hydroxylation()
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let c = &cands[0];
        let atoms = [c.site];
        let goal = residual_cost_after_site_cast(&parent_diff, &c.pattern.effect, &atoms, &[]);
        // Casting hydroxyl onto a needs-oxygen site should claim the O gap.
        assert!(
            goal < parent_diff.cost(),
            "goal={goal} parent={}",
            parent_diff.cost()
        );
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let lifted =
            try_atom_diff_for_child_goal(&parent, &parent_diff, &child, &target, Some(goal))
                .expect("lift");
        assert!(
            lifted.cost() <= goal,
            "lifted={} goal={goal} parent_diff={parent_diff:?}",
            lifted.cost()
        );
        assert_eq!(lifted.cost(), 0);
    }

    #[test]
    fn site_cast_residual_cleavage_credits_leave() {
        let reactant = parse_mol("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let diff = atom_diff(&reactant, &target);
        let effect = Effect {
            adds: Some("O".into()),
            cleaves: true,
            leave_count: Some(1),
            leave_formula: crate::pattern::leave_me(),
            delta_formula: crate::pattern::compose_delta_formula(
                Some("O"),
                None,
                &crate::pattern::leave_me(),
            ),
            ..Effect::default()
        };
        // Site = cleaved bond atoms (O–Me): use any cleavage bond from the diff.
        let (a, b) = *diff
            .cleavage_bonds
            .iter()
            .next()
            .expect("anisole→phenol cleavage bond");
        let goal = residual_cost_after_site_cast(&diff, &effect, &[a, b], &[]);
        assert!(
            goal < diff.cost(),
            "cleavage cast should drop leave/bond cost; goal={goal} parent={}",
            diff.cost()
        );
    }

    #[test]
    fn tag_lift_shrinks_on_dealkylation() {
        use crate::forest_mol::ForestMol;

        // Anisole → phenol: methyl carbon removed, no heavy add.
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let parent_diff = atom_diff(parent.mol(), &target);
        let parent_cost = parent_diff.cost();
        let cands = dealkylation()
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!cands.is_empty());
        let phenol = canon_of("Oc1ccccc1").unwrap();
        let mut child = None;
        for c in &cands {
            let pieces = c.materialize_mols(parent.mol()).unwrap();
            for piece in pieces {
                let adopted = parent.adopt_product(piece);
                if adopted.csmi().as_ref() == phenol.as_str() {
                    child = Some(adopted);
                    break;
                }
            }
            if child.is_some() {
                break;
            }
        }
        let child = child.expect("phenol from anisole dealkylation");
        assert!(
            added_heavy_atoms(&parent, &child).is_empty(),
            "dealkylation keeps no new heavies on the kept fragment"
        );
        assert!(
            try_atom_diff_for_child(&parent, &parent_diff, &child, &target).is_none(),
            "closer still refuses shrink lift"
        );
        let lifted = try_lift_cleaved_child(&parent, &parent_diff, &child, &target)
            .expect("cleavage shrink should tag-lift");
        let via = atom_diff_after_cleavage(&parent, &parent_diff, &child, &target);
        assert!(
            via.cost() <= parent_cost,
            "toward phenol: parent={parent_cost} child={}",
            via.cost()
        );
        // When lift already drops cost, after_cleavage should not need MCS.
        if lifted.cost() < parent_cost {
            assert_eq!(via.cost(), lifted.cost());
        }
    }

    #[test]
    fn dh_product_ends_match_after_hydroquinone_dh() {
        use crate::forest_mol::ForestMol;
        use crate::rules::dehydrogenation;

        let parent = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let pairs = dehydrogenation()
            .pair_candidates_leaf(parent.mol())
            .unwrap();
        assert!(!pairs.is_empty());
        let pair = &pairs[0];
        let (end_a, end_b) = pair.end_atoms().expect("DH pair ends");
        // Pre-application reactant ends also match here — but the gate is
        // defined on the product.
        let pieces = pair.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let tag_a = parent.tag_of(end_a).expect("tagged");
        let tag_b = parent.tag_of(end_b).expect("tagged");
        let product_ends = [
            child.index_of(tag_a).expect("product end a"),
            child.index_of(tag_b).expect("product end b"),
        ];
        assert!(dh_product_ends_match(child.mol(), &product_ends, &target));
    }

    #[test]
    fn dh_product_ends_refuse_wrong_connectivity() {
        // Identity "product" still carrying methoxy vs quinone target.
        let product = parse_mol("COc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let oxygens: Vec<usize> = (0..product.atom_count())
            .filter(|&i| product.atom(atom_idx(i)).element.atomic_number() == 8)
            .collect();
        assert_eq!(oxygens.len(), 2);
        assert!(!dh_product_ends_match(&product, &oxygens, &target));
    }
}
