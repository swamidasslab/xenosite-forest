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
    pub needs_oxygen: BTreeSet<usize>,
    pub needs_carbonyl: BTreeSet<usize>,
    pub needs_alcohol: BTreeSet<usize>,
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
        let h_off = self.h_delta.values().filter(|&&d| d != 0).count();
        3 * self.cleaved.len()
            + 3 * self.n_extra
            + 2 * self.needs_oxygen.len()
            + self.loses_aromaticity.len()
            + h_off
            + 3 * self.cleavage_bonds.len()
            + self.bond_order_mismatches
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

    let mut needs_oxygen = BTreeSet::new();
    let mut needs_carbonyl = BTreeSet::new();
    let mut needs_alcohol = BTreeSet::new();
    for (t_idx, atom) in target.atoms() {
        let t = atom_usize(t_idx);
        if atom.element.atomic_number() != 8 || image.contains(&t) {
            continue;
        }
        for (nbr, bond_idx) in target.neighbors(t_idx) {
            let n = atom_usize(nbr);
            let order = order_value(target.bond(bond_idx).order);
            for (&r_idx, &t_mapped) in mapping {
                if t_mapped != n {
                    continue;
                }
                needs_oxygen.insert(r_idx);
                if order >= 1.5 {
                    needs_carbonyl.insert(r_idx);
                } else {
                    needs_alcohol.insert(r_idx);
                }
            }
        }
    }

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

    let needs_oxygen = expand_by_rank(reactant, needs_oxygen);
    let needs_carbonyl = expand_by_rank(reactant, needs_carbonyl);
    let needs_alcohol = expand_by_rank(reactant, needs_alcohol);
    let loses_aromaticity = expand_by_rank(reactant, loses_aromaticity);
    let h_delta = expand_h_delta(reactant, h_delta);

    let mut diff = AtomDiff {
        mapping: mapping.clone(),
        mappings: vec![mapping.clone()],
        needs_oxygen,
        needs_carbonyl,
        needs_alcohol,
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
        primary.needs_oxygen.extend(&view.needs_oxygen);
        primary.needs_carbonyl.extend(&view.needs_carbonyl);
        primary.needs_alcohol.extend(&view.needs_alcohol);
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
/// unmapped — callers may [`extend_mapping_for_added`] or fall back to full MCS.
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

/// Child [`AtomDiff`] via tag-lifted parent MCS when possible.
///
/// **Safe lift only:** parent and child share the same heavy-atom tag set
/// (no add/remove — typically DH / bond-order edits). Add or remove can
/// produce a cheaper/dearer cost than true MCS and poison the closer; those
/// return `None` so the caller runs full [`atom_diff`].
///
/// Cleavage shrinks: use [`try_lift_cleaved_child`] / [`atom_diff_after_cleavage`]
/// (lift is allowed there; closer still refuses shrink via this function).
///
/// Never runs MCS itself.
pub fn try_atom_diff_for_child(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> Option<AtomDiff> {
    if !added_heavy_atoms(parent, child).is_empty() {
        return None;
    }
    // Any heavy parent atom missing on the child → shrink; lift cost unreliable
    // for the closer. Cleavage expand uses [`try_lift_cleaved_child`] instead.
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
    Some(atom_diff_from_mappings(child.mol(), target, lifted))
}

/// Tag-lift after a cleavage shrink (no heavy adds). Removed atoms drop out of
/// the parent MCS map via [`lift_mappings`].
///
/// Lifted cost may overshoot true MCS; callers that need a hard closer bound
/// should fall back with [`atom_diff_after_cleavage`].
pub fn try_lift_cleaved_child(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> Option<AtomDiff> {
    if !added_heavy_atoms(parent, child).is_empty() {
        return None;
    }
    let parent_maps = if parent_diff.mappings.is_empty() {
        vec![parent_diff.mapping.clone()]
    } else {
        parent_diff.mappings.clone()
    };
    let lifted = lift_mappings(parent, child, &parent_maps)?;
    Some(atom_diff_from_mappings(child.mol(), target, lifted))
}

/// Child diff after cleavage: prefer tag-lift when it already shows a strict
/// cost drop vs the parent; otherwise full MCS.
pub fn atom_diff_after_cleavage(
    parent: &crate::forest_mol::ForestMol,
    parent_diff: &AtomDiff,
    child: &crate::forest_mol::ForestMol,
    target: &Molecule,
) -> AtomDiff {
    if let Some(lifted) = try_lift_cleaved_child(parent, parent_diff, child, target) {
        if lifted.cost() < parent_diff.cost() {
            return lifted;
        }
    }
    atom_diff(child.mol(), target)
}

/// Child [`AtomDiff`] via tag-lift when possible; else full MCS.
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

/// Pattern-level gate (Python `_pattern_could_help`).
pub fn pattern_could_help(effect: &Effect, diff: &AtomDiff) -> bool {
    pattern_could_help_on(effect, diff, None)
}

/// Same as [`pattern_could_help`], with optional live mol for formula-O gate.
pub fn pattern_could_help_on(effect: &Effect, diff: &AtomDiff, mol: Option<&Molecule>) -> bool {
    let can_cleave = effect.cleaves;
    if diff.target_smaller() && !can_cleave {
        return false;
    }
    if effect_adds_oxygen(effect) && !can_cleave && diff.needs_oxygen.is_empty() {
        return false;
    }
    if effect_adds_oxygen(effect) && !can_cleave && !effect.dearomatizes {
        if let Some(m) = mol {
            // Target O count is not on AtomDiff; approximate via needs + extras.
            // Prefer live formula vs target heavy oxygen when available.
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

/// Live-mol oxygen gate (Python formula check).
pub fn pattern_could_help_mol(
    effect: &Effect,
    diff: &AtomDiff,
    mol: &Molecule,
    target: &Molecule,
) -> bool {
    if !pattern_could_help_on(effect, diff, Some(mol)) {
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

/// Site-level gate (Python `_site_could_help`) for a deferred candidate.
pub fn candidate_could_help(candidate: &Candidate, diff: &AtomDiff) -> bool {
    candidate_could_help_on(candidate, diff, None, None)
}

fn scope_could_help(
    effect: &Effect,
    atoms: &[usize],
    path_ends: &[usize],
    diff: &AtomDiff,
) -> bool {
    let mut scope: HashSet<usize> = atoms.iter().copied().collect();
    scope.extend(path_ends.iter().copied());
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

/// Full site gate with live mol (leave_count / methide partner).
pub fn candidate_could_help_on(
    candidate: &Candidate,
    diff: &AtomDiff,
    mol: Option<&Molecule>,
    target: Option<&Molecule>,
) -> bool {
    let effect = &candidate.pattern.effect;
    let ok = match (mol, target) {
        (Some(m), Some(t)) => pattern_could_help_mol(effect, diff, m, t),
        _ => pattern_could_help(effect, diff),
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

    if effect_adds_oxygen(effect) && !effect.dearomatizes {
        let oxygen_sites: Vec<_> = atoms
            .iter()
            .copied()
            .filter(|a| diff.needs_oxygen.contains(a))
            .collect();
        if oxygen_sites.is_empty() {
            return false;
        }
        if oxygen_sites
            .iter()
            .all(|a| diff.needs_carbonyl.contains(a) && diff.loses_aromaticity.contains(a))
        {
            return false;
        }
    }

    // Single-site methide: alkyl partner needs an exocyclic C–C bond raise.
    // Pair ends use [`pair_could_help`] (per-end partner), not this merge.
    if effect.partner.as_deref() == Some("C") {
        if let Some(m) = mol {
            if !atoms.iter().any(|&a| alkyl_bond_raises(m, a, diff)) {
                return false;
            }
        }
    }

    scope_could_help(effect, &atoms, &[], diff)
}

/// Pair-site gate (Python `_site_could_help` when ``"ends"`` is on the info).
///
/// Per-end oxygen and ``partner == "C"`` (methide) read each end's effect, not
/// the merged span. Dearomatize / H gates use end atoms ∪ path anchors.
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

    let ends = [(&pair.left.effect, end_a), (&pair.right.effect, end_b)];
    for (end, atom) in ends {
        if effect_adds_oxygen(end) && !diff.needs_oxygen.contains(&atom) {
            return false;
        }
        if end.partner.as_deref() == Some("C") && !alkyl_bond_raises(mol, atom, diff) {
            return false;
        }
    }

    let (p0, p1) = pair.path_ends();
    scope_could_help(effect, &atoms, &[p0, p1], diff)
}

/// Sort key for expand: cleavage / dearom / oxygen first (Python `order_key`).
pub fn candidate_order_key(candidate: &Candidate, diff: &AtomDiff) -> (u8, u8, u8, String) {
    let effect = &candidate.pattern.effect;
    let want_cleave = diff.target_smaller() || diff.has_cleavage();
    let want_dear = !diff.loses_aromaticity.is_empty();
    let want_oxy = !diff.needs_oxygen.is_empty();
    let cleave = if effect.cleaves { 0 } else { 1 };
    let dear = if effect.dearomatizes { 0 } else { 1 };
    let oxy = if effect_adds_oxygen(effect) { 0 } else { 1 };
    let primary = if want_cleave { cleave } else { 0 };
    let secondary = if want_dear { dear } else { 0 };
    let tertiary = if want_oxy { oxy } else { 0 };
    (primary, secondary, tertiary, candidate.pattern.name.clone())
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
        assert!(!diff.needs_oxygen.is_empty(), "{diff:?}");
        assert!(!diff.has_cleavage());
        let set = hydroxylation();
        let cands = set
            .candidates(&reactant)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(cands.iter().any(|c| candidate_could_help(c, &diff)));
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
            !diff.loses_aromaticity.is_empty() || !diff.needs_oxygen.is_empty(),
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
        // Add/remove lifts are unsafe for closer cost → try is None; MCS fallback.
        assert!(
            try_atom_diff_for_child(&parent, &parent_diff, &child, &target).is_none(),
            "add should fall back to full MCS"
        );
        let via = atom_diff_for_child(&parent, &parent_diff, &child, &target);
        assert_eq!(via.cost(), 0, "{via:?}");
        assert!(child.shares_tag_gen(&parent));
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
}
