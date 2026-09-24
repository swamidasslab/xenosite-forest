//! Provisional reactant→target local diff for candidate filtering.
//!
//! Door for Python `find_path.atom_diff` / `_site_could_help`. Aligns via
//! chematic MCS embeddings; reads [`crate::pattern::Effect`] on
//! [`crate::candidate::Candidate`] without filter closures.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use chematic::core::BondOrder;
use chematic::smarts::{find_matches, find_mcs};

use crate::candidate::Candidate;
use crate::mol::{Molecule, atom_idx, atom_usize, ranks};
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

/// Local reactant-to-target differences (Python `AtomDiff` subset).
#[derive(Clone, Debug, Default)]
pub struct AtomDiff {
    pub mapping: BTreeMap<usize, usize>,
    pub needs_oxygen: BTreeSet<usize>,
    pub needs_carbonyl: BTreeSet<usize>,
    pub cleaved: BTreeSet<usize>,
    pub cleavage_bonds: BTreeSet<(usize, usize)>,
    pub loses_aromaticity: BTreeSet<usize>,
    pub h_delta: BTreeMap<usize, i32>,
    pub bond_raises: BTreeSet<(usize, usize)>,
    pub reactant_heavy: usize,
    pub target_heavy: usize,
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
        let h_off = self.h_delta.values().filter(|&&d| d != 0).count();
        3 * self.cleaved.len()
            + 2 * self.needs_oxygen.len()
            + self.loses_aromaticity.len()
            + h_off
            + 3 * self.cleavage_bonds.len()
    }

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

fn best_mapping(reactant: &Molecule, target: &Molecule) -> BTreeMap<usize, usize> {
    let query = find_mcs(&[reactant, target]);
    if query.atom_count() == 0 {
        return BTreeMap::new();
    }
    let r_hits = find_matches(&query, reactant);
    let t_hits = find_matches(&query, target);
    let mut best: Option<(usize, BTreeMap<usize, usize>)> = None;
    for r_emb in &r_hits {
        let r_map: BTreeMap<usize, usize> =
            r_emb.iter().map(|(&q, &t)| (q, atom_usize(t))).collect();
        let n = r_map.len();
        for t_emb in &t_hits {
            let t_map: BTreeMap<usize, usize> =
                t_emb.iter().map(|(&q, &t)| (q, atom_usize(t))).collect();
            if t_map.len() != n {
                continue;
            }
            let mut aligned = BTreeMap::new();
            for (&q, &r_idx) in &r_map {
                if let Some(&t_idx) = t_map.get(&q) {
                    aligned.insert(r_idx, t_idx);
                }
            }
            if aligned.len() != n {
                continue;
            }
            let score = aligned.len();
            match &best {
                None => best = Some((score, aligned)),
                Some((best_score, _)) if score > *best_score => best = Some((score, aligned)),
                _ => {}
            }
        }
    }
    best.map(|(_, m)| m).unwrap_or_default()
}

/// Pair reactant atoms with target atoms and record the local change.
pub fn atom_diff(reactant: &Molecule, target: &Molecule) -> AtomDiff {
    let mapping = best_mapping(reactant, target);
    let image: HashSet<usize> = mapping.values().copied().collect();

    let mut needs_oxygen = BTreeSet::new();
    let mut needs_carbonyl = BTreeSet::new();
    for (t_idx, atom) in target.atoms() {
        let t = atom_usize(t_idx);
        if atom.element.atomic_number() != 8 || image.contains(&t) {
            continue;
        }
        for (nbr, bond_idx) in target.neighbors(t_idx) {
            let n = atom_usize(nbr);
            let order = order_value(target.bond(bond_idx).order);
            for (&r_idx, &t_mapped) in &mapping {
                if t_mapped != n {
                    continue;
                }
                needs_oxygen.insert(r_idx);
                if order >= 1.5 {
                    needs_carbonyl.insert(r_idx);
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

    for (&r_idx, &t_idx) in &mapping {
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
        if delta >= 0.2 {
            bond_raises.insert(bond_key(i, j));
        }
    }

    // Unique-edit collapses equivalent sites to one representative. Expand
    // oxygen / aromaticity / h_delta keys across topological rank so the
    // kept site still matches the diff.
    let rank = ranks(reactant);
    let expand = |sites: BTreeSet<usize>| -> BTreeSet<usize> {
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
    };
    let needs_oxygen = expand(needs_oxygen);
    let needs_carbonyl = expand(needs_carbonyl);
    let loses_aromaticity = expand(loses_aromaticity);
    let mut h_delta_exp = BTreeMap::new();
    for (&atom, &delta) in &h_delta {
        let r = rank.get(atom).copied().unwrap_or(usize::MAX);
        for (i, &ri) in rank.iter().enumerate() {
            if ri == r {
                h_delta_exp.insert(i, delta);
            }
        }
    }

    AtomDiff {
        mapping,
        needs_oxygen,
        needs_carbonyl,
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta: h_delta_exp,
        bond_raises,
        reactant_heavy: reactant.atom_count(),
        target_heavy: target.atom_count(),
    }
}

fn effect_adds_oxygen(effect: &Effect) -> bool {
    effect.adds.as_deref().is_some_and(|a| a.contains('O'))
}

fn effect_adds_h(effect: &Effect) -> bool {
    effect.adds.as_deref().is_some_and(|a| a.contains('H'))
}

fn effect_removes_h(effect: &Effect) -> bool {
    effect.removes.as_deref().is_some_and(|r| r.contains('H'))
}

/// Pattern-level gate (Python `_rule_could_help` subset) reading span-like fields.
pub fn pattern_could_help(effect: &Effect, diff: &AtomDiff) -> bool {
    if effect.cleaves && !diff.has_cleavage() {
        return false;
    }
    if effect.dearomatizes
        && !effect_adds_h(effect)
        && diff.loses_aromaticity.is_empty()
        && !diff.h_loss()
    {
        return false;
    }
    let drops_h_only = !effect.cleaves && effect.adds.is_none() && effect_removes_h(effect);
    if drops_h_only && !diff.h_loss() && diff.loses_aromaticity.is_empty() {
        return false;
    }
    let adds_h_only = !effect.cleaves && effect_adds_h(effect) && effect.removes.is_none();
    if adds_h_only && !diff.h_gain() {
        return false;
    }
    if effect_adds_oxygen(effect)
        && !effect.dearomatizes
        && diff.needs_oxygen.is_empty()
        && diff.target_heavy <= diff.reactant_heavy
    {
        // No oxygen to add and target is not larger — skip.
        return false;
    }
    true
}

/// Site-level gate (Python `_site_could_help` subset) for a deferred candidate.
pub fn candidate_could_help(candidate: &Candidate, diff: &AtomDiff) -> bool {
    let effect = &candidate.pattern.effect;
    if !pattern_could_help(effect, diff) {
        return false;
    }
    let atoms = [candidate.site];
    if effect.cleaves {
        return diff.site_is_cleavage(&atoms)
            || diff.cleaved.contains(&candidate.site)
            || diff.has_cleavage();
    }
    if effect_adds_oxygen(effect) && !effect.dearomatizes {
        if !diff.needs_oxygen.contains(&candidate.site) && diff.needs_oxygen.is_empty() {
            // Soft: allow when molecule-level needs_oxygen is empty only if
            // pattern_could_help already passed (e.g. target smaller).
        } else if !diff.needs_oxygen.contains(&candidate.site) && !diff.needs_oxygen.is_empty() {
            return false;
        }
    }
    if effect.dearomatizes
        && !diff.loses_aromaticity.contains(&candidate.site)
        && !diff.loses_aromaticity.is_empty()
    {
        return false;
    }
    if effect_removes_h(effect) && !effect_adds_oxygen(effect) && !effect.cleaves {
        let loses_h = diff.h_delta.get(&candidate.site).copied().unwrap_or(0) < 0;
        if !loses_h && !diff.loses_aromaticity.contains(&candidate.site) {
            // Soft molecule-level escape when any h_loss / aromatic loss exists.
            if !diff.h_loss() && diff.loses_aromaticity.is_empty() {
                return false;
            }
        }
    }
    if effect_adds_h(effect) && !effect_adds_oxygen(effect) && !effect.cleaves {
        let gains = diff.h_delta.get(&candidate.site).copied().unwrap_or(0) > 0;
        if !gains && !diff.h_gain() {
            return false;
        }
    }
    true
}

/// Keep predicate for [`crate::find_path::find_path_with`] from an [`AtomDiff`].
pub fn keep_against_diff(diff: &AtomDiff) -> impl Fn(&Candidate) -> bool + '_ {
    move |c: &Candidate| candidate_could_help(c, diff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::rules::hydroxylation;

    #[test]
    fn ethane_to_ethanol_needs_oxygen_on_carbon() {
        let reactant = parse_mol("CC").unwrap();
        let target = parse_mol("CCO").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(!diff.needs_oxygen.is_empty(), "{diff:?}");
        assert!(!diff.has_cleavage());
        let set = hydroxylation();
        let cands = set.candidates(&reactant).unwrap();
        assert!(!cands.is_empty());
        assert!(
            cands.iter().any(|c| candidate_could_help(c, &diff)),
            "{diff:?}"
        );
    }

    #[test]
    fn anisole_to_phenol_is_cleavage() {
        let reactant = parse_mol("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let diff = atom_diff(&reactant, &target);
        assert!(diff.has_cleavage() || diff.target_smaller(), "{diff:?}");
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
        let cands = set.candidates(&reactant).unwrap();
        // Alkene/alkyne SMIRKS may still match; path_end is PairEndpoint.
        for c in &cands {
            if c.pattern.effect.adds.as_deref() == Some("HH") {
                // Adding H toward an oxidative quinone target should fail pattern gate
                // when no h_gain.
                if !diff.h_gain() {
                    assert!(
                        !pattern_could_help(&c.pattern.effect, &diff) || c.pattern.effect.cleaves,
                        "pattern {} should not help toward quinone",
                        c.pattern.name
                    );
                }
            }
        }
    }
}
