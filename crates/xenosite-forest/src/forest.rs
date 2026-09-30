//! Formula and structure-answer bags held by [`crate::forest_mol::ForestMol`].
//!
//! Not a Python `_forest` dict. [`ForestMol`] owns these as fields.

use std::collections::BTreeMap;
use std::rc::Rc;

use crate::mol::Molecule;

/// Charge-unit keys in [`Formula::counts`]: magnitudes only (never zero, never
/// negative). Net formal charge is still [`Formula::charge`] (= `+` − `-`).
pub const CHARGE_PLUS: &str = "+";
pub const CHARGE_MINUS: &str = "-";

/// Element → count (hydrogens included) plus formal charge.
///
/// Cached on each [`crate::forest_mol::ForestMol`]. ``counts`` is the map
/// filters and pattern ``delta_formula`` compare against. Formal charge units
/// live in ``counts`` under [`CHARGE_PLUS`] / [`CHARGE_MINUS`] (strictly
/// positive when present); [`Self::charge`] is the signed net.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Formula {
    /// Map of element symbol → count, including hydrogens, plus ``"+"`` /
    /// ``"-"`` formal-charge unit counts when nonzero.
    pub counts: BTreeMap<String, i32>,
    pub charge: i32,
}

impl Formula {
    /// Insert or clear ``"+"`` / ``"-"`` from absolute positive / negative
    /// formal-charge unit totals and set [`Self::charge`] = pos − neg.
    pub fn with_charge_units(mut self, positive: i32, negative: i32) -> Self {
        debug_assert!(positive >= 0 && negative >= 0);
        if positive > 0 {
            self.counts.insert(CHARGE_PLUS.into(), positive);
        } else {
            self.counts.remove(CHARGE_PLUS);
        }
        if negative > 0 {
            self.counts.insert(CHARGE_MINUS.into(), negative);
        } else {
            self.counts.remove(CHARGE_MINUS);
        }
        self.charge = positive - negative;
        self
    }
}

/// True for formula keys that are formal-charge units, not elements.
pub fn is_charge_key(key: &str) -> bool {
    key == CHARGE_PLUS || key == CHARGE_MINUS
}

/// Structure-dependent answers, filled on first read. Shared by `Rc` when
/// the chemistry is the same; a new bag after an edit.
#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub csmi: Option<Rc<str>>,
    /// Fail-closed dedup key ([`crate::mol::stable_csmi_key`]).
    /// Outer `None` = not filled; inner `None` = proven unstable (do not index).
    pub stable_csmi: Option<Option<Rc<str>>>,
    pub formula: Option<Rc<Formula>>,
    pub topol_equiv: Option<Rc<Vec<usize>>>,
    /// Atom+bond automorphism generators (canonaut). Needed for site and
    /// higher-order (plan) orbits; filled once per structure bag.
    pub atom_bond_generators: Option<Rc<Vec<crate::orbits::AtomBondGenerator>>>,
    pub smarts_matches: BTreeMap<String, Rc<Vec<BTreeMap<u16, usize>>>>,
}

/// Heavy-atom counts, explicit + implicit hydrogens, and formal charge.
///
/// Positive formal charges sum under ``"+"``, absolute negative under ``"-"``
/// (strictly positive counts). [`Formula::charge`] is the signed net.
pub fn molecule_formula(mol: &Molecule) -> Formula {
    let mut counts = BTreeMap::new();
    let mut positive = 0i32;
    let mut negative = 0i32;
    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        let c = i32::from(atom.charge);
        if c > 0 {
            positive += c;
        } else if c < 0 {
            negative += -c;
        }
        if z == 1 {
            *counts.entry("H".to_string()).or_insert(0) += 1;
            continue;
        }
        *counts.entry(atom.element.symbol().to_string()).or_insert(0) += 1;
        let hydrogens = mol.implicit_hydrogen_count(idx) as i32;
        if hydrogens > 0 {
            *counts.entry("H".to_string()).or_insert(0) += hydrogens;
        }
    }
    Formula { counts, charge: 0 }.with_charge_units(positive, negative)
}

/// Change in element counts (and charge) from ``before`` to ``after``.
///
/// Zero-count keys are omitted from ``counts``. Charge units appear as
/// ``"+"`` / ``"-"`` deltas (may be negative when units are lost).
/// [`Formula::charge`] is the signed net difference.
pub fn formula_delta(before: &Formula, after: &Formula) -> Formula {
    let mut counts = BTreeMap::new();
    let mut keys: Vec<&str> = before
        .counts
        .keys()
        .chain(after.counts.keys())
        .map(String::as_str)
        .collect();
    keys.sort_unstable();
    keys.dedup();
    for key in keys {
        let delta = after.counts.get(key).copied().unwrap_or(0)
            - before.counts.get(key).copied().unwrap_or(0);
        if delta != 0 {
            counts.insert(key.to_string(), delta);
        }
    }
    Formula {
        counts,
        charge: after.charge - before.charge,
    }
}

/// Element-count L1 distance between two formulas (**includes H and
/// formal-charge units** ``"+"`` / ``"-"``).
///
/// Used by [`crate::find_path::MatchScoreSpec`] closeness / improvement scoring
/// and by expand formula distance (same door as atom_diff residual scoring).
pub fn formula_l1(a: &Formula, b: &Formula) -> usize {
    let mut keys: Vec<&str> = a
        .counts
        .keys()
        .chain(b.counts.keys())
        .map(String::as_str)
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let mut dist = 0usize;
    for key in keys {
        let da = a.counts.get(key).copied().unwrap_or(0);
        let db = b.counts.get(key).copied().unwrap_or(0);
        dist += da.abs_diff(db) as usize;
    }
    dist
}

/// Heavy-atom L1 distance between two formulas (H and charge units ignored).
pub fn formula_heavy_l1(a: &Formula, b: &Formula) -> usize {
    let mut keys: Vec<&str> = a
        .counts
        .keys()
        .chain(b.counts.keys())
        .map(String::as_str)
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let mut dist = 0usize;
    for key in keys {
        if key == "H" || is_charge_key(key) {
            continue;
        }
        let da = a.counts.get(key).copied().unwrap_or(0);
        let db = b.counts.get(key).copied().unwrap_or(0);
        dist += da.abs_diff(db) as usize;
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn formula_l1_counts_hydrogen_heavy_skips_it() {
        let ethane = molecule_formula(&parse_mol("CC").unwrap());
        let ethene = molecule_formula(&parse_mol("C=C").unwrap());
        // Same heavy (C2); H differs 6 vs 4.
        assert_eq!(formula_heavy_l1(&ethane, &ethene), 0);
        assert_eq!(formula_l1(&ethane, &ethene), 2);
    }

    #[test]
    fn formula_charge_keys_are_strictly_positive_units() {
        let cation = molecule_formula(&parse_mol("C[N+](C)(C)C").unwrap());
        assert_eq!(cation.counts.get(CHARGE_PLUS), Some(&1));
        assert!(!cation.counts.contains_key(CHARGE_MINUS));
        assert_eq!(cation.charge, 1);

        let anion = molecule_formula(&parse_mol("CC(=O)[O-]").unwrap());
        assert_eq!(anion.counts.get(CHARGE_MINUS), Some(&1));
        assert!(!anion.counts.contains_key(CHARGE_PLUS));
        assert_eq!(anion.charge, -1);

        let zwitter = molecule_formula(&parse_mol("[NH3+]CC(=O)[O-]").unwrap());
        assert_eq!(zwitter.counts.get(CHARGE_PLUS), Some(&1));
        assert_eq!(zwitter.counts.get(CHARGE_MINUS), Some(&1));
        assert_eq!(zwitter.charge, 0);

        let neutral = molecule_formula(&parse_mol("CCO").unwrap());
        assert!(!neutral.counts.contains_key(CHARGE_PLUS));
        assert!(!neutral.counts.contains_key(CHARGE_MINUS));
        assert_eq!(neutral.charge, 0);
    }

    #[test]
    fn formula_l1_and_delta_include_charge_units() {
        let neutral = molecule_formula(&parse_mol("N").unwrap());
        let cation = molecule_formula(&parse_mol("[NH4+]").unwrap());
        // H also changes ammonia→ammonium; charge unit +1.
        assert!(formula_l1(&neutral, &cation) >= 1);
        assert_eq!(formula_heavy_l1(&neutral, &cation), 0);
        let delta = formula_delta(&neutral, &cation);
        assert_eq!(delta.counts.get(CHARGE_PLUS), Some(&1));
        assert_eq!(delta.charge, 1);
        assert!(!delta.counts.contains_key(CHARGE_MINUS));
    }
}
