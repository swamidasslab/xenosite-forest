//! Formula and structure-answer bags held by [`crate::forest_mol::ForestMol`].
//!
//! Not a Python `_forest` dict. [`ForestMol`] owns these as fields.

use std::collections::BTreeMap;
use std::rc::Rc;

use crate::mol::Molecule;

/// Element → count (hydrogens included) plus formal charge.
///
/// Cached on each [`crate::forest_mol::ForestMol`]. ``counts`` is the map
/// filters and pattern ``delta_formula`` compare against.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Formula {
    /// Map of element symbol → count, including hydrogens.
    pub counts: BTreeMap<String, i32>,
    pub charge: i32,
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
pub fn molecule_formula(mol: &Molecule) -> Formula {
    let mut counts = BTreeMap::new();
    let mut charge = 0i32;
    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        if z == 1 {
            *counts.entry("H".to_string()).or_insert(0) += 1;
            charge += i32::from(atom.charge);
            continue;
        }
        *counts.entry(atom.element.symbol().to_string()).or_insert(0) += 1;
        let hydrogens = mol.implicit_hydrogen_count(idx) as i32;
        if hydrogens > 0 {
            *counts.entry("H".to_string()).or_insert(0) += hydrogens;
        }
        charge += i32::from(atom.charge);
    }
    Formula { counts, charge }
}

/// Change in element counts (and charge) from ``before`` to ``after``.
///
/// Zero-count keys are omitted from ``counts``. Charge is always the signed
/// difference (may be zero).
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

/// Heavy-atom L1 distance between two formulas (H ignored; charge ignored).
///
/// Used by [`crate::find_path::MatchScoreSpec`] closeness / improvement scoring.
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
        if key == "H" {
            continue;
        }
        let da = a.counts.get(key).copied().unwrap_or(0);
        let db = b.counts.get(key).copied().unwrap_or(0);
        dist += da.abs_diff(db) as usize;
    }
    dist
}
