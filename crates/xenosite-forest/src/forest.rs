//! Formula and structure-answer bags held by [`crate::forest_mol::ForestMol`].
//!
//! Not a Python `_forest` dict. [`ForestMol`] owns these as fields.

use std::collections::BTreeMap;
use std::rc::Rc;

use crate::mol::Molecule;

/// Heavy-atom counts, hydrogens included, and formal charge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Formula {
    pub counts: BTreeMap<String, i32>,
    pub charge: i32,
}

/// Structure-dependent answers, filled on first read. Shared by `Rc` when
/// the chemistry is the same; a new bag after an edit.
#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub csmi: Option<Rc<str>>,
    pub formula: Option<Rc<Formula>>,
    pub topol_equiv: Option<Rc<Vec<usize>>>,
    pub smarts_matches: BTreeMap<String, Rc<Vec<BTreeMap<u16, usize>>>>,
}

pub fn molecule_formula(mol: &Molecule) -> Formula {
    let mut counts = BTreeMap::new();
    let mut charge = 0i32;
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 1 {
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
