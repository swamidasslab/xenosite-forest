//! `_forest` schema: three layers, like `records.Forest`.
//!
//! - `cache` — structure-dependent answers (`records.Structure`)
//! - `is_terminal_product` — mutable flag
//! - `immutable` / `atom_trace` wait on a full port
//!
//! `forest_copy(same_structure=true)` keeps `cache` by identity.
//! A different structure drops it.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::mol::Molecule;

/// Heavy-atom counts, hydrogens included, and formal charge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Formula {
    pub counts: BTreeMap<String, i32>,
    pub charge: i32,
}

/// Schema of `_forest["cache"]`. Keys fill in on first read.
#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub csmi: Option<Rc<str>>,
    pub formula: Option<Rc<Formula>>,
    pub topol_equiv: Option<Rc<Vec<usize>>>,
    pub smarts_matches: BTreeMap<String, Rc<Vec<BTreeMap<u16, usize>>>>,
}

/// The molecule's `_forest`. `cache` is the structure bag.
#[derive(Clone, Debug)]
pub struct Forest {
    pub cache: Rc<RefCell<Structure>>,
    pub is_terminal_product: bool,
}

/// Fresh forest shell: empty mutable cache.
pub fn empty_forest() -> Forest {
    Forest {
        cache: Rc::new(RefCell::new(Structure::default())),
        is_terminal_product: false,
    }
}

/// Copy a forest with the three-layer policy.
///
/// `same_structure`: keep `cache` by identity. Otherwise drop it.
pub fn forest_copy(forest: &Forest, same_structure: bool) -> Forest {
    Forest {
        cache: if same_structure {
            Rc::clone(&forest.cache)
        } else {
            Rc::new(RefCell::new(Structure::default()))
        },
        is_terminal_product: forest.is_terminal_product,
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_structure_copy_shares_cache() {
        let forest = empty_forest();
        forest.cache.borrow_mut().csmi = Some(Rc::from("CC"));
        let kept = forest_copy(&forest, true);
        assert!(Rc::ptr_eq(&forest.cache, &kept.cache));
        assert_eq!(kept.cache.borrow().csmi.as_deref(), Some("CC"));
    }

    #[test]
    fn different_structure_copy_drops_cache() {
        let forest = empty_forest();
        forest.cache.borrow_mut().csmi = Some(Rc::from("stale"));
        let dropped = forest_copy(&forest, false);
        assert!(!Rc::ptr_eq(&forest.cache, &dropped.cache));
        assert!(dropped.cache.borrow().csmi.is_none());
    }
}
