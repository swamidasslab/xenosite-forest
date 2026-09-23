//! [`ForestMol`]: owned chematic mol plus caches.
//!
//! Python attaches `_forest` and mints `xf` on a foreign RDKit `Mol`. Rust
//! owns the type, so those are not a layer. Caches live on this object:
//! structure answers (`csmi`, formula, …) and a shared per-system kekulé bag.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::forest::{Formula, Structure, molecule_formula};
use crate::kekule::{KekuleCache, ensure_kekule_parents};
use crate::mol::{ForestError, Molecule, canon_smiles, parse_mol, ranks};
use crate::smarts::smarts_matches;

/// Owned molecule. Parse installs empty caches, filled on demand.
#[derive(Clone)]
pub struct ForestMol {
    mol: Molecule,
    structure: Rc<RefCell<Structure>>,
    kekule: Rc<RefCell<KekuleCache>>,
    pub is_terminal_product: Cell<bool>,
}

impl ForestMol {
    pub fn parse(smiles: &str) -> Result<Self, ForestError> {
        Ok(Self::new(parse_mol(smiles)?))
    }

    /// Wrap chemistry with **new** caches (disconnected from any parent tree).
    pub fn new(mol: Molecule) -> Self {
        Self {
            mol,
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::new(RefCell::new(KekuleCache::default())),
            is_terminal_product: Cell::new(false),
        }
    }

    /// Product of an edit: new structure bag, **same** kekulé cache `Rc`.
    ///
    /// Unmodified systems still hit. An edit that changes a system's shape
    /// is a new [`crate::kekule::SystemKey`] and starts an empty bag.
    pub fn product(mol: Molecule, parent: &Self) -> Self {
        Self {
            mol,
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::clone(&parent.kekule),
            is_terminal_product: Cell::new(false),
        }
    }

    pub fn mol(&self) -> &Molecule {
        &self.mol
    }

    /// Same structure: share structure + kekulé caches by identity.
    pub fn copy_mol(&self) -> Self {
        Self {
            mol: self.mol.clone(),
            structure: Rc::clone(&self.structure),
            kekule: Rc::clone(&self.kekule),
            is_terminal_product: Cell::new(self.is_terminal_product.get()),
        }
    }

    /// Chemistry clone ready to edit: new structure bag, shared kekulé `Rc`.
    pub fn edit_copy(&self) -> Self {
        Self::product(self.mol.clone(), self)
    }

    pub fn shares_structure(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.structure, &other.structure)
    }

    pub fn shares_kekule(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.kekule, &other.kekule)
    }

    pub fn kekule(&self) -> Rc<RefCell<KekuleCache>> {
        Rc::clone(&self.kekule)
    }

    pub fn clear_structure(&self) {
        *self.structure.borrow_mut() = Structure::default();
    }

    pub fn csmi(&self) -> Rc<str> {
        let mut cache = self.structure.borrow_mut();
        if cache.csmi.is_none() {
            cache.csmi = Some(Rc::from(canon_smiles(&self.mol)));
        }
        Rc::clone(cache.csmi.as_ref().expect("csmi filled"))
    }

    pub fn formula(&self) -> Rc<Formula> {
        let mut cache = self.structure.borrow_mut();
        if cache.formula.is_none() {
            cache.formula = Some(Rc::new(molecule_formula(&self.mol)));
        }
        Rc::clone(cache.formula.as_ref().expect("formula filled"))
    }

    pub fn topol_equiv(&self) -> Rc<Vec<usize>> {
        let mut cache = self.structure.borrow_mut();
        if cache.topol_equiv.is_none() {
            cache.topol_equiv = Some(Rc::new(ranks(&self.mol)));
        }
        Rc::clone(cache.topol_equiv.as_ref().expect("topol_equiv filled"))
    }

    pub fn smarts_matches(
        &self,
        smarts: &str,
    ) -> Result<Rc<Vec<BTreeMap<u16, usize>>>, ForestError> {
        let mut cache = self.structure.borrow_mut();
        if !cache.smarts_matches.contains_key(smarts) {
            let hits = smarts_matches(&self.mol, smarts)?;
            cache
                .smarts_matches
                .insert(smarts.to_string(), Rc::new(hits));
        }
        Ok(Rc::clone(
            cache
                .smarts_matches
                .get(smarts)
                .expect("smarts_matches filled"),
        ))
    }

    /// Kekulize the system containing this bond. Shared with relatives.
    pub fn ensure_kekule(&self, left: usize, right: usize) {
        let mut cache = self.kekule.borrow_mut();
        ensure_kekule_parents(&self.mol, left, right, &mut cache);
    }

    pub fn reactant_parent(
        &self,
        mapped: &BTreeMap<u16, usize>,
        smirks: &str,
    ) -> Result<Molecule, ForestError> {
        let mut cache = self.kekule.borrow_mut();
        crate::kekule::reactant_parent(&self.mol, mapped, smirks, &mut cache)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::canon_of;

    #[test]
    fn parse_has_caches_empty_until_read() {
        let mol = ForestMol::parse("CC").unwrap();
        assert!(mol.structure.borrow().csmi.is_none());
        let _ = mol.csmi();
        assert!(mol.structure.borrow().csmi.is_some());
    }

    #[test]
    fn cached_answers_are_the_same_object() {
        let mol = ForestMol::parse("CCC").unwrap();
        assert!(Rc::ptr_eq(&mol.topol_equiv(), &mol.topol_equiv()));
        assert!(Rc::ptr_eq(&mol.csmi(), &mol.csmi()));
        let smarts = "[C:1][C:2]";
        assert!(Rc::ptr_eq(
            &mol.smarts_matches(smarts).unwrap(),
            &mol.smarts_matches(smarts).unwrap()
        ));
        assert!(Rc::ptr_eq(&mol.formula(), &mol.formula()));
    }

    #[test]
    fn writings_of_ethanol_share_cached_identity_not_a_spelling() {
        let mol = ForestMol::parse("CCO").unwrap();
        assert_eq!(mol.csmi().as_ref(), canon_of("C(C)O").unwrap());
        let formula = mol.formula();
        assert_eq!(formula.counts.get("C"), Some(&2));
        assert_eq!(formula.counts.get("O"), Some(&1));
        assert_eq!(formula.counts.get("H"), Some(&6));
        assert_eq!(formula.charge, 0);
    }

    #[test]
    fn disconnected_wrap_does_not_share_caches() {
        let parent = ForestMol::parse("c1ccccc1O").unwrap();
        let parent_csmi = parent.csmi();
        let product = ForestMol::new(parent.mol().clone());
        assert!(!product.shares_structure(&parent));
        assert!(!product.shares_kekule(&parent));
        assert!(product.structure.borrow().csmi.is_none());
        assert_eq!(product.csmi().as_ref(), parent_csmi.as_ref());
        assert!(!Rc::ptr_eq(&product.csmi(), &parent_csmi));
    }

    #[test]
    fn copy_shares_both_caches_edit_copy_shares_only_kekule() {
        let parent = ForestMol::parse("CCO").unwrap();
        let _ = parent.csmi();
        let copied = parent.copy_mol();
        assert!(copied.shares_structure(&parent));
        assert!(copied.shares_kekule(&parent));
        assert!(Rc::ptr_eq(&copied.csmi(), &parent.csmi()));
        let edited = parent.edit_copy();
        assert!(!edited.shares_structure(&parent));
        assert!(edited.shares_kekule(&parent));
        let _ = edited.csmi();
        assert!(!Rc::ptr_eq(&edited.csmi(), &parent.csmi()));
    }

    #[test]
    fn clear_structure_drops_csmi_keeps_kekule_rc() {
        let mol = ForestMol::parse("CC").unwrap();
        let first = mol.csmi();
        let kekule = mol.kekule();
        mol.clear_structure();
        assert!(mol.structure.borrow().csmi.is_none());
        assert!(Rc::ptr_eq(&kekule, &mol.kekule()));
        let second = mol.csmi();
        assert_eq!(first.as_ref(), second.as_ref());
        assert!(!Rc::ptr_eq(&first, &second));
    }
}
