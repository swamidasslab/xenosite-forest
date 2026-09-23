//! Held molecule + mint-on-read `xf` facade.
//!
//! Python: `mol.xf` is a short-lived facade; answers live on `mol._forest`.
//! Rust: [`ForestMol`] owns the chematic mol and an optional [`Forest`].
//! [`ForestMol::xf`] mints a facade; nothing is stored under the name `xf`.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::forest::{Forest, Formula, Structure, empty_forest, forest_copy, molecule_formula};
use crate::mol::{ForestError, Molecule, canon_smiles, parse_mol, ranks};
use crate::smarts::smarts_matches;

/// A molecule that may carry a `_forest`. Parse does not install one.
#[derive(Clone)]
pub struct ForestMol {
    mol: Molecule,
    forest: RefCell<Option<Forest>>,
}

impl ForestMol {
    pub fn parse(smiles: &str) -> Result<Self, ForestError> {
        Ok(Self {
            mol: parse_mol(smiles)?,
            forest: RefCell::new(None),
        })
    }

    /// Chemistry copy with no forest (Python `Chem.Mol(parent)` / `rw_copy`).
    pub fn from_molecule(mol: Molecule) -> Self {
        Self {
            mol,
            forest: RefCell::new(None),
        }
    }

    pub fn mol(&self) -> &Molecule {
        &self.mol
    }

    /// True when `_forest` is already present. Does not install one.
    pub fn has_forest(&self) -> bool {
        self.forest.borrow().is_some()
    }

    pub fn ensure_forest(&self) -> Forest {
        let mut slot = self.forest.borrow_mut();
        if slot.is_none() {
            *slot = Some(empty_forest());
        }
        slot.as_ref().expect("forest installed").clone()
    }

    pub fn wipe_forest(&self) {
        *self.forest.borrow_mut() = None;
    }

    /// Mint a facade (Python `mol.xf`).
    pub fn xf(&self) -> Xf<'_> {
        Xf { held: self }
    }

    /// `copy_mol`: same structure, keep cache by identity.
    pub fn copy_mol(&self) -> Self {
        let forest = self
            .forest
            .borrow()
            .as_ref()
            .map(|forest| forest_copy(forest, true));
        Self {
            mol: self.mol.clone(),
            forest: RefCell::new(forest),
        }
    }

    /// `rw_copy`: editable chemistry copy, forest not carried.
    pub fn rw_copy(&self) -> Self {
        Self::from_molecule(self.mol.clone())
    }

    pub fn shares_cache(&self, other: &Self) -> bool {
        match (&*self.forest.borrow(), &*other.forest.borrow()) {
            (Some(left), Some(right)) => Rc::ptr_eq(&left.cache, &right.cache),
            _ => false,
        }
    }
}

/// Ephemeral facade for one [`ForestMol`]. Strong borrow only.
pub struct Xf<'a> {
    held: &'a ForestMol,
}

impl<'a> Xf<'a> {
    pub fn mol(&self) -> &'a Molecule {
        self.held.mol()
    }

    pub fn has_forest(&self) -> bool {
        self.held.has_forest()
    }

    /// Ensures forest; returns the parent (Python `forestmol`).
    pub fn forestmol(&self) -> &'a ForestMol {
        let _ = self.held.ensure_forest();
        self.held
    }

    pub fn forest(&self) -> Forest {
        self.held.ensure_forest()
    }

    pub fn is_terminal(&self) -> bool {
        self.held
            .forest
            .borrow()
            .as_ref()
            .is_some_and(|forest| forest.is_terminal_product)
    }

    pub fn clear_structure(&self) {
        let forest = self.held.ensure_forest();
        *forest.cache.borrow_mut() = Structure::default();
    }

    pub fn csmi(&self) -> Rc<str> {
        let forest = self.held.ensure_forest();
        let mut cache = forest.cache.borrow_mut();
        if cache.csmi.is_none() {
            cache.csmi = Some(Rc::from(canon_smiles(&self.held.mol)));
        }
        Rc::clone(cache.csmi.as_ref().expect("csmi filled"))
    }

    pub fn formula(&self) -> Rc<Formula> {
        let forest = self.held.ensure_forest();
        let mut cache = forest.cache.borrow_mut();
        if cache.formula.is_none() {
            cache.formula = Some(Rc::new(molecule_formula(&self.held.mol)));
        }
        Rc::clone(cache.formula.as_ref().expect("formula filled"))
    }

    pub fn topol_equiv(&self) -> Rc<Vec<usize>> {
        let forest = self.held.ensure_forest();
        let mut cache = forest.cache.borrow_mut();
        if cache.topol_equiv.is_none() {
            cache.topol_equiv = Some(Rc::new(ranks(&self.held.mol)));
        }
        Rc::clone(cache.topol_equiv.as_ref().expect("topol_equiv filled"))
    }

    pub fn smarts_matches(
        &self,
        smarts: &str,
    ) -> Result<Rc<Vec<BTreeMap<u16, usize>>>, ForestError> {
        let forest = self.held.ensure_forest();
        let mut cache = forest.cache.borrow_mut();
        if !cache.smarts_matches.contains_key(smarts) {
            let hits = smarts_matches(&self.held.mol, smarts)?;
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::canon_of;

    #[test]
    fn parse_does_not_install_forest() {
        let mol = ForestMol::parse("CC").unwrap();
        assert!(!mol.xf().has_forest());
        let _ = mol.xf().csmi();
        assert!(mol.xf().has_forest());
    }

    #[test]
    fn cached_answers_are_the_same_object() {
        let mol = ForestMol::parse("CCC").unwrap();
        assert!(Rc::ptr_eq(&mol.xf().topol_equiv(), &mol.xf().topol_equiv()));
        assert!(Rc::ptr_eq(&mol.xf().csmi(), &mol.xf().csmi()));
        let smarts = "[C:1][C:2]";
        assert!(Rc::ptr_eq(
            &mol.xf().smarts_matches(smarts).unwrap(),
            &mol.xf().smarts_matches(smarts).unwrap()
        ));
        assert!(Rc::ptr_eq(&mol.xf().formula(), &mol.xf().formula()));
    }

    #[test]
    fn writings_of_ethanol_share_cached_identity_not_a_spelling() {
        let mol = ForestMol::parse("CCO").unwrap();
        assert_eq!(mol.xf().csmi().as_ref(), canon_of("C(C)O").unwrap());
        let formula = mol.xf().formula();
        assert_eq!(formula.counts.get("C"), Some(&2));
        assert_eq!(formula.counts.get("O"), Some(&1));
        assert_eq!(formula.counts.get("H"), Some(&6));
        assert_eq!(formula.charge, 0);
    }

    #[test]
    fn fresh_mol_does_not_reuse_parent_caches() {
        let parent = ForestMol::parse("c1ccccc1O").unwrap();
        let parent_csmi = parent.xf().csmi();
        let product = ForestMol::from_molecule(parent.mol().clone());
        assert!(!product.has_forest());
        let _ = product.ensure_forest();
        assert!(!product.shares_cache(&parent));
        assert!(product.xf().forest().cache.borrow().csmi.is_none());
        assert_eq!(product.xf().csmi().as_ref(), parent_csmi.as_ref());
        assert!(!Rc::ptr_eq(&product.xf().csmi(), &parent_csmi));
    }

    #[test]
    fn copy_mol_shares_cache_rw_copy_does_not() {
        let parent = ForestMol::parse("CCO").unwrap();
        let _ = parent.xf().csmi();
        let copied = parent.copy_mol();
        assert!(copied.has_forest());
        assert!(copied.shares_cache(&parent));
        assert!(Rc::ptr_eq(&copied.xf().csmi(), &parent.xf().csmi()));
        let rw = parent.rw_copy();
        assert!(!rw.has_forest());
        let _ = rw.xf().csmi();
        assert!(!rw.shares_cache(&parent));
    }

    #[test]
    fn clear_structure_drops_answers_forest_stays() {
        let mol = ForestMol::parse("CC").unwrap();
        let first = mol.xf().csmi();
        mol.xf().clear_structure();
        assert!(mol.has_forest());
        assert!(mol.xf().forest().cache.borrow().csmi.is_none());
        let second = mol.xf().csmi();
        assert_eq!(first.as_ref(), second.as_ref());
        assert!(!Rc::ptr_eq(&first, &second));
    }

    #[test]
    fn wipe_forest_is_honest() {
        let mol = ForestMol::parse("CC").unwrap();
        assert!(!mol.xf().has_forest());
        let held = mol.xf().forestmol();
        assert!(std::ptr::eq(held, &mol) && mol.xf().has_forest());
        mol.wipe_forest();
        assert!(!mol.xf().has_forest());
    }
}
