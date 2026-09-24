//! [`ForestMol`]: owned chematic mol plus caches.
//!
//! Python attaches `_forest` and mints `xf` on a foreign RDKit `Mol`. Rust
//! owns the type, so those are not a layer. Caches live on this object:
//! structure answers (`csmi`, formula, …) and a shared per-system kekulé bag.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::atom_tracker::AtomTracker;
use crate::forest::{Formula, Structure, molecule_formula};
use crate::kekule::{KekuleCache, ensure_kekule_parents};
use crate::labels::{self, Tag};
use crate::mol::{ForestError, Molecule, atom_idx, canon_smiles, parse_mol, ranks};
use crate::smarts::smarts_matches;
use chematic::smiles::canonical_smiles_with_order;

/// Write sidecar labels onto chematic `Atom.tag` so apply/fragments copy them.
fn sync_tags_to_mol(mol: &mut Molecule, labels: &[Option<Tag>]) {
    let n = mol.atom_count().min(labels.len());
    for i in 0..n {
        mol.set_tag(atom_idx(i), labels[i].map(|t| t.0));
    }
}

/// Canonical SMILES → `parse_mol` (aromatize), tags remapped by visit order.
fn normalize_tagged_product(product: &Molecule) -> Result<Molecule, ForestError> {
    if product.atom_count() == 0 {
        return Ok(product.clone());
    }
    let (smi, order) = canonical_smiles_with_order(product);
    let mut fresh = parse_mol(&smi)?;
    if fresh.atom_count() != order.len() {
        // Fall back to plain canon_smiles round-trip without tags.
        return parse_mol(&canon_smiles(product));
    }
    for (new_i, &old_idx) in order.iter().enumerate() {
        let tag = product.atom(old_idx).tag;
        fresh.set_tag(atom_idx(new_i), tag);
    }
    Ok(fresh)
}

/// Owned molecule. Parse installs empty caches, filled on demand.
#[derive(Clone)]
pub struct ForestMol {
    mol: Molecule,
    labels: Vec<Option<Tag>>,
    tag_gen: Rc<Cell<u32>>,
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
        let mut mol = mol;
        let (labels, next) = labels::stamp(mol.atom_count());
        sync_tags_to_mol(&mut mol, &labels);
        Self {
            mol,
            labels,
            tag_gen: Rc::new(Cell::new(next)),
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::new(RefCell::new(KekuleCache::default())),
            is_terminal_product: Cell::new(false),
        }
    }

    /// Product of an edit: new structure bag, **same** kekulé cache `Rc`.
    ///
    /// Indexes of atoms that still exist are assumed stable (clone / append).
    /// SMIRKS apply that rewrites indexes must use [`Self::from_apply`].
    ///
    /// Unmodified systems still hit. An edit that changes a system's shape
    /// is a new [`crate::kekule::SystemKey`] and starts an empty bag.
    pub fn product(mol: Molecule, parent: &Self) -> Self {
        let mut mol = mol;
        let next = parent.tag_gen.get();
        let (labels, next) = labels::remap_index_stable(&parent.labels, mol.atom_count(), next);
        parent.tag_gen.set(next);
        sync_tags_to_mol(&mut mol, &labels);
        Self {
            mol,
            labels,
            tag_gen: Rc::clone(&parent.tag_gen),
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::clone(&parent.kekule),
            is_terminal_product: Cell::new(false),
        }
    }

    /// Product of a reindexing apply. `src_to_new[src] = Some(dst)` or `None`.
    pub fn from_apply(&self, mol: Molecule, src_to_new: &[Option<usize>]) -> Self {
        let mut mol = mol;
        let next = self.tag_gen.get();
        let (labels, next) = labels::remap_apply(&self.labels, src_to_new, mol.atom_count(), next);
        self.tag_gen.set(next);
        sync_tags_to_mol(&mut mol, &labels);
        Self {
            mol,
            labels,
            tag_gen: Rc::clone(&self.tag_gen),
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::clone(&self.kekule),
            is_terminal_product: Cell::new(false),
        }
    }

    /// Adopt a chematic product: re-parse via [`canon_smiles`] + [`parse_mol`]
    /// (aromaticity parity with the old string walk) while remapping `Atom.tag`
    /// by canonical visit order.
    pub fn adopt_product(&self, product: Molecule) -> Self {
        let normalized = normalize_tagged_product(&product).unwrap_or(product);
        let src_to_new = AtomTracker::src_to_new(self.mol(), &normalized);
        self.from_apply(normalized, &src_to_new)
    }

    /// Same tags after a permutation: `old_at_new[new] = old`.
    pub fn after_permute(&self, mol: &Molecule, old_at_new: &[usize]) -> Self {
        let mut mol = mol.clone();
        let labels = labels::remap_permute(&self.labels, old_at_new);
        sync_tags_to_mol(&mut mol, &labels);
        Self {
            mol,
            labels,
            tag_gen: Rc::clone(&self.tag_gen),
            structure: Rc::new(RefCell::new(Structure::default())),
            kekule: Rc::clone(&self.kekule),
            is_terminal_product: Cell::new(false),
        }
    }

    pub fn heavy_atom_count(&self) -> usize {
        self.mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() > 1)
            .count()
    }

    pub fn tag_of(&self, idx: usize) -> Option<Tag> {
        self.labels.get(idx).copied().flatten()
    }

    pub fn index_of(&self, tag: Tag) -> Option<usize> {
        self.labels.iter().position(|&held| held == Some(tag))
    }

    pub fn shares_tag_gen(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.tag_gen, &other.tag_gen)
    }

    pub fn mol(&self) -> &Molecule {
        &self.mol
    }

    /// Same structure: share structure + kekulé caches by identity.
    pub fn copy_mol(&self) -> Self {
        Self {
            mol: self.mol.clone(),
            labels: self.labels.clone(),
            tag_gen: Rc::clone(&self.tag_gen),
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
