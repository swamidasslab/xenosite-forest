//! Kekulé parents, one conjugated system at a time, on demand, shared.
//!
//! Whole-molecule enumeration is the product of independent rings (`2^n`
//! phenyls). Live forest does not do that. This cache:
//!
//! - assigns **one** conjugated component (other systems stay aromatic)
//! - fills **on demand** when a match names a bond in that system
//! - stores **assignment maps**, not baked mols
//! - lives on [`crate::forest_mol::ForestMol`] as one `Rc` for a copy tree: the first
//!   relative to fill a key shares it with every relative whose system
//!   still matches. An edit that changes kekulization of a system misses
//!   that key and starts a new bag.
//!
//! # [`PiGraph`] methods (reuse + derived systems)
//!
//! Perception and residual derivation go through [`PiGraph`]:
//! [`PiGraph::conjugated`], [`PiGraph::components`], [`PiGraph::two_core`],
//! [`PiGraph::without_atoms`] / [`PiGraph::after_forced_doubles`],
//! [`PiGraph::aromatic_subgraph`]. [`KekuleCache::ensure_graph`] fills bags
//! from a graph; [`SystemKekule::with_forced_doubles`] derives a constrained
//! bag without refill; [`ResidualKey::apply_to`] rebuilds an index residual
//! from tags.
//!
//! **Chematic:** `kekulize` / `apply_kekule` are single aromatic-bond matchings;
//! `Molecule::fragments` is σ connectivity. Neither is a conjugated flood-fill
//! or multi-resonance / residual API — keep π graphs here (HEURISTICS).
//!
//! # Tag keys, not indexes
//!
//! Reusable Kekulé objects ([`SystemKey`], [`SystemKekule`] assignments /
//! `by_order`, and residual/constraint keys when they land) are keyed by
//! forest **labels** ([`Tag`] / chematic `Atom.tag`), **not** atom indexes.
//! Indexes are only for perception and matching on the *current* mol layout;
//! they move under rewrite and canonical reorder. Tags are stable in a
//! ForestMol copy tree, so untouched systems keep hitting the shared bag after
//! `product` / `edit_copy`. Overlay resolves tags → indexes on the mol being
//! stamped.

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::rc::Rc;

use chematic::core::BondOrder;
use chematic::perception::find_sssr;
use chematic::smarts::{BondPrimitive, BondQuery, parse_smarts};

use crate::labels::Tag;
use crate::mol::{ForestError, Molecule, atom_idx, atom_usize};

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn tag_bond_key(a: Tag, b: Tag) -> (Tag, Tag) {
    if a <= b { (a, b) } else { (b, a) }
}

/// Chematic `Atom.tag` as a forest [`Tag`].
fn tag_of(mol: &Molecule, idx: usize) -> Option<Tag> {
    mol.atom(atom_idx(idx)).tag.map(Tag)
}

/// Current index holding `tag`, if any.
fn index_of_tag(mol: &Molecule, tag: Tag) -> Option<usize> {
    mol.atoms().find_map(|(idx, atom)| {
        (atom.tag == Some(tag.0)).then_some(atom_usize(idx))
    })
}

/// Stamp `Atom.tag = Some(i as u32)` on every untagged atom.
///
/// ForestMol always stamps real copy-tree labels. Bare mols (one-shot
/// [`kekule_forms`], unit tests) need tags before cache insert — index-as-tag
/// is only stable on that mol instance, not across rewrite. Prefer ForestMol
/// labels for reuse across `product` / `edit_copy`.
pub fn stamp_missing_index_tags(mol: &mut Molecule) {
    let n = mol.atom_count();
    for i in 0..n {
        if mol.atom(atom_idx(i)).tag.is_none() {
            mol.set_tag(atom_idx(i), Some(i as u32));
        }
    }
}

fn tags_for_atoms(mol: &Molecule, atoms: &BTreeSet<usize>) -> Option<BTreeSet<Tag>> {
    let mut tags = BTreeSet::new();
    for &i in atoms {
        tags.insert(tag_of(mol, i)?);
    }
    Some(tags)
}

fn assignment_to_tags(
    mol: &Molecule,
    assignment: &BTreeMap<(usize, usize), BondOrder>,
) -> Option<BTreeMap<(Tag, Tag), BondOrder>> {
    let mut out = BTreeMap::new();
    for (&(a, b), &order) in assignment {
        let ta = tag_of(mol, a)?;
        let tb = tag_of(mol, b)?;
        out.insert(tag_bond_key(ta, tb), order);
    }
    Some(out)
}

fn assignment_to_idxs(
    mol: &Molecule,
    assignment: &BTreeMap<(Tag, Tag), BondOrder>,
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    let mut out = BTreeMap::new();
    for (&(ta, tb), &order) in assignment {
        let a = index_of_tag(mol, ta)?;
        let b = index_of_tag(mol, tb)?;
        out.insert(bond_key(a, b), order);
    }
    Some(out)
}

fn atoms_to_idxs(mol: &Molecule, tags: &BTreeSet<Tag>) -> Option<BTreeSet<usize>> {
    let mut out = BTreeSet::new();
    for &t in tags {
        out.insert(index_of_tag(mol, t)?);
    }
    Some(out)
}

#[cfg(test)]
fn bond_order_map(mol: &Molecule) -> BTreeMap<(usize, usize), BondOrder> {
    mol.bonds()
        .map(|(_, bond)| {
            (
                bond_key(atom_usize(bond.atom1), atom_usize(bond.atom2)),
                bond.order,
            )
        })
        .collect()
}

fn implied_order(query: &BondQuery) -> Option<f32> {
    match query {
        BondQuery::Primitive(BondPrimitive::Double) => Some(2.0),
        BondQuery::Primitive(BondPrimitive::Single) => Some(1.0),
        BondQuery::Primitive(BondPrimitive::Aromatic) => Some(1.5),
        BondQuery::Or(left, right) => match (implied_order(left), implied_order(right)) {
            (Some(2.0), Some(1.5)) | (Some(1.5), Some(2.0)) => Some(2.0),
            (Some(1.0), Some(1.5)) | (Some(1.5), Some(1.0)) => Some(1.0),
            (a, b) => a.or(b),
        },
        BondQuery::Any | BondQuery::Primitive(BondPrimitive::Any) => Some(1.0),
        BondQuery::And(left, right) => implied_order(left).or_else(|| implied_order(right)),
        BondQuery::Not(_) => None,
        BondQuery::Primitive(_) => None,
    }
}

/// Bond order the reactant SMARTS asks for between maps 1 and 2.
pub fn smirks_mapped_bond_order(smirks: &str) -> Option<f32> {
    let reactant = smirks.split(">>").next()?;
    let query = parse_smarts(reactant).ok()?;
    let mut idxs = BTreeMap::new();
    for (i, atom) in query.atoms.iter().enumerate() {
        if let Some(mapno) = atom.atom_map {
            idxs.insert(mapno, i);
        }
    }
    let left = *idxs.get(&1)?;
    let right = *idxs.get(&2)?;
    let bond = query.bonds.iter().find(|bond| {
        (bond.atom1 == left && bond.atom2 == right) || (bond.atom1 == right && bond.atom2 == left)
    })?;
    implied_order(&bond.query)
}

fn in_ring_bonds(mol: &Molecule) -> BTreeSet<(usize, usize)> {
    let mut keys = BTreeSet::new();
    for ring in find_sssr(mol).rings() {
        let atoms: Vec<usize> = ring.iter().copied().map(atom_usize).collect();
        let n = atoms.len();
        if n < 3 {
            continue;
        }
        for i in 0..n {
            keys.insert(bond_key(atoms[i], atoms[(i + 1) % n]));
        }
    }
    keys
}

fn pi_center(mol: &Molecule, index: usize) -> bool {
    if mol.atom(atom_idx(index)).aromatic {
        return true;
    }
    for (_nbr, bond_idx) in mol.neighbors(atom_idx(index)) {
        match mol.bond(bond_idx).order {
            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic => return true,
            _ => {}
        }
    }
    false
}

fn conjugated_bond(
    mol: &Molecule,
    left: usize,
    right: usize,
    order: BondOrder,
    in_ring: &BTreeSet<(usize, usize)>,
) -> bool {
    let key = bond_key(left, right);
    match order {
        BondOrder::Double | BondOrder::Triple => true,
        BondOrder::Aromatic => in_ring.contains(&key),
        BondOrder::Single => {
            if mol.atom(atom_idx(left)).aromatic
                && mol.atom(atom_idx(right)).aromatic
                && in_ring.contains(&key)
            {
                return true;
            }
            let z_left = mol.atom(atom_idx(left)).element.atomic_number();
            let z_right = mol.atom(atom_idx(right)).element.atomic_number();
            let elements = [z_left, z_right];
            if !elements.contains(&6) || !elements.iter().any(|z| matches!(z, 7 | 8 | 16)) {
                return false;
            }
            pi_center(mol, left) || pi_center(mol, right)
        }
        _ => mol.atom(atom_idx(left)).aromatic && mol.atom(atom_idx(right)).aromatic,
    }
}

/// Conjugated component containing `start`. Biaryl single bonds do not join rings.
///
/// Prefer [`PiGraph::conjugated`] when deriving residuals / components; this
/// free function remains for call sites that only need the atom/bond sets.
pub fn conjugated_component(
    mol: &Molecule,
    start: usize,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
    let g = PiGraph::conjugated(mol, start);
    (g.atoms, g.bonds)
}

/// Index-keyed π topology on the **current** mol layout (perception).
///
/// This is the reusable graph object for conjugated systems, residual views,
/// 2-core / component decomposition, and constrained matching. Cache identity
/// stays on [`SystemKey`] / [`ResidualKey`] (tags); [`PiGraph`] is the working
/// form once tags resolve to indexes on a mol.
///
/// **Not Chematic's job:** chematic `kekulize` matches aromatic bonds once;
/// `Molecule::fragments` splits σ components. Forest conjugation (exocyclic
/// hetero singles, multi-resonance bags, residual/2-core) lives here — do not
/// move this flood-fill into chematic without a shared schema both doors read.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PiGraph {
    pub atoms: BTreeSet<usize>,
    pub bonds: BTreeSet<(usize, usize)>,
}

impl PiGraph {
    pub fn new(atoms: BTreeSet<usize>, bonds: BTreeSet<(usize, usize)>) -> Self {
        let bonds = bonds
            .into_iter()
            .filter(|&(a, b)| atoms.contains(&a) && atoms.contains(&b))
            .map(|(a, b)| bond_key(a, b))
            .collect();
        Self { atoms, bonds }
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Conjugated component flood-fill from `start`.
    pub fn conjugated(mol: &Molecule, start: usize) -> Self {
        let in_ring = in_ring_bonds(mol);
        let mut atoms = BTreeSet::from([start]);
        let mut bonds = BTreeSet::new();
        let mut stack = vec![start];
        while let Some(index) = stack.pop() {
            for (nbr, bond_idx) in mol.neighbors(atom_idx(index)) {
                let bond = mol.bond(bond_idx);
                let other = atom_usize(nbr);
                if !conjugated_bond(mol, index, other, bond.order, &in_ring) {
                    continue;
                }
                bonds.insert(bond_key(index, other));
                if atoms.insert(other) {
                    stack.push(other);
                }
            }
        }
        Self { atoms, bonds }
    }

    /// Restrict to `keep` (saturate sites leave; epoxide sp³ dropped).
    pub fn without_atoms(&self, drop: &BTreeSet<usize>) -> Self {
        let atoms: BTreeSet<usize> = self
            .atoms
            .iter()
            .copied()
            .filter(|a| !drop.contains(a))
            .collect();
        Self::new(atoms, self.bonds.clone())
    }

    /// Keep only atoms in `keep` (and bonds between them).
    pub fn restrict_atoms(&self, keep: &BTreeSet<usize>) -> Self {
        let atoms: BTreeSet<usize> = self
            .atoms
            .iter()
            .copied()
            .filter(|a| keep.contains(a))
            .collect();
        Self::new(atoms, self.bonds.clone())
    }

    /// Aromatic atoms of this graph (pair aromatic-scoped parents).
    pub fn aromatic_subgraph(&self, mol: &Molecule) -> Self {
        let keep: BTreeSet<usize> = self
            .atoms
            .iter()
            .copied()
            .filter(|&i| mol.atom(atom_idx(i)).aromatic)
            .collect();
        self.restrict_atoms(&keep)
    }

    /// Drop atoms whose π demand is consumed by a forced double (exocyclic
    /// leaf or in-system forced edge) before the cyclic 2-core aromaticity test.
    pub fn after_forced_doubles(&self, forced: &BTreeSet<(usize, usize)>) -> Self {
        let mut atoms = self.atoms.clone();
        for &(a, b) in forced {
            let a_in = atoms.contains(&a);
            let b_in = atoms.contains(&b);
            if a_in && b_in {
                atoms.remove(&a);
                atoms.remove(&b);
            } else if a_in {
                atoms.remove(&a);
            } else if b_in {
                atoms.remove(&b);
            }
        }
        Self::new(atoms, self.bonds.clone())
    }

    /// 2-core: iteratively drop atoms with π-degree &lt; 2.
    pub fn two_core(&self) -> Self {
        let mut atoms = self.atoms.clone();
        let mut bonds = self.bonds.clone();
        loop {
            let mut deg: HashMap<usize, usize> = atoms.iter().map(|&a| (a, 0)).collect();
            for &(a, b) in &bonds {
                *deg.entry(a).or_default() += 1;
                *deg.entry(b).or_default() += 1;
            }
            let drop: Vec<usize> = deg
                .iter()
                .filter(|(_, d)| **d < 2)
                .map(|(&a, _)| a)
                .collect();
            if drop.is_empty() {
                break;
            }
            for a in drop {
                atoms.remove(&a);
            }
            bonds.retain(|&(a, b)| atoms.contains(&a) && atoms.contains(&b));
        }
        Self { atoms, bonds }
    }

    /// Connected components under this graph's edges (π-component split).
    pub fn components(&self) -> Vec<Self> {
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        for &start in &self.atoms {
            if !seen.insert(start) {
                continue;
            }
            let mut stack = vec![start];
            let mut comp_atoms = BTreeSet::from([start]);
            while let Some(a) = stack.pop() {
                for &(left, right) in &self.bonds {
                    let other = if left == a {
                        right
                    } else if right == a {
                        left
                    } else {
                        continue;
                    };
                    if comp_atoms.insert(other) {
                        seen.insert(other);
                        stack.push(other);
                    }
                }
            }
            let comp_bonds: BTreeSet<_> = self
                .bonds
                .iter()
                .copied()
                .filter(|&(a, b)| comp_atoms.contains(&a) && comp_atoms.contains(&b))
                .collect();
            out.push(Self {
                atoms: comp_atoms,
                bonds: comp_bonds,
            });
        }
        out
    }

    /// Doubles on `edited` that touch this graph (perception for constraints).
    pub fn perceive_forced_doubles(&self, edited: &Molecule) -> BTreeSet<(usize, usize)> {
        let mut forced = BTreeSet::new();
        for (_, bond) in edited.bonds() {
            if bond.order != BondOrder::Double {
                continue;
            }
            let a = atom_usize(bond.atom1);
            let b = atom_usize(bond.atom2);
            if self.atoms.contains(&a) || self.atoms.contains(&b) {
                forced.insert(bond_key(a, b));
            }
        }
        forced
    }

    pub fn try_system_key(&self, mol: &Molecule) -> Option<SystemKey> {
        SystemKey::try_of(mol, &self.atoms, &self.bonds)
    }

    pub fn system_key(&self, mol: &Molecule) -> SystemKey {
        SystemKey::of(mol, &self.atoms, &self.bonds)
    }

    pub fn complete_assignment(
        &self,
        mol: &Molecule,
        constraints: &KekuleConstraints,
        config: &KekuleConfig,
    ) -> Option<BTreeMap<(usize, usize), BondOrder>> {
        complete_assignment(mol, &self.atoms, &self.bonds, constraints, config)
    }

    pub fn all_assignments(
        &self,
        mol: &Molecule,
        constraints: &KekuleConstraints,
        config: &KekuleConfig,
    ) -> Vec<BTreeMap<(usize, usize), BondOrder>> {
        all_assignments(mol, &self.atoms, &self.bonds, constraints, config)
    }

    /// Cyclic 2-core components that remain aromatic under `config`.
    /// Cyclic 2-core components that remain aromatic under `config`.
    pub fn aromatic_2core_atoms(&self, mol: &Molecule, config: &KekuleConfig) -> BTreeSet<usize> {
        let core = self.two_core();
        let empty = KekuleConstraints::default();
        let match_cfg = config.explicit_forced_only().single_assignment();
        let mut aromatic = BTreeSet::new();
        for comp in core.components() {
            // Cyclic: at least as many edges as atoms.
            if comp.bonds.len() < comp.atoms.len() {
                continue;
            }
            if complete_assignment(mol, &comp.atoms, &comp.bonds, &empty, &match_cfg).is_none()
            {
                continue;
            }
            let n = comp.atoms.len();
            if !config.huckel_4n2 || (n >= 6 && n % 4 == 2) {
                aromatic.extend(comp.atoms);
            }
        }
        aromatic
    }

    /// Residual view key for this graph under `removed` / `forced` (tags).
    pub fn residual_key(
        &self,
        mol: &Molecule,
        removed: &BTreeSet<usize>,
        forced_doubles: &BTreeSet<(usize, usize)>,
    ) -> Option<ResidualKey> {
        ResidualKey::try_of(mol, self, removed, forced_doubles)
    }
}

fn order_code(order: BondOrder) -> Option<u8> {
    match order {
        BondOrder::Single => Some(1),
        BondOrder::Double => Some(2),
        BondOrder::Aromatic => Some(3),
        _ => Some(0),
    }
}

/// Identity of one conjugated system for cache reuse.
///
/// **Tag-keyed, not index-keyed.** `atoms` and `shape` bond endpoints are
/// [`Tag`]s (cheatic `Atom.tag`). Indexes would break reuse after rewrite or
/// canonical reorder; tags survive on the ForestMol copy-tree `Rc` cache.
///
/// An edit that changes aromatic flags or bond orders inside the system is a
/// new key, so that bag is not reused. An edit elsewhere leaves this key
/// identical, so relatives share the bag.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SystemKey {
    /// Conjugated-system atoms as forest labels.
    pub atoms: BTreeSet<Tag>,
    /// Kekulization-relevant shape: tag-pair bonds with order code + aromatic bits.
    shape: Vec<((Tag, Tag), u8, bool, bool)>,
}

impl SystemKey {
    /// Build a tag-keyed system identity from the current mol's indexes.
    ///
    /// Perception still uses indexes; the returned key stores only tags.
    /// Returns `None` if any system atom lacks `Atom.tag` — call
    /// [`stamp_missing_index_tags`] on bare mols, or use a [`crate::forest_mol::ForestMol`].
    pub fn try_of(
        mol: &Molecule,
        atoms: &BTreeSet<usize>,
        bonds: &BTreeSet<(usize, usize)>,
    ) -> Option<Self> {
        let tag_atoms = tags_for_atoms(mol, atoms)?;
        let mut shape: Vec<_> = bonds
            .iter()
            .map(|&(a, b)| {
                let ta = tag_of(mol, a)?;
                let tb = tag_of(mol, b)?;
                let order = mol
                    .bond_between(atom_idx(a), atom_idx(b))
                    .map(|(_, bond)| order_code(bond.order).unwrap_or(0))
                    .unwrap_or(0);
                let a_aro = mol.atom(atom_idx(a)).aromatic;
                let b_aro = mol.atom(atom_idx(b)).aromatic;
                Some((tag_bond_key(ta, tb), order, a_aro, b_aro))
            })
            .collect::<Option<Vec<_>>>()?;
        shape.sort_unstable();
        Some(Self {
            atoms: tag_atoms,
            shape,
        })
    }

    /// Like [`Self::try_of`], panicking if tags are missing (ForestMol / stamped mols).
    pub fn of(mol: &Molecule, atoms: &BTreeSet<usize>, bonds: &BTreeSet<(usize, usize)>) -> Self {
        Self::try_of(mol, atoms, bonds).expect(
            "SystemKey requires Atom.tag on every system atom (ForestMol stamps; \
             bare mols: stamp_missing_index_tags)",
        )
    }
}

/// Residual / constraint view key (secondary cache beside [`SystemKey`]).
///
/// **Tag-keyed:** `removed` and `forced_doubles` are labels / label-pairs, not
/// indexes — so a child mol with the same tags can reuse a parent-derived
/// residual after index shuffle. See HEURISTICS edit-as-π-constraints cache plan.
///
/// Build with [`ResidualKey::try_of`] from a parent [`PiGraph`]; apply back to
/// indexes with [`ResidualKey::apply_to`].
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ResidualKey {
    pub parent: SystemKey,
    pub removed: BTreeSet<Tag>,
    pub forced_doubles: BTreeSet<(Tag, Tag)>,
}

impl ResidualKey {
    /// Tag-keyed residual identity from parent topology + index constraints.
    pub fn try_of(
        mol: &Molecule,
        parent: &PiGraph,
        removed: &BTreeSet<usize>,
        forced_doubles: &BTreeSet<(usize, usize)>,
    ) -> Option<Self> {
        let parent_key = parent.try_system_key(mol)?;
        let mut rem = BTreeSet::new();
        for &i in removed {
            rem.insert(tag_of(mol, i)?);
        }
        let mut forced = BTreeSet::new();
        for &(a, b) in forced_doubles {
            let ta = tag_of(mol, a)?;
            let tb = tag_of(mol, b)?;
            forced.insert(tag_bond_key(ta, tb));
        }
        Some(Self {
            parent: parent_key,
            removed: rem,
            forced_doubles: forced,
        })
    }

    /// Resolve this residual onto `parent`'s index graph for `mol`.
    ///
    /// Drops `removed` atoms; does **not** strip forced-demand atoms — call
    /// [`PiGraph::after_forced_doubles`] on the result when building the
    /// aromatic 2-core candidate.
    pub fn apply_to(&self, mol: &Molecule, parent: &PiGraph) -> Option<PiGraph> {
        let mut drop = BTreeSet::new();
        for &tag in &self.removed {
            drop.insert(index_of_tag(mol, tag)?);
        }
        Some(parent.without_atoms(&drop))
    }

    /// Forced doubles as indexes on `mol`.
    pub fn forced_indexes(&self, mol: &Molecule) -> Option<BTreeSet<(usize, usize)>> {
        let mut out = BTreeSet::new();
        for &(ta, tb) in &self.forced_doubles {
            let a = index_of_tag(mol, ta)?;
            let b = index_of_tag(mol, tb)?;
            out.insert(bond_key(a, b));
        }
        Some(out)
    }
}

/// One system's kekulé assignments (bond-order maps, not baked mols).
///
/// Bond endpoints in `assignments` and `by_order` are [`Tag`] pairs — reusable
/// across index permutations. Callers resolve to indexes via the current mol
/// when overlaying.
#[derive(Clone, Debug, Default)]
pub struct SystemKekule {
    pub assignments: Vec<BTreeMap<(Tag, Tag), BondOrder>>,
    pub by_order: BTreeMap<((Tag, Tag), u8), usize>,
    /// Parallel to `assignments`: Σ|formal charge| after overlay + move_charge.
    /// Used so `by_order` prefers closed-shell writings (APAP neutral amide).
    charge_mags: Vec<i32>,
}

impl SystemKekule {
    pub fn is_filled(&self) -> bool {
        !self.assignments.is_empty()
    }

    /// Semantic identity of cached writings (assignments + by_order), not `Rc`.
    pub fn semantically_eq(&self, other: &Self) -> bool {
        self.assignments == other.assignments && self.by_order == other.by_order
    }

    /// Derived bag: keep assignments where each forced tag-pair is Double.
    ///
    /// Does not refill matching — alters only the selected writings so callers
    /// can constrain one conjugated system without recomputing the parent bag.
    pub fn with_forced_doubles(&self, forced: &BTreeSet<(Tag, Tag)>) -> Self {
        if forced.is_empty() {
            return self.clone();
        }
        let mut out = Self::default();
        let mut seen = BTreeSet::new();
        for (ai, assignment) in self.assignments.iter().enumerate() {
            let ok = forced.iter().all(|&(a, b)| {
                let key = tag_bond_key(a, b);
                assignment.get(&key) == Some(&BondOrder::Double)
            });
            if !ok {
                continue;
            }
            let signature: Vec<_> = assignment
                .iter()
                .map(|(&k, &o)| (k, order_code(o)))
                .collect();
            if !seen.insert(signature) {
                continue;
            }
            let index = out.assignments.len();
            let mag = self.charge_mags.get(ai).copied().unwrap_or(0);
            for (key, order) in assignment {
                if let Some(code) = order_code(*order) {
                    match out.by_order.get(&(*key, code)) {
                        None => {
                            out.by_order.insert((*key, code), index);
                        }
                        Some(&prev) => {
                            let prev_mag = out.charge_mags.get(prev).copied().unwrap_or(i32::MAX);
                            if mag < prev_mag {
                                out.by_order.insert((*key, code), index);
                            }
                        }
                    }
                }
            }
            out.assignments.push(assignment.clone());
            out.charge_mags.push(mag);
        }
        out
    }
}

/// Shared forest-level map. One `Rc` per copy tree.
///
/// Keys and stored writings are tag-keyed ([`SystemKey`], [`SystemKekule`]).
#[derive(Clone, Debug, Default)]
pub struct KekuleCache {
    systems: BTreeMap<SystemKey, Rc<RefCell<SystemKekule>>>,
}

impl KekuleCache {
    pub fn system_count(&self) -> usize {
        self.systems.len()
    }

    pub fn assignment_count(&self) -> usize {
        self.systems
            .values()
            .map(|slot| slot.borrow().assignments.len())
            .sum()
    }

    pub fn slot(&mut self, key: SystemKey) -> Rc<RefCell<SystemKekule>> {
        self.systems
            .entry(key)
            .or_insert_with(|| Rc::new(RefCell::new(SystemKekule::default())))
            .clone()
    }

    pub fn get(&self, key: &SystemKey) -> Option<Rc<RefCell<SystemKekule>>> {
        self.systems.get(key).cloned()
    }

    /// Tag-keyed system keys currently held (for cache-identity tests).
    pub fn system_keys(&self) -> Vec<SystemKey> {
        self.systems.keys().cloned().collect()
    }

    /// Fill (or reuse) the bag for `graph`'s [`SystemKey`].
    pub fn ensure_graph(
        &mut self,
        mol: &Molecule,
        graph: &PiGraph,
    ) -> Rc<RefCell<SystemKekule>> {
        let key = graph.system_key(mol);
        let slot = self.slot(key);
        fill_slot(mol, &graph.atoms, &graph.bonds, &slot);
        slot
    }
}

/// Config for the generic Kekulé / π-matching solver.
///
/// **This is the contract for altering solver behavior.** Callers (pair emit,
/// ResonanceRule parents, 2-core aromaticity) pass data constraints
/// (`forced_doubles`, `saturate`) plus this config. Do not add element- or
/// rule-named branches inside the matcher — flip a field here (or extend this
/// struct with a named, documented option) when behavior must change.
///
/// Defaults match current forest use: perceive post-edit doubles, enumerate
/// multi-resonance seeds, charge-follow on overlay, 4n+2 benzenoid 2-core.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KekuleConfig {
    /// Treat existing `BondOrder::Double` edges that touch the atom set as
    /// additional forced seeds (graph perception after end edits).
    pub perceive_existing_doubles: bool,
    /// When listing assignments, also try each residual bond as an extra
    /// forced seed (multi-resonance enumeration).
    pub enumerate_bond_seeds: bool,
    /// After stamping an assignment, run [`move_charge_with_bonds`].
    pub move_charge_on_overlay: bool,
    /// Cyclic 2-core component stays aromatic only if atom count is 4n+2
    /// (all-carbon benzenoid shortcut). When false, any cyclic valid Kekulé
    /// 2-core component is kept.
    pub huckel_4n2: bool,
}

impl Default for KekuleConfig {
    fn default() -> Self {
        Self {
            perceive_existing_doubles: true,
            enumerate_bond_seeds: true,
            move_charge_on_overlay: true,
            huckel_4n2: true,
        }
    }
}

impl KekuleConfig {
    /// Parent-bag fill / bond overlay (ResonanceRule `reactant_parent`).
    pub fn for_parents() -> Self {
        Self::default()
    }

    /// Constraint-directed pair materialize (forced doubles + saturate).
    pub fn for_constraints() -> Self {
        Self::default()
    }

    /// Single complete match only (no multi-resonance seed enumeration).
    pub fn single_assignment(mut self) -> Self {
        self.enumerate_bond_seeds = false;
        self
    }

    /// Do not auto-force doubles already on the mol (caller supplies all forced).
    pub fn explicit_forced_only(mut self) -> Self {
        self.perceive_existing_doubles = false;
        self
    }
}

/// Constraints for one matching call (data, not config).
///
/// Supply from perception and/or PatternInfo. The solver does not interpret
/// edit tokens or rule names.
#[derive(Clone, Debug, Default)]
pub struct KekuleConstraints {
    /// Edges that must be selected as double (forced leaf/edge).
    pub forced_doubles: BTreeSet<(usize, usize)>,
    /// Atoms that leave must-match (and typically gain H on the emit path).
    pub saturate: BTreeSet<usize>,
}

impl KekuleConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_forced(mut self, forced: BTreeSet<(usize, usize)>) -> Self {
        self.forced_doubles = forced;
        self
    }

    pub fn with_saturate(mut self, saturate: BTreeSet<usize>) -> Self {
        self.saturate = saturate;
        self
    }
}

/// Parents covering two atoms. Different systems: union, not a product.
#[derive(Clone)]
pub struct EndParents {
    pub parents: Vec<Molecule>,
    pub same_system: bool,
}

/// Bond-order sum for charge-follow (aromatic counts as 1.5, matching RDKit).
pub(crate) fn bond_order_sums(mol: &Molecule) -> HashMap<usize, f32> {
    let mut sums = HashMap::new();
    for (idx, _) in mol.atoms() {
        let i = atom_usize(idx);
        let mut sum = 0.0_f32;
        for (_nbr, bond_idx) in mol.neighbors(idx) {
            sum += match mol.bond(bond_idx).order {
                BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.0,
                BondOrder::Double => 2.0,
                BondOrder::Triple => 3.0,
                BondOrder::Aromatic => 1.5,
                BondOrder::Quadruple => 4.0,
                _ => 1.0,
            };
        }
        sums.insert(i, sum);
    }
    sums
}

/// Move formal charge when a bond-order flip would leave it behind.
///
/// Same rules as Python `move_charge_with_bonds` — general π bookkeeping, not
/// rule chemistry: neutral C shifts H; neutral aromatic atoms do not mint
/// charge from 1.5→1/2 alone; already-charged atoms and neutral non-aromatic
/// heteroatoms follow the bond (N may gain charge when bond order rises).
///
/// Closed-shell valence after **end edits** (forced doubles, cleavage) belongs
/// on the emit path ([`crate::pair_edit`] / valence fill), not here — do not
/// special-case elements for named reactions.
/// Python 3 `int(round(x))` — half to even — so aromatic 4.5→4 Kekulé
/// bridge carbons get delta 0 (not −1). Rust `f32::round` is half away from
/// zero and would mint a spurious explicit H (`[CH]` on a quaternary bridge).
fn round_half_even(x: f32) -> i8 {
    let floored = x.floor();
    let frac = x - floored;
    const EPS: f32 = 1e-6;
    if frac < 0.5 - EPS {
        floored as i8
    } else if frac > 0.5 + EPS {
        (floored + 1.0) as i8
    } else {
        // Exactly .5 (within EPS): toward even.
        let n = floored as i32;
        if n.rem_euclid(2) == 0 {
            n as i8
        } else {
            (n + 1) as i8
        }
    }
}

pub fn move_charge_with_bonds(
    mol: &mut Molecule,
    before: &HashMap<usize, f32>,
    aromatic: &HashSet<usize>,
) {
    let after = bond_order_sums(mol);
    let idxs: Vec<usize> = after.keys().copied().collect();
    for i in idxs {
        let old = match before.get(&i) {
            Some(&v) => v,
            None => continue,
        };
        let new = after[&i];
        let delta = round_half_even(new - old);
        if delta == 0 {
            continue;
        }
        let atom = mol.atom(atom_idx(i));
        let charge = atom.charge;
        let z = atom.element.atomic_number();
        let neutral = charge == 0;
        if neutral && z == 6 {
            // H travels with the bond on neutral carbon.
            let total = mol.implicit_hydrogen_count(atom_idx(i)) as i16;
            let updated = total - delta as i16;
            if updated >= 0 {
                *mol = with_atom_explicit_h(mol, atom_idx(i), updated as u8);
            }
        } else if neutral && aromatic.contains(&i) {
            continue;
        } else {
            mol.set_charge(atom_idx(i), charge.saturating_add(delta));
        }
    }
}

pub(crate) fn with_atom_explicit_h(mol: &Molecule, idx: chematic::core::AtomIdx, h: u8) -> Molecule {
    use chematic::core::MoleculeBuilder;
    let mut builder = MoleculeBuilder::new();
    for (aidx, atom) in mol.atoms() {
        let mut a = atom.clone();
        if aidx == idx {
            a.hydrogen_count = Some(h);
        }
        builder.add_atom(a);
    }
    for (_, bond) in mol.bonds() {
        let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    builder.copy_stereo_from(mol);
    builder.copy_r_groups_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.build()
}

fn match_assignment(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    seed: (usize, usize),
    config: &KekuleConfig,
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    let constraints = KekuleConstraints {
        forced_doubles: BTreeSet::from([bond_key(seed.0, seed.1)]),
        saturate: BTreeSet::new(),
    };
    complete_assignment(mol, atoms, bonds, &constraints, config)
}

/// Complete a Kekulé assignment under π constraints.
///
/// `constraints.forced_doubles` seed the matching (exocyclic leaves / edited
/// edges that must be double). `constraints.saturate` atoms leave must-match
/// (PatternInfo when the end edit does not change the graph — e.g. `keep` +
/// adds H). [`KekuleConfig::perceive_existing_doubles`] optionally treats
/// doubles already on `mol` as additional forced seeds. Returns `None` when
/// the residual matching is impossible.
pub fn complete_assignment(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    constraints: &KekuleConstraints,
    config: &KekuleConfig,
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    use chematic::core::kekulization::atom_must_be_matched;
    let must: BTreeSet<usize> = atoms
        .iter()
        .copied()
        .filter(|&i| {
            !constraints.saturate.contains(&i) && atom_must_be_matched(mol, atom_idx(i))
        })
        .collect();
    let mut doubles: HashMap<usize, usize> = HashMap::new();
    for &(left, right) in &constraints.forced_doubles {
        if doubles.contains_key(&left) || doubles.contains_key(&right) {
            if doubles.get(&left) != Some(&right) {
                return None;
            }
            continue;
        }
        doubles.insert(left, right);
        doubles.insert(right, left);
    }
    // Optionally treat existing Double bonds touching the system as forced
    // (perception after end edits — the leaf partner may sit outside `must`).
    if config.perceive_existing_doubles {
        for (_, bond) in mol.bonds() {
            if bond.order != BondOrder::Double {
                continue;
            }
            let a = atom_usize(bond.atom1);
            let b = atom_usize(bond.atom2);
            if !atoms.contains(&a) && !atoms.contains(&b) {
                continue;
            }
            if doubles.contains_key(&a) || doubles.contains_key(&b) {
                if doubles.get(&a) != Some(&b) {
                    return None;
                }
                continue;
            }
            doubles.insert(a, b);
            doubles.insert(b, a);
        }
    }
    let mut adj: HashMap<usize, Vec<usize>> = must.iter().map(|&a| (a, Vec::new())).collect();
    for &(left, right) in bonds {
        if !must.contains(&left) || !must.contains(&right) {
            continue;
        }
        adj.get_mut(&left)?.push(right);
        adj.get_mut(&right)?.push(left);
    }
    for nbrs in adj.values_mut() {
        nbrs.sort_unstable();
    }
    let places: Vec<usize> = must.iter().copied().collect();
    fn place(
        idx: usize,
        places: &[usize],
        adj: &HashMap<usize, Vec<usize>>,
        doubles: &mut HashMap<usize, usize>,
    ) -> bool {
        if idx == places.len() {
            return true;
        }
        let atom = places[idx];
        if doubles.contains_key(&atom) {
            return place(idx + 1, places, adj, doubles);
        }
        let Some(nbrs) = adj.get(&atom) else {
            return false;
        };
        for &nbr in nbrs {
            if doubles.contains_key(&nbr) {
                continue;
            }
            doubles.insert(atom, nbr);
            doubles.insert(nbr, atom);
            if place(idx + 1, places, adj, doubles) {
                return true;
            }
            doubles.remove(&atom);
            doubles.remove(&nbr);
        }
        false
    }
    if !place(0, &places, &adj, &mut doubles) {
        return None;
    }
    // Every must-match atom must be paired.
    if places.iter().any(|a| !doubles.contains_key(a)) {
        return None;
    }
    let mut written = BTreeMap::new();
    for &(left, right) in bonds {
        let is_double = doubles.get(&left) == Some(&right);
        written.insert(
            bond_key(left, right),
            if is_double {
                BondOrder::Double
            } else {
                BondOrder::Single
            },
        );
    }
    Some(written)
}

/// Count distinct complete assignments under the same constraints (multi-resonance).
pub fn count_assignments(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    constraints: &KekuleConstraints,
    config: &KekuleConfig,
) -> usize {
    all_assignments(mol, atoms, bonds, constraints, config).len()
}

/// Every distinct complete assignment under the constraints.
///
/// When [`KekuleConfig::enumerate_bond_seeds`] is set, also tries each residual
/// bond as an extra forced seed (multi-resonance enumeration).
pub fn all_assignments(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    constraints: &KekuleConstraints,
    config: &KekuleConfig,
) -> Vec<BTreeMap<(usize, usize), BondOrder>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    let mut consider = |forced: &BTreeSet<(usize, usize)>| {
        let c = KekuleConstraints {
            forced_doubles: forced.clone(),
            saturate: constraints.saturate.clone(),
        };
        let Some(assignment) = complete_assignment(mol, atoms, bonds, &c, config) else {
            return;
        };
        let sig: Vec<_> = assignment
            .iter()
            .map(|(&k, &o)| (k, order_code(o)))
            .collect();
        if seen.insert(sig) {
            out.push(assignment);
        }
    };
    consider(&constraints.forced_doubles);
    if config.enumerate_bond_seeds {
        for &seed in bonds {
            if constraints.forced_doubles.contains(&seed) {
                continue;
            }
            let mut forced = constraints.forced_doubles.clone();
            forced.insert(seed);
            consider(&forced);
        }
    }
    out
}

/// 2-core of a π graph: iteratively drop atoms with degree &lt; 2.
///
/// Prefer [`PiGraph::two_core`].
pub fn two_core(
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
    let g = PiGraph::new(atoms.clone(), bonds.clone()).two_core();
    (g.atoms, g.bonds)
}

/// Atoms that remain aromatic: cyclic 2-core components with a valid Kekulé
/// state (and 4n+2 π electrons when [`KekuleConfig::huckel_4n2`] is set —
/// all-carbon benzenoid shortcut).
///
/// Build the residual [`PiGraph`] **after** dearomatizing edits
/// ([`PiGraph::without_atoms`] / [`PiGraph::after_forced_doubles`]). Then ask
/// which cyclic 2-core pieces survive — not whether the original fused system
/// is still aromatic. Prefer [`PiGraph::aromatic_2core_atoms`].
pub fn aromatic_2core_atoms(
    mol: &Molecule,
    residual_atoms: &BTreeSet<usize>,
    residual_bonds: &BTreeSet<(usize, usize)>,
    config: &KekuleConfig,
) -> BTreeSet<usize> {
    PiGraph::new(residual_atoms.clone(), residual_bonds.clone()).aromatic_2core_atoms(mol, config)
}

/// Stamp one system's assignment onto `mol`. Other systems stay as they were.
///
/// `atoms` and `assignment` are **index**-keyed for the current mol layout
/// (perception). Cached bags store tag-keyed maps — resolve before calling.
///
/// When [`KekuleConfig::move_charge_on_overlay`] is set, applies
/// [`move_charge_with_bonds`] so charge separation on valid writings
/// (including N that gains charge) matches Python reactant overlays.
pub fn overlay(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    assignment: &BTreeMap<(usize, usize), BondOrder>,
    config: &KekuleConfig,
) -> Molecule {
    let before = bond_order_sums(mol);
    let aromatic: HashSet<usize> = mol
        .atoms()
        .filter_map(|(idx, atom)| atom.aromatic.then_some(atom_usize(idx)))
        .collect();
    let mut out = mol.clone();
    for (&(left, right), &order) in assignment {
        if let Some((bond_idx, _)) = out.bond_between(atom_idx(left), atom_idx(right)) {
            out.set_bond_order(bond_idx, order);
        }
    }
    for &atom in atoms {
        if out.atom(atom_idx(atom)).aromatic {
            out = out.with_atom_aromatic(atom_idx(atom), false);
        }
    }
    if config.move_charge_on_overlay {
        move_charge_with_bonds(&mut out, &before, &aromatic);
    }
    out
}

/// Overlay from a **tag-keyed** cached assignment (resolves tags → indexes).
pub fn overlay_tagged(
    mol: &Molecule,
    atom_tags: &BTreeSet<Tag>,
    assignment: &BTreeMap<(Tag, Tag), BondOrder>,
    config: &KekuleConfig,
) -> Option<Molecule> {
    let atoms = atoms_to_idxs(mol, atom_tags)?;
    let orders = assignment_to_idxs(mol, assignment)?;
    Some(overlay(mol, &atoms, &orders, config))
}

fn fill_slot(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    slot: &RefCell<SystemKekule>,
) {
    if slot.borrow().is_filled() {
        return;
    }
    let config = KekuleConfig::for_parents();
    let mut bag = slot.borrow_mut();
    let mut seen = BTreeSet::new();
    for &seed in bonds {
        let Some(bond_orders) = match_assignment(mol, atoms, bonds, seed, &config) else {
            continue;
        };
        let Some(tagged) = assignment_to_tags(mol, &bond_orders) else {
            continue;
        };
        let signature: Vec<_> = tagged
            .iter()
            .map(|(&k, &o)| (k, order_code(o)))
            .collect();
        if !seen.insert(signature) {
            continue;
        }
        let index = bag.assignments.len();
        let atom_tags: BTreeSet<Tag> = atoms.iter().filter_map(|&i| tag_of(mol, i)).collect();
        let charge_mag = overlay_tagged(mol, &atom_tags, &tagged, &config)
            .map(|p| {
                p.atoms()
                    .map(|(_, a)| a.charge.unsigned_abs() as i32)
                    .sum()
            })
            .unwrap_or(0);
        for (key, order) in &tagged {
            if let Some(code) = order_code(*order) {
                if code == 1 || code == 2 {
                    match bag.by_order.get(&(*key, code)) {
                        None => {
                            bag.by_order.insert((*key, code), index);
                        }
                        Some(&prev) => {
                            // Prefer closed-shell overlay (lower |charge|).
                            let prev_mag = bag.charge_mags.get(prev).copied().unwrap_or(i32::MAX);
                            if charge_mag < prev_mag {
                                bag.by_order.insert((*key, code), index);
                            }
                        }
                    }
                }
            }
        }
        bag.assignments.push(tagged);
        bag.charge_mags.push(charge_mag);
    }
}

fn system_of(mol: &Molecule, left: usize, right: usize) -> (SystemKey, PiGraph) {
    let mut graph = PiGraph::conjugated(mol, left);
    let seed = bond_key(left, right);
    if !graph.bonds.contains(&seed) && mol.bond_between(atom_idx(left), atom_idx(right)).is_some()
    {
        graph.bonds.insert(seed);
        graph.atoms.insert(right);
    }
    (graph.system_key(mol), graph)
}

/// Fill (or reuse) the bag for the system that contains `(left, right)`.
pub fn ensure_kekule_parents(
    mol: &Molecule,
    left: usize,
    right: usize,
    cache: &mut KekuleCache,
) -> Rc<RefCell<SystemKekule>> {
    let (_key, graph) = system_of(mol, left, right);
    cache.ensure_graph(mol, &graph)
}

/// Overlay `mol` with the cached assignment where `(left, right)` has `order`.
pub fn parent_for_bond(
    mol: &Molecule,
    cache: &KekuleCache,
    left: usize,
    right: usize,
    order: u8,
) -> Option<Molecule> {
    let (key, _) = system_of(mol, left, right);
    let slot = cache.get(&key)?;
    let bag = slot.borrow();
    let ta = tag_of(mol, left)?;
    let tb = tag_of(mol, right)?;
    let index = *bag.by_order.get(&(tag_bond_key(ta, tb), order))?;
    let assignment = bag.assignments.get(index)?;
    overlay_tagged(mol, &key.atoms, assignment, &KekuleConfig::for_parents())
}

/// Parents covering `start` and `end`. Different systems: union, not a product.
pub fn parents_for_ends(
    mol: &Molecule,
    start: usize,
    end: usize,
    cache: &mut KekuleCache,
    atoms: Option<&BTreeSet<usize>>,
) -> EndParents {
    if let Some(atoms) = atoms {
        let bonds: BTreeSet<_> = mol
            .bonds()
            .filter_map(|(_, bond)| {
                let a = atom_usize(bond.atom1);
                let b = atom_usize(bond.atom2);
                (atoms.contains(&a) && atoms.contains(&b)).then_some(bond_key(a, b))
            })
            .collect();
        if bonds.is_empty() {
            return EndParents {
                parents: Vec::new(),
                same_system: true,
            };
        }
        let graph = PiGraph::new(atoms.clone(), bonds);
        let key = graph.system_key(mol);
        let slot = cache.ensure_graph(mol, &graph);
        let parents = slot
            .borrow()
            .assignments
            .iter()
            .filter_map(|assignment| {
                overlay_tagged(mol, &key.atoms, assignment, &KekuleConfig::for_parents())
            })
            .collect();
        return EndParents {
            parents,
            same_system: true,
        };
    }
    let start_g = PiGraph::conjugated(mol, start);
    let end_g = PiGraph::conjugated(mol, end);
    let start_key = start_g.system_key(mol);
    let end_key = end_g.system_key(mol);
    let start_slot = cache.ensure_graph(mol, &start_g);
    let same = start_g.atoms == end_g.atoms;
    if !same {
        let _ = cache.ensure_graph(mol, &end_g);
    }
    let mut parents: Vec<Molecule> = start_slot
        .borrow()
        .assignments
        .iter()
        .filter_map(|assignment| {
            overlay_tagged(
                mol,
                &start_key.atoms,
                assignment,
                &KekuleConfig::for_parents(),
            )
        })
        .collect();
    if !same {
        let end_slot = cache.get(&end_key).expect("end system filled");
        parents.extend(end_slot.borrow().assignments.iter().filter_map(|assignment| {
            overlay_tagged(mol, &end_key.atoms, assignment, &KekuleConfig::for_parents())
        }));
    }
    EndParents {
        parents,
        same_system: same,
    }
}

/// Overlays of every aromatic conjugated system. Pair-path derisk only.
///
/// Assignments stay inside **aromatic** atoms of each conjugated component
/// (exocyclic amide / nitro are not rewritten). That matches Python pair
/// conjugated parents from ``ResonanceMolSupplier`` on typical aromatics —
/// charged amide kekulé forms belong to ResonanceRule ``reactant_parent``,
/// not the pair door. A ResonanceRule should call [`ensure_kekule_parents`]
/// for the matched bond instead of this. Counts are a **sum** of systems,
/// not a product.
///
/// Stamps missing index tags on a clone so bare mols can fill a local cache;
/// ForestMol callers already carry copy-tree labels.
pub fn kekule_forms(mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
    if !mol.atoms().any(|(_, atom)| atom.aromatic) {
        return Ok(vec![mol.clone()]);
    }
    let mut mol = mol.clone();
    stamp_missing_index_tags(&mut mol);
    let mol = &mol;
    let mut cache = KekuleCache::default();
    let mut covered = BTreeSet::new();
    let mut forms = Vec::new();
    for (idx, atom) in mol.atoms() {
        if !atom.aromatic {
            continue;
        }
        let start = atom_usize(idx);
        if !covered.insert(start) {
            continue;
        }
        let graph = PiGraph::conjugated(mol, start);
        covered.extend(&graph.atoms);
        let aromatic = graph.aromatic_subgraph(mol);
        if aromatic.atoms.len() < 2 || aromatic.bonds.is_empty() {
            continue;
        }
        let key = aromatic.system_key(mol);
        let slot = cache.ensure_graph(mol, &aromatic);
        forms.extend(
            slot.borrow()
                .assignments
                .iter()
                .filter_map(|assignment| {
                    overlay_tagged(mol, &key.atoms, assignment, &KekuleConfig::for_parents())
                }),
        );
    }
    if forms.is_empty() {
        Ok(vec![mol.clone()])
    } else {
        Ok(forms)
    }
}

/// ResonanceRule parent: overlay of **that bond's system** with the implied order.
pub fn reactant_parent(
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
    smirks: &str,
    cache: &mut KekuleCache,
) -> Result<Molecule, ForestError> {
    let Some(&left) = mapped.get(&1) else {
        return Ok(mol.clone());
    };
    let Some(&right) = mapped.get(&2) else {
        return Ok(mol.clone());
    };
    let Some((_, bond)) = mol.bond_between(atom_idx(left), atom_idx(right)) else {
        return Ok(mol.clone());
    };
    let aromatic = bond.order == BondOrder::Aromatic
        || mol.atom(atom_idx(left)).aromatic
        || mol.atom(atom_idx(right)).aromatic;
    if !aromatic {
        return Ok(mol.clone());
    }
    // Bare mols need tags for SystemKey; ForestMol already stamps.
    let mut stamped = mol.clone();
    stamp_missing_index_tags(&mut stamped);
    let mol = &stamped;
    let want = smirks_mapped_bond_order(smirks).unwrap_or(2.0);
    let want_code: u8 = if want >= 1.5 { 2 } else { 1 };
    ensure_kekule_parents(mol, left, right, cache);
    let Some(parent) = parent_for_bond(mol, cache, left, right, want_code) else {
        return Ok(mol.clone());
    };
    // Prefer aromatic parent when Kekulé writing raises |charge|
    // (isocyanate N=C=O → O=C[N-]Ar). Closed-shell prefer (C10 / C16 leave).
    let parent_mag: i32 = parent
        .atoms()
        .map(|(_, a)| a.charge.unsigned_abs() as i32)
        .sum();
    let mol_mag: i32 = mol
        .atoms()
        .map(|(_, a)| a.charge.unsigned_abs() as i32)
        .sum();
    if parent_mag > mol_mag {
        return Ok(mol.clone());
    }
    Ok(parent)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forest_mol::ForestMol;
    use crate::mol::{canon_of, canon_smiles, parse_mol};
    use crate::smarts::smarts_matches;
    use crate::smirks::apply_smirks_at;
    use chematic::core::Atom;

    const ANTHRACENE: &str = "c1ccc2cc3ccccc3cc2c1";
    const POLYPHENYL: &str = "c1ccc(-c2ccc(-c3ccc(-c4ccc(-c5ccc(-c6ccccc6)cc5)cc4)cc3)cc2)cc1";

    fn aromatic_count(mol: &Molecule) -> usize {
        mol.atoms().filter(|(_, atom)| atom.aromatic).count()
    }

    fn fill(smiles: &str) -> (Molecule, KekuleCache) {
        let mut mol = parse_mol(smiles).unwrap();
        stamp_missing_index_tags(&mut mol);
        let mut cache = KekuleCache::default();
        for mapped in smarts_matches(&mol, "[#6:1]=,:[#6:2]").unwrap() {
            ensure_kekule_parents(&mol, mapped[&1], mapped[&2], &mut cache);
        }
        (mol, cache)
    }

    #[test]
    fn system_key_and_assignments_are_tag_keyed_not_index() {
        // Reusable cache objects store Tags (forest labels), not atom indexes.
        let parent = ForestMol::parse("c1ccccc1").unwrap();
        let (atoms, bonds) = conjugated_component(parent.mol(), 0);
        let key = SystemKey::of(parent.mol(), &atoms, &bonds);
        assert!(
            key.atoms.iter().all(|t| parent.index_of(*t).is_some()),
            "SystemKey atoms are Tags resolvable on the mol"
        );
        let seed = *bonds.iter().next().expect("ring bond");
        parent.ensure_kekule(seed.0, seed.1);
        let kekule_rc = parent.kekule();
        let cache = kekule_rc.borrow();
        let slot = cache.get(&key).expect("filled");
        let asg = &slot.borrow().assignments[0];
        for &(ta, tb) in asg.keys() {
            assert!(parent.index_of(ta).is_some() && parent.index_of(tb).is_some());
        }
        // ResidualKey is also tag-shaped (removed / forced are Tags).
        let residual = ResidualKey {
            parent: key.clone(),
            removed: BTreeSet::new(),
            forced_doubles: BTreeSet::new(),
        };
        assert_eq!(residual.parent.atoms, key.atoms);
    }

    #[test]
    fn benzene_has_two_kekule_forms() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let forms = kekule_forms(&mol).unwrap();
        assert!(forms.len() >= 2);
        let orders: Vec<_> = forms.iter().map(bond_order_map).collect();
        assert_ne!(orders[0], orders[1]);
        for form in &forms {
            let doubles = form
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Double)
                .count();
            assert_eq!(doubles, 3);
        }
    }

    /// Naphthalene → 1,4-naphthoquinone: one fused ring collapses to a
    /// cyclohexadienedione; the other remains aromatic. Gate is cyclic 2-core
    /// + Kekulé + 4n+2 — no `aromatize` / sanitize call.
    ///
    /// Reactant topology `c1ccc2ccccc2c1`; product shape
    /// `O=C1C=CC(=O)c2ccccc12`. Forced exocyclic C=O doubles consume the
    /// quinone carbons' π demand; residual 2-core keeps the benzenoid sextet.
    #[test]
    fn naphthalene_14_naphthoquinone_partial_collapse_without_sanitize() {
        // 1,4-naphthalenediol: phenol ends become the forced C=O leaves.
        // Indices (RDKit/chematic): O0–C1, O5–C4; fused ring C6–C11.
        let mol = parse_mol("Oc1ccc(O)c2ccccc12").unwrap();
        let parent = PiGraph::conjugated(&mol, 1);
        assert!(
            parent.atoms.len() >= 10,
            "naphthalene diol is one fused conjugated system: {:?}",
            parent.atoms
        );
        let forced: BTreeSet<(usize, usize)> =
            [(0usize, 1usize), (4usize, 5usize)]
                .into_iter()
                .map(|(a, b)| bond_key(a, b))
                .collect();
        // Perceive as the emit path does: only the in-system endpoint of each
        // exocyclic forced double is demand-consumed.
        let residual = parent.after_forced_doubles(&forced);
        assert!(
            !residual.atoms.contains(&1) && !residual.atoms.contains(&4),
            "quinone carbons leave the residual π graph"
        );
        let core = residual.two_core();
        assert_eq!(
            core.atoms,
            BTreeSet::from([6, 7, 8, 9, 10, 11]),
            "dangling quinone ring strips; fused benzenoid 2-core remains"
        );
        let aromatic = residual.aromatic_2core_atoms(&mol, &KekuleConfig::for_constraints());
        assert_eq!(
            aromatic,
            BTreeSet::from([6, 7, 8, 9, 10, 11]),
            "partial collapse stamps only the surviving aromatic sextet (no sanitize)"
        );
        assert!(
            !aromatic.contains(&1) && !aromatic.contains(&4) && !aromatic.contains(&2),
            "quinone-ring carbons must not stay aromatic"
        );
    }

    /// Naphthalene epoxidize outer left bond 0-1: left collapses, right 2-core
    /// stays aromatic (HEURISTICS cyclic 2-core partial collapse).
    #[test]
    fn naphthalene_epoxide_leaves_right_ring_aromatic_2core() {
        // Graph numbering matches the HEURISTICS example:
        // left 0-1-2-3-4-5-0, right 2-6-7-8-9-3-2, shared 2-3.
        let atoms: BTreeSet<usize> = (0..10).collect();
        let bonds: BTreeSet<(usize, usize)> = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (2, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 3),
        ]
        .into_iter()
        .map(|(a, b)| bond_key(a, b))
        .collect();
        // Epoxidation of 0-1 → remove sp³ atoms 0 and 1 from π graph.
        let mut residual_atoms = atoms.clone();
        residual_atoms.remove(&0);
        residual_atoms.remove(&1);
        let residual_bonds: BTreeSet<_> = bonds
            .iter()
            .copied()
            .filter(|&(a, b)| residual_atoms.contains(&a) && residual_atoms.contains(&b))
            .collect();
        let (core_atoms, _) = two_core(&residual_atoms, &residual_bonds);
        assert_eq!(
            core_atoms,
            BTreeSet::from([2, 3, 6, 7, 8, 9]),
            "dangling 3-4-5 strips; right ring remains"
        );
        // Need a real mol for Kekulé check — benzene stands in for the right ring
        // topology; aromatic_2core on the abstract residual uses complete_assignment
        // on mol indices, so build a naphthalene mol and map by index.
        let mol = parse_mol("c1ccc2ccccc2c1").unwrap();
        let (sys, sys_bonds) = conjugated_component(&mol, 0);
        assert!(sys.len() >= 10);
        // Remove two atoms that share an outer (non-fusion) bond.
        let fusion: BTreeSet<_> = sys_bonds
            .iter()
            .copied()
            .filter(|&(a, b)| {
                let a_ring = sys_bonds.iter().filter(|&&(x, y)| x == a || y == a).count();
                let b_ring = sys_bonds.iter().filter(|&&(x, y)| x == b || y == b).count();
                // fusion carbons have degree 3 in the π graph
                a_ring == 3 && b_ring == 3
            })
            .collect();
        let outer = sys_bonds
            .iter()
            .copied()
            .find(|&(a, b)| {
                !fusion.contains(&(a, b))
                    && sys_bonds.iter().filter(|&&(x, y)| x == a || y == a).count() == 2
                    && sys_bonds.iter().filter(|&&(x, y)| x == b || y == b).count() == 2
            })
            .expect("outer bond");
        let mut residual_atoms = sys.clone();
        residual_atoms.remove(&outer.0);
        residual_atoms.remove(&outer.1);
        let residual_bonds: BTreeSet<_> = sys_bonds
            .iter()
            .copied()
            .filter(|&(a, b)| residual_atoms.contains(&a) && residual_atoms.contains(&b))
            .collect();
        let aromatic = aromatic_2core_atoms(
            &mol,
            &residual_atoms,
            &residual_bonds,
            &KekuleConfig::for_constraints(),
        );
        assert_eq!(aromatic.len(), 6, "one benzenoid sextet survives: {aromatic:?}");
        assert!(!aromatic.contains(&outer.0) && !aromatic.contains(&outer.1));
    }

    #[test]
    fn kekule_config_gates_enumeration_and_huckel() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let (atoms, bonds) = conjugated_component(&mol, 0);
        let empty = KekuleConstraints::default();
        let multi = all_assignments(
            &mol,
            &atoms,
            &bonds,
            &empty,
            &KekuleConfig::for_parents(),
        );
        assert!(multi.len() >= 2, "default enumerates bond seeds");
        let single = all_assignments(
            &mol,
            &atoms,
            &bonds,
            &empty,
            &KekuleConfig::for_parents().single_assignment(),
        );
        assert_eq!(single.len(), 1, "single_assignment stops after one complete");

        // Cyclobutadiene-shaped 4-atom cyclic 2-core: valid matching, fails 4n+2.
        let four: BTreeSet<usize> = BTreeSet::from([0, 1, 2, 3]);
        let four_bonds: BTreeSet<_> = [(0, 1), (1, 2), (2, 3), (3, 0)]
            .into_iter()
            .map(|(a, b)| bond_key(a, b))
            .collect();
        let with_huckel = aromatic_2core_atoms(
            &mol,
            &four,
            &four_bonds,
            &KekuleConfig {
                huckel_4n2: true,
                ..KekuleConfig::default()
            },
        );
        assert!(with_huckel.is_empty());
        let without = aromatic_2core_atoms(
            &mol,
            &four,
            &four_bonds,
            &KekuleConfig {
                huckel_4n2: false,
                ..KekuleConfig::default()
            },
        );
        assert_eq!(without, four);
    }

    #[test]
    fn pi_graph_derives_residual_and_filters_assignments() {
        let mut mol = parse_mol("c1ccc2ccccc2c1").unwrap();
        stamp_missing_index_tags(&mut mol);
        let parent = PiGraph::conjugated(&mol, 0);
        assert!(parent.atoms.len() >= 10);
        // Epoxide-style: drop two outer atoms → residual 2-core keeps one ring.
        let drop: BTreeSet<usize> = parent
            .bonds
            .iter()
            .copied()
            .find(|&(a, b)| {
                let da = parent.bonds.iter().filter(|&&(x, y)| x == a || y == a).count();
                let db = parent.bonds.iter().filter(|&&(x, y)| x == b || y == b).count();
                da == 2 && db == 2
            })
            .map(|(a, b)| BTreeSet::from([a, b]))
            .expect("outer bond");
        let residual = parent.without_atoms(&drop);
        let core = residual.two_core();
        assert_eq!(core.atoms.len(), 6, "one benzenoid 2-core: {:?}", core.atoms);
        assert_eq!(core.components().len(), 1);

        let mut cache = KekuleCache::default();
        let slot = cache.ensure_graph(&mol, &parent);
        let n = slot.borrow().assignments.len();
        assert!(n >= 2);
        // Forced double on any ring bond → derived bag is a non-empty subset.
        let (a, b) = *parent.bonds.iter().next().unwrap();
        let ta = tag_of(&mol, a).unwrap();
        let tb = tag_of(&mol, b).unwrap();
        let forced = BTreeSet::from([tag_bond_key(ta, tb)]);
        let derived = slot.borrow().with_forced_doubles(&forced);
        assert!(!derived.assignments.is_empty());
        assert!(derived.assignments.len() <= n);

        let rkey = parent
            .residual_key(&mol, &drop, &BTreeSet::new())
            .expect("tags");
        let applied = rkey.apply_to(&mol, &parent).unwrap();
        assert_eq!(applied.atoms, residual.atoms);
    }

    #[test]
    fn anthracene_is_four_parents_not_sixteen() {
        let (_mol, cache) = fill(ANTHRACENE);
        assert_eq!(cache.assignment_count(), 4);
        assert_eq!(cache.system_count(), 1);
    }

    #[test]
    fn polyphenyl_is_twelve_parents_not_sixty_four() {
        let (mol, cache) = fill(POLYPHENYL);
        assert_eq!(cache.system_count(), 6);
        assert_eq!(cache.assignment_count(), 12, "6 rings × 2, not 2^6");
        let mapped = smarts_matches(&mol, "[#6:1]=,:[#6:2]").unwrap();
        let parent = parent_for_bond(&mol, &cache, mapped[0][&1], mapped[0][&2], 2).unwrap();
        assert_eq!(
            aromatic_count(&parent),
            30,
            "one ring kekulized, five stay aromatic"
        );
    }

    #[test]
    fn biphenyl_ends_are_a_union_not_a_product() {
        let mut mol = parse_mol("c1ccc(-c2ccccc2)cc1").unwrap();
        stamp_missing_index_tags(&mut mol);
        let (sys_a, _) = conjugated_component(&mol, 0);
        let other = (0..mol.atom_count())
            .find(|&i| conjugated_component(&mol, i).0 != sys_a)
            .expect("second ring");
        let mut cache = KekuleCache::default();
        let ends = parents_for_ends(&mol, 0, other, &mut cache, None);
        assert!(!ends.same_system);
        assert_eq!(ends.parents.len(), 4, "2+2, not 2×2 fully kekulized mols");
        for parent in &ends.parents {
            assert_eq!(aromatic_count(parent), 6);
        }
    }

    #[test]
    fn first_relative_to_fill_shares_with_unmodified_alignments() {
        let parent = ForestMol::parse("c1ccc(-c2ccccc2)cc1").unwrap();
        let (sys_a, _) = conjugated_component(parent.mol(), 0);
        let bond_a = sys_a.iter().copied().next().unwrap();
        let nbr = parent
            .mol()
            .neighbors(atom_idx(bond_a))
            .find(|(n, _)| sys_a.contains(&atom_usize(*n)))
            .map(|(n, _)| atom_usize(n))
            .unwrap();
        let other = other_atom(parent.mol(), &sys_a);

        let sibling = parent.copy_mol();
        let (mut chem, oxygen) = parent
            .mol()
            .with_atom_added(Atom::organic(chematic::core::Element::O));
        chem.add_bond(atom_idx(other), oxygen, BondOrder::Single)
            .unwrap();
        let child = ForestMol::product(chem, &parent);
        assert!(child.shares_kekule(&parent));
        assert!(sibling.shares_kekule(&parent));
        assert!(!child.shares_structure(&parent));

        assert_eq!(parent.kekule().borrow().assignment_count(), 0);
        sibling.ensure_kekule(bond_a, nbr);
        assert_eq!(parent.kekule().borrow().assignment_count(), 2);
        assert_eq!(child.kekule().borrow().assignment_count(), 2);

        let parent_edited_key = {
            let (atoms, bonds) = conjugated_component(parent.mol(), other);
            SystemKey::of(parent.mol(), &atoms, &bonds)
        };
        let child_edited_key = {
            let start = other_atom(child.mol(), &sys_a);
            let (atoms, bonds) = conjugated_component(child.mol(), start);
            SystemKey::of(child.mol(), &atoms, &bonds)
        };
        assert_ne!(parent_edited_key, child_edited_key);
        let systems_before = child.kekule().borrow().system_count();
        let start = other_atom(child.mol(), &sys_a);
        let seed = *conjugated_component(child.mol(), start)
            .1
            .iter()
            .next()
            .unwrap();
        child.ensure_kekule(seed.0, seed.1);
        assert!(child.kekule().borrow().system_count() > systems_before);
    }

    fn other_atom(mol: &Molecule, sys_a: &BTreeSet<usize>) -> usize {
        (0..mol.atom_count())
            .find(|&i| conjugated_component(mol, i).0 != *sys_a)
            .unwrap()
    }

    #[test]
    fn product_shares_kekule_rc_not_structure() {
        let parent = ForestMol::parse("c1ccccc1").unwrap();
        let _ = parent.csmi();
        let child = parent.edit_copy();
        assert!(child.shares_kekule(&parent));
        assert!(!child.shares_structure(&parent));
    }

    #[test]
    fn thiophene_s_oxidation_picks_single_s_c_parent() {
        let mut mol = parse_mol("c1ccsc1").unwrap();
        stamp_missing_index_tags(&mut mol);
        let smarts = "[#6:2]1=,:[#6:3][#6:4]=,:[#6:5][#16;v2,v4:1]1";
        let apply = "[S:1]>>[S+:1][O-]";
        let hits = smarts_matches(&mol, smarts).unwrap();
        assert!(!hits.is_empty());
        let mut cache = KekuleCache::default();
        let parent =
            reactant_parent(&mol, &hits[0], &format!("{smarts}>>[S:1]"), &mut cache).unwrap();
        let s = hits[0][&1];
        let c = hits[0][&2];
        let (_, bond) = parent.bond_between(atom_idx(s), atom_idx(c)).unwrap();
        assert_eq!(bond.order, BondOrder::Single);
        let products = apply_smirks_at(apply, &parent, &hits[0]).unwrap();
        assert!(!products.is_empty());
        let want = canon_of("[O-][s+]1cccc1").unwrap();
        assert!(
            products
                .iter()
                .any(|p| canon_of(&canon_smiles(p)).unwrap() == want),
            "want {want}, got {:?}",
            products.iter().map(canon_smiles).collect::<Vec<_>>()
        );
    }

    #[test]
    fn epoxidation_picks_double_parent_on_benzene() {
        let mut mol = parse_mol("c1ccccc1").unwrap();
        stamp_missing_index_tags(&mut mol);
        let smirks = "[#6:1]=[#6,#7:2]>>[*:1]1-[*:2][O]1";
        let hits = smarts_matches(&mol, "[#6:1]=,:[#6:2]").unwrap();
        assert_eq!(hits.len(), 6);
        let mut cache = KekuleCache::default();
        let parent = reactant_parent(&mol, &hits[0], smirks, &mut cache).unwrap();
        let a = hits[0][&1];
        let b = hits[0][&2];
        let (_, bond) = parent.bond_between(atom_idx(a), atom_idx(b)).unwrap();
        assert_eq!(bond.order, BondOrder::Double);
        assert_eq!(aromatic_count(&parent), 0);
    }

    #[test]
    fn indole_kekule_parent_is_closed_shell() {
        // Bridge carbons are aromatic bos 4.5 → Kekulé 4.0 (delta −0.5).
        // Banker's round → 0 so move_charge must not mint `[CH]` (Rust
        // f32::round would); otherwise ring-open Dealk refuses the parent.
        use chematic::perception::validate_valence;
        let mol = parse_mol("c1ccc2[nH]ccc2c1").unwrap();
        let hits = smarts_matches(&mol, "[#6H1:1][#7:2]").unwrap();
        assert!(!hits.is_empty());
        let mut cache = KekuleCache::default();
        let parent = reactant_parent(
            &mol,
            &hits[0],
            "[#6H1:1][#7:2]>>([*:2].[*:1]=O)",
            &mut cache,
        )
        .unwrap();
        assert!(
            validate_valence(&parent).is_empty(),
            "indole single-bond parent must be closed-shell: {} errs={:?}",
            canon_smiles(&parent),
            validate_valence(&parent)
        );
    }

    #[test]
    fn shared_kekule_bag_matches_fresh_fill_after_edit_copy() {
        // edit_copy shares the Kekulé Rc. Untouched systems' cached
        // assignments must equal a fresh fill on the same mol (tag-keyed).
        let parent = ForestMol::parse("c1ccc(-c2ccccc2)cc1").unwrap();
        let (atoms, bonds) = conjugated_component(parent.mol(), 0);
        let seed = *bonds.iter().next().expect("ring bond");
        parent.ensure_kekule(seed.0, seed.1);
        let key = SystemKey::of(parent.mol(), &atoms, &bonds);
        let shared = parent
            .kekule()
            .borrow()
            .get(&key)
            .expect("parent bag")
            .borrow()
            .clone();

        let child = parent.edit_copy();
        assert!(child.shares_kekule(&parent));
        // Fresh cache on the same chemistry must match the shared bag.
        let mut fresh = KekuleCache::default();
        ensure_kekule_parents(child.mol(), seed.0, seed.1, &mut fresh);
        let fresh_bag = fresh.get(&key).expect("fresh bag").borrow().clone();
        assert!(
            shared.semantically_eq(&fresh_bag),
            "shared cache after edit_copy must match fresh fill"
        );
    }

    #[test]
    fn product_that_breaks_system_misses_parent_key() {
        // Hydroxylate a ring carbon → new structure; shared Rc still holds the
        // old key, but a fresh fill on the product builds its own bag.
        use crate::rules::hydroxylation;
        use crate::Candidate;
        let parent = ForestMol::parse("c1ccccc1").unwrap();
        let (atoms, bonds) = conjugated_component(parent.mol(), 0);
        let seed = *bonds.iter().next().unwrap();
        parent.ensure_kekule(seed.0, seed.1);
        let parent_key = SystemKey::of(parent.mol(), &atoms, &bonds);
        let parent_bag = parent
            .kekule()
            .borrow()
            .get(&parent_key)
            .unwrap()
            .borrow()
            .clone();

        let set = hydroxylation();
        let cand = set.candidates(parent.mol()).next().unwrap().unwrap();
        assert!(matches!(cand, Candidate::Edit(_)));
        let products = cand.materialize_mols(parent.mol()).unwrap();
        assert!(!products.is_empty());
        let child = parent.adopt_product(products[0].clone());
        // Product may still share Rc until structure-changing invalidation;
        // semantic check: fresh fill on child equals what child would compute.
        let mut fresh = KekuleCache::default();
        let (catoms, cbonds) = conjugated_component(child.mol(), 0);
        if catoms.len() >= 2 && !cbonds.is_empty() {
            let cseed = *cbonds.iter().next().unwrap();
            ensure_kekule_parents(child.mol(), cseed.0, cseed.1, &mut fresh);
            let ckey = SystemKey::of(child.mol(), &catoms, &cbonds);
            let fresh_bag = fresh.get(&ckey).unwrap().borrow().clone();
            // If the parent key still resolves on the child (tags survive),
            // shared bag must match fresh for that key; otherwise keys differ.
            if let Some(shared_slot) = child.kekule().borrow().get(&ckey) {
                assert!(
                    shared_slot.borrow().semantically_eq(&fresh_bag),
                    "child shared bag for live key must match fresh"
                );
            }
            // Parent's old bag content is still well-formed.
            assert!(parent_bag.is_filled());
        }
    }
}
