//! Kekulé parents, one conjugated system at a time, on demand, shared.
//!
//! Whole-molecule enumeration is the product of independent rings (`2^n`
//! phenyls). Live forest does not do that. This cache:
//!
//! - assigns **one** conjugated component (other systems stay aromatic)
//! - fills **on demand** when a match names a bond in that system
//! - stores **assignment maps**, not baked mols, keyed by the system's
//!   atom set plus a fingerprint of its aromatic/bond shape
//! - lives on [`crate::forest_mol::ForestMol`] as one `Rc` for a copy tree: the first
//!   relative to fill a key shares it with every relative whose system
//!   still matches. An edit that changes kekulization of a system misses
//!   that key and starts a new bag.

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::rc::Rc;

use chematic::core::BondOrder;
use chematic::perception::find_sssr;
use chematic::smarts::{BondPrimitive, BondQuery, parse_smarts};

use crate::mol::{ForestError, Molecule, atom_idx, atom_usize};

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
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
pub fn conjugated_component(
    mol: &Molecule,
    start: usize,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
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
    (atoms, bonds)
}

fn order_code(order: BondOrder) -> Option<u8> {
    match order {
        BondOrder::Single => Some(1),
        BondOrder::Double => Some(2),
        BondOrder::Aromatic => Some(3),
        _ => Some(0),
    }
}

/// Identity of one conjugated system: indexes plus kekulization-relevant shape.
///
/// An edit that changes aromatic flags or bond orders inside the system is a
/// new key, so that bag is not reused. An edit elsewhere leaves this key
/// identical, so relatives share the bag.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SystemKey {
    pub atoms: BTreeSet<usize>,
    shape: Vec<((usize, usize), u8, bool, bool)>,
}

impl SystemKey {
    pub fn of(mol: &Molecule, atoms: &BTreeSet<usize>, bonds: &BTreeSet<(usize, usize)>) -> Self {
        let mut shape: Vec<_> = bonds
            .iter()
            .map(|&(a, b)| {
                let order = mol
                    .bond_between(atom_idx(a), atom_idx(b))
                    .map(|(_, bond)| order_code(bond.order).unwrap_or(0))
                    .unwrap_or(0);
                let a_aro = mol.atom(atom_idx(a)).aromatic;
                let b_aro = mol.atom(atom_idx(b)).aromatic;
                ((a, b), order, a_aro, b_aro)
            })
            .collect();
        shape.sort_unstable();
        Self {
            atoms: atoms.clone(),
            shape,
        }
    }
}

/// One system's kekulé assignments (bond-order maps, not baked mols).
#[derive(Clone, Debug, Default)]
pub struct SystemKekule {
    pub assignments: Vec<BTreeMap<(usize, usize), BondOrder>>,
    pub by_order: BTreeMap<((usize, usize), u8), usize>,
}

impl SystemKekule {
    pub fn is_filled(&self) -> bool {
        !self.assignments.is_empty()
    }
}

/// Shared forest-level map. One `Rc` per copy tree.
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
/// Same rules as Python `move_charge_with_bonds`: neutral C shifts H; neutral
/// aromatic atoms do not mint charge from 1.5→1/2 alone; already-charged atoms
/// and neutral non-aromatic heteroatoms follow the bond (N may gain charge when
/// bond order rises — iminium / amide resonance / pyridinium bookkeeping).
///
/// Closed-shell exception for oxygen: bond-order rise on neutral OH consumes H
/// (phenol/enol → carbonyl) instead of minting `[OH+]`. Radicals and protonated
/// carbonyls are not products — emit-path correctness, not sanitize rescue (C10).
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
        let delta = (new - old).round() as i8;
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
        } else if neutral && z == 8 && delta > 0 {
            // Phenol/enol → carbonyl: closed-shell O (H=0, charge=0). Do not
            // mint [OH+] — after the double is written, valence-inferred H is
            // already 0 so a "consume H else charge" branch wrongly charges.
            *mol = with_atom_explicit_h(mol, atom_idx(i), 0);
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
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    complete_assignment(
        mol,
        atoms,
        bonds,
        &BTreeSet::from([bond_key(seed.0, seed.1)]),
        &BTreeSet::new(),
    )
}

/// Complete a Kekulé assignment under π constraints.
///
/// `forced_doubles` seed the matching (exocyclic carbonyl/imine leaves, edited
/// single→double edges). `saturate` atoms leave must-match (Hydrogenation
/// `keep` + adds H — PatternInfo when the edit does not change the graph).
/// Returns `None` when the residual matching is impossible.
pub fn complete_assignment(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    forced_doubles: &BTreeSet<(usize, usize)>,
    saturate: &BTreeSet<usize>,
) -> Option<BTreeMap<(usize, usize), BondOrder>> {
    use chematic::core::kekulization::atom_must_be_matched;
    let must: BTreeSet<usize> = atoms
        .iter()
        .copied()
        .filter(|&i| !saturate.contains(&i) && atom_must_be_matched(mol, atom_idx(i)))
        .collect();
    let mut doubles: HashMap<usize, usize> = HashMap::new();
    for &(left, right) in forced_doubles {
        if doubles.contains_key(&left) || doubles.contains_key(&right) {
            if doubles.get(&left) != Some(&right) {
                return None;
            }
            continue;
        }
        doubles.insert(left, right);
        doubles.insert(right, left);
    }
    // Also treat existing Double bonds touching the system as forced (perception
    // after end edits — carbonyl O may sit outside `must`).
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
    forced_doubles: &BTreeSet<(usize, usize)>,
    saturate: &BTreeSet<usize>,
) -> usize {
    all_assignments(mol, atoms, bonds, forced_doubles, saturate).len()
}

/// Every distinct complete assignment under the constraints.
pub fn all_assignments(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    forced_doubles: &BTreeSet<(usize, usize)>,
    saturate: &BTreeSet<usize>,
) -> Vec<BTreeMap<(usize, usize), BondOrder>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    let mut consider = |forced: &BTreeSet<(usize, usize)>| {
        let Some(assignment) = complete_assignment(mol, atoms, bonds, forced, saturate) else {
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
    consider(forced_doubles);
    for &seed in bonds {
        if forced_doubles.contains(&seed) {
            continue;
        }
        let mut forced = forced_doubles.clone();
        forced.insert(seed);
        consider(&forced);
    }
    out
}

/// 2-core of a π graph: iteratively drop atoms with degree &lt; 2.
pub fn two_core(
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
    let mut atoms = atoms.clone();
    let mut bonds: BTreeSet<(usize, usize)> = bonds
        .iter()
        .copied()
        .filter(|&(a, b)| atoms.contains(&a) && atoms.contains(&b))
        .collect();
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
    (atoms, bonds)
}

/// Atoms that remain aromatic: cyclic 2-core components with a valid Kekulé
/// state and 4n+2 π electrons (all-carbon benzenoid shortcut).
///
/// Build `residual_atoms` / `residual_bonds` **after** dearomatizing edits
/// (epoxidation: drop the two sp³ atoms; quinone: drop the carbonyl carbon
/// from the aromatic candidate, or equivalently consume its demand with
/// forced C=O so it is absent here). Then ask which cyclic 2-core pieces
/// survive — not whether the original fused system is still aromatic.
pub fn aromatic_2core_atoms(
    mol: &Molecule,
    residual_atoms: &BTreeSet<usize>,
    residual_bonds: &BTreeSet<(usize, usize)>,
) -> BTreeSet<usize> {
    let (core_atoms, core_bonds) = two_core(residual_atoms, residual_bonds);
    let mut aromatic = BTreeSet::new();
    let mut seen = BTreeSet::new();
    for &start in &core_atoms {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut comp_atoms = BTreeSet::from([start]);
        while let Some(a) = stack.pop() {
            for &(left, right) in &core_bonds {
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
        let comp_bonds: BTreeSet<_> = core_bonds
            .iter()
            .copied()
            .filter(|&(a, b)| comp_atoms.contains(&a) && comp_atoms.contains(&b))
            .collect();
        // Cyclic: connected component with at least as many edges as atoms.
        if comp_bonds.len() < comp_atoms.len() {
            continue;
        }
        if complete_assignment(
            mol,
            &comp_atoms,
            &comp_bonds,
            &BTreeSet::new(),
            &BTreeSet::new(),
        )
        .is_none()
        {
            continue;
        }
        let n = comp_atoms.len();
        if n >= 6 && n % 4 == 2 {
            aromatic.extend(comp_atoms);
        }
    }
    aromatic
}

/// Stamp one system's assignment onto `mol`. Other systems stay as they were.
///
/// Applies [`move_charge_with_bonds`] so charge separation on valid writings
/// (including N that gains charge) matches Python reactant overlays.
pub fn overlay(
    mol: &Molecule,
    atoms: &BTreeSet<usize>,
    assignment: &BTreeMap<(usize, usize), BondOrder>,
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
    move_charge_with_bonds(&mut out, &before, &aromatic);
    out
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
    let mut bag = slot.borrow_mut();
    let mut seen = BTreeSet::new();
    for &seed in bonds {
        let Some(bond_orders) = match_assignment(mol, atoms, bonds, seed) else {
            continue;
        };
        let signature: Vec<_> = bond_orders
            .iter()
            .map(|(&k, &o)| (k, order_code(o)))
            .collect();
        if !seen.insert(signature) {
            continue;
        }
        let index = bag.assignments.len();
        for (key, order) in &bond_orders {
            if let Some(code) = order_code(*order) {
                if code == 1 || code == 2 {
                    bag.by_order.entry((*key, code)).or_insert(index);
                }
            }
        }
        bag.assignments.push(bond_orders);
    }
}

fn system_of(mol: &Molecule, left: usize, right: usize) -> (SystemKey, BTreeSet<(usize, usize)>) {
    let (mut atoms, mut bonds) = conjugated_component(mol, left);
    let seed = bond_key(left, right);
    if !bonds.contains(&seed) && mol.bond_between(atom_idx(left), atom_idx(right)).is_some() {
        bonds.insert(seed);
        atoms.insert(right);
    }
    (SystemKey::of(mol, &atoms, &bonds), bonds)
}

/// Fill (or reuse) the bag for the system that contains `(left, right)`.
pub fn ensure_kekule_parents(
    mol: &Molecule,
    left: usize,
    right: usize,
    cache: &mut KekuleCache,
) -> Rc<RefCell<SystemKekule>> {
    let (key, bonds) = system_of(mol, left, right);
    let atoms = key.atoms.clone();
    let slot = cache.slot(key);
    fill_slot(mol, &atoms, &bonds, &slot);
    slot
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
    let index = *bag.by_order.get(&(bond_key(left, right), order))?;
    let assignment = bag.assignments.get(index)?;
    Some(overlay(mol, &key.atoms, assignment))
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
        let key = SystemKey::of(mol, atoms, &bonds);
        let slot = cache.slot(key);
        fill_slot(mol, atoms, &bonds, &slot);
        let parents = slot
            .borrow()
            .assignments
            .iter()
            .map(|assignment| overlay(mol, atoms, assignment))
            .collect();
        return EndParents {
            parents,
            same_system: true,
        };
    }
    let (start_atoms, start_bonds) = conjugated_component(mol, start);
    let (end_atoms, end_bonds) = conjugated_component(mol, end);
    let start_key = SystemKey::of(mol, &start_atoms, &start_bonds);
    let end_key = SystemKey::of(mol, &end_atoms, &end_bonds);
    let start_slot = cache.slot(start_key);
    fill_slot(mol, &start_atoms, &start_bonds, &start_slot);
    let same = start_atoms == end_atoms;
    if !same {
        let end_slot = cache.slot(end_key);
        fill_slot(mol, &end_atoms, &end_bonds, &end_slot);
    }
    let mut parents: Vec<Molecule> = start_slot
        .borrow()
        .assignments
        .iter()
        .map(|assignment| overlay(mol, &start_atoms, assignment))
        .collect();
    if !same {
        let end_slot = cache
            .get(&SystemKey::of(mol, &end_atoms, &end_bonds))
            .expect("end system filled");
        parents.extend(
            end_slot
                .borrow()
                .assignments
                .iter()
                .map(|assignment| overlay(mol, &end_atoms, assignment)),
        );
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
pub fn kekule_forms(mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
    if !mol.atoms().any(|(_, atom)| atom.aromatic) {
        return Ok(vec![mol.clone()]);
    }
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
        let (atoms, bonds) = conjugated_component(mol, start);
        covered.extend(&atoms);
        let aromatic: BTreeSet<usize> = atoms
            .iter()
            .copied()
            .filter(|&i| mol.atom(atom_idx(i)).aromatic)
            .collect();
        if aromatic.len() < 2 {
            continue;
        }
        let arbonds: BTreeSet<_> = bonds
            .iter()
            .copied()
            .filter(|&(a, b)| aromatic.contains(&a) && aromatic.contains(&b))
            .collect();
        if arbonds.is_empty() {
            continue;
        }
        let key = SystemKey::of(mol, &aromatic, &arbonds);
        let slot = cache.slot(key);
        fill_slot(mol, &aromatic, &arbonds, &slot);
        forms.extend(
            slot.borrow()
                .assignments
                .iter()
                .map(|assignment| overlay(mol, &aromatic, assignment)),
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
    let want = smirks_mapped_bond_order(smirks).unwrap_or(2.0);
    let want_code: u8 = if want >= 1.5 { 2 } else { 1 };
    ensure_kekule_parents(mol, left, right, cache);
    Ok(parent_for_bond(mol, cache, left, right, want_code).unwrap_or_else(|| mol.clone()))
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
        let mol = parse_mol(smiles).unwrap();
        let mut cache = KekuleCache::default();
        for mapped in smarts_matches(&mol, "[#6:1]=,:[#6:2]").unwrap() {
            ensure_kekule_parents(&mol, mapped[&1], mapped[&2], &mut cache);
        }
        (mol, cache)
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
        let aromatic = aromatic_2core_atoms(&mol, &residual_atoms, &residual_bonds);
        assert_eq!(aromatic.len(), 6, "one benzenoid sextet survives: {aromatic:?}");
        assert!(!aromatic.contains(&outer.0) && !aromatic.contains(&outer.1));
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
        let mol = parse_mol("c1ccc(-c2ccccc2)cc1").unwrap();
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
        let mol = parse_mol("c1ccsc1").unwrap();
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
        let mol = parse_mol("c1ccccc1").unwrap();
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
}
