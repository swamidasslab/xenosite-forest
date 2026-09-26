//! ResonancePair edits as π-matching constraints (no path flip).
//!
//! Match each [`Edit::PairEndpoint`] SMARTS, apply end edits, perceive forced
//! doubles / saturate sites, complete the residual Kekulé matching, then mark
//! surviving aromatic atoms via the cyclic 2-core test (HEURISTICS). Methide
//! is an effect field — two methide ends are allowed.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use chematic::core::{Atom, BondOrder, Element};
use chematic::perception::find_sssr;

use crate::kekule::{
    all_assignments, aromatic_2core_atoms, bond_order_sums, conjugated_component,
    move_charge_with_bonds, with_atom_explicit_h,
};
use crate::mol::{ForestError, Molecule, aromatize, atom_idx, atom_usize, canon_smiles};
use crate::pattern::{Edit, PatternInfo};
use crate::smarts::smarts_matches;
use crate::valence::accept_product;

fn bond_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn clear_aromatic(mol: &mut Molecule) {
    for atom in mol.atoms().map(|(idx, _)| idx).collect::<Vec<_>>() {
        if mol.atom(atom).aromatic {
            *mol = mol.with_atom_aromatic(atom, false);
        }
    }
}

/// Demote aromatic bonds touching `system` to single so matching can assign.
fn demote_aromatic_bonds(mol: &mut Molecule, system: &HashSet<usize>) {
    let idxs: Vec<_> = mol
        .bonds()
        .filter_map(|(idx, bond)| {
            if bond.order != BondOrder::Aromatic {
                return None;
            }
            let a = atom_usize(bond.atom1);
            let b = atom_usize(bond.atom2);
            if system.contains(&a) || system.contains(&b) {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    for idx in idxs {
        mol.set_bond_order(idx, BondOrder::Single);
    }
}

/// `keep` + adds H → saturate site atoms (leave must-match, +H later).
fn saturate_sites(pattern: &PatternInfo, mapped: &BTreeMap<u16, usize>) -> BTreeSet<usize> {
    let Edit::PairEndpoint(edit) = &pattern.edit else {
        return BTreeSet::new();
    };
    if edit.as_str() != "keep" {
        return BTreeSet::new();
    }
    let adds = pattern.effect.adds.as_deref().unwrap_or("");
    if !adds.contains('H') {
        return BTreeSet::new();
    }
    mapped.get(&1).copied().into_iter().collect()
}

/// Residual π graph: conjugated edges among system atoms minus saturate sites.
///
/// Topology from the parent (pre-edit) conjugated component. Saturate atoms
/// leave must-match; demand-consumed atoms for the aromatic 2-core gate are
/// stripped after forced doubles are perceived.
fn residual_pi_graph(
    parent: &Molecule,
    system: &HashSet<usize>,
    saturate: &BTreeSet<usize>,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
    let atoms: BTreeSet<usize> = system
        .iter()
        .copied()
        .filter(|a| !saturate.contains(a))
        .collect();
    let mut bonds = BTreeSet::new();
    if let Some(&seed) = atoms.iter().next().or_else(|| system.iter().next()) {
        let (_, comp_bonds) = conjugated_component(parent, seed);
        for (a, b) in comp_bonds {
            if atoms.contains(&a) && atoms.contains(&b) {
                bonds.insert(bond_key(a, b));
            }
        }
    }
    // Fallback: any bond between residual atoms (full-mol system path).
    if bonds.is_empty() {
        for (_, bond) in parent.bonds() {
            let a = atom_usize(bond.atom1);
            let b = atom_usize(bond.atom2);
            if atoms.contains(&a) && atoms.contains(&b) {
                match bond.order {
                    BondOrder::Single
                    | BondOrder::Double
                    | BondOrder::Triple
                    | BondOrder::Aromatic => {
                        bonds.insert(bond_key(a, b));
                    }
                    _ => {}
                }
            }
        }
    }
    (atoms, bonds)
}

/// Forced doubles from the post-edit graph (perception; PatternInfo fills gaps
/// via saturate_sites). Any double touching the residual atoms is a constraint.
fn perceive_forced_doubles(
    edited: &Molecule,
    atoms: &BTreeSet<usize>,
) -> BTreeSet<(usize, usize)> {
    let mut forced = BTreeSet::new();
    for (_, bond) in edited.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let a = atom_usize(bond.atom1);
        let b = atom_usize(bond.atom2);
        if atoms.contains(&a) || atoms.contains(&b) {
            forced.insert(bond_key(a, b));
        }
    }
    forced
}

/// Drop atoms whose π demand is consumed by a forced double (exocyclic leaf or
/// in-system forced edge) before the cyclic 2-core aromaticity test.
fn aromatic_residual_after_forced(
    atoms: &BTreeSet<usize>,
    bonds: &BTreeSet<(usize, usize)>,
    forced: &BTreeSet<(usize, usize)>,
) -> (BTreeSet<usize>, BTreeSet<(usize, usize)>) {
    let mut arom_atoms = atoms.clone();
    for &(a, b) in forced {
        let a_in = arom_atoms.contains(&a);
        let b_in = arom_atoms.contains(&b);
        if a_in && b_in {
            // Forced edge inside the residual (phenol C=O, iminium, …): both
            // demands consumed — drop both from the aromatic candidate.
            arom_atoms.remove(&a);
            arom_atoms.remove(&b);
        } else if a_in {
            arom_atoms.remove(&a);
        } else if b_in {
            arom_atoms.remove(&b);
        }
    }
    let arom_bonds: BTreeSet<_> = bonds
        .iter()
        .copied()
        .filter(|&(a, b)| arom_atoms.contains(&a) && arom_atoms.contains(&b))
        .collect();
    (arom_atoms, arom_bonds)
}

fn bump_h(mol: &mut Molecule, atom: usize, delta: i8) {
    let idx = atom_idx(atom);
    // implicit_hydrogen_count already returns stored hydrogen_count when set.
    let total = mol.implicit_hydrogen_count(idx) as i16 + delta as i16;
    if total >= 0 {
        *mol = with_atom_explicit_h(mol, idx, total as u8);
    }
}

/// Capability OR, resolved against whether the conjugated system is aromatic.
fn merge_dearomatizes(left: &PatternInfo, right: &PatternInfo, system_aromatic: bool) -> bool {
    (left.effect.dearomatizes || right.effect.dearomatizes) && system_aromatic
}

fn ring_sets(mol: &Molecule) -> HashMap<usize, BTreeSet<usize>> {
    let mut out: HashMap<usize, BTreeSet<usize>> = HashMap::new();
    for (ring_i, ring) in find_sssr(mol).rings().iter().enumerate() {
        for &atom in ring {
            out.entry(atom_usize(atom)).or_default().insert(ring_i);
        }
    }
    out
}

fn same_rings(rings: &HashMap<usize, BTreeSet<usize>>, a: usize, b: usize) -> bool {
    rings.get(&a) == rings.get(&b)
}

fn site_atom(mapped: &BTreeMap<u16, usize>, pattern: &PatternInfo) -> Option<usize> {
    mapped.get(&pattern.primary_map()).copied()
}

/// Non-site mapped atoms when [`Effect::exclusive_partner`] is set.
fn exclusive_partner_atoms(
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
) -> BTreeSet<usize> {
    if !pattern.effect.exclusive_partner {
        return BTreeSet::new();
    }
    let Some(site) = site_atom(mapped, pattern) else {
        return BTreeSet::new();
    };
    mapped
        .values()
        .copied()
        .filter(|&idx| idx != site)
        .collect()
}

/// True when an exclusive partner atom appears on the other end's map.
fn shared_exclusive_partner(
    map1: &BTreeMap<u16, usize>,
    info1: &PatternInfo,
    map2: &BTreeMap<u16, usize>,
    info2: &PatternInfo,
) -> bool {
    let exclusive1 = exclusive_partner_atoms(map1, info1);
    let exclusive2 = exclusive_partner_atoms(map2, info2);
    if exclusive1.is_empty() && exclusive2.is_empty() {
        return false;
    }
    let other1: BTreeSet<usize> = map2.values().copied().collect();
    let other2: BTreeSet<usize> = map1.values().copied().collect();
    !exclusive1.is_disjoint(&other1) || !exclusive2.is_disjoint(&other2)
}

fn edit_end(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let Edit::PairEndpoint(edit) = &pattern.edit else {
        return false;
    };
    match edit.as_str() {
        "single_to_double" => edit_single_to_double(mol, mapped, pattern, rings),
        "add_carbonyl_o" => edit_add_carbonyl_o(mol, mapped),
        "replace_halogen" => edit_replace_halogen(mol, mapped),
        "iminium" => {
            if !edit_single_to_double(mol, mapped, pattern, rings) {
                return false;
            }
            let Some(&n) = mapped.get(&2) else {
                return false;
            };
            *mol = mol.with_atom_charge(atom_idx(n), 1);
            true
        }
        "dealkylate" => edit_dealkylate(mol, mapped, pattern, rings),
        "keep" => mapped.contains_key(&1),
        _ => false,
    }
}

/// Atoms whose σ skeleton changed in an end edit (cleavage / atom add). Their
/// bond-order baseline for [`move_charge_with_bonds`] is the post-edit graph so
/// H is not double-counted (methyl leave must be CH4, not [CH5]).
fn skeleton_changed_atoms(
    left: &PatternInfo,
    map1: &BTreeMap<u16, usize>,
    right: &PatternInfo,
    map2: &BTreeMap<u16, usize>,
) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    for (pattern, mapped) in [(left, map1), (right, map2)] {
        let Edit::PairEndpoint(edit) = &pattern.edit else {
            continue;
        };
        match edit.as_str() {
            "dealkylate" => {
                if let (Some(&hetero), Some(&alkyl)) = (mapped.get(&2), mapped.get(&3)) {
                    out.insert(hetero);
                    out.insert(alkyl);
                }
            }
            "add_carbonyl_o" => {
                // New oxygen is absent from parent `before`; ring carbon’s new
                // C=O is a π constraint — keep parent baseline for that carbon.
            }
            _ => {}
        }
    }
    out
}

fn edit_single_to_double(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let (Some(&a), Some(&b)) = (mapped.get(&1), mapped.get(&2)) else {
        return false;
    };
    if pattern.skip_same_rings && same_rings(rings, a, b) {
        return false;
    }
    let Some((bond_idx, _)) = mol.bond_between(atom_idx(a), atom_idx(b)) else {
        return false;
    };
    mol.set_bond_order(bond_idx, BondOrder::Double);
    true
}

fn edit_add_carbonyl_o(mol: &mut Molecule, mapped: &BTreeMap<u16, usize>) -> bool {
    let Some(&carbon) = mapped.get(&1) else {
        return false;
    };
    let (mut next, oxygen) = mol.with_atom_added(Atom::organic(Element::O));
    if next
        .add_bond(atom_idx(carbon), oxygen, BondOrder::Double)
        .is_err()
    {
        return false;
    }
    *mol = next;
    true
}

fn edit_replace_halogen(mol: &mut Molecule, mapped: &BTreeMap<u16, usize>) -> bool {
    let (Some(&carbon), Some(&halogen)) = (mapped.get(&1), mapped.get(&2)) else {
        return false;
    };
    *mol = mol.with_atom_element(atom_idx(halogen), Element::O);
    *mol = mol.with_atom_charge(atom_idx(halogen), 0);
    let Some((bond_idx, _)) = mol.bond_between(atom_idx(carbon), atom_idx(halogen)) else {
        return false;
    };
    mol.set_bond_order(bond_idx, BondOrder::Double);
    true
}

fn edit_dealkylate(
    mol: &mut Molecule,
    mapped: &BTreeMap<u16, usize>,
    pattern: &PatternInfo,
    rings: &HashMap<usize, BTreeSet<usize>>,
) -> bool {
    let (Some(&hetero), Some(&alkyl)) = (mapped.get(&2), mapped.get(&3)) else {
        return false;
    };
    if !edit_single_to_double(mol, mapped, pattern, rings) {
        return false;
    }
    if let Some((bi, _)) = mol.bond_between(atom_idx(hetero), atom_idx(alkyl)) {
        *mol = mol.with_bond_removed(bi);
    } else {
        return false;
    }
    // Closed-shell H after cleavage (no radicals, no [CH5]/[CH6]). Fill from
    // remaining bond orders + charge; do not stack bump on valence-inferred H.
    fill_closed_shell_h(mol, alkyl);
    fill_closed_shell_h(mol, hetero);
    true
}

/// Set explicit H so atom valence is complete (closed shell). Used after
/// cleavage so leave fragments are CH4 / HN= / O= without sanitize rescue.
fn fill_closed_shell_h(mol: &mut Molecule, atom: usize) {
    let idx = atom_idx(atom);
    let a = mol.atom(idx);
    let z = a.element.atomic_number();
    let charge = a.charge as i16;
    let mut bond_sum = 0.0_f32;
    for (_nbr, bidx) in mol.neighbors(idx) {
        bond_sum += match mol.bond(bidx).order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::Quadruple => 4.0,
            _ => 1.0,
        };
    }
    // Organic valence targets (neutral): C 4, N 3, O 2. Charge adjusts.
    let target = match z {
        6 => 4 + charge,
        7 => 3 + charge,
        8 => 2 + charge,
        _ => return,
    };
    let need = target - bond_sum.round() as i16;
    if need >= 0 {
        *mol = with_atom_explicit_h(mol, idx, need as u8);
    }
}

/// One pair emission before RuleSet packaging.
#[derive(Clone, Debug)]
pub struct PairEmission {
    pub site: usize,
    pub pattern_name: String,
    pub products: Vec<String>,
}

/// Discovered pair site before constraint materialize / product CSMI.
///
/// Carries merged [`crate::pattern::Effect`] so a search can filter without
/// running the edit. [`PairCandidate::materialize`] applies end edits as π
/// constraints and completes the residual matching (no path flip).
///
/// [`Self::rule_path`] is the same leaf-first namespace as [`crate::candidate::Candidate`].
/// Prefer [`crate::ruleset::RuleSet::candidates`] /
/// [`crate::ruleset::RuleSet::metabolites`] so the leaf name is stamped
/// automatically. Bare [`pair_candidates`] is `pub(crate)` and leaves
/// `rule_path` empty — RuleSet doors stamp via `stamp_pair_paths`.
#[derive(Clone, Debug)]
pub struct PairCandidate {
    pub site: usize,
    pub pattern_name: String,
    pub left: PatternInfo,
    pub right: PatternInfo,
    /// Merged end effects (dearomatizes resolved against system aromaticity).
    pub effect: crate::pattern::Effect,
    /// Leaf-first rule namespace (emitting set, then containers). Empty when
    /// discovered via bare [`pair_candidates`] without a
    /// [`crate::ruleset::RuleSet`] stamp.
    pub rule_path: Vec<Option<String>>,
    map1: BTreeMap<u16, usize>,
    map2: BTreeMap<u16, usize>,
    start: usize,
    end: usize,
    system: HashSet<usize>,
}

impl PairCandidate {
    /// Emitting (leaf) rule name when discovered under a named set.
    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }

    /// Named segments of [`Self::rule_path`] (unnamed sets omitted).
    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    /// Hop / emission rule label: leaf [`RuleSet`] name, else [`PatternInfo::name`]
    /// on the first end (same fallback as [`crate::candidate::Candidate::rule_name`]).
    pub fn rule_name(&self) -> &str {
        self.leaf_rule().unwrap_or(self.left.name.as_str())
    }

    /// Stamp leaf-first namespace (replaces any prior path).
    ///
    /// Prefer [`crate::ruleset::RuleSet::stamp_pair_paths`] /
    /// [`crate::ruleset::RuleSet::with_outer_path`] at discovery time so call
    /// sites do not mint leaf names by hand.
    pub fn with_rule_path(mut self, rule_path: Vec<Option<String>>) -> Self {
        self.rule_path = rule_path;
        self
    }

    /// Discovery site atoms for each end (Python `end_atoms`).
    pub fn end_atoms(&self) -> Option<(usize, usize)> {
        let a = site_atom(&self.map1, &self.left)?;
        let b = site_atom(&self.map2, &self.right)?;
        Some((a, b))
    }

    /// Conjugated-system anchors (legacy path_ends; materialize uses constraints).
    pub fn path_ends(&self) -> (usize, usize) {
        (self.start, self.end)
    }

    pub fn materialize_mols(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        let rings = ring_sets(mol);
        let mut rw = mol.clone();
        clear_aromatic(&mut rw);
        demote_aromatic_bonds(&mut rw, &self.system);
        if !edit_end(&mut rw, &self.map1, &self.left, &rings) {
            return Ok(Vec::new());
        }
        if !edit_end(&mut rw, &self.map2, &self.right, &rings) {
            return Ok(Vec::new());
        }

        let saturate = saturate_sites(&self.left, &self.map1)
            .union(&saturate_sites(&self.right, &self.map2))
            .copied()
            .collect::<BTreeSet<_>>();
        let (atoms, bonds) = residual_pi_graph(mol, &self.system, &saturate);
        let forced = perceive_forced_doubles(&rw, &atoms);
        let assignments = all_assignments(&rw, &atoms, &bonds, &forced, &saturate);
        if assignments.is_empty() {
            return Ok(Vec::new());
        }

        // Residual aromaticity: cyclic 2-core after dropping demand-consumed
        // atoms (forced doubles / saturate) — HEURISTICS, not sanitize flags.
        let (arom_atoms, arom_bonds) = aromatic_residual_after_forced(&atoms, &bonds, &forced);
        let aromatic_core = aromatic_2core_atoms(&rw, &arom_atoms, &arom_bonds);
        if self.effect.dearomatizes {
            let edited: HashSet<usize> = self
                .system
                .iter()
                .copied()
                .filter(|i| mol.atom(atom_idx(*i)).aromatic)
                .collect();
            if !edited.is_empty() && edited.iter().all(|i| aromatic_core.contains(i)) {
                // Capability claimed dearomatization; whole edited system still aromatic.
                return Ok(Vec::new());
            }
        }

        let mut products = Vec::new();
        let mut local_csmi = BTreeSet::new();
        // Charge/H follow parent π bond sums (aromatic = 1.5). Cleaved atoms
        // use the post-edit baseline so leave fragments stay closed-shell
        // (CH4 not [CH5]) — no sanitize rescue, no Rust regression on C16.
        let mut before = bond_order_sums(mol);
        let post_edit = bond_order_sums(&rw);
        for a in skeleton_changed_atoms(&self.left, &self.map1, &self.right, &self.map2) {
            if let Some(&v) = post_edit.get(&a) {
                before.insert(a, v);
            }
        }
        let was_aromatic: HashSet<usize> = mol
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(atom_usize(idx)))
            .collect();
        for assignment in assignments {
            let mut product = rw.clone();
            for (&(left, right), &order) in &assignment {
                if let Some((bond_idx, _)) =
                    product.bond_between(atom_idx(left), atom_idx(right))
                {
                    product.set_bond_order(bond_idx, order);
                }
            }
            for &atom in &atoms {
                if product.atom(atom_idx(atom)).aromatic {
                    product = product.with_atom_aromatic(atom_idx(atom), false);
                }
            }
            for &atom in &saturate {
                bump_h(&mut product, atom, 1);
            }
            move_charge_with_bonds(&mut product, &before, &was_aromatic);
            // Stamp surviving aromatic 2-core; leave the rest localized.
            for &atom in &aromatic_core {
                product = product.with_atom_aromatic(atom_idx(atom), true);
            }
            // Perception finish (RDKit-parity aromaticity). Dearomatize refuse
            // already used the 2-core gate above — not sanitize flags.
            let checked = aromatize(&product);
            for frag in checked.fragments() {
                if !accept_product(&frag) {
                    continue;
                }
                let smiles = canon_smiles(&frag);
                if local_csmi.insert(smiles) {
                    products.push(frag);
                }
            }
        }
        Ok(products)
    }

    pub fn materialize(&self, mol: &Molecule) -> Result<Vec<String>, ForestError> {
        Ok(self
            .materialize_mols(mol)?
            .iter()
            .map(canon_smiles)
            .collect())
    }

    /// Elementary plan site atoms (pair ends), falling back to discovery site.
    pub fn plan_site_atoms(&self) -> Vec<usize> {
        match self.end_atoms() {
            Some((a, b)) => vec![a, b],
            None => vec![self.site],
        }
    }

    pub fn emit(&self, mol: &Molecule) -> Result<Option<PairEmission>, ForestError> {
        let mols = self.materialize_mols(mol)?;
        if mols.is_empty() {
            return Ok(None);
        }
        crate::formula_check::check_effect_delta_formula(
            mol,
            &self.effect,
            &mols,
            &self.pattern_name,
        );
        let products = mols
            .iter()
            .map(crate::mol::canon_smiles)
            .collect::<Vec<_>>();
        Ok(Some(PairEmission {
            site: self.site,
            pattern_name: self.pattern_name.clone(),
            products,
        }))
    }
}

type EndpointHit<'a> = (BTreeMap<u16, usize>, &'a PatternInfo);

fn merge_effect_fields(
    left: &PatternInfo,
    right: &PatternInfo,
    system_aromatic: bool,
) -> crate::pattern::Effect {
    let adds = match (&left.effect.adds, &right.effect.adds) {
        (None, None) => None,
        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
        (Some(a), Some(b)) => Some(format!("{a}{b}")),
    };
    let removes = match (&left.effect.removes, &right.effect.removes) {
        (None, None) => None,
        (Some(a), None) | (None, Some(a)) => Some(a.clone()),
        (Some(a), Some(b)) => Some(format!("{a}{b}")),
    };
    crate::pattern::Effect {
        adds,
        removes,
        delta_formula: crate::pattern::merge_delta_formula(
            &left.effect.resolved_delta_formula(),
            &right.effect.resolved_delta_formula(),
        ),
        leave_formula: {
            let mut leave = left.effect.leave_formula.clone();
            for (el, n) in &right.effect.leave_formula {
                *leave.entry(el.clone()).or_insert(0) += n;
            }
            leave.retain(|_, n| *n != 0);
            leave
        },
        cleaves: left.effect.cleaves || right.effect.cleaves,
        leave_count: left.effect.leave_count.or(right.effect.leave_count),
        methide: left.effect.methide || right.effect.methide,
        exclusive_partner: left.effect.exclusive_partner || right.effect.exclusive_partner,
        dearomatizes: merge_dearomatizes(left, right, system_aromatic),
        partner: left
            .effect
            .partner
            .clone()
            .or_else(|| right.effect.partner.clone()),
        when: None,
    }
}

/// Discover pair sites without materializing products.
///
/// Crate-internal discovery primitive: returns pairs with empty
/// [`PairCandidate::rule_path`]. Prefer [`crate::ruleset::RuleSet::candidates`]
/// / [`crate::ruleset::RuleSet::metabolites`].
pub(crate) fn pair_candidates(
    mol: &Molecule,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairCandidate>, ForestError> {
    let endpoints: Vec<&PatternInfo> = endpoints
        .iter()
        .filter(|p| matches!(p.edit, Edit::PairEndpoint(_)))
        .collect();
    if endpoints.is_empty() {
        return Ok(Vec::new());
    }

    let mut hits: HashMap<usize, Vec<EndpointHit<'_>>> = HashMap::new();
    for pattern in &endpoints {
        for mapped in smarts_matches(mol, &pattern.smarts)? {
            let Some(&anchor) = mapped.get(&1) else {
                continue;
            };
            hits.entry(anchor).or_default().push((mapped, *pattern));
        }
    }
    if hits.len() < 2 {
        return Ok(Vec::new());
    }

    let mut systems: Vec<HashSet<usize>> = Vec::new();
    let mut covered = HashSet::new();
    for &anchor in hits.keys() {
        if !covered.insert(anchor) {
            continue;
        }
        let (atoms, _) = conjugated_component(mol, anchor);
        let set: HashSet<usize> = atoms.iter().copied().collect();
        covered.extend(&set);
        if hits.keys().filter(|a| set.contains(a)).count() >= 2 {
            systems.push(set);
        }
    }
    if systems.is_empty() {
        let all: HashSet<usize> = mol.atoms().map(|(i, _)| atom_usize(i)).collect();
        if hits.keys().filter(|a| all.contains(a)).count() >= 2 {
            systems.push(all);
        }
    }

    let mut out = Vec::new();
    let mut seen_sig: BTreeSet<(usize, usize, String, String, usize)> = BTreeSet::new();
    let gens = crate::orbits::atom_bond_generators(mol);
    let n_atoms = mol.atom_count();

    for system in &systems {
        let anchors: Vec<usize> = hits
            .keys()
            .copied()
            .filter(|a| system.contains(a))
            .collect();
        if anchors.len() < 2 {
            continue;
        }
        let system_aromatic = system.iter().any(|&i| mol.atom(atom_idx(i)).aromatic);
        let mut sorted_anchors = anchors;
        sorted_anchors.sort_unstable();
        for i in 0..sorted_anchors.len() {
            for j in (i + 1)..sorted_anchors.len() {
                let start = sorted_anchors[i];
                let end = sorted_anchors[j];
                let Some(hits_a) = hits.get(&start) else {
                    continue;
                };
                let Some(hits_b) = hits.get(&end) else {
                    continue;
                };
                for (map1, info1) in hits_a {
                    for (map2, info2) in hits_b {
                        let (Some(site_a), Some(site_b)) =
                            (site_atom(map1, info1), site_atom(map2, info2))
                        else {
                            continue;
                        };
                        if site_a == site_b {
                            continue;
                        }
                        if shared_exclusive_partner(map1, info1, map2, info2) {
                            continue;
                        }
                        let (n1, n2) = if info1.name <= info2.name {
                            (info1.name.clone(), info2.name.clone())
                        } else {
                            (info2.name.clone(), info1.name.clone())
                        };
                        let sa = site_a.min(site_b);
                        let sb = site_a.max(site_b);
                        let orbit =
                            crate::orbits::atom_pair_orbit_id_with_gens(&gens, n_atoms, sa, sb);
                        if !seen_sig.insert((sa, sb, n1.clone(), n2.clone(), orbit)) {
                            continue;
                        }
                        out.push(PairCandidate {
                            site: sa,
                            pattern_name: format!("{n1}+{n2}"),
                            left: (*info1).clone(),
                            right: (*info2).clone(),
                            effect: merge_effect_fields(info1, info2, system_aromatic),
                            rule_path: Vec::new(),
                            map1: map1.clone(),
                            map2: map2.clone(),
                            start,
                            end,
                            system: system.clone(),
                        });
                    }
                }
            }
        }
    }
    Ok(out)
}

/// Run ResonancePair metabolize for the given endpoint patterns.
pub fn pair_metabolize(
    mol: &Molecule,
    endpoints: &[PatternInfo],
) -> Result<Vec<PairEmission>, ForestError> {
    let mut emissions = Vec::new();
    let mut seen_csmi: BTreeSet<BTreeSet<String>> = BTreeSet::new();
    for candidate in pair_candidates(mol, endpoints)? {
        let Some(emission) = candidate.emit(mol)? else {
            continue;
        };
        let key: BTreeSet<String> = emission.products.iter().cloned().collect();
        if !seen_csmi.insert(key) {
            continue;
        }
        emissions.push(emission);
    }
    Ok(emissions)
}

/// Path dehydrogenation of para-hydroquinone (door smoke for pair metabolize).
pub fn dehydrogenate_hydroquinone(mol: &Molecule) -> Result<Vec<String>, ForestError> {
    let phenol_end = PatternInfo {
        name: "phenol_end".into(),
        smarts: "[#6:1]-[#8H:2]".into(),
        site_kind: crate::pattern::SiteKind::AtomPair,
        site_map: vec![2],
        edit: Edit::PairEndpoint("single_to_double".into()),
        effect: crate::pattern::Effect {
            removes: Some("H".into()),
            dearomatizes: true,
            ..Default::default()
        }
        .sealed(),
        possibilities: Vec::new(),
        skip_same_rings: false,
        cleave_side_group: None,
        search_bias: 0,
    };
    Ok(pair_metabolize(mol, &[phenol_end])?
        .into_iter()
        .flat_map(|e| e.products)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use crate::rules::{dehydrogenation, quinone_formation};

    #[test]
    fn hydroquinone_yields_benzoquinone() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let products = dehydrogenate_hydroquinone(&mol).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            products.iter().any(|s| canon_of(s).unwrap() == want),
            "want {want}, got {products:?}"
        );
    }

    #[test]
    fn dehydrogenation_rule_pair_door_on_hydroquinone() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        assert!(!endpoints.is_empty());
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| { e.products.iter().any(|p| canon_of(p).unwrap() == want) }),
            "{emissions:?}"
        );
    }

    #[test]
    fn quinone_formation_add_carbonyl_on_benzene() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let endpoints: Vec<_> = quinone_formation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let emissions = pair_metabolize(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions.iter().any(|e| {
                e.pattern_name.contains("add_carbonyl_o")
                    && e.products.iter().any(|p| canon_of(p).unwrap() == want)
            }),
            "expected para-quinone from two add_carbonyl_o; got {emissions:?}"
        );
    }

    #[test]
    fn pair_candidates_defer_materialize() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let endpoints: Vec<_> = dehydrogenation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect();
        let cands = pair_candidates(&mol, &endpoints).unwrap();
        assert!(!cands.is_empty());
        assert!(
            cands.iter().any(|c| c.effect.dearomatizes),
            "aromatic hydroquinone pair should resolve dearomatizes"
        );
        // Refuse before materialize.
        let kept: Vec<_> = cands
            .into_iter()
            .filter(|c| !c.effect.dearomatizes)
            .collect();
        assert!(kept.is_empty());
        // Accept and materialize.
        let cands = pair_candidates(&mol, &endpoints).unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(cands.iter().any(|c| {
            c.materialize(&mol)
                .unwrap()
                .iter()
                .any(|p| canon_of(p).unwrap() == want)
        }));
    }

    #[test]
    fn quinone_dealkylate_splits_fragments_like_find_path() {
        // Python split_fragments: ['C', 'O=C1C=CC(=O)C=C1'] as two products.
        // find_path bifurcation needs n_products >= 2 (not one C.quinone mol).
        use crate::forest_mol::ForestMol;
        use crate::product_graph::{ProductGraphConfig, product_layer};
        use crate::rules::quinone_formation;

        let set = quinone_formation();
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let layer = ProductGraphConfig {
            target: None,
            max_nodes: usize::MAX,
            max_depth: usize::MAX,
        };
        let children = product_layer(&parent, &set, &layer).unwrap();
        let want_q = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        let want_me = canon_of("C").unwrap();
        let csmi: Vec<String> = children
            .iter()
            .map(|c| c.child.csmi().as_ref().to_string())
            .collect();
        assert!(
            csmi.iter().any(|s| canon_of(s).unwrap() == want_q),
            "quinone fragment missing: {csmi:?}"
        );
        assert!(
            csmi.iter().any(|s| canon_of(s).unwrap() == want_me),
            "methyl fragment missing: {csmi:?}"
        );
        assert!(
            csmi.iter().all(|s| !s.contains('.')),
            "disconnected CSMI should be split before yield: {csmi:?}"
        );
        let with_both = children.iter().any(|c| {
            c.hop.rule == "QuinoneFormation"
                && c.hop.cleaves
                && c.hop.products.len() >= 2
                && c.hop.products.iter().any(|p| canon_of(p).unwrap() == want_q)
                && c.hop.products.iter().any(|p| canon_of(p).unwrap() == want_me)
        });
        assert!(
            with_both,
            "cleaving hop must name QuinoneFormation and list both fragments"
        );
    }

    /// C16: Rust QF dealkylates phenyl–N=C=X (isocyanate / isothiocyanate /
    /// carbodiimide). Python misses the leave + quinone-imine pair — do not
    /// "fix" parity by dropping these Rust products.
    #[test]
    fn quinone_formation_dealkylates_phenyl_ncx_ahead_of_python() {
        let endpoints = qf_pair_endpoints();
        for (smi, leave, imines) in [
            ("O=C=Nc1ccccc1", "C=O", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
            ("S=C=Nc1ccccc1", "C=S", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
            ("N=C=Nc1ccccc1", "C=N", ["N=C1C=CC(=O)C=C1", "N=C1C=CC=CC1=O"]),
        ] {
            let mol = parse_mol(smi).unwrap();
            let emissions = pair_metabolize(&mol, &endpoints).unwrap();
            let products: std::collections::BTreeSet<String> = emissions
                .iter()
                .flat_map(|e| e.products.iter().cloned())
                .map(|s| canon_of(&s).unwrap())
                .collect();
            let want_leave = canon_of(leave).unwrap();
            assert!(
                products.contains(&want_leave),
                "{smi}: missing leave {leave}; got {products:?}"
            );
            for imine in imines {
                let want = canon_of(imine).unwrap();
                assert!(
                    products.contains(&want),
                    "{smi}: missing quinone-imine {imine}; got {products:?}"
                );
            }
        }
    }

    fn qf_pair_endpoints() -> Vec<PatternInfo> {
        quinone_formation()
            .patterns()
            .into_iter()
            .filter(|&p| matches!(p.edit, Edit::PairEndpoint(_)))
            .cloned()
            .collect()
    }

    #[test]
    fn quinone_formation_declares_exclusive_partner_on_hetero_arms() {
        let set = quinone_formation();
        let by_name: std::collections::BTreeMap<_, _> = set
            .patterns()
            .into_iter()
            .map(|p| (p.name.as_str(), p))
            .collect();
        for name in [
            "single_to_double",
            "replace_halogen",
            "iminium",
            "dealkylate",
        ] {
            let p = by_name.get(name).unwrap_or_else(|| panic!("missing {name}"));
            assert!(
                p.effect.exclusive_partner,
                "{name} should set exclusive_partner"
            );
        }
        let methide = by_name.get("methide_end").unwrap();
        assert!(
            !methide.effect.exclusive_partner,
            "methide alkyl partner stays off exclusive_partner"
        );
        assert_eq!(methide.effect.partner.as_deref(), Some("C"));
    }

    #[test]
    fn shared_exclusive_partner_detects_bridging_atom() {
        let endpoints = qf_pair_endpoints();
        let std = endpoints
            .iter()
            .find(|p| p.name == "single_to_double")
            .unwrap();
        assert!(std.effect.exclusive_partner);
        // Site map 1 = ring C; map 2 = partner heteroatom. Same partner atom
        // on both ends → refuse.
        let mut map1 = BTreeMap::new();
        map1.insert(1, 0);
        map1.insert(2, 10);
        let mut map2 = BTreeMap::new();
        map2.insert(1, 5);
        map2.insert(2, 10);
        assert!(shared_exclusive_partner(&map1, std, &map2, std));
        // Distinct partners (catechol-style) → allow.
        map2.insert(2, 11);
        assert!(!shared_exclusive_partner(&map1, std, &map2, std));
    }

    #[test]
    fn bridging_n_pair_candidates_do_not_share_exclusive_partner() {
        for smiles in [
            "c1ccc(N(C)c2ccccc2)cc1",
            "c1ccc2c(c1)Nc1ccccc1C2",
            "c1ccc2c(c1)Nc1ccccc1O2",
        ] {
            let mol = parse_mol(smiles).unwrap();
            let endpoints = qf_pair_endpoints();
            let cands = pair_candidates(&mol, &endpoints).unwrap();
            for c in &cands {
                assert!(
                    !shared_exclusive_partner(&c.map1, &c.left, &c.map2, &c.right),
                    "{smiles}: survivor {} still shares exclusive partner",
                    c.pattern_name
                );
            }
        }
    }

    #[test]
    fn catechol_identical_o_partners_still_emit() {
        let mol = parse_mol("Oc1ccccc1O").unwrap();
        let endpoints = qf_pair_endpoints();
        let cands = pair_candidates(&mol, &endpoints).unwrap();
        let saw_phenol_pair = cands.iter().any(|c| {
            c.left.effect.partner.as_deref() == Some("O")
                && c.right.effect.partner.as_deref() == Some("O")
                && c.left.effect.exclusive_partner
                && c.right.effect.exclusive_partner
        });
        // Partner string may live only after resolve_effect — check PatternInfo
        // on QF single_to_double: partner is None in Rust catalog (Python sets
        // via possibilities). Distinct O atoms: maps differ on map 2.
        let saw_distinct_o = cands.iter().any(|c| {
            c.left.name == "single_to_double"
                && c.right.name == "single_to_double"
                && c.map1.get(&2) != c.map2.get(&2)
                && c.map1.get(&2).is_some()
                && c.map2.get(&2).is_some()
        });
        assert!(
            saw_phenol_pair || saw_distinct_o,
            "expected ortho catechol pair with distinct O partners; got {:?}",
            cands
                .iter()
                .map(|c| (
                    c.pattern_name.as_str(),
                    c.map1.clone(),
                    c.map2.clone(),
                    c.left.effect.partner.clone(),
                    c.right.effect.partner.clone(),
                ))
                .collect::<Vec<_>>()
        );
        assert!(!cands.is_empty(), "catechol should still emit pair candidates");
    }

    #[test]
    fn constraint_products_are_closed_shell() {
        // No radicals / overfilled leaves — emit-path correctness (C10), and no
        // regression on Rust-ahead QF dealkylate (C16).
        let endpoints = qf_pair_endpoints();
        for smi in ["Oc1ccc(O)cc1", "COc1ccccc1", "O=C=Nc1ccccc1"] {
            let mol = parse_mol(smi).unwrap();
            for c in pair_candidates(&mol, &endpoints).unwrap() {
                for p in c.materialize(&mol).unwrap() {
                    assert!(!p.contains("[C]"), "{smi}: radical carbon in {p}");
                    assert!(!p.contains("[CH5]") && !p.contains("[CH6]"), "{smi}: bad methyl {p}");
                    assert!(!p.contains("[OH+]"), "{smi}: protonated carbonyl in {p}");
                }
            }
        }
    }
}
