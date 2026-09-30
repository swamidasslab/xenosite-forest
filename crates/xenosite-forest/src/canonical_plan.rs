//! Elementary plans: the plan **is** [`Deps`] (steps + precedes + maybe).
//!
//! No parallel `CanonicalStep` dialect. A hop emits elementary [`Step`]s;
//! [`Deps::bind`] rewrites [`PlanAtom::WillAdd`] → [`PlanAtom::AddedBy`] and
//! builds precedes. Composite leaves own a [`CanonicalPlanFn`] (Python
//! `canonical_plan`) that returns steps named after existing elementary
//! rules — not a `PlanKind` enum in search.
//!
//! Cleavage fragments discarded by the walk live on the plan as [`Maybe`]
//! (not a sibling on the path outcome, and not searched).
//!
//! Apply: [`Deps::linearizations`] → [`StepSequence::apply`] through named
//! elementary rules at resolved sites.

use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
use std::ops::Deref;

use crate::ForestError;
use crate::ForestMol;
use crate::labels::Tag;
use crate::mol::{Molecule, atom_idx, atom_usize, canon_of, canon_smiles};
use crate::pattern::Effect;

/// One atom note in a [`Step`] site.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum PlanAtom {
    /// Known index on the mol at emit / resolve time.
    Index(usize),
    /// Atom a prep step will add: element at this anchor index.
    /// Becomes [`PlanAtom::AddedBy`] when the plan is bound.
    WillAdd { element: String, at: usize },
    /// Atom created by an earlier elementary step (rule + that step's anchors).
    AddedBy { rule: String, anchors: Vec<usize> },
}

impl PlanAtom {
    pub fn index(idx: usize) -> Self {
        Self::Index(idx)
    }

    pub fn will_add(element: impl Into<String>, at: usize) -> Self {
        Self::WillAdd {
            element: element.into(),
            at,
        }
    }

    pub fn oxygen_at(at: usize) -> Self {
        Self::will_add("O", at)
    }

    /// Anchor index this note depends on (known atom or will-add site).
    pub fn anchor(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::WillAdd { at, .. } => Some(*at),
            Self::AddedBy { anchors, .. } => anchors.first().copied(),
        }
    }
}

fn plan_atom_sort_key(a: &PlanAtom) -> (u8, String, Vec<usize>) {
    match a {
        PlanAtom::Index(i) => (0, String::new(), vec![*i]),
        PlanAtom::WillAdd { element, at } => (1, element.clone(), vec![*at]),
        PlanAtom::AddedBy { rule, anchors } => (2, rule.clone(), anchors.clone()),
    }
}

/// One elementary reaction at a site. The plan language — not a search-only note.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Step {
    pub rule: String,
    pub site: Vec<PlanAtom>,
    /// Automorphism orbit of the site atoms (generator closure), sorted.
    /// Empty means unknown / not filled; treat as the resolved site indexes.
    /// Filled from ForestMol-cached gens when the plan is emitted.
    pub orbit: Vec<usize>,
}

/// Compat alias while call sites migrate.
pub type CanonicalStep = Step;

impl Step {
    pub fn new(rule: impl Into<String>, site: impl IntoIterator<Item = PlanAtom>) -> Self {
        let mut site: Vec<_> = site.into_iter().collect();
        site.sort_by_key(plan_atom_sort_key);
        site.dedup();
        Self {
            rule: rule.into(),
            site,
            orbit: Vec::new(),
        }
    }

    pub fn with_orbit(mut self, orbit: impl IntoIterator<Item = usize>) -> Self {
        let mut orbit: Vec<_> = orbit.into_iter().collect();
        orbit.sort_unstable();
        orbit.dedup();
        self.orbit = orbit;
        self
    }

    /// Origin / will-add anchors named by this step's site notes.
    pub fn anchors(&self) -> HashSet<usize> {
        self.site.iter().filter_map(PlanAtom::anchor).collect()
    }

    /// Orbit for equivalence checks: filled orbit, or resolved index anchors.
    pub fn site_orbit(&self) -> Vec<usize> {
        if !self.orbit.is_empty() {
            return self.orbit.clone();
        }
        let mut atoms: Vec<_> = self.anchors().into_iter().collect();
        atoms.sort_unstable();
        atoms
    }

    /// Same rule and site classes under passed / filled orbits.
    pub fn same_site_class(&self, other: &Self) -> bool {
        if self.rule != other.rule {
            return false;
        }
        let a = self.site_orbit();
        let b = other.site_orbit();
        match (a.first(), b.first()) {
            (Some(&ai), Some(&bi)) => crate::same_site_orbit(ai, &a, bi, &b),
            (None, None) => true,
            _ => false,
        }
    }

    /// Resolve site notes to current heavy-atom indices on `mol`.
    pub fn resolve_site(&self, mol: &Molecule) -> Result<HashSet<usize>, ForestError> {
        let mut out = HashSet::new();
        for note in &self.site {
            out.insert(resolve_atom(note, mol)?);
        }
        Ok(out)
    }

    /// Run this elementary rule at the resolved site; return product mols.
    ///
    /// Products keep atom indices from the edit (no SMILES round-trip) so
    /// later steps' anchors still resolve.
    ///
    /// Unique-edit emits one orbit representative. When the plan names a
    /// different atom in that class, remap a non-pair [`DeferredSite`] onto the
    /// **wanted** site so the edit lands on the planned atom (not the rep) —
    /// otherwise free multi-Index plans and [`Deps::reaches`] miss orbit mates.
    pub fn apply(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        let wanted = self.resolve_site(mol)?;
        let Some(rule) = crate::rules::leaf_rule(&self.rule) else {
            return Err(ForestError::Plan(format!(
                "unknown plan rule {}",
                self.rule
            )));
        };
        // Discovery only — preserve existing chematic tags so multi-hop Index
        // replay can compose src_to_new (ForestMol::new would restamp).
        let forest = ForestMol::wrap_preserving_labels(mol.clone());
        let gens = forest.atom_bond_generators();
        let n = mol.atom_count();
        let mut accepted = wanted.clone();
        for &w in &wanted {
            for i in crate::orbits::atom_orbit_with_gens(gens.as_ref(), n, w) {
                accepted.insert(i);
            }
        }
        let mut products = Vec::new();
        let mut seen = HashSet::new();
        for c in rule.candidates(&forest) {
            let c = c?;
            let c = if c.is_pair() {
                if !pair_matches_wanted(mol, &c, &accepted) {
                    continue;
                }
                c
            } else {
                match remap_deferred_site_to_wanted(&c, &wanted, gens.as_ref(), n, mol) {
                    Some(remapped) => remapped,
                    None => continue,
                }
            };
            for p in c.materialize_mols()? {
                let smi = canon_smiles(&p);
                if seen.insert(smi) {
                    products.push(p);
                }
            }
        }
        let endpoints = rule.leaf_pair_endpoints();
        if !endpoints.is_empty() {
            for pair in crate::pair_edit::compose_candidates_from_endpoints(
                std::rc::Rc::new(forest.copy_mol()),
                &endpoints,
            )? {
                if !pair_matches_wanted(mol, &pair, &accepted) {
                    continue;
                }
                for p in pair.materialize_mols()? {
                    let smi = canon_smiles(&p);
                    if seen.insert(smi) {
                        products.push(p);
                    }
                }
            }
        }
        Ok(products)
    }
}

/// Prefer exact wanted site; else accept when every mapped site_map atom is
/// wanted. Else remap a unique-edit embedding onto a wanted atom in the same
/// automorphism orbit (singleton sites, or multi-map ends retargeted together).
///
/// Pairs are not remapped here — callers match them with [`pair_matches_wanted`].
fn remap_deferred_site_to_wanted(
    c: &crate::candidate::DeferredSite,
    wanted: &HashSet<usize>,
    gens: &[crate::orbits::AtomBondGenerator],
    n_atoms: usize,
    mol: &Molecule,
) -> Option<crate::candidate::DeferredSite> {
    if c.is_pair() {
        return None;
    }
    if wanted.contains(&c.site) {
        return Some(c.clone());
    }
    if !c.mapped.is_empty() && c.mapped.values().all(|i| wanted.contains(i)) {
        return Some(c.clone());
    }
    let orbit = crate::orbits::atom_orbit_with_gens(gens, n_atoms, c.site);
    let &want = wanted.iter().find(|w| orbit.contains(w))?;
    if c.mapped.len() <= 1 {
        let mut e2 = c.clone();
        e2.site = want;
        e2.info.site = want;
        for v in e2.mapped.values_mut() {
            *v = want;
        }
        return Some(e2);
    }
    // Multi-atom SMIRKS map: move each mapped end onto the wanted atom of the
    // same element in the orbit (C→wanted carbon; O/N/S→hetero bonded to it).
    let want_z = mol.atom(atom_idx(want)).element.atomic_number();
    let mut e2 = c.clone();
    e2.site = want;
    e2.info.site = want;
    for v in e2.mapped.values_mut() {
        let z = mol.atom(atom_idx(*v)).element.atomic_number();
        if z == want_z {
            *v = want;
        } else if matches!(z, 7 | 8 | 16) {
            let hetero = wanted
                .iter()
                .copied()
                .find(|&w| {
                    mol.atom(atom_idx(w)).element.atomic_number() == z
                        && mol
                            .neighbors(atom_idx(want))
                            .any(|(n, _)| atom_usize(n) == w)
                })
                .or_else(|| {
                    mol.neighbors(atom_idx(want)).find_map(|(n, _)| {
                        let j = atom_usize(n);
                        (mol.atom(atom_idx(j)).element.atomic_number() == z).then_some(j)
                    })
                })?;
            *v = hetero;
        } else {
            return None;
        }
    }
    Some(e2)
}

fn pair_matches_wanted(
    mol: &Molecule,
    pair: &crate::candidate::DeferredSite,
    wanted: &HashSet<usize>,
) -> bool {
    let Some((a, b)) = pair.end_atoms() else {
        return false;
    };
    if wanted.len() == 2 && wanted.contains(&a) && wanted.contains(&b) {
        return true;
    }
    // Plan named partner heteroatoms (phenol O, etc.).
    let partners: HashSet<usize> = [a, b]
        .into_iter()
        .filter_map(|carbon| {
            mol.neighbors(atom_idx(carbon)).find_map(|(n, _)| {
                let z = mol.atom(n).element.atomic_number();
                (z == 8 || z == 7 || z == 16).then_some(atom_usize(n))
            })
        })
        .collect();
    partners.len() == wanted.len() && partners.iter().all(|p| wanted.contains(p))
}

fn resolve_atom(note: &PlanAtom, mol: &Molecule) -> Result<usize, ForestError> {
    match note {
        PlanAtom::Index(i) => {
            if *i >= mol.atom_count() {
                return Err(ForestError::Plan(format!("site index {i} out of range")));
            }
            Ok(*i)
        }
        PlanAtom::WillAdd { element, at } => resolve_added_element(mol, *at, element),
        PlanAtom::AddedBy { anchors, .. } => {
            let Some(&at) = anchors.first() else {
                return Err(ForestError::Plan("AddedBy with empty anchors".into()));
            };
            resolve_added_element(mol, at, "O")
                .or_else(|_| resolve_added_element(mol, at, "N"))
                .or_else(|_| resolve_added_element(mol, at, "S"))
        }
    }
}

fn resolve_added_element(mol: &Molecule, at: usize, element: &str) -> Result<usize, ForestError> {
    let z = match element {
        "O" => 8u8,
        "N" => 7,
        "S" => 16,
        "C" => 6,
        _ => {
            return Err(ForestError::Plan(format!(
                "cannot resolve added element {element}"
            )));
        }
    };
    bonded(mol, at, z).ok_or_else(|| {
        ForestError::Plan(format!(
            "no {element} bonded to atom {at} (will-add / added-by unresolved)"
        ))
    })
}

/// Fragment discarded at a cleavage. Not searched.
///
/// `site` is the cleavage site on the parent. `opens` are earlier ring-open
/// sites on the same walk, oldest first. `side` is the discarded fragment CSMI.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageSide {
    pub site: BTreeSet<usize>,
    pub side: String,
    pub opens: Vec<BTreeSet<usize>>,
}

impl CleavageSide {
    pub fn new(
        site: impl IntoIterator<Item = usize>,
        side: impl Into<String>,
        opens: impl IntoIterator<Item = BTreeSet<usize>>,
    ) -> Self {
        Self {
            site: site.into_iter().collect(),
            side: side.into(),
            opens: opens.into_iter().collect(),
        }
    }

    /// Formation site plus prior ring-opens (Python `span_sites`).
    pub fn span_sites(&self) -> Vec<&BTreeSet<usize>> {
        let mut out: Vec<&BTreeSet<usize>> = self.opens.iter().collect();
        out.push(&self.site);
        out
    }
}

/// Uncleared cleavage fragments carried on a [`Deps`] plan. Not a step.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Maybe {
    pub entries: Vec<CleavageSide>,
}

impl Maybe {
    pub fn new(entries: impl IntoIterator<Item = CleavageSide>) -> Self {
        Self {
            entries: entries.into_iter().collect(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn sides(&self) -> Vec<&str> {
        self.entries.iter().map(|e| e.side.as_str()).collect()
    }

    /// True when `side` is a discarded fragment, or `site` overlaps a bag span.
    ///
    /// The bifurcating cleavage site itself does not pass: that step is already
    /// in the required plan. A different reaction on overlapping atoms does.
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        if self.entries.is_empty() {
            return false;
        }
        if let Some(want) = side {
            let want = canon_of(want).unwrap_or_else(|_| want.to_string());
            return self.entries.iter().any(|e| e.side == want);
        }
        let Some(keys) = site else {
            return true;
        };
        if keys.is_empty() {
            return false;
        }
        for entry in &self.entries {
            if keys == &entry.site {
                continue;
            }
            for span in entry.span_sites() {
                if !keys.is_disjoint(span) {
                    return true;
                }
            }
        }
        false
    }
}

/// MS1 example: many hydroxylation / O-placing leaves as arms, `count = 3`
/// when the spectrum expects three oxygenations. Arms are elementary rule
/// names (catalog leaf names). Pattern-level OR inside a leaf stays on
/// [`crate::pattern::PatternInfo`] — not duplicated here.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ApplyN {
    /// Allowed elementary rule names (OR). Order is not significance.
    pub arms: Vec<String>,
    /// Exact number of applications required from `arms`.
    pub count: u16,
}

impl ApplyN {
    pub fn new(arms: impl IntoIterator<Item = impl Into<String>>, count: u16) -> Self {
        let mut arms: Vec<_> = arms.into_iter().map(Into::into).collect();
        arms.sort();
        arms.dedup();
        Self { arms, count }
    }

    pub fn allows(&self, rule: &str) -> bool {
        self.arms.iter().any(|a| a == rule)
    }

    /// Orbit-deduped unordered site combinations of length [`Self::count`].
    ///
    /// `eligible` must already include unique-edit orbit atoms (not just
    /// representatives) — see [`eligible_sites_for_apply_n`]. Counts and
    /// dedup for ApplyN read this set, not sequential unique-edit alone.
    pub fn site_combinations(&self, mol: &Molecule, eligible: &[usize]) -> Vec<Vec<usize>> {
        crate::orbits::unordered_site_combinations(mol, eligible, self.count as usize)
    }

    /// Number of orbit-deduped site combinations ([`Self::site_combinations`]).
    pub fn n_combinations(&self, mol: &Molecule, eligible: &[usize]) -> usize {
        self.site_combinations(mol, eligible).len()
    }
}

/// Eligible site atoms for an [`ApplyN`] pool: union of unique-edit / candidate
/// orbits for arms that fire on `mol` under `ruleset`.
pub fn eligible_sites_for_apply_n(
    mol: &ForestMol,
    ruleset: &crate::ruleset::RuleSet,
    pool: &ApplyN,
) -> Result<Vec<usize>, ForestError> {
    let mut eligible = BTreeSet::new();
    for cand in ruleset.candidates(mol) {
        let cand = cand?;
        let rule = cand.leaf_rule().unwrap_or_else(|| cand.pattern_name());
        if !pool.allows(rule) {
            continue;
        }
        if cand.orbit.is_empty() {
            eligible.insert(cand.site);
        } else {
            eligible.extend(cand.orbit.iter().copied());
        }
    }
    Ok(eligible.into_iter().collect())
}

/// One distinct ApplyN product plus the step plans that cover paths to it.
///
/// Each [`Deps`] has free elementary steps (no precedes) for one Aut-deduped
/// site combination. [`Deps::linearizations`] are the ordered paths; several
/// plans appear only when distinct site combos collapse to the same CSMI.
#[derive(Clone, Debug)]
pub struct ApplyNProduct {
    pub smiles: String,
    pub plans: Vec<Deps>,
}

impl ApplyNProduct {
    /// Sum of [`Deps::n_linearizations`] over covering plans.
    pub fn n_covering_linearizations(&self) -> usize {
        self.plans.iter().map(Deps::n_linearizations).sum()
    }
}

/// Counts from [`apply_n_emit_products`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ApplyNEmitStats {
    pub n_eligible_sites: usize,
    pub n_combinations: usize,
    pub n_products: usize,
    pub n_plans: usize,
    /// Σ plan linearizations (paths covered across all products).
    pub n_covering_linearizations: usize,
}

/// Emit distinct products for `pool`, each with covering step plans.
///
/// **Efficient:** one materialize pass per Aut-deduped site combination (sorted
/// tag order), not `k!` apply orders. Paths to a product are the free-plan
/// linearizations on its covering [`Deps`] (and any extra combos that share
/// the CSMI).
pub fn apply_n_emit_products(
    reactant: &str,
    ruleset: &crate::ruleset::RuleSet,
    pool: &ApplyN,
) -> Result<(Vec<ApplyNProduct>, ApplyNEmitStats), ForestError> {
    let start = ForestMol::parse(reactant)?;
    let eligible = eligible_sites_for_apply_n(&start, ruleset, pool)?;
    let combos = pool.site_combinations(start.mol(), &eligible);
    let mut by_csmi: BTreeMap<String, Vec<Deps>> = BTreeMap::new();

    for combo in &combos {
        let mut tags: Vec<Tag> = combo
            .iter()
            .map(|&i| {
                start
                    .tag_of(i)
                    .ok_or_else(|| ForestError::Plan(format!("no tag on eligible site index {i}")))
            })
            .collect::<Result<_, _>>()?;
        tags.sort_unstable();
        let Some((smiles, steps)) = apply_combo_sorted(&start, ruleset, pool, &tags)? else {
            continue;
        };
        // Free plan: no precedes — linearizations cover all step orders.
        let plan = Deps::bind(steps);
        let entry = by_csmi.entry(smiles).or_default();
        if !entry.iter().any(|p| p.same_linearizations(&plan)) {
            entry.push(plan);
        }
    }

    let products: Vec<ApplyNProduct> = by_csmi
        .into_iter()
        .map(|(smiles, plans)| ApplyNProduct { smiles, plans })
        .collect();
    let stats = ApplyNEmitStats {
        n_eligible_sites: eligible.len(),
        n_combinations: combos.len(),
        n_products: products.len(),
        n_plans: products.iter().map(|p| p.plans.len()).sum(),
        n_covering_linearizations: products
            .iter()
            .map(ApplyNProduct::n_covering_linearizations)
            .sum(),
    };
    Ok((products, stats))
}

/// Apply `tags` in given order; return product CSMI and elementary steps.
///
/// Sites are tracked across [`DeferredSite::apply`] / [`ForestMol::from_edit_product`]
/// via tags. Emitted plans use stamp-origin Indices so replay is index-stable.
/// [`PlanAtom::Index`] notes use **start** indices (carbons keep indices when
/// O is appended). Step orbits are atom indices from
/// [`crate::orbits::atom_orbit_with_gens`] on the start mol — not tag ids.
fn apply_combo_sorted(
    start: &ForestMol,
    ruleset: &crate::ruleset::RuleSet,
    pool: &ApplyN,
    tags: &[Tag],
) -> Result<Option<(String, Vec<Step>)>, ForestError> {
    let mut cur = start.clone();
    let mut steps = Vec::with_capacity(tags.len());
    let start_gens = start.atom_bond_generators();
    let start_n = start.mol().atom_count();
    for &tag in tags {
        let Some(idx) = cur.index_of(tag) else {
            return Ok(None);
        };
        let Some(start_idx) = start.index_of(tag) else {
            return Ok(None);
        };
        let gens = cur.atom_bond_generators();
        let n = cur.mol().atom_count();
        let mut applied = false;
        for cand in ruleset.candidates(&cur) {
            let cand = cand?;
            let rule = cand.leaf_rule().unwrap_or_else(|| cand.pattern_name());
            if !pool.allows(rule) {
                continue;
            }
            let wanted: HashSet<usize> = [idx].into_iter().collect();
            let Some(cand) =
                remap_deferred_site_to_wanted(&cand, &wanted, gens.as_ref(), n, cur.mol())
            else {
                continue;
            };
            let Some(em) = cand.apply()? else {
                continue;
            };
            let piece = if cand.effect.cleaves && em.products.len() > 1 {
                // Keep the largest fragment (MS1 continue side); leave goes to Maybe.
                em.products
                    .into_iter()
                    .max_by_key(|p| p.mol().atom_count())
                    .expect("non-empty")
            } else {
                em.products.into_iter().next().expect("non-empty")
            };
            let orbit =
                crate::orbits::atom_orbit_with_gens(start_gens.as_ref(), start_n, start_idx);
            steps.push(Step::new(rule, [PlanAtom::index(start_idx)]).with_orbit(orbit));
            cur = piece;
            applied = true;
            break;
        }
        if !applied {
            return Ok(None);
        }
    }
    Ok(Some((cur.csmi().as_ref().to_string(), steps)))
}

/// Sorted unique product CSMIs for `pool` (see [`apply_n_emit_products`]).
pub fn apply_n_distinct_products(
    reactant: &str,
    ruleset: &crate::ruleset::RuleSet,
    pool: &ApplyN,
) -> Result<Vec<String>, ForestError> {
    let (products, _) = apply_n_emit_products(reactant, ruleset, pool)?;
    Ok(products.into_iter().map(|p| p.smiles).collect())
}

/// Number of distinct products from [`apply_n_distinct_products`].
pub fn apply_n_n_distinct_products(
    reactant: &str,
    ruleset: &crate::ruleset::RuleSet,
    pool: &ApplyN,
) -> Result<usize, ForestError> {
    Ok(apply_n_emit_products(reactant, ruleset, pool)?.1.n_products)
}

/// Compact set of linearizations under step + precedes constraints.
///
/// Python calls this `StepPlan`. Each linear extension is a [`StepSequence`]
/// (ordered steps that act like one step). Enumerate with [`Self::linearizations`]
/// / [`Self::n_linearizations`] without always listing them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Deps {
    steps: Vec<Step>,
    precedes: Vec<(usize, usize)>,
    /// Discarded cleavage fragments (Python `PathOutcome.maybe`, on the plan).
    maybe: Maybe,
    /// OR-pools with exact apply counts (MS1 / spectrum constraints).
    apply_n: Vec<ApplyN>,
}

impl Deps {
    /// Build from already-bound steps + raw precedes; stores the reduction.
    pub fn new(
        steps: impl IntoIterator<Item = Step>,
        precedes: impl IntoIterator<Item = (usize, usize)>,
    ) -> Self {
        let steps: Vec<_> = steps.into_iter().collect();
        let n = steps.len();
        let raw: Vec<_> = precedes.into_iter().collect();
        for &(a, b) in &raw {
            assert!(
                a < n && b < n && a != b,
                "invalid precedes ({a}, {b}) for n={n}"
            );
        }
        let precedes = if n == 0 {
            Vec::new()
        } else {
            canonical_dependency_edges(n, &raw).expect("cycle in precedes")
        };
        Self {
            steps,
            precedes,
            maybe: Maybe::default(),
            apply_n: Vec::new(),
        }
    }

    /// Bind will-add notes → added-by, then build precedes (the plan identity).
    pub fn bind(steps: impl IntoIterator<Item = Step>) -> Self {
        bind_deps(steps.into_iter().collect())
    }

    /// Attach cleavage-side bags (Python `Maybe` on the outcome, here on Deps).
    pub fn with_maybe(mut self, maybe: Maybe) -> Self {
        self.maybe = maybe;
        self
    }

    /// Attach OR-apply pools (MS1: N hydroxylations from a long arm list).
    pub fn with_apply_n(mut self, pools: impl IntoIterator<Item = ApplyN>) -> Self {
        self.apply_n = pools.into_iter().collect();
        self
    }

    pub fn steps(&self) -> &[Step] {
        &self.steps
    }

    pub fn precedes(&self) -> &[(usize, usize)] {
        &self.precedes
    }

    pub fn maybe(&self) -> &Maybe {
        &self.maybe
    }

    pub fn apply_n(&self) -> &[ApplyN] {
        &self.apply_n
    }

    /// Delegate to [`Maybe::allows`] (site overlap or discarded side SMILES).
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        self.maybe.allows(site, side)
    }

    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Every topological sort under precedes (tiny graphs).
    pub fn linearizations(&self) -> Vec<StepSequence> {
        let n = self.steps.len();
        if n == 0 {
            return vec![StepSequence { steps: Vec::new() }];
        }
        let mut outgoing = vec![Vec::new(); n];
        let mut indeg = vec![0usize; n];
        for &(a, b) in &self.precedes {
            outgoing[a].push(b);
            indeg[b] += 1;
        }
        let mut out = Vec::new();
        let mut path = Vec::new();
        let mut indeg_work = indeg.clone();
        fn rec(
            steps: &[Step],
            outgoing: &[Vec<usize>],
            indeg: &mut [usize],
            path: &mut Vec<usize>,
            out: &mut Vec<StepSequence>,
        ) {
            if path.len() == steps.len() {
                out.push(StepSequence {
                    steps: path.iter().map(|&i| steps[i].clone()).collect(),
                });
                return;
            }
            let ready: Vec<_> = (0..steps.len()).filter(|&i| indeg[i] == 0).collect();
            for i in ready {
                indeg[i] = usize::MAX; // mark used
                path.push(i);
                for &b in &outgoing[i] {
                    indeg[b] -= 1;
                }
                rec(steps, outgoing, indeg, path, out);
                for &b in &outgoing[i] {
                    indeg[b] += 1;
                }
                path.pop();
                indeg[i] = 0;
            }
        }
        rec(&self.steps, &outgoing, &mut indeg_work, &mut path, &mut out);
        out
    }

    /// True if some linearization apply reaches `target` CSMI (or canon spelling).
    ///
    /// Reactant is parsed as a [`ForestMol`] so atom tags survive edits;
    /// [`StepSequence::apply_forest`] remaps Index notes through `src_to_new`
    /// (same adopt path as find_path).
    pub fn reaches(&self, reactant: &str, target: &str) -> Result<bool, ForestError> {
        let want = canon_of(target)?;
        let start = ForestMol::parse(reactant)?;
        for lin in self.linearizations() {
            for product in lin.apply_forest(&start)? {
                if product.csmi().as_ref() == want {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// Sorted unique product CSMIs from every linearization apply on `reactant`.
    ///
    /// Empty apply results (sites that fail to rebind) contribute nothing.
    /// Pair with [`Self::n_linearizations`]: `n_lin` can exceed distinct products
    /// when several orders yield the same CSMI, or exceed zero while products are
    /// empty when replay fails.
    pub fn distinct_products(&self, reactant: &str) -> Result<Vec<String>, ForestError> {
        let start = ForestMol::parse(reactant)?;
        let mut seen = BTreeSet::new();
        for lin in self.linearizations() {
            for product in lin.apply_forest(&start)? {
                seen.insert(product.csmi().as_ref().to_string());
            }
        }
        Ok(seen.into_iter().collect())
    }

    /// [`Self::distinct_products`] length.
    pub fn n_distinct_products(&self, reactant: &str) -> Result<usize, ForestError> {
        Ok(self.distinct_products(reactant)?.len())
    }

    /// `(n_linearizations, n_distinct_products)` for reporting / fuzz checks.
    pub fn replay_stats(&self, reactant: &str) -> Result<(usize, usize), ForestError> {
        Ok((self.n_linearizations(), self.n_distinct_products(reactant)?))
    }

    /// Number of topological sorts under precedes.
    ///
    /// Counts via bitmask DP (no materializing [`StepSequence`]s). `n > 20`
    /// falls back to enumerating [`Self::linearizations`] (plans that large are
    /// not expected in PhaseOne multipath).
    pub fn n_linearizations(&self) -> usize {
        let n = self.steps.len();
        if n == 0 {
            return 1;
        }
        if n > 20 {
            return self.linearizations().len();
        }
        count_topological_sorts(n, &self.precedes)
    }

    /// True iff `other` admits exactly the same total orders.
    ///
    /// Requires the same multiset of leaf [`Step`]s (exact rule + site notes —
    /// unique-edit already emits one canonical site per class, so orbit-aware
    /// align is not needed), then compares canonical precedes after aligning
    /// indices. Prefer this over `==` when declaration order of steps may
    /// differ; construction already stores the transitive reduction, so `==`
    /// also sees reduced edges when step order matches.
    pub fn same_linearizations(&self, other: &Deps) -> bool {
        let Some(aligned) = align_deps_indices(&self.steps, &other.steps) else {
            return false;
        };
        let mut edges_b: Vec<_> = other
            .precedes
            .iter()
            .map(|&(a, b)| (aligned[a], aligned[b]))
            .collect();
        edges_b.sort_unstable();
        self.precedes == edges_b
    }

    /// Yield-key when site indices remapped across free-step reorderings: same
    /// rule multiset, same Maybe side CSMI multiset, and precedes isomorphic
    /// under rule-name alignment. Complements exact [`Self::same_linearizations`]
    /// when Index notes differ only because intermediates renumbered.
    pub fn same_rule_maybe_skeleton(&self, other: &Deps) -> bool {
        if self.steps.len() != other.steps.len() {
            return false;
        }
        let mut sides_a: Vec<_> = self.maybe.sides();
        let mut sides_b: Vec<_> = other.maybe.sides();
        sides_a.sort_unstable();
        sides_b.sort_unstable();
        if sides_a != sides_b {
            return false;
        }
        let Some(aligned) = align_deps_indices_by_rule(&self.steps, &other.steps) else {
            return false;
        };
        let mut edges_b: Vec<_> = other
            .precedes
            .iter()
            .map(|&(a, b)| (aligned[a], aligned[b]))
            .collect();
        edges_b.sort_unstable();
        self.precedes == edges_b
    }

    /// True when `other` is a longer walk that only adds steps beyond `self`
    /// (Maybe sides of `self` ⊆ `other`; each of `self`'s steps matches a
    /// distinct step of `other` by rule name). Drops dominated multipath hits.
    pub fn dominates_extension_of(&self, other: &Deps) -> bool {
        if self.steps.len() >= other.steps.len() {
            return false;
        }
        let mut sides_a: Vec<_> = self.maybe.sides();
        let mut sides_b: Vec<_> = other.maybe.sides();
        sides_a.sort_unstable();
        sides_b.sort_unstable();
        if !sides_a.iter().all(|s| sides_b.contains(s)) {
            return false;
        }
        let mut used = vec![false; other.steps.len()];
        for step in &self.steps {
            let found = other
                .steps
                .iter()
                .enumerate()
                .find_map(|(j, o)| (!used[j] && step.rule == o.rule).then_some(j));
            let Some(j) = found else {
                return false;
            };
            used[j] = true;
        }
        true
    }

    /// Count of shared total orders: `|L(self) ∩ L(other)|`.
    ///
    /// Same step multiset required (else `0`). After aligning indices, the
    /// intersection of topological sorts is the sorts of the **edge union**;
    /// a cycle in that union means empty intersection (`0`). Equal to
    /// [`Self::n_linearizations`] on both sides iff [`Self::same_linearizations`].
    pub fn linearization_overlap(&self, other: &Deps) -> usize {
        let Some(aligned) = align_deps_indices(&self.steps, &other.steps) else {
            return 0;
        };
        let n = self.steps.len();
        if n == 0 {
            return 1;
        }
        let mut edges = self.precedes.clone();
        for &(a, b) in &other.precedes {
            edges.push((aligned[a], aligned[b]));
        }
        let Ok(reduced) = canonical_dependency_edges(n, &edges) else {
            return 0;
        };
        // Same nodes + union edges; maybe / apply_n do not affect required orders.
        Deps {
            steps: self.steps.clone(),
            precedes: reduced,
            maybe: Maybe::default(),
            apply_n: Vec::new(),
        }
        .n_linearizations()
    }
}

/// Count topological sorts of a DAG on `0..n` (bitmask DP).
fn count_topological_sorts(n: usize, precedes: &[(usize, usize)]) -> usize {
    debug_assert!(n <= 20);
    let mut preds = vec![0u32; n];
    for &(a, b) in precedes {
        if a < n && b < n {
            preds[b] |= 1u32 << a;
        }
    }
    let full = 1usize << n;
    let mut dp = vec![0u128; full];
    dp[0] = 1;
    for mask in 0..full {
        let ways = dp[mask];
        if ways == 0 {
            continue;
        }
        for (v, pred) in preds.iter().enumerate() {
            let bit = 1usize << v;
            if mask & bit != 0 {
                continue;
            }
            if (mask & *pred as usize) != *pred as usize {
                continue;
            }
            dp[mask | bit] = dp[mask | bit].saturating_add(ways);
        }
    }
    usize::try_from(dp[full - 1]).unwrap_or(usize::MAX)
}

/// Map indices in `steps_b` → indices in `steps_a` by [`Step`] equality.
///
/// `None` if the leaf multisets differ. Duplicate equal steps match greedily.
/// Unique-edit canonical sites make exact equality the right identity — do not
/// widen to [`Step::same_site_class`] here.
pub fn align_deps_indices(steps_a: &[Step], steps_b: &[Step]) -> Option<Vec<usize>> {
    if steps_a.len() != steps_b.len() {
        return None;
    }
    let mut used = vec![false; steps_b.len()];
    // remap[j_in_b] = i_in_a
    let mut remap = vec![0usize; steps_b.len()];
    for (i, step) in steps_a.iter().enumerate() {
        let found = steps_b
            .iter()
            .enumerate()
            .find_map(|(j, other)| (!used[j] && other == step).then_some(j));
        let j = found?;
        used[j] = true;
        remap[j] = i;
    }
    Some(remap)
}

/// Align by rule name only (greedy). Used for remapped-index yield keys.
fn align_deps_indices_by_rule(steps_a: &[Step], steps_b: &[Step]) -> Option<Vec<usize>> {
    if steps_a.len() != steps_b.len() {
        return None;
    }
    let mut used = vec![false; steps_b.len()];
    let mut remap = vec![0usize; steps_b.len()];
    for (i, step) in steps_a.iter().enumerate() {
        let found = steps_b
            .iter()
            .enumerate()
            .find_map(|(j, other)| (!used[j] && other.rule == step.rule).then_some(j));
        let j = found?;
        used[j] = true;
        remap[j] = i;
    }
    Some(remap)
}

impl Deref for Deps {
    type Target = [Step];

    fn deref(&self) -> &Self::Target {
        &self.steps
    }
}

/// Ordered elementary steps that act like a single step (composite).
///
/// One total order from a [`Deps`] / StepPlan. Enumerate via
/// [`Deps::linearizations`] (method keeps the poset “linear extension” name;
/// this type is the concrete sequence).
///
/// `apply` is the same door as for one hop: length 1 is not a special case.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepSequence {
    pub steps: Vec<Step>,
}

impl StepSequence {
    /// Apply steps in order. Returns product molecules after the last step.
    ///
    /// Soft-continues when a current fragment cannot resolve a site (cleavage
    /// leave sides). Index notes are in stamp-origin / reactant frame and are
    /// remapped through composed `src_to_new` after each hop.
    pub fn apply(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        if self.steps.is_empty() {
            return Ok(vec![mol.clone()]);
        }
        let mut currents: Vec<(Molecule, Vec<Option<usize>>)> = {
            let map: Vec<Option<usize>> = (0..mol.atom_count()).map(Some).collect();
            vec![(mol.clone(), map)]
        };
        for step in &self.steps {
            let mut next = Vec::new();
            let mut seen = HashSet::new();
            for (cur, map) in &currents {
                let Ok(mapped_step) = remap_step_through(step, map) else {
                    continue;
                };
                let Ok(prods) = mapped_step.apply(cur) else {
                    continue;
                };
                for product in prods {
                    let smi = canon_smiles(&product);
                    if !seen.insert(smi) {
                        continue;
                    }
                    let src_to_new = crate::atom_tracker::AtomTracker::src_to_new(cur, &product);
                    next.push((product, compose_index_map(map, &src_to_new)));
                }
            }
            if next.is_empty() {
                return Ok(Vec::new());
            }
            currents = next;
        }
        Ok(currents.into_iter().map(|(m, _)| m).collect())
    }

    /// Apply steps in order, adopting each product onto a [`ForestMol`].
    ///
    /// Soft-continues on site resolve failures (cleavage leaves). Stamp-origin
    /// Index notes are remapped through composed `src_to_new` after each adopt.
    pub fn apply_forest(&self, start: &ForestMol) -> Result<Vec<ForestMol>, ForestError> {
        if self.steps.is_empty() {
            return Ok(vec![start.clone()]);
        }
        let mut currents: Vec<(ForestMol, Vec<Option<usize>>)> = {
            let map: Vec<Option<usize>> = (0..start.mol().atom_count()).map(Some).collect();
            vec![(start.clone(), map)]
        };
        for step in &self.steps {
            let mut next = Vec::new();
            let mut seen = HashSet::new();
            for (cur, map) in &currents {
                let Ok(mapped_step) = remap_step_through(step, map) else {
                    continue;
                };
                let Ok(prods) = mapped_step.apply(cur.mol()) else {
                    continue;
                };
                for product in prods {
                    let child = cur.from_edit_product(product);
                    let smi = child.csmi().as_ref().to_string();
                    if !seen.insert(smi) {
                        continue;
                    }
                    let src_to_new =
                        crate::atom_tracker::AtomTracker::src_to_new(cur.mol(), child.mol());
                    next.push((child, compose_index_map(map, &src_to_new)));
                }
            }
            if next.is_empty() {
                return Ok(Vec::new());
            }
            currents = next;
        }
        Ok(currents.into_iter().map(|(m, _)| m).collect())
    }
}

fn compose_index_map(map: &[Option<usize>], src_to_new: &[Option<usize>]) -> Vec<Option<usize>> {
    map.iter()
        .map(|cur| cur.and_then(|i| src_to_new.get(i).copied().flatten()))
        .collect()
}

fn remap_step_through(step: &Step, map: &[Option<usize>]) -> Result<Step, ForestError> {
    let mut site = Vec::with_capacity(step.site.len());
    for note in &step.site {
        site.push(match note {
            PlanAtom::Index(i) => {
                let j = map.get(*i).copied().flatten().ok_or_else(|| {
                    ForestError::Plan(format!("origin index {i} left this piece"))
                })?;
                PlanAtom::index(j)
            }
            PlanAtom::WillAdd { element, at } => {
                let j = map.get(*at).copied().flatten().ok_or_else(|| {
                    ForestError::Plan(format!("will-add anchor {at} left this piece"))
                })?;
                PlanAtom::will_add(element.clone(), j)
            }
            PlanAtom::AddedBy { rule, anchors } => {
                let mut out = Vec::with_capacity(anchors.len());
                for &a in anchors {
                    let j = map.get(a).copied().flatten().ok_or_else(|| {
                        ForestError::Plan(format!("added-by anchor {a} left this piece"))
                    })?;
                    out.push(j);
                }
                PlanAtom::AddedBy {
                    rule: rule.clone(),
                    anchors: out,
                }
            }
        });
    }
    let orbit: Vec<_> = step
        .orbit
        .iter()
        .filter_map(|&i| map.get(i).copied().flatten())
        .collect();
    Ok(Step::new(step.rule.clone(), site).with_orbit(orbit))
}

/// Former name of [`StepSequence`]. Prefer `StepSequence`.
pub type Linearization = StepSequence;

/// One bitmask per node: bit `b` set iff `a` must precede `b`.
///
/// Port of Python `transitive_closure_masks`. Nodes are opaque index labels
/// `0..n-1` — this does **not** check that two plans' steps are the same
/// reactions/sites. Align node identity first (see [`Deps::same_linearizations`]).
///
/// `# Errors`
/// Invalid edge, self-loop, or cycle. `n > 128` (bitmask width).
pub fn transitive_closure_masks(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<u128>, &'static str> {
    if n > 128 {
        return Err("n > 128");
    }
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut indegree = vec![0usize; n];
    let mut seen = HashSet::new();
    for &(a, b) in edges {
        if a >= n || b >= n {
            return Err("invalid edge");
        }
        if a == b {
            return Err("cycle in precedes");
        }
        if seen.insert((a, b)) {
            outgoing[a].push(b);
            indegree[b] += 1;
        }
    }
    let mut queue: VecDeque<usize> = (0..n).filter(|&v| indegree[v] == 0).collect();
    let mut topo = Vec::with_capacity(n);
    while let Some(a) = queue.pop_front() {
        topo.push(a);
        for &b in &outgoing[a] {
            indegree[b] -= 1;
            if indegree[b] == 0 {
                queue.push_back(b);
            }
        }
    }
    if topo.len() != n {
        return Err("cycle in precedes");
    }
    let mut closure = vec![0u128; n];
    for &a in topo.iter().rev() {
        for &b in &outgoing[a] {
            // closure[a] |= (1 << b) | closure[b]
            closure[a] |= (1u128 << b) | closure[b];
        }
    }
    Ok(closure)
}

/// Unique minimal dependency edge set (transitive reduction) for this DAG.
///
/// Port of Python `canonical_dependency_edges`. Two DAGs over the **same
/// labeled nodes** `0..n-1` have exactly the same valid orderings iff this
/// returns the same edge list. [`Deps`] stores this form on construction.
pub fn canonical_dependency_edges(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<(usize, usize)>, &'static str> {
    let closure = transitive_closure_masks(n, edges)?;
    let mut canonical_edges = Vec::new();
    for a in 0..n {
        let descendants = closure[a];
        let mut through = 0u128;
        let mut remaining = descendants;
        while remaining != 0 {
            let bit = remaining & remaining.wrapping_neg();
            let b = bit.trailing_zeros() as usize;
            remaining -= bit;
            through |= closure[b];
        }
        let mut direct = descendants & !through;
        while direct != 0 {
            let bit = direct & direct.wrapping_neg();
            let b = bit.trailing_zeros() as usize;
            direct -= bit;
            canonical_edges.push((a, b));
        }
    }
    canonical_edges.sort_unstable();
    Ok(canonical_edges)
}

/// Bind will-add → added-by and collect precedes from those notes.
///
/// A [`PlanAtom::WillAdd`] depends only on earlier **prep** steps that supply
/// that element at the anchor (e.g. Hydroxylation for oxygen) — not on every
/// prior edit whose site happens to include the same index (that overstates
/// free dealk ≺ DH).
pub fn bind_deps(steps: Vec<Step>) -> Deps {
    let mut edges = Vec::new();
    let mut bound = Vec::with_capacity(steps.len());
    for (later, step) in steps.iter().enumerate() {
        let mut site = Vec::with_capacity(step.site.len());
        for item in &step.site {
            match item {
                PlanAtom::WillAdd { element, at } => {
                    let mut bound_note = PlanAtom::Index(*at);
                    for (earlier, previous) in steps[..later].iter().enumerate() {
                        if !previous.anchors().contains(at) {
                            continue;
                        }
                        if !step_supplies_will_add(previous, element) {
                            continue;
                        }
                        edges.push((earlier, later));
                        let mut anchors: Vec<_> = previous.anchors().into_iter().collect();
                        anchors.sort_unstable();
                        bound_note = PlanAtom::AddedBy {
                            rule: previous.rule.clone(),
                            anchors,
                        };
                    }
                    site.push(bound_note);
                }
                PlanAtom::AddedBy { rule, anchors } => {
                    let wanted: HashSet<_> = anchors.iter().copied().collect();
                    for (earlier, previous) in steps[..later].iter().enumerate() {
                        if previous.rule == *rule && previous.anchors() == wanted {
                            edges.push((earlier, later));
                        }
                    }
                    let mut anchors = anchors.clone();
                    anchors.sort_unstable();
                    site.push(PlanAtom::AddedBy {
                        rule: rule.clone(),
                        anchors,
                    });
                }
                PlanAtom::Index(_) => site.push(item.clone()),
            }
        }
        bound.push(Step::new(step.rule.clone(), site).with_orbit(step.orbit.iter().copied()));
    }
    Deps::new(bound, edges)
}

/// Whether `step` is a prep that can supply `element` for a [`PlanAtom::WillAdd`].
///
/// Named after the elementary rules quinone / DH / epoxide-hydration plans
/// emit as preps. Matching any earlier site that merely *touches* the anchor
/// (e.g. Dealkylation) would invent precedes between free cleaves and later DH.
fn step_supplies_will_add(step: &Step, element: &str) -> bool {
    match element {
        "O" => matches!(
            step.rule.as_str(),
            "Hydroxylation" | "OxidativeDehalogenation" | "Epoxidation"
        ),
        _ => false,
    }
}

/// [`Deps::bind`] alias (Python `as_deps`).
pub fn as_deps(steps: impl IntoIterator<Item = Step>) -> Deps {
    Deps::bind(steps)
}

/// Leaf-owned plan expander (Python `ReactionRule.canonical_plan`).
///
/// Returns elementary [`Step`]s named after catalog rules (`Hydroxylation`,
/// `Dehydrogenation`, …). `None` on a [`crate::ruleset::RuleSet`] means
/// identity: one step at the discovery site.
pub type CanonicalPlanFn = fn(
    mol: &Molecule,
    rule_name: &str,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step>;

/// Identity plan: `rule` at the given site atoms.
pub fn identity_plan(rule: impl Into<String>, site: impl IntoIterator<Item = usize>) -> Vec<Step> {
    vec![Step::new(rule, site.into_iter().map(PlanAtom::index))]
}

/// Identity plan with automorphism orbit of the site (from generators).
pub fn identity_plan_with_orbit(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
    orbit: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    let site: Vec<_> = site.into_iter().collect();
    vec![Step::new(rule, site.into_iter().map(PlanAtom::index)).with_orbit(orbit)]
}

/// Identity plan on a tagged [`ForestMol`]: root-survivor sites as stamp-origin
/// [`PlanAtom::Index`] (`Tag(k)` ↔ index `k-1` on the stamped reactant), born
/// O/N/S as [`PlanAtom::AddedBy`] at a carbon anchor (not a raw Index of the
/// born atom). Replay remaps origin Indices through `src_to_new`.
pub fn identity_plan_on_forest(
    rule: impl Into<String>,
    forest: &ForestMol,
    site_atoms: impl IntoIterator<Item = usize>,
    orbit: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    let site_atoms: Vec<usize> = site_atoms.into_iter().collect();
    let mol = forest.mol();
    let n = mol.atom_count();
    let mut carbons = Vec::new();
    let mut heteros = Vec::new();
    for &i in &site_atoms {
        if i >= n {
            continue;
        }
        let z = mol.atom(atom_idx(i)).element.atomic_number();
        match z {
            7 | 8 | 16 => heteros.push(i),
            _ => carbons.push(i),
        }
    }
    let origin_of = |i: usize| -> usize {
        forest
            .tag_of(i)
            .and_then(|t| forest.stamp_origin_index(t))
            .unwrap_or(i)
    };
    let mut site = Vec::new();
    for &c in &carbons {
        site.push(PlanAtom::index(origin_of(c)));
    }
    for &h in &heteros {
        let z = mol.atom(atom_idx(h)).element.atomic_number();
        let element = match z {
            7 => "N",
            16 => "S",
            _ => "O",
        };
        let born = forest
            .tag_of(h)
            .map(|t| forest.stamp_origin_index(t).is_none())
            .unwrap_or(true);
        if born {
            let anchor = carbons.first().copied().unwrap_or(h);
            // WillAdd so Deps::bind can wire Hydroxylation ≺ DH (AddedBy
            // rule:"O" never matches a prep step name).
            site.push(PlanAtom::will_add(element, origin_of(anchor)));
        } else {
            site.push(PlanAtom::index(origin_of(h)));
        }
    }
    if site.is_empty() {
        for &i in &site_atoms {
            if i < n {
                site.push(PlanAtom::index(origin_of(i)));
            }
        }
    }
    let orbit: Vec<usize> = orbit
        .into_iter()
        .map(origin_of)
        .filter(|&i| i < n)
        .collect();
    vec![Step::new(rule, site).with_orbit(orbit)]
}

/// Compat name.
pub fn identity_canonical_plan(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    identity_plan(rule, site)
}

/// Resolve a leaf's plan hook (or identity).
pub fn steps_for_leaf(
    plan: Option<CanonicalPlanFn>,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    match plan {
        Some(f) => f(mol, rule_name, site_atoms, end_effects),
        None => identity_plan(rule_name, site_atoms.iter().copied()),
    }
}

/// Bound [`Deps`] for one hop.
pub fn plan_for_leaf(
    plan: Option<CanonicalPlanFn>,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Deps {
    Deps::bind(steps_for_leaf(
        plan,
        rule_name,
        mol,
        site_atoms,
        end_effects,
    ))
}

fn end_needs_oxygen(effect: &Effect) -> bool {
    let partner = effect.partner.as_deref().unwrap_or("");
    let needs_o = effect.adds.as_deref().is_some_and(|a| a.contains('O'));
    needs_o && partner != "O"
}

fn bonded(mol: &Molecule, idx: usize, atomic_num: u8) -> Option<usize> {
    mol.neighbors(atom_idx(idx)).find_map(|(nbr, _)| {
        let n = atom_usize(nbr);
        (mol.atom(nbr).element.atomic_number() == atomic_num).then_some(n)
    })
}

const HALOGEN: &[&str] = &["F", "Cl", "Br", "I", "At"];

fn halogen_z(partner: &str) -> Option<u8> {
    match partner {
        "F" => Some(9),
        "Cl" => Some(17),
        "Br" => Some(35),
        "I" => Some(53),
        "At" => Some(85),
        _ => None,
    }
}

/// Preps that supply missing oxygens, then one dehydrogenation.
pub fn hydroxylation_then_dehydrogenation(
    mol: &Molecule,
    ends: &[&Effect],
    end_atoms: &[usize],
) -> Vec<Step> {
    let mut preps = Vec::new();
    let mut dh_refs = Vec::new();
    for (end, &atom) in ends.iter().zip(end_atoms.iter()) {
        let partner = end.partner.as_deref().unwrap_or("");
        if end_needs_oxygen(end) {
            let anchor = PlanAtom::index(atom);
            if HALOGEN.contains(&partner) {
                let Some(z) = halogen_z(partner) else {
                    continue;
                };
                let Some(halo) = bonded(mol, atom, z) else {
                    continue;
                };
                preps.push(Step::new(
                    "OxidativeDehalogenation",
                    [anchor.clone(), PlanAtom::index(halo)],
                ));
            } else {
                preps.push(Step::new("Hydroxylation", [anchor.clone()]));
            }
            dh_refs.push(PlanAtom::oxygen_at(atom));
            continue;
        }
        let atomic_num = match partner {
            "O" => Some(8),
            "N" => Some(7),
            "C" => Some(6),
            "S" => Some(16),
            _ => None,
        };
        if let Some(z) = atomic_num
            && let Some(hetero) = bonded(mol, atom, z)
        {
            dh_refs.push(PlanAtom::index(hetero));
        }
    }
    if dh_refs.is_empty() {
        return Vec::new();
    }
    preps.push(Step::new("Dehydrogenation", dh_refs));
    preps
}

/// Python `QuinoneFormation.canonical_plan`: prep missing oxygens, then DH.
///
/// Steps name existing elementary rules (`Hydroxylation`,
/// `OxidativeDehalogenation`, `Dehydrogenation`). Wired on the QF leaf via
/// [`crate::ruleset::RuleSet::with_canonical_plan`].
pub fn quinone_canonical_plan(
    mol: &Molecule,
    _rule_name: &str,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    if let (Some(ends), true) = (end_effects, site_atoms.len() >= 2) {
        let plan = hydroxylation_then_dehydrogenation(mol, ends, site_atoms);
        if !plan.is_empty() {
            return plan;
        }
    }
    identity_plan("Dehydrogenation", site_atoms.iter().copied())
}

/// Stable oxygenation (epoxidation) then hydrolysis (epoxide opening).
///
/// Wired on [`crate::rules::epoxide_hydration`] — one metabolize hop to the
/// vicinal diol; the plan is the elementary split for search / replay.
pub fn epoxide_hydration_canonical_plan(
    _mol: &Molecule,
    _rule_name: &str,
    site_atoms: &[usize],
    _end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    if site_atoms.len() < 2 {
        return identity_plan("EpoxideHydration", site_atoms.iter().copied());
    }
    let c0 = site_atoms[0];
    let c1 = site_atoms[1];
    vec![
        Step::new("Epoxidation", [PlanAtom::index(c0), PlanAtom::index(c1)]),
        // Both epoxide carbons plus WillAdd O: one-carbon Opening misses after
        // SMIRKS reorder on substituted alkenes; AddedBy resolves the born O.
        Step::new(
            "EpoxideOpening",
            [
                PlanAtom::index(c0),
                PlanAtom::index(c1),
                PlanAtom::oxygen_at(c0),
            ],
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::pattern::Effect;

    #[test]
    fn identity_is_one_step() {
        let plan = identity_plan("Hydroxylation", [0]);
        assert_eq!(plan.len(), 1);
        assert_eq!(plan[0].rule, "Hydroxylation");
        assert_eq!(plan[0].site, vec![PlanAtom::index(0)]);
    }

    #[test]
    fn hydroquinone_pair_needs_no_prep() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let carbons: Vec<usize> = (0..mol.atom_count())
            .filter(|&i| {
                let a = atom_idx(i);
                mol.atom(a).element.atomic_number() == 6
                    && mol
                        .neighbors(a)
                        .any(|(n, _)| mol.atom(n).element.atomic_number() == 8)
            })
            .collect();
        assert_eq!(carbons.len(), 2);
        let phenol = Effect {
            adds: None,
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: Some("O".into()),
            ..Default::default()
        };
        let ends = [&phenol, &phenol];
        let plan = hydroxylation_then_dehydrogenation(&mol, &ends, &carbons);
        assert_eq!(plan.len(), 1);
        assert_eq!(plan[0].rule, "Dehydrogenation");
        assert_eq!(plan[0].site.len(), 2);
    }

    #[test]
    fn bare_carbon_end_preps_hydroxylation() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let end = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: None,
            ..Default::default()
        };
        let ends = [&end, &end];
        let plan = hydroxylation_then_dehydrogenation(&mol, &ends, &[0, 3]);
        assert_eq!(plan.len(), 3);
        assert_eq!(plan[0].rule, "Hydroxylation");
        assert_eq!(plan[1].rule, "Hydroxylation");
        assert_eq!(plan[2].rule, "Dehydrogenation");
        assert!(
            plan[2]
                .site
                .iter()
                .any(|a| matches!(a, PlanAtom::WillAdd { .. }))
        );
    }

    #[test]
    fn epoxide_hydration_plan_is_epoxidation_then_opening() {
        let plan = epoxide_hydration_canonical_plan(
            &parse_mol("C=C").unwrap(),
            "EpoxideHydration",
            &[0, 1],
            None,
        );
        assert_eq!(plan.len(), 2);
        assert_eq!(plan[0].rule, "Epoxidation");
        assert_eq!(plan[1].rule, "EpoxideOpening");
        let deps = Deps::bind(plan);
        assert_eq!(deps.precedes(), &[(0, 1)]);
    }

    #[test]
    fn bind_benzene_qf_precedes_both_oh_before_dh() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let end = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: None,
            ..Default::default()
        };
        let ends = [&end, &end];
        let plan = hydroxylation_then_dehydrogenation(&mol, &ends, &[0, 3]);
        let deps = Deps::bind(plan);
        assert_eq!(deps.len(), 3);
        let edges: HashSet<_> = deps.precedes().iter().copied().collect();
        assert!(!edges.contains(&(0, 1)));
        assert!(!edges.contains(&(1, 0)));
        assert!(edges.contains(&(0, 2)));
        assert!(edges.contains(&(1, 2)));
        assert!(
            deps[2]
                .site
                .iter()
                .any(|a| matches!(a, PlanAtom::AddedBy { .. }))
        );
        assert_eq!(deps.linearizations().len(), 2); // OH arms commute
    }

    #[test]
    fn bind_phenol_one_oh_precedes_dh() {
        let plan = vec![
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
            Step::new(
                "Dehydrogenation",
                [PlanAtom::index(0), PlanAtom::oxygen_at(1)],
            ),
        ];
        let deps = Deps::bind(plan);
        assert_eq!(deps.precedes(), &[(0, 1)]);
    }

    #[test]
    fn bind_will_add_ignores_dealk_at_same_carbon() {
        // Free dealk at the carbon must not invent Dealk ≺ DH; only the OH prep.
        let plan = vec![
            Step::new("Dealkylation", [PlanAtom::index(0), PlanAtom::index(1)]),
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
            Step::new(
                "Dehydrogenation",
                [PlanAtom::index(2), PlanAtom::oxygen_at(1)],
            ),
        ];
        let deps = Deps::bind(plan);
        let edges: HashSet<_> = deps.precedes().iter().copied().collect();
        assert!(
            !edges.contains(&(0, 2)),
            "Dealkylation must not precede DH: {edges:?}"
        );
        assert!(edges.contains(&(1, 2)), "OH must precede DH: {edges:?}");
        assert_eq!(edges.len(), 1);
    }

    #[test]
    fn replay_benzene_qf_plan_reaches_quinone() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let end = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: None,
            ..Default::default()
        };
        let ends = [&end, &end];
        let deps = Deps::bind(hydroxylation_then_dehydrogenation(&mol, &ends, &[0, 3]));
        assert!(
            deps.reaches("c1ccccc1", "O=C1C=CC(=O)C=C1").unwrap(),
            "plan={deps:?}"
        );
    }

    /// Every QF prep-spine shape: OH/DH cases must [`Deps::reaches`]; OxDehal
    /// spines check structure + aryl/aliphatic leaf apply (aromatic specialize).
    #[test]
    fn qf_canonical_plan_shapes_all_replay() {
        fn rules_of(deps: &Deps) -> Vec<&str> {
            deps.steps().iter().map(|s| s.rule.as_str()).collect()
        }
        fn ring_c_bonded_to(mol: &Molecule, z: u8) -> usize {
            (0..mol.atom_count())
                .find(|&i| {
                    let a = atom_idx(i);
                    mol.atom(a).element.atomic_number() == 6
                        && mol.atom(a).aromatic
                        && mol
                            .neighbors(a)
                            .any(|(n, _)| mol.atom(n).element.atomic_number() == z)
                })
                .expect("ring C bonded to Z")
        }
        fn bare_para_to(mol: &Molecule, tagged: usize) -> usize {
            (0..mol.atom_count())
                .filter(|&i| {
                    if i == tagged {
                        return false;
                    }
                    let a = atom_idx(i);
                    if mol.atom(a).element.atomic_number() != 6 || !mol.atom(a).aromatic {
                        return false;
                    }
                    !mol.neighbors(a).any(|(n, _)| {
                        let z = mol.atom(n).element.atomic_number();
                        z != 6 && z != 1
                    })
                })
                .max_by_key(|&i| {
                    let mut dist = vec![usize::MAX; mol.atom_count()];
                    dist[tagged] = 0;
                    let mut q = std::collections::VecDeque::from([tagged]);
                    while let Some(u) = q.pop_front() {
                        for (n, _) in mol.neighbors(atom_idx(u)) {
                            let v = atom_usize(n);
                            if dist[v] == usize::MAX {
                                dist[v] = dist[u] + 1;
                                q.push_back(v);
                            }
                        }
                    }
                    dist[i]
                })
                .expect("para CH")
        }

        let need_o = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            dearomatizes: true,
            ..Default::default()
        };
        let phenol = Effect {
            removes: Some("H".into()),
            dearomatizes: true,
            partner: Some("O".into()),
            ..Default::default()
        };
        let cl_end = Effect {
            adds: Some("O".into()),
            dearomatizes: true,
            partner: Some("Cl".into()),
            ..Default::default()
        };

        // 1) DH only — hydroquinone.
        {
            let smi = "Oc1ccc(O)cc1";
            let mol = parse_mol(smi).unwrap();
            let c0 = ring_c_bonded_to(&mol, 8);
            let c1 = (0..mol.atom_count())
                .find(|&i| {
                    i != c0
                        && mol.atom(atom_idx(i)).element.atomic_number() == 6
                        && mol.atom(atom_idx(i)).aromatic
                        && mol
                            .neighbors(atom_idx(i))
                            .any(|(n, _)| mol.atom(n).element.atomic_number() == 8)
                })
                .unwrap();
            let deps = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&phenol, &phenol],
                &[c0, c1],
            ));
            assert_eq!(rules_of(&deps), ["Dehydrogenation"]);
            assert!(deps.reaches(smi, "O=C1C=CC(=O)C=C1").unwrap(), "{deps:?}");
        }

        // 2) OH → DH — phenol.
        {
            let smi = "Oc1ccccc1";
            let mol = parse_mol(smi).unwrap();
            let c_oh = ring_c_bonded_to(&mol, 8);
            let c_h = bare_para_to(&mol, c_oh);
            let deps = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&need_o, &phenol],
                &[c_h, c_oh],
            ));
            assert_eq!(rules_of(&deps), ["Hydroxylation", "Dehydrogenation"]);
            assert!(deps.reaches(smi, "O=C1C=CC(=O)C=C1").unwrap(), "{deps:?}");
        }

        // 3) OH → OH → DH — benzene.
        {
            let smi = "c1ccccc1";
            let mol = parse_mol(smi).unwrap();
            let deps = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&need_o, &need_o],
                &[0, 3],
            ));
            assert_eq!(
                rules_of(&deps),
                ["Hydroxylation", "Hydroxylation", "Dehydrogenation"]
            );
            assert!(deps.reaches(smi, "O=C1C=CC(=O)C=C1").unwrap());
        }

        // 4–7) OxDehal spines — structure only (cleavage remaps Indices for reaches).
        {
            let smi = "Oc1ccc(Cl)cc1";
            let mol = parse_mol(smi).unwrap();
            let c_cl = ring_c_bonded_to(&mol, 17);
            let c_oh = ring_c_bonded_to(&mol, 8);
            let deps = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&cl_end, &phenol],
                &[c_cl, c_oh],
            ));
            assert_eq!(
                rules_of(&deps),
                ["OxidativeDehalogenation", "Dehydrogenation"]
            );
        }
        {
            let smi = "Clc1ccc(Cl)cc1";
            let mol = parse_mol(smi).unwrap();
            let mut cls: Vec<usize> = (0..mol.atom_count())
                .filter(|&i| {
                    mol.atom(atom_idx(i)).element.atomic_number() == 6
                        && mol.atom(atom_idx(i)).aromatic
                        && mol
                            .neighbors(atom_idx(i))
                            .any(|(n, _)| mol.atom(n).element.atomic_number() == 17)
                })
                .collect();
            cls.sort_unstable();
            let deps = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&cl_end, &cl_end],
                &[cls[0], cls[1]],
            ));
            assert_eq!(
                rules_of(&deps),
                [
                    "OxidativeDehalogenation",
                    "OxidativeDehalogenation",
                    "Dehydrogenation"
                ]
            );
        }
        {
            let smi = "Clc1ccccc1";
            let mol = parse_mol(smi).unwrap();
            let c_cl = ring_c_bonded_to(&mol, 17);
            let c_h = bare_para_to(&mol, c_cl);
            let oh_first = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&need_o, &cl_end],
                &[c_h, c_cl],
            ));
            assert_eq!(
                rules_of(&oh_first),
                [
                    "Hydroxylation",
                    "OxidativeDehalogenation",
                    "Dehydrogenation"
                ]
            );
            let ox_first = Deps::bind(hydroxylation_then_dehydrogenation(
                &mol,
                &[&cl_end, &need_o],
                &[c_cl, c_h],
            ));
            assert_eq!(
                rules_of(&ox_first),
                [
                    "OxidativeDehalogenation",
                    "Hydroxylation",
                    "Dehydrogenation"
                ]
            );
        }

        // Aryl OxDehal leaf applies via aromatic specialize; aliphatic too.
        {
            let mol = parse_mol("Clc1ccccc1").unwrap();
            let hits =
                crate::smarts::smarts_matches(&mol, "[#9,#17,#35,#53,#85:1]-[#6:2]").unwrap();
            assert_eq!(hits.len(), 1);
            let c = *hits[0].get(&2).unwrap();
            let step = Step::new("OxidativeDehalogenation", [PlanAtom::index(c)]);
            let products = step.apply(&mol).unwrap();
            let phenol = canon_of("Oc1ccccc1").unwrap();
            assert!(
                products
                    .iter()
                    .any(|p| canon_of(&canon_smiles(p)).unwrap() == phenol),
                "aryl OxDehal got {:?}",
                products.iter().map(canon_smiles).collect::<Vec<_>>()
            );
        }
        {
            let mol = parse_mol("CCCl").unwrap();
            let hits =
                crate::smarts::smarts_matches(&mol, "[#9,#17,#35,#53,#85:1]-[#6:2]").unwrap();
            assert_eq!(hits.len(), 1);
            let c = *hits[0].get(&2).unwrap();
            let x = *hits[0].get(&1).unwrap();
            let step = Step::new(
                "OxidativeDehalogenation",
                [PlanAtom::index(c), PlanAtom::index(x)],
            );
            let products = step.apply(&mol).unwrap();
            assert!(
                products.iter().any(|p| {
                    let got = canon_of(&canon_smiles(p)).unwrap();
                    got == canon_of("C(C)O").unwrap() || got == canon_of("CCO").unwrap()
                }),
                "got {:?}",
                products.iter().map(canon_smiles).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn replay_hydroxylation_ethane() {
        let mol = crate::as_forest_mol("CC").unwrap();
        let rule = crate::rules::hydroxylation();
        let em = rule
            .metabolize(&mol, |_, _, _| true, |_, _, _| true, true)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!em.is_empty());
        let deps = Deps::bind(identity_plan("Hydroxylation", [em[0].site]));
        assert!(deps.reaches("CC", "CCO").unwrap());
    }

    #[test]
    fn canonical_edges_drop_transitive() {
        let reduced = canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert_eq!(reduced, vec![(0, 1), (1, 2)]);
    }

    #[test]
    fn transitive_closure_and_canonical() {
        // Python test_transitive_closure_and_canonical
        let closure = transitive_closure_masks(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert!(closure[0] & (1 << 1) != 0 && closure[0] & (1 << 2) != 0);
        assert!(closure[1] & (1 << 2) != 0);
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap(),
            vec![(0, 1), (1, 2)]
        );
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2)]).unwrap(),
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap()
        );
    }

    #[test]
    fn transitive_closure_cycle_raises() {
        assert!(
            transitive_closure_masks(2, &[(0, 1), (1, 0)])
                .unwrap_err()
                .contains("cycle")
        );
        assert!(
            transitive_closure_masks(2, &[(0, 0)])
                .unwrap_err()
                .contains("cycle")
        );
    }

    #[test]
    fn deps_same_linearizations_via_transitive_closure() {
        // Python test_deps_same_linearizations_via_transitive_closure:
        // lin-set identity is canonical edges — not == / raw precedes.
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        let chain = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2)]);
        let with_transitive =
            Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2), (0, 2)]);
        assert_eq!(chain, with_transitive);
        assert_eq!(chain.precedes(), &[(0, 1), (1, 2)]);
        assert_eq!(with_transitive.precedes(), &[(0, 1), (1, 2)]);
        assert!(chain.same_linearizations(&with_transitive));
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap(),
            vec![(0, 1), (1, 2)]
        );

        let flipped = Deps::new([c.clone(), a.clone(), b.clone()], [(1, 2), (2, 0)]); // a≺b≺c
        assert!(chain.same_linearizations(&flipped));
        assert_ne!(chain, flipped); // == is order-sensitive

        let layered = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 2), (1, 2)]);
        assert!(layered.same_linearizations(&Deps::new(
            [a.clone(), b.clone(), c.clone()],
            [(0, 2), (1, 2)]
        )));

        let free_dealk = Deps::new([a.clone(), b.clone(), c.clone()], [(1, 2)]); // only b≺c
        assert!(!free_dealk.same_linearizations(&chain));

        let other_nodes = Deps::new(
            [
                Step::new("X", [PlanAtom::index(0)]),
                Step::new("Y", [PlanAtom::index(1)]),
                Step::new("Z", [PlanAtom::index(2)]),
            ],
            [(0, 1), (1, 2)],
        );
        assert!(!chain.same_linearizations(&other_nodes));
    }

    #[test]
    fn deps_stores_canonical_precedes() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        let d = Deps::new([a, b, c], [(0, 1), (1, 2), (0, 2)]);
        assert_eq!(d.precedes(), &[(0, 1), (1, 2)]);
    }

    #[test]
    fn same_linearizations_two_prep_then_final() {
        let h0 = Step::new("Hydroxylation", [PlanAtom::index(0)]);
        let h3 = Step::new("Hydroxylation", [PlanAtom::index(3)]);
        let dh = Step::new(
            "Dehydrogenation",
            [
                PlanAtom::AddedBy {
                    rule: "Hydroxylation".into(),
                    anchors: vec![0],
                },
                PlanAtom::AddedBy {
                    rule: "Hydroxylation".into(),
                    anchors: vec![3],
                },
            ],
        );
        let layered = Deps::new([h0.clone(), h3.clone(), dh.clone()], [(0, 2), (1, 2)]);
        let swapped = Deps::new([h3, h0, dh], [(0, 2), (1, 2)]);
        assert!(layered.same_linearizations(&swapped));
        assert_eq!(layered.n_linearizations(), 2);
        assert_eq!(layered.linearization_overlap(&swapped), 2);
    }

    #[test]
    fn linearization_overlap_partial_and_conflicting() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        // Free A∥B ≺ C → 2 orders
        let free = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 2), (1, 2)]);
        // Chain A≺B≺C → 1 order (subset of free)
        let chain = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2)]);
        assert_eq!(free.linearization_overlap(&chain), 1);
        assert_eq!(chain.linearization_overlap(&free), 1);
        assert!(!free.same_linearizations(&chain));

        // Opposite A/B orders → union cycles → 0
        let ab = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1)]);
        let ba = Deps::new([a.clone(), b.clone(), c.clone()], [(1, 0)]);
        assert_eq!(ab.linearization_overlap(&ba), 0);

        // Different step multiset → 0
        let other = Deps::new(
            [a, b, Step::new("D", [PlanAtom::index(2)])],
            [(0, 2), (1, 2)],
        );
        assert_eq!(free.linearization_overlap(&other), 0);
        assert_eq!(
            Deps::new([], []).linearization_overlap(&Deps::new([], [])),
            1
        );
    }

    #[test]
    fn same_linearizations_duplicate_equal_steps() {
        // Python test_align_duplicate_steps_greedy
        let a1 = Step::new("A", [PlanAtom::index(0)]);
        let a2 = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let d1 = Deps::new([a1.clone(), a2.clone(), b.clone()], [(0, 2), (1, 2)]);
        let d2 = Deps::new([a2, b, a1], [(0, 1), (2, 1)]);
        assert!(d1.same_linearizations(&d2));
    }

    #[test]
    fn free_nodes_n_linearizations_is_factorial() {
        let steps: Vec<_> = (0..4)
            .map(|i| Step::new(format!("S{i}"), [PlanAtom::index(i)]))
            .collect();
        assert_eq!(Deps::new(steps, []).n_linearizations(), 24);
        assert_eq!(Deps::new([], []).n_linearizations(), 1);
        assert_eq!(
            Deps::new([Step::new("A", [PlanAtom::index(0)])], []).n_linearizations(),
            1
        );
    }

    #[test]
    fn same_rule_maybe_skeleton_collapses_remapped_free_dealks() {
        let a = Step::new("Dealkylation", [PlanAtom::index(0)]);
        let b = Step::new("Dealkylation", [PlanAtom::index(9)]);
        let c = Step::new("Dealkylation", [PlanAtom::index(12)]);
        let d = Step::new("Dealkylation", [PlanAtom::index(0)]);
        let maybe = Maybe::new([
            CleavageSide::new([0], "OC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new([9], "OC", std::iter::empty::<BTreeSet<usize>>()),
        ]);
        let p0 = Deps::new([a, b], []).with_maybe(maybe.clone());
        let p1 = Deps::new([c, d], []).with_maybe(maybe);
        assert!(!p0.same_linearizations(&p1)); // exact sites differ
        assert!(p0.same_rule_maybe_skeleton(&p1));
    }

    #[test]
    fn dominates_extension_drops_longer_same_maybe_walk() {
        let short = Deps::new(
            [
                Step::new("Dealkylation", [PlanAtom::index(8)]),
                Step::new("Dealkylation", [PlanAtom::index(1)]),
            ],
            [],
        )
        .with_maybe(Maybe::new([
            CleavageSide::new([8], "CNC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new(
                [1],
                "c1ccc(C(C)(C)C)cc1",
                std::iter::empty::<BTreeSet<usize>>(),
            ),
        ]));
        let longer = Deps::new(
            [
                Step::new("Dealkylation", [PlanAtom::index(8)]),
                Step::new("Dealkylation", [PlanAtom::index(1)]),
                Step::new("Dehydrogenation", [PlanAtom::index(6)]),
            ],
            [],
        )
        .with_maybe(Maybe::new([
            CleavageSide::new([8], "CNC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new(
                [1],
                "c1ccc(C(C)(C)C)cc1",
                std::iter::empty::<BTreeSet<usize>>(),
            ),
        ]));
        assert!(short.dominates_extension_of(&longer));
        assert!(!longer.dominates_extension_of(&short));
        assert!(!short.same_linearizations(&longer));
    }

    #[test]
    fn n_linearizations_matches_enumeration_on_layered() {
        let h0 = Step::new("Hydroxylation", [PlanAtom::index(0)]);
        let h3 = Step::new("Hydroxylation", [PlanAtom::index(3)]);
        let dh = Step::new("Dehydrogenation", [PlanAtom::index(0)]);
        let layered = Deps::new([h0, h3, dh], [(0, 2), (1, 2)]);
        assert_eq!(layered.n_linearizations(), layered.linearizations().len());
        assert_eq!(layered.n_linearizations(), 2);
    }

    #[test]
    fn replay_stats_free_hydroxylations_share_product() {
        // Two free OH steps (no precedes): 2 linearizations, one product when
        // both orders reach the same CSMI on ethane.
        let plan = Deps::bind([
            Step::new("Hydroxylation", [PlanAtom::index(0)]),
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
        ]);
        assert_eq!(plan.n_linearizations(), 2);
        let (n_lin, n_prod) = plan.replay_stats("CC").unwrap();
        assert_eq!(n_lin, 2);
        assert_eq!(n_prod, plan.n_distinct_products("CC").unwrap());
        assert_eq!(
            n_prod,
            1,
            "both orders → same ethane diol: {:?}",
            plan.distinct_products("CC")
        );
    }

    #[test]
    fn apply_n_benzene_oh2_three_combinations_three_products() {
        let set = crate::rules::hydroxylation();
        let pool = ApplyN::new(["Hydroxylation"], 2);
        let mol = ForestMol::parse("c1ccccc1").unwrap();
        let eligible = eligible_sites_for_apply_n(&mol, &set, &pool).unwrap();
        assert_eq!(eligible.len(), 6);
        assert_eq!(pool.n_combinations(mol.mol(), &eligible), 3);
        let (products, stats) = apply_n_emit_products("c1ccccc1", &set, &pool).unwrap();
        assert_eq!(stats.n_combinations, 3);
        assert_eq!(stats.n_products, 3, "{products:?}");
        assert_eq!(products.len(), 3);
        // Each product: one free 2-step plan → 2 linearizations (path orders).
        for p in &products {
            assert_eq!(p.plans.len(), 1, "{}", p.smiles);
            assert_eq!(p.n_covering_linearizations(), 2, "{}", p.smiles);
            assert!(
                p.plans[0].reaches("c1ccccc1", &p.smiles).unwrap(),
                "{} plan={:?}",
                p.smiles,
                p.plans[0]
            );
        }
        assert_eq!(stats.n_covering_linearizations, 6);
        assert_eq!(
            apply_n_n_distinct_products("c1ccccc1", &set, &pool).unwrap(),
            3
        );
    }

    #[test]
    fn apply_n_emit_covers_paths_without_permuting_applies() {
        // Ethane OH×2: one combo, one product; free plan has 2 lins.
        let set = crate::rules::hydroxylation();
        let pool = ApplyN::new(["Hydroxylation"], 2);
        let (products, stats) = apply_n_emit_products("CC", &set, &pool).unwrap();
        assert_eq!(stats.n_combinations, 1);
        assert_eq!(
            stats.n_products,
            1,
            "{:?}",
            products.iter().map(|p| &p.smiles).collect::<Vec<_>>()
        );
        assert_eq!(products[0].n_covering_linearizations(), 2);
        let (n_lin, n_prod) = products[0].plans[0].replay_stats("CC").unwrap();
        assert_eq!(n_lin, 2);
        assert_eq!(n_prod, 1);
    }

    #[test]
    fn deps_invalid_edge_panics() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        assert!(
            std::panic::catch_unwind(|| {
                let _ = Deps::new([a.clone(), b.clone()], [(0, 5)]);
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                let _ = Deps::new([a, b], [(0, 1), (1, 0)]);
            })
            .is_err()
        );
    }

    #[test]
    fn apply_n_composes_with_maybe_on_deps() {
        let steps = Deps::new([Step::new("Dealkylation", [PlanAtom::index(0)])], [])
            .with_maybe(Maybe::new([CleavageSide::new(
                [0],
                "C",
                std::iter::empty::<BTreeSet<usize>>(),
            )]))
            .with_apply_n([ApplyN::new(["Hydroxylation", "Epoxidation"], 2)]);
        assert_eq!(steps.apply_n().len(), 1);
        assert_eq!(steps.apply_n()[0].count, 2);
        assert!(steps.apply_n()[0].allows("Hydroxylation"));
        assert!(steps.apply_n()[0].allows("Epoxidation"));
        assert!(!steps.apply_n()[0].allows("Dealkylation"));
        assert!(!steps.maybe().is_empty());
        assert_eq!(steps.len(), 1);
    }

    /// Key ops: ApplyN beside precedes (WillAdd bind) and Maybe — linearizations,
    /// reaches, and bag attachment stay consistent.
    #[test]
    fn apply_n_composes_with_precedes_maybe_reaches() {
        let mol = parse_mol("C=C").unwrap();
        let plan = epoxide_hydration_canonical_plan(&mol, "EpoxideHydration", &[0, 1], None);
        let deps = Deps::bind(plan)
            .with_maybe(Maybe::new([CleavageSide::new(
                [0],
                "leave",
                std::iter::empty::<BTreeSet<usize>>(),
            )]))
            .with_apply_n([ApplyN::new(["EpoxideHydration"], 1)]);
        assert_eq!(deps.precedes(), &[(0, 1)]);
        assert_eq!(deps.n_linearizations(), 1);
        assert!(!deps.maybe().is_empty());
        assert_eq!(deps.apply_n().len(), 1);
        assert_eq!(deps.apply_n()[0].count, 1);
        assert!(deps.apply_n()[0].allows("EpoxideHydration"));
        assert!(deps.reaches("C=C", "OCCO").unwrap());
    }

    #[test]
    fn apply_n_epoxide_hydration_plan_replays_after_hydrogenation() {
        // Multi-hop substrate: alkyne → ene → diol. Opening names both carbons
        // so Step.apply Index layout still resolves.
        let h2_then_hyd = Deps::bind([
            Step::new("Hydrogenation", [PlanAtom::index(0), PlanAtom::index(1)]),
            Step::new("Epoxidation", [PlanAtom::index(0), PlanAtom::index(1)]),
            Step::new(
                "EpoxideOpening",
                [
                    PlanAtom::index(0),
                    PlanAtom::index(1),
                    PlanAtom::oxygen_at(0),
                ],
            ),
        ])
        .with_apply_n([ApplyN::new(
            ["Hydrogenation", "EpoxideHydration", "Epoxidation"],
            2,
        )]);
        assert_eq!(h2_then_hyd.precedes(), &[(1, 2)]);
        assert!(h2_then_hyd.reaches("C#C", "OCCO").unwrap());
        assert_eq!(h2_then_hyd.apply_n()[0].count, 2);
    }

    /// Alcohol DH after epoxide hydration: identity from site_atoms (C+O when
    /// the deferred site maps both ends) must replay to the aldehyde.
    #[test]
    fn ethene_hyd_then_alcohol_dh_pair_identity_replays() {
        let start = ForestMol::parse("C=C").unwrap();
        let hyd = crate::rules::epoxide_hydration();
        let c = hyd.candidates(&start).next().unwrap().unwrap();
        let glycol = c.apply().unwrap().unwrap().products.remove(0);
        let dh = crate::rules::dehydrogenation();
        let alcohol = dh
            .candidates(&glycol)
            .find(|c| c.as_ref().is_ok_and(|c| c.pattern_name() == "alcohol"))
            .unwrap()
            .unwrap();
        let mut steps =
            epoxide_hydration_canonical_plan(start.mol(), "EpoxideHydration", &[0, 1], None);
        // Stamp-origin Indices + AddedBy for born alcohol O (not raw Index of O).
        let atoms = alcohol.site_atoms();
        steps.extend(identity_plan_on_forest(
            alcohol
                .leaf_rule()
                .unwrap_or_else(|| alcohol.pattern_name()),
            &glycol,
            atoms.clone(),
            if alcohol.is_pair() {
                atoms
            } else {
                alcohol.orbit.clone()
            },
        ));
        let deps = Deps::bind(steps);
        assert_eq!(deps.steps()[2].rule, "Dehydrogenation");
        assert!(
            !deps.steps()[2].site.is_empty(),
            "alcohol identity site: {:?}",
            deps.steps()[2]
        );
        assert!(deps.reaches("C=C", "OCC=O").unwrap(), "plan={deps:?}");
    }

    #[test]
    fn apply_n_dedups_and_sorts_arms() {
        let pool = ApplyN::new(["Epoxidation", "Hydroxylation", "Hydroxylation"], 3);
        assert_eq!(
            pool.arms,
            vec!["Epoxidation".to_string(), "Hydroxylation".to_string()]
        );
    }
}
