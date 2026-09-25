//! First-run plan-guided search over a [`crate::ruleset::RuleSet`].
//!
//! Walks carry tagged [`crate::forest_mol::ForestMol`] (structure/`csmi` cache,
//! atom tags through edits). Expands via [`RuleSet::candidates`]: site–pattern–
//! parent triples. A search reads [`PatternInfo`] / effect fields to decide,
//! then materializes only for survivors. Nested sets stay namespaces on each
//! step's leaf-first `rule_path`.
//!
//! Outcomes carry [`crate::canonical_plan::Deps`] plans (elementary steps +
//! precedes + [`crate::canonical_plan::Maybe`] cleavage bags). Composite leaves
//! own a `canonical_plan` hook (Python) that returns steps named after existing
//! rules. Closer uses atom-diff cost. Lazy: try tag-lift at enqueue (same-heavy-
//! tag edits only); full MCS on pop when lift is `None`. Eager:
//! [`crate::atom_diff::atom_diff_for_child`] (lift else MCS) at enqueue.

use std::cmp::Ordering;
use std::collections::{BTreeSet, BinaryHeap, HashSet};
use std::rc::Rc;

use crate::ForestError;
use crate::candidate::Candidate;
use crate::canonical_plan::{CanonicalStep, CleavageSide, Deps, Maybe, as_deps};
use crate::forest_mol::ForestMol;
use crate::mol::{canon_of, parse_mol};
use crate::pattern::{CleaveFoldKey, CleaveSideSig, PatternInfo, SiteInfo};
use crate::rules::default_ruleset;
use crate::ruleset::RuleSet;

/// Billed work for one search (Python `PathCounters` subset).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PathCounters {
    pub nodes: usize,
    pub mol_edits: usize,
    pub expansions: usize,
    /// Hit not emitted: same linearizations or same rule/Maybe skeleton.
    pub dropped_duplicate_plan: usize,
    /// Hit that *would* match [`Deps::dominates_extension_of`] against a yield.
    /// Counted only — not used to drop or abort (that lever over-collapsed
    /// multipath; HEURISTICS).
    pub signal_contained_plan: usize,
}

impl PathCounters {
    pub fn billed(&self) -> usize {
        self.mol_edits + self.nodes
    }

    /// Optimization signal: duplicate yield drops + contained-extension matches.
    /// Do not use to abort walks; more careful heuristics might reduce later.
    pub fn plan_drops(&self) -> usize {
        self.dropped_duplicate_plan + self.signal_contained_plan
    }
}

/// One accepted edit on a walk: RuleSet namespace path + pattern + kept product.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathStep {
    pub rule_path: Vec<Option<String>>,
    pub pattern_name: String,
    pub site: usize,
    /// Unique-edit primary-map orbit (includes `site`). Plan equivalence under
    /// automorphism compares sites via [`crate::same_site_orbit`].
    pub site_orbit: Vec<usize>,
    /// Kept fragment CSMI (the search node).
    pub product: String,
    /// Cleaved-off fragment CSMIs (not expanded).
    pub sides: Vec<String>,
}

impl PathStep {
    /// Named segments of [`Self::rule_path`] (unnamed sets omitted).
    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }

    /// Same rule/pattern and sites in one unique-edit orbit.
    pub fn same_site_class(&self, other: &Self) -> bool {
        self.pattern_name == other.pattern_name
            && self.leaf_rule() == other.leaf_rule()
            && crate::same_site_orbit(self.site, &self.site_orbit, other.site, &other.site_orbit)
    }

    fn from_emission(emission: &ForestEmission, product: String, sides: Vec<String>) -> Self {
        Self {
            rule_path: emission.rule_path.clone(),
            pattern_name: emission.pattern_name.clone(),
            site: emission.site,
            site_orbit: emission.site_orbit.clone(),
            product,
            sides,
        }
    }
}

/// One reactant→target hit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathOutcome {
    pub steps: Vec<PathStep>,
    /// Elementary plan with precedes and cleavage [`Maybe`] (on the plan).
    pub plan: Deps,
    pub smiles: String,
}

impl PathOutcome {
    /// Cleavage-side bags on the plan (Python `PathOutcome.maybe`).
    pub fn maybe(&self) -> &Maybe {
        self.plan.maybe()
    }

    /// [`Deps::allows`] / [`Maybe::allows`] — discarded side or overlapping site.
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        self.plan.allows(site, side)
    }
}

#[derive(Clone)]
struct Walk {
    mol: ForestMol,
    steps: Vec<PathStep>,
    plan: Vec<CanonicalStep>,
    /// Accumulated cleavage bags (attached to [`Deps`] at hit yield).
    maybe: Vec<CleavageSide>,
    /// Uncleared ring-open sites (Python `_Walk.opens`).
    opens: Vec<BTreeSet<usize>>,
    /// Parent's [`crate::atom_diff::AtomDiff::cost`] when this walk was
    /// enqueued. `None` = root (always expand).
    parent_cost: Option<usize>,
    /// Diff of **this** mol vs target. Filled at enqueue by
    /// [`try_atom_diff_for_child`] / cleavage lift when safe; else `None`
    /// and pop runs full MCS.
    diff: Option<crate::atom_diff::AtomDiff>,
}

/// Heap entry: hits first, then FIFO (`seq`). Lower priority value pops first.
#[derive(Clone)]
struct HeapItem {
    target_hit: bool,
    seq: usize,
    walk: Walk,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.target_hit == other.target_hit && self.seq == other.seq
    }
}

impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.target_hit
            .cmp(&other.target_hit)
            .then_with(|| other.seq.cmp(&self.seq))
    }
}

fn ha_distance(ha: usize, target_ha: usize) -> usize {
    ha.abs_diff(target_ha)
}

/// Tagged emission: ForestMol products + elementary plan + cleavage site data.
#[derive(Clone)]
struct ForestEmission {
    site: usize,
    site_orbit: Vec<usize>,
    /// Discovery site atoms (Python frozenset site for CleavageSide).
    site_atoms: BTreeSet<usize>,
    cleaves: bool,
    /// Cross-rule Or fold signature (from [`PatternInfo::cleave_side_group`]).
    cleave_side_sig: CleaveSideSig,
    pattern_name: String,
    rule_path: Vec<Option<String>>,
    products: Vec<ForestMol>,
    plan: Vec<CanonicalStep>,
}

impl ForestEmission {
    fn product_csmis(&self) -> Vec<String> {
        self.products
            .iter()
            .map(|p| p.csmi().as_ref().to_string())
            .collect()
    }

    fn cleave_fold_key(&self) -> CleaveFoldKey {
        self.cleave_side_sig.fold_key(&self.product_csmis())
    }
}

/// Keep fragments that are worth continuing toward ``target``.
///
/// On a bifurcation with a parent diff, every product that passes the cleavage
/// expand gate (strict cost drop via tag-lift when possible, ha ≥ target) is
/// kept — **both** sides if both match. Uses tagged [`ForestMol`] products;
/// never re-parses SMILES.
fn keep_fragments(
    parent: &ForestMol,
    products: &[ForestMol],
    target_csmi: &str,
    target_ha: usize,
    parent_diff: Option<&crate::atom_diff::AtomDiff>,
    target_mol: Option<&crate::Molecule>,
) -> Vec<(ForestMol, Vec<String>, Option<crate::atom_diff::AtomDiff>)> {
    if products.is_empty() {
        return Vec::new();
    }
    if products.len() >= 2 {
        if let (Some(pdiff), Some(tmol)) = (parent_diff, target_mol) {
            let mut kept = Vec::new();
            for (i, child) in products.iter().enumerate() {
                let csmi = child.csmi();
                let is_hit = csmi.as_ref() == target_csmi;
                if !is_hit && child.heavy_atom_count() < target_ha {
                    continue;
                }
                let child_diff =
                    crate::atom_diff::atom_diff_after_cleavage(parent, pdiff, child, tmol);
                if is_hit || child_diff.cost() < pdiff.cost() {
                    let mut sides = Vec::new();
                    for (j, other) in products.iter().enumerate() {
                        if j != i && other.csmi().as_ref() != csmi.as_ref() {
                            sides.push(other.csmi().as_ref().to_string());
                        }
                    }
                    kept.push((child.clone(), sides, Some(child_diff)));
                }
            }
            if !kept.is_empty() {
                return kept;
            }
        }
    }
    keep_fragment(products, target_csmi, target_ha)
        .into_iter()
        .map(|(mol, sides)| (mol, sides, None))
        .collect()
}

/// Keep the fragment closest to ``target``. Prefer exact CSMI hit.
fn keep_fragment(
    products: &[ForestMol],
    target_csmi: &str,
    target_ha: usize,
) -> Option<(ForestMol, Vec<String>)> {
    if products.is_empty() {
        return None;
    }
    let mut best: Option<(usize, Rc<str>, usize)> = None;
    for (i, mol) in products.iter().enumerate() {
        let csmi = mol.csmi();
        let cost = if csmi.as_ref() == target_csmi {
            0
        } else {
            1 + ha_distance(mol.heavy_atom_count(), target_ha)
        };
        match &best {
            None => best = Some((i, Rc::clone(&csmi), cost)),
            Some((_, _, best_cost)) if cost < *best_cost => {
                best = Some((i, Rc::clone(&csmi), cost));
            }
            Some((_, kept, best_cost)) if cost == *best_cost && csmi.as_ref() < kept.as_ref() => {
                best = Some((i, Rc::clone(&csmi), cost));
            }
            _ => {}
        }
    }
    let (best_i, kept_csmi, _) = best?;
    let kept = products[best_i].clone();
    let mut sides = Vec::new();
    for (i, mol) in products.iter().enumerate() {
        if i == best_i {
            continue;
        }
        let csmi = mol.csmi();
        if csmi.as_ref() != kept_csmi.as_ref() {
            sides.push(csmi.as_ref().to_string());
        }
    }
    Some((kept, sides))
}

fn closer(parent_ha: usize, child_ha: usize, target_ha: usize, target_hit: bool) -> bool {
    // Refuse only walks that grow more distant in heavy-atom count. Equal
    // distance (same-size DH / oxidation hops) must stay open — a strict `<`
    // drops propane→propene→epoxide. When atom_diff is on, this is the cheap
    // enqueue gate; full cost closer runs on pop (lazy).
    target_hit || ha_distance(child_ha, target_ha) <= ha_distance(parent_ha, target_ha)
}

fn cost_closer(parent_cost: usize, child_cost: usize, target_hit: bool) -> bool {
    target_hit || child_cost < parent_cost
}

/// Search bounds. Defaults match Python `find_path` knobs.
#[derive(Clone, Copy, Debug)]
pub struct FindPathConfig {
    pub max_paths: usize,
    pub max_nodes: usize,
    /// When true, build an [`crate::atom_diff::AtomDiff`] once per expansion
    /// and refuse candidates that cannot help (no filter closures).
    pub use_atom_diff: bool,
    /// When true with `use_atom_diff`, defer full cost closer to pop: enqueue
    /// every kept fragment, verify `cost() < parent_cost` before expanding.
    /// Skips MCS on siblings never popped. Do not HA-gate at enqueue —
    /// oxidation can raise HA distance while lowering cost. Default on.
    pub lazy_closer: bool,
}

impl Default for FindPathConfig {
    fn default() -> Self {
        Self {
            max_paths: 1,
            max_nodes: 800,
            // Match Python live `use_filters=True`.
            use_atom_diff: true,
            lazy_closer: true,
        }
    }
}

fn accept_all_candidates(_c: &Candidate) -> bool {
    true
}

/// [`FindPath`] with the built-in accept-all keep predicate.
pub type OpenFindPath<'a, 'b> = FindPath<'a, 'b, fn(&Candidate) -> bool>;

/// HEURISTICS: a later walk that is only a reordering of an already-yielded
/// [`Deps`] is not a new path. Also drop remapped-index free-step twins
/// ([`Deps::same_rule_maybe_skeleton`]). Dominated-extension is **not** a
/// yield drop (over-collapsed multipath); it only bumps
/// [`PathCounters::signal_contained_plan`].
fn plan_already_yielded(found: &[PathOutcome], plan: &Deps) -> bool {
    found
        .iter()
        .any(|h| h.plan.same_linearizations(plan) || h.plan.same_rule_maybe_skeleton(plan))
}

fn record_yield_plan_signals(counters: &mut PathCounters, found: &[PathOutcome], plan: &Deps) {
    if plan_already_yielded(found, plan) {
        counters.dropped_duplicate_plan += 1;
        return;
    }
    if found.iter().any(|h| h.plan.dominates_extension_of(plan)) {
        counters.signal_contained_plan += 1;
    }
}

/// Yield walks that turn ``reactant`` into ``target``.
///
/// Pull iterator (Python generator parity): each [`Iterator::next`] resumes the
/// search until the next [`PathOutcome`]. Nested sets keep step namespaces on
/// each `rule_path`. Callers that want a list use `.collect()`.
pub fn find_path<'a, 'b>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'a, 'b>, ForestError> {
    find_path_with(
        reactant,
        target,
        ruleset,
        counters,
        FindPathConfig::default(),
        accept_all_candidates,
    )
}

/// [`find_path`] with the Phase-I default catalog.
pub fn find_path_default<'b>(
    reactant: &str,
    target: &str,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'static, 'b>, ForestError> {
    // Leak-free: borrow the process-wide default via once / static ruleset.
    find_path(reactant, target, default_ruleset_ref(), counters)
}

fn default_ruleset_ref() -> &'static RuleSet {
    use std::sync::OnceLock;
    static DEFAULT: OnceLock<RuleSet> = OnceLock::new();
    DEFAULT.get_or_init(default_ruleset)
}

/// [`find_path`] with atom-diff candidate gating (no filter closures).
pub fn find_path_diff<'a, 'b>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'a, 'b>, ForestError> {
    find_path_with(
        reactant,
        target,
        ruleset,
        counters,
        FindPathConfig {
            use_atom_diff: true,
            ..FindPathConfig::default()
        },
        accept_all_candidates,
    )
}

/// Same as [`find_path`], with caller bounds and a candidate predicate.
///
/// `keep` reads the triple (site, pattern, parent) **before** materialize.
/// Dropped candidates do not count as `mol_edits`.
pub fn find_path_with<'a, 'b, K>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    config: FindPathConfig,
    keep: K,
) -> Result<FindPath<'a, 'b, K>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let start = ForestMol::parse(reactant)?;
    let start_csmi = start.csmi();
    let target_csmi = canon_of(target)?;
    let target_mol = parse_mol(&target_csmi)?;
    let target_ha = target_mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi.as_ref() == target_csmi.as_str(),
        seq,
        walk: Walk {
            mol: start,
            steps: Vec::new(),
            plan: Vec::new(),
            maybe: Vec::new(),
            opens: Vec::new(),
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi.as_ref().to_string());

    Ok(FindPath {
        ruleset,
        counters,
        keep,
        config,
        target_csmi,
        target_mol,
        target_ha,
        heap,
        seq,
        seen,
        yielded: Vec::new(),
        done: false,
    })
}

/// Pull search over a [`RuleSet`]: yields one [`PathOutcome`] per `next`.
///
/// Expand pulls candidates, buffers survivors only for `order_key` sort, then
/// materializes one emission at a time (never a full emission list).
pub struct FindPath<'a, 'b, K> {
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    keep: K,
    config: FindPathConfig,
    target_csmi: String,
    target_mol: crate::Molecule,
    target_ha: usize,
    heap: BinaryHeap<HeapItem>,
    seq: usize,
    seen: HashSet<String>,
    /// Already-yielded hits (for [`plan_already_yielded`] only).
    yielded: Vec<PathOutcome>,
    done: bool,
}

impl<K> FindPath<'_, '_, K>
where
    K: Fn(&Candidate) -> bool,
{
    /// Drain the search into a list (Python `list(find_path(...))`).
    pub fn collect_all(self) -> Result<Vec<PathOutcome>, ForestError> {
        self.collect()
    }
}

impl<K> Iterator for FindPath<'_, '_, K>
where
    K: Fn(&Candidate) -> bool,
{
    type Item = Result<PathOutcome, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let FindPathConfig {
            max_paths,
            max_nodes,
            use_atom_diff,
            lazy_closer,
        } = self.config;

        while let Some(item) = self.heap.pop() {
            if self.yielded.len() >= max_paths || self.counters.nodes >= max_nodes {
                break;
            }
            let walk = item.walk;
            let here = walk.mol.csmi();
            if here.as_ref() == self.target_csmi.as_str() {
                self.counters.nodes += 1;
                let plan = as_deps(walk.plan).with_maybe(Maybe::new(walk.maybe));
                record_yield_plan_signals(self.counters, &self.yielded, &plan);
                if plan_already_yielded(&self.yielded, &plan) {
                    continue;
                }
                let outcome = PathOutcome {
                    steps: walk.steps,
                    plan,
                    smiles: here.as_ref().to_string(),
                };
                self.yielded.push(outcome.clone());
                return Some(Ok(outcome));
            }

            let diff = if use_atom_diff {
                let d = match walk.diff {
                    Some(d) => d,
                    None => crate::atom_diff::atom_diff(walk.mol.mol(), &self.target_mol),
                };
                if lazy_closer {
                    if let Some(pc) = walk.parent_cost {
                        if !cost_closer(pc, d.cost(), false) {
                            continue;
                        }
                    }
                }
                Some(d)
            } else {
                None
            };
            self.counters.nodes += 1;
            self.counters.expansions += 1;
            let parent_cost = diff.as_ref().map(|d| d.cost());
            let parent_ha = walk.mol.heavy_atom_count();
            let mut hits_from_here = 0usize;
            // Cross-rule cleave Or: enqueue each (fold_key, continue_csmi) once.
            let mut seen_cleave_continues: HashSet<(CleaveFoldKey, String)> = HashSet::new();

            let expand = match Expand::new(
                self.ruleset,
                &walk.mol,
                &self.target_mol,
                self.counters,
                &self.keep,
                diff.as_ref(),
            ) {
                Ok(e) => e,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            };

            for emission in expand {
                let emission = match emission {
                    Ok(e) => e,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                let keeps = keep_fragments(
                    &walk.mol,
                    &emission.products,
                    &self.target_csmi,
                    self.target_ha,
                    diff.as_ref(),
                    Some(&self.target_mol),
                );
                let cleave_key = if emission.cleaves && emission.products.len() >= 2 {
                    Some(emission.cleave_fold_key())
                } else {
                    None
                };
                for (kept, sides, lifted_diff) in keeps {
                    let kept_csmi = kept.csmi().as_ref().to_string();
                    if let Some(key) = &cleave_key {
                        if !seen_cleave_continues.insert((key.clone(), kept_csmi.clone())) {
                            // Another arm (possibly another rule) already queued
                            // this continuation for the same Or bag.
                            continue;
                        }
                    }
                    let child_ha = kept.heavy_atom_count();
                    let target_hit = kept_csmi == self.target_csmi;

                    let mut child_diff = if use_atom_diff {
                        lifted_diff.or_else(|| {
                            diff.as_ref().and_then(|parent_d| {
                                crate::atom_diff::try_atom_diff_for_child(
                                    &walk.mol,
                                    parent_d,
                                    &kept,
                                    &self.target_mol,
                                )
                            })
                        })
                    } else {
                        None
                    };

                    let allow = if use_atom_diff {
                        if lazy_closer {
                            true
                        } else if let Some(pc) = parent_cost {
                            if child_diff.is_none() {
                                child_diff =
                                    Some(crate::atom_diff::atom_diff(kept.mol(), &self.target_mol));
                            }
                            cost_closer(pc, child_diff.as_ref().unwrap().cost(), target_hit)
                        } else {
                            true
                        }
                    } else {
                        closer(parent_ha, child_ha, self.target_ha, target_hit)
                    };
                    if !allow {
                        continue;
                    }
                    if self.seen.contains(&kept_csmi) && !target_hit {
                        continue;
                    }
                    self.seen.insert(kept_csmi.clone());

                    let (child_maybe, child_opens) = accumulate_maybe(
                        &walk.maybe,
                        &walk.opens,
                        &emission.site_atoms,
                        emission.cleaves,
                        emission.products.len(),
                        &sides,
                    );
                    let mut steps = walk.steps.clone();
                    steps.push(PathStep::from_emission(&emission, kept_csmi.clone(), sides));
                    let mut plan = walk.plan.clone();
                    plan.extend(emission.plan.iter().cloned());
                    self.heap.push(HeapItem {
                        target_hit,
                        seq: self.seq,
                        walk: Walk {
                            mol: kept,
                            steps,
                            plan,
                            maybe: child_maybe,
                            opens: child_opens,
                            parent_cost,
                            diff: child_diff,
                        },
                    });
                    self.seq += 1;
                    if target_hit {
                        hits_from_here += 1;
                        if self.yielded.len() + hits_from_here >= max_paths {
                            break;
                        }
                    }
                }
                if self.yielded.len() + hits_from_here >= max_paths {
                    break;
                }
            }
        }

        self.done = true;
        None
    }
}

/// Python walk update: ring-open accumulates `opens`; bifurcation appends bags.
fn accumulate_maybe(
    parent_maybe: &[CleavageSide],
    parent_opens: &[BTreeSet<usize>],
    site_atoms: &BTreeSet<usize>,
    cleaves: bool,
    n_products: usize,
    discarded: &[String],
) -> (Vec<CleavageSide>, Vec<BTreeSet<usize>>) {
    if n_products == 1 && cleaves {
        let mut opens = parent_opens.to_vec();
        opens.push(site_atoms.clone());
        (parent_maybe.to_vec(), opens)
    } else if n_products > 1 {
        let mut maybe = parent_maybe.to_vec();
        for side in discarded {
            maybe.push(CleavageSide::new(
                site_atoms.iter().copied(),
                side.clone(),
                parent_opens.iter().cloned(),
            ));
        }
        (maybe, parent_opens.to_vec())
    } else {
        (parent_maybe.to_vec(), parent_opens.to_vec())
    }
}

fn candidate_site_atoms(candidate: &Candidate) -> BTreeSet<usize> {
    let mut atoms: BTreeSet<usize> = candidate
        .pattern
        .site_map
        .iter()
        .filter_map(|m| candidate.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.insert(candidate.site);
    }
    atoms
}

/// One DFS frame while streaming ResonancePair emissions.
struct PairFrame<'a> {
    set: &'a RuleSet,
    member_i: usize,
}

/// Pending pair after filter + order_key sort (discovery only).
struct PendingPair<'a> {
    pair: crate::pair_edit::PairCandidate,
    set: &'a RuleSet,
    rule_path: Vec<Option<String>>,
}

/// Pull tagged emissions one at a time.
///
/// Candidate survivors may be buffered for `order_key` sort (discovery, no
/// edit). Materialize runs per `next`. Pair leaves are walked depth-first the
/// same way — never a full emission `Vec`.
struct Expand<'a, K> {
    parent: &'a ForestMol,
    target: &'a crate::Molecule,
    counters: &'a mut PathCounters,
    keep: &'a K,
    diff: Option<&'a crate::atom_diff::AtomDiff>,
    deferred: std::vec::IntoIter<Candidate>,
    pair_stack: Vec<PairFrame<'a>>,
    pair_pending: std::vec::IntoIter<PendingPair<'a>>,
    done: bool,
}

impl<'a, K> Expand<'a, K>
where
    K: Fn(&Candidate) -> bool,
{
    fn new(
        ruleset: &'a RuleSet,
        parent: &'a ForestMol,
        target: &'a crate::Molecule,
        counters: &'a mut PathCounters,
        keep: &'a K,
        diff: Option<&'a crate::atom_diff::AtomDiff>,
    ) -> Result<Self, ForestError> {
        let mol = parent.mol();
        // Pull candidates; buffer only survivors for order_key sort.
        let mut deferred = Vec::new();
        for c in ruleset.candidates(mol) {
            let c = c?;
            let ok = if let Some(d) = diff {
                keep(&c)
                    && crate::atom_diff::candidate_could_help_on(&c, d, Some(mol), Some(target))
            } else {
                keep(&c)
            };
            if ok {
                deferred.push(c);
            }
        }
        if let Some(d) = diff {
            deferred.sort_by_key(|c| crate::atom_diff::candidate_order_key(c, d));
        }
        Ok(Self {
            parent,
            target,
            counters,
            keep,
            diff,
            deferred: deferred.into_iter(),
            pair_stack: vec![PairFrame {
                set: ruleset,
                member_i: 0,
            }],
            pair_pending: Vec::new().into_iter(),
            done: false,
        })
    }

    fn emit_candidate(
        &mut self,
        candidate: &Candidate,
    ) -> Result<Option<ForestEmission>, ForestError> {
        let pieces = candidate.materialize_mols(self.parent.mol())?;
        if pieces.is_empty() {
            return Ok(None);
        }
        let products: Vec<ForestMol> = pieces
            .into_iter()
            .map(|piece| self.parent.adopt_product(piece))
            .collect();
        Ok(Some(ForestEmission {
            site: candidate.site,
            site_orbit: candidate.orbit.clone(),
            site_atoms: candidate_site_atoms(candidate),
            cleaves: candidate.pattern.effect.cleaves,
            cleave_side_sig: candidate.pattern.cleave_side_sig(),
            pattern_name: candidate.pattern.name.clone(),
            rule_path: candidate.rule_path.clone(),
            products,
            plan: candidate.identity_plan_with_gens(
                &self.parent.atom_bond_generators(),
                self.parent.mol().atom_count(),
            ),
        }))
    }

    fn emit_pair(
        &mut self,
        pending: &PendingPair<'a>,
    ) -> Result<Option<ForestEmission>, ForestError> {
        let mol = self.parent.mol();
        let pair = &pending.pair;
        // Python ResonancePair bumps mol_edits only after alternating paths exist.
        let pieces = pair.materialize_mols(mol)?;
        if pieces.is_empty() {
            return Ok(None);
        }
        self.counters.mol_edits += 1;
        let products: Vec<ForestMol> = pieces
            .into_iter()
            .map(|piece| self.parent.adopt_product(piece))
            .collect();
        let site_atoms = pair.plan_site_atoms();
        let ends = [&pair.left.effect, &pair.right.effect];
        let plan = pending.set.canonical_plan(mol, &site_atoms, Some(&ends));
        Ok(Some(ForestEmission {
            site: pair.site,
            site_orbit: vec![pair.site],
            site_atoms: pair.plan_site_atoms().into_iter().collect(),
            cleaves: pair.effect.cleaves,
            cleave_side_sig: {
                let left = pair.left.cleave_side_sig();
                let right = pair.right.cleave_side_sig();
                if left == right {
                    left
                } else {
                    CleaveSideSig::Ungrouped
                }
            },
            pattern_name: pair.pattern_name.clone(),
            rule_path: pending.rule_path.clone(),
            products,
            plan,
        }))
    }

    fn load_leaf_pairs(&mut self, set: &'a RuleSet) -> Result<(), ForestError> {
        let mol = self.parent.mol();
        let mut pairs = set.pair_candidates_leaf(mol)?;
        if let Some(d) = self.diff {
            pairs.retain(|pair| {
                if !keep_pair(pair, self.keep) {
                    return false;
                }
                crate::atom_diff::pair_could_help(pair, d, mol, self.target)
            });
            pairs.sort_by_key(|p| {
                let cleave = if p.effect.cleaves { 0u8 } else { 1 };
                let dear = if p.effect.dearomatizes { 0u8 } else { 1 };
                let oxy = if p.effect.adds.as_deref().is_some_and(|a| a.contains('O')) {
                    0u8
                } else {
                    1
                };
                let want_cleave = d.target_smaller() || d.has_cleavage();
                let want_dear = !d.loses_aromaticity.is_empty();
                let want_oxy = !d.needs_oxygen.is_empty();
                (
                    if want_cleave { cleave } else { 0 },
                    if want_dear { dear } else { 0 },
                    if want_oxy { oxy } else { 0 },
                    p.pattern_name.clone(),
                )
            });
        } else {
            pairs.retain(|p| keep_pair(p, self.keep));
        }
        // Leaf-first path: this set, then ancestors still on the stack (root last).
        let mut rule_path = vec![set.name.clone()];
        for frame in self.pair_stack.iter().rev() {
            rule_path.push(frame.set.name.clone());
        }
        let pending: Vec<PendingPair<'a>> = pairs
            .into_iter()
            .map(|pair| PendingPair {
                pair,
                set,
                rule_path: rule_path.clone(),
            })
            .collect();
        self.pair_pending = pending.into_iter();
        Ok(())
    }

    fn advance_pairs(&mut self) -> Result<(), ForestError> {
        loop {
            if !self.pair_pending.as_slice().is_empty() {
                return Ok(());
            }
            let action = {
                let Some(frame) = self.pair_stack.last_mut() else {
                    return Ok(());
                };
                let members = frame.set.members();
                if frame.member_i < members.len() {
                    let member = &members[frame.member_i];
                    frame.member_i += 1;
                    match member {
                        crate::ruleset::RuleMember::Set(child) => {
                            // `'a` outlives the frame; child lives in the ruleset tree.
                            let child: &'a RuleSet = child;
                            Some(PairAction::Push(child))
                        }
                        crate::ruleset::RuleMember::Pattern(_) => Some(PairAction::Continue),
                    }
                } else {
                    let set = frame.set;
                    Some(PairAction::LoadLeaf(set))
                }
            };
            match action {
                None => return Ok(()),
                Some(PairAction::Continue) => {}
                Some(PairAction::Push(child)) => {
                    self.pair_stack.push(PairFrame {
                        set: child,
                        member_i: 0,
                    });
                }
                Some(PairAction::LoadLeaf(set)) => {
                    self.pair_stack.pop();
                    self.load_leaf_pairs(set)?;
                }
            }
        }
    }
}

enum PairAction<'a> {
    Continue,
    Push(&'a RuleSet),
    LoadLeaf(&'a RuleSet),
}

impl<K> Iterator for Expand<'_, K>
where
    K: Fn(&Candidate) -> bool,
{
    type Item = Result<ForestEmission, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            while let Some(candidate) = self.deferred.next() {
                self.counters.mol_edits += 1;
                match self.emit_candidate(&candidate) {
                    Ok(Some(emission)) => return Some(Ok(emission)),
                    Ok(None) => {}
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            if let Some(pending) = self.pair_pending.next() {
                match self.emit_pair(&pending) {
                    Ok(Some(emission)) => return Some(Ok(emission)),
                    Ok(None) => continue,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            match self.advance_pairs() {
                Ok(()) => {
                    if self.pair_pending.as_slice().is_empty() && self.pair_stack.is_empty() {
                        self.done = true;
                        return None;
                    }
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

fn keep_pair<K>(pair: &crate::pair_edit::PairCandidate, keep: &K) -> bool
where
    K: Fn(&Candidate) -> bool,
{
    let mut stand_in = Candidate {
        site: pair.site,
        orbit: vec![pair.site],
        pattern: pair.left.clone(),
        rule_path: Vec::new(),
        mapped: Default::default(),
        parent: crate::candidate::ParentRef::Context,
    };
    stand_in.pattern.effect = pair.effect.clone();
    stand_in.pattern.name = pair.pattern_name.clone();
    keep(&stand_in)
}

/// Python-style filter closures (optional). Prefer [`find_path_with`] + reading
/// [`Candidate`] fields when the predicate does not need the leaf `RuleSet`.
///
/// Pull iterator: metabolize emissions are consumed one at a time (no full
/// emission list). Callers that want a list use `.collect()`.
pub fn find_path_with_filters<'a, 'b, R, S>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    config: FindPathConfig,
    filter_rules: R,
    filter_sites: S,
) -> Result<FindPathFilters<'a, 'b, R, S>, ForestError>
where
    R: Fn(&crate::Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&crate::Molecule, usize, &SiteInfo) -> bool,
{
    let start = ForestMol::parse(reactant)?;
    let start_csmi = start.csmi();
    let target_csmi = canon_of(target)?;
    let target_ha = ForestMol::parse(&target_csmi)?.heavy_atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi.as_ref() == target_csmi.as_str(),
        seq,
        walk: Walk {
            mol: start,
            steps: Vec::new(),
            plan: Vec::new(),
            maybe: Vec::new(),
            opens: Vec::new(),
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi.as_ref().to_string());

    Ok(FindPathFilters {
        ruleset,
        counters,
        filter_rules,
        filter_sites,
        config,
        target_csmi,
        target_ha,
        heap,
        seq,
        seen,
        yielded: Vec::new(),
        done: false,
    })
}

/// Filter-closure search: pulls [`RuleSet::metabolize`] one emission at a time.
pub struct FindPathFilters<'a, 'b, R, S> {
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    filter_rules: R,
    filter_sites: S,
    config: FindPathConfig,
    target_csmi: String,
    target_ha: usize,
    heap: BinaryHeap<HeapItem>,
    seq: usize,
    seen: HashSet<String>,
    yielded: Vec<PathOutcome>,
    done: bool,
}

impl<R, S> Iterator for FindPathFilters<'_, '_, R, S>
where
    R: Fn(&crate::Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&crate::Molecule, usize, &SiteInfo) -> bool,
{
    type Item = Result<PathOutcome, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let FindPathConfig {
            max_paths,
            max_nodes,
            ..
        } = self.config;

        while let Some(item) = self.heap.pop() {
            if self.yielded.len() >= max_paths || self.counters.nodes >= max_nodes {
                break;
            }
            self.counters.nodes += 1;
            let walk = item.walk;
            let here = walk.mol.csmi();
            if here.as_ref() == self.target_csmi.as_str() {
                let plan = as_deps(walk.plan).with_maybe(Maybe::new(walk.maybe));
                record_yield_plan_signals(self.counters, &self.yielded, &plan);
                if plan_already_yielded(&self.yielded, &plan) {
                    continue;
                }
                let outcome = PathOutcome {
                    steps: walk.steps,
                    plan,
                    smiles: here.as_ref().to_string(),
                };
                self.yielded.push(outcome.clone());
                return Some(Ok(outcome));
            }

            let mol = walk.mol.mol();
            self.counters.expansions += 1;
            let mut hits_from_here = 0usize;

            // Pull metabolize one emission at a time — no full list.
            let emissions =
                self.ruleset
                    .metabolize(mol, &self.filter_rules, &self.filter_sites, true);
            for emission in emissions {
                let emission = match emission {
                    Ok(e) => e,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                self.counters.mol_edits += 1;
                // Filter-path emissions are CSMI strings (no tag continuity). Prefer
                // find_path_with + Candidate materialize for tagged walks.
                let products: Result<Vec<_>, _> = emission
                    .products
                    .iter()
                    .map(|s| ForestMol::parse(s))
                    .collect();
                let products = match products {
                    Ok(p) => p,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                let keeps = keep_fragments(
                    &walk.mol,
                    &products,
                    &self.target_csmi,
                    self.target_ha,
                    None,
                    None,
                );
                for (kept, sides, _) in keeps {
                    let kept_csmi = kept.csmi().as_ref().to_string();
                    let child_ha = kept.heavy_atom_count();
                    let target_hit = kept_csmi == self.target_csmi;
                    if !closer(
                        walk.mol.heavy_atom_count(),
                        child_ha,
                        self.target_ha,
                        target_hit,
                    ) {
                        continue;
                    }
                    if self.seen.contains(&kept_csmi) && !target_hit {
                        continue;
                    }
                    self.seen.insert(kept_csmi.clone());

                    let site_atoms: BTreeSet<usize> = emission.site_atoms.iter().copied().collect();
                    let (child_maybe, child_opens) = accumulate_maybe(
                        &walk.maybe,
                        &walk.opens,
                        &site_atoms,
                        emission.cleaves,
                        products.len(),
                        &sides,
                    );
                    let mut steps = walk.steps.clone();
                    steps.push(PathStep {
                        rule_path: emission.rule_path.clone(),
                        pattern_name: emission.pattern_name.clone(),
                        site: emission.site,
                        site_orbit: emission.site_orbit.clone(),
                        product: kept_csmi.clone(),
                        sides,
                    });
                    let mut plan = walk.plan.clone();
                    plan.extend(emission.plan.iter().cloned());
                    self.heap.push(HeapItem {
                        target_hit,
                        seq: self.seq,
                        walk: Walk {
                            mol: kept,
                            steps,
                            plan,
                            maybe: child_maybe,
                            opens: child_opens,
                            parent_cost: None,
                            diff: None,
                        },
                    });
                    self.seq += 1;
                    if target_hit {
                        hits_from_here += 1;
                        if self.yielded.len() + hits_from_here >= max_paths {
                            break;
                        }
                    }
                }
                if self.yielded.len() + hits_from_here >= max_paths {
                    break;
                }
            }
        }

        self.done = true;
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydroxylation::hydroxylation;
    use crate::mol::canon_of;
    use crate::ruleset::o_dealkylation;

    #[test]
    fn ethane_to_ethanol_is_one_hydroxylation() {
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let outcome = &hits[0];
        assert_eq!(outcome.smiles, canon_of("CCO").unwrap());
        assert_eq!(outcome.steps.len(), 1);
        assert_eq!(outcome.steps[0].leaf_rule(), Some("Hydroxylation"));
        assert_eq!(outcome.steps[0].namespace(), vec!["Hydroxylation"]);
        assert!(outcome.steps[0].sides.is_empty());
        assert_eq!(outcome.plan.len(), 1);
        assert_eq!(outcome.plan[0].rule, "Hydroxylation");
        assert_eq!(counters.mol_edits, 1);
    }

    #[test]
    fn anisole_to_phenol_cleaves_and_records_side() {
        let mut counters = PathCounters::default();
        let hits = find_path("COc1ccccc1", "Oc1ccccc1", &o_dealkylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let outcome = &hits[0];
        assert_eq!(outcome.smiles, canon_of("Oc1ccccc1").unwrap());
        assert_eq!(outcome.steps.len(), 1);
        assert_eq!(outcome.steps[0].leaf_rule(), Some("Dealkylation"));
        assert_eq!(outcome.steps[0].pattern_name, "O-Me");
        assert!(
            !outcome.steps[0].sides.is_empty(),
            "cleavage should leave a side fragment"
        );
        assert!(!outcome.maybe().is_empty(), "maybe bags on the plan");
        assert!(
            outcome.allows(None, Some(outcome.maybe().entries[0].side.as_str())),
            "allows discarded side"
        );
        assert_eq!(outcome.plan[0].rule, "Dealkylation");
    }

    #[test]
    fn composed_ruleset_path_keeps_leaf_and_outer_namespace() {
        let set = RuleSet::compose(Some("Forest".into()), [hydroxylation(), o_dealkylation()]);
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &set, &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        let step = &hits[0].steps[0];
        assert_eq!(step.namespace(), vec!["Hydroxylation", "Forest"]);
        assert_eq!(step.leaf_rule(), Some("Hydroxylation"));
    }

    #[test]
    fn nested_ruleset_path_appends_each_set() {
        let inner = RuleSet::compose(Some("Inner".into()), [hydroxylation()]);
        let outer = RuleSet::compose(Some("Outer".into()), [inner]);
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &outer, &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        assert_eq!(
            hits[0].steps[0].namespace(),
            vec!["Hydroxylation", "Inner", "Outer"]
        );
    }

    #[test]
    fn already_at_target_yields_empty_plan() {
        let mut counters = PathCounters::default();
        let hits = find_path("CCO", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert_eq!(hits.len(), 1);
        assert!(hits[0].steps.is_empty());
        assert!(hits[0].plan.is_empty());
        assert!(hits[0].maybe().is_empty());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        assert_eq!(counters.nodes, 1);
        assert_eq!(counters.mol_edits, 0);
    }

    #[test]
    fn keep_predicate_skips_materialize() {
        let mut counters = PathCounters::default();
        let refuse_h2 = |c: &Candidate| c.pattern.name != "h2";
        let hits = find_path_with(
            "CC",
            "CCO",
            &hydroxylation(),
            &mut counters,
            FindPathConfig::default(),
            refuse_h2,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        // Ethane only matches h2; refusing it yields no path and no edits.
        assert!(hits.is_empty());
        assert_eq!(counters.mol_edits, 0);
    }

    #[test]
    fn default_ruleset_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("CC", "CCO", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Hydroxylation"));
        assert!(!hits[0].plan.is_empty());
    }

    #[test]
    fn default_ruleset_anisole_to_phenol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("COc1ccccc1", "Oc1ccccc1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("Oc1ccccc1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Dealkylation"));
        assert!(!hits[0].steps[0].sides.is_empty());
    }

    #[test]
    fn phase_one_ethene_to_epoxide() {
        let mut counters = PathCounters::default();
        let hits = find_path("C=C", "C1CO1", &crate::rules::phase_one(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("C1CO1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Epoxidation"));
    }

    #[test]
    fn default_ruleset_hydroquinone_to_quinone() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("O=C1C=CC(=O)C=C1").unwrap());
        let leaf = hits[0].steps[0].leaf_rule();
        assert!(
            leaf == Some("Dehydrogenation") || leaf == Some("QuinoneFormation"),
            "leaf={leaf:?} steps={:?}",
            hits[0].steps
        );
        assert!(
            !hits[0].plan.is_empty(),
            "expected elementary plan, got {:?}",
            hits[0].plan
        );
    }

    #[test]
    fn benzene_to_quinone_plan_precedes_both_oh_before_dh() {
        // Python test_find_path_uses_atom_diff_filters: precedes (0,2),(1,2).
        let mut counters = PathCounters::default();
        let hits = find_path_default("c1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let names: Vec<_> = hits[0].plan.iter().map(|s| s.rule.as_str()).collect();
        assert_eq!(
            names.iter().filter(|&&n| n == "Hydroxylation").count(),
            2,
            "plan={names:?}"
        );
        assert_eq!(
            names.iter().filter(|&&n| n == "Dehydrogenation").count(),
            1,
            "plan={names:?}"
        );
        assert!(!names.contains(&"QuinoneFormation"));
        let hydroxyl: Vec<_> = names
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n == "Hydroxylation").then_some(i))
            .collect();
        let dh = names.iter().position(|&n| n == "Dehydrogenation").unwrap();
        let edges: HashSet<_> = hits[0].plan.precedes().iter().copied().collect();
        assert!(!edges.contains(&(hydroxyl[0], hydroxyl[1])));
        assert!(!edges.contains(&(hydroxyl[1], hydroxyl[0])));
        assert!(edges.contains(&(hydroxyl[0], dh)));
        assert!(edges.contains(&(hydroxyl[1], dh)));
    }

    #[test]
    fn phenol_to_quinone_plan_oh_precedes_dh() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("Oc1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let names: Vec<_> = hits[0].plan.iter().map(|s| s.rule.as_str()).collect();
        assert!(!names.contains(&"QuinoneFormation"));
        assert_eq!(names.iter().filter(|&&n| n == "Hydroxylation").count(), 1);
        assert_eq!(names.iter().filter(|&&n| n == "Dehydrogenation").count(), 1);
        let oh = names.iter().position(|&n| n == "Hydroxylation").unwrap();
        let dh = names.iter().position(|&n| n == "Dehydrogenation").unwrap();
        let edges: HashSet<_> = hits[0].plan.precedes().iter().copied().collect();
        assert!(edges.contains(&(oh, dh)));
    }

    #[test]
    fn path_hit_plan_replays_to_hit() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("c1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert!(
            hits[0].plan.reaches("c1ccccc1", &hits[0].smiles).unwrap(),
            "plan={:?} smiles={}",
            hits[0].plan,
            hits[0].smiles
        );
    }

    #[test]
    fn ethane_plan_replays_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        assert!(hits[0].plan.reaches("CC", "CCO").unwrap());
    }

    #[test]
    fn atom_diff_gates_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_diff("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        // Diff gating should not inflate edits beyond the one helpful site.
        assert_eq!(counters.mol_edits, 1);
    }

    #[test]
    fn atom_diff_default_ruleset_hydroquinone() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "Oc1ccc(O)cc1",
            "O=C1C=CC(=O)C=C1",
            &default_ruleset(),
            &mut counters,
            FindPathConfig {
                use_atom_diff: true,
                ..FindPathConfig::default()
            },
            accept_all_candidates,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("O=C1C=CC(=O)C=C1").unwrap());
    }

    #[test]
    fn walk_preserves_tags_through_hydroxylation() {
        let parent = ForestMol::parse("CC").unwrap();
        let t0 = parent.tag_of(0).expect("stamped");
        let t1 = parent.tag_of(1).expect("stamped");
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        // Product mol is not on the outcome; check adopt via a fresh emit.
        let set = hydroxylation();
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let pieces = cands[0].materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        assert!(child.shares_tag_gen(&parent));
        assert_eq!(child.tag_of(child.index_of(t0).unwrap()), Some(t0));
        assert_eq!(child.tag_of(child.index_of(t1).unwrap()), Some(t1));
        assert_eq!(child.mol().atom_count(), 3);
    }

    #[test]
    fn cleavage_side_and_maybe_allows_on_deps() {
        use crate::canonical_plan::{CleavageSide, Maybe};
        use std::collections::BTreeSet;

        let side = CleavageSide::new([1, 2], canon_of("CCO").unwrap(), []);
        let empty = Maybe::default();
        assert!(empty.is_empty());

        let filled = Maybe::new([side]);
        assert!(!filled.is_empty());
        let demethyl: BTreeSet<_> = [0usize, 1].into_iter().collect();
        let formation: BTreeSet<_> = [1usize, 2].into_iter().collect();
        assert!(filled.allows(Some(&demethyl), None));
        assert!(!filled.allows(Some(&formation), None));
        assert!(filled.allows(None, Some("CCO")));
        let far: BTreeSet<_> = [9usize, 10].into_iter().collect();
        assert!(!filled.allows(Some(&far), None));

        let plan = Deps::bind([crate::canonical_plan::Step::new(
            "Dealkylation",
            [crate::canonical_plan::PlanAtom::index(1)],
        )])
        .with_maybe(filled);
        assert!(plan.allows(Some(&demethyl), None));
        assert!(!plan.allows(Some(&formation), None));
    }

    #[test]
    fn hydrolysis_emits_maybe_bag_on_plan() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            &crate::rules::hydrolysis(),
            &mut counters,
            FindPathConfig {
                max_paths: 3,
                max_nodes: 40,
                use_atom_diff: false,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let with_bag: Vec<_> = hits.iter().filter(|h| !h.maybe().is_empty()).collect();
        assert!(!with_bag.is_empty(), "expected CleavageSide on hydrolysis");
        let outcome = with_bag[0];
        assert_eq!(
            outcome
                .plan
                .iter()
                .map(|s| s.rule.as_str())
                .collect::<Vec<_>>(),
            vec!["Hydrolysis"]
        );
        let want_side = canon_of("CC(C)(C)O").unwrap();
        assert!(
            outcome.maybe().entries.iter().any(|e| e.side == want_side),
            "sides={:?}",
            outcome.maybe().sides()
        );
        assert!(outcome.allows(None, Some(&want_side)));
    }

    #[test]
    fn dealkylation_emits_maybe_and_multipath() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::dealkylation(),
            &mut counters,
            FindPathConfig {
                max_paths: 4,
                max_nodes: 80,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert!(
            hits.iter().any(|h| !h.maybe().is_empty()),
            "expected maybe bag"
        );
        let parent_ha = ForestMol::parse("CN(C)Cc1ccccc1")
            .unwrap()
            .heavy_atom_count();
        for h in &hits {
            for e in &h.maybe().entries {
                let side_ha = ForestMol::parse(&e.side).unwrap().heavy_atom_count();
                assert!(side_ha < parent_ha, "side={} ha={side_ha}", e.side);
            }
        }
    }

    #[test]
    fn multipath_emits_up_to_max_paths() {
        let mut counters = PathCounters::default();
        let one = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 1,
                max_nodes: 100,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert_eq!(one.len(), 1);

        let mut counters = PathCounters::default();
        let many = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 4,
                max_nodes: 200,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(
            many.len() > 1 && many.len() <= 4,
            "max_paths=1 → {} hits; max_paths=4 → {} (want >1)",
            one.len(),
            many.len()
        );
        assert!(many.iter().all(|h| h.smiles == one[0].smiles));
        // HEURISTICS: a reordering of the same Deps is not a new path.
        // Stronger observational check: multipath hits share no total orders.
        for (i, a) in many.iter().enumerate() {
            for b in many.iter().skip(i + 1) {
                assert!(
                    !a.plan.same_linearizations(&b.plan),
                    "hit{i} same_linearizations as later hit"
                );
                assert_eq!(
                    a.plan.linearization_overlap(&b.plan),
                    0,
                    "hit{i} linearization_overlap={} with later hit",
                    a.plan.linearization_overlap(&b.plan)
                );
            }
        }
    }

    #[test]
    fn plan_drop_counters_signal_only_no_contained_prune() {
        // Duplicate → count + drop. Contained-extension → count only, keep yield.
        use crate::canonical_plan::{PlanAtom, Step};
        let short = as_deps([Step::new("Dealkylation", [PlanAtom::index(0)])]);
        let longer = as_deps([
            Step::new("Dealkylation", [PlanAtom::index(0)]),
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
        ]);
        let twin = short.clone();
        let found = vec![PathOutcome {
            steps: vec![],
            plan: short,
            smiles: "C".into(),
        }];
        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &twin);
        assert_eq!(counters.dropped_duplicate_plan, 1);
        assert_eq!(counters.signal_contained_plan, 0);
        assert!(plan_already_yielded(&found, &twin));

        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &longer);
        assert_eq!(counters.dropped_duplicate_plan, 0);
        assert_eq!(counters.signal_contained_plan, 1);
        assert!(!plan_already_yielded(&found, &longer));
        assert_eq!(counters.plan_drops(), 1);
    }

    #[test]
    fn tba_cleavage_maybe_allows_overlapping_site() {
        // Terbinafine → TBA: formation Dealkylation leaves naphthyl-amine bag.
        const TERB: &str = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12";
        const TBA: &str = r"C(#C/C=C/C=O)C(C)(C)C";
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            TERB,
            TBA,
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 3,
                max_nodes: 80,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let direct: Vec<_> = hits
            .iter()
            .filter(|h| h.plan.len() == 1 && h.plan[0].rule == "Dealkylation")
            .collect();
        assert!(!direct.is_empty(), "expected Dealkylation → TBA");
        let outcome = direct[0];
        assert!(!outcome.maybe().is_empty());
        assert!(
            outcome
                .maybe()
                .entries
                .iter()
                .any(|e| e.side.contains("cccc")),
            "sides={:?}",
            outcome.maybe().sides()
        );
        let formation = outcome.maybe().entries[0].site.clone();
        assert!(!outcome.allows(Some(&formation), None));
        // Overlapping different site (share one atom) passes Maybe.
        let mut demethyl = BTreeSet::new();
        demethyl.insert(*formation.iter().next().unwrap());
        demethyl.insert(formation.iter().next().unwrap() + 1000);
        assert!(outcome.allows(Some(&demethyl), None));
        assert!(outcome.allows(None, Some(outcome.maybe().entries[0].side.as_str())));
    }
}
