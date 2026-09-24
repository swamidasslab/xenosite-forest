//! First-run plan-guided search over a [`crate::ruleset::RuleSet`].
//!
//! Walks carry tagged [`crate::forest_mol::ForestMol`] (structure/`csmi` cache,
//! atom tags through edits). Expands via [`RuleSet::candidates`]: site–pattern–
//! parent triples. A search reads [`PatternInfo`] / effect fields to decide,
//! then materializes only for survivors. Nested sets stay namespaces on each
//! step's leaf-first `rule_path`.
//!
//! Outcomes carry elementary [`crate::canonical_plan::CanonicalStep`] plans
//! (identity, or quinone-shaped prep-then-DH from [`PlanKind`] on the leaf).
//! Closer uses atom-diff cost; child diffs prefer tag-lifted parent MCS when
//! tags allow ([`crate::atom_diff::try_atom_diff_for_child`]) — default lazy
//! and eager both lift at enqueue; full MCS only on lift miss.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};
use std::rc::Rc;

use crate::ForestError;
use crate::candidate::Candidate;
use crate::canonical_plan::{CanonicalStep, steps_for_kind};
use crate::forest_mol::ForestMol;
use crate::mol::{canon_of, parse_mol};
use crate::pattern::{PatternInfo, SiteInfo};
use crate::rules::default_ruleset;
use crate::ruleset::RuleSet;

/// Billed work for one search (Python `PathCounters` subset).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PathCounters {
    pub nodes: usize,
    pub mol_edits: usize,
    pub expansions: usize,
}

impl PathCounters {
    pub fn billed(&self) -> usize {
        self.mol_edits + self.nodes
    }
}

/// One accepted edit on a walk: RuleSet namespace path + pattern + kept product.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathStep {
    pub rule_path: Vec<Option<String>>,
    pub pattern_name: String,
    pub site: usize,
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

    fn from_emission(emission: &ForestEmission, product: String, sides: Vec<String>) -> Self {
        Self {
            rule_path: emission.rule_path.clone(),
            pattern_name: emission.pattern_name.clone(),
            site: emission.site,
            product,
            sides,
        }
    }
}

/// One reactant→target hit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathOutcome {
    pub steps: Vec<PathStep>,
    /// Elementary plan (Python `CanonicalStep` / `Deps` shape).
    pub plan: Vec<CanonicalStep>,
    pub smiles: String,
}

#[derive(Clone)]
struct Walk {
    mol: ForestMol,
    steps: Vec<PathStep>,
    plan: Vec<CanonicalStep>,
    /// Parent's [`crate::atom_diff::AtomDiff::cost`] when this walk was
    /// enqueued. `None` = root (always expand).
    parent_cost: Option<usize>,
    /// Diff of **this** mol vs target. Filled at enqueue by tag-lift (+ OH
    /// extend) when possible; on pop, full MCS only if still `None`.
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

/// Tagged emission: ForestMol products + elementary plan.
#[derive(Clone)]
struct ForestEmission {
    site: usize,
    pattern_name: String,
    rule_path: Vec<Option<String>>,
    products: Vec<ForestMol>,
    plan: Vec<CanonicalStep>,
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

/// Yield walks that turn ``reactant`` into ``target``.
///
/// Discovers [`Candidate`]s, keeps all of them, materializes, plus ResonancePair
/// leaf emissions. Nested sets keep step namespaces on each `rule_path`.
pub fn find_path(
    reactant: &str,
    target: &str,
    ruleset: &RuleSet,
    counters: &mut PathCounters,
) -> Result<Vec<PathOutcome>, ForestError> {
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
pub fn find_path_default(
    reactant: &str,
    target: &str,
    counters: &mut PathCounters,
) -> Result<Vec<PathOutcome>, ForestError> {
    find_path(reactant, target, &default_ruleset(), counters)
}

/// [`find_path`] with atom-diff candidate gating (no filter closures).
pub fn find_path_diff(
    reactant: &str,
    target: &str,
    ruleset: &RuleSet,
    counters: &mut PathCounters,
) -> Result<Vec<PathOutcome>, ForestError> {
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
pub fn find_path_with<K>(
    reactant: &str,
    target: &str,
    ruleset: &RuleSet,
    counters: &mut PathCounters,
    config: FindPathConfig,
    keep: K,
) -> Result<Vec<PathOutcome>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let FindPathConfig {
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
    } = config;
    let start = ForestMol::parse(reactant)?;
    let start_csmi = start.csmi();
    let target_csmi = canon_of(target)?;
    let target_mol = parse_mol(&target_csmi)?;
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
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi.as_ref().to_string());
    let mut found = Vec::new();

    while let Some(item) = heap.pop() {
        if found.len() >= max_paths || counters.nodes >= max_nodes {
            break;
        }
        let walk = item.walk;
        let here = walk.mol.csmi();
        if here.as_ref() == target_csmi.as_str() {
            counters.nodes += 1;
            found.push(PathOutcome {
                steps: walk.steps,
                plan: walk.plan,
                smiles: here.as_ref().to_string(),
            });
            continue;
        }

        let diff = if use_atom_diff {
            // Prefer tag-lifted diff from enqueue; full MCS only when missing.
            let d = match walk.diff {
                Some(d) => d,
                None => crate::atom_diff::atom_diff(walk.mol.mol(), &target_mol),
            };
            // Lazy closer: verify cost against parent before expanding.
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
        counters.nodes += 1;
        counters.expansions += 1;
        let parent_cost = diff.as_ref().map(|d| d.cost());
        let emissions = expand(
            ruleset,
            &walk.mol,
            &target_mol,
            counters,
            &keep,
            diff.as_ref(),
        )?;
        let mut hits_from_here = 0usize;
        let parent_ha = walk.mol.heavy_atom_count();

        for emission in emissions {
            let Some((kept, sides)) = keep_fragment(&emission.products, &target_csmi, target_ha)
            else {
                continue;
            };
            let kept_csmi = kept.csmi().as_ref().to_string();
            let child_ha = kept.heavy_atom_count();
            let target_hit = kept_csmi == target_csmi;

            // Tag-lift (+ local add extend / remove shrink). `None` → full MCS
            // later (eager: now; lazy: on pop).
            let mut child_diff = if use_atom_diff {
                diff.as_ref().and_then(|parent_d| {
                    crate::atom_diff::try_atom_diff_for_child(
                        &walk.mol,
                        parent_d,
                        &kept,
                        &target_mol,
                    )
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
                            Some(crate::atom_diff::atom_diff(kept.mol(), &target_mol));
                    }
                    cost_closer(pc, child_diff.as_ref().unwrap().cost(), target_hit)
                } else {
                    true
                }
            } else {
                closer(parent_ha, child_ha, target_ha, target_hit)
            };
            if !allow {
                continue;
            }
            if seen.contains(&kept_csmi) && !target_hit {
                continue;
            }
            seen.insert(kept_csmi.clone());

            let mut steps = walk.steps.clone();
            steps.push(PathStep::from_emission(&emission, kept_csmi.clone(), sides));
            let mut plan = walk.plan.clone();
            plan.extend(emission.plan.iter().cloned());
            heap.push(HeapItem {
                target_hit,
                seq,
                walk: Walk {
                    mol: kept,
                    steps,
                    plan,
                    parent_cost,
                    diff: child_diff,
                },
            });
            seq += 1;
            if target_hit {
                hits_from_here += 1;
                if found.len() + hits_from_here >= max_paths {
                    break;
                }
            }
        }
    }

    Ok(found)
}

fn expand<K>(
    ruleset: &RuleSet,
    parent: &ForestMol,
    target: &crate::Molecule,
    counters: &mut PathCounters,
    keep: &K,
    diff: Option<&crate::atom_diff::AtomDiff>,
) -> Result<Vec<ForestEmission>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let mol = parent.mol();
    let mut candidates = ruleset.candidates(mol)?;
    if let Some(d) = diff {
        candidates.retain(|c| {
            keep(c) && crate::atom_diff::candidate_could_help_on(c, d, Some(mol), Some(target))
        });
        candidates.sort_by_key(|c| crate::atom_diff::candidate_order_key(c, d));
    } else {
        candidates.retain(|c| keep(c));
    }
    let mut emissions = Vec::new();
    for candidate in candidates {
        counters.mol_edits += 1;
        if let Some(emission) = emit_candidate(&candidate, parent)? {
            emissions.push(emission);
        }
    }
    emissions.extend(expand_pairs(ruleset, parent, target, counters, keep, diff)?);
    Ok(emissions)
}

fn emit_candidate(
    candidate: &Candidate,
    parent: &ForestMol,
) -> Result<Option<ForestEmission>, ForestError> {
    let pieces = candidate.materialize_mols(parent.mol())?;
    if pieces.is_empty() {
        return Ok(None);
    }
    let products: Vec<ForestMol> = pieces
        .into_iter()
        .map(|piece| parent.adopt_product(piece))
        .collect();
    Ok(Some(ForestEmission {
        site: candidate.site,
        pattern_name: candidate.pattern.name.clone(),
        rule_path: candidate.rule_path.clone(),
        products,
        plan: candidate.identity_plan(),
    }))
}

fn expand_pairs<K>(
    ruleset: &RuleSet,
    parent: &ForestMol,
    target: &crate::Molecule,
    counters: &mut PathCounters,
    keep: &K,
    diff: Option<&crate::atom_diff::AtomDiff>,
) -> Result<Vec<ForestEmission>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let mol = parent.mol();
    let mut out = Vec::new();
    for member in ruleset.members() {
        if let crate::ruleset::RuleMember::Set(child) = member {
            for mut emission in expand_pairs(child, parent, target, counters, keep, diff)? {
                emission.rule_path.push(ruleset.name.clone());
                out.push(emission);
            }
        }
    }
    let mut pairs = ruleset.pair_candidates_leaf(mol)?;
    if let Some(d) = diff {
        pairs.retain(|pair| {
            if !keep_pair(pair, keep) {
                return false;
            }
            crate::atom_diff::pair_could_help(pair, d, mol, target)
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
        pairs.retain(|p| keep_pair(p, keep));
    }
    for pair in pairs {
        // Python ResonancePair bumps mol_edits only after alternating paths exist
        // (emit yields products). Empty materializations are not billed.
        let pieces = pair.materialize_mols(mol)?;
        if pieces.is_empty() {
            continue;
        }
        counters.mol_edits += 1;
        let products: Vec<ForestMol> = pieces
            .into_iter()
            .map(|piece| parent.adopt_product(piece))
            .collect();
        let site_atoms = pair.plan_site_atoms();
        let ends = [&pair.left.effect, &pair.right.effect];
        let leaf = ruleset
            .name
            .as_deref()
            .unwrap_or(pair.pattern_name.as_str());
        let plan = steps_for_kind(ruleset.plan_kind, leaf, mol, &site_atoms, Some(&ends));
        out.push(ForestEmission {
            site: pair.site,
            pattern_name: pair.pattern_name.clone(),
            rule_path: vec![ruleset.name.clone()],
            products,
            plan,
        });
    }
    Ok(out)
}

fn keep_pair<K>(pair: &crate::pair_edit::PairCandidate, keep: &K) -> bool
where
    K: Fn(&Candidate) -> bool,
{
    let mut stand_in = Candidate {
        site: pair.site,
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
pub fn find_path_with_filters<R, S>(
    reactant: &str,
    target: &str,
    ruleset: &RuleSet,
    counters: &mut PathCounters,
    config: FindPathConfig,
    filter_rules: R,
    filter_sites: S,
) -> Result<Vec<PathOutcome>, ForestError>
where
    R: Fn(&crate::Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&crate::Molecule, usize, &SiteInfo) -> bool,
{
    let FindPathConfig {
        max_paths,
        max_nodes,
        use_atom_diff: _,
        lazy_closer: _,
    } = config;
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
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi.as_ref().to_string());
    let mut found = Vec::new();

    while let Some(item) = heap.pop() {
        if found.len() >= max_paths || counters.nodes >= max_nodes {
            break;
        }
        counters.nodes += 1;
        let walk = item.walk;
        let here = walk.mol.csmi();
        if here.as_ref() == target_csmi.as_str() {
            found.push(PathOutcome {
                steps: walk.steps,
                plan: walk.plan,
                smiles: here.as_ref().to_string(),
            });
            continue;
        }

        let mol = walk.mol.mol();
        counters.expansions += 1;
        let string_emissions = ruleset.metabolize(mol, &filter_rules, &filter_sites, true)?;
        let mut hits_from_here = 0usize;

        for emission in string_emissions {
            counters.mol_edits += 1;
            // Re-parse products into ForestMol (filter path; no tag continuity).
            let products: Result<Vec<_>, _> = emission
                .products
                .iter()
                .map(|s| ForestMol::parse(s))
                .collect();
            let products = products?;
            let Some((kept, sides)) = keep_fragment(&products, &target_csmi, target_ha) else {
                continue;
            };
            let kept_csmi = kept.csmi().as_ref().to_string();
            let child_ha = kept.heavy_atom_count();
            let target_hit = kept_csmi == target_csmi;
            if !closer(walk.mol.heavy_atom_count(), child_ha, target_ha, target_hit) {
                continue;
            }
            if seen.contains(&kept_csmi) && !target_hit {
                continue;
            }
            seen.insert(kept_csmi.clone());

            let mut steps = walk.steps.clone();
            steps.push(PathStep {
                rule_path: emission.rule_path.clone(),
                pattern_name: emission.pattern_name.clone(),
                site: emission.site,
                product: kept_csmi.clone(),
                sides,
            });
            let mut plan = walk.plan.clone();
            plan.extend(emission.plan.iter().cloned());
            heap.push(HeapItem {
                target_hit,
                seq,
                walk: Walk {
                    mol: kept,
                    steps,
                    plan,
                    parent_cost: None,
                    // Filter path re-parses; no tag continuity → no lift.
                    diff: None,
                },
            });
            seq += 1;
            if target_hit {
                hits_from_here += 1;
                if found.len() + hits_from_here >= max_paths {
                    break;
                }
            }
        }
    }

    Ok(found)
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
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters).unwrap();
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
        let hits = find_path("COc1ccccc1", "Oc1ccccc1", &o_dealkylation(), &mut counters).unwrap();
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
        assert_eq!(outcome.plan[0].rule, "Dealkylation");
    }

    #[test]
    fn composed_ruleset_path_keeps_leaf_and_outer_namespace() {
        let set = RuleSet::compose(Some("Forest".into()), [hydroxylation(), o_dealkylation()]);
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &set, &mut counters).unwrap();
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
        let hits = find_path("CC", "CCO", &outer, &mut counters).unwrap();
        assert!(!hits.is_empty());
        assert_eq!(
            hits[0].steps[0].namespace(),
            vec!["Hydroxylation", "Inner", "Outer"]
        );
    }

    #[test]
    fn already_at_target_yields_empty_plan() {
        let mut counters = PathCounters::default();
        let hits = find_path("CCO", "CCO", &hydroxylation(), &mut counters).unwrap();
        assert_eq!(hits.len(), 1);
        assert!(hits[0].steps.is_empty());
        assert!(hits[0].plan.is_empty());
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
        .unwrap();
        // Ethane only matches h2; refusing it yields no path and no edits.
        assert!(hits.is_empty());
        assert_eq!(counters.mol_edits, 0);
    }

    #[test]
    fn default_ruleset_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("CC", "CCO", &mut counters).unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Hydroxylation"));
        assert!(!hits[0].plan.is_empty());
    }

    #[test]
    fn default_ruleset_anisole_to_phenol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("COc1ccccc1", "Oc1ccccc1", &mut counters).unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("Oc1ccccc1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Dealkylation"));
        assert!(!hits[0].steps[0].sides.is_empty());
    }

    #[test]
    fn phase_one_ethene_to_epoxide() {
        let mut counters = PathCounters::default();
        let hits = find_path("C=C", "C1CO1", &crate::rules::phase_one(), &mut counters).unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("C1CO1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Epoxidation"));
    }

    #[test]
    fn default_ruleset_hydroquinone_to_quinone() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1", &mut counters).unwrap();
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
    fn atom_diff_gates_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_diff("CC", "CCO", &hydroxylation(), &mut counters).unwrap();
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
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters).unwrap();
        assert!(!hits.is_empty());
        // Product mol is not on the outcome; check adopt via a fresh emit.
        let set = hydroxylation();
        let cands = set.candidates(parent.mol()).unwrap();
        let pieces = cands[0].materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        assert!(child.shares_tag_gen(&parent));
        assert_eq!(child.tag_of(child.index_of(t0).unwrap()), Some(t0));
        assert_eq!(child.tag_of(child.index_of(t1).unwrap()), Some(t1));
        assert_eq!(child.mol().atom_count(), 3);
    }
}
