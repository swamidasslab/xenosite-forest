//! First-run plan-guided search over a [`crate::ruleset::RuleSet`].
//!
//! Expands via [`RuleSet::candidates`]: site–pattern–parent triples. A search
//! reads [`PatternInfo`] / effect fields to decide, then
//! [`crate::candidate::Candidate::materialize`] only for survivors — no filter
//! closures required. Nested sets stay namespaces on each step's leaf-first
//! `rule_path`.
//!
//! This is not full Python `find_path`: no atom-diff filters, no `CanonicalStep` /
//! `Deps` plan, no stale-priority heap rescore. Closer is a provisional heavy-atom
//! distance (refuse a child that moves away from the target size). Cleavage keeps
//! the fragment that hits the target CSMI, else the one nearest in heavy-atom count.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use crate::ForestError;
use crate::candidate::Candidate;
use crate::mol::{canon_of, canon_smiles, parse_mol};
use crate::pattern::{Emission, PatternInfo, SiteInfo};
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

    fn from_emission(emission: &Emission, product: String, sides: Vec<String>) -> Self {
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
    pub smiles: String,
}

#[derive(Clone, Debug)]
struct Walk {
    smiles: String,
    heavy: usize,
    steps: Vec<PathStep>,
}

/// Heap entry: hits first, then FIFO (`seq`). Lower priority value pops first.
#[derive(Clone, Debug)]
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

fn heavy_atoms(smiles: &str) -> Result<usize, ForestError> {
    Ok(parse_mol(smiles)?.atom_count())
}

fn ha_distance(ha: usize, target_ha: usize) -> usize {
    ha.abs_diff(target_ha)
}

/// Keep the fragment closest to ``target``. Prefer exact CSMI hit.
fn keep_fragment(
    products: &[String],
    target_csmi: &str,
    target_ha: usize,
) -> Result<Option<(String, Vec<String>)>, ForestError> {
    if products.is_empty() {
        return Ok(None);
    }
    let mut best: Option<(String, usize)> = None;
    for raw in products {
        let csmi = canon_of(raw)?;
        let cost = if csmi == target_csmi {
            0
        } else {
            1 + ha_distance(heavy_atoms(&csmi)?, target_ha)
        };
        match &best {
            None => best = Some((csmi, cost)),
            Some((_, best_cost)) if cost < *best_cost => best = Some((csmi, cost)),
            Some((kept, best_cost)) if cost == *best_cost && csmi < *kept => {
                best = Some((csmi, cost));
            }
            _ => {}
        }
    }
    let Some((kept, _)) = best else {
        return Ok(None);
    };
    let mut sides = Vec::new();
    for raw in products {
        let csmi = canon_of(raw)?;
        if csmi != kept {
            sides.push(csmi);
        }
    }
    Ok(Some((kept, sides)))
}

fn closer(parent_ha: usize, child_ha: usize, target_ha: usize, target_hit: bool) -> bool {
    target_hit || ha_distance(child_ha, target_ha) < ha_distance(parent_ha, target_ha)
}

/// Search bounds. Defaults match Python `find_path` knobs.
#[derive(Clone, Copy, Debug)]
pub struct FindPathConfig {
    pub max_paths: usize,
    pub max_nodes: usize,
}

impl Default for FindPathConfig {
    fn default() -> Self {
        Self {
            max_paths: 1,
            max_nodes: 800,
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
    } = config;
    let start = parse_mol(reactant)?;
    let start_csmi = canon_smiles(&start);
    let target_csmi = canon_of(target)?;
    let target_ha = heavy_atoms(&target_csmi)?;
    let start_ha = start.atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi == target_csmi,
        seq,
        walk: Walk {
            smiles: start_csmi.clone(),
            heavy: start_ha,
            steps: Vec::new(),
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi);
    let mut found = Vec::new();

    while let Some(item) = heap.pop() {
        if found.len() >= max_paths || counters.nodes >= max_nodes {
            break;
        }
        counters.nodes += 1;
        let walk = item.walk;
        if walk.smiles == target_csmi {
            found.push(PathOutcome {
                steps: walk.steps,
                smiles: walk.smiles,
            });
            continue;
        }

        let mol = parse_mol(&walk.smiles)?;
        counters.expansions += 1;
        let emissions = expand(ruleset, &mol, counters, &keep)?;
        let mut hits_from_here = 0usize;

        for emission in emissions {
            let Some((kept, sides)) = keep_fragment(&emission.products, &target_csmi, target_ha)?
            else {
                continue;
            };
            let child_ha = heavy_atoms(&kept)?;
            let target_hit = kept == target_csmi;
            if !closer(walk.heavy, child_ha, target_ha, target_hit) {
                continue;
            }
            if seen.contains(&kept) && !target_hit {
                continue;
            }
            seen.insert(kept.clone());

            let mut steps = walk.steps.clone();
            steps.push(PathStep::from_emission(&emission, kept.clone(), sides));
            heap.push(HeapItem {
                target_hit,
                seq,
                walk: Walk {
                    smiles: kept,
                    heavy: child_ha,
                    steps,
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
    mol: &crate::Molecule,
    counters: &mut PathCounters,
    keep: &K,
) -> Result<Vec<Emission>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let mut emissions = Vec::new();
    for candidate in ruleset.candidates(mol)? {
        if !keep(&candidate) {
            continue;
        }
        counters.mol_edits += 1;
        if let Some(emission) = candidate.emit(mol)? {
            emissions.push(emission);
        }
    }
    emissions.extend(expand_pairs(ruleset, mol, counters, keep)?);
    Ok(emissions)
}

fn expand_pairs<K>(
    ruleset: &RuleSet,
    mol: &crate::Molecule,
    counters: &mut PathCounters,
    keep: &K,
) -> Result<Vec<Emission>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let mut out = Vec::new();
    for member in ruleset.members() {
        if let crate::ruleset::RuleMember::Set(child) = member {
            for mut emission in expand_pairs(child, mol, counters, keep)? {
                emission.rule_path.push(ruleset.name.clone());
                out.push(emission);
            }
        }
    }
    for pair in ruleset.pair_candidates_leaf(mol)? {
        if !keep_pair(&pair, keep) {
            continue;
        }
        counters.mol_edits += 1;
        if let Some(emission) = pair.emit(mol)? {
            out.push(Emission {
                site: emission.site,
                pattern_name: emission.pattern_name,
                rule_path: vec![ruleset.name.clone()],
                products: emission.products,
            });
        }
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
    } = config;
    let start = parse_mol(reactant)?;
    let start_csmi = canon_smiles(&start);
    let target_csmi = canon_of(target)?;
    let target_ha = heavy_atoms(&target_csmi)?;
    let start_ha = start.atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    heap.push(HeapItem {
        target_hit: start_csmi == target_csmi,
        seq,
        walk: Walk {
            smiles: start_csmi.clone(),
            heavy: start_ha,
            steps: Vec::new(),
        },
    });
    seq += 1;

    let mut seen = HashSet::new();
    seen.insert(start_csmi);
    let mut found = Vec::new();

    while let Some(item) = heap.pop() {
        if found.len() >= max_paths || counters.nodes >= max_nodes {
            break;
        }
        counters.nodes += 1;
        let walk = item.walk;
        if walk.smiles == target_csmi {
            found.push(PathOutcome {
                steps: walk.steps,
                smiles: walk.smiles,
            });
            continue;
        }

        let mol = parse_mol(&walk.smiles)?;
        counters.expansions += 1;
        let emissions = ruleset.metabolize(&mol, &filter_rules, &filter_sites, true)?;
        let mut hits_from_here = 0usize;

        for emission in emissions {
            counters.mol_edits += 1;
            let Some((kept, sides)) = keep_fragment(&emission.products, &target_csmi, target_ha)?
            else {
                continue;
            };
            let child_ha = heavy_atoms(&kept)?;
            let target_hit = kept == target_csmi;
            if !closer(walk.heavy, child_ha, target_ha, target_hit) {
                continue;
            }
            if seen.contains(&kept) && !target_hit {
                continue;
            }
            seen.insert(kept.clone());

            let mut steps = walk.steps.clone();
            steps.push(PathStep::from_emission(&emission, kept.clone(), sides));
            heap.push(HeapItem {
                target_hit,
                seq,
                walk: Walk {
                    smiles: kept,
                    heavy: child_ha,
                    steps,
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
    }

    #[test]
    fn candidates_filter_without_closures() {
        let mol = parse_mol("CCC").unwrap();
        let set = hydroxylation();
        let kept: Vec<_> = set
            .candidates(&mol)
            .unwrap()
            .into_iter()
            .filter(|c| c.pattern.name == "h")
            .collect();
        assert!(kept.is_empty(), "propane has no #6h1");
        let h2: Vec<_> = set
            .candidates(&mol)
            .unwrap()
            .into_iter()
            .filter(|c| c.pattern.name == "h2")
            .collect();
        assert_eq!(h2.len(), 2);
        let products = h2[0].materialize(&mol).unwrap();
        assert!(!products.is_empty());
    }
}
