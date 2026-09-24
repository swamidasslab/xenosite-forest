//! Cleavage product graph: Or of arms that share the same fragment multiset.
//!
//! Both sides of a bifurcation are first-class products (nodes). Arms that
//! produce the same sorted fragment CSMI multiset fold into one [`CleavageOr`].
//! Choosing an arm **and** a continuation fragment builds the linked
//! [`Maybe`] from the other side(s).
//!
//! MCS / [`crate::atom_diff`] gates which fragments stay expandable — the graph
//! does not bypass diff.
//!
//! This is **not** archive And/Or plan trees: Deps stay the plan language.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::ForestError;
use crate::atom_diff::{AtomDiff, atom_diff, candidate_could_help, pair_could_help};
use crate::canonical_plan::{CleavageSide, Maybe};
use crate::forest_mol::ForestMol;
use crate::mol::{Molecule, canon_of};
use crate::pair_edit::PairCandidate;
use crate::ruleset::RuleSet;

/// One way to produce a fragment multiset: rule/site on the parent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageArm {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    pub site_atoms: Vec<usize>,
    /// Both (all) bifurcation fragments, sorted CSMIs.
    pub products: Vec<String>,
}

impl CleavageArm {
    /// Linked Maybe when continuing the walk on `continue_csmi`.
    ///
    /// The other product CSMIs become [`CleavageSide`] entries for this arm.
    pub fn maybe_for(&self, continue_csmi: &str) -> Option<Maybe> {
        if !self.products.iter().any(|p| p == continue_csmi) {
            return None;
        }
        let entries: Vec<_> = self
            .products
            .iter()
            .filter(|p| p.as_str() != continue_csmi)
            .map(|side| {
                CleavageSide::new(
                    self.site_atoms.iter().copied(),
                    side.clone(),
                    std::iter::empty::<BTreeSet<usize>>(),
                )
            })
            .collect();
        Some(Maybe::new(entries))
    }
}

/// Or over arms that share the same sorted fragment multiset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageOr {
    /// Canonical key: sorted fragment CSMIs (both sides kept).
    pub fragments: Vec<String>,
    pub arms: Vec<CleavageArm>,
}

impl CleavageOr {
    pub fn n_arms(&self) -> usize {
        self.arms.len()
    }

    /// Choose an arm and which fragment to continue on; Maybe = the other side(s).
    pub fn choose(&self, arm_index: usize, continue_csmi: &str) -> Option<(&CleavageArm, Maybe)> {
        let arm = self.arms.get(arm_index)?;
        let maybe = arm.maybe_for(continue_csmi)?;
        Some((arm, maybe))
    }
}

/// One-hop cleavage layer: folded Or products from a single mol.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CleavageLayer {
    pub products: Vec<CleavageOr>,
}

impl CleavageLayer {
    pub fn n_products(&self) -> usize {
        self.products.len()
    }

    pub fn n_arms(&self) -> usize {
        self.products.iter().map(|p| p.arms.len()).sum()
    }

    pub fn max_or_fanin(&self) -> usize {
        self.products
            .iter()
            .map(|p| p.arms.len())
            .max()
            .unwrap_or(0)
    }

    /// Distinct fragment CSMIs appearing on either side of any Or.
    pub fn distinct_fragments(&self) -> BTreeSet<&str> {
        self.products
            .iter()
            .flat_map(|o| o.fragments.iter().map(|s| s.as_str()))
            .collect()
    }
}

/// Node in the multi-hop cleavage product graph (one fragment CSMI).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageNode {
    pub csmi: String,
    /// Incoming Or folds that produced this fragment (among others).
    /// `(parent CSMI, Or, this fragment was one of Or.fragments)`.
    pub via: Vec<(String, CleavageOr)>,
}

/// Multi-hop cleavage product graph (BFS of cleaving edits only).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CleavageGraph {
    pub nodes: Vec<CleavageNode>,
}

impl CleavageGraph {
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn n_ors(&self) -> usize {
        self.nodes.iter().map(|n| n.via.len()).sum()
    }

    pub fn n_arms(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|n| n.via.iter())
            .map(|(_, o)| o.arms.len())
            .sum()
    }

    pub fn max_or_fanin(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|n| n.via.iter())
            .map(|(_, o)| o.arms.len())
            .max()
            .unwrap_or(0)
    }
}

/// Options for building a cleavage layer / graph.
#[derive(Clone, Debug)]
pub struct CleavageGraphConfig {
    /// When set, MCS-gate fragments toward this target.
    pub target: Option<String>,
    /// Max distinct fragment CSMIs in a BFS graph (including the root).
    pub max_nodes: usize,
    /// Max BFS depth (root = 0).
    pub max_depth: usize,
}

impl Default for CleavageGraphConfig {
    fn default() -> Self {
        Self {
            target: None,
            max_nodes: 64,
            max_depth: 6,
        }
    }
}

fn sorted_fragments(products: &[String]) -> Vec<String> {
    let mut v = products.to_vec();
    v.sort();
    v
}

fn fold_layer(arms: Vec<CleavageArm>) -> CleavageLayer {
    let mut buckets: BTreeMap<Vec<String>, Vec<CleavageArm>> = BTreeMap::new();
    for arm in arms {
        let key = arm.products.clone();
        buckets.entry(key).or_default().push(arm);
    }
    let products = buckets
        .into_iter()
        .map(|(fragments, arms)| CleavageOr { fragments, arms })
        .collect();
    CleavageLayer { products }
}

/// True when this fragment is allowed to expand toward the target.
fn mcs_allows_fragment(
    parent_diff: Option<&AtomDiff>,
    fragment_csmi: &str,
    target: &Molecule,
) -> bool {
    let Ok(frag) = ForestMol::parse(fragment_csmi) else {
        return false;
    };
    let child_diff = atom_diff(frag.mol(), target);
    match parent_diff {
        Some(parent) => child_diff.cost() <= parent.cost(),
        None => true,
    }
}

/// A split is kept if it has ≥2 products and at least one fragment MCS-gates
/// (when a target is set). Both sides remain in the Or record either way.
fn split_usable(
    products: &[String],
    parent_diff: Option<&AtomDiff>,
    target: Option<&Molecule>,
) -> bool {
    if products.len() < 2 {
        return false;
    }
    match (parent_diff, target) {
        (Some(diff), Some(t)) => products
            .iter()
            .any(|p| mcs_allows_fragment(Some(diff), p, t)),
        _ => true,
    }
}

/// One-hop cleavage layer from `mol`, folding arms by sorted fragment multiset.
///
/// Both sides of every bifurcation are retained on [`CleavageArm::products`] /
/// [`CleavageOr::fragments`]. When `config.target` is set, a split is kept only
/// if at least one fragment does not worsen MCS cost vs the parent.
pub fn cleavage_layer(
    mol: &Molecule,
    ruleset: &RuleSet,
    config: &CleavageGraphConfig,
) -> Result<CleavageLayer, ForestError> {
    let target_mol = match &config.target {
        Some(t) => Some(crate::mol::parse_mol(t)?),
        None => None,
    };
    let parent_diff = target_mol.as_ref().map(|t| atom_diff(mol, t));

    let mut raw: Vec<CleavageArm> = Vec::new();

    for c in ruleset.candidates(mol)? {
        if !c.pattern.effect.cleaves {
            continue;
        }
        if let Some(diff) = &parent_diff {
            if !candidate_could_help(&c, diff) {
                continue;
            }
        }
        let Some(emission) = c.emit(mol)? else {
            continue;
        };
        if !emission.cleaves {
            continue;
        }
        let products = sorted_fragments(&emission.products);
        if !split_usable(&products, parent_diff.as_ref(), target_mol.as_ref()) {
            continue;
        }
        let rule = emission
            .leaf_rule()
            .unwrap_or(emission.pattern_name.as_str())
            .to_string();
        raw.push(CleavageArm {
            rule,
            pattern_name: emission.pattern_name,
            site: emission.site,
            site_orbit: emission.site_orbit,
            site_atoms: emission.site_atoms,
            products,
        });
    }

    for pair in ruleset.pair_candidates(mol)? {
        if !pair.effect.cleaves {
            continue;
        }
        if let (Some(diff), Some(t)) = (&parent_diff, &target_mol) {
            if !pair_could_help(&pair, diff, mol, t) {
                continue;
            }
        }
        push_pair_arm(
            mol,
            &pair,
            parent_diff.as_ref(),
            target_mol.as_ref(),
            &mut raw,
        )?;
    }

    Ok(fold_layer(raw))
}

fn push_pair_arm(
    mol: &Molecule,
    pair: &PairCandidate,
    parent_diff: Option<&AtomDiff>,
    target_mol: Option<&Molecule>,
    raw: &mut Vec<CleavageArm>,
) -> Result<(), ForestError> {
    let Some(emission) = pair.emit(mol)? else {
        return Ok(());
    };
    let products = sorted_fragments(&emission.products);
    if !split_usable(&products, parent_diff, target_mol) {
        return Ok(());
    }
    raw.push(CleavageArm {
        rule: pair
            .pattern_name
            .split('+')
            .next()
            .unwrap_or("Pair")
            .to_string(),
        pattern_name: emission.pattern_name,
        site: emission.site,
        site_orbit: vec![emission.site],
        site_atoms: pair.plan_site_atoms(),
        products,
    });
    Ok(())
}

/// BFS cleavage product graph from `start`.
///
/// Every fragment CSMI from each Or is a node (both sides). Non-cleaving
/// chemistry is not applied here — callers run that across nodes.
pub fn cleavage_product_graph(
    start: &str,
    ruleset: &RuleSet,
    config: &CleavageGraphConfig,
) -> Result<CleavageGraph, ForestError> {
    let root_csmi = canon_of(start)?;
    let mut nodes = vec![CleavageNode {
        csmi: root_csmi.clone(),
        via: Vec::new(),
    }];
    let mut index: BTreeMap<String, usize> = BTreeMap::new();
    index.insert(root_csmi.clone(), 0);

    let mut queue: VecDeque<(usize, usize)> = VecDeque::new();
    queue.push_back((0, 0));

    let target_mol = match &config.target {
        Some(t) => Some(crate::mol::parse_mol(t)?),
        None => None,
    };

    while let Some((ni, depth)) = queue.pop_front() {
        if depth >= config.max_depth || nodes.len() >= config.max_nodes {
            continue;
        }
        let parent_csmi = nodes[ni].csmi.clone();
        let parent_mol = crate::mol::parse_mol(&parent_csmi)?;
        let parent_diff = target_mol.as_ref().map(|t| atom_diff(&parent_mol, t));
        let layer = cleavage_layer(&parent_mol, ruleset, config)?;

        for or in layer.products {
            for frag in &or.fragments {
                if let (Some(diff), Some(t)) = (&parent_diff, &target_mol) {
                    if !mcs_allows_fragment(Some(diff), frag, t) {
                        continue;
                    }
                }
                let child_i = match index.entry(frag.clone()) {
                    std::collections::btree_map::Entry::Occupied(o) => *o.get(),
                    std::collections::btree_map::Entry::Vacant(v) => {
                        if nodes.len() >= config.max_nodes {
                            break;
                        }
                        let i = nodes.len();
                        nodes.push(CleavageNode {
                            csmi: frag.clone(),
                            via: Vec::new(),
                        });
                        v.insert(i);
                        queue.push_back((i, depth + 1));
                        i
                    }
                };
                // Avoid duplicate identical Or records on the same via list.
                let already = nodes[child_i]
                    .via
                    .iter()
                    .any(|(p, o)| p == &parent_csmi && o.fragments == or.fragments);
                if !already {
                    nodes[child_i].via.push((parent_csmi.clone(), or.clone()));
                }
            }
        }
    }

    Ok(CleavageGraph { nodes })
}

/// Summary counts for benches / tests.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageGraphStats {
    pub nodes: usize,
    pub ors: usize,
    pub arms: usize,
    pub max_or_fanin: usize,
    pub one_hop_products: usize,
    pub one_hop_arms: usize,
    pub one_hop_fragments: usize,
}

pub fn cleavage_graph_stats(
    start: &str,
    target: Option<&str>,
    ruleset: &RuleSet,
) -> Result<(CleavageGraphStats, CleavageGraph), ForestError> {
    let config = CleavageGraphConfig {
        target: target.map(|t| t.to_string()),
        ..CleavageGraphConfig::default()
    };
    let start_mol = crate::mol::parse_mol(start)?;
    let layer = cleavage_layer(&start_mol, ruleset, &config)?;
    let graph = cleavage_product_graph(start, ruleset, &config)?;
    Ok((
        CleavageGraphStats {
            nodes: graph.n_nodes(),
            ors: graph.n_ors(),
            arms: graph.n_arms(),
            max_or_fanin: graph.max_or_fanin(),
            one_hop_products: layer.n_products(),
            one_hop_arms: layer.n_arms(),
            one_hop_fragments: layer.distinct_fragments().len(),
        },
        graph,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::{leaf_rule, phase_one};

    #[test]
    fn anisole_keeps_both_cleavage_fragments() {
        let mol = crate::mol::parse_mol("COc1ccccc1").unwrap();
        let dealk = leaf_rule("Dealkylation").expect("Dealkylation");
        let layer = cleavage_layer(
            &mol,
            &dealk,
            &CleavageGraphConfig {
                target: Some("Oc1ccccc1".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(layer.n_products() >= 1, "products={:?}", layer.products);
        let or = layer
            .products
            .iter()
            .find(|o| o.fragments.len() >= 2)
            .expect("expected a bifurcation Or");
        assert_eq!(or.fragments.len(), 2, "fragments={:?}", or.fragments);
        // Both sides retained — not a single kept/side pick.
        assert_eq!(or.arms[0].products, or.fragments);
    }

    #[test]
    fn choose_continuation_puts_other_side_in_maybe() {
        let mol = crate::mol::parse_mol("COc1ccc(OC)cc1").unwrap();
        let layer = cleavage_layer(
            &mol,
            &phase_one(),
            &CleavageGraphConfig {
                target: Some("Oc1ccc(O)cc1".into()),
                ..Default::default()
            },
        )
        .unwrap();
        let or = layer
            .products
            .iter()
            .find(|o| o.fragments.len() >= 2)
            .unwrap();
        let continue_on = or.fragments[0].as_str();
        let other = or.fragments[1].as_str();
        let (arm, maybe) = or.choose(0, continue_on).unwrap();
        assert!(arm.products.contains(&continue_on.to_string()));
        assert_eq!(maybe.sides(), vec![other]);
    }

    #[test]
    fn tetramethoxy_biphenyl_or_fanin_and_both_sides() {
        let start = "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC";
        let target = "Oc1ccc(-c2ccc(O)c(O)c2)cc1O";
        let (stats, graph) = cleavage_graph_stats(start, Some(target), &phase_one()).unwrap();
        assert!(
            stats.one_hop_arms >= stats.one_hop_products,
            "stats={stats:?}"
        );
        assert!(
            stats.max_or_fanin >= 2 || stats.one_hop_arms >= 2,
            "stats={stats:?}"
        );
        assert!(
            stats.one_hop_fragments >= 2,
            "both sides should appear, stats={stats:?}"
        );
        assert!(stats.nodes >= 2, "graph should keep multiple fragments");
        assert!(stats.nodes < 40, "nodes={}", stats.nodes);
        // At least one via Or on a non-root node carries ≥2 fragments.
        let both_sides = graph
            .nodes
            .iter()
            .any(|n| n.via.iter().any(|(_, o)| o.fragments.len() >= 2));
        assert!(both_sides, "no Or retained both cleavage fragments");
    }

    #[test]
    fn skip_non_cleaving_patterns() {
        let mol = crate::mol::parse_mol("c1ccccc1").unwrap();
        let layer = cleavage_layer(
            &mol,
            &crate::hydroxylation::hydroxylation(),
            &CleavageGraphConfig::default(),
        )
        .unwrap();
        assert!(layer.products.is_empty());
    }
}
