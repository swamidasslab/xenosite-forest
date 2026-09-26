//! Stream metabolites up to a depth (BFS or DFS) with path info.
//!
//! Sibling of Python `bfs` / `dfs` in `xenosite.forest.find_path`: frontier
//! expand via [`crate::product_graph::product_layer`] (no target MCS gate).
//! Root is not yielded; each child is emitted once (CSMI dedup) with the hop
//! path from the reactant.

use std::collections::{HashSet, VecDeque};

use crate::ForestError;
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::product_graph::{ProductGraphConfig, ProductHop, product_layer};
use crate::ruleset::RuleSet;

/// Frontier order for [`enumerate_metabolites`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EnumOrder {
    /// Queue: finish shallower depths before deeper ones.
    Bfs,
    /// Stack: expand the newest child before the rest of that generation.
    Dfs,
}

/// One reaction hop on an enumeration path.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathHop {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    pub site_tags: Vec<Tag>,
    pub added_tags: Vec<Tag>,
    /// All fragment CSMIs from this emission (sorted; len≥2 when cleaving).
    pub products: Vec<String>,
    pub cleaves: bool,
}

impl From<&ProductHop> for PathHop {
    fn from(hop: &ProductHop) -> Self {
        Self {
            rule: hop.rule.clone(),
            pattern_name: hop.pattern_name.clone(),
            site: hop.site,
            site_orbit: hop.site_orbit.clone(),
            site_tags: hop.site_tags.clone(),
            added_tags: hop.added_tags.clone(),
            products: hop.products.clone(),
            cleaves: hop.cleaves,
        }
    }
}

/// Path from the reactant to a yielded metabolite.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathInfo {
    /// Hops from the reactant (length == [`Self::depth`]).
    pub hops: Vec<PathHop>,
}

impl PathInfo {
    pub fn depth(&self) -> usize {
        self.hops.len()
    }

    pub fn last(&self) -> Option<&PathHop> {
        self.hops.last()
    }

    pub fn rules(&self) -> impl Iterator<Item = &str> {
        self.hops.iter().map(|h| h.rule.as_str())
    }
}

/// One streamed metabolite with its path from the reactant.
///
/// No `Debug` derive: [`ForestMol`] is not `Debug`.
#[derive(Clone)]
pub struct Metabolite {
    pub mol: ForestMol,
    pub path: PathInfo,
}

impl Metabolite {
    pub fn smiles(&self) -> String {
        self.mol.csmi().as_ref().to_string()
    }

    pub fn depth(&self) -> usize {
        self.path.depth()
    }
}

/// Bounds for streaming enumeration.
#[derive(Clone, Debug)]
pub struct EnumConfig {
    pub order: EnumOrder,
    /// Max reaction depth (root = 0; depth-1 children are first products).
    pub max_depth: usize,
    /// Cap on distinct CSMIs enqueued (including the root). `0` = unlimited.
    pub max_nodes: usize,
}

impl Default for EnumConfig {
    fn default() -> Self {
        Self {
            order: EnumOrder::Bfs,
            max_depth: 1,
            max_nodes: 0,
        }
    }
}

impl EnumConfig {
    pub fn bfs(max_depth: usize) -> Self {
        Self {
            order: EnumOrder::Bfs,
            max_depth,
            ..Self::default()
        }
    }

    pub fn dfs(max_depth: usize) -> Self {
        Self {
            order: EnumOrder::Dfs,
            max_depth,
            ..Self::default()
        }
    }

    pub fn with_max_nodes(mut self, max_nodes: usize) -> Self {
        self.max_nodes = max_nodes;
        self
    }
}

struct Expand {
    mol: ForestMol,
    /// Empty for the root (not yielded).
    path: PathInfo,
}

/// Pull iterator: metabolites of `ruleset` up to [`EnumConfig::max_depth`].
pub struct MetaboliteEnum<'a> {
    ruleset: &'a RuleSet,
    config: EnumConfig,
    layer: ProductGraphConfig,
    frontier: VecDeque<Expand>,
    seen: HashSet<String>,
    /// Distinct CSMIs enqueued (root counted).
    n_nodes: usize,
    err: Option<ForestError>,
    done: bool,
}

impl<'a> MetaboliteEnum<'a> {
    fn new(reactant: &str, ruleset: &'a RuleSet, config: EnumConfig) -> Result<Self, ForestError> {
        let start = ForestMol::parse(reactant)?;
        let mut seen = HashSet::new();
        if let Some(k) = start.stable_csmi_key() {
            seen.insert(k.as_ref().to_string());
        } else {
            seen.insert(start.csmi().as_ref().to_string());
        }
        let mut frontier = VecDeque::new();
        frontier.push_back(Expand {
            mol: start,
            path: PathInfo { hops: Vec::new() },
        });
        Ok(Self {
            ruleset,
            config,
            layer: ProductGraphConfig {
                target: None,
                max_nodes: usize::MAX,
                max_depth: usize::MAX,
            },
            frontier,
            seen,
            n_nodes: 1,
            err: None,
            done: false,
        })
    }

    fn child_key(child: &ForestMol) -> String {
        child
            .stable_csmi_key()
            .map(|k| k.as_ref().to_string())
            .unwrap_or_else(|| child.csmi().as_ref().to_string())
    }

    fn enqueue_children(&mut self, parent: &Expand) -> Result<(), ForestError> {
        if parent.path.depth() >= self.config.max_depth {
            return Ok(());
        }
        if self.config.max_nodes > 0 && self.n_nodes >= self.config.max_nodes {
            return Ok(());
        }
        let children = product_layer(&parent.mol, self.ruleset, &self.layer)?;
        let mut batch = Vec::new();
        for ch in children {
            if self.config.max_nodes > 0 && self.n_nodes >= self.config.max_nodes {
                break;
            }
            let key = Self::child_key(&ch.child);
            if !self.seen.insert(key) {
                continue;
            }
            self.n_nodes += 1;
            let mut hops = parent.path.hops.clone();
            hops.push(PathHop::from(&ch.hop));
            batch.push(Expand {
                mol: ch.child,
                path: PathInfo { hops },
            });
        }
        match self.config.order {
            EnumOrder::Bfs => {
                for node in batch {
                    self.frontier.push_back(node);
                }
            }
            EnumOrder::Dfs => {
                // LIFO: reverse so the first child is expanded first (Python dfs).
                for node in batch.into_iter().rev() {
                    self.frontier.push_back(node);
                }
            }
        }
        Ok(())
    }
}

impl Iterator for MetaboliteEnum<'_> {
    type Item = Result<Metabolite, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if let Some(e) = self.err.take() {
            self.done = true;
            return Some(Err(e));
        }
        while let Some(node) = match self.config.order {
            EnumOrder::Bfs => self.frontier.pop_front(),
            EnumOrder::Dfs => self.frontier.pop_back(),
        } {
            let yield_me = !node.path.hops.is_empty();
            if let Err(e) = self.enqueue_children(&node) {
                self.done = true;
                return Some(Err(e));
            }
            if yield_me {
                return Some(Ok(Metabolite {
                    mol: node.mol,
                    path: node.path,
                }));
            }
        }
        self.done = true;
        None
    }
}

/// Stream metabolites of `ruleset` from `reactant` under `config`.
pub fn enumerate_metabolites<'a>(
    reactant: &str,
    ruleset: &'a RuleSet,
    config: EnumConfig,
) -> Result<MetaboliteEnum<'a>, ForestError> {
    MetaboliteEnum::new(reactant, ruleset, config)
}

/// Breadth-first metabolites up to `max_depth` (Python `bfs`).
pub fn bfs<'a>(
    reactant: &str,
    ruleset: &'a RuleSet,
    max_depth: usize,
) -> Result<MetaboliteEnum<'a>, ForestError> {
    enumerate_metabolites(reactant, ruleset, EnumConfig::bfs(max_depth))
}

/// Depth-first metabolites up to `max_depth` (Python `dfs`).
pub fn dfs<'a>(
    reactant: &str,
    ruleset: &'a RuleSet,
    max_depth: usize,
) -> Result<MetaboliteEnum<'a>, ForestError> {
    enumerate_metabolites(reactant, ruleset, EnumConfig::dfs(max_depth))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::{dealkylation, hydroxylation, phase_one};
    use crate::ruleset::o_dealkylation;

    fn smiles(hit: &Metabolite) -> String {
        hit.smiles().to_string()
    }

    #[test]
    fn bfs_ethane_depth1_yields_ethanol() {
        let set = hydroxylation();
        let hits: Vec<_> = bfs("CC", &set, 1)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        assert!(!hits.is_empty());
        assert!(hits.iter().all(|h| h.depth() == 1));
        assert!(hits.iter().any(|h| {
            let s = h.smiles();
            s == "CCO" || s == "C(C)O"
        }));
        for h in &hits {
            assert_eq!(h.path.hops.len(), 1);
            assert_eq!(h.path.last().unwrap().rule, "Hydroxylation");
        }
    }

    #[test]
    fn dfs_and_bfs_same_products_depth1() {
        let set = hydroxylation();
        let mut bfs_s: Vec<_> = bfs("c1ccccc1", &set, 1)
            .unwrap()
            .map(|h| smiles(&h.unwrap()))
            .collect();
        let mut dfs_s: Vec<_> = dfs("c1ccccc1", &set, 1)
            .unwrap()
            .map(|h| smiles(&h.unwrap()))
            .collect();
        bfs_s.sort();
        dfs_s.sort();
        assert_eq!(bfs_s, dfs_s);
        assert_eq!(bfs_s.len(), 1); // phenol (orbit)
    }

    #[test]
    fn bfs_depth2_path_info_accumulates() {
        let set = hydroxylation();
        let hits: Vec<_> = bfs("CC", &set, 2)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        let d1: Vec<_> = hits.iter().filter(|h| h.depth() == 1).collect();
        let d2: Vec<_> = hits.iter().filter(|h| h.depth() == 2).collect();
        assert!(!d1.is_empty());
        assert!(!d2.is_empty(), "expected depth-2 diols; got {:?}", hits.iter().map(smiles).collect::<Vec<_>>());
        for h in &d2 {
            assert_eq!(h.path.hops.len(), 2);
            assert!(h.path.rules().all(|r| r == "Hydroxylation"));
        }
        // BFS: all depth-1 before any depth-2.
        let mut saw_d2 = false;
        for h in &hits {
            if h.depth() == 2 {
                saw_d2 = true;
            } else if saw_d2 {
                panic!("BFS yielded depth-1 after depth-2");
            }
        }
    }

    #[test]
    fn dfs_can_reach_depth2_before_finishing_depth1() {
        let set = hydroxylation();
        let hits: Vec<_> = dfs("CCC", &set, 2)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        assert!(hits.iter().any(|h| h.depth() == 2));
        // With a stack, a depth-2 hit may appear before every depth-1 is done.
        let first_d2 = hits.iter().position(|h| h.depth() == 2);
        let last_d1 = hits.iter().rposition(|h| h.depth() == 1);
        if let (Some(i2), Some(i1)) = (first_d2, last_d1) {
            assert!(
                i2 < i1,
                "expected DFS to interleave depths; first_d2={i2} last_d1={i1} order={:?}",
                hits.iter().map(|h| h.depth()).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn anisole_dealk_streams_phenol_with_path() {
        let set = o_dealkylation();
        let hits: Vec<_> = bfs("COc1ccccc1", &set, 1)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        let phenol = hits.iter().find(|h| {
            let s = h.smiles();
            s == "Oc1ccccc1" || s == "c1(O)ccccc1"
        });
        assert!(
            phenol.is_some(),
            "{:?}",
            hits.iter().map(smiles).collect::<Vec<_>>()
        );
        let p = phenol.unwrap();
        assert_eq!(p.depth(), 1);
        let last = p.path.last().unwrap();
        assert!(
            last.cleaves || last.rule.contains("Dealk") || last.rule == "ODealkylation",
            "rule={}",
            last.rule
        );
    }

    #[test]
    fn max_nodes_caps_enumeration() {
        let set = phase_one();
        let hits: Vec<_> = enumerate_metabolites(
            "COc1ccccc1",
            &set,
            EnumConfig::bfs(3).with_max_nodes(8),
        )
        .unwrap()
        .map(|h| h.unwrap())
        .collect();
        // Root counts as 1 → at most 7 yielded metabolites.
        assert!(hits.len() <= 7, "len={}", hits.len());
    }

    #[test]
    fn csmi_dedup_across_depths() {
        let set = dealkylation();
        let hits: Vec<_> = bfs("COc1ccc(OC)cc1", &set, 2)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        let mut seen = HashSet::new();
        for h in &hits {
            assert!(seen.insert(h.smiles().to_string()), "duplicate {}", h.smiles());
        }
    }
}
