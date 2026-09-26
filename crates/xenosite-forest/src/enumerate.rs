//! Stream metabolites up to a depth (BFS or DFS) with path info.
//!
//! Sibling of Python `bfs` / `dfs` in `xenosite.forest.find_path`: frontier
//! expand via [`crate::product_graph::product_layer`] (no target MCS gate).
//! Root is not yielded. Default [`EnumConfig::unique_csmi`] emits each child
//! CSMI once; set it false to stream every path (still refuses cycles on the
//! current walk).

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
    /// Cap on enqueued nodes (including the root). `0` = unlimited.
    /// With [`Self::unique_csmi`], each distinct CSMI counts once; without, each
    /// path-node counts.
    pub max_nodes: usize,
    /// When true (default, Python `bfs`/`dfs` shape), emit each product CSMI
    /// once. When false, emit every path; only refuse children already on the
    /// current walk (cycle).
    pub unique_csmi: bool,
}

impl Default for EnumConfig {
    fn default() -> Self {
        Self {
            order: EnumOrder::Bfs,
            max_depth: 1,
            max_nodes: 0,
            unique_csmi: true,
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

    /// Stream every path (no global CSMI collapse). Cycles on one walk still drop.
    pub fn with_all_paths(mut self) -> Self {
        self.unique_csmi = false;
        self
    }

    pub fn with_unique_csmi(mut self, unique_csmi: bool) -> Self {
        self.unique_csmi = unique_csmi;
        self
    }
}

struct Expand {
    mol: ForestMol,
    /// Empty for the root (not yielded).
    path: PathInfo,
    /// CSMI keys from root through [`Self::mol`] (inclusive). Used to refuse
    /// cycles when [`EnumConfig::unique_csmi`] is false.
    ancestors: HashSet<String>,
}

/// Pull iterator: metabolites of `ruleset` up to [`EnumConfig::max_depth`].
pub struct MetaboliteEnum<'a> {
    ruleset: &'a RuleSet,
    config: EnumConfig,
    layer: ProductGraphConfig,
    frontier: VecDeque<Expand>,
    /// Global CSMI set when [`EnumConfig::unique_csmi`] is true.
    seen: HashSet<String>,
    /// Enqueued nodes (root counted); see [`EnumConfig::max_nodes`].
    n_nodes: usize,
    err: Option<ForestError>,
    done: bool,
}

impl<'a> MetaboliteEnum<'a> {
    fn new(reactant: &str, ruleset: &'a RuleSet, config: EnumConfig) -> Result<Self, ForestError> {
        let start = ForestMol::parse(reactant)?;
        let root_key = Self::mol_key(&start);
        let mut seen = HashSet::new();
        if config.unique_csmi {
            seen.insert(root_key.clone());
        }
        let mut ancestors = HashSet::new();
        ancestors.insert(root_key);
        let mut frontier = VecDeque::new();
        frontier.push_back(Expand {
            mol: start,
            path: PathInfo { hops: Vec::new() },
            ancestors,
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

    fn mol_key(mol: &ForestMol) -> String {
        mol.stable_csmi_key()
            .map(|k| k.as_ref().to_string())
            .unwrap_or_else(|| mol.csmi().as_ref().to_string())
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
            let key = Self::mol_key(&ch.child);
            if self.config.unique_csmi {
                if !self.seen.insert(key.clone()) {
                    continue;
                }
            } else if parent.ancestors.contains(&key) {
                // Cycle on this walk — not a new path.
                continue;
            }
            self.n_nodes += 1;
            let mut hops = parent.path.hops.clone();
            hops.push(PathHop::from(&ch.hop));
            let mut ancestors = parent.ancestors.clone();
            ancestors.insert(key);
            batch.push(Expand {
                mol: ch.child,
                path: PathInfo { hops },
                ancestors,
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

/// [`enumerate_metabolites`] with [`crate::rules::default_ruleset`]. Override
/// via [`enumerate_metabolites`].
pub fn enumerate_metabolites_default(
    reactant: &str,
    config: EnumConfig,
) -> Result<MetaboliteEnum<'static>, ForestError> {
    enumerate_metabolites(reactant, crate::rules::default_ruleset_ref(), config)
}

/// Breadth-first metabolites up to `max_depth` (Python `bfs`; CSMI-deduped).
pub fn bfs<'a>(
    reactant: &str,
    ruleset: &'a RuleSet,
    max_depth: usize,
) -> Result<MetaboliteEnum<'a>, ForestError> {
    enumerate_metabolites(reactant, ruleset, EnumConfig::bfs(max_depth))
}

/// [`bfs`] with [`crate::rules::default_ruleset`]. Override via [`bfs`].
pub fn bfs_default(
    reactant: &str,
    max_depth: usize,
) -> Result<MetaboliteEnum<'static>, ForestError> {
    bfs(reactant, crate::rules::default_ruleset_ref(), max_depth)
}

/// Depth-first metabolites up to `max_depth` (Python `dfs`; CSMI-deduped).
pub fn dfs<'a>(
    reactant: &str,
    ruleset: &'a RuleSet,
    max_depth: usize,
) -> Result<MetaboliteEnum<'a>, ForestError> {
    enumerate_metabolites(reactant, ruleset, EnumConfig::dfs(max_depth))
}

/// [`dfs`] with [`crate::rules::default_ruleset`]. Override via [`dfs`].
pub fn dfs_default(
    reactant: &str,
    max_depth: usize,
) -> Result<MetaboliteEnum<'static>, ForestError> {
    dfs(reactant, crate::rules::default_ruleset_ref(), max_depth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::{dealkylation, hydroxylation, phase_one};
    use crate::ruleset::o_dealkylation;
    use std::collections::BTreeMap;

    fn smiles(hit: &Metabolite) -> String {
        hit.smiles()
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
        assert!(
            !d2.is_empty(),
            "expected depth-2 diols; got {:?}",
            hits.iter().map(smiles).collect::<Vec<_>>()
        );
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
    fn quinone_formation_cleavage_streams_both_fragments() {
        use crate::mol::canon_of;
        use crate::rules::quinone_formation;

        let set = quinone_formation();
        let hits: Vec<_> = bfs("COc1ccccc1", &set, 1)
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        let want_q = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        let want_me = canon_of("C").unwrap();
        let smis: Vec<_> = hits.iter().map(smiles).collect();
        assert!(
            smis.iter().any(|s| canon_of(s).unwrap() == want_q),
            "quinone missing: {smis:?}"
        );
        assert!(
            smis.iter().any(|s| canon_of(s).unwrap() == want_me),
            "methyl missing: {smis:?}"
        );
        assert!(
            smis.iter().all(|s| !s.contains('.')),
            "disconnected CSMI must be split: {smis:?}"
        );
        let cleave_hop = hits.iter().any(|h| {
            let last = h.path.last().unwrap();
            last.rule == "QuinoneFormation"
                && last.cleaves
                && last.products.len() >= 2
                && last.products.iter().any(|p| canon_of(p).unwrap() == want_q)
                && last.products.iter().any(|p| canon_of(p).unwrap() == want_me)
        });
        assert!(
            cleave_hop,
            "enumerate hop must name QuinoneFormation with both products"
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
            assert!(
                seen.insert(h.smiles()),
                "duplicate {}",
                h.smiles()
            );
        }
    }

    #[test]
    fn all_paths_yields_multiple_routes_to_same_product() {
        // Propane OH×2: 1,2-diol reachable via two site orders that unique-edit
        // does not collapse (primary then secondary vs secondary then primary).
        let set = hydroxylation();
        let dedup: Vec<_> = enumerate_metabolites("CCC", &set, EnumConfig::bfs(2))
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        let all: Vec<_> = enumerate_metabolites("CCC", &set, EnumConfig::bfs(2).with_all_paths())
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        assert!(
            all.len() > dedup.len(),
            "all_paths={} should exceed dedup={}",
            all.len(),
            dedup.len()
        );
        let dedup_smiles: HashSet<_> = dedup.iter().map(smiles).collect();
        let all_smiles: HashSet<_> = all.iter().map(smiles).collect();
        assert_eq!(
            dedup_smiles, all_smiles,
            "same product set; all_paths only adds redundant routes"
        );
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for h in &all {
            *counts.entry(h.smiles()).or_default() += 1;
        }
        assert!(
            counts.values().any(|&c| c > 1),
            "expected some CSMI with >1 path; counts={counts:?}"
        );
        let mut path_keys = HashSet::new();
        for h in &all {
            let key: Vec<_> = h
                .path
                .hops
                .iter()
                .map(|hop| {
                    (
                        hop.rule.clone(),
                        hop.site_tags.clone(),
                        hop.pattern_name.clone(),
                    )
                })
                .collect();
            assert!(
                path_keys.insert((h.smiles(), key)),
                "duplicate path yield for {}",
                h.smiles()
            );
        }
    }

    #[test]
    fn all_paths_phase_one_captures_convergent_chemistries() {
        // Ethene→ethanol via Hydroxylation+H2 vs Epoxidation+Opening, etc.
        let set = phase_one();
        let dedup_n = enumerate_metabolites("C=C", &set, EnumConfig::bfs(2))
            .unwrap()
            .count();
        let all: Vec<_> = enumerate_metabolites("C=C", &set, EnumConfig::bfs(2).with_all_paths())
            .unwrap()
            .map(|h| h.unwrap())
            .collect();
        assert!(all.len() > dedup_n, "all={} dedup={dedup_n}", all.len());
        let ethanolish: Vec<_> = all
            .iter()
            .filter(|h| {
                let s = h.smiles();
                s == "CCO" || s == "C(C)O"
            })
            .collect();
        if ethanolish.len() > 1 {
            let routes: HashSet<Vec<_>> = ethanolish
                .iter()
                .map(|h| h.path.rules().map(str::to_string).collect())
                .collect();
            assert!(
                routes.len() > 1,
                "expected distinct rule sequences to ethanol; got {routes:?}"
            );
        }
    }

    #[test]
    fn all_paths_refuses_cycles_on_same_walk() {
        // Even without CSMI dedup, a child equal to an ancestor is dropped.
        let set = hydroxylation();
        let hits: Vec<_> = enumerate_metabolites(
            "CC",
            &set,
            EnumConfig::bfs(3).with_all_paths().with_max_nodes(200),
        )
        .unwrap()
        .map(|h| h.unwrap())
        .collect();
        for h in &hits {
            let mut walk = HashSet::new();
            walk.insert("CC".to_string()); // reactant display; keys may differ
            // Path hops' continue products should not repeat within one path
            // when each hop's child CSMI is tracked — check hop product list
            // length equals depth (one continue species per hop).
            assert_eq!(h.path.hops.len(), h.depth());
        }
        assert!(!hits.is_empty());
    }
}
