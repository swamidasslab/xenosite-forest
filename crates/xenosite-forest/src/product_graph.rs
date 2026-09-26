//! Product graph: molecules as nodes, reaction applications as edges.
//!
//! Generalizes the cleavage-only [`crate::cleavage_graph`] to **all** PhaseOne
//! edits (hydroxylation, DH/QF, dealkylation, …). Cleaving bifurcations still
//! keep both fragment CSMIs as nodes when the expand gate passes.
//!
//! This is a measurement / structure door for multipath redundancy and for a
//! later diversity term that can read the graph — not a find_path phase and
//! not a revival of archive NetworkX [`MetaboliteNetwork`] (DROPPED).
//!
//! MCS / [`crate::atom_diff`] gates expansion toward an optional target: strict
//! cost drop and ha ≥ target (same idea as the cleavage product graph).

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::ForestError;
use crate::atom_diff::atom_diff;
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::mol::{Molecule, canon_of};
use crate::ruleset::RuleSet;

/// One parent→child hop recorded on the product graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProductHop {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    /// Parent forest Tags of the discovery site atoms (sorted).
    pub site_tags: Vec<Tag>,
    /// Tags minted for atoms this product gained (sorted).
    pub added_tags: Vec<Tag>,
    /// All fragment CSMIs from this emission (sorted; len≥2 when cleaving).
    pub products: Vec<String>,
    pub cleaves: bool,
}

/// One molecule node in the product graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProductNode {
    pub csmi: String,
    /// Incoming hops: `(parent CSMI, hop)`.
    pub via: Vec<(String, ProductHop)>,
}

impl ProductNode {
    pub fn n_inbound(&self) -> usize {
        self.via.len()
    }
}

/// BFS product graph over all kept reaction products.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ProductGraph {
    pub nodes: Vec<ProductNode>,
}

impl ProductGraph {
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn n_edges(&self) -> usize {
        self.nodes.iter().map(|n| n.via.len()).sum()
    }

    /// Distinct `(rule, pattern_name)` pairs appearing on any inbound hop.
    pub fn n_rule_patterns(&self) -> usize {
        let mut keys = BTreeSet::new();
        for n in &self.nodes {
            for (_, hop) in &n.via {
                keys.insert((hop.rule.as_str(), hop.pattern_name.as_str()));
            }
        }
        keys.len()
    }

    pub fn index_of(&self, csmi: &str) -> Option<usize> {
        self.nodes.iter().position(|n| n.csmi == csmi)
    }

    /// True when `target` CSMI appears as a node.
    pub fn reaches(&self, target_csmi: &str) -> bool {
        self.index_of(target_csmi).is_some()
    }
}

/// Options for building a product layer / graph.
#[derive(Clone, Debug)]
pub struct ProductGraphConfig {
    /// When set, MCS-gate children toward this target.
    pub target: Option<String>,
    /// Max distinct product CSMIs (including the root).
    pub max_nodes: usize,
    /// Max BFS depth (root = 0).
    pub max_depth: usize,
}

impl Default for ProductGraphConfig {
    fn default() -> Self {
        Self {
            target: None,
            max_nodes: 256,
            max_depth: 6,
        }
    }
}

impl ProductGraphConfig {
    pub fn toward(target: impl Into<String>) -> Self {
        Self {
            target: Some(target.into()),
            ..Self::default()
        }
    }
}

/// One kept child from a single emission on `parent`.
///
/// No `Debug` derive: [`ForestMol`] is not `Debug` (same as [`crate::cleavage_graph::CleavageSeed`]).
#[derive(Clone)]
pub struct ProductChild {
    pub hop: ProductHop,
    pub child: ForestMol,
    /// True when BFS should expand this child (target gate).
    pub expand: bool,
}

/// One-hop product layer from a tagged parent.
///
/// One walk: [`RuleSet::metabolites`] (SMIRKS + pairs). Rule names come from
/// emission `rule_path` (leaf RuleSet stamped by that door), same as find_path —
/// no pair branch and no pattern-name aliases.
pub fn product_layer(
    parent: &ForestMol,
    ruleset: &RuleSet,
    config: &ProductGraphConfig,
) -> Result<Vec<ProductChild>, ForestError> {
    let mol = parent.mol();
    let target_mol = match &config.target {
        Some(t) => Some(crate::mol::parse_mol(t)?),
        None => None,
    };
    let target_csmi = match &config.target {
        Some(t) => Some(canon_of(t)?),
        None => None,
    };
    let target_ha = target_mol.as_ref().map(|t| {
        t.atoms()
            .filter(|(_, a)| a.element.atomic_number() > 1)
            .count()
    });
    let parent_diff = target_mol.as_ref().map(|t| atom_diff(mol, t));
    let parent_ha = parent.heavy_atom_count();
    let have_target = config.target.is_some();

    let mut out = Vec::new();

    for emission in ruleset.metabolites(mol, false) {
        let emission = emission?;
        if emission.mols.is_empty() {
            continue;
        }
        let mut products_sorted = emission.products.clone();
        products_sorted.sort();
        let cleaves = emission.cleaves && emission.mols.len() >= 2;
        let rule = emission.rule_name().to_string();
        let site_tags = site_tags_of(parent, &emission.site_atoms);

        for piece in emission.mols {
            let child = parent.adopt_product(piece);
            let child_csmi = child.csmi().as_ref().to_string();
            let expand = child_worth_expanding(
                parent_diff.as_ref(),
                parent_ha,
                &child,
                &child_csmi,
                target_mol.as_ref(),
                target_csmi.as_deref(),
                target_ha,
            );
            if have_target && !expand && !cleaves {
                continue;
            }
            if cleaves && have_target && !expand {
                continue;
            }
            let added_tags = added_tags_of(parent, &child);
            out.push(ProductChild {
                hop: ProductHop {
                    rule: rule.clone(),
                    pattern_name: emission.pattern_name.clone(),
                    site: emission.site,
                    site_orbit: emission.site_orbit.clone(),
                    site_tags: site_tags.clone(),
                    added_tags,
                    products: products_sorted.clone(),
                    cleaves,
                },
                child,
                expand: !have_target || expand,
            });
        }
    }

    Ok(out)
}

fn site_tags_of(parent: &ForestMol, site_atoms: &[usize]) -> Vec<Tag> {
    let mut tags: Vec<Tag> = site_atoms
        .iter()
        .filter_map(|&i| parent.tag_of(i))
        .collect();
    tags.sort();
    tags.dedup();
    tags
}

fn added_tags_of(parent: &ForestMol, child: &ForestMol) -> Vec<Tag> {
    let mut tags: Vec<Tag> = (0..child.mol().atom_count())
        .filter_map(|i| child.tag_of(i))
        .filter(|t| parent.index_of(*t).is_none())
        .collect();
    tags.sort();
    tags.dedup();
    tags
}

fn child_worth_expanding(
    parent_diff: Option<&crate::atom_diff::AtomDiff>,
    parent_ha: usize,
    child: &ForestMol,
    child_csmi: &str,
    target: Option<&Molecule>,
    target_csmi: Option<&str>,
    target_ha: Option<usize>,
) -> bool {
    let (Some(t), Some(tc), Some(tha)) = (target, target_csmi, target_ha) else {
        return true;
    };
    if child_csmi == tc {
        return true;
    }
    // Cleavage shrinks; a piece smaller than the target cannot reach it.
    if child.heavy_atom_count() < tha {
        return false;
    }
    match parent_diff {
        Some(parent) => {
            let child_diff = atom_diff(child.mol(), t);
            child_diff.cost() < parent.cost()
        }
        None => {
            // HA non-worsening when no MCS yet (same spirit as find_path closer).
            let phd = parent_ha.abs_diff(tha);
            let chd = child.heavy_atom_count().abs_diff(tha);
            chd <= phd
        }
    }
}

/// BFS product graph from `start` under `config`.
pub fn product_graph(
    start: &str,
    ruleset: &RuleSet,
    config: &ProductGraphConfig,
) -> Result<ProductGraph, ForestError> {
    let root = ForestMol::parse(start)?;
    let root_csmi = root.csmi().as_ref().to_string();
    let mut nodes = vec![ProductNode {
        csmi: root_csmi.clone(),
        via: Vec::new(),
    }];
    let mut index: BTreeMap<String, usize> = BTreeMap::new();
    index.insert(root_csmi, 0);

    // Queue carries ForestMol for tag continuity through adopt_product.
    let mut queue: VecDeque<(usize, usize, ForestMol)> = VecDeque::new();
    let mut scheduled: BTreeSet<usize> = BTreeSet::new();
    queue.push_back((0, 0, root));
    scheduled.insert(0);

    while let Some((ni, depth, parent)) = queue.pop_front() {
        if depth >= config.max_depth || nodes.len() >= config.max_nodes {
            continue;
        }
        let parent_csmi = nodes[ni].csmi.clone();
        let layer = product_layer(&parent, ruleset, config)?;

        // Expandable children first so dead-ends do not eat the node budget.
        let mut ranked = layer;
        ranked.sort_by_key(|c| !c.expand);

        for child in ranked {
            let child_csmi = child.child.csmi().as_ref().to_string();
            let child_i = match index.entry(child_csmi.clone()) {
                std::collections::btree_map::Entry::Occupied(o) => *o.get(),
                std::collections::btree_map::Entry::Vacant(v) => {
                    if nodes.len() >= config.max_nodes {
                        continue;
                    }
                    let i = nodes.len();
                    nodes.push(ProductNode {
                        csmi: child_csmi,
                        via: Vec::new(),
                    });
                    v.insert(i);
                    i
                }
            };

            let already = nodes[child_i].via.iter().any(|(p, h)| {
                p == &parent_csmi
                    && h.pattern_name == child.hop.pattern_name
                    && h.site == child.hop.site
                    && h.products == child.hop.products
            });
            if !already {
                nodes[child_i]
                    .via
                    .push((parent_csmi.clone(), child.hop.clone()));
            }

            if child.expand && depth + 1 < config.max_depth && scheduled.insert(child_i) {
                queue.push_back((child_i, depth + 1, child.child));
            }
        }
    }

    Ok(ProductGraph { nodes })
}

/// [`product_graph`] with [`crate::rules::default_ruleset`]. Override via
/// [`product_graph`].
pub fn product_graph_default(
    start: &str,
    config: &ProductGraphConfig,
) -> Result<ProductGraph, ForestError> {
    product_graph(start, crate::rules::default_ruleset_ref(), config)
}

/// Compact stats for benches.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ProductGraphStats {
    pub n_nodes: usize,
    pub n_edges: usize,
    pub n_rule_patterns: usize,
    pub reaches_target: bool,
}

pub fn product_graph_stats(
    start: &str,
    target: Option<&str>,
    ruleset: &RuleSet,
    max_nodes: usize,
    max_depth: usize,
) -> Result<(ProductGraphStats, ProductGraph), ForestError> {
    let config = ProductGraphConfig {
        target: target.map(str::to_string),
        max_nodes,
        max_depth,
    };
    let graph = product_graph(start, ruleset, &config)?;
    let reaches = match target {
        Some(t) => {
            let tc = canon_of(t)?;
            graph.reaches(&tc)
        }
        None => false,
    };
    Ok((
        ProductGraphStats {
            n_nodes: graph.n_nodes(),
            n_edges: graph.n_edges(),
            n_rule_patterns: graph.n_rule_patterns(),
            reaches_target: reaches,
        },
        graph,
    ))
}

/// [`product_graph_stats`] with [`crate::rules::default_ruleset`]. Override via
/// [`product_graph_stats`].
pub fn product_graph_stats_default(
    start: &str,
    target: Option<&str>,
    max_nodes: usize,
    max_depth: usize,
) -> Result<(ProductGraphStats, ProductGraph), ForestError> {
    product_graph_stats(
        start,
        target,
        crate::rules::default_ruleset_ref(),
        max_nodes,
        max_depth,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cleavage_graph::fragment_worth_expanding;
    use crate::rules::{hydroxylation, phase_one};
    use crate::ruleset::o_dealkylation;

    #[test]
    fn ethane_hydroxylation_layer_has_ethanol() {
        let parent = ForestMol::parse("CC").unwrap();
        let layer =
            product_layer(&parent, &hydroxylation(), &ProductGraphConfig::default()).unwrap();
        assert!(!layer.is_empty());
        let ethanol = canon_of("CCO").unwrap();
        assert!(
            layer.iter().any(|c| c.child.csmi().as_ref() == ethanol),
            "layer CSMIs: {:?}",
            layer
                .iter()
                .map(|c| c.child.csmi().as_ref().to_string())
                .collect::<Vec<_>>()
        );
        let hop = &layer[0].hop;
        assert!(!hop.site_tags.is_empty());
        // Hydroxylation adds O → one born tag.
        assert_eq!(hop.added_tags.len(), 1);
    }

    #[test]
    fn anisole_product_graph_reaches_phenol() {
        let (stats, graph) =
            product_graph_stats("COc1ccccc1", Some("Oc1ccccc1"), &o_dealkylation(), 64, 4).unwrap();
        assert!(
            stats.reaches_target,
            "nodes={:?}",
            graph.nodes.iter().map(|n| &n.csmi).collect::<Vec<_>>()
        );
        assert!(stats.n_nodes >= 2);
        assert!(stats.n_edges >= 1);
    }

    #[test]
    fn ethane_to_ethanol_graph_one_hop() {
        let graph = product_graph(
            "CC",
            &hydroxylation(),
            &ProductGraphConfig {
                target: Some("CCO".into()),
                max_nodes: 16,
                max_depth: 2,
            },
        )
        .unwrap();
        let ethanol = canon_of("CCO").unwrap();
        assert!(graph.reaches(&ethanol));
        let ei = graph.index_of(&ethanol).unwrap();
        assert_eq!(graph.nodes[ei].via.len(), 1);
        assert_eq!(graph.nodes[ei].via[0].1.rule, "Hydroxylation");
    }

    #[test]
    fn phase_one_toward_target_stays_bounded() {
        let (stats, _) = product_graph_stats(
            "COc1ccc(O)cc1",
            Some("O=C1C=C(O)C(=O)C(O)=C1"),
            &phase_one(),
            64,
            4,
        )
        .unwrap();
        assert!(stats.n_nodes <= 64);
        assert!(stats.n_nodes >= 1);
    }

    #[test]
    fn fragment_gate_shared_with_cleavage_graph() {
        // Sanity: product_graph reuses cleavage expand helper for shrinks.
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let target = crate::mol::parse_mol("Oc1ccccc1").unwrap();
        let tc = canon_of("Oc1ccccc1").unwrap();
        let diff = atom_diff(parent.mol(), &target);
        assert!(fragment_worth_expanding(Some(&diff), &tc, &target, &tc, 7));
    }
}
