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

use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;

use crate::ForestError;
use crate::atom_diff::atom_diff;
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::metabolic_network::{MetabolicHop, MetabolicNetwork, NodeIdx};
use crate::mol::{Molecule, canon_of};
use crate::ruleset::RuleSet;

/// Compat alias — hops are stored on [`MetabolicNetwork`].
pub type ProductHop = MetabolicHop;

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

/// One-hop product layer from a tagged parent (candidates + pairs).
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

    let mut out = Vec::new();

    for c in ruleset.candidates(parent) {
        let c = c?;
        if let (Some(diff), Some(t)) = (&parent_diff, &target_mol)
            && !c.could_help_on(diff, Some(t))
        {
            continue;
        }
        let pieces = c.materialize_mols()?;
        if pieces.is_empty() {
            continue;
        }
        let products: Vec<String> = pieces.iter().map(crate::mol::canon_smiles).collect();
        let mut products_sorted = products.clone();
        products_sorted.sort();
        let site_atoms: Vec<usize> = {
            let mut atoms: BTreeSet<usize> = c
                .pattern
                .site_map
                .iter()
                .filter_map(|m| c.mapped.get(m).copied())
                .collect();
            if atoms.is_empty() {
                atoms.insert(c.site);
            }
            atoms.into_iter().collect()
        };
        let site_tags = site_tags_of(parent, &site_atoms);
        let rule = c.leaf_rule().unwrap_or(c.pattern.name.as_str()).to_string();
        let cleaves = c.pattern.effect.cleaves && pieces.len() >= 2;

        for piece in pieces {
            let child = parent.from_edit_product(piece);
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
            if config.target.is_some() && !expand && !cleaves {
                // Non-cleaving dead-ends toward the target are not recorded.
                continue;
            }
            if cleaves {
                // Cleavage: keep fragment as a node only if expandable or no target.
                if config.target.is_some() && !expand {
                    continue;
                }
            }
            let added_tags = added_tags_of(parent, &child);
            let plan = c.elementary_plan();
            out.push(ProductChild {
                hop: MetabolicHop {
                    rule: rule.clone(),
                    pattern_name: c.pattern.name.clone(),
                    site: c.site,
                    site_orbit: c.orbit.clone(),
                    site_tags: site_tags.clone(),
                    added_tags,
                    products: products_sorted.clone(),
                    cleaves,
                    plan,
                    discarded_sides: Vec::new(),
                    site_atoms: site_atoms.clone(),
                },
                child,
                expand: config.target.is_none() || expand,
            });
        }
    }

    for pair in ruleset
        .candidates(parent)
        .filter(|c| matches!(c, Ok(s) if s.is_pair()))
    {
        let pair = pair?;
        if let (Some(diff), Some(t)) = (&parent_diff, &target_mol)
            && !pair.could_help_on(diff, Some(t))
        {
            continue;
        }
        push_pair_children(
            parent,
            &pair,
            parent_diff.as_ref(),
            parent_ha,
            target_mol.as_ref(),
            target_csmi.as_deref(),
            target_ha,
            config.target.is_some(),
            &mut out,
        )?;
    }

    Ok(out)
}

#[allow(clippy::too_many_arguments)] // target gate + parent walk; keep flat
fn push_pair_children(
    parent: &ForestMol,
    pair: &crate::candidate::DeferredSite,
    parent_diff: Option<&crate::atom_diff::AtomDiff>,
    parent_ha: usize,
    target_mol: Option<&Molecule>,
    target_csmi: Option<&str>,
    target_ha: Option<usize>,
    have_target: bool,
    out: &mut Vec<ProductChild>,
) -> Result<(), ForestError> {
    let pieces = pair.materialize_mols()?;
    if pieces.is_empty() {
        return Ok(());
    }
    let products: Vec<String> = pieces.iter().map(crate::mol::canon_smiles).collect();
    let mut products_sorted = products;
    products_sorted.sort();
    let site_atoms = pair.plan_site_atoms();
    let site_tags = site_tags_of(parent, &site_atoms);
    let rule = pair
        .pattern_name
        .split('+')
        .next()
        .unwrap_or("Pair")
        .to_string();
    let cleaves = pair.effect.cleaves && pieces.len() >= 2;

    for piece in pieces {
        let child = parent.from_edit_product(piece);
        let child_csmi = child.csmi().as_ref().to_string();
        let expand = child_worth_expanding(
            parent_diff,
            parent_ha,
            &child,
            &child_csmi,
            target_mol,
            target_csmi,
            target_ha,
        );
        if have_target && !expand && !cleaves {
            continue;
        }
        if cleaves && have_target && !expand {
            continue;
        }
        let added_tags = added_tags_of(parent, &child);
        let plan = pair.elementary_plan();
        out.push(ProductChild {
            hop: MetabolicHop {
                rule: rule.clone(),
                pattern_name: pair.pattern_name.clone(),
                site: pair.site,
                site_orbit: vec![pair.site],
                site_tags: site_tags.clone(),
                added_tags,
                products: products_sorted.clone(),
                cleaves,
                plan,
                discarded_sides: Vec::new(),
                site_atoms: site_atoms.clone(),
            },
            child,
            expand: !have_target || expand,
        });
    }
    Ok(())
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

/// BFS product exploration from `start` under `config` (same type as search ``network=``).
pub fn product_graph(
    start: &str,
    ruleset: &RuleSet,
    config: &ProductGraphConfig,
) -> Result<MetabolicNetwork, ForestError> {
    let mut net = MetabolicNetwork::new();
    product_graph_into(&mut net, start, ruleset, config)?;
    Ok(net)
}

/// Grow an existing [`MetabolicNetwork`] via BFS from `start` (multipath edges, tagged mols).
pub fn product_graph_into(
    net: &mut MetabolicNetwork,
    start: &str,
    ruleset: &RuleSet,
    config: &ProductGraphConfig,
) -> Result<(), ForestError> {
    let root = ForestMol::parse(start)?;
    let root_i = if net.n_nodes() == 0 {
        net.ensure_root_mol(root.clone())
    } else {
        net.ensure_node_mol(root.clone())
    };

    let mut queue: VecDeque<(NodeIdx, usize, ForestMol)> = VecDeque::new();
    let mut scheduled: BTreeSet<NodeIdx> = BTreeSet::new();
    queue.push_back((root_i, 0, root));
    scheduled.insert(root_i);

    while let Some((parent_i, depth, parent)) = queue.pop_front() {
        if depth >= config.max_depth || net.n_nodes() >= config.max_nodes {
            continue;
        }
        let layer = product_layer(&parent, ruleset, config)?;

        let mut ranked = layer;
        ranked.sort_by_key(|c| !c.expand);

        for item in ranked {
            let child_csmi = item.child.csmi().as_ref().to_string();
            if net.index_of(&child_csmi).is_none() && net.n_nodes() >= config.max_nodes {
                continue;
            }
            let child_i = net.ensure_node_mol(item.child.clone());
            let kept = Arc::new(item.child.clone());
            net.record_hop_idx(parent_i, child_i, item.hop, Some(kept));

            if item.expand && depth + 1 < config.max_depth && scheduled.insert(child_i) {
                queue.push_back((child_i, depth + 1, item.child));
            }
        }
    }
    Ok(())
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
) -> Result<(ProductGraphStats, MetabolicNetwork), ForestError> {
    let config = ProductGraphConfig {
        target: target.map(str::to_string),
        max_nodes,
        max_depth,
    };
    let net = product_graph(start, ruleset, &config)?;
    let reaches = match target {
        Some(t) => {
            let tc = canon_of(t)?;
            net.reaches(&tc)
        }
        None => false,
    };
    Ok((
        ProductGraphStats {
            n_nodes: net.n_nodes(),
            n_edges: net.n_edges(),
            n_rule_patterns: net.n_rule_patterns(),
            reaches_target: reaches,
        },
        net,
    ))
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
        let (stats, _net) =
            product_graph_stats("COc1ccccc1", Some("Oc1ccccc1"), &o_dealkylation(), 64, 4).unwrap();
        assert!(stats.reaches_target, "n_nodes={}", stats.n_nodes);
        assert!(stats.n_nodes >= 2);
        assert!(stats.n_edges >= 1);
    }

    #[test]
    fn ethane_to_ethanol_graph_one_hop() {
        let net = product_graph(
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
        assert!(net.reaches(&ethanol));
        let ei = net.index_of(&ethanol).unwrap();
        assert_eq!(net.nodes[ei].inbound.len(), 1);
        assert_eq!(net.nodes[ei].inbound[0].hop.rule, "Hydroxylation");
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
