//! Explored metabolic network: reactant root + metabolites + hops.
//!
//! User-facing explore object for [`crate::find_path`] / [`crate::find_path_partial`].
//! Evolves the former [`crate::product_graph::ProductGraph`] schema (nodes + hops).
//! Not a revival of archive NetworkX `MetaboliteNetwork` (DROPPED).

use std::collections::{BTreeMap, BTreeSet};

use crate::ForestError;
use crate::atom_diff::{AtomDiffResidual, atom_diff, residual_from_diff};
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::mol::parse_mol;

/// One parent→child hop recorded on the network.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetabolicHop {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    pub site_tags: Vec<Tag>,
    pub added_tags: Vec<Tag>,
    pub products: Vec<String>,
    pub cleaves: bool,
}

/// One molecule node (reactant or metabolite).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetabolicNode {
    pub csmi: String,
    /// Incoming hops: `(parent CSMI, hop)`.
    pub via: Vec<(String, MetabolicHop)>,
    /// Marked sterile / sealed local-min during partial search.
    pub sealed: bool,
    /// Expanded at least once during search.
    pub expanded: bool,
}

impl MetabolicNode {
    pub fn n_inbound(&self) -> usize {
        self.via.len()
    }
}

/// Explored network: always includes the reactant root when seeded.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MetabolicNetwork {
    pub root_csmi: Option<String>,
    pub nodes: Vec<MetabolicNode>,
    index: BTreeMap<String, usize>,
}

impl MetabolicNetwork {
    pub fn new() -> Self {
        Self::default()
    }

    /// Ensure root exists (call once at search start).
    pub fn ensure_root(&mut self, csmi: impl Into<String>) -> usize {
        let csmi = csmi.into();
        if self.root_csmi.is_none() {
            self.root_csmi = Some(csmi.clone());
        }
        self.ensure_node(csmi)
    }

    pub fn ensure_node(&mut self, csmi: impl Into<String>) -> usize {
        let csmi = csmi.into();
        if let Some(&i) = self.index.get(&csmi) {
            return i;
        }
        let i = self.nodes.len();
        self.nodes.push(MetabolicNode {
            csmi: csmi.clone(),
            via: Vec::new(),
            sealed: false,
            expanded: false,
        });
        self.index.insert(csmi, i);
        i
    }

    pub fn record_hop(&mut self, parent_csmi: &str, child_csmi: &str, hop: MetabolicHop) {
        let _ = self.ensure_node(parent_csmi);
        let child_i = self.ensure_node(child_csmi);
        let already = self.nodes[child_i].via.iter().any(|(p, h)| {
            p == parent_csmi
                && h.pattern_name == hop.pattern_name
                && h.site == hop.site
                && h.products == hop.products
        });
        if !already {
            self.nodes[child_i].via.push((parent_csmi.to_string(), hop));
        }
    }

    pub fn mark_expanded(&mut self, csmi: &str) {
        if let Some(&i) = self.index.get(csmi) {
            self.nodes[i].expanded = true;
        }
    }

    pub fn mark_sealed(&mut self, csmi: &str) {
        if let Some(&i) = self.index.get(csmi) {
            self.nodes[i].sealed = true;
        }
    }

    pub fn is_sealed(&self, csmi: &str) -> bool {
        self.index
            .get(csmi)
            .map(|&i| self.nodes[i].sealed)
            .unwrap_or(false)
    }

    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn n_edges(&self) -> usize {
        self.nodes.iter().map(|n| n.via.len()).sum()
    }

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
        self.index.get(csmi).copied()
    }

    pub fn reaches(&self, target_csmi: &str) -> bool {
        self.index_of(target_csmi).is_some()
    }

    /// Parent CSMIs of a node.
    pub fn parents(&self, csmi: &str) -> Vec<&str> {
        self.index_of(csmi)
            .map(|i| self.nodes[i].via.iter().map(|(p, _)| p.as_str()).collect())
            .unwrap_or_default()
    }

    /// Child CSMIs that list `csmi` as a parent.
    pub fn children(&self, csmi: &str) -> Vec<&str> {
        self.nodes
            .iter()
            .filter(|n| n.via.iter().any(|(p, _)| p == csmi))
            .map(|n| n.csmi.as_str())
            .collect()
    }

    /// Atom-diff cost of each node vs `target` SMILES; lowest first (up to `k`).
    pub fn closest(&self, target: &str, k: usize) -> Result<Vec<(String, usize)>, ForestError> {
        let tmol = parse_mol(target)?;
        let mut scored: Vec<(String, usize)> = Vec::new();
        for n in &self.nodes {
            let Ok(rmol) = parse_mol(&n.csmi) else {
                continue;
            };
            let d = atom_diff(&rmol, &tmol);
            scored.push((n.csmi.clone(), d.cost()));
        }
        scored.sort_by_key(|(_, c)| *c);
        scored.truncate(k);
        Ok(scored)
    }

    /// Structured residual of `csmi` vs `target`.
    pub fn missed(
        &self,
        csmi: &str,
        target: &str,
    ) -> Result<Option<AtomDiffResidual>, ForestError> {
        if self.index_of(csmi).is_none() {
            return Ok(None);
        }
        let rmol = parse_mol(csmi)?;
        let tmol = parse_mol(target)?;
        let d = atom_diff(&rmol, &tmol);
        Ok(Some(residual_from_diff(&d, Some(&rmol), Some(&tmol))))
    }

    /// Cost of one node vs target (parse both).
    pub fn cost_to(&self, csmi: &str, target: &str) -> Result<Option<usize>, ForestError> {
        if self.index_of(csmi).is_none() {
            return Ok(None);
        }
        let rmol = parse_mol(csmi)?;
        let tmol = parse_mol(target)?;
        Ok(Some(atom_diff(&rmol, &tmol).cost()))
    }
}

/// Rebuild index after deserializing nodes-only (internal).
#[allow(dead_code)]
pub(crate) fn network_from_nodes(
    nodes: Vec<MetabolicNode>,
    root_csmi: Option<String>,
) -> MetabolicNetwork {
    let mut index = BTreeMap::new();
    for (i, n) in nodes.iter().enumerate() {
        index.insert(n.csmi.clone(), i);
    }
    MetabolicNetwork {
        root_csmi,
        nodes,
        index,
    }
}

/// Convert a search emission hop into a network hop.
#[allow(clippy::too_many_arguments)]
pub fn hop_from_parts(
    rule: impl Into<String>,
    pattern_name: impl Into<String>,
    site: usize,
    site_orbit: Vec<usize>,
    site_tags: Vec<Tag>,
    added_tags: Vec<Tag>,
    products: Vec<String>,
    cleaves: bool,
) -> MetabolicHop {
    MetabolicHop {
        rule: rule.into(),
        pattern_name: pattern_name.into(),
        site,
        site_orbit,
        site_tags,
        added_tags,
        products,
        cleaves,
    }
}

/// Tags on `mol` for the given atom indices (sorted).
pub fn tags_for_atoms(mol: &ForestMol, atoms: &std::collections::BTreeSet<usize>) -> Vec<Tag> {
    let mut tags: Vec<Tag> = atoms.iter().filter_map(|&a| mol.tag_of(a)).collect();
    tags.sort_unstable();
    tags.dedup();
    tags
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ensure_root_and_hop() {
        let mut net = MetabolicNetwork::new();
        net.ensure_root("CCO");
        net.record_hop(
            "CCO",
            "CCO",
            hop_from_parts(
                "Hydroxylation",
                "h",
                0,
                vec![0],
                vec![],
                vec![],
                vec!["CCO".into()],
                false,
            ),
        );
        assert!(net.reaches("CCO"));
        assert_eq!(net.n_nodes(), 1);
        assert_eq!(net.root_csmi.as_deref(), Some("CCO"));
    }
}
