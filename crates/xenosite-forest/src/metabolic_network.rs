//! Explored metabolic network: reactant root + metabolites + multipath hops.
//!
//! Product graph door: bidirectional adjacency, optional node/edge attrs,
//! summarization queries ([`Self::step_plan_between`]), and prune ([`Self::prune_to_targets`]).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::sync::Arc;

use crate::ForestError;
use crate::atom_diff::{AtomDiffResidual, atom_diff, residual_from_diff};
use crate::canonical_plan::{CleavageSide, Deps, Maybe, Step, as_deps, identity_plan_on_forest};
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::mol::{parse_mol, stable_csmi_key, stable_csmi_key_of};

pub type NodeIdx = usize;

/// Mutable attribute value on nodes / edges (NetworkX-style extension).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphValue {
    Bool(bool),
    I64(i64),
    String(String),
}

pub type AttrMap = BTreeMap<String, GraphValue>;

/// Well-known node attribute keys.
pub const NODE_SEALED: &str = "sealed";
pub const NODE_EXPANDED: &str = "expanded";

/// Well-known edge attribute keys (hop / emission evidence).
pub const EDGE_RULE: &str = "rule";
pub const EDGE_PATTERN: &str = "pattern_name";
pub const EDGE_SITE: &str = "site";
pub const EDGE_CLEAVES: &str = "cleaves";

/// One parent→child hop (legacy flat record + attrs).
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
    /// Elementary plan steps for this hop (composite leaves may emit >1).
    /// Empty → [`MetabolicNetwork::step_plan_between`] synthesizes identity /
    /// catalog plan from parent mol + site.
    pub plan: Vec<Step>,
    /// Discarded cleavage fragment CSMIs (kept side is the child node).
    pub discarded_sides: Vec<String>,
    /// Parent site atoms for this hop (Maybe bags / plan synthesis).
    pub site_atoms: Vec<usize>,
}

impl MetabolicHop {
    /// Elementary steps for summarization: stored plan, else catalog / identity.
    pub fn elementary_steps(&self, parent: &ForestMol) -> Vec<Step> {
        if !self.plan.is_empty() {
            return self.plan.clone();
        }
        let atoms: Vec<usize> = if self.site_atoms.is_empty() {
            vec![self.site]
        } else {
            self.site_atoms.clone()
        };
        let orbit = if self.site_orbit.is_empty() {
            atoms.clone()
        } else {
            self.site_orbit.clone()
        };
        match crate::rules::leaf_rule(&self.rule).filter(|leaf| leaf.has_plan_hook()) {
            Some(leaf) => {
                let steps = leaf.canonical_plan(parent.mol(), &atoms, None);
                if steps.is_empty() {
                    identity_plan_on_forest(self.rule.as_str(), parent, atoms, orbit)
                } else {
                    steps
                }
            }
            None => identity_plan_on_forest(self.rule.as_str(), parent, atoms, orbit),
        }
    }
}

#[derive(Clone, Debug)]
pub struct InboundEdge {
    pub parent_idx: NodeIdx,
    pub kept: Arc<ForestMol>,
    pub hop: MetabolicHop,
    pub attrs: AttrMap,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OutboundRef {
    pub child_idx: NodeIdx,
    pub inbound_slot: usize,
}

/// One molecule node.
#[derive(Clone, Debug)]
pub struct MetabolicNode {
    pub csmi: String,
    pub mol: Arc<ForestMol>,
    pub inbound: Vec<InboundEdge>,
    pub outbound: Vec<OutboundRef>,
    pub attrs: AttrMap,
}

impl MetabolicNode {
    pub fn n_inbound(&self) -> usize {
        self.inbound.len()
    }

    pub fn sealed(&self) -> bool {
        self.attrs
            .get(NODE_SEALED)
            .is_some_and(|v| matches!(v, GraphValue::Bool(true)))
    }

    pub fn expanded(&self) -> bool {
        self.attrs
            .get(NODE_EXPANDED)
            .is_some_and(|v| matches!(v, GraphValue::Bool(true)))
    }

    /// Legacy view: `(parent CSMI, hop)` for callers not yet on indices.
    pub fn via(&self, nodes: &[MetabolicNode]) -> Vec<(String, MetabolicHop)> {
        self.inbound
            .iter()
            .map(|e| {
                let p_csmi = nodes[e.parent_idx].csmi.clone();
                (p_csmi, e.hop.clone())
            })
            .collect()
    }
}

/// Explored network: reactant root + recorded multipath hops.
#[derive(Clone, Debug, Default)]
pub struct MetabolicNetwork {
    pub root_csmi: Option<String>,
    root: Option<NodeIdx>,
    pub targets: HashSet<NodeIdx>,
    pub nodes: Vec<MetabolicNode>,
    index: HashMap<Arc<str>, NodeIdx>,
}

impl MetabolicNetwork {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn root_idx(&self) -> Option<NodeIdx> {
        self.root
    }

    fn vertex_key(mol: &ForestMol) -> Option<Arc<str>> {
        mol.stable_csmi_key()
            .map(|k| Arc::from(k.as_ref()))
            .or_else(|| Some(Arc::from(mol.csmi().as_ref())))
    }

    pub fn ensure_root(&mut self, csmi: impl Into<String>) -> NodeIdx {
        let csmi = csmi.into();
        if self.root_csmi.is_none() {
            self.root_csmi = Some(csmi.clone());
        }
        let idx = self.ensure_node_csmi(csmi);
        self.root = Some(idx);
        idx
    }

    /// Root reactant with preserved ForestMol tags (BFS / product exploration).
    pub fn ensure_root_mol(&mut self, mol: ForestMol) -> NodeIdx {
        let csmi = mol.csmi().as_ref().to_string();
        if self.root_csmi.is_none() {
            self.root_csmi = Some(csmi);
        }
        let idx = self.ensure_node_mol(mol);
        self.root = Some(idx);
        idx
    }

    pub fn ensure_node(&mut self, csmi: impl Into<String>) -> NodeIdx {
        self.ensure_node_csmi(csmi.into())
    }

    pub fn ensure_node_mol(&mut self, mol: ForestMol) -> NodeIdx {
        if let Some(key) = Self::vertex_key(&mol) {
            if let Some(&i) = self.index.get(key.as_ref()) {
                return i;
            }
            let csmi = mol.csmi().as_ref().to_string();
            let i = self.push_node(csmi, Arc::new(mol));
            self.index.insert(key, i);
            return i;
        }
        let csmi = mol.csmi().as_ref().to_string();
        self.push_node(csmi, Arc::new(mol))
    }

    fn ensure_node_csmi(&mut self, csmi: String) -> NodeIdx {
        if let Some(&i) = self.index.get(csmi.as_str()) {
            return i;
        }
        let mol = Arc::new(
            ForestMol::parse(&csmi).unwrap_or_else(|_| ForestMol::parse("C").expect("fallback")),
        );
        if let Some(key) = stable_csmi_key(mol.mol()).map(|k| Arc::<str>::from(k.as_str())) {
            if let Some(&i) = self.index.get(key.as_ref()) {
                return i;
            }
            let i = self.push_node(csmi, mol);
            self.index.insert(key, i);
            return i;
        }
        self.push_node(csmi, mol)
    }

    fn push_node(&mut self, csmi: String, mol: Arc<ForestMol>) -> NodeIdx {
        let i = self.nodes.len();
        self.nodes.push(MetabolicNode {
            csmi,
            mol,
            inbound: Vec::new(),
            outbound: Vec::new(),
            attrs: AttrMap::new(),
        });
        i
    }

    pub fn mark_target(&mut self, idx: NodeIdx) {
        if idx < self.nodes.len() {
            self.targets.insert(idx);
        }
    }

    pub fn record_hop(&mut self, parent_csmi: &str, child_csmi: &str, hop: MetabolicHop) {
        let parent_i = self.ensure_node(parent_csmi.to_string());
        let child_i = self.ensure_node(child_csmi.to_string());
        self.record_hop_idx(parent_i, child_i, hop, None);
    }

    pub fn record_hop_idx(
        &mut self,
        parent_i: NodeIdx,
        child_i: NodeIdx,
        hop: MetabolicHop,
        kept: Option<Arc<ForestMol>>,
    ) {
        let kept = kept.unwrap_or_else(|| self.nodes[child_i].mol.clone());
        let slot = self.nodes[child_i].inbound.len();
        let mut attrs = AttrMap::new();
        attrs.insert(
            EDGE_RULE.into(),
            GraphValue::String(hop.rule.clone()),
        );
        attrs.insert(
            EDGE_PATTERN.into(),
            GraphValue::String(hop.pattern_name.clone()),
        );
        attrs.insert(EDGE_SITE.into(), GraphValue::I64(hop.site as i64));
        attrs.insert(
            EDGE_CLEAVES.into(),
            GraphValue::Bool(hop.cleaves),
        );
        self.nodes[child_i].inbound.push(InboundEdge {
            parent_idx: parent_i,
            kept,
            hop,
            attrs,
        });
        self.nodes[parent_i]
            .outbound
            .push(OutboundRef {
                child_idx: child_i,
                inbound_slot: slot,
            });
    }

    pub fn mark_expanded(&mut self, csmi: &str) {
        if let Some(i) = self.index_of(csmi) {
            self.nodes[i]
                .attrs
                .insert(NODE_EXPANDED.into(), GraphValue::Bool(true));
        }
    }

    pub fn mark_sealed(&mut self, csmi: &str) {
        if let Some(i) = self.index_of(csmi) {
            self.nodes[i]
                .attrs
                .insert(NODE_SEALED.into(), GraphValue::Bool(true));
        }
    }

    pub fn is_sealed(&self, csmi: &str) -> bool {
        self.index_of(csmi)
            .map(|i| self.nodes[i].sealed())
            .unwrap_or(false)
    }

    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn n_edges(&self) -> usize {
        self.nodes.iter().map(|n| n.inbound.len()).sum()
    }

    pub fn n_rule_patterns(&self) -> usize {
        let mut keys = BTreeSet::new();
        for n in &self.nodes {
            for e in &n.inbound {
                keys.insert((e.hop.rule.as_str(), e.hop.pattern_name.as_str()));
            }
        }
        keys.len()
    }

    pub fn index_of(&self, csmi: &str) -> Option<NodeIdx> {
        self.index.get(csmi).copied().or_else(|| {
            stable_csmi_key_of(csmi).and_then(|k| self.index.get(k.as_str()).copied())
        })
    }

    pub fn reaches(&self, target_csmi: &str) -> bool {
        self.index_of(target_csmi).is_some()
    }

    pub fn parents(&self, csmi: &str) -> Vec<&str> {
        self.index_of(csmi)
            .map(|i| {
                self.nodes[i]
                    .inbound
                    .iter()
                    .map(|e| self.nodes[e.parent_idx].csmi.as_str())
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn children(&self, csmi: &str) -> Vec<&str> {
        let Some(i) = self.index_of(csmi) else {
            return Vec::new();
        };
        self.nodes[i]
            .outbound
            .iter()
            .map(|o| self.nodes[o.child_idx].csmi.as_str())
            .collect()
    }

    /// Remove forward dead-ends that are not marked targets (§3j).
    pub fn prune_to_targets(&mut self) -> usize {
        let root = match self.root {
            Some(r) => r,
            None => return 0,
        };
        let mut queue: VecDeque<NodeIdx> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(i, n)| {
                *i != root && n.outbound.is_empty() && !self.targets.contains(i)
            })
            .map(|(i, _)| i)
            .collect();
        let mut removed = 0usize;
        while let Some(n) = queue.pop_front() {
            if n >= self.nodes.len() || self.targets.contains(&n) || self.root == Some(n) {
                continue;
            }
            if !self.nodes[n].outbound.is_empty() {
                continue;
            }
            self.remove_node(n);
            removed += 1;
            for (pi, parent) in self.nodes.iter().enumerate() {
                if parent.outbound.is_empty()
                    && pi != root
                    && !self.targets.contains(&pi)
                {
                    queue.push_back(pi);
                }
            }
        }
        removed
    }

    fn remove_node(&mut self, n: NodeIdx) {
        if n >= self.nodes.len() {
            return;
        }
        for e in self.nodes[n].inbound.clone() {
            if e.parent_idx < self.nodes.len() {
                self.nodes[e.parent_idx].outbound.retain(|o| o.child_idx != n);
            }
        }
        for o in self.nodes[n].outbound.clone() {
            if o.child_idx < self.nodes.len() {
                self.nodes[o.child_idx]
                    .inbound
                    .retain(|e| e.parent_idx != n);
            }
        }
        self.compact_remove(&[n]);
    }

    fn compact_remove(&mut self, to_drop: &[NodeIdx]) {
        let drop: HashSet<_> = to_drop.iter().copied().collect();
        if drop.is_empty() {
            return;
        }
        let mut new_nodes = Vec::new();
        let mut remap = HashMap::new();
        for (i, node) in self.nodes.iter().enumerate() {
            if drop.contains(&i) {
                continue;
            }
            remap.insert(i, new_nodes.len());
            new_nodes.push(node.clone());
        }
        for node in &mut new_nodes {
            node.inbound.retain(|e| !drop.contains(&e.parent_idx));
            for e in &mut node.inbound {
                e.parent_idx = remap[&e.parent_idx];
            }
            node.outbound.retain(|o| !drop.contains(&o.child_idx));
            for o in &mut node.outbound {
                o.child_idx = remap[&o.child_idx];
            }
        }
        self.nodes = new_nodes;
        self.index.clear();
        for (i, n) in self.nodes.iter().enumerate() {
            self.index.insert(Arc::from(n.csmi.as_str()), i);
            if let Some(k) = stable_csmi_key(n.mol.mol()) {
                self.index.insert(Arc::from(k.as_str()), i);
            }
        }
        self.root = self
            .root_csmi
            .as_ref()
            .and_then(|r| self.index.get(r.as_str()).copied());
        self.targets = self.targets.iter().filter_map(|&i| remap.get(&i).copied()).collect();
    }

    /// Summarization query: collapse recorded routes between two nodes (v1).
    ///
    /// Reconstructs one `from → to` path via **backward BFS** on inbound edges
    /// (forward-greedy outbound fails under multipath). Concatenates each hop's
    /// elementary steps; cleavage `discarded_sides` become [`Maybe`] bags
    /// (opens empty in v1 — callers may overlay richer walk Maybe).
    pub fn step_plan_between(&self, from: NodeIdx, to: NodeIdx) -> Result<Deps, ForestError> {
        if from >= self.nodes.len() || to >= self.nodes.len() {
            return Err(ForestError::Plan("invalid node index".into()));
        }
        if from == to {
            return Ok(as_deps([]));
        }
        // back_edge[parent] = (child, inbound_slot on child): forward hop parent→child.
        let mut back_edge: HashMap<NodeIdx, (NodeIdx, usize)> = HashMap::new();
        let mut seen = HashSet::from([to]);
        let mut q = VecDeque::from([to]);
        let mut found = false;
        'bfs: while let Some(cur) = q.pop_front() {
            if cur >= self.nodes.len() {
                continue;
            }
            for (slot, e) in self.nodes[cur].inbound.iter().enumerate() {
                let parent = e.parent_idx;
                if !seen.insert(parent) {
                    continue;
                }
                back_edge.insert(parent, (cur, slot));
                if parent == from {
                    found = true;
                    break 'bfs;
                }
                q.push_back(parent);
            }
        }
        if !found {
            return Err(ForestError::Plan("no path between nodes".into()));
        }
        let mut steps = Vec::new();
        let mut maybe_entries = Vec::new();
        let mut cur = from;
        while cur != to {
            let Some(&(child, slot)) = back_edge.get(&cur) else {
                return Err(ForestError::Plan("could not synthesize path steps".into()));
            };
            if child >= self.nodes.len() || slot >= self.nodes[child].inbound.len() {
                return Err(ForestError::Plan("could not synthesize path steps".into()));
            }
            let e = &self.nodes[child].inbound[slot];
            let parent_mol = &self.nodes[e.parent_idx].mol;
            steps.extend(e.hop.elementary_steps(parent_mol));
            let site: Vec<usize> = if e.hop.site_atoms.is_empty() {
                vec![e.hop.site]
            } else {
                e.hop.site_atoms.clone()
            };
            for side in &e.hop.discarded_sides {
                maybe_entries.push(CleavageSide::new(
                    site.iter().copied(),
                    side.clone(),
                    std::iter::empty::<BTreeSet<usize>>(),
                ));
            }
            cur = child;
        }
        Ok(as_deps(steps).with_maybe(Maybe::new(maybe_entries)))
    }

    pub fn closest(&self, target: &str, k: usize) -> Result<Vec<(String, usize)>, ForestError> {
        let tmol = parse_mol(target)?;
        let mut scored: Vec<(String, usize)> = Vec::new();
        for n in &self.nodes {
            if n.inbound.is_empty() && Some(n.csmi.as_str()) != self.root_csmi.as_deref() {
                continue;
            }
            let d = atom_diff(n.mol.mol(), &tmol);
            scored.push((n.csmi.clone(), d.cost()));
        }
        scored.sort_by_key(|(_, c)| *c);
        scored.truncate(k);
        Ok(scored)
    }

    pub fn missed(
        &self,
        csmi: &str,
        target: &str,
    ) -> Result<Option<AtomDiffResidual>, ForestError> {
        let Some(i) = self.index_of(csmi) else {
            return Ok(None);
        };
        let rmol = self.nodes[i].mol.mol();
        let tmol = parse_mol(target)?;
        let d = atom_diff(rmol, &tmol);
        Ok(Some(residual_from_diff(&d, Some(rmol), Some(&tmol))))
    }

    pub fn cost_to(&self, csmi: &str, target: &str) -> Result<Option<usize>, ForestError> {
        if self.index_of(csmi).is_none() {
            return Ok(None);
        }
        let rmol = parse_mol(csmi)?;
        let tmol = parse_mol(target)?;
        Ok(Some(atom_diff(&rmol, &tmol).cost()))
    }
}

#[allow(dead_code)]
pub(crate) fn network_from_nodes(
    nodes: Vec<MetabolicNode>,
    root_csmi: Option<String>,
) -> MetabolicNetwork {
    let mut index = HashMap::new();
    for (i, n) in nodes.iter().enumerate() {
        if let Some(k) = stable_csmi_key(n.mol.mol()) {
            index.insert(Arc::from(k.as_str()), i);
        }
        index.insert(Arc::from(n.csmi.as_str()), i);
    }
    MetabolicNetwork {
        root_csmi: root_csmi.clone(),
        root: root_csmi
            .as_ref()
            .and_then(|r| index.get(r.as_str()).copied()),
        targets: HashSet::new(),
        nodes,
        index,
    }
}

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
        plan: Vec::new(),
        discarded_sides: Vec::new(),
        site_atoms: vec![site],
    }
}

pub fn tags_for_atoms(mol: &ForestMol, atoms: &BTreeSet<usize>) -> Vec<Tag> {
    let mut tags: Vec<Tag> = atoms.iter().filter_map(|&a| mol.tag_of(a)).collect();
    tags.sort_unstable();
    tags.dedup();
    tags
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ensure_root_and_multipath_hop() {
        let mut net = MetabolicNetwork::new();
        net.ensure_root("CCO");
        let hop = hop_from_parts(
            "Hydroxylation",
            "h",
            0,
            vec![0],
            vec![],
            vec![],
            vec!["CCO".into()],
            false,
        );
        net.record_hop("CCO", "CCO", hop.clone());
        net.record_hop("CCO", "CCO", hop);
        assert!(net.reaches("CCO"));
        assert_eq!(net.n_nodes(), 1);
        assert_eq!(net.n_edges(), 2);
    }

    #[test]
    fn prune_removes_forward_leaf() {
        let mut net = MetabolicNetwork::new();
        let root = net.ensure_root("A");
        let _ = net.ensure_node("B");
        net.record_hop("A", "B", hop_from_parts("r", "p", 0, vec![0], vec![], vec![], vec!["B".into()], false));
        net.mark_target(root);
        let removed = net.prune_to_targets();
        assert_eq!(removed, 1);
        assert_eq!(net.n_nodes(), 1);
        assert!(net.reaches("A"));
        assert!(!net.reaches("B"));
    }

    #[test]
    fn step_plan_between_ethane_hydroxylation() {
        let mut net = MetabolicNetwork::new();
        let root = ForestMol::parse("CC").unwrap();
        let child = ForestMol::parse("CCO").unwrap();
        let ri = net.ensure_root_mol(root.copy_mol());
        let ci = net.ensure_node_mol(child.copy_mol());
        let mut hop = hop_from_parts(
            "Hydroxylation",
            "h",
            0,
            vec![0, 1],
            vec![],
            vec![],
            vec![child.csmi().as_ref().to_string()],
            false,
        );
        hop.plan = crate::canonical_plan::identity_plan_on_forest(
            "Hydroxylation",
            &root,
            [0],
            [0, 1],
        );
        hop.site_atoms = vec![0];
        net.record_hop_idx(ri, ci, hop, Some(std::sync::Arc::new(child)));
        let plan = net.step_plan_between(ri, ci).unwrap();
        assert_eq!(plan.len(), 1);
        assert_eq!(plan.steps()[0].rule, "Hydroxylation");
        assert!(net.step_plan_between(ri, ri).unwrap().is_empty());
    }
}
