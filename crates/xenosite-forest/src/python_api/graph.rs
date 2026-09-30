//! [`MetabolicNetwork`] handle, graph views, and BFS product exploration.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::metabolic_network::{AttrMap, GraphValue, MetabolicNetwork, NODE_EXPANDED, NODE_SEALED};
use crate::product_graph::{
    ProductGraphConfig, product_graph as product_graph_rs, product_graph_into,
};
use crate::rules::product_graph_ruleset as product_graph_ruleset_rs;
use crate::ruleset::RuleSet;

use super::common::{NETWORK_LOCK_POISONED, network_lock_err, py_err};
use super::mol::PyForestMol;
use super::plan::PyStepPlan;
use super::rules::PyRuleSet;

fn graph_value_to_py(py: Python<'_>, v: &GraphValue) -> PyResult<Py<PyAny>> {
    match v {
        GraphValue::Bool(b) => Ok(pyo3::types::PyBool::new(py, *b)
            .to_owned()
            .into_any()
            .unbind()),
        GraphValue::I64(i) => Ok(pyo3::types::PyInt::new(py, *i)
            .to_owned()
            .into_any()
            .unbind()),
        GraphValue::String(s) => Ok(pyo3::types::PyString::new(py, s)
            .to_owned()
            .into_any()
            .unbind()),
    }
}

fn py_to_graph_value(value: &Bound<'_, PyAny>) -> PyResult<GraphValue> {
    if let Ok(b) = value.extract::<bool>() {
        return Ok(GraphValue::Bool(b));
    }
    if let Ok(i) = value.extract::<i64>() {
        return Ok(GraphValue::I64(i));
    }
    if let Ok(s) = value.extract::<String>() {
        return Ok(GraphValue::String(s));
    }
    Err(PyValueError::new_err(
        "attr value must be bool, int, or str",
    ))
}

fn attrs_to_dict(py: Python<'_>, attrs: &AttrMap) -> PyResult<Py<PyAny>> {
    let d = pyo3::types::PyDict::new(py);
    for (k, v) in attrs {
        d.set_item(k, graph_value_to_py(py, v)?)?;
    }
    Ok(d.unbind().into_any())
}

/// View of one molecule row in a [`MetabolicNetwork`].
#[pyclass(name = "GraphNode", unsendable)]
pub struct PyGraphNode {
    net: Arc<Mutex<MetabolicNetwork>>,
    idx: usize,
}

#[pymethods]
impl PyGraphNode {
    #[getter]
    fn index(&self) -> usize {
        self.idx
    }

    #[getter]
    fn csmi(&self) -> PyResult<String> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        net.nodes
            .get(self.idx)
            .map(|n| n.csmi.clone())
            .ok_or_else(|| PyValueError::new_err(format!("node index {} out of range", self.idx)))
    }

    #[getter]
    fn mol(&self) -> PyResult<PyForestMol> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let node = net.nodes.get(self.idx).ok_or_else(|| {
            PyValueError::new_err(format!("node index {} out of range", self.idx))
        })?;
        Ok(PyForestMol::wrap(node.mol.as_ref().copy_mol()))
    }

    #[getter]
    fn sealed(&self) -> PyResult<bool> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        Ok(net.nodes.get(self.idx).map(|n| n.sealed()).unwrap_or(false))
    }

    #[getter]
    fn expanded(&self) -> PyResult<bool> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        Ok(net
            .nodes
            .get(self.idx)
            .map(|n| n.expanded())
            .unwrap_or(false))
    }

    /// Read a node attr (well-known keys include ``sealed``, ``expanded``).
    fn get_attr(&self, py: Python<'_>, key: &str) -> PyResult<Option<Py<PyAny>>> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let node = net.nodes.get(self.idx).ok_or_else(|| {
            PyValueError::new_err(format!("node index {} out of range", self.idx))
        })?;
        match node.attrs.get(key) {
            Some(v) => Ok(Some(graph_value_to_py(py, v)?)),
            None => Ok(None),
        }
    }

    /// Set a node attr (bool / int / str). Mutates the live network.
    fn set_attr(&self, key: &str, value: Bound<'_, PyAny>) -> PyResult<()> {
        let mut net = self.net.lock().map_err(|_| network_lock_err())?;
        let node = net.nodes.get_mut(self.idx).ok_or_else(|| {
            PyValueError::new_err(format!("node index {} out of range", self.idx))
        })?;
        node.attrs
            .insert(key.to_string(), py_to_graph_value(&value)?);
        Ok(())
    }

    fn __contains__(&self, key: &str) -> PyResult<bool> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        Ok(net
            .nodes
            .get(self.idx)
            .is_some_and(|n| n.attrs.contains_key(key)))
    }

    fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        self.get_attr(py, key)?
            .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(key.to_string()))
    }

    fn __setitem__(&self, key: &str, value: Bound<'_, PyAny>) -> PyResult<()> {
        self.set_attr(key, value)
    }

    /// Spine + attrs for notebooks / JSON (opt-in marshal).
    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let node = net.nodes.get(self.idx).ok_or_else(|| {
            PyValueError::new_err(format!("node index {} out of range", self.idx))
        })?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("index", self.idx)?;
        d.set_item("csmi", &node.csmi)?;
        d.set_item(NODE_SEALED, node.sealed())?;
        d.set_item(NODE_EXPANDED, node.expanded())?;
        d.set_item("attrs", attrs_to_dict(py, &node.attrs)?)?;
        Ok(d.unbind().into_any())
    }

    fn n_inbound(&self) -> PyResult<usize> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        Ok(net
            .nodes
            .get(self.idx)
            .map(|n| n.inbound.len())
            .unwrap_or(0))
    }

    fn inbound_edge(&self, slot: usize) -> PyResult<PyGraphEdge> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let n_in = net
            .nodes
            .get(self.idx)
            .map(|n| n.inbound.len())
            .unwrap_or(0);
        if slot >= n_in {
            return Err(PyValueError::new_err(format!(
                "inbound slot {slot} out of range (n_inbound={n_in})"
            )));
        }
        Ok(PyGraphEdge {
            net: self.net.clone(),
            child_idx: self.idx,
            slot,
        })
    }

    fn __str__(&self) -> String {
        match (
            self.csmi(),
            self.sealed(),
            self.expanded(),
            self.n_inbound(),
        ) {
            (Ok(csmi), Ok(sealed), Ok(expanded), Ok(n_in)) => {
                super::display::truncate_display(&format!(
                    "GraphNode[{}] {}\n  sealed={} expanded={} inbound={}",
                    self.idx, csmi, sealed, expanded, n_in
                ))
            }
            _ => format!("GraphNode[{}]", self.idx),
        }
    }

    fn __repr__(&self) -> String {
        match self.csmi() {
            Ok(csmi) => format!("GraphNode({}, {:?})", self.idx, csmi),
            Err(_) => format!("GraphNode({})", self.idx),
        }
    }
}

/// One recorded parent→child hop (multipath slot on the child).
#[pyclass(name = "GraphEdge", unsendable)]
pub struct PyGraphEdge {
    net: Arc<Mutex<MetabolicNetwork>>,
    child_idx: usize,
    slot: usize,
}

#[pymethods]
impl PyGraphEdge {
    #[getter]
    fn child_index(&self) -> usize {
        self.child_idx
    }

    #[getter]
    fn inbound_slot(&self) -> usize {
        self.slot
    }

    #[getter]
    fn rule(&self) -> PyResult<String> {
        Ok(self.hop()?.rule.clone())
    }

    #[getter]
    fn pattern_name(&self) -> PyResult<String> {
        Ok(self.hop()?.pattern_name.clone())
    }

    #[getter]
    fn site(&self) -> PyResult<usize> {
        Ok(self.hop()?.site)
    }

    #[getter]
    fn parent_index(&self) -> PyResult<usize> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let edge = net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))?;
        Ok(edge.parent_idx)
    }

    #[getter]
    fn cleaves(&self) -> PyResult<bool> {
        Ok(self.hop()?.cleaves)
    }

    #[getter]
    fn kept_mol(&self) -> PyResult<PyForestMol> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let edge = net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))?;
        Ok(PyForestMol::wrap(edge.kept.as_ref().copy_mol()))
    }

    /// Parent node mol (``site`` indexes this reactant, not the kept child).
    #[getter]
    fn parent_mol(&self) -> PyResult<PyForestMol> {
        let parent_idx = self.parent_index()?;
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let node = net.nodes.get(parent_idx).ok_or_else(|| {
            PyValueError::new_err(format!("parent index {parent_idx} out of range"))
        })?;
        Ok(PyForestMol::wrap(node.mol.as_ref().copy_mol()))
    }

    fn products_csmi(&self) -> PyResult<Vec<String>> {
        Ok(self.hop()?.products.clone())
    }

    fn get_attr(&self, py: Python<'_>, key: &str) -> PyResult<Option<Py<PyAny>>> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let edge = net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))?;
        match edge.attrs.get(key) {
            Some(v) => Ok(Some(graph_value_to_py(py, v)?)),
            None => Ok(None),
        }
    }

    fn set_attr(&self, key: &str, value: Bound<'_, PyAny>) -> PyResult<()> {
        let mut net = self.net.lock().map_err(|_| network_lock_err())?;
        let edge = net.nodes[self.child_idx]
            .inbound
            .get_mut(self.slot)
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))?;
        edge.attrs
            .insert(key.to_string(), py_to_graph_value(&value)?);
        Ok(())
    }

    fn __contains__(&self, key: &str) -> PyResult<bool> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        Ok(net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .is_some_and(|e| e.attrs.contains_key(key)))
    }

    fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        self.get_attr(py, key)?
            .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(key.to_string()))
    }

    fn __setitem__(&self, key: &str, value: Bound<'_, PyAny>) -> PyResult<()> {
        self.set_attr(key, value)
    }

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        let edge = net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("child_index", self.child_idx)?;
        d.set_item("inbound_slot", self.slot)?;
        d.set_item("parent_index", edge.parent_idx)?;
        d.set_item("rule", &edge.hop.rule)?;
        d.set_item("pattern_name", &edge.hop.pattern_name)?;
        d.set_item("site", edge.hop.site)?;
        d.set_item("cleaves", edge.hop.cleaves)?;
        d.set_item("attrs", attrs_to_dict(py, &edge.attrs)?)?;
        Ok(d.unbind().into_any())
    }

    fn __str__(&self) -> String {
        match (
            self.rule(),
            self.site(),
            self.parent_index(),
            self.products_csmi(),
        ) {
            (Ok(rule), Ok(site), Ok(parent), Ok(prods)) => {
                let kept = prods.first().map(String::as_str).unwrap_or("?");
                super::display::truncate_display(&format!(
                    "GraphEdge  {} @ site={}\n  parent={} → child={}\n  kept: {}",
                    rule, site, parent, self.child_idx, kept
                ))
            }
            _ => format!("GraphEdge(child={}, slot={})", self.child_idx, self.slot),
        }
    }

    fn __repr__(&self) -> String {
        match (self.rule(), self.site(), self.parent_index()) {
            (Ok(rule), Ok(site), Ok(parent)) => {
                format!("GraphEdge({rule:?}@{site}, {parent}→{})", self.child_idx)
            }
            _ => format!("GraphEdge(child={}, slot={})", self.child_idx, self.slot),
        }
    }
}

impl PyGraphEdge {
    fn hop(&self) -> PyResult<crate::metabolic_network::MetabolicHop> {
        let net = self.net.lock().map_err(|_| network_lock_err())?;
        net.nodes[self.child_idx]
            .inbound
            .get(self.slot)
            .map(|e| e.hop.clone())
            .ok_or_else(|| PyValueError::new_err("inbound slot out of range"))
    }
}

/// Explored metabolic network (reactant root + hops). Mutated by search when
/// passed as ``network=``.
#[pyclass(name = "MetabolicNetwork", unsendable)]
pub struct PyMetabolicNetwork {
    pub(crate) inner: Arc<Mutex<MetabolicNetwork>>,
}

#[pymethods]
impl PyMetabolicNetwork {
    #[new]
    fn new() -> Self {
        Self {
            inner: Arc::new(Mutex::new(MetabolicNetwork::new())),
        }
    }

    #[getter]
    fn root_csmi(&self) -> PyResult<Option<String>> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .root_csmi
            .clone())
    }

    fn n_nodes(&self) -> PyResult<usize> {
        Ok(self.inner.lock().map_err(|_| network_lock_err())?.n_nodes())
    }

    fn n_edges(&self) -> PyResult<usize> {
        Ok(self.inner.lock().map_err(|_| network_lock_err())?.n_edges())
    }

    fn reaches(&self, target_csmi: &str) -> PyResult<bool> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .reaches(target_csmi))
    }

    fn prune_to_targets(&self) -> PyResult<usize> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .prune_to_targets())
    }

    fn step_plan_between(&self, from_idx: usize, to_idx: usize) -> PyResult<PyStepPlan> {
        let plan = self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .step_plan_between(from_idx, to_idx)
            .map_err(py_err)?;
        Ok(PyStepPlan { inner: plan })
    }

    fn closest(&self, target: &str, k: usize) -> PyResult<Vec<(String, usize)>> {
        self.inner
            .lock()
            .map_err(|_| network_lock_err())?
            .closest(target, k)
            .map_err(py_err)
    }

    fn node_csmi(&self, idx: usize) -> PyResult<String> {
        let net = self.inner.lock().map_err(|_| network_lock_err())?;
        net.nodes
            .get(idx)
            .map(|n| n.csmi.clone())
            .ok_or_else(|| PyValueError::new_err(format!("node index {idx} out of range")))
    }

    fn index_of(&self, csmi: &str) -> PyResult<Option<usize>> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .index_of(csmi))
    }

    fn root_idx(&self) -> PyResult<Option<usize>> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .root_idx())
    }

    fn n_rule_patterns(&self) -> PyResult<usize> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .n_rule_patterns())
    }

    fn mark_target(&self, idx: usize) -> PyResult<()> {
        self.inner
            .lock()
            .map_err(|_| network_lock_err())?
            .mark_target(idx);
        Ok(())
    }

    fn parents(&self, csmi: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .parents(csmi)
            .into_iter()
            .map(str::to_string)
            .collect())
    }

    fn children(&self, csmi: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .children(csmi)
            .into_iter()
            .map(str::to_string)
            .collect())
    }

    /// Indexed node view (same graph handle; no clone of mols).
    fn node(&self, idx: usize) -> PyResult<PyGraphNode> {
        let n = self.n_nodes()?;
        if idx >= n {
            return Err(PyValueError::new_err(format!(
                "node index {idx} out of range"
            )));
        }
        Ok(PyGraphNode {
            net: self.inner.clone(),
            idx,
        })
    }

    fn __len__(&self) -> PyResult<usize> {
        self.n_nodes()
    }

    fn __getitem__(&self, idx: usize) -> PyResult<PyGraphNode> {
        self.node(idx)
    }

    /// Marked target CSMIs (sorted), for notebooks / summaries.
    fn target_csmis(&self) -> PyResult<Vec<String>> {
        let net = self.inner.lock().map_err(|_| network_lock_err())?;
        let mut out: Vec<String> = net
            .targets
            .iter()
            .filter_map(|&i| net.nodes.get(i).map(|n| n.csmi.clone()))
            .collect();
        out.sort();
        Ok(out)
    }

    fn __str__(&self) -> String {
        match self.inner.lock() {
            Ok(net) => super::display::format_network(&net),
            Err(_) => "MetabolicNetwork(<lock poisoned>)".into(),
        }
    }

    fn __repr__(&self) -> String {
        match self.inner.lock() {
            Ok(net) => format!(
                "MetabolicNetwork(nodes={}, edges={})",
                net.n_nodes(),
                net.n_edges()
            ),
            Err(_) => "MetabolicNetwork(<lock poisoned>)".into(),
        }
    }

    fn missed(&self, py: Python<'_>, csmi: &str, target: &str) -> PyResult<Option<Py<PyAny>>> {
        let residual = self
            .inner
            .lock()
            .map_err(|_| network_lock_err())?
            .missed(csmi, target)
            .map_err(py_err)?;
        Ok(match residual {
            None => None,
            Some(r) => {
                let d = pyo3::types::PyDict::new(py);
                d.set_item("cost", r.cost)?;
                d.set_item("n_extra", r.n_extra)?;
                d.set_item("categories", r.categories.clone())?;
                d.set_item("unresolvable", r.unresolvable)?;
                Some(d.unbind().into_any())
            }
        })
    }
}

fn product_graph_config(
    target: Option<&str>,
    max_nodes: usize,
    max_depth: usize,
) -> ProductGraphConfig {
    ProductGraphConfig {
        target: target.map(str::to_string),
        max_nodes,
        max_depth,
    }
}

fn rules_for_product_graph(ruleset: Option<&Bound<'_, PyRuleSet>>) -> RuleSet {
    match ruleset {
        Some(rs) => rs.borrow().inner.clone(),
        None => product_graph_ruleset_rs(),
    }
}

fn mark_target_csmi(net: &mut MetabolicNetwork, target: Option<&str>) {
    if let Some(t) = target
        && let Some(i) = net.index_of(t)
    {
        net.mark_target(i);
    }
}

/// BFS product exploration into a [`MetabolicNetwork`].
///
/// Releases the GIL for the Rust BFS (including ``network=``).
#[pyfunction]
#[pyo3(name = "product_graph_bfs", signature = (
    start,
    *,
    target=None,
    max_nodes=256,
    max_depth=6,
    ruleset=None,
    network=None,
))]
pub fn product_graph_bfs(
    py: Python<'_>,
    start: &str,
    target: Option<&str>,
    max_nodes: usize,
    max_depth: usize,
    ruleset: Option<&Bound<'_, PyRuleSet>>,
    network: Option<PyRefMut<'_, PyMetabolicNetwork>>,
) -> PyResult<PyMetabolicNetwork> {
    let rules = rules_for_product_graph(ruleset);
    let config = product_graph_config(target, max_nodes, max_depth);
    let start = start.to_owned();
    let target = target.map(str::to_string);

    if let Some(net_py) = network {
        let arc = net_py.inner.clone();
        drop(net_py);
        py.detach(|| {
            let mut guard = arc.lock().map_err(|_| NETWORK_LOCK_POISONED.to_string())?;
            product_graph_into(&mut guard, &start, &rules, &config).map_err(|e| e.to_string())?;
            mark_target_csmi(&mut guard, target.as_deref());
            Ok::<_, String>(())
        })
        .map_err(PyValueError::new_err)?;
        return Ok(PyMetabolicNetwork { inner: arc });
    }

    let net = py
        .detach(|| {
            let mut net = product_graph_rs(&start, &rules, &config).map_err(|e| e.to_string())?;
            mark_target_csmi(&mut net, target.as_deref());
            Ok::<_, String>(net)
        })
        .map_err(PyValueError::new_err)?;
    Ok(PyMetabolicNetwork {
        inner: Arc::new(Mutex::new(net)),
    })
}

#[pyfunction]
#[pyo3(name = "product_graph_into", signature = (
    network,
    start,
    *,
    target=None,
    max_nodes=256,
    max_depth=6,
    ruleset=None,
))]
pub fn product_graph_into_py(
    py: Python<'_>,
    network: PyRefMut<'_, PyMetabolicNetwork>,
    start: &str,
    target: Option<&str>,
    max_nodes: usize,
    max_depth: usize,
    ruleset: Option<&Bound<'_, PyRuleSet>>,
) -> PyResult<()> {
    let rules = rules_for_product_graph(ruleset);
    let config = product_graph_config(target, max_nodes, max_depth);
    let arc = network.inner.clone();
    drop(network);
    let start = start.to_owned();
    let target = target.map(str::to_string);
    py.detach(|| {
        let mut guard = arc.lock().map_err(|_| NETWORK_LOCK_POISONED.to_string())?;
        product_graph_into(&mut guard, &start, &rules, &config).map_err(|e| e.to_string())?;
        mark_target_csmi(&mut guard, target.as_deref());
        Ok::<_, String>(())
    })
    .map_err(PyValueError::new_err)?;
    Ok(())
}
