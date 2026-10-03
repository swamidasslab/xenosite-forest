//! find_path / find_path_partial and the shared ``PathOutcome`` pyclass.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubs")]
use pyo3_stub_gen::derive::*;

use crate::export::{
    HopDict, MaybeEntryDict, PathCountersDict, PathOutcomeDict, ResidualDict, hops, maybe_entries,
};
use crate::find_path::{
    FindPathConfig, HeapScoreMode, PathCounters, PathOutcome, find_path_partial,
    find_path_with_network,
};
use crate::metabolic_network::MetabolicNetwork;
use crate::rules::default_ruleset as default_ruleset_rs;

use super::common::NETWORK_LOCK_POISONED;
use super::graph::PyMetabolicNetwork;
use super::mol::PyForestMol;
use super::plan::PyStepPlan;

fn lock_network(
    arc: &Arc<Mutex<MetabolicNetwork>>,
) -> Result<std::sync::MutexGuard<'_, MetabolicNetwork>, String> {
    arc.lock().map_err(|_| NETWORK_LOCK_POISONED.to_string())
}

/// ``(outcomes, counters, network)`` — shared by ``find_path`` and
/// ``find_path_partial``.
type SearchResult = (Vec<PyPathOutcome>, PyPathCounters, PyMetabolicNetwork);

/// Run the default-ruleset search; `partial` adds the closest-reach flush.
fn run_search(
    py: Python<'_>,
    reactant: &str,
    target: &str,
    config: FindPathConfig,
    network: Option<&Bound<'_, PyMetabolicNetwork>>,
    partial: bool,
) -> PyResult<SearchResult> {
    let rules = default_ruleset_rs();
    let reactant = reactant.to_owned();
    let target = target.to_owned();
    let net_arc = match network {
        Some(n) => n.borrow().inner.clone(),
        None => Arc::new(Mutex::new(MetabolicNetwork::new())),
    };

    let (outcomes, counters) = py
        .detach({
            let net_arc = net_arc.clone();
            move || -> Result<_, String> {
                let mut counters = PathCounters::default();
                let mut guard = lock_network(&net_arc)?;
                let net = Some(&mut *guard);
                let outcomes = if partial {
                    find_path_partial(
                        &reactant,
                        &target,
                        &rules,
                        &mut counters,
                        config,
                        net,
                        |_| true,
                    )
                } else {
                    find_path_with_network(
                        &reactant,
                        &target,
                        &rules,
                        &mut counters,
                        config,
                        net,
                        |_| true,
                    )
                    .and_then(|it| it.collect_all())
                }
                .map_err(|e| e.to_string())?;
                Ok((outcomes, counters))
            }
        })
        .map_err(PyValueError::new_err)?;

    Ok((
        outcomes
            .into_iter()
            .map(|inner| PyPathOutcome { inner })
            .collect(),
        PyPathCounters { inner: counters },
        PyMetabolicNetwork { inner: net_arc },
    ))
}

/// Native chematic ``find_path`` (default ruleset).
///
/// Returns ``(hits, counters, network)``; every hit is exact
/// (``residual_cost == 0``). Search always records hops on a
/// [`MetabolicNetwork`]: pass ``network=`` to extend a live graph, or omit it
/// to get a fresh one back. Default ruleset is QuinoneFormation +
/// EpoxideHydration + Tautomerization + PhaseOne core.
///
/// Releases the GIL for the Rust search (including ``network=``).
#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (
    reactant,
    target,
    *,
    max_paths=1,
    max_nodes=800,
    use_atom_diff=true,
    lazy_closer=false,
    diversity=false,
    drop_skeleton_twins=true,
    score="log-neg-pc",
    timeout=None,
    network=None,
    normalize_tautomer=false,
    invert_target_tautomer=false,
))]
#[allow(clippy::too_many_arguments)]
pub fn find_path(
    py: Python<'_>,
    reactant: &str,
    target: &str,
    max_paths: usize,
    max_nodes: usize,
    use_atom_diff: bool,
    lazy_closer: bool,
    diversity: bool,
    drop_skeleton_twins: bool,
    score: &str,
    timeout: Option<f64>,
    network: Option<&Bound<'_, PyMetabolicNetwork>>,
    normalize_tautomer: bool,
    invert_target_tautomer: bool,
) -> PyResult<SearchResult> {
    let mut config = parse_find_path_config(
        score,
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
        diversity,
        drop_skeleton_twins,
        timeout,
    )?;
    config.normalize_tautomer = normalize_tautomer;
    config.invert_target_tautomer = invert_target_tautomer;
    run_search(py, reactant, target, config, network, false)
}

#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "PathOutcome", unsendable)]
/// One reactant→target reach from ``find_path`` / ``find_path_partial``.
///
/// Exact hits have ``residual_cost == 0``; ``find_path_partial`` may also
/// return closest non-exact reaches (``residual_cost > 0``).
pub struct PyPathOutcome {
    inner: PathOutcome,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
#[pymethods]
impl PyPathOutcome {
    #[getter]
    fn smiles(&self) -> &str {
        &self.inner.smiles
    }

    #[getter]
    fn mol(&self) -> PyForestMol {
        PyForestMol::wrap(self.inner.mol.clone())
    }

    #[getter]
    fn plan(&self, py: Python<'_>) -> PyResult<Py<PyStepPlan>> {
        Py::new(
            py,
            PyStepPlan {
                inner: self.inner.plan.clone(),
            },
        )
    }

    /// Cleavage Maybe bags (same as ``plan.maybe()``).
    fn maybe(&self) -> Vec<MaybeEntryDict> {
        maybe_entries(self.inner.plan.maybe())
    }

    /// Search hops (pattern / site on reactant / product) for pathway displays.
    fn hops(&self) -> Vec<HopDict> {
        hops(&self.inner.steps)
    }

    /// Leftover disagreement with the target; ``0`` for exact hits.
    #[getter]
    fn residual_cost(&self) -> usize {
        self.inner.residual.cost
    }

    /// ``residual_cost == 0``: reached the target.
    #[getter]
    fn is_exact(&self) -> bool {
        self.inner.is_exact()
    }

    /// Full residual (cost, extra atoms, category hints); all-zero when exact.
    fn residual(&self) -> ResidualDict {
        ResidualDict::from(&self.inner.residual)
    }

    fn to_dict(&self) -> PathOutcomeDict {
        PathOutcomeDict::from(&self.inner)
    }

    fn __str__(&self) -> String {
        super::display::format_path_outcome(
            &self.inner.smiles,
            self.inner.residual.cost,
            self.inner.steps.len(),
            self.inner.plan.steps().len(),
            self.inner.plan.n_linearizations(),
        )
    }

    fn __repr__(&self) -> String {
        if self.inner.is_exact() {
            format!("PathOutcome({:?})", self.inner.smiles)
        } else {
            format!(
                "PathOutcome({:?}, residual_cost={})",
                self.inner.smiles, self.inner.residual.cost
            )
        }
    }
}

#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "PathCounters")]
pub struct PyPathCounters {
    inner: PathCounters,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
#[pymethods]
impl PyPathCounters {
    fn to_dict(&self) -> PathCountersDict {
        PathCountersDict::from(&self.inner)
    }

    #[getter]
    fn billed(&self) -> usize {
        self.inner.billed()
    }

    #[getter]
    fn timed_out(&self) -> bool {
        self.inner.timed_out
    }

    #[getter]
    fn diversity_repush(&self) -> usize {
        self.inner.diversity_repush
    }

    fn __str__(&self) -> String {
        super::display::format_path_counters(&self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "PathCounters(billed={}, nodes={}, timed_out={})",
            self.inner.billed(),
            self.inner.nodes,
            self.inner.timed_out
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn parse_find_path_config(
    score: &str,
    max_paths: usize,
    max_nodes: usize,
    use_atom_diff: bool,
    lazy_closer: bool,
    diversity: bool,
    drop_skeleton_twins: bool,
    timeout: Option<f64>,
) -> PyResult<FindPathConfig> {
    let heap_score = HeapScoreMode::from_label(score).ok_or_else(|| {
        PyValueError::new_err(format!(
            "unknown score {score:?}; try log-neg-pc, soft, add-both, …"
        ))
    })?;
    let timeout = match timeout {
        None => None,
        Some(secs) if secs.is_finite() && secs >= 0.0 => {
            Some(std::time::Duration::from_secs_f64(secs))
        }
        Some(secs) => {
            return Err(PyValueError::new_err(format!(
                "timeout must be a non-negative finite number of seconds; got {secs}"
            )));
        }
    };
    Ok(FindPathConfig {
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
        heap_score,
        drop_skeleton_twins,
        diversity,
        timeout,
        ..FindPathConfig::default()
    })
}

/// Like ``find_path``, but also returns the closest non-exact reaches.
///
/// Same return type as ``find_path``: ``(outcomes, counters, network)``.
/// Exact hits come first (``residual_cost == 0``), then — only when they do
/// not fill ``max_paths`` — closest partials by ascending ``residual_cost``;
/// at most ``max_paths`` in total. Search always records hops on a
/// [`MetabolicNetwork`] (pass ``network=`` to extend one).
///
/// Releases the GIL for the Rust search (including ``network=``).
#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "find_path_partial", signature = (
    reactant,
    target,
    *,
    max_paths=1,
    max_nodes=800,
    use_atom_diff=true,
    lazy_closer=false,
    diversity=false,
    drop_skeleton_twins=true,
    score="log-neg-pc",
    timeout=None,
    network=None,
    normalize_tautomer=false,
    invert_target_tautomer=false,
))]
#[allow(clippy::too_many_arguments)]
pub fn find_path_partial_py(
    py: Python<'_>,
    reactant: &str,
    target: &str,
    max_paths: usize,
    max_nodes: usize,
    use_atom_diff: bool,
    lazy_closer: bool,
    diversity: bool,
    drop_skeleton_twins: bool,
    score: &str,
    timeout: Option<f64>,
    network: Option<&Bound<'_, PyMetabolicNetwork>>,
    normalize_tautomer: bool,
    invert_target_tautomer: bool,
) -> PyResult<SearchResult> {
    let mut config = parse_find_path_config(
        score,
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
        diversity,
        drop_skeleton_twins,
        timeout,
    )?;
    config.normalize_tautomer = normalize_tautomer;
    config.invert_target_tautomer = invert_target_tautomer;
    run_search(py, reactant, target, config, network, true)
}
