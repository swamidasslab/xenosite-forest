//! find_path / find_path_partial and path outcome pyclasses.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::export::{PathCountersView, PathOutcomeView};
use crate::find_path::{
    FindPathConfig, FindPathPartialResult, HeapScoreMode, PartialOutcome, PathCounters,
    PathOutcome, find_path_partial, find_path_with, find_path_with_network,
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

/// Native chematic ``find_path`` (default ruleset). Returns ``(hits, counters)``.
///
/// Default ruleset is QuinoneFormation + EpoxideHydration + Tautomerization +
/// PhaseOne core. Each hit is ``{"smiles": str, "steps": [{"rule": str,
/// "site": [...]}]``. Counters is a plain dict of the billed fields. Separate
/// from the Python RDKit ``xenosite.forest.find_path`` walk.
///
/// Releases the GIL for the Rust search (including ``network=``).
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
) -> PyResult<(Vec<PyPathOutcome>, PyPathCounters)> {
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
    let rules = default_ruleset_rs();
    let reactant = reactant.to_owned();
    let target = target.to_owned();
    let net_arc = network.map(|n| n.borrow().inner.clone());

    let (hits, counters) = py
        .detach(move || -> Result<_, String> {
            let mut counters = PathCounters::default();
            let hits = match &net_arc {
                Some(arc) => {
                    let mut guard = lock_network(arc)?;
                    find_path_with_network(
                        &reactant,
                        &target,
                        &rules,
                        &mut counters,
                        config,
                        Some(&mut *guard),
                        |_| true,
                    )
                    .and_then(|it| it.collect_all())
                    .map_err(|e| e.to_string())?
                }
                None => find_path_with(
                    &reactant,
                    &target,
                    &rules,
                    &mut counters,
                    config,
                    |_| true,
                )
                .and_then(|it| it.collect_all())
                .map_err(|e| e.to_string())?,
            };
            Ok((hits, counters))
        })
        .map_err(PyValueError::new_err)?;

    Ok((
        hits.into_iter()
            .map(|h| PyPathOutcome { inner: h })
            .collect(),
        PyPathCounters { inner: counters },
    ))
}

#[pyclass(name = "PathOutcome", unsendable)]
pub struct PyPathOutcome {
    inner: PathOutcome,
}

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

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let view = PathOutcomeView::from(&self.inner);
        Ok(pythonize::pythonize(py, &view)?.unbind().into_any())
    }
}

#[pyclass(name = "PartialOutcome", unsendable)]
pub struct PyPartialPathOutcome {
    inner: PartialOutcome,
}

#[pymethods]
impl PyPartialPathOutcome {
    #[getter]
    fn smiles(&self) -> &str {
        &self.inner.smiles
    }

    #[getter]
    fn mol(&self) -> PyForestMol {
        PyForestMol::wrap(self.inner.mol.clone())
    }

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let hit = PathOutcome {
            steps: self.inner.steps.clone(),
            plan: self.inner.plan.clone(),
            mol: self.inner.mol.clone(),
            smiles: self.inner.smiles.clone(),
        };
        let view = PathOutcomeView::from(&hit);
        let d = pyo3::types::PyDict::new(py);
        d.set_item("smiles", view.smiles)?;
        d.set_item("steps", pythonize::pythonize(py, &view.steps)?)?;
        let residual = pyo3::types::PyDict::new(py);
        residual.set_item("cost", self.inner.residual.cost)?;
        residual.set_item("n_extra", self.inner.residual.n_extra)?;
        residual.set_item("categories", self.inner.residual.categories.clone())?;
        residual.set_item("unresolvable", self.inner.residual.unresolvable)?;
        d.set_item("residual", residual)?;
        Ok(d.unbind().into_any())
    }
}

#[pyclass(name = "PathCounters")]
pub struct PyPathCounters {
    inner: PathCounters,
}

#[pymethods]
impl PyPathCounters {
    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let view = PathCountersView::from(&self.inner);
        Ok(pythonize::pythonize(py, &view)?.unbind().into_any())
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

/// Like ``find_path``, but also returns closest unreachable reaches.
///
/// Returns ``(exact_hits, partials, counters)``. Each partial is
/// ``{"smiles", "steps", "residual": {"cost", "categories", ...}}``.
/// Pass ``network=`` to record hops on a [`MetabolicNetwork`].
///
/// Releases the GIL for the Rust search (including ``network=``).
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
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
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
) -> PyResult<(Vec<PyPathOutcome>, Vec<PyPartialPathOutcome>, PyPathCounters)> {
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
    let rules = default_ruleset_rs();
    let reactant = reactant.to_owned();
    let target = target.to_owned();
    let net_arc = network.map(|n| n.borrow().inner.clone());

    let (exact, partials, counters) = py
        .detach(move || -> Result<_, String> {
            let mut counters = PathCounters::default();
            let result = match &net_arc {
                Some(arc) => {
                    let mut guard = lock_network(arc)?;
                    find_path_partial(
                        &reactant,
                        &target,
                        &rules,
                        &mut counters,
                        config,
                        Some(&mut *guard),
                        |_| true,
                    )
                    .map_err(|e| e.to_string())?
                }
                None => find_path_partial(
                    &reactant,
                    &target,
                    &rules,
                    &mut counters,
                    config,
                    None,
                    |_| true,
                )
                .map_err(|e| e.to_string())?,
            };
            let FindPathPartialResult { exact, partials } = result;
            Ok((exact, partials, counters))
        })
        .map_err(PyValueError::new_err)?;

    Ok((
        exact
            .into_iter()
            .map(|h| PyPathOutcome { inner: h })
            .collect(),
        partials
            .into_iter()
            .map(|p| PyPartialPathOutcome { inner: p })
            .collect(),
        PyPathCounters { inner: counters },
    ))
}
