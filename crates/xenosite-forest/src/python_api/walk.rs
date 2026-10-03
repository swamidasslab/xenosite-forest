//! random_path and RandomPathOutcome.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubs")]
use pyo3_stub_gen::derive::*;

use crate::export::{PatternSummaryDict, RandomPathOutcomeDict, RandomPathStepDict};
use crate::pathway::PathwayOptions;
use crate::random_path::{RandomPathOutcome, random_path as random_path_rs, random_path_with};
use crate::rules::phase_one as phase_one_rs;
use crate::ruleset::RuleSet;

use super::rules::PyRuleSet;

#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "RandomPathOutcome")]
pub struct PyRandomPathOutcome {
    inner: RandomPathOutcome,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
#[pymethods]
impl PyRandomPathOutcome {
    #[getter]
    fn smiles(&self) -> &str {
        &self.inner.smiles
    }

    #[getter]
    fn path(&self) -> Vec<String> {
        self.inner.path.clone()
    }

    fn to_dict(&self) -> RandomPathOutcomeDict {
        RandomPathOutcomeDict::from(&self.inner)
    }

    fn steps(&self) -> Vec<RandomPathStepDict> {
        self.inner
            .steps
            .iter()
            .map(RandomPathStepDict::from)
            .collect()
    }

    fn patterns(&self) -> Vec<PatternSummaryDict> {
        self.inner
            .patterns
            .iter()
            .map(PatternSummaryDict::from)
            .collect()
    }

    fn __eq__(&self, other: Bound<'_, PyAny>) -> PyResult<bool> {
        if let Ok(o) = other.extract::<PyRef<Self>>() {
            return Ok(self.inner == o.inner);
        }
        Ok(false)
    }

    fn __str__(&self) -> String {
        let path = &self.inner.path;
        let shown = if path.len() <= 5 {
            path.join(" → ")
        } else {
            format!(
                "{} → … → {} ({} hops)",
                path.first().map(String::as_str).unwrap_or("?"),
                path.last().map(String::as_str).unwrap_or("?"),
                path.len().saturating_sub(1)
            )
        };
        format!("RandomPathOutcome  {}\n  {}", self.inner.smiles, shown)
    }

    fn __repr__(&self) -> String {
        format!("RandomPathOutcome({:?})", self.inner.smiles)
    }
}

/// Seeded random walk: apply up to ``max_steps`` rules. Returns
/// [`RandomPathOutcome`].
///
/// Default ruleset is PhaseOne. Pass a ``RuleSet`` to override.
/// ``skip_multicomponent`` / ``skip_seen`` map to [`PathwayOptions`] (off by
/// default; same knobs for StepSequence / PathOutcome ``apply``).
///
/// Releases the GIL for the Rust walk.
#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(signature = (
    reactant,
    seed,
    *,
    max_steps=1,
    ruleset=None,
    skip_multicomponent=false,
    skip_seen=false,
))]
pub fn random_path(
    py: Python<'_>,
    reactant: &str,
    seed: u64,
    max_steps: usize,
    ruleset: Option<&Bound<'_, PyRuleSet>>,
    skip_multicomponent: bool,
    skip_seen: bool,
) -> PyResult<PyRandomPathOutcome> {
    let rules: RuleSet = match ruleset {
        Some(rs) => rs.borrow().inner.clone(),
        None => phase_one_rs(),
    };
    let options = PathwayOptions {
        skip_multicomponent,
        skip_seen,
    };
    let reactant = reactant.to_owned();
    let outcome = py
        .detach(move || {
            if options == PathwayOptions::default() {
                random_path_rs(&reactant, seed, &rules, max_steps)
            } else {
                random_path_with(&reactant, seed, &rules, max_steps, options)
            }
            .map_err(|e| e.to_string())
        })
        .map_err(PyValueError::new_err)?;
    Ok(PyRandomPathOutcome { inner: outcome })
}
