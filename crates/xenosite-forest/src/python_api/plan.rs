//! Elementary plan handle ([`Deps`](crate::canonical_plan::Deps) → Python ``StepPlan``).

use std::collections::BTreeSet;

use pyo3::prelude::*;

use crate::canonical_plan::{Deps, Maybe};

/// Serialize [`Maybe`] bags for Python (``site`` / ``side`` / ``opens`` / ``span_sites``).
pub(crate) fn maybe_to_py(py: Python<'_>, maybe: &Maybe) -> PyResult<Py<PyAny>> {
    let rows: Vec<Py<PyAny>> = maybe
        .entries
        .iter()
        .map(|e| {
            let d = pyo3::types::PyDict::new(py);
            let site: Vec<usize> = e.site.iter().copied().collect();
            d.set_item("site", site)?;
            d.set_item("side", &e.side)?;
            let opens: Vec<Vec<usize>> = e
                .opens
                .iter()
                .map(|o| o.iter().copied().collect())
                .collect();
            d.set_item("opens", opens)?;
            let span: Vec<Vec<usize>> = e
                .span_sites()
                .into_iter()
                .map(|s| s.iter().copied().collect())
                .collect();
            d.set_item("span_sites", span)?;
            Ok(d.unbind().into_any())
        })
        .collect::<PyResult<_>>()?;
    Ok(pyo3::types::PyList::new(py, rows)?.unbind().into_any())
}

/// Elementary plan with precedes (Python ``StepPlan``).
#[pyclass(name = "StepPlan", unsendable)]
pub struct PyStepPlan {
    pub(crate) inner: Deps,
}

#[pymethods]
impl PyStepPlan {
    fn __len__(&self) -> usize {
        self.inner.steps().len()
    }

    fn n_linearizations(&self) -> usize {
        self.inner.n_linearizations()
    }

    /// Enumerate topological sorts under precedes (tiny plans).
    ///
    /// Each row is a list of ``{"rule", "site"}`` like plan steps. Caps at
    /// 64 rows so accidental factorial explosions stay bounded.
    fn linearizations(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        const MAX: usize = 64;
        let lins = self.inner.linearizations();
        let rows: Vec<Py<PyAny>> = lins
            .iter()
            .take(MAX)
            .map(|lin| {
                let steps: Vec<Py<PyAny>> = lin
                    .steps
                    .iter()
                    .map(|step| {
                        let view = crate::export::PlanStepView::from(step);
                        Ok(pythonize::pythonize(py, &view)?.unbind().into_any())
                    })
                    .collect::<PyResult<_>>()?;
                Ok(pyo3::types::PyList::new(py, steps)?.unbind().into_any())
            })
            .collect::<PyResult<_>>()?;
        Ok(pyo3::types::PyList::new(py, rows)?.unbind().into_any())
    }

    fn allows(&self, site: Option<Vec<usize>>, side: Option<&str>) -> bool {
        let set: Option<BTreeSet<usize>> = site.map(|v| v.into_iter().collect());
        self.inner.allows(set.as_ref(), side)
    }

    /// Cleavage-side bags on this plan (`site` / `side` / `opens`).
    fn maybe(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        maybe_to_py(py, self.inner.maybe())
    }

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let steps: Vec<Py<PyAny>> = self
            .inner
            .steps()
            .iter()
            .map(|step| {
                let view = crate::export::PlanStepView::from(step);
                Ok(pythonize::pythonize(py, &view)?.unbind().into_any())
            })
            .collect::<PyResult<_>>()?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("steps", steps)?;
        d.set_item("n_linearizations", self.inner.n_linearizations())?;
        let precedes: Vec<(usize, usize)> = self.inner.precedes().to_vec();
        d.set_item("precedes", precedes)?;
        d.set_item("maybe", maybe_to_py(py, self.inner.maybe())?)?;
        Ok(d.unbind().into_any())
    }

    fn __str__(&self) -> String {
        super::display::format_step_plan(&self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "StepPlan({} steps, ~{} lin)",
            self.inner.steps().len(),
            self.inner.n_linearizations()
        )
    }
}
