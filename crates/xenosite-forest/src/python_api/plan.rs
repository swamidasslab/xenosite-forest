//! Elementary plan handle ([`Deps`](crate::canonical_plan::Deps) → Python ``StepPlan``).

use std::collections::BTreeSet;

use pyo3::prelude::*;

use crate::canonical_plan::Deps;

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
