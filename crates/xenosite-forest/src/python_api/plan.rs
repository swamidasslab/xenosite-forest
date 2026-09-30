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
        Ok(d.unbind().into_any())
    }
}
