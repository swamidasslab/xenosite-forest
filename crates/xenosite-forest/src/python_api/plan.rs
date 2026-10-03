//! Elementary plan handle ([`Deps`](crate::canonical_plan::Deps) → Python ``StepPlan``).

use std::collections::BTreeSet;

use pyo3::prelude::*;
#[cfg(feature = "stubs")]
use pyo3_stub_gen::derive::*;

use crate::canonical_plan::Deps;
use crate::export::{MaybeEntryDict, PlanStepDict, StepPlanDict, maybe_entries, plan_steps};

/// Elementary plan with precedes (Python ``StepPlan``).
#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "StepPlan", unsendable)]
pub struct PyStepPlan {
    pub(crate) inner: Deps,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
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
    fn linearizations(&self) -> Vec<Vec<PlanStepDict>> {
        const MAX: usize = 64;
        self.inner
            .linearizations()
            .iter()
            .take(MAX)
            .map(|lin| plan_steps(&lin.steps))
            .collect()
    }

    fn allows(&self, site: Option<Vec<usize>>, side: Option<&str>) -> bool {
        let set: Option<BTreeSet<usize>> = site.map(|v| v.into_iter().collect());
        self.inner.allows(set.as_ref(), side)
    }

    /// Cleavage-side bags on this plan (`site` / `side` / `opens`).
    fn maybe(&self) -> Vec<MaybeEntryDict> {
        maybe_entries(self.inner.maybe())
    }

    fn to_dict(&self) -> StepPlanDict {
        StepPlanDict {
            steps: plan_steps(self.inner.steps()),
            n_linearizations: self.inner.n_linearizations(),
            precedes: self.inner.precedes().to_vec(),
            maybe: maybe_entries(self.inner.maybe()),
        }
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
