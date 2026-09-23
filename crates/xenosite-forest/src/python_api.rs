//! PyO3 class wrap of [`crate::xf::ForestMol`].
//!
//! `#[pyclass]` stores the Rust struct as the Python instance payload. One
//! Python object ↔ one `ForestMol`. Getters are methods on that payload;
//! they fill `_forest["cache"]` the same way native `xf` does.
//!
//! Native-only: CPython C-API. Not compiled for `wasm32-unknown-unknown`.

use std::collections::HashMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::forest::Formula;
use crate::xf::ForestMol;

fn py_err(err: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// Python-visible formula. Nested `#[pyclass]` wrap of [`Formula`].
#[pyclass(name = "Formula", frozen)]
#[derive(Clone)]
pub struct PyFormula {
    #[pyo3(get)]
    pub counts: HashMap<String, i32>,
    #[pyo3(get)]
    pub charge: i32,
}

impl From<&Formula> for PyFormula {
    fn from(formula: &Formula) -> Self {
        Self {
            counts: formula.counts.clone().into_iter().collect(),
            charge: formula.charge,
        }
    }
}

/// Python class wrapping [`ForestMol`].
///
/// `unsendable`: the payload holds `RefCell` cache state and a chematic mol.
/// `#[new]` is `__init__`. Getters become Python properties.
#[pyclass(name = "ForestMol", unsendable)]
pub struct PyForestMol {
    inner: ForestMol,
    csmi: Option<Py<PyString>>,
    formula: Option<Py<PyFormula>>,
}

impl PyForestMol {
    fn wrap(inner: ForestMol) -> Self {
        Self {
            inner,
            csmi: None,
            formula: None,
        }
    }
}

#[pymethods]
impl PyForestMol {
    #[new]
    fn new(smiles: &str) -> PyResult<Self> {
        Ok(Self::wrap(ForestMol::parse(smiles).map_err(py_err)?))
    }

    #[getter]
    fn has_forest(&self) -> bool {
        self.inner.xf().has_forest()
    }

    /// Cached canonical SMILES. Interned so `mol.csmi is mol.csmi`.
    #[getter]
    fn csmi(&mut self, py: Python<'_>) -> Py<PyString> {
        if let Some(held) = &self.csmi {
            return held.clone_ref(py);
        }
        let s = self.inner.xf().csmi();
        let interned = PyString::intern(py, s.as_ref()).unbind();
        self.csmi = Some(interned.clone_ref(py));
        interned
    }

    #[getter]
    fn formula(&mut self, py: Python<'_>) -> PyResult<Py<PyFormula>> {
        if let Some(held) = &self.formula {
            return Ok(held.clone_ref(py));
        }
        let obj = Py::new(py, PyFormula::from(self.inner.xf().formula().as_ref()))?;
        self.formula = Some(obj.clone_ref(py));
        Ok(obj)
    }

    fn clear_structure(&mut self) {
        self.inner.xf().clear_structure();
        self.csmi = None;
        self.formula = None;
    }

    fn copy(&self, py: Python<'_>) -> Self {
        let mut out = Self::wrap(self.inner.copy_mol());
        out.csmi = self.csmi.as_ref().map(|s| s.clone_ref(py));
        out.formula = self.formula.as_ref().map(|f| f.clone_ref(py));
        out
    }

    fn rw_copy(&self) -> Self {
        Self::wrap(self.inner.rw_copy())
    }

    fn wipe_forest(&mut self) {
        self.inner.wipe_forest();
        self.csmi = None;
        self.formula = None;
    }

    fn smarts_matches(&self, smarts: &str) -> PyResult<Vec<HashMap<u16, usize>>> {
        let hits = self.inner.xf().smarts_matches(smarts).map_err(py_err)?;
        Ok(hits
            .iter()
            .map(|mapped| mapped.iter().map(|(&k, &v)| (k, v)).collect())
            .collect())
    }

    fn __repr__(&self) -> String {
        if self.inner.has_forest() {
            format!("ForestMol({:?})", self.inner.xf().csmi())
        } else {
            "ForestMol(<no forest>)".to_string()
        }
    }
}

#[pymodule]
fn xenosite_forest(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForestMol>()?;
    m.add_class::<PyFormula>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pyclass_wraps_forest_mol_and_keeps_csmi_identity() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "xenosite_forest").unwrap();
            module.add_class::<PyForestMol>().unwrap();
            module.add_class::<PyFormula>().unwrap();
            let class = module.getattr("ForestMol").unwrap();
            let mol = class.call1(("CCO",)).unwrap();
            assert_eq!(mol.get_type().name().unwrap(), "ForestMol");
            assert!(
                !mol.getattr("has_forest")
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );
            let csmi = mol.getattr("csmi").unwrap();
            assert!(
                mol.getattr("has_forest")
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );
            let again = mol.getattr("csmi").unwrap();
            assert!(
                csmi.is(&again),
                "cached csmi must be the same Python object"
            );
            let formula = mol.getattr("formula").unwrap();
            assert_eq!(formula.get_type().name().unwrap(), "Formula");
            let again_formula = mol.getattr("formula").unwrap();
            assert!(
                formula.is(&again_formula),
                "cached formula must be the same Python object"
            );
            let charge: i32 = formula.getattr("charge").unwrap().extract().unwrap();
            assert_eq!(charge, 0);
            let copied = mol.call_method0("copy").unwrap();
            assert!(
                copied
                    .getattr("has_forest")
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );
            let rw = mol.call_method0("rw_copy").unwrap();
            assert!(!rw.getattr("has_forest").unwrap().extract::<bool>().unwrap());
            mol.call_method0("wipe_forest").unwrap();
            assert!(
                !mol.getattr("has_forest")
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );
        });
    }
}
