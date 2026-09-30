//! ForestMol / Formula pyclasses and tautomer normalize.

use std::collections::HashMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::forest::Formula;
use crate::forest_mol::ForestMol;

use super::common::py_err;

/// Rust ``ForestMol`` pyclass or the public ``xenosite.forest.mol.ForestMol`` wrapper (``_rust``).
pub fn py_forest_mol_ref<'py>(mol: &Bound<'py, PyAny>) -> PyResult<PyRef<'py, PyForestMol>> {
    if let Ok(r) = mol.extract::<PyRef<'_, PyForestMol>>() {
        return Ok(r);
    }
    if mol.hasattr("_rust")? {
        let inner = mol.getattr("_rust")?;
        return inner.extract::<PyRef<'py, PyForestMol>>();
    }
    Err(PyValueError::new_err(
        "expected ForestMol (extension or xenosite.forest.mol wrapper)",
    ))
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
/// `unsendable`: pyclass stays on one Python thread; payload `ForestMol` is Send
/// (Arc caches) so pure-Rust doors can release the GIL.
/// `#[new]` is `__init__`. Getters become Python properties.
#[pyclass(name = "ForestMol", unsendable)]
pub struct PyForestMol {
    pub(crate) inner: ForestMol,
    csmi: Option<Py<PyString>>,
    formula: Option<Py<PyFormula>>,
}

impl PyForestMol {
    pub(crate) fn wrap(inner: ForestMol) -> Self {
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

    /// Cached canonical SMILES. Interned so `mol.csmi is mol.csmi`.
    #[getter]
    fn csmi(&mut self, py: Python<'_>) -> Py<PyString> {
        if let Some(held) = &self.csmi {
            return held.clone_ref(py);
        }
        let s = self.inner.csmi();
        let interned = PyString::intern(py, s.as_ref()).unbind();
        self.csmi = Some(interned.clone_ref(py));
        interned
    }

    /// Fail-closed dedup key (Chematic `canonical_smiles_stable_key`).
    ///
    /// **Can return `None`.** Do not fall back to [`Self::csmi`] for HashSet /
    /// yield identity — skip CSMI dedup instead.
    #[getter]
    fn stable_csmi_key(&self, py: Python<'_>) -> Option<Py<PyString>> {
        self.inner
            .stable_csmi_key()
            .map(|s| PyString::intern(py, s.as_ref()).unbind())
    }

    #[getter]
    fn formula(&mut self, py: Python<'_>) -> PyResult<Py<PyFormula>> {
        if let Some(held) = &self.formula {
            return Ok(held.clone_ref(py));
        }
        let obj = Py::new(py, PyFormula::from(self.inner.formula().as_ref()))?;
        self.formula = Some(obj.clone_ref(py));
        Ok(obj)
    }

    fn clear_structure(&mut self) {
        self.inner.clear_structure();
        self.csmi = None;
        self.formula = None;
    }

    fn copy(&self, py: Python<'_>) -> Self {
        let mut out = Self::wrap(self.inner.copy_mol());
        out.csmi = self.csmi.as_ref().map(|s| s.clone_ref(py));
        out.formula = self.formula.as_ref().map(|f| f.clone_ref(py));
        out
    }

    fn edit_copy(&self) -> Self {
        Self::wrap(self.inner.edit_copy())
    }

    /// Chematic tautomer pick adopted with Forest tracing.
    ///
    /// Returns ``(ForestMol, changed)``.
    fn normalize_tautomer(&self) -> PyResult<(Self, bool)> {
        let out = self.inner.normalize_tautomer().map_err(py_err)?;
        Ok((Self::wrap(out.mol), out.changed))
    }

    fn smarts_matches(&self, smarts: &str) -> PyResult<Vec<HashMap<u16, usize>>> {
        let hits = self.inner.smarts_matches(smarts).map_err(py_err)?;
        Ok(hits
            .iter()
            .map(|mapped| mapped.iter().map(|(&k, &v)| (k, v)).collect())
            .collect())
    }

    fn __repr__(&self) -> String {
        format!("ForestMol({:?})", self.inner.csmi())
    }
}

/// Chematic tautomer pick adopted as a tagged [`ForestMol`].
///
/// Returns ``(mol, changed)`` where ``mol`` is a :class:`ForestMol` and
/// ``changed`` is whether the form differed from the input.
#[pyfunction]
pub fn normalize_tautomer(smiles: &str) -> PyResult<(PyForestMol, bool)> {
    let out = crate::normalize_tautomer(smiles).map_err(py_err)?;
    Ok((PyForestMol::wrap(out.mol), out.changed))
}
