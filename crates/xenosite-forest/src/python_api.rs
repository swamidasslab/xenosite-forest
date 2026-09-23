//! PyO3 class wrap of [`crate::xf::ForestMol`].
//!
//! `#[pyclass]` stores the Rust struct as the Python instance payload. One
//! Python object ↔ one `ForestMol`. Getters are methods on that payload;
//! they fill `_forest["cache"]` the same way native `xf` does.
//!
//! Native-only: CPython C-API. Not compiled for `wasm32-unknown-unknown`.

use std::cell::RefCell;
use std::collections::HashMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::forest::Formula;
use crate::mol::Molecule;
use crate::pattern::{Edit, Effect, PatternInfo, SiteInfo};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};
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

fn parse_edit(edit: &str) -> Edit {
    if edit.eq_ignore_ascii_case("hydroxyl") {
        Edit::Hydroxyl
    } else {
        Edit::Smirks(edit.to_string())
    }
}

fn edit_label(edit: &Edit) -> String {
    match edit {
        Edit::Hydroxyl => "hydroxyl".into(),
        Edit::Smirks(smirks) => smirks.clone(),
    }
}

/// Python wrap of [`PatternInfo`]. Frozen data; `RuleSet` clones it in.
#[pyclass(name = "PatternInfo", frozen)]
#[derive(Clone)]
pub struct PyPatternInfo {
    inner: PatternInfo,
}

#[pymethods]
impl PyPatternInfo {
    #[new]
    #[pyo3(signature = (name, smarts, edit, adds=None, removes=None, cleaves=false, methide=false))]
    fn new(
        name: String,
        smarts: String,
        edit: String,
        adds: Option<String>,
        removes: Option<String>,
        cleaves: bool,
        methide: bool,
    ) -> Self {
        Self {
            inner: PatternInfo::new(
                name,
                smarts,
                parse_edit(&edit),
                Effect {
                    adds,
                    removes,
                    cleaves,
                    methide,
                },
            ),
        }
    }

    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn smarts(&self) -> &str {
        &self.inner.smarts
    }

    #[getter]
    fn edit(&self) -> String {
        edit_label(&self.inner.edit)
    }

    #[getter]
    fn adds(&self) -> Option<String> {
        self.inner.effect.adds.clone()
    }

    #[getter]
    fn removes(&self) -> Option<String> {
        self.inner.effect.removes.clone()
    }

    #[getter]
    fn cleaves(&self) -> bool {
        self.inner.effect.cleaves
    }

    #[getter]
    fn methide(&self) -> bool {
        self.inner.effect.methide
    }

    fn __repr__(&self) -> String {
        format!(
            "PatternInfo({:?}, {:?})",
            self.inner.name, self.inner.smarts
        )
    }
}

/// Python wrap of [`RuleSet`]. Patterns are copied into Rust at construction.
#[pyclass(name = "RuleSet", unsendable)]
pub struct PyRuleSet {
    inner: RuleSet,
}

#[pymethods]
impl PyRuleSet {
    #[new]
    #[pyo3(signature = (patterns=None, name=None))]
    fn new(patterns: Option<Vec<PyPatternInfo>>, name: Option<String>) -> Self {
        Self {
            inner: RuleSet::new(
                name,
                patterns
                    .unwrap_or_default()
                    .into_iter()
                    .map(|pattern| pattern.inner),
            ),
        }
    }

    #[staticmethod]
    fn hydroxylation() -> Self {
        Self {
            inner: crate::hydroxylation::hydroxylation(),
        }
    }

    #[staticmethod]
    fn o_dealkylation() -> Self {
        Self {
            inner: crate::ruleset::o_dealkylation(),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (sets, name=None))]
    fn compose(sets: Vec<PyRef<'_, PyRuleSet>>, name: Option<String>) -> Self {
        Self {
            inner: RuleSet::compose(name, sets.into_iter().map(|set| set.inner.clone())),
        }
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.patterns().len()
    }

    fn patterns(&self) -> Vec<PyPatternInfo> {
        self.inner
            .patterns()
            .iter()
            .cloned()
            .map(|inner| PyPatternInfo { inner })
            .collect()
    }

    /// Run the owned patterns. `filter_rules` / `filter_sites` are optional Python
    /// callables; omit them and Rust `accept_all_*` runs with no GIL per site.
    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyForestMol>,
        filter_rules: Option<Bound<'_, PyAny>>,
        filter_sites: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Vec<(String, usize, Vec<String>)>> {
        let chemistry = mol.borrow().inner.mol().clone();
        let set = slf.borrow().inner.clone();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            set.metabolize(&chemistry, accept_all_rules, accept_all_sites, true)
                .map_err(py_err)?
        } else {
            metabolize_with_python(slf, mol, &set, &chemistry, filter_rules, filter_sites)?
        };
        Ok(emissions
            .into_iter()
            .map(|e| (e.pattern_name, e.site, e.products))
            .collect())
    }

    fn __repr__(&self) -> String {
        match &self.inner.name {
            Some(name) => format!(
                "RuleSet({name:?}, {} patterns)",
                self.inner.patterns().len()
            ),
            None => format!("RuleSet({} patterns)", self.inner.patterns().len()),
        }
    }
}

fn metabolize_with_python(
    slf: &Bound<'_, PyRuleSet>,
    mol: &Bound<'_, PyForestMol>,
    set: &RuleSet,
    chemistry: &Molecule,
    filter_rules: Option<Bound<'_, PyAny>>,
    filter_sites: Option<Bound<'_, PyAny>>,
) -> PyResult<Vec<crate::pattern::Emission>> {
    let py_rules = filter_rules.map(|cb| cb.unbind());
    let py_sites = filter_sites.map(|cb| cb.unbind());
    let py_mol = mol.clone().unbind();
    let py_set = slf.clone().unbind();
    let err: RefCell<Option<PyErr>> = RefCell::new(None);
    let take_bool = |result: PyResult<bool>, slot: &RefCell<Option<PyErr>>| match result {
        Ok(keep) => keep,
        Err(e) => {
            *slot.borrow_mut() = Some(e);
            false
        }
    };
    let rules = |_: &Molecule, _: &RuleSet, pattern: &PatternInfo| {
        if err.borrow().is_some() {
            return false;
        }
        let Some(cb) = &py_rules else {
            return true;
        };
        take_bool(
            Python::attach(|py| {
                let info = Py::new(
                    py,
                    PyPatternInfo {
                        inner: pattern.clone(),
                    },
                )?;
                cb.bind(py)
                    .call1((py_mol.bind(py), py_set.bind(py), info))?
                    .extract::<bool>()
            }),
            &err,
        )
    };
    let sites = |_: &Molecule, site: usize, info: &SiteInfo| {
        if err.borrow().is_some() {
            return false;
        }
        let Some(cb) = &py_sites else {
            return true;
        };
        take_bool(
            Python::attach(|py| {
                let bag = Py::new(
                    py,
                    PyPatternInfo {
                        inner: info.pattern.clone(),
                    },
                )?;
                cb.bind(py)
                    .call1((py_mol.bind(py), site, bag))?
                    .extract::<bool>()
            }),
            &err,
        )
    };
    let emissions = set
        .metabolize(chemistry, rules, sites, true)
        .map_err(py_err)?;
    if let Some(e) = err.into_inner() {
        return Err(e);
    }
    Ok(emissions)
}

#[pymodule]
fn xenosite_forest(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForestMol>()?;
    m.add_class::<PyFormula>()?;
    m.add_class::<PyPatternInfo>()?;
    m.add_class::<PyRuleSet>()?;
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

    #[test]
    fn python_composes_ruleset_once_then_metabolize_is_a_handle() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "xenosite_forest").unwrap();
            module.add_class::<PyForestMol>().unwrap();
            module.add_class::<PyPatternInfo>().unwrap();
            module.add_class::<PyRuleSet>().unwrap();
            let pattern = module.getattr("PatternInfo").unwrap();
            let ruleset = module.getattr("RuleSet").unwrap();
            let mol_cls = module.getattr("ForestMol").unwrap();
            let mol = mol_cls.call1(("c1ccccc1",)).unwrap();
            let h = pattern
                .call1(("h", "[#6h1:1]", "hydroxyl", "O", "H"))
                .unwrap();
            let h2 = pattern
                .call1(("h2", "[#6h2,#6h3:1]", "hydroxyl", "O", "H"))
                .unwrap();
            let rs = ruleset.call1((vec![h, h2], "Hydroxylation")).unwrap();
            assert_eq!(
                rs.call_method0("__len__")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                2
            );
            let products = rs.call_method1("metabolize", (&mol,)).unwrap();
            let products: Vec<(String, usize, Vec<String>)> = products.extract().unwrap();
            assert_eq!(products.len(), 1);
            assert_eq!(products[0].0, "h");
            let filt = py
                .eval(c"lambda m, rule, p: p.name == 'h2'", None, None)
                .unwrap();
            let kept: Vec<(String, usize, Vec<String>)> = rs
                .call_method1("metabolize", (&mol, filt))
                .unwrap()
                .extract()
                .unwrap();
            assert!(kept.is_empty());
            let built_in = ruleset.call_method0("hydroxylation").unwrap();
            assert_eq!(
                built_in
                    .call_method0("__len__")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                2
            );
            let dealkyl = ruleset.call_method0("o_dealkylation").unwrap();
            let composed = ruleset
                .call_method1("compose", (vec![built_in, dealkyl], "probe"))
                .unwrap();
            assert_eq!(
                composed
                    .call_method0("__len__")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                3
            );
        });
    }
}
