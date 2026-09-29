//! PyO3 class wrap of [`crate::forest_mol::ForestMol`].
//!
//! This is the **product** Python door (`xenosite_forest` / `xenosite.forest`).
//! New public APIs are added here — not in `xenosite.forest.native` (frozen
//! RDKit reference; see `docs/forest/NATIVE.md`).
//!
//! `#[pyclass]` stores the Rust struct as the Python instance payload. One
//! Python object ↔ one `ForestMol`. Getters are methods on that payload.
//!
//! Native-only: CPython C-API. Not compiled for `wasm32-unknown-unknown`.

use std::cell::RefCell;
use std::collections::HashMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::find_path::{
    FindPathConfig, FindPathPartialResult, HeapScoreMode, PathCounters, find_path_partial,
    find_path_with, find_path_with_network,
};
use crate::forest::Formula;
use crate::forest_mol::ForestMol;
use crate::metabolic_network::MetabolicNetwork;
use crate::mol::Molecule;
use crate::pathway::PathwayOptions;
use crate::pattern::{Edit, Effect, PatternInfo, SiteInfo};
use crate::random_path::{random_path as random_path_rs, random_path_with};
use crate::rules::{
    catalog_names, dealkylation as dealkylation_rs, default_ruleset as default_ruleset_rs,
    dehydrogenation as dehydrogenation_rs, epoxidation as epoxidation_rs,
    epoxide_opening as epoxide_opening_rs, hydrolysis as hydrolysis_rs,
    hydroxylation as hydroxylation_rs, leaf_rule as leaf_rule_rs,
    n_dealkylation as n_dealkylation_rs, phase_one as phase_one_rs,
    quinone_formation as quinone_formation_rs,
};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};

fn py_err(err: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// One `RuleSet.metabolize` row: `(pattern_name, site, products, rule_path)`.
type MetabolizeRow = (String, usize, Vec<String>, Vec<Option<String>>);

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
        Edit::PairEndpoint(name) => format!("pair:{name}"),
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
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
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

/// Python wrap of [`crate::bound_pattern::BoundPattern`].
#[pyclass(name = "BoundPattern", unsendable)]
#[derive(Clone)]
pub struct PyBoundPattern {
    inner: crate::bound_pattern::BoundPattern,
}

#[pymethods]
impl PyBoundPattern {
    #[getter]
    fn rule(&self) -> PyRuleSet {
        PyRuleSet {
            inner: self.inner.rule().clone(),
        }
    }

    #[getter]
    fn pattern(&self) -> PyPatternInfo {
        PyPatternInfo {
            inner: self.inner.pattern().clone(),
        }
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn rule_name(&self) -> Option<String> {
        self.inner.rule_name().map(str::to_string)
    }

    #[getter]
    fn curie(&self) -> String {
        self.inner.curie()
    }

    #[getter]
    fn iri(&self) -> String {
        self.inner.iri()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __contains__(&self, key: &str) -> bool {
        self.inner.contains_name(key)
    }

    fn __getitem__(&self, key: Bound<'_, PyAny>) -> PyResult<PyBoundPattern> {
        if let Ok(index) = key.extract::<isize>() {
            let index = if index < 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                    "BoundPattern index out of range",
                ));
            } else {
                index as usize
            };
            return self
                .inner
                .get(index)
                .map(|inner| PyBoundPattern { inner })
                .ok_or_else(|| {
                    PyErr::new::<pyo3::exceptions::PyIndexError, _>("BoundPattern index out of range")
                });
        }
        let name: String = key.extract()?;
        self.inner
            .get_str(&name)
            .map(|inner| PyBoundPattern { inner })
            .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>(name))
    }

    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyForestMol>,
        filter_rules: Option<Bound<'_, PyAny>>,
        filter_sites: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Vec<MetabolizeRow>> {
        let forest = mol.borrow().inner.copy_mol();
        let bp = slf.borrow().inner.clone();
        let set = bp.rule().clone();
        let pattern_name = bp.name().to_string();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            bp.metabolize_default(&forest, true)
                .collect::<Result<Vec<_>, _>>()
                .map_err(py_err)?
        } else {
            // Compose BoundPattern name filter with optional Python callables.
            let name = pattern_name.clone();
            metabolize_with_python_bound(mol, &set, &forest, &name, filter_rules, filter_sites)?
        };
        let _ = pattern_name;
        Ok(emissions
            .into_iter()
            .map(|e| {
                let products = e.product_csmis();
                (e.pattern_name, e.site, products, e.rule_path)
            })
            .collect())
    }

    fn __repr__(&self) -> String {
        format!("BoundPattern({:?})", self.inner.curie())
    }
}

fn metabolize_with_python_bound(
    mol: &Bound<'_, PyForestMol>,
    set: &RuleSet,
    forest: &crate::ForestMol,
    pattern_name: &str,
    filter_rules: Option<Bound<'_, PyAny>>,
    filter_sites: Option<Bound<'_, PyAny>>,
) -> PyResult<Vec<crate::pattern::Emission>> {
    let pname = pattern_name.to_string();
    let py_rules = filter_rules.map(|cb| cb.unbind());
    let py_sites = filter_sites.map(|cb| cb.unbind());
    let py_mol = mol.clone().unbind();
    let err: RefCell<Option<PyErr>> = RefCell::new(None);
    let take_bool = |result: PyResult<bool>, slot: &RefCell<Option<PyErr>>| match result {
        Ok(keep) => keep,
        Err(e) => {
            *slot.borrow_mut() = Some(e);
            false
        }
    };
    let rules = |_: &Molecule, leaf: &RuleSet, pattern: &PatternInfo| {
        if pattern.name != pname {
            return false;
        }
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
                let rule = Py::new(
                    py,
                    PyRuleSet {
                        inner: leaf.clone(),
                    },
                )?;
                let mol = py_mol.bind(py);
                cb.bind(py)
                    .call1((mol, rule, info))?
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
                let mol = py_mol.bind(py);
                cb.bind(py)
                    .call1((mol, site, info.site))?
                    .extract::<bool>()
            }),
            &err,
        )
    };
    let out = set
        .metabolize(forest, rules, sites, true)
        .collect::<Result<Vec<_>, _>>()
        .map_err(py_err)?;
    if let Some(e) = err.into_inner() {
        return Err(e);
    }
    Ok(out)
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

    /// Named leaf catalog rule (`Hydroxylation`, `Dealkylation`, …).
    #[staticmethod]
    fn leaf(name: &str) -> PyResult<Self> {
        leaf_rule_rs(name)
            .map(|inner| Self { inner })
            .ok_or_else(|| PyValueError::new_err(format!("unknown leaf rule: {name:?}")))
    }

    /// Leaf names in catalog order.
    #[staticmethod]
    fn catalog_names() -> Vec<String> {
        catalog_names().iter().map(|s| (*s).to_string()).collect()
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

    /// Cross-language parity excuse, if any (see Rust `RuleSet::parity_exception`).
    #[getter]
    fn parity_exception(&self) -> Option<String> {
        self.inner.parity_exception.clone()
    }

    /// Direct member count (nested sets count as one member each).
    fn __len__(&self) -> usize {
        self.inner.members().len()
    }

    fn __contains__(&self, key: &str) -> bool {
        self.inner.contains_name(key)
    }

    /// Index by int or name: catalog → child ``RuleSet``; leaf → ``BoundPattern``.
    fn __getitem__(&self, key: Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            if let Ok(mut index) = key.extract::<isize>() {
                let len = self.inner.len() as isize;
                if index < 0 {
                    index += len;
                }
                if index < 0 || index >= len {
                    return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                        "RuleSet index out of range",
                    ));
                }
                let index = index as usize;
                if self.inner.is_catalog() {
                    let child = self.inner.get(index).ok_or_else(|| {
                        PyErr::new::<pyo3::exceptions::PyIndexError, _>("RuleSet index out of range")
                    })?;
                    return Ok(Py::new(py, PyRuleSet { inner: child })?.into_any());
                }
                let bp = self.inner.bound_pattern_at(index).ok_or_else(|| {
                    PyErr::new::<pyo3::exceptions::PyIndexError, _>("RuleSet index out of range")
                })?;
                return Ok(Py::new(py, PyBoundPattern { inner: bp })?.into_any());
            }
            let name: String = key.extract()?;
            if self.inner.is_catalog() {
                if let Some(child) = self.inner.get_str(&name) {
                    return Ok(Py::new(py, PyRuleSet { inner: child })?.into_any());
                }
            }
            if let Some(bp) = self.inner.bound_pattern(&name) {
                return Ok(Py::new(py, PyBoundPattern { inner: bp })?.into_any());
            }
            Err(PyErr::new::<pyo3::exceptions::PyKeyError, _>(name))
        })
    }

    /// Flat leaf patterns under this set (including nested members).
    fn patterns(&self) -> Vec<PyPatternInfo> {
        self.inner
            .patterns()
            .into_iter()
            .cloned()
            .map(|inner| PyPatternInfo { inner })
            .collect()
    }

    /// Run the owned members. `filter_rules` / `filter_sites` are optional Python
    /// callables; omit them and Rust `accept_all_*` runs with no GIL per site.
    ///
    /// Each row is `(pattern_name, site, products, rule_path)` where `rule_path`
    /// is leaf-first namespace names (`None` for an unnamed set).
    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyForestMol>,
        filter_rules: Option<Bound<'_, PyAny>>,
        filter_sites: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Vec<MetabolizeRow>> {
        let forest = mol.borrow().inner.copy_mol();
        let set = slf.borrow().inner.clone();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            set.metabolize(&forest, accept_all_rules, accept_all_sites, true)
                .collect::<Result<Vec<_>, _>>()
                .map_err(py_err)?
        } else {
            metabolize_with_python(mol, &set, &forest, filter_rules, filter_sites)?
        };
        // Explicit CSMI downgrade at the Python string-row boundary.
        Ok(emissions
            .into_iter()
            .map(|e| {
                let products = e.product_csmis();
                (e.pattern_name, e.site, products, e.rule_path)
            })
            .collect())
    }

    fn __repr__(&self) -> String {
        match &self.inner.name {
            Some(name) => format!(
                "RuleSet({name:?}, {} members, {} patterns)",
                self.inner.members().len(),
                self.inner.patterns().len()
            ),
            None => format!(
                "RuleSet({} members, {} patterns)",
                self.inner.members().len(),
                self.inner.patterns().len()
            ),
        }
    }
}

fn metabolize_with_python(
    mol: &Bound<'_, PyForestMol>,
    set: &RuleSet,
    forest: &crate::ForestMol,
    filter_rules: Option<Bound<'_, PyAny>>,
    filter_sites: Option<Bound<'_, PyAny>>,
) -> PyResult<Vec<crate::pattern::Emission>> {
    let py_rules = filter_rules.map(|cb| cb.unbind());
    let py_sites = filter_sites.map(|cb| cb.unbind());
    let py_mol = mol.clone().unbind();
    let err: RefCell<Option<PyErr>> = RefCell::new(None);
    let take_bool = |result: PyResult<bool>, slot: &RefCell<Option<PyErr>>| match result {
        Ok(keep) => keep,
        Err(e) => {
            *slot.borrow_mut() = Some(e);
            false
        }
    };
    // Filters see the leaf RuleSet (namespace), not the outer compose container.
    let rules = |_: &Molecule, leaf: &RuleSet, pattern: &PatternInfo| {
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
                let leaf_py = Py::new(
                    py,
                    PyRuleSet {
                        inner: leaf.clone(),
                    },
                )?;
                cb.bind(py)
                    .call1((py_mol.bind(py), leaf_py, info))?
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
        .metabolize(forest, rules, sites, true)
        .collect::<Result<Vec<_>, _>>()
        .map_err(py_err)?;
    if let Some(e) = err.into_inner() {
        return Err(e);
    }
    Ok(emissions)
}

/// Chematic tautomer pick adopted as a tagged [`ForestMol`].
///
/// Returns ``(mol, changed)`` where ``mol`` is a :class:`ForestMol` and
/// ``changed`` is whether the form differed from the input.
#[pyfunction]
fn normalize_tautomer(smiles: &str) -> PyResult<(PyForestMol, bool)> {
    let out = crate::normalize_tautomer(smiles).map_err(py_err)?;
    Ok((PyForestMol::wrap(out.mol), out.changed))
}

/// Native chematic ``find_path`` (default ruleset). Returns ``(hits, counters)``.
///
/// Default ruleset is QuinoneFormation + EpoxideHydration + Tautomerization +
/// PhaseOne core. Each hit is ``{"smiles": str, "steps": [{"rule": str,
/// "site": [...]}]``. Counters is a plain dict of the billed fields. Separate
/// from the Python RDKit ``xenosite.forest.find_path`` walk.
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
fn find_path(
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
) -> PyResult<(Vec<Py<PyAny>>, Py<PyAny>)> {
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
    let mut counters = PathCounters::default();
    let hits = match network {
        Some(net) => {
            let mut py_net = net.borrow_mut();
            find_path_with_network(
                reactant,
                target,
                &rules,
                &mut counters,
                config,
                Some(&mut py_net.inner),
                |_| true,
            )
            .map_err(py_err)?
            .collect_all()
            .map_err(py_err)?
        }
        None => find_path_with(reactant, target, &rules, &mut counters, config, |_| true)
            .map_err(py_err)?
            .collect_all()
            .map_err(py_err)?,
    };
    Ok((path_outcome_dicts(py, &hits)?, counters_dict(py, &counters)?))
}

/// Explored metabolic network (reactant root + hops). Mutated by search when
/// passed as ``network=``.
#[pyclass(name = "MetabolicNetwork", unsendable)]
struct PyMetabolicNetwork {
    inner: MetabolicNetwork,
}

#[pymethods]
impl PyMetabolicNetwork {
    #[new]
    fn new() -> Self {
        Self {
            inner: MetabolicNetwork::new(),
        }
    }

    #[getter]
    fn root_csmi(&self) -> Option<String> {
        self.inner.root_csmi.clone()
    }

    fn n_nodes(&self) -> usize {
        self.inner.n_nodes()
    }

    fn n_edges(&self) -> usize {
        self.inner.n_edges()
    }

    fn reaches(&self, target_csmi: &str) -> bool {
        self.inner.reaches(target_csmi)
    }

    fn closest(&self, target: &str, k: usize) -> PyResult<Vec<(String, usize)>> {
        self.inner.closest(target, k).map_err(py_err)
    }

    fn missed(&self, py: Python<'_>, csmi: &str, target: &str) -> PyResult<Option<Py<PyAny>>> {
        let residual = self.inner.missed(csmi, target).map_err(py_err)?;
        Ok(match residual {
            None => None,
            Some(r) => {
                let d = pyo3::types::PyDict::new(py);
                d.set_item("cost", r.cost)?;
                d.set_item("n_extra", r.n_extra)?;
                d.set_item("categories", r.categories.clone())?;
                d.set_item("unresolvable", r.unresolvable)?;
                Some(d.unbind().into_any())
            }
        })
    }
}

fn path_outcome_dicts(
    py: Python<'_>,
    hits: &[crate::find_path::PathOutcome],
) -> PyResult<Vec<Py<PyAny>>> {
    let mut out = Vec::with_capacity(hits.len());
    for hit in hits {
        let steps: Vec<Py<PyAny>> = hit
            .plan
            .iter()
            .map(|step| {
                let site: Vec<String> = step
                    .site
                    .iter()
                    .map(|a| match a {
                        crate::PlanAtom::Index(i) => i.to_string(),
                        other => format!("{other:?}"),
                    })
                    .collect();
                let d = pyo3::types::PyDict::new(py);
                d.set_item("rule", step.rule.as_str())?;
                d.set_item("site", site)?;
                Ok::<_, PyErr>(d.unbind().into_any())
            })
            .collect::<PyResult<_>>()?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("smiles", hit.smiles.as_str())?;
        d.set_item("steps", steps)?;
        out.push(d.unbind().into_any());
    }
    Ok(out)
}

fn counters_dict(py: Python<'_>, counters: &PathCounters) -> PyResult<Py<PyAny>> {
    let c = pyo3::types::PyDict::new(py);
    c.set_item("nodes", counters.nodes)?;
    c.set_item("mol_edits", counters.mol_edits)?;
    c.set_item("expansions", counters.expansions)?;
    c.set_item("billed", counters.billed())?;
    c.set_item("dropped_duplicate_plan", counters.dropped_duplicate_plan)?;
    c.set_item("dropped_exact_plan", counters.dropped_exact_plan)?;
    c.set_item("dropped_skeleton_twin", counters.dropped_skeleton_twin)?;
    c.set_item("diversity_repush", counters.diversity_repush)?;
    c.set_item("unstable_csmi_key", counters.unstable_csmi_key)?;
    c.set_item("timed_out", counters.timed_out)?;
    Ok(c.unbind().into_any())
}

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
fn find_path_partial_py(
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
) -> PyResult<(Vec<Py<PyAny>>, Vec<Py<PyAny>>, Py<PyAny>)> {
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
    let mut counters = PathCounters::default();
    let result = match network {
        Some(net) => {
            let mut py_net = net.borrow_mut();
            find_path_partial(
                reactant,
                target,
                &rules,
                &mut counters,
                config,
                Some(&mut py_net.inner),
                |_| true,
            )
            .map_err(py_err)?
        }
        None => find_path_partial(
            reactant,
            target,
            &rules,
            &mut counters,
            config,
            None,
            |_| true,
        )
        .map_err(py_err)?,
    };
    let FindPathPartialResult { exact, partials } = result;
    let exact_out = path_outcome_dicts(py, &exact)?;
    let mut partial_out = Vec::with_capacity(partials.len());
    for p in partials {
        let steps: Vec<Py<PyAny>> = p
            .plan
            .iter()
            .map(|step| {
                let site: Vec<String> = step
                    .site
                    .iter()
                    .map(|a| match a {
                        crate::PlanAtom::Index(i) => i.to_string(),
                        other => format!("{other:?}"),
                    })
                    .collect();
                let d = pyo3::types::PyDict::new(py);
                d.set_item("rule", step.rule.as_str())?;
                d.set_item("site", site)?;
                Ok::<_, PyErr>(d.unbind().into_any())
            })
            .collect::<PyResult<_>>()?;
        let residual = pyo3::types::PyDict::new(py);
        residual.set_item("cost", p.residual.cost)?;
        residual.set_item("n_extra", p.residual.n_extra)?;
        residual.set_item("categories", p.residual.categories.clone())?;
        residual.set_item("unresolvable", p.residual.unresolvable)?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("smiles", p.smiles.as_str())?;
        d.set_item("steps", steps)?;
        d.set_item("residual", residual)?;
        partial_out.push(d.unbind().into_any());
    }
    Ok((exact_out, partial_out, counters_dict(py, &counters)?))
}

/// Seeded random walk: apply up to ``max_steps`` rules. Returns a dict with
/// ``smiles``, ``path``, ``steps``, and ``patterns``.
///
/// Default ruleset is PhaseOne. Pass a ``RuleSet`` to override.
/// ``skip_multicomponent`` / ``skip_seen`` map to [`PathwayOptions`] (off by
/// default; same knobs for StepSequence / PathOutcome ``apply``).
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
fn random_path(
    py: Python<'_>,
    reactant: &str,
    seed: u64,
    max_steps: usize,
    ruleset: Option<&Bound<'_, PyRuleSet>>,
    skip_multicomponent: bool,
    skip_seen: bool,
) -> PyResult<Py<PyAny>> {
    let owned;
    let rules = match ruleset {
        Some(rs) => {
            owned = rs.borrow().inner.clone();
            &owned
        }
        None => {
            owned = phase_one_rs();
            &owned
        }
    };
    let options = PathwayOptions {
        skip_multicomponent,
        skip_seen,
    };
    let outcome = if options == PathwayOptions::default() {
        random_path_rs(reactant, seed, rules, max_steps)
    } else {
        random_path_with(reactant, seed, rules, max_steps, options)
    }
    .map_err(py_err)?;

    let steps: Vec<Py<PyAny>> = outcome
        .steps
        .iter()
        .map(|step| {
            let d = pyo3::types::PyDict::new(py);
            d.set_item("rule", step.rule.as_str())?;
            d.set_item("pattern", step.pattern_name.as_str())?;
            d.set_item("site", step.site.clone())?;
            d.set_item("products", step.products.clone())?;
            d.set_item("chosen", step.chosen)?;
            Ok::<_, PyErr>(d.unbind().into_any())
        })
        .collect::<PyResult<_>>()?;

    let patterns: Vec<Py<PyAny>> = outcome
        .patterns
        .iter()
        .map(|p| {
            let d = pyo3::types::PyDict::new(py);
            d.set_item("name", p.name.as_str())?;
            d.set_item("smarts", p.smarts.as_str())?;
            d.set_item("cleaves", p.effect.cleaves)?;
            d.set_item("search_bias", p.search_bias)?;
            Ok::<_, PyErr>(d.unbind().into_any())
        })
        .collect::<PyResult<_>>()?;

    let d = pyo3::types::PyDict::new(py);
    d.set_item("smiles", outcome.smiles.as_str())?;
    d.set_item("path", outcome.path.clone())?;
    d.set_item("steps", steps)?;
    d.set_item("patterns", patterns)?;
    Ok(d.unbind().into_any())
}

fn wrap_ruleset(inner: RuleSet) -> PyRuleSet {
    PyRuleSet { inner }
}

#[pyfunction]
fn phase_one() -> PyRuleSet {
    wrap_ruleset(phase_one_rs())
}

/// Look up a sealed leaf by catalog name (`LEAF_CTORS`).
///
/// Used by native↔Rust product parity over the shared coverage substrate pool.
#[pyfunction]
fn leaf_rule(name: &str) -> PyResult<PyRuleSet> {
    leaf_rule_rs(name)
        .map(wrap_ruleset)
        .ok_or_else(|| PyValueError::new_err(format!("unknown leaf rule: {name}")))
}

#[pyfunction]
fn epoxidation() -> PyRuleSet {
    wrap_ruleset(epoxidation_rs())
}

#[pyfunction]
fn quinone_formation() -> PyRuleSet {
    wrap_ruleset(quinone_formation_rs())
}

#[pyfunction]
fn epoxide_opening() -> PyRuleSet {
    wrap_ruleset(epoxide_opening_rs())
}

#[pyfunction]
fn n_dealkylation() -> PyRuleSet {
    wrap_ruleset(n_dealkylation_rs())
}

#[pyfunction]
fn hydroxylation() -> PyRuleSet {
    wrap_ruleset(hydroxylation_rs())
}

#[pyfunction]
fn dehydrogenation() -> PyRuleSet {
    wrap_ruleset(dehydrogenation_rs())
}

#[pyfunction]
fn dealkylation() -> PyRuleSet {
    wrap_ruleset(dealkylation_rs())
}

#[pyfunction]
fn hydrolysis() -> PyRuleSet {
    wrap_ruleset(hydrolysis_rs())
}

#[pyfunction]
fn default_ruleset() -> PyRuleSet {
    wrap_ruleset(default_ruleset_rs())
}

#[pyfunction]
fn forest_xmet_sssom() -> String {
    crate::mapping::forest_xmet_sssom().to_string()
}

#[pyfunction]
#[pyo3(name = "resolve")]
fn resolve_py(id: &str) -> PyResult<Py<PyAny>> {
    Python::attach(|py| match crate::mapping::resolve(id).map_err(py_err)? {
        crate::mapping::Resolved::Rule(inner) => {
            Ok(Py::new(py, PyRuleSet { inner })?.into_any())
        }
        crate::mapping::Resolved::Pattern(inner) => {
            Ok(Py::new(py, PyBoundPattern { inner })?.into_any())
        }
    })
}

#[pyfunction]
#[pyo3(name = "expand_iri")]
fn expand_iri_py(curie_or_iri: &str) -> String {
    crate::mapping::expand_iri(curie_or_iri)
}

#[pyfunction]
#[pyo3(name = "to_curie")]
fn to_curie_py(iri: &str) -> String {
    crate::mapping::to_curie(iri)
}

#[pymodule]
#[pyo3(name = "_rust")]
fn xenosite_forest(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForestMol>()?;
    m.add_class::<PyFormula>()?;
    m.add_class::<PyPatternInfo>()?;
    m.add_class::<PyRuleSet>()?;
    m.add_class::<PyBoundPattern>()?;
    m.add_class::<PyMetabolicNetwork>()?;
    m.add_function(wrap_pyfunction!(find_path, m)?)?;
    m.add_function(wrap_pyfunction!(find_path_partial_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_tautomer, m)?)?;
    m.add_function(wrap_pyfunction!(random_path, m)?)?;
    m.add_function(wrap_pyfunction!(phase_one, m)?)?;
    m.add_function(wrap_pyfunction!(leaf_rule, m)?)?;
    m.add_function(wrap_pyfunction!(epoxidation, m)?)?;
    m.add_function(wrap_pyfunction!(quinone_formation, m)?)?;
    m.add_function(wrap_pyfunction!(epoxide_opening, m)?)?;
    m.add_function(wrap_pyfunction!(n_dealkylation, m)?)?;
    m.add_function(wrap_pyfunction!(hydroxylation, m)?)?;
    m.add_function(wrap_pyfunction!(dehydrogenation, m)?)?;
    m.add_function(wrap_pyfunction!(dealkylation, m)?)?;
    m.add_function(wrap_pyfunction!(hydrolysis, m)?)?;
    m.add_function(wrap_pyfunction!(default_ruleset, m)?)?;
    m.add_function(wrap_pyfunction!(forest_xmet_sssom, m)?)?;
    m.add_function(wrap_pyfunction!(resolve_py, m)?)?;
    m.add_function(wrap_pyfunction!(expand_iri_py, m)?)?;
    m.add_function(wrap_pyfunction!(to_curie_py, m)?)?;
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
            let csmi = mol.getattr("csmi").unwrap();
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
            assert!(copied.getattr("csmi").unwrap().is(&csmi));
            let _edited = mol.call_method0("edit_copy").unwrap();
            mol.call_method0("clear_structure").unwrap();
            let after = mol.getattr("csmi").unwrap();
            assert_eq!(
                after.extract::<String>().unwrap(),
                csmi.extract::<String>().unwrap()
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
            let products: Vec<MetabolizeRow> = products.extract().unwrap();
            assert_eq!(products.len(), 1);
            assert_eq!(products[0].0, "h");
            assert_eq!(products[0].3, vec![Some("Hydroxylation".into())]);
            let filt = py
                .eval(c"lambda m, rule, p: p.name == 'h2'", None, None)
                .unwrap();
            let kept: Vec<MetabolizeRow> = rs
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
            // Nested namespaces: two members, three flat patterns.
            assert_eq!(
                composed
                    .call_method0("__len__")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                2
            );
            assert_eq!(composed.call_method0("patterns").unwrap().len().unwrap(), 3);
        });
    }
}
