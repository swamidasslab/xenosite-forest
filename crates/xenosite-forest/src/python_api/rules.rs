//! RuleSet, BoundPattern, Emission, PatternInfo, and catalog factories.

use std::cell::RefCell;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubs")]
use pyo3_stub_gen::derive::*;

use crate::export::EmissionDict;
use crate::mol::Molecule;
use crate::pattern::{Edit, Effect, Emission, PatternInfo, SiteInfo, SiteKind};
use crate::rules::{
    catalog_names, dealkylation as dealkylation_rs, default_ruleset as default_ruleset_rs,
    dehydrogenation as dehydrogenation_rs, epoxidation as epoxidation_rs,
    epoxide_opening as epoxide_opening_rs, hydrolysis as hydrolysis_rs,
    hydroxylation as hydroxylation_rs, leaf_rule as leaf_rule_rs,
    n_dealkylation as n_dealkylation_rs, phase_one as phase_one_rs,
    product_graph_ruleset as product_graph_ruleset_rs, quinone_formation as quinone_formation_rs,
    reactivity as reactivity_rs,
};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};

#[cfg(feature = "stubs")]
use super::args::literal_str;
use super::args::{
    BoundSiteFilter, Callback, IndexKey, MolArg, RuleFilter, RuleOrPattern, SiteFilter,
};
use super::common::py_err;
use super::mol::PyForestMol;

fn wrap_emissions(emissions: Vec<Emission>) -> Vec<PyEmission> {
    emissions
        .into_iter()
        .map(|inner| PyEmission { inner })
        .collect()
}

#[cfg(feature = "stubs")]
const SITE_KINDS: [SiteKind; 4] = [
    SiteKind::Atom,
    SiteKind::Bond,
    SiteKind::DirectedBond,
    SiteKind::AtomPair,
];

/// ``site_kind`` string; stub type is the ``Literal`` of [`site_kind_str`].
pub struct SiteKindName(&'static str);

#[cfg(feature = "stubs")]
impl pyo3_stub_gen::PyStubType for SiteKindName {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        literal_str(&SITE_KINDS.map(site_kind_str))
    }
}

impl<'py> IntoPyObject<'py> for SiteKindName {
    type Target = pyo3::types::PyString;
    type Output = Bound<'py, pyo3::types::PyString>;
    type Error = std::convert::Infallible;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        self.0.into_pyobject(py)
    }
}

fn site_kind_str(kind: SiteKind) -> &'static str {
    match kind {
        SiteKind::Atom => "atom",
        SiteKind::Bond => "bond",
        SiteKind::DirectedBond => "directed_bond",
        SiteKind::AtomPair => "atom_pair",
    }
}

/// Look up ``site_kind`` on the leaf rule pattern named by the emission.
fn site_kind_for_emission(em: &Emission) -> &'static str {
    let Some(leaf_name) = em.leaf_rule() else {
        return "atom";
    };
    let Some(set) = leaf_rule_rs(leaf_name) else {
        return "atom";
    };
    for pattern in set.patterns() {
        if pattern.name == em.pattern_name {
            return site_kind_str(pattern.site_kind);
        }
    }
    // Composite pair names (``left+right``) still come from an AtomPair leaf.
    if em.site_atoms.len() >= 2 {
        return "atom_pair";
    }
    "atom"
}

/// One applied metabolize row (Python ``Emission``).
#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "Emission", unsendable)]
pub struct PyEmission {
    pub(crate) inner: Emission,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
#[pymethods]
impl PyEmission {
    #[getter]
    fn pattern_name(&self) -> &str {
        self.inner.pattern_name.as_str()
    }

    #[getter]
    fn site(&self) -> usize {
        self.inner.site
    }

    /// Discovery site atoms on the reactant (bond ends / pair ends / atom).
    #[getter]
    fn site_atoms(&self) -> Vec<usize> {
        self.inner.site_atoms.clone()
    }

    /// From the leaf rule's pattern ``site_kind`` (``atom`` / ``bond`` /
    /// ``directed_bond`` / ``atom_pair``) — how notebook SOM is drawn.
    #[getter]
    fn site_kind(&self) -> SiteKindName {
        SiteKindName(site_kind_for_emission(&self.inner))
    }

    /// Trace labels for [`Self::site_atoms`] on the reactant (notebook SOM).
    #[getter]
    fn site_tags(&self) -> Vec<u16> {
        self.inner
            .site_atoms
            .iter()
            .filter_map(|&i| self.inner.reactant.tag_of(i).map(|t| t.get()))
            .collect()
    }

    fn product_csmis(&self) -> Vec<String> {
        self.inner.product_csmis()
    }

    fn products(&self) -> Vec<PyForestMol> {
        self.inner
            .products
            .iter()
            .map(|m| PyForestMol::wrap(m.clone()))
            .collect()
    }

    /// Mol the pattern was applied to (`site` indexes this mol).
    #[getter]
    fn reactant(&self) -> PyForestMol {
        PyForestMol::wrap(self.inner.reactant.clone())
    }

    fn rule_path(&self) -> Vec<Option<String>> {
        self.inner.rule_path.clone()
    }

    fn to_dict(&self) -> EmissionDict {
        EmissionDict::from(&self.inner)
    }

    fn __str__(&self) -> String {
        super::display::format_emission(&self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "Emission({:?}, site={})",
            self.inner.pattern_name, self.inner.site
        )
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
#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "PatternInfo", frozen)]
#[derive(Clone)]
pub struct PyPatternInfo {
    pub(crate) inner: PatternInfo,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
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

    /// ``atom`` / ``bond`` / ``directed_bond`` / ``atom_pair``.
    #[getter]
    fn site_kind(&self) -> SiteKindName {
        SiteKindName(site_kind_str(self.inner.site_kind))
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

    fn __str__(&self) -> String {
        super::display::format_pattern_info(&self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "PatternInfo({:?}, {:?})",
            self.inner.name, self.inner.smarts
        )
    }
}

/// Python wrap of [`RuleSet`]. Patterns are copied into Rust at construction.
#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "RuleSet", unsendable)]
pub struct PyRuleSet {
    pub(crate) inner: RuleSet,
}

/// Python wrap of [`crate::bound_pattern::BoundPattern`].
#[cfg_attr(feature = "stubs", gen_stub_pyclass)]
#[pyclass(name = "BoundPattern", unsendable)]
#[derive(Clone)]
pub struct PyBoundPattern {
    inner: crate::bound_pattern::BoundPattern,
}

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
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

    fn __getitem__(&self, key: IndexKey) -> PyResult<PyBoundPattern> {
        let index_err =
            || PyErr::new::<pyo3::exceptions::PyIndexError, _>("BoundPattern index out of range");
        match key {
            IndexKey::Index(index) => {
                let index = usize::try_from(index).map_err(|_| index_err())?;
                self.inner
                    .get(index)
                    .map(|inner| PyBoundPattern { inner })
                    .ok_or_else(index_err)
            }
            IndexKey::Name(name) => self
                .inner
                .get_str(&name)
                .map(|inner| PyBoundPattern { inner })
                .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>(name)),
        }
    }

    /// Releases the GIL when no Python filters are passed.
    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        py: Python<'_>,
        mol: MolArg<'_>,
        filter_rules: Option<Callback<'_, RuleFilter>>,
        filter_sites: Option<Callback<'_, BoundSiteFilter>>,
    ) -> PyResult<Vec<PyEmission>> {
        let forest = mol.mol.inner.copy_mol();
        let filter_rules = filter_rules.map(|c| c.func);
        let filter_sites = filter_sites.map(|c| c.func);
        let mol = &mol.obj;
        let bp = slf.borrow().inner.clone();
        let set = bp.rule().clone();
        let pattern_name = bp.name().to_string();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            py.detach(move || {
                bp.metabolize_default(&forest, true)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| e.to_string())
            })
            .map_err(PyValueError::new_err)?
        } else {
            // Compose BoundPattern name filter with optional Python callables.
            let name = pattern_name.clone();
            metabolize_with_python_bound(mol, &set, &forest, &name, filter_rules, filter_sites)?
        };
        let _ = pattern_name;
        Ok(wrap_emissions(emissions))
    }

    fn __repr__(&self) -> String {
        format!("BoundPattern({:?})", self.inner.curie())
    }

    fn __str__(&self) -> String {
        format!(
            "BoundPattern {}  ({} patterns under leaf)",
            self.inner.curie(),
            self.inner.len()
        )
    }
}

fn metabolize_with_python_bound(
    mol: &Bound<'_, PyAny>,
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
                cb.bind(py).call1((mol, rule, info))?.extract::<bool>()
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
                cb.bind(py).call1((mol, site, info.site))?.extract::<bool>()
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

#[cfg_attr(feature = "stubs", gen_stub_pymethods)]
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
    fn __getitem__(&self, key: IndexKey) -> PyResult<RuleOrPattern> {
        let index_err =
            || PyErr::new::<pyo3::exceptions::PyIndexError, _>("RuleSet index out of range");
        match key {
            IndexKey::Index(mut index) => {
                let len = self.inner.len() as isize;
                if index < 0 {
                    index += len;
                }
                if index < 0 || index >= len {
                    return Err(index_err());
                }
                let index = index as usize;
                if self.inner.is_catalog() {
                    let child = self.inner.get(index).ok_or_else(index_err)?;
                    return Ok(RuleOrPattern::Rule(PyRuleSet { inner: child }));
                }
                let bp = self.inner.bound_pattern_at(index).ok_or_else(index_err)?;
                Ok(RuleOrPattern::Pattern(PyBoundPattern { inner: bp }))
            }
            IndexKey::Name(name) => {
                if self.inner.is_catalog()
                    && let Some(child) = self.inner.get_str(&name)
                {
                    return Ok(RuleOrPattern::Rule(PyRuleSet { inner: child }));
                }
                if let Some(bp) = self.inner.bound_pattern(&name) {
                    return Ok(RuleOrPattern::Pattern(PyBoundPattern { inner: bp }));
                }
                Err(PyErr::new::<pyo3::exceptions::PyKeyError, _>(name))
            }
        }
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
    ///
    /// Releases the GIL when no Python filters are passed.
    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        py: Python<'_>,
        mol: MolArg<'_>,
        filter_rules: Option<Callback<'_, RuleFilter>>,
        filter_sites: Option<Callback<'_, SiteFilter>>,
    ) -> PyResult<Vec<PyEmission>> {
        let forest = mol.mol.inner.copy_mol();
        let filter_rules = filter_rules.map(|c| c.func);
        let filter_sites = filter_sites.map(|c| c.func);
        let mol = &mol.obj;
        let set = slf.borrow().inner.clone();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            py.detach(move || {
                set.metabolize(&forest, accept_all_rules, accept_all_sites, true)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| e.to_string())
            })
            .map_err(PyValueError::new_err)?
        } else {
            metabolize_with_python(mol, &set, &forest, filter_rules, filter_sites)?
        };
        Ok(wrap_emissions(emissions))
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

    fn __str__(&self) -> String {
        super::display::format_ruleset(&self.inner)
    }
}

fn metabolize_with_python(
    mol: &Bound<'_, PyAny>,
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

fn wrap_ruleset(inner: RuleSet) -> PyRuleSet {
    PyRuleSet { inner }
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn phase_one() -> PyRuleSet {
    wrap_ruleset(phase_one_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn reactivity() -> PyRuleSet {
    wrap_ruleset(reactivity_rs())
}

/// Look up a sealed leaf by catalog name (`LEAF_CTORS`).
///
/// Used by native↔Rust product parity over the shared coverage substrate pool.
#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn leaf_rule(name: &str) -> PyResult<PyRuleSet> {
    leaf_rule_rs(name)
        .map(wrap_ruleset)
        .ok_or_else(|| PyValueError::new_err(format!("unknown leaf rule: {name}")))
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn epoxidation() -> PyRuleSet {
    wrap_ruleset(epoxidation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn quinone_formation() -> PyRuleSet {
    wrap_ruleset(quinone_formation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn epoxide_opening() -> PyRuleSet {
    wrap_ruleset(epoxide_opening_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn n_dealkylation() -> PyRuleSet {
    wrap_ruleset(n_dealkylation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn hydroxylation() -> PyRuleSet {
    wrap_ruleset(hydroxylation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn dehydrogenation() -> PyRuleSet {
    wrap_ruleset(dehydrogenation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn dealkylation() -> PyRuleSet {
    wrap_ruleset(dealkylation_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn hydrolysis() -> PyRuleSet {
    wrap_ruleset(hydrolysis_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn default_ruleset() -> PyRuleSet {
    wrap_ruleset(default_ruleset_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn product_graph_ruleset() -> PyRuleSet {
    wrap_ruleset(product_graph_ruleset_rs())
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
pub fn forest_xmet_sssom() -> String {
    crate::mapping::forest_xmet_sssom().to_string()
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "resolve")]
pub fn resolve_py(id: &str) -> PyResult<RuleOrPattern> {
    Ok(match crate::mapping::resolve(id).map_err(py_err)? {
        crate::mapping::Resolved::Rule(inner) => RuleOrPattern::Rule(PyRuleSet { inner }),
        crate::mapping::Resolved::Pattern(inner) => {
            RuleOrPattern::Pattern(PyBoundPattern { inner })
        }
    })
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "expand_iri")]
pub fn expand_iri_py(curie_or_iri: &str) -> String {
    crate::mapping::expand_iri(curie_or_iri)
}

#[cfg_attr(feature = "stubs", gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "to_curie")]
pub fn to_curie_py(iri: &str) -> String {
    crate::mapping::to_curie(iri)
}
