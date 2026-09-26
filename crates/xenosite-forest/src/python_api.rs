//! PyO3 class wrap of [`crate::forest_mol::ForestMol`].
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

use crate::enumerate::{EnumConfig, bfs as enum_bfs, dfs as enum_dfs, enumerate_metabolites};
use crate::find_path::{FindPathConfig, HeapScoreMode, PathCounters, find_path_with};
use crate::forest::Formula;
use crate::forest_mol::ForestMol;
use crate::mol::{Molecule, ranks};
use crate::pattern::{Edit, Effect, PatternInfo, SiteInfo};
use crate::rules::{catalog_names, default_ruleset, leaf_rule, phase_one};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};

fn py_err(err: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// One `RuleSet.metabolize` row:
/// `(pattern_name, site, site_atoms, site_orbit, products, rule_path)`.
type MetabolizeRow = (
    String,
    usize,
    Vec<usize>,
    Vec<usize>,
    Vec<String>,
    Vec<Option<String>>,
);

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

    fn smarts_matches(&self, smarts: &str) -> PyResult<Vec<HashMap<u16, usize>>> {
        let hits = self.inner.smarts_matches(smarts).map_err(py_err)?;
        Ok(hits
            .iter()
            .map(|mapped| mapped.iter().map(|(&k, &v)| (k, v)).collect())
            .collect())
    }

    /// Topological equivalence ranks (one per atom index). Same role as
    /// RDKit ``CanonicalRankAtoms(..., breakTies=False)`` for site identity.
    fn ranks(&self) -> Vec<usize> {
        ranks(self.inner.mol())
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
        leaf_rule(name)
            .map(|inner| Self { inner })
            .ok_or_else(|| PyValueError::new_err(format!("unknown leaf rule {name:?}")))
    }

    /// Leaf names in catalog order.
    #[staticmethod]
    fn catalog_names() -> Vec<String> {
        catalog_names().iter().map(|s| (*s).to_string()).collect()
    }

    #[staticmethod]
    fn phase_one() -> Self {
        Self { inner: phase_one() }
    }

    #[staticmethod]
    fn default_ruleset() -> Self {
        Self {
            inner: default_ruleset(),
        }
    }

    #[staticmethod]
    fn all_rules() -> Self {
        Self {
            inner: crate::rules::all_rules(),
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

    /// Cross-language parity excuse, if any (see Rust `RuleSet::parity_exception`).
    #[getter]
    fn parity_exception(&self) -> Option<String> {
        self.inner.parity_exception.clone()
    }

    /// Direct member count (nested sets count as one member each).
    fn __len__(&self) -> usize {
        self.inner.members().len()
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
    /// Each row is
    /// `(pattern_name, site, site_atoms, site_orbit, products, rule_path)`
    /// where `rule_path` is leaf-first namespace names (`None` for an unnamed set).
    #[pyo3(signature = (mol, filter_rules=None, filter_sites=None))]
    fn metabolize(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyForestMol>,
        filter_rules: Option<Bound<'_, PyAny>>,
        filter_sites: Option<Bound<'_, PyAny>>,
    ) -> PyResult<Vec<MetabolizeRow>> {
        let chemistry = mol.borrow().inner.mol().clone();
        let set = slf.borrow().inner.clone();
        let emissions = if filter_rules.is_none() && filter_sites.is_none() {
            set.metabolize(&chemistry, accept_all_rules, accept_all_sites, true)
                .collect::<Result<Vec<_>, _>>()
                .map_err(py_err)?
        } else {
            metabolize_with_python(mol, &set, &chemistry, filter_rules, filter_sites)?
        };
        Ok(emissions
            .into_iter()
            .map(|e| {
                (
                    e.pattern_name,
                    e.site,
                    e.site_atoms,
                    e.site_orbit,
                    e.products,
                    e.rule_path,
                )
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
    chemistry: &Molecule,
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
        .metabolize(chemistry, rules, sites, true)
        .collect::<Result<Vec<_>, _>>()
        .map_err(py_err)?;
    if let Some(e) = err.into_inner() {
        return Err(e);
    }
    Ok(emissions)
}

/// Native chematic ``find_path``. Returns ``(hits, counters)``.
///
/// Default catalog is PhaseOne; pass ``ruleset=`` to override. Each hit is
/// ``{"smiles": str, "steps": [{"rule": str, "site": [...]}]}``.
#[pyfunction]
#[pyo3(signature = (
    reactant,
    target,
    *,
    ruleset=None,
    max_paths=1,
    max_nodes=800,
    use_atom_diff=true,
    lazy_closer=false,
    diversity=false,
    drop_skeleton_twins=true,
    score="log-neg-pc",
))]
fn find_path(
    py: Python<'_>,
    reactant: &str,
    target: &str,
    ruleset: Option<PyRef<'_, PyRuleSet>>,
    max_paths: usize,
    max_nodes: usize,
    use_atom_diff: bool,
    lazy_closer: bool,
    diversity: bool,
    drop_skeleton_twins: bool,
    score: &str,
) -> PyResult<(Vec<Py<PyAny>>, Py<PyAny>)> {
    let heap_score = HeapScoreMode::from_label(score).ok_or_else(|| {
        PyValueError::new_err(format!(
            "unknown score {score:?}; try log-neg-pc, soft, add-both, …"
        ))
    })?;
    let config = FindPathConfig {
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
        heap_score,
        drop_skeleton_twins,
        diversity,
    };
    let owned;
    let rules: &RuleSet = match &ruleset {
        Some(rs) => &rs.inner,
        None => {
            owned = phase_one();
            &owned
        }
    };
    let mut counters = PathCounters::default();
    let hits = find_path_with(reactant, target, rules, &mut counters, config, |_| true)
        .map_err(py_err)?
        .collect_all()
        .map_err(py_err)?;

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
                        crate::PlanAtom::Label(t) => t.0.to_string(),
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
    Ok((out, c.unbind().into_any()))
}

fn enum_hits_to_py(
    py: Python<'_>,
    iter: impl Iterator<Item = Result<crate::enumerate::Metabolite, crate::ForestError>>,
) -> PyResult<Vec<Py<PyAny>>> {
    let mut out = Vec::new();
    for hit in iter {
        let hit = hit.map_err(py_err)?;
        let hops: Vec<Py<PyAny>> = hit
            .path
            .hops
            .iter()
            .map(|hop| {
                let d = pyo3::types::PyDict::new(py);
                d.set_item("rule", hop.rule.as_str())?;
                d.set_item("pattern", hop.pattern_name.as_str())?;
                d.set_item("site", hop.site)?;
                d.set_item("products", hop.products.clone())?;
                d.set_item("cleaves", hop.cleaves)?;
                Ok::<_, PyErr>(d.unbind().into_any())
            })
            .collect::<PyResult<_>>()?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("smiles", hit.smiles())?;
        d.set_item("depth", hit.path.depth())?;
        d.set_item("hops", hops)?;
        out.push(d.unbind().into_any());
    }
    Ok(out)
}

/// Breadth-first metabolites (default catalog: ``default_ruleset``).
#[pyfunction]
#[pyo3(signature = (reactant, max_depth, *, ruleset=None))]
fn bfs(
    py: Python<'_>,
    reactant: &str,
    max_depth: usize,
    ruleset: Option<PyRef<'_, PyRuleSet>>,
) -> PyResult<Vec<Py<PyAny>>> {
    let owned;
    let rules: &RuleSet = match &ruleset {
        Some(rs) => &rs.inner,
        None => {
            owned = default_ruleset();
            &owned
        }
    };
    let stream = enum_bfs(reactant, rules, max_depth).map_err(py_err)?;
    enum_hits_to_py(py, stream)
}

/// Depth-first metabolites (default catalog: ``default_ruleset``).
#[pyfunction]
#[pyo3(signature = (reactant, max_depth, *, ruleset=None))]
fn dfs(
    py: Python<'_>,
    reactant: &str,
    max_depth: usize,
    ruleset: Option<PyRef<'_, PyRuleSet>>,
) -> PyResult<Vec<Py<PyAny>>> {
    let owned;
    let rules: &RuleSet = match &ruleset {
        Some(rs) => &rs.inner,
        None => {
            owned = default_ruleset();
            &owned
        }
    };
    let stream = enum_dfs(reactant, rules, max_depth).map_err(py_err)?;
    enum_hits_to_py(py, stream)
}

/// Enumerate metabolites with explicit depth / order / dedup knobs.
#[pyfunction]
#[pyo3(signature = (
    reactant,
    *,
    ruleset=None,
    max_depth=1,
    order="bfs",
    max_nodes=0,
    unique_csmi=true,
))]
fn enumerate(
    py: Python<'_>,
    reactant: &str,
    ruleset: Option<PyRef<'_, PyRuleSet>>,
    max_depth: usize,
    order: &str,
    max_nodes: usize,
    unique_csmi: bool,
) -> PyResult<Vec<Py<PyAny>>> {
    let owned;
    let rules: &RuleSet = match &ruleset {
        Some(rs) => &rs.inner,
        None => {
            owned = default_ruleset();
            &owned
        }
    };
    let mut config = match order {
        "bfs" => EnumConfig::bfs(max_depth),
        "dfs" => EnumConfig::dfs(max_depth),
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown order {other:?}; use bfs or dfs"
            )));
        }
    };
    config = config.with_max_nodes(max_nodes).with_unique_csmi(unique_csmi);
    let stream = enumerate_metabolites(reactant, rules, config).map_err(py_err)?;
    enum_hits_to_py(py, stream)
}

#[pymodule]
fn xenosite_forest(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForestMol>()?;
    m.add_class::<PyFormula>()?;
    m.add_class::<PyPatternInfo>()?;
    m.add_class::<PyRuleSet>()?;
    m.add_function(wrap_pyfunction!(find_path, m)?)?;
    m.add_function(wrap_pyfunction!(bfs, m)?)?;
    m.add_function(wrap_pyfunction!(dfs, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate, m)?)?;
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
            assert_eq!(products[0].5, vec![Some("Hydroxylation".into())]);
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
