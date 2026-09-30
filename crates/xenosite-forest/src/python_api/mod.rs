//! PyO3 product door (`xenosite_forest` / `xenosite.forest`).
//!
//! New public APIs land here — not in `xenosite.forest.native` (frozen RDKit
//! reference; see `docs/forest/NATIVE.md`).
//!
//! Native-only: CPython C-API. Not compiled for `wasm32-unknown-unknown`.

mod common;
mod display;
mod graph;
mod macros;
mod mol;
mod path;
mod plan;
mod rules;
mod walk;

pub use graph::{
    PyGraphEdge, PyGraphNode, PyMetabolicNetwork, product_graph_bfs, product_graph_into_py,
};
pub use mol::{PyForestMol, PyFormula, normalize_tautomer};
pub use path::{
    PyPartialPathOutcome, PyPathCounters, PyPathOutcome, find_path, find_path_partial_py,
};
pub use plan::PyStepPlan;
pub use rules::{
    PyBoundPattern, PyEmission, PyPatternInfo, PyRuleSet, dealkylation, default_ruleset,
    dehydrogenation, epoxidation, epoxide_opening, expand_iri_py, forest_xmet_sssom, hydrolysis,
    hydroxylation, leaf_rule, n_dealkylation, phase_one, product_graph_ruleset, quinone_formation,
    reactivity, resolve_py, to_curie_py,
};
pub use walk::{PyRandomPathOutcome, random_path};

use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "_rust")]
fn xenosite_forest(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForestMol>()?;
    m.add_class::<PyFormula>()?;
    m.add_class::<PyPatternInfo>()?;
    m.add_class::<PyRuleSet>()?;
    m.add_class::<PyBoundPattern>()?;
    m.add_class::<PyMetabolicNetwork>()?;
    m.add_class::<PyGraphNode>()?;
    m.add_class::<PyGraphEdge>()?;
    m.add_class::<PyPathOutcome>()?;
    m.add_class::<PyPartialPathOutcome>()?;
    m.add_class::<PyPathCounters>()?;
    m.add_class::<PyStepPlan>()?;
    m.add_class::<PyEmission>()?;
    m.add_class::<PyRandomPathOutcome>()?;
    m.add_function(wrap_pyfunction!(find_path, m)?)?;
    m.add_function(wrap_pyfunction!(product_graph_bfs, m)?)?;
    m.add_function(wrap_pyfunction!(product_graph_into_py, m)?)?;
    m.add_function(wrap_pyfunction!(find_path_partial_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_tautomer, m)?)?;
    m.add_function(wrap_pyfunction!(random_path, m)?)?;
    m.add_function(wrap_pyfunction!(phase_one, m)?)?;
    m.add_function(wrap_pyfunction!(reactivity, m)?)?;
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
    m.add_function(wrap_pyfunction!(product_graph_ruleset, m)?)?;
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
            let emissions: Vec<Py<PyEmission>> = products.extract().unwrap();
            assert_eq!(emissions.len(), 1);
            let e0 = emissions[0].bind(py);
            assert_eq!(
                e0.getattr("pattern_name")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "h"
            );
            assert_eq!(
                e0.call_method0("rule_path")
                    .unwrap()
                    .extract::<Vec<Option<String>>>()
                    .unwrap(),
                vec![Some("Hydroxylation".into())]
            );
            let filt = py
                .eval(c"lambda m, rule, p: p.name == 'h2'", None, None)
                .unwrap();
            let kept: Vec<Py<PyEmission>> = rs
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
