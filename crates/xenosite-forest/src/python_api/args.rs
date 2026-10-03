//! Typed argument / return wrappers for `Bound<PyAny>` edges of the API.
//!
//! Runtime conversions live here; their `PyStubType`s (in [`stub_types`],
//! `stubs` feature only) derive the `_rust.pyi` signatures from these Rust
//! types rather than string overrides.

use std::marker::PhantomData;

use pyo3::prelude::*;

use super::mol::{PyForestMol, py_forest_mol_ref};
use super::rules::{PyBoundPattern, PyRuleSet};

#[cfg(feature = "stubs")]
pub use stub_types::literal_str;

/// A Rust `ForestMol` pyclass or the `xenosite.forest.mol.ForestMol` wrapper.
///
/// Keeps the original object (filters are called with what the caller passed).
pub struct MolArg<'py> {
    pub obj: Bound<'py, PyAny>,
    pub mol: PyRef<'py, PyForestMol>,
}

impl<'py> FromPyObject<'py> for MolArg<'py> {
    fn extract_bound(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        Ok(Self {
            obj: obj.clone(),
            mol: py_forest_mol_ref(obj)?,
        })
    }
}

/// Python callable; `S` names its stub signature (see `CallbackSig`).
pub struct Callback<'py, S> {
    pub func: Bound<'py, PyAny>,
    sig: PhantomData<S>,
}

impl<'py, S> FromPyObject<'py> for Callback<'py, S> {
    fn extract_bound(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        Ok(Self {
            func: obj.clone(),
            sig: PhantomData,
        })
    }
}

/// `filter_rules(mol, leaf_rule, pattern) -> bool`.
pub struct RuleFilter;

/// `RuleSet.metabolize` `filter_sites(mol, site, pattern) -> bool`.
pub struct SiteFilter;

/// `BoundPattern.metabolize` `filter_sites(mol, site, primary_site) -> bool`.
pub struct BoundSiteFilter;

/// `__getitem__` key: position or name.
#[derive(FromPyObject)]
pub enum IndexKey {
    Index(isize),
    Name(String),
}

/// Catalog child (`RuleSet`) or leaf pattern handle (`BoundPattern`).
// Return-only value converted straight into a Python object; boxing buys nothing.
#[allow(clippy::large_enum_variant)]
#[derive(IntoPyObject)]
pub enum RuleOrPattern {
    Rule(PyRuleSet),
    Pattern(PyBoundPattern),
}

#[cfg(feature = "stubs")]
mod stub_types {
    use pyo3_stub_gen::{PyStubType, TypeInfo};

    use super::super::mol::PyForestMol;
    use super::super::rules::{PyBoundPattern, PyPatternInfo, PyRuleSet};
    use super::*;
    use crate::export::AttrValue;

    /// `typing.Literal[...]` over string values.
    pub fn literal_str(values: &[&str]) -> TypeInfo {
        let body = values
            .iter()
            .map(|v| format!("{v:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        TypeInfo::with_module(&format!("typing.Literal[{body}]"), "typing".into())
    }

    /// `collections.abc.Callable[[args...], ret]`, merging argument imports.
    fn callable(args: Vec<TypeInfo>, ret: TypeInfo) -> TypeInfo {
        let mut out = TypeInfo::with_module("", "collections.abc".into());
        let names: Vec<String> = args.iter().map(|a| a.name.clone()).collect();
        out.name = format!(
            "collections.abc.Callable[[{}], {}]",
            names.join(", "),
            ret.name
        );
        for t in args.into_iter().chain(std::iter::once(ret)) {
            out.import.extend(t.import);
            out.type_refs.extend(t.type_refs);
        }
        out
    }

    impl PyStubType for MolArg<'_> {
        fn type_output() -> TypeInfo {
            // ``from xenosite.forest import mol`` → ``mol.ForestMol``. No
            // ``type_refs`` entry: it is keyed by bare name and would collide
            // with the extension's own ``ForestMol``.
            let mut wrapper = TypeInfo::locally_defined("ForestMol", "xenosite.forest.mol".into());
            wrapper.type_refs.clear();
            PyForestMol::type_output() | wrapper
        }
    }

    /// Python signature of a [`Callback`].
    pub trait CallbackSig {
        fn args() -> Vec<TypeInfo>;
        fn ret() -> TypeInfo {
            bool::type_output()
        }
    }

    impl<S: CallbackSig> PyStubType for Callback<'_, S> {
        fn type_output() -> TypeInfo {
            callable(S::args(), S::ret())
        }
    }

    impl CallbackSig for RuleFilter {
        fn args() -> Vec<TypeInfo> {
            vec![
                MolArg::type_input(),
                PyRuleSet::type_output(),
                PyPatternInfo::type_output(),
            ]
        }
    }

    impl CallbackSig for SiteFilter {
        fn args() -> Vec<TypeInfo> {
            vec![
                MolArg::type_input(),
                usize::type_output(),
                PyPatternInfo::type_output(),
            ]
        }
    }

    impl CallbackSig for BoundSiteFilter {
        fn args() -> Vec<TypeInfo> {
            vec![
                MolArg::type_input(),
                usize::type_output(),
                usize::type_output(),
            ]
        }
    }

    impl PyStubType for IndexKey {
        fn type_output() -> TypeInfo {
            isize::type_output() | String::type_output()
        }
    }

    impl PyStubType for RuleOrPattern {
        fn type_output() -> TypeInfo {
            PyRuleSet::type_output() | PyBoundPattern::type_output()
        }
    }

    impl PyStubType for AttrValue {
        fn type_output() -> TypeInfo {
            bool::type_output() | i64::type_output() | String::type_output()
        }
    }
}
