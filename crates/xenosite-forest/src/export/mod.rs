//! Serde export views for explicit `.to_dict()` / JSON — not the primary API.
//!
//! Every view is declared with [`typed_dict!`]: one Rust struct is both the
//! runtime dict (serde → `pythonize`) and a `typing.TypedDict`: defined on
//! `_rust` at import (`python`) and written to `xenosite/forest/_rust.pyi`
//! (`stubs`, dev-only). Add fields here, then `make stubs` — never edit the
//! `.pyi` by hand.

use std::collections::BTreeMap;

use serde::Serialize;

use crate::atom_diff::AtomDiffResidual;
use crate::canonical_plan::{CleavageSide, Maybe, PlanAtom, Step};
use crate::find_path::{PathCounters, PathOutcome, PathStep};
use crate::metabolic_network::{AttrMap, GraphValue};
use crate::pattern::Emission;
use crate::pattern::PatternInfo;
use crate::random_path::{RandomPathOutcome, RandomPathStep};

/// One field of a [`typed_dict!`] struct.
#[cfg(feature = "python")]
pub struct TypedDictField {
    pub name: &'static str,
    pub doc: &'static str,
    #[cfg(feature = "stubs")]
    pub r#type: fn() -> pyo3_stub_gen::TypeInfo,
}

/// Registered [`typed_dict!`] struct (runtime TypedDicts + stub generation).
#[cfg(feature = "python")]
pub struct TypedDictInfo {
    pub name: &'static str,
    pub doc: &'static str,
    pub type_id: fn() -> std::any::TypeId,
    pub fields: &'static [TypedDictField],
}

#[cfg(feature = "python")]
inventory::collect!(TypedDictInfo);

/// Declare a serde view that is also a Python `TypedDict`.
///
/// Expands to the struct (`Serialize`). With `python`: `IntoPyObject` via
/// `pythonize` (so `#[pymethods]` can return it directly) and an inventory
/// entry that defines the runtime `TypedDict`. With `stubs`: a `PyStubType`
/// naming it, and field types for `class Name(typing.TypedDict)` in the stub.
macro_rules! typed_dict {
    (
        $(#[doc = $doc:literal])*
        pub struct $name:ident {
            $(
                $(#[doc = $fdoc:literal])*
                pub $field:ident : $ty:ty
            ),* $(,)?
        }
    ) => {
        $(#[doc = $doc])*
        #[derive(Clone, Debug, PartialEq, Serialize)]
        pub struct $name {
            $(
                $(#[doc = $fdoc])*
                pub $field: $ty,
            )*
        }

        #[cfg(feature = "stubs")]
        impl pyo3_stub_gen::PyStubType for $name {
            fn type_output() -> pyo3_stub_gen::TypeInfo {
                pyo3_stub_gen::TypeInfo::locally_defined(
                    stringify!($name),
                    pyo3_stub_gen::ModuleRef::Default,
                )
            }
        }

        #[cfg(feature = "python")]
        impl<'py> pyo3::IntoPyObject<'py> for $name {
            type Target = pyo3::PyAny;
            type Output = pyo3::Bound<'py, pyo3::PyAny>;
            type Error = pyo3::PyErr;

            fn into_pyobject(self, py: pyo3::Python<'py>) -> Result<Self::Output, Self::Error> {
                Ok(pythonize::pythonize(py, &self)?)
            }
        }

        #[cfg(feature = "python")]
        inventory::submit! {
            TypedDictInfo {
                name: stringify!($name),
                doc: concat!($($doc, "\n"),*),
                type_id: std::any::TypeId::of::<$name>,
                fields: &[$(
                    TypedDictField {
                        name: stringify!($field),
                        doc: concat!($($fdoc, "\n"),*),
                        #[cfg(feature = "stubs")]
                        r#type: <$ty as pyo3_stub_gen::PyStubType>::type_output,
                    },
                )*],
            }
        }
    };
}

/// Node / edge attribute value (`bool` / `int` / `str`).
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(untagged)]
pub enum AttrValue {
    Bool(bool),
    Int(i64),
    Str(String),
}

impl From<&GraphValue> for AttrValue {
    fn from(v: &GraphValue) -> Self {
        match v {
            GraphValue::Bool(b) => Self::Bool(*b),
            GraphValue::I64(i) => Self::Int(*i),
            GraphValue::String(s) => Self::Str(s.clone()),
        }
    }
}

impl From<AttrValue> for GraphValue {
    fn from(v: AttrValue) -> Self {
        match v {
            AttrValue::Bool(b) => Self::Bool(b),
            AttrValue::Int(i) => Self::I64(i),
            AttrValue::Str(s) => Self::String(s),
        }
    }
}

pub fn attrs_view(attrs: &AttrMap) -> BTreeMap<String, AttrValue> {
    attrs
        .iter()
        .map(|(k, v)| (k.clone(), AttrValue::from(v)))
        .collect()
}

typed_dict! {
    /// One elementary plan step (``rule`` + plan-atom ``site`` labels).
    pub struct PlanStepDict {
        pub rule: String,
        /// Reactant indexes, or ``WillAdd(...)`` / ``AddedBy(...)`` labels.
        pub site: Vec<String>,
    }
}

typed_dict! {
    /// ``PathOutcome.to_dict()``.
    pub struct PathOutcomeDict {
        pub smiles: String,
        pub steps: Vec<PlanStepDict>,
        /// All-zero for exact hits; ``cost > 0`` only from ``find_path_partial``.
        pub residual: ResidualDict,
    }
}

typed_dict! {
    /// ``PathCounters.to_dict()``.
    pub struct PathCountersDict {
        pub nodes: usize,
        pub mol_edits: usize,
        pub expansions: usize,
        pub billed: usize,
        pub dropped_duplicate_plan: usize,
        pub dropped_exact_plan: usize,
        pub dropped_skeleton_twin: usize,
        pub signal_contained_plan: usize,
        pub diversity_repush: usize,
        pub unstable_csmi_key: usize,
        pub timed_out: bool,
    }
}

typed_dict! {
    /// ``Emission.to_dict()``.
    pub struct EmissionDict {
        pub pattern_name: String,
        pub site: usize,
        /// Reactant CSMI.
        pub reactant: String,
        /// Product CSMIs.
        pub products: Vec<String>,
        /// Leaf-first namespace names (``None`` for an unnamed set).
        pub rule_path: Vec<Option<String>>,
    }
}

typed_dict! {
    /// One hop of ``RandomPathOutcome.steps()``.
    pub struct RandomPathStepDict {
        pub rule: String,
        pub pattern: String,
        pub site: Vec<usize>,
        pub products: Vec<String>,
        /// Index into ``products`` the walk continued from.
        pub chosen: usize,
    }
}

typed_dict! {
    /// Pattern summary in ``RandomPathOutcome.patterns()``.
    pub struct PatternSummaryDict {
        pub name: String,
        pub smarts: String,
        pub cleaves: bool,
        pub search_bias: i8,
    }
}

typed_dict! {
    /// ``RandomPathOutcome.to_dict()``.
    pub struct RandomPathOutcomeDict {
        pub smiles: String,
        pub path: Vec<String>,
        pub steps: Vec<RandomPathStepDict>,
        pub patterns: Vec<PatternSummaryDict>,
    }
}

typed_dict! {
    /// One search hop of ``PathOutcome.hops()``.
    pub struct HopDict {
        pub pattern_name: String,
        /// Site on ``reactant``.
        pub site: usize,
        /// Unique-edit orbit (includes ``site``).
        pub site_orbit: Vec<usize>,
        /// CSMI the rule was applied to.
        pub reactant: String,
        /// Kept fragment CSMI.
        pub product: String,
        /// Cleaved-off fragment CSMIs.
        pub sides: Vec<String>,
        /// Leaf rule name (``""`` if unnamed).
        pub rule: String,
        /// Outer→leaf namespace names (``Default`` root omitted).
        pub rule_path: Vec<Option<String>>,
    }
}

typed_dict! {
    /// One cleavage-side bag of ``StepPlan.maybe()``.
    pub struct MaybeEntryDict {
        pub site: Vec<usize>,
        pub side: String,
        pub opens: Vec<Vec<usize>>,
        pub span_sites: Vec<Vec<usize>>,
    }
}

typed_dict! {
    /// ``StepPlan.to_dict()``.
    pub struct StepPlanDict {
        pub steps: Vec<PlanStepDict>,
        pub n_linearizations: usize,
        /// ``(before, after)`` step-index pairs.
        pub precedes: Vec<(usize, usize)>,
        pub maybe: Vec<MaybeEntryDict>,
    }
}

typed_dict! {
    /// Leftover disagreement (``PathOutcome.residual()``, ``MetabolicNetwork.missed()``).
    pub struct ResidualDict {
        pub cost: usize,
        pub n_extra: usize,
        /// Soft category hints (conjugation-sized add, needs oxygen, …).
        pub categories: Vec<String>,
        /// True when no helper candidate resolves this residual.
        pub unresolvable: bool,
    }
}

typed_dict! {
    /// ``GraphNode.to_dict()``.
    pub struct GraphNodeDict {
        pub index: usize,
        pub csmi: String,
        pub sealed: bool,
        pub expanded: bool,
        pub attrs: BTreeMap<String, AttrValue>,
    }
}

typed_dict! {
    /// ``GraphEdge.to_dict()``.
    pub struct GraphEdgeDict {
        pub child_index: usize,
        pub inbound_slot: usize,
        pub parent_index: usize,
        pub rule: String,
        pub pattern_name: String,
        pub site: usize,
        pub cleaves: bool,
        pub attrs: BTreeMap<String, AttrValue>,
    }
}

fn plan_atom_label(atom: &PlanAtom) -> String {
    match atom {
        PlanAtom::Index(i) => i.to_string(),
        PlanAtom::WillAdd { element, at } => format!("WillAdd({element}@{at})"),
        PlanAtom::AddedBy { rule, anchors } => format!("AddedBy({rule:?},{anchors:?})"),
    }
}

impl From<&Step> for PlanStepDict {
    fn from(step: &Step) -> Self {
        Self {
            rule: step.rule.clone(),
            site: step.site.iter().map(plan_atom_label).collect(),
        }
    }
}

pub fn plan_steps(steps: &[Step]) -> Vec<PlanStepDict> {
    steps.iter().map(PlanStepDict::from).collect()
}

impl From<&PathOutcome> for PathOutcomeDict {
    fn from(hit: &PathOutcome) -> Self {
        Self {
            smiles: hit.smiles.clone(),
            steps: plan_steps(hit.plan.steps()),
            residual: ResidualDict::from(&hit.residual),
        }
    }
}

impl From<&Emission> for EmissionDict {
    fn from(e: &Emission) -> Self {
        Self {
            pattern_name: e.pattern_name.clone(),
            site: e.site,
            reactant: e.reactant.csmi().as_ref().to_string(),
            products: e.product_csmis(),
            rule_path: e.rule_path.clone(),
        }
    }
}

impl From<&RandomPathStep> for RandomPathStepDict {
    fn from(s: &RandomPathStep) -> Self {
        Self {
            rule: s.rule.clone(),
            pattern: s.pattern_name.clone(),
            site: s.site.clone(),
            products: s.products.clone(),
            chosen: s.chosen,
        }
    }
}

impl From<&PatternInfo> for PatternSummaryDict {
    fn from(p: &PatternInfo) -> Self {
        Self {
            name: p.name.clone(),
            smarts: p.smarts.clone(),
            cleaves: p.effect.cleaves,
            search_bias: p.search_bias,
        }
    }
}

impl From<&RandomPathOutcome> for RandomPathOutcomeDict {
    fn from(o: &RandomPathOutcome) -> Self {
        Self {
            smiles: o.smiles.clone(),
            path: o.path.clone(),
            steps: o.steps.iter().map(RandomPathStepDict::from).collect(),
            patterns: o.patterns.iter().map(PatternSummaryDict::from).collect(),
        }
    }
}

impl From<&PathCounters> for PathCountersDict {
    fn from(c: &PathCounters) -> Self {
        Self {
            nodes: c.nodes,
            mol_edits: c.mol_edits,
            expansions: c.expansions,
            billed: c.billed(),
            dropped_duplicate_plan: c.dropped_duplicate_plan,
            dropped_exact_plan: c.dropped_exact_plan,
            dropped_skeleton_twin: c.dropped_skeleton_twin,
            signal_contained_plan: c.signal_contained_plan,
            diversity_repush: c.diversity_repush,
            unstable_csmi_key: c.unstable_csmi_key,
            timed_out: c.timed_out,
        }
    }
}

impl From<&PathStep> for HopDict {
    fn from(s: &PathStep) -> Self {
        Self {
            pattern_name: s.pattern_name.clone(),
            site: s.site,
            site_orbit: s.site_orbit.clone(),
            reactant: s.reactant.clone(),
            product: s.product.clone(),
            sides: s.sides.clone(),
            rule: s.leaf_rule().unwrap_or("").to_string(),
            // Outer→leaf for display; omit Default catalog root.
            rule_path: s
                .rule_path
                .iter()
                .rev()
                .filter(|p| p.as_deref() != Some("Default"))
                .cloned()
                .collect(),
        }
    }
}

pub fn hops(steps: &[PathStep]) -> Vec<HopDict> {
    steps.iter().map(HopDict::from).collect()
}

impl From<&CleavageSide> for MaybeEntryDict {
    fn from(e: &CleavageSide) -> Self {
        Self {
            site: e.site.iter().copied().collect(),
            side: e.side.clone(),
            opens: e
                .opens
                .iter()
                .map(|o| o.iter().copied().collect())
                .collect(),
            span_sites: e
                .span_sites()
                .into_iter()
                .map(|s| s.iter().copied().collect())
                .collect(),
        }
    }
}

pub fn maybe_entries(maybe: &Maybe) -> Vec<MaybeEntryDict> {
    maybe.entries.iter().map(MaybeEntryDict::from).collect()
}

impl From<&AtomDiffResidual> for ResidualDict {
    fn from(r: &AtomDiffResidual) -> Self {
        Self {
            cost: r.cost,
            n_extra: r.n_extra,
            categories: r.categories.clone(),
            unresolvable: r.unresolvable,
        }
    }
}
