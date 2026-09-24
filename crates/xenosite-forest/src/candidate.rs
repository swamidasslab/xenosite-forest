//! Deferred site hits: site + pattern + parent, materialize on demand.
//!
//! Calling a [`crate::ruleset::RuleSet`] can iterate these triples instead of
//! products. A search reads [`PatternInfo`] / [`Effect`] to filter, then calls
//! [`Candidate::materialize`] only for survivors — no filter closures inside
//! the rule walk.

use std::collections::BTreeMap;

use crate::ForestError;
use crate::mol::Molecule;
use crate::pattern::{Edit, Emission, PatternInfo};
use crate::ruleset::apply_edit_for_candidate;

/// Which mol to apply the edit on.
///
/// [`ParentRef::Context`] is the aromatic (or input) parent — ResonanceRule
/// SMIRKS still go through Kekulé `reactant_parent` at materialize time.
/// [`ParentRef::Form`] is an exclusive-double match already on a Kekulé writing
/// (Epoxidation); apply runs on that form directly.
#[derive(Clone)]
pub enum ParentRef {
    Context,
    Form(Box<Molecule>),
}

impl std::fmt::Debug for ParentRef {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Context => f.write_str("Context"),
            Self::Form(mol) => write!(f, "Form(atoms={})", mol.atom_count()),
        }
    }
}

/// One discovered site that has not yet been edited.
#[derive(Clone, Debug)]
pub struct Candidate {
    pub site: usize,
    pub pattern: PatternInfo,
    /// Leaf-first rule namespace (emitting set, then containers).
    pub rule_path: Vec<Option<String>>,
    pub mapped: BTreeMap<u16, usize>,
    pub parent: ParentRef,
}

impl Candidate {
    pub fn pattern_name(&self) -> &str {
        &self.pattern.name
    }

    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }

    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    /// Run the edit and return product CSMIs. Call only after filtering.
    pub fn materialize(&self, context: &Molecule) -> Result<Vec<String>, ForestError> {
        let work = match &self.parent {
            ParentRef::Context => context,
            ParentRef::Form(form) => form,
        };
        apply_edit_for_candidate(work, &self.pattern, &self.mapped)
    }

    /// Materialize into an [`Emission`] (site + namespace + products).
    pub fn emit(&self, context: &Molecule) -> Result<Option<Emission>, ForestError> {
        let products = self.materialize(context)?;
        if products.is_empty() {
            return Ok(None);
        }
        Ok(Some(Emission {
            site: self.site,
            pattern_name: self.pattern.name.clone(),
            rule_path: self.rule_path.clone(),
            products,
        }))
    }

    pub fn is_pair_endpoint(&self) -> bool {
        matches!(self.pattern.edit, Edit::PairEndpoint(_))
    }
}
