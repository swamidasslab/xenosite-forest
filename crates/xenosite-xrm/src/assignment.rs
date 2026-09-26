//! JSONL assignment rules. Match = all declared constraints pass.

use crate::chemistry::ReactionChemistry;
use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// One JSONL assignment rule. Only fields you set are checked.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AssignmentRule {
    pub id: String,
    pub emit: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reactant_smarts: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub product_smarts: Option<String>,
    #[serde(default)]
    pub require_map: Vec<u16>,
    #[serde(default)]
    pub delta: BTreeMap<String, i32>,
    /// Opaque tags; match if any is present (ignored when empty).
    #[serde(default)]
    pub tags_any: Vec<String>,
    /// If set, first reactant SMARTS hit must be aromatic / not.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site_aromatic: Option<bool>,
}

#[derive(Clone, Debug, Default)]
pub struct Assignments {
    pub rules: Vec<AssignmentRule>,
}

impl Assignments {
    pub fn load_jsonl(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_jsonl_str(&fs::read_to_string(path.as_ref())?)
    }

    pub fn from_jsonl_str(text: &str) -> Result<Self> {
        let mut rules = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            let rule: AssignmentRule = serde_json::from_str(t)
                .map_err(|e| Error::Config(format!("assignment line {}: {e}", i + 1)))?;
            if rule.emit.is_empty() {
                return Err(Error::Config(format!("{}: empty emit", rule.id)));
            }
            if !rule.has_criterion() {
                return Err(Error::Config(format!(
                    "{}: need at least one constraint",
                    rule.id
                )));
            }
            rules.push(rule);
        }
        Ok(Self { rules })
    }

    pub fn merge(&mut self, other: Assignments) {
        self.rules.extend(other.rules);
    }

    pub fn matching<'a>(&'a self, chem: &ReactionChemistry) -> Result<Vec<&'a AssignmentRule>> {
        let mut hits = Vec::new();
        for rule in &self.rules {
            if rule.matches(chem)? {
                hits.push(rule);
            }
        }
        Ok(hits)
    }
}

impl AssignmentRule {
    fn has_criterion(&self) -> bool {
        self.reactant_smarts.is_some()
            || self.product_smarts.is_some()
            || !self.delta.is_empty()
            || !self.tags_any.is_empty()
            || self.site_aromatic.is_some()
    }

    fn matches(&self, chem: &ReactionChemistry) -> Result<bool> {
        if !self.tags_any.is_empty() && !self.tags_any.iter().any(|t| chem.tags.contains(t)) {
            return Ok(false);
        }
        for (k, want) in &self.delta {
            if chem.delta.get(k).copied().unwrap_or(0) != *want {
                return Ok(false);
            }
        }
        if let Some(s) = &self.reactant_smarts {
            if !chem.smarts_hits(false, s, &self.require_map)? {
                return Ok(false);
            }
        }
        if let Some(s) = &self.product_smarts {
            if !chem.smarts_hits(true, s, &self.require_map)? {
                return Ok(false);
            }
        }
        if let Some(want) = self.site_aromatic {
            let Some(probe) = &self.reactant_smarts else {
                return Err(Error::Config(format!(
                    "{}: site_aromatic needs reactant_smarts",
                    self.id
                )));
            };
            match chem.any_mapped_reactant_aromatic(probe)? {
                Some(is_ar) if is_ar == want => {}
                _ => return Ok(false),
            }
        }
        Ok(true)
    }
}
