//! JSONL structural assignment rules (config only; no forest imports).

use crate::chemistry::ReactionChemistry;
use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// One JSONL assignment rule.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AssignmentRule {
    pub id: String,
    /// Concept CURIEs to emit when this rule matches.
    pub emit: Vec<String>,
    /// Also emit every `skos:broader` ancestor of each emit id.
    #[serde(default = "default_true")]
    pub include_ancestors: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reactant_smarts: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub product_smarts: Option<String>,
    /// Atom-map numbers that must appear in the SMARTS hits.
    #[serde(default)]
    pub require_map: Vec<u16>,
    /// Required formula delta (product − reactant), e.g. `{"O": 1}`.
    #[serde(default)]
    pub delta: BTreeMap<String, i32>,
    /// All of these opaque tags must be present on the query.
    #[serde(default)]
    pub tags_all: Vec<String>,
    /// At least one of these opaque tags must be present (ignored if empty).
    #[serde(default)]
    pub tags_any: Vec<String>,
    /// When set, require a reactant SMARTS hit whose first mapped atom is aromatic / not.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site_aromatic: Option<bool>,
    /// SMARTS used for the aromaticity probe (defaults to reactant_smarts).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aromatic_probe_smarts: Option<String>,
    /// Priority: higher wins when sorting emitted evidence (specificity tie-break).
    #[serde(default)]
    pub priority: i32,
}

fn default_true() -> bool {
    true
}

#[derive(Clone, Debug, Default)]
pub struct Assignments {
    pub rules: Vec<AssignmentRule>,
}

impl Assignments {
    pub fn load_jsonl(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())?;
        Self::from_jsonl_str(&text)
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
                return Err(Error::Config(format!(
                    "assignment {} has empty emit",
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

    /// Matching rules, highest priority first.
    pub fn matching<'a>(&'a self, chem: &ReactionChemistry) -> Result<Vec<&'a AssignmentRule>> {
        let mut hits = Vec::new();
        for rule in &self.rules {
            if self.rule_matches(rule, chem)? {
                hits.push(rule);
            }
        }
        hits.sort_by(|a, b| b.priority.cmp(&a.priority).then_with(|| a.id.cmp(&b.id)));
        Ok(hits)
    }

    fn rule_matches(&self, rule: &AssignmentRule, chem: &ReactionChemistry) -> Result<bool> {
        if !rule.tags_all.is_empty() && !rule.tags_all.iter().all(|t| chem.tags.contains(t)) {
            return Ok(false);
        }
        if !rule.tags_any.is_empty() && !rule.tags_any.iter().any(|t| chem.tags.contains(t)) {
            return Ok(false);
        }
        for (k, want) in &rule.delta {
            if chem.delta.get(k).copied().unwrap_or(0) != *want {
                return Ok(false);
            }
        }
        if let Some(s) = &rule.reactant_smarts {
            if !chem.smarts_hits(false, s, &rule.require_map)? {
                return Ok(false);
            }
        }
        if let Some(s) = &rule.product_smarts {
            if !chem.smarts_hits(true, s, &rule.require_map)? {
                return Ok(false);
            }
        }
        if let Some(want) = rule.site_aromatic {
            let probe = rule
                .aromatic_probe_smarts
                .as_deref()
                .or(rule.reactant_smarts.as_deref())
                .ok_or_else(|| {
                    Error::Config(format!(
                        "assignment {} sets site_aromatic without SMARTS probe",
                        rule.id
                    ))
                })?;
            match chem.any_mapped_reactant_aromatic(probe)? {
                Some(is_ar) if is_ar == want => {}
                _ => return Ok(false),
            }
        }
        // A rule with only empty constraints would match everything — require
        // at least one positive criterion.
        let has_criterion = rule.reactant_smarts.is_some()
            || rule.product_smarts.is_some()
            || !rule.delta.is_empty()
            || !rule.tags_all.is_empty()
            || !rule.tags_any.is_empty()
            || rule.site_aromatic.is_some();
        Ok(has_criterion)
    }
}
