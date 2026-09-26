//! JSONL assignment rules. Match = all declared constraints pass.
//! Hits are site-localized when SMARTS or `@map` tags provide a site.

use crate::chemistry::ReactionChemistry;
use crate::error::{Error, Result};
use crate::term::SiteRef;
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
    /// Opaque tags (base names); match if any is present (ignored when empty).
    /// Callers may localize with `@map` suffixes; bases are compared without `@…`.
    #[serde(default)]
    pub tags_any: Vec<String>,
    /// If set, first reactant SMARTS hit must be aromatic / not.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site_aromatic: Option<bool>,
}

/// One rule hit, optionally localized to a site.
#[derive(Clone, Debug)]
pub struct AssignmentHit<'a> {
    pub rule: &'a AssignmentRule,
    pub site: SiteRef,
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

    pub fn matching<'a>(&'a self, chem: &ReactionChemistry) -> Result<Vec<AssignmentHit<'a>>> {
        let mut hits = Vec::new();
        for rule in &self.rules {
            hits.extend(rule.hits(chem)?);
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

    fn hits<'a>(&'a self, chem: &ReactionChemistry) -> Result<Vec<AssignmentHit<'a>>> {
        // Tags: base-match; collect localized sites from caller tags.
        let mut tag_sites = Vec::new();
        if !self.tags_any.is_empty() {
            tag_sites = chem.tag_sites(&self.tags_any);
            if tag_sites.is_empty() {
                return Ok(Vec::new());
            }
        }

        for (k, want) in &self.delta {
            if chem.delta.get(k).copied().unwrap_or(0) != *want {
                return Ok(Vec::new());
            }
        }

        let mut smarts_sites = Vec::new();
        let has_smarts = self.reactant_smarts.is_some() || self.product_smarts.is_some();
        if has_smarts {
            if let Some(s) = &self.reactant_smarts {
                let sites = chem.smarts_sites(false, s, &self.require_map)?;
                if sites.is_empty() {
                    return Ok(Vec::new());
                }
                smarts_sites = sites;
            }
            if let Some(s) = &self.product_smarts {
                let psites = chem.smarts_sites(true, s, &self.require_map)?;
                if psites.is_empty() {
                    return Ok(Vec::new());
                }
                if smarts_sites.is_empty() {
                    smarts_sites = psites;
                } else {
                    // Keep reactant sites that share a map number with a product hit,
                    // or keep all reactant sites when maps are absent.
                    let product_maps: BTreeMap<(), ()> = BTreeMap::new();
                    let _ = product_maps;
                    let pmaps: std::collections::BTreeSet<u16> = psites
                        .iter()
                        .flat_map(|s| s.map_nums.iter().copied())
                        .collect();
                    if !pmaps.is_empty() {
                        smarts_sites.retain(|s| {
                            s.map_nums.is_empty() || s.map_nums.iter().any(|m| pmaps.contains(m))
                        });
                        if smarts_sites.is_empty() {
                            return Ok(Vec::new());
                        }
                    }
                }
            }
            if let Some(want) = self.site_aromatic {
                let Some(probe) = &self.reactant_smarts else {
                    return Err(Error::Config(format!(
                        "{}: site_aromatic needs reactant_smarts",
                        self.id
                    )));
                };
                // Filter sites by aromaticity of first reactant atom in each site.
                let mut kept = Vec::new();
                for site in smarts_sites {
                    let is_ar = if let Some(&a) = site.reactant_atoms.first() {
                        chem.reactant.atom(chematic::core::AtomIdx(a as u32)).aromatic
                    } else {
                        match chem.any_mapped_reactant_aromatic(probe)? {
                            Some(v) => v,
                            None => continue,
                        }
                    };
                    if is_ar == want {
                        kept.push(site);
                    }
                }
                if kept.is_empty() {
                    return Ok(Vec::new());
                }
                smarts_sites = kept;
            }
        } else if self.site_aromatic.is_some() {
            return Err(Error::Config(format!(
                "{}: site_aromatic needs reactant_smarts",
                self.id
            )));
        }

        // Prefer SMARTS sites when present (structure over reaction-tool tags).
        let sites = if !smarts_sites.is_empty() {
            // If tags also localized, intersect by map number when both have maps.
            if !tag_sites.is_empty()
                && tag_sites.iter().any(|s| !s.map_nums.is_empty())
                && smarts_sites.iter().any(|s| !s.map_nums.is_empty())
            {
                let tag_maps: std::collections::BTreeSet<u16> = tag_sites
                    .iter()
                    .flat_map(|s| s.map_nums.iter().copied())
                    .collect();
                smarts_sites
                    .into_iter()
                    .filter(|s| s.map_nums.iter().any(|m| tag_maps.contains(m)))
                    .collect()
            } else {
                smarts_sites
            }
        } else if !tag_sites.is_empty() {
            tag_sites
        } else {
            // Delta-only / unconstrained molecule-level hit.
            vec![SiteRef::default()]
        };

        if sites.is_empty() {
            return Ok(Vec::new());
        }
        Ok(sites
            .into_iter()
            .map(|site| AssignmentHit { rule: self, site })
            .collect())
    }
}
