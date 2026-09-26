//! Rich named terms returned by the namer.

use crate::curie::Curie;
use serde::{Deserialize, Serialize};

/// Which concept scheme / ontology a term belongs to.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OntologyRef {
    /// Scheme id (e.g. `xrm:scheme`).
    pub id: Curie,
    /// Human label (e.g. "Xenobiotic Reaction Metabolism").
    pub label: String,
    /// Optional IRI for the scheme.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub iri: Option<String>,
}

/// How specific a term is relative to the thesaurus root(s).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Specificity {
    /// Depth from a scheme root via `skos:broader` (0 = top concept).
    pub depth: u32,
    /// Number of `skos:narrower` children in the loaded thesaurus.
    pub child_count: usize,
    /// True when the term has no narrower concepts in this thesaurus.
    pub is_leaf: bool,
}

/// A typed link to another concept (intra- or inter-ontology).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TermLink {
    pub target: Curie,
    /// SKOS predicate local name: `broader`, `narrower`, `exactMatch`, …
    pub predicate: String,
    /// Ontology prefix of the target, when parseable.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_ontology: Option<String>,
}

/// Localized site of metabolism on the mapped reactant/product.
///
/// Prefer atom-map numbers from mapped SMILES (`[C:1]`). When maps are absent,
/// reactant atom indices are retained so multi-change cases can still be split.
/// An empty [`SiteRef`] means the term is molecule-level (not site-localized).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct SiteRef {
    /// Stable atom-map numbers from the mapped SMILES (sorted, unique).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub map_nums: Vec<u16>,
    /// Reactant atom indices for the matched site (sorted, unique).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reactant_atoms: Vec<usize>,
}

impl SiteRef {
    pub fn is_empty(&self) -> bool {
        self.map_nums.is_empty() && self.reactant_atoms.is_empty()
    }

    pub fn from_map_nums(nums: impl IntoIterator<Item = u16>) -> Self {
        let mut map_nums: Vec<u16> = nums.into_iter().filter(|n| *n != 0).collect();
        map_nums.sort_unstable();
        map_nums.dedup();
        Self {
            map_nums,
            reactant_atoms: Vec::new(),
        }
    }

    /// Compact display: `@1,2` (map nums) or `@a3,a5` (atom indices) or empty.
    pub fn display_suffix(&self) -> String {
        if !self.map_nums.is_empty() {
            format!(
                "@{}",
                self.map_nums
                    .iter()
                    .map(|n| n.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        } else if !self.reactant_atoms.is_empty() {
            format!(
                "@{}",
                self.reactant_atoms
                    .iter()
                    .map(|a| format!("a{a}"))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        } else {
            String::new()
        }
    }
}

/// One systematic name hit for a mapped reaction (optionally site-localized).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Term {
    pub id: Curie,
    pub ontology: OntologyRef,
    pub pref_label: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub alt_labels: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub definition: Option<String>,
    pub specificity: Specificity,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub broader: Vec<TermLink>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub narrower: Vec<TermLink>,
    /// Intra-ontology synonyms (`skos:altLabel` already listed; exact/close within scheme).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub intra_matches: Vec<TermLink>,
    /// Inter-ontology mappings (MeSH, MOP, opaque `forest.rule:…`, …).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inter_matches: Vec<TermLink>,
    /// Root→leaf prefLabel path for this concept in its scheme.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub path_labels: Vec<String>,
    /// Assignment rule ids from config that emitted this term.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub evidence: Vec<String>,
    /// Site localization for this emit; empty = molecule-level.
    #[serde(default, skip_serializing_if = "SiteRef::is_empty")]
    pub site: SiteRef,
}

impl Term {
    /// True when `self` is a descendant of `other` via broader links.
    pub fn is_more_specific_than(&self, other: &Term) -> bool {
        self.specificity.depth > other.specificity.depth
            && self
                .path_labels
                .windows(other.path_labels.len())
                .any(|w| w == other.path_labels.as_slice())
            || (self.path_labels.len() > other.path_labels.len()
                && self.path_labels[..other.path_labels.len()] == other.path_labels[..])
    }
}
