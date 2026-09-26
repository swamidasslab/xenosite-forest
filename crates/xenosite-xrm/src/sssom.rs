//! SSSOM (Simple Standard for Sharing Ontology Mappings) TSV loader.

use crate::curie::Curie;
use crate::error::{Error, Result};
use crate::term::TermLink;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// SSSOM mapping predicate (SKOS / SSSOM vocabulary).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum MappingPredicate {
    ExactMatch,
    CloseMatch,
    BroadMatch,
    NarrowMatch,
    RelatedMatch,
    Other(String),
}

impl MappingPredicate {
    pub fn parse(s: &str) -> Self {
        let local = s.rsplit_once(':').map(|(_, l)| l).unwrap_or(s);
        match local {
            "exactMatch" | "skosExactMatch" => Self::ExactMatch,
            "closeMatch" | "skosCloseMatch" => Self::CloseMatch,
            "broadMatch" | "skosBroadMatch" => Self::BroadMatch,
            "narrowMatch" | "skosNarrowMatch" => Self::NarrowMatch,
            "relatedMatch" | "skosRelatedMatch" => Self::RelatedMatch,
            other => Self::Other(other.to_string()),
        }
    }

    pub fn skos_name(&self) -> &str {
        match self {
            Self::ExactMatch => "exactMatch",
            Self::CloseMatch => "closeMatch",
            Self::BroadMatch => "broadMatch",
            Self::NarrowMatch => "narrowMatch",
            Self::RelatedMatch => "relatedMatch",
            Self::Other(s) => s.as_str(),
        }
    }
}

/// One SSSOM mapping row.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SssomMapping {
    pub subject_id: Curie,
    pub predicate_id: MappingPredicate,
    pub object_id: Curie,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mapping_justification: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subject_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub object_label: Option<String>,
}

/// Index of SSSOM mappings by subject.
#[derive(Clone, Debug, Default)]
pub struct SssomTable {
    pub by_subject: BTreeMap<String, Vec<SssomMapping>>,
}

impl SssomTable {
    pub fn load_tsv(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())?;
        Self::from_tsv_str(&text)
    }

    pub fn from_tsv_str(text: &str) -> Result<Self> {
        let mut lines = text.lines().filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        });
        let header = lines
            .next()
            .ok_or_else(|| Error::Config("SSSOM TSV empty".into()))?;
        let cols: Vec<&str> = header.split('\t').collect();
        let idx = |name: &str| cols.iter().position(|c| *c == name);
        let i_sub = idx("subject_id").ok_or_else(|| Error::Config("SSSOM missing subject_id".into()))?;
        let i_pred = idx("predicate_id")
            .ok_or_else(|| Error::Config("SSSOM missing predicate_id".into()))?;
        let i_obj = idx("object_id").ok_or_else(|| Error::Config("SSSOM missing object_id".into()))?;
        let i_just = idx("mapping_justification");
        let i_slab = idx("subject_label");
        let i_olab = idx("object_label");

        let mut by_subject: BTreeMap<String, Vec<SssomMapping>> = BTreeMap::new();
        for (lineno, line) in lines.enumerate() {
            let fields: Vec<&str> = line.split('\t').collect();
            let get = |i: usize| -> Result<&str> {
                fields.get(i).copied().ok_or_else(|| {
                    Error::Config(format!("SSSOM row {} short (need col {i})", lineno + 2))
                })
            };
            let mapping = SssomMapping {
                subject_id: Curie::new(get(i_sub)?),
                predicate_id: MappingPredicate::parse(get(i_pred)?),
                object_id: Curie::new(get(i_obj)?),
                mapping_justification: i_just.and_then(|i| fields.get(i).map(|s| s.to_string())),
                subject_label: i_slab.and_then(|i| {
                    fields
                        .get(i)
                        .filter(|s| !s.is_empty())
                        .map(|s| s.to_string())
                }),
                object_label: i_olab.and_then(|i| {
                    fields
                        .get(i)
                        .filter(|s| !s.is_empty())
                        .map(|s| s.to_string())
                }),
            };
            by_subject
                .entry(mapping.subject_id.as_str().to_string())
                .or_default()
                .push(mapping);
        }
        Ok(Self { by_subject })
    }

    pub fn merge(&mut self, other: SssomTable) {
        for (k, mut v) in other.by_subject {
            self.by_subject.entry(k).or_default().append(&mut v);
        }
    }

    pub fn links_for(&self, subject: &str, scheme_prefix: &str) -> (Vec<TermLink>, Vec<TermLink>) {
        let mut intra = Vec::new();
        let mut inter = Vec::new();
        let Some(rows) = self.by_subject.get(subject) else {
            return (intra, inter);
        };
        for row in rows {
            let link = TermLink {
                target_ontology: row.object_id.prefix().map(str::to_string),
                target: row.object_id.clone(),
                predicate: row.predicate_id.skos_name().to_string(),
            };
            let obj_prefix = row.object_id.prefix().unwrap_or("");
            if obj_prefix == scheme_prefix {
                intra.push(link);
            } else {
                inter.push(link);
            }
        }
        (intra, inter)
    }
}
