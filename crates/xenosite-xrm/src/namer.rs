//! Load config; name a mapped reactant→product.

use crate::assignment::Assignments;
use crate::chemistry::{MappedReaction, ReactionChemistry};
use crate::error::{Error, Result};
use crate::skos::Thesaurus;
use crate::sssom::SssomTable;
use crate::term::Term;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

#[derive(Clone, Debug, Deserialize)]
pub struct NamerConfig {
    #[serde(default)]
    pub ontologies: Vec<String>,
    #[serde(default)]
    pub mappings: Vec<String>,
    #[serde(default)]
    pub assignments: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct Namer {
    pub thesaurus: Thesaurus,
    pub mappings: SssomTable,
    pub assignments: Assignments,
    pub scheme_prefix: String,
}

impl Namer {
    pub fn from_manifest(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let base = path.parent().unwrap_or_else(|| Path::new("."));
        let cfg: NamerConfig = serde_json::from_str(&fs::read_to_string(path)?)?;
        Self::from_config(&cfg, base)
    }

    pub fn from_config(cfg: &NamerConfig, base: &Path) -> Result<Self> {
        if cfg.ontologies.is_empty() {
            return Err(Error::Config("manifest.ontologies is empty".into()));
        }
        let mut thesaurus = Thesaurus::default();
        for rel in &cfg.ontologies {
            thesaurus.merge(Thesaurus::load_jsonld(base.join(rel))?);
        }
        let mut mappings = SssomTable::default();
        for rel in &cfg.mappings {
            mappings.merge(SssomTable::load_tsv(base.join(rel))?);
        }
        let mut assignments = Assignments::default();
        for rel in &cfg.assignments {
            assignments.merge(Assignments::load_jsonl(base.join(rel))?);
        }
        thesaurus.validate()?;
        let scheme_prefix = thesaurus
            .schemes
            .values()
            .next()
            .and_then(|s| s.id.prefix())
            .unwrap_or("xrm")
            .to_string();
        Ok(Self {
            thesaurus,
            mappings,
            assignments,
            scheme_prefix,
        })
    }

    /// Name a mapped reactant→product. Most specific terms first; ancestors included.
    pub fn name(&self, query: &MappedReaction) -> Result<Vec<Term>> {
        let chem = ReactionChemistry::prepare(query)?;
        let hits = self.assignments.matching(&chem)?;

        let mut evidence: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut ids = BTreeSet::new();
        for rule in hits {
            for id in &rule.emit {
                if self.thesaurus.get(id).is_none() {
                    return Err(Error::UnknownConcept(format!("{id} ({})", rule.id)));
                }
                evidence.entry(id.clone()).or_default().push(rule.id.clone());
                ids.insert(id.clone());
                for anc in self.thesaurus.ancestors(id) {
                    ids.insert(anc.as_str().to_string());
                }
            }
        }

        let mut terms = Vec::with_capacity(ids.len());
        for id in ids {
            let mut term = self
                .thesaurus
                .term_shell(&id, evidence.get(&id).cloned().unwrap_or_default())?;
            let (intra, inter) = self.mappings.links_for(&id, &self.scheme_prefix);
            term.intra_matches.extend(intra);
            term.inter_matches.extend(inter);
            terms.push(term);
        }
        terms.sort_by(|a, b| {
            b.specificity
                .depth
                .cmp(&a.specificity.depth)
                .then_with(|| a.pref_label.cmp(&b.pref_label))
        });
        Ok(terms)
    }

    pub fn name_smiles(
        &self,
        reactant: &str,
        product: &str,
        tags: &[&str],
    ) -> Result<Vec<Term>> {
        self.name(&MappedReaction::new(reactant, product).with_tags(tags.iter().copied()))
    }

    /// Compact lines for feedback dumps: `prefLabel [id] ← path`.
    pub fn format_sample_lines(&self, terms: &[Term]) -> Vec<String> {
        terms
            .iter()
            .map(|t| {
                format!(
                    "{} [{}] ← {}",
                    t.pref_label,
                    t.id,
                    t.path_labels.join(" > ")
                )
            })
            .collect()
    }
}
