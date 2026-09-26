//! Namer: load config, name a mapped reactant→product.

use crate::assignment::Assignments;
use crate::chemistry::{MappedReaction, ReactionChemistry};
use crate::error::{Error, Result};
use crate::skos::Thesaurus;
use crate::sssom::SssomTable;
use crate::term::Term;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

/// Manifest listing ontology / mapping / assignment files (paths relative to manifest dir).
#[derive(Clone, Debug, Deserialize)]
pub struct NamerConfig {
    #[serde(default)]
    pub ontologies: Vec<String>,
    #[serde(default)]
    pub mappings: Vec<String>,
    #[serde(default)]
    pub assignments: Vec<String>,
}

/// Config-driven reaction namer.
#[derive(Clone, Debug)]
pub struct Namer {
    pub thesaurus: Thesaurus,
    pub mappings: SssomTable,
    pub assignments: Assignments,
    /// Scheme prefix used for intra- vs inter-ontology classification (e.g. `xrm`).
    pub scheme_prefix: String,
}

impl Namer {
    pub fn from_manifest(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let base = path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let text = fs::read_to_string(path)?;
        let cfg: NamerConfig = serde_json::from_str(&text)?;
        Self::from_config(&cfg, &base)
    }

    pub fn from_config(cfg: &NamerConfig, base: &Path) -> Result<Self> {
        if cfg.ontologies.is_empty() {
            return Err(Error::Config("manifest.ontologies is empty".into()));
        }
        let mut thesaurus = Thesaurus::default();
        for rel in &cfg.ontologies {
            let t = Thesaurus::load_jsonld(base.join(rel))?;
            thesaurus.merge(t);
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

    /// Name a mapped reactant→product. Returns rich terms, most specific first.
    pub fn name(&self, query: &MappedReaction) -> Result<Vec<Term>> {
        let chem = ReactionChemistry::prepare(query)?;
        let hits = self.assignments.matching(&chem)?;
        let mut evidence: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut order: Vec<String> = Vec::new();
        let mut seen = BTreeSet::new();

        for rule in hits {
            let mut emit_ids = rule.emit.clone();
            if rule.include_ancestors {
                for id in rule.emit.clone() {
                    for anc in self.thesaurus.ancestors(&id) {
                        emit_ids.push(anc.as_str().to_string());
                    }
                }
            }
            for id in emit_ids {
                if self.thesaurus.get(&id).is_none() {
                    return Err(Error::UnknownConcept(format!(
                        "{id} (from assignment {})",
                        rule.id
                    )));
                }
                evidence.entry(id.clone()).or_default().push(rule.id.clone());
                if seen.insert(id.clone()) {
                    order.push(id);
                }
            }
        }

        let mut terms = Vec::with_capacity(order.len());
        for id in &order {
            let mut term = self
                .thesaurus
                .term_shell(id, evidence.get(id).cloned().unwrap_or_default())?;
            let (intra, inter) = self.mappings.links_for(id, &self.scheme_prefix);
            // Merge SKOS-declared matches already on the shell with SSSOM rows.
            term.intra_matches.extend(intra);
            term.inter_matches.extend(inter);
            // Drop enzyme-facet primary labels: never emit terms whose prefLabel
            // looks like an enzyme family (config should not put them in chemical
            // emit paths; this is a hard guard).
            if looks_like_enzyme_label(&term.pref_label) {
                continue;
            }
            terms.push(term);
        }

        // Most specific first (deeper path), then pref_label.
        terms.sort_by(|a, b| {
            b.specificity
                .depth
                .cmp(&a.specificity.depth)
                .then_with(|| a.pref_label.cmp(&b.pref_label))
        });
        Ok(terms)
    }

    /// Convenience: SMILES pair + optional opaque tags.
    pub fn name_smiles(
        &self,
        reactant: &str,
        product: &str,
        tags: &[&str],
    ) -> Result<Vec<Term>> {
        let q = MappedReaction::new(reactant, product).with_tags(tags.iter().copied());
        self.name(&q)
    }
}

fn looks_like_enzyme_label(label: &str) -> bool {
    let l = label.to_ascii_lowercase();
    l.contains("cytochrome")
        || l.contains("cyp")
        || l.starts_with("ugt")
        || l.contains("transferase enzyme")
        || l.contains("ec ")
}
