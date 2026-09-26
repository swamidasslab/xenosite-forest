//! Load config; name a mapped reactant→product.

use crate::assignment::Assignments;
use crate::bundle::AnnotationBundle;
use crate::chemistry::{MappedReaction, ReactionChemistry};
use crate::error::{Error, Result};
use crate::skos::Thesaurus;
use crate::sssom::SssomTable;
use crate::term::{SiteRef, Term};
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
    ///
    /// Terms are site-localized when assignment hits carry a [`SiteRef`]. The same
    /// concept may appear once per distinct site so multi-change cases disambiguate.
    pub fn name(&self, query: &MappedReaction) -> Result<Vec<Term>> {
        let chem = ReactionChemistry::prepare(query)?;
        let hits = self.assignments.matching(&chem)?;

        // Key: (concept id, site) → evidence rule ids
        let mut evidence: BTreeMap<(String, SiteRef), Vec<String>> = BTreeMap::new();
        let mut keys: BTreeSet<(String, SiteRef)> = BTreeSet::new();

        // Companion spines reached via skos:relatedMatch (data), not ad-hoc emit lists.
        const RELATED_EMIT_SPINES: &[&str] = &[
            "xrm:1500000", // reactive metabolite / product class
            "xrm:1400000", // medchem liability
            "xrm:2200000", // leaving group
        ];

        for hit in hits {
            for id in &hit.rule.emit {
                if self.thesaurus.get(id).is_none() {
                    return Err(Error::UnknownConcept(format!("{id} ({})", hit.rule.id)));
                }
                let key = (id.clone(), hit.site.clone());
                evidence
                    .entry(key.clone())
                    .or_default()
                    .push(hit.rule.id.clone());
                keys.insert(key);
                for anc in self.thesaurus.ancestors(id) {
                    let akey = (anc.as_str().to_string(), hit.site.clone());
                    keys.insert(akey);
                }
                // General gap coverage: follow relatedMatch into companion spines.
                for rel in self
                    .thesaurus
                    .related_under_spines(id, RELATED_EMIT_SPINES)
                {
                    let rkey = (rel.as_str().to_string(), hit.site.clone());
                    evidence
                        .entry(rkey.clone())
                        .or_default()
                        .push(hit.rule.id.clone());
                    keys.insert(rkey);
                }
            }
        }

        let mut terms = Vec::with_capacity(keys.len());
        for (id, site) in keys {
            let ev = evidence.get(&(id.clone(), site.clone())).cloned().unwrap_or_default();
            let mut term = self.thesaurus.term_shell(&id, ev)?;
            term.site = site;
            let (intra, inter) = self.mappings.links_for(&id, &self.scheme_prefix);
            term.intra_matches.extend(intra);
            term.inter_matches.extend(inter);
            terms.push(term);
        }
        terms.sort_by(|a, b| {
            b.specificity
                .depth
                .cmp(&a.specificity.depth)
                .then_with(|| a.site.cmp(&b.site))
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

    /// Facet flat terms into per-site annotation bundles (`xmet:*` combination layer).
    ///
    /// Combinatorial display names live in [`AnnotationBundle::site_label`], not as
    /// ontology concepts.
    pub fn annotate(&self, query: &MappedReaction) -> Result<Vec<AnnotationBundle>> {
        let terms = self.name(query)?;
        Ok(AnnotationBundle::from_terms(&self.thesaurus, &terms))
    }

    pub fn annotate_smiles(
        &self,
        reactant: &str,
        product: &str,
        tags: &[&str],
    ) -> Result<Vec<AnnotationBundle>> {
        self.annotate(&MappedReaction::new(reactant, product).with_tags(tags.iter().copied()))
    }

    /// Compact lines for feedback dumps: `prefLabel [id]@site ← path`.
    pub fn format_sample_lines(&self, terms: &[Term]) -> Vec<String> {
        terms
            .iter()
            .map(|t| {
                format!(
                    "{} [{}]{} ← {}",
                    t.pref_label,
                    t.id,
                    t.site.display_suffix(),
                    t.path_labels.join(" > ")
                )
            })
            .collect()
    }
}
