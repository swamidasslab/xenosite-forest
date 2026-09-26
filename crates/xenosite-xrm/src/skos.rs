//! SKOS thesaurus load (JSON-LD subset).

use crate::curie::Curie;
use crate::error::{Error, Result};
use crate::term::{OntologyRef, Specificity, Term, TermLink};
use serde::Deserialize;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

/// One SKOS concept as stored after load.
#[derive(Clone, Debug)]
pub struct SkosConcept {
    pub id: Curie,
    pub pref_label: String,
    pub alt_labels: Vec<String>,
    pub definition: Option<String>,
    pub in_scheme: Curie,
    pub broader: Vec<Curie>,
    pub narrower: Vec<Curie>,
    pub exact_match: Vec<Curie>,
    pub close_match: Vec<Curie>,
    pub broad_match: Vec<Curie>,
    pub related_match: Vec<Curie>,
    pub top_concept: bool,
}

/// Concept scheme metadata.
#[derive(Clone, Debug)]
pub struct ConceptScheme {
    pub id: Curie,
    pub label: String,
    pub iri: Option<String>,
}

/// In-memory SKOS thesaurus.
#[derive(Clone, Debug, Default)]
pub struct Thesaurus {
    pub schemes: BTreeMap<String, ConceptScheme>,
    pub concepts: BTreeMap<String, SkosConcept>,
}

impl Thesaurus {
    pub fn load_jsonld(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path.as_ref())?;
        Self::from_jsonld_str(&text)
    }

    pub fn from_jsonld_str(text: &str) -> Result<Self> {
        let root: Value = serde_json::from_str(text)?;
        let graph = root
            .get("@graph")
            .and_then(|g| g.as_array())
            .ok_or_else(|| Error::Config("SKOS JSON-LD missing @graph".into()))?;

        let mut schemes = BTreeMap::new();
        let mut concepts = BTreeMap::new();

        for node in graph {
            let id = node_id(node)?;
            let types = node_types(node);
            if types.iter().any(|t| t.ends_with("ConceptScheme") || t == "skos:ConceptScheme")
            {
                let label = lang_string(node, "prefLabel")
                    .or_else(|| lang_string(node, "skos:prefLabel"))
                    .or_else(|| lang_string(node, "label"))
                    .unwrap_or_else(|| id.clone());
                schemes.insert(
                    id.clone(),
                    ConceptScheme {
                        id: Curie::new(id),
                        label,
                        iri: node
                            .get("iri")
                            .and_then(|v| v.as_str())
                            .map(str::to_string),
                    },
                );
                continue;
            }
            if types
                .iter()
                .any(|t| t.ends_with("Concept") || t == "skos:Concept")
            {
                let concept = parse_concept(node, &id)?;
                concepts.insert(id, concept);
            }
        }

        // Infer narrower from broader when not asserted.
        let broader_edges: Vec<(String, String)> = concepts
            .values()
            .flat_map(|c| {
                c.broader
                    .iter()
                    .map(|b| (c.id.as_str().to_string(), b.as_str().to_string()))
                    .collect::<Vec<_>>()
            })
            .collect();
        for (child, parent) in broader_edges {
            if let Some(p) = concepts.get_mut(&parent) {
                let curie = Curie::new(child);
                if !p.narrower.iter().any(|n| n == &curie) {
                    p.narrower.push(curie);
                }
            }
        }

        if schemes.is_empty() {
            return Err(Error::Config("SKOS file has no ConceptScheme".into()));
        }
        Ok(Self { schemes, concepts })
    }

    pub fn get(&self, id: &str) -> Option<&SkosConcept> {
        self.concepts.get(id)
    }

    pub fn ontology_ref(&self, scheme_id: &str) -> OntologyRef {
        self.schemes
            .get(scheme_id)
            .map(|s| OntologyRef {
                id: s.id.clone(),
                label: s.label.clone(),
                iri: s.iri.clone(),
            })
            .unwrap_or_else(|| OntologyRef {
                id: Curie::new(scheme_id),
                label: scheme_id.to_string(),
                iri: None,
            })
    }

    pub fn depth(&self, id: &str) -> u32 {
        let mut depth = 0u32;
        let mut current = id.to_string();
        let mut seen = BTreeSet::new();
        while let Some(c) = self.concepts.get(&current) {
            if !seen.insert(current.clone()) {
                break;
            }
            match c.broader.first() {
                Some(b) => {
                    depth += 1;
                    current = b.as_str().to_string();
                }
                None => break,
            }
        }
        depth
    }

    pub fn path_labels(&self, id: &str) -> Vec<String> {
        let mut chain = Vec::new();
        let mut current = id.to_string();
        let mut seen = BTreeSet::new();
        while let Some(c) = self.concepts.get(&current) {
            if !seen.insert(current.clone()) {
                break;
            }
            chain.push(c.pref_label.clone());
            match c.broader.first() {
                Some(b) => current = b.as_str().to_string(),
                None => break,
            }
        }
        chain.reverse();
        chain
    }

    pub fn specificity(&self, id: &str) -> Specificity {
        let child_count = self
            .concepts
            .get(id)
            .map(|c| c.narrower.len())
            .unwrap_or(0);
        Specificity {
            depth: self.depth(id),
            child_count,
            is_leaf: child_count == 0,
        }
    }

    /// Build a [`Term`] shell from thesaurus data (matches filled by caller).
    pub fn term_shell(&self, id: &str, evidence: Vec<String>) -> Result<Term> {
        let c = self
            .get(id)
            .ok_or_else(|| Error::UnknownConcept(id.to_string()))?;
        let ontology = self.ontology_ref(c.in_scheme.as_str());
        let link = |targets: &[Curie], pred: &str| -> Vec<TermLink> {
            targets
                .iter()
                .map(|t| TermLink {
                    target_ontology: t.prefix().map(str::to_string),
                    target: t.clone(),
                    predicate: pred.to_string(),
                })
                .collect()
        };
        Ok(Term {
            id: c.id.clone(),
            ontology,
            pref_label: c.pref_label.clone(),
            alt_labels: c.alt_labels.clone(),
            definition: c.definition.clone(),
            specificity: self.specificity(id),
            broader: link(&c.broader, "broader"),
            narrower: link(&c.narrower, "narrower"),
            intra_matches: Vec::new(),
            inter_matches: link(&c.exact_match, "exactMatch")
                .into_iter()
                .chain(link(&c.close_match, "closeMatch"))
                .chain(link(&c.broad_match, "broadMatch"))
                .chain(link(&c.related_match, "relatedMatch"))
                .collect(),
            path_labels: self.path_labels(id),
            evidence,
        })
    }

    pub fn ancestors(&self, id: &str) -> Vec<Curie> {
        let mut out = Vec::new();
        let mut current = id.to_string();
        let mut seen = BTreeSet::new();
        while let Some(c) = self.concepts.get(&current) {
            if !seen.insert(current.clone()) {
                break;
            }
            for b in &c.broader {
                out.push(b.clone());
                current = b.as_str().to_string();
                break;
            }
            if c.broader.is_empty() {
                break;
            }
        }
        out
    }
}

fn node_id(node: &Value) -> Result<String> {
    node.get("id")
        .or_else(|| node.get("@id"))
        .and_then(|v| v.as_str())
        .map(str::to_string)
        .ok_or_else(|| Error::Config("SKOS node missing id".into()))
}

fn node_types(node: &Value) -> Vec<String> {
    match node.get("type").or_else(|| node.get("@type")) {
        Some(Value::String(s)) => vec![s.clone()],
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|v| v.as_str().map(str::to_string))
            .collect(),
        _ => Vec::new(),
    }
}

fn lang_string(node: &Value, key: &str) -> Option<String> {
    match node.get(key)? {
        Value::String(s) => Some(s.clone()),
        Value::Object(m) => m
            .get("@value")
            .and_then(|v| v.as_str())
            .map(str::to_string),
        Value::Array(a) => a.iter().find_map(|v| match v {
            Value::String(s) => Some(s.clone()),
            Value::Object(m) => m
                .get("@value")
                .and_then(|x| x.as_str())
                .map(str::to_string),
            _ => None,
        }),
        _ => None,
    }
}

fn id_list(node: &Value, key: &str) -> Vec<Curie> {
    match node.get(key) {
        Some(Value::String(s)) => vec![Curie::new(s.clone())],
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|v| match v {
                Value::String(s) => Some(Curie::new(s.clone())),
                Value::Object(m) => m
                    .get("@id")
                    .or_else(|| m.get("id"))
                    .and_then(|x| x.as_str())
                    .map(|s| Curie::new(s)),
                _ => None,
            })
            .collect(),
        Some(Value::Object(m)) => m
            .get("@id")
            .or_else(|| m.get("id"))
            .and_then(|x| x.as_str())
            .map(|s| vec![Curie::new(s)])
            .unwrap_or_default(),
        _ => Vec::new(),
    }
}

fn parse_concept(node: &Value, id: &str) -> Result<SkosConcept> {
    let pref = lang_string(node, "prefLabel")
        .or_else(|| lang_string(node, "skos:prefLabel"))
        .ok_or_else(|| Error::Config(format!("concept {id} missing prefLabel")))?;
    let alt = match node.get("altLabel").or_else(|| node.get("skos:altLabel")) {
        Some(Value::String(s)) => vec![s.clone()],
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|v| match v {
                Value::String(s) => Some(s.clone()),
                Value::Object(m) => m
                    .get("@value")
                    .and_then(|x| x.as_str())
                    .map(str::to_string),
                _ => None,
            })
            .collect(),
        Some(Value::Object(m)) => m
            .get("@value")
            .and_then(|x| x.as_str())
            .map(|s| vec![s.to_string()])
            .unwrap_or_default(),
        _ => Vec::new(),
    };
    let in_scheme = id_list(node, "inScheme")
        .into_iter()
        .next()
        .or_else(|| id_list(node, "skos:inScheme").into_iter().next())
        .unwrap_or_else(|| Curie::new("xrm:scheme"));
    Ok(SkosConcept {
        id: Curie::new(id),
        pref_label: pref,
        alt_labels: alt,
        definition: lang_string(node, "definition")
            .or_else(|| lang_string(node, "skos:definition")),
        in_scheme,
        broader: id_list(node, "broader")
            .into_iter()
            .chain(id_list(node, "skos:broader"))
            .collect(),
        narrower: id_list(node, "narrower")
            .into_iter()
            .chain(id_list(node, "skos:narrower"))
            .collect(),
        exact_match: id_list(node, "exactMatch")
            .into_iter()
            .chain(id_list(node, "skos:exactMatch"))
            .collect(),
        close_match: id_list(node, "closeMatch")
            .into_iter()
            .chain(id_list(node, "skos:closeMatch"))
            .collect(),
        broad_match: id_list(node, "broadMatch")
            .into_iter()
            .chain(id_list(node, "skos:broadMatch"))
            .collect(),
        related_match: id_list(node, "relatedMatch")
            .into_iter()
            .chain(id_list(node, "skos:relatedMatch"))
            .collect(),
        top_concept: node
            .get("topConcept")
            .and_then(|v| v.as_bool())
            .unwrap_or(false),
    })
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct UnusedGuard;
