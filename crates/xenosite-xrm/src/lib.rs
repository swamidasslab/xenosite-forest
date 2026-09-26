//! Med-chem-oriented SKOS thesaurus for xenobiotic metabolism reaction naming.
//!
//! The missing layer between reaction-generating SMARTS/rules and human-useful
//! med-chem metabolism names — **not** an enzyme, pathway, compound, or
//! named-synthesis ontology. Separate from RXNO/MOP/GO/… with SSSOM mappings
//! (see `data/ontology/RELATED_ONTOLOGIES.md`).
//!
//! # Boundary
//!
//! This crate **must not** depend on `xenosite-forest` (or any forest rule
//! implementation). Naming is driven only by text / JSON-LD / JSONL / SSSOM
//! config. Opaque CURIEs such as `forest.rule:Hydroxylation` may appear in
//! config and in caller-supplied tags; they are strings, never resolved by
//! importing forest code. Forest SMIRKS may be harvested offline for validation.
//!
//! # Formats
//!
//! - **SKOS** (JSON-LD) — concept thesaurus
//! - **SSSOM** (TSV) — inter- and intra-ontology mappings
//! - **JSONL** — structural assignment rules over a mapped reactant→product
//! - **manifest.json** — paths to the above

mod assignment;
mod bundle;
mod chemistry;
mod curie;
mod error;
mod namer;
mod skos;
mod sssom;
mod term;

pub use assignment::{AssignmentHit, AssignmentRule, Assignments};
pub use bundle::{spines, AnnotationBundle};
pub use chemistry::{AtomMap, LocalizedTag, MappedReaction}; // LocalizedTag for site-localized caller tags
pub use curie::Curie;
pub use error::{Error, Result};
pub use namer::{Namer, NamerConfig};
pub use skos::{ConceptScheme, SkosConcept, Thesaurus};
pub use sssom::{MappingPredicate, SssomMapping, SssomTable};
pub use term::{OntologyRef, SiteRef, Specificity, Term, TermLink};

/// Default bundled manifest (relative to this crate's `data/` directory).
pub const DEFAULT_MANIFEST: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/manifest.json");
