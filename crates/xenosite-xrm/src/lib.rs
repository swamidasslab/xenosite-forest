//! Config-driven xenobiotic reaction naming.
//!
//! # Boundary
//!
//! This crate **must not** depend on `xenosite-forest` (or any forest rule
//! implementation). Naming is driven only by text / JSON-LD / JSONL / SSSOM
//! config. Opaque CURIEs such as `forest.rule:Hydroxylation` may appear in
//! config and in caller-supplied tags; they are strings, never resolved by
//! importing forest code.
//!
//! # Formats
//!
//! - **SKOS** (JSON-LD) — concept thesaurus
//! - **SSSOM** (TSV) — inter- and intra-ontology mappings
//! - **JSONL** — structural assignment rules over a mapped reactant→product
//! - **manifest.json** — paths to the above

mod assignment;
mod chemistry;
mod curie;
mod error;
mod namer;
mod skos;
mod sssom;
mod term;

pub use assignment::{AssignmentRule, Assignments};
pub use chemistry::{AtomMap, MappedReaction};
pub use curie::Curie;
pub use error::{Error, Result};
pub use namer::{Namer, NamerConfig};
pub use skos::{ConceptScheme, SkosConcept, Thesaurus};
pub use sssom::{MappingPredicate, SssomMapping, SssomTable};
pub use term::{OntologyRef, Specificity, Term, TermLink};

/// Default bundled manifest (relative to this crate's `data/` directory).
pub const DEFAULT_MANIFEST: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/manifest.json");
