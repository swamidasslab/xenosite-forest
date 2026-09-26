//! Compact URI (CURIE) identifiers for concepts and external mappings.

use serde::{Deserialize, Serialize};
use std::fmt;

/// `prefix:local` identifier (e.g. `xrm:0000100`, `mesh:D050216`).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Curie(pub String);

impl Curie {
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Ontology / prefix part before `:`, if present.
    pub fn prefix(&self) -> Option<&str> {
        self.0.split_once(':').map(|(p, _)| p)
    }

    pub fn local(&self) -> &str {
        self.0
            .split_once(':')
            .map(|(_, l)| l)
            .unwrap_or(self.0.as_str())
    }
}

impl fmt::Display for Curie {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for Curie {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}

impl From<String> for Curie {
    fn from(s: String) -> Self {
        Self::new(s)
    }
}
