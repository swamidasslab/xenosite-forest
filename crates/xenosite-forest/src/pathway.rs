//! Shared options for applying / walking a pathway on a molecule.
//!
//! A [`crate::StepSequence`] is an ordered container of steps that acts like
//! a step. Public **`apply`** (and [`crate::random_path`]) take these options
//! and return one outcome shape for length 1..n.
//!
//! Filtering and label ensure-missing (never overwrite) belong on that shared
//! pathway layer, not on low-level materialize helpers.

/// Filters for choosing among products when applying steps or sampling a walk.
///
/// Defaults are **off**: every materialized product is eligible (including
/// multi-component SMILES and revisiting a CSMI already on the path).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PathwayOptions {
    /// Refuse products whose CSMI contains `.` (cleavage bags / fragments).
    pub skip_multicomponent: bool,
    /// Refuse products whose CSMI already appeared on the walk (no loops).
    pub skip_seen: bool,
}

impl PathwayOptions {
    /// Common test / fuzz preset: single-component, no revisits.
    pub const fn no_loops_or_fragments() -> Self {
        Self {
            skip_multicomponent: true,
            skip_seen: true,
        }
    }

    /// True when `smi` is allowed under these options given `seen` so far.
    pub fn allows(&self, smi: &str, seen: &std::collections::HashSet<String>) -> bool {
        if smi.is_empty() {
            return false;
        }
        if self.skip_multicomponent && smi.contains('.') {
            return false;
        }
        if self.skip_seen && seen.contains(smi) {
            return false;
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn default_allows_multicomponent_and_seen() {
        let opts = PathwayOptions::default();
        let mut seen = HashSet::new();
        seen.insert("CC".into());
        assert!(opts.allows("CCO.O", &seen));
        assert!(opts.allows("CC", &seen));
    }

    #[test]
    fn no_loops_preset_refuses_both() {
        let opts = PathwayOptions::no_loops_or_fragments();
        let mut seen = HashSet::new();
        seen.insert("CC".into());
        assert!(!opts.allows("CCO.O", &seen));
        assert!(!opts.allows("CC", &seen));
        assert!(opts.allows("CCO", &seen));
    }
}
