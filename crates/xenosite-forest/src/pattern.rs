//! Pattern records. The algorithm reads these; it does not subclass them.
//!
//! `SiteKind`, `Edit`, and `Effect` are the categories. Methide is an effect
//! field, not a pathway flag.

/// What kind of site this pattern names. Discovery indexes follow this.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SiteKind {
    Atom,
}

/// How the pattern edits the matched atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    /// Graph-add a hydroxyl oxygen at map 1 (chematic SMIRKS dialect aside).
    Hydroxyl,
    /// Apply this SMIRKS at the unique-edit match.
    Smirks(String),
}

/// One concrete outcome. Filters read these fields.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Effect {
    pub adds: Option<String>,
    pub removes: Option<String>,
    pub cleaves: bool,
    /// Effect bit, not a `pathways=("methide",)` switch.
    pub methide: bool,
}

/// What a SMARTS pattern can do, before a match.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatternInfo {
    pub name: String,
    pub smarts: String,
    pub site_kind: SiteKind,
    pub edit: Edit,
    pub effect: Effect,
}

impl PatternInfo {
    pub fn new(
        name: impl Into<String>,
        smarts: impl Into<String>,
        edit: Edit,
        effect: Effect,
    ) -> Self {
        Self {
            name: name.into(),
            smarts: smarts.into(),
            site_kind: SiteKind::Atom,
            edit,
            effect,
        }
    }

    pub fn hydroxyl(name: impl Into<String>, smarts: impl Into<String>) -> Self {
        Self::new(
            name,
            smarts,
            Edit::Hydroxyl,
            Effect {
                adds: Some("O".into()),
                removes: Some("H".into()),
                cleaves: false,
                methide: false,
            },
        )
    }
}

/// Bag `filter_sites` sees after unique-edit, before the edit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteInfo {
    pub site: usize,
    pub pattern: PatternInfo,
}

/// One metabolize emission: discovery site, pattern, rule namespace, product CSMIs.
///
/// `rule_path` is leaf-first (emitting rule, then each containing [`crate::ruleset::RuleSet`]),
/// matching Python `info["rule"]` / addition chain order. Unnamed sets stay on the
/// chain as `None`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Emission {
    pub site: usize,
    pub pattern_name: String,
    pub rule_path: Vec<Option<String>>,
    pub products: Vec<String>,
}

impl Emission {
    /// Named segments of [`Self::rule_path`] (unnamed sets omitted).
    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    /// Emitting (leaf) rule name, if the leaf was named.
    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }
}
