//! Pattern records. The algorithm reads these; it does not subclass them.
//!
//! `SiteKind`, `Edit`, and `Effect` are the categories. Methide is an effect
//! field, not a pathway flag.

/// What kind of site this pattern names. Discovery indexes follow this.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SiteKind {
    Atom,
    Bond,
    DirectedBond,
    AtomPair,
}

/// How the pattern edits the matched atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    /// Graph-add a hydroxyl oxygen at map 1 (chematic SMIRKS dialect aside).
    Hydroxyl,
    /// Apply this SMIRKS at the unique-edit match.
    Smirks(String),
    /// Resonance-pair endpoint (path flip via [`crate::pair_edit`]).
    PairEndpoint(String),
}

/// One concrete outcome. Filters read these fields.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Effect {
    pub adds: Option<String>,
    pub removes: Option<String>,
    pub cleaves: bool,
    /// Named leaving heavy-atom count (methyl dealkylation = 1). `None` = open.
    pub leave_count: Option<u16>,
    /// Effect bit, not a `pathways=("methide",)` switch.
    pub methide: bool,
    /// Capability: pair/path may dearomatize. Resolved against system aromaticity.
    pub dearomatizes: bool,
    /// Methide / alkyl partner element hint (`"C"`). Filters read this.
    pub partner: Option<String>,
}

/// What a SMARTS pattern can do, before a match.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatternInfo {
    pub name: String,
    pub smarts: String,
    pub site_kind: SiteKind,
    /// Atom-map numbers that form the discovery site (Python `site_map`).
    pub site_map: Vec<u16>,
    pub edit: Edit,
    pub effect: Effect,
    /// Refuse single-to-double when maps 1 and 2 share the same ring set.
    pub skip_same_rings: bool,
    /// Cleavage side groups `(leave, keep)` aligned to [`Self::site_map`] order.
    ///
    /// Cleavage analogue of pair ``swap_group``, but pooled **across rules** at
    /// expand (not within one ResonancePair). `None` = ungrouped: fold key is
    /// fragment CSMI multiset only. Equal non-empty strings = swappable sides
    /// (normalize order in the key). Shared leave labels (e.g. `"Me"` on O- and
    /// N-methyl dealks) let distinct rules fold when fragments match.
    pub cleave_side_group: Option<(String, String)>,
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
            site_map: vec![1],
            edit,
            effect,
            skip_same_rings: false,
            cleave_side_group: None,
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
                leave_count: None,
                methide: false,
                dearomatizes: false,
                partner: None,
            },
        )
    }

    /// Set cleavage side groups (leave, keep). Equal labels ⇒ swappable.
    pub fn with_cleave_side_group(
        mut self,
        leave: impl Into<String>,
        keep: impl Into<String>,
    ) -> Self {
        self.cleave_side_group = Some((leave.into(), keep.into()));
        self
    }

    /// First map in [`Self::site_map`], or 1.
    pub fn primary_map(&self) -> u16 {
        self.site_map.first().copied().unwrap_or(1)
    }

    /// Resolved side-group signature for cross-rule cleavage fold.
    pub fn cleave_side_sig(&self) -> CleaveSideSig {
        CleaveSideSig::resolve(self.cleave_side_group.as_ref())
    }
}

/// How cleavage sides participate in cross-rule Or fold keys.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CleaveSideSig {
    /// No side-group data — fold by fragment multiset only.
    Ungrouped,
    /// Directed leave vs keep (unequal labels).
    Directed(String, String),
    /// Swappable sides (equal non-empty labels).
    Swap(String),
}

impl CleaveSideSig {
    pub fn resolve(groups: Option<&(String, String)>) -> Self {
        match groups {
            None => Self::Ungrouped,
            Some((a, b)) if a == b && !a.is_empty() => Self::Swap(a.clone()),
            Some((a, b)) => Self::Directed(a.clone(), b.clone()),
        }
    }

    /// Fold key with sorted fragment CSMIs (both sides first-class).
    pub fn fold_key(&self, fragments: &[String]) -> CleaveFoldKey {
        let mut fragments = fragments.to_vec();
        fragments.sort();
        CleaveFoldKey {
            side: self.clone(),
            fragments,
        }
    }
}

/// Cross-rule cleavage Or bucket: side signature + fragment CSMI multiset.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CleaveFoldKey {
    pub side: CleaveSideSig,
    pub fragments: Vec<String>,
}

/// Bag `filter_sites` sees after unique-edit, before the edit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteInfo {
    pub site: usize,
    /// Primary-map atoms that share this site's unique-edit class (sorted).
    /// Always includes `site`. Length > 1 when topology collapses equivalents.
    pub orbit: Vec<usize>,
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
    /// Primary-map orbit passed down from unique-edit (see [`SiteInfo::orbit`]).
    pub site_orbit: Vec<usize>,
    /// Discovery site atoms for cleavage bookkeeping ([`crate::canonical_plan::CleavageSide`]).
    pub site_atoms: Vec<usize>,
    pub cleaves: bool,
    pub pattern_name: String,
    pub rule_path: Vec<Option<String>>,
    pub products: Vec<String>,
    /// Elementary steps for this hop (identity or quinone-shaped expansion).
    /// Bind with [`crate::canonical_plan::Deps::bind`] for precedes / replay.
    pub plan: Vec<crate::canonical_plan::Step>,
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
