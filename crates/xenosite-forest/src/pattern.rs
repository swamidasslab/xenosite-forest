//! Pattern records. The algorithm reads these; it does not subclass them.
//!
//! `SiteKind`, `Edit`, and `Effect` are the categories. Methide is an effect
//! field, not a pathway flag.

use std::collections::BTreeMap;

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

/// Constraint that picks one branch of a SMARTS OR once atoms are known.
///
/// Same role as Python ``When``. Used when ``delta_formula`` (or other effect
/// fields) disagree across OR arms — annotate each arm, do not invent a search
/// branch.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct When {
    pub map: u16,
    pub z: Option<u8>,
    pub h: Option<u8>,
}

impl When {
    pub fn atomic(map: u16, z: u8) -> Self {
        Self {
            map,
            z: Some(z),
            h: None,
        }
    }

    pub fn atomic_h(map: u16, z: u8, h: u8) -> Self {
        Self {
            map,
            z: Some(z),
            h: Some(h),
        }
    }
}

/// Known element symbols, longest first (for bag strings like ``Cl``, ``Br``).
const ELEMENT_SYMBOLS: &[&str] = &["At", "Br", "Cl", "I", "F", "O", "N", "S", "P", "C", "H"];

/// Parse an ``adds`` / ``removes`` bag (``"OO"``, ``"HH"``, ``"Cl"``, ``"OH"``)
/// into element → count. Unknown characters are skipped.
pub fn bag_counts(bag: &str) -> BTreeMap<String, i32> {
    let mut counts = BTreeMap::new();
    let bytes = bag.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let rest = &bag[i..];
        let mut matched = None;
        for sym in ELEMENT_SYMBOLS {
            if rest.starts_with(sym) {
                matched = Some(*sym);
                break;
            }
        }
        match matched {
            Some(sym) => {
                *counts.entry(sym.to_string()).or_insert(0) += 1;
                i += sym.len();
            }
            None => i += 1,
        }
    }
    counts
}

/// Net formula change from ``adds`` / ``removes`` bags. Zero-count keys omitted.
pub fn bag_delta_formula(adds: Option<&str>, removes: Option<&str>) -> BTreeMap<String, i32> {
    let mut delta = BTreeMap::new();
    if let Some(bag) = adds {
        for (el, n) in bag_counts(bag) {
            *delta.entry(el).or_insert(0) += n;
        }
    }
    if let Some(bag) = removes {
        for (el, n) in bag_counts(bag) {
            *delta.entry(el).or_insert(0) -= n;
        }
    }
    delta.retain(|_, n| *n != 0);
    delta
}

/// Junction bags ± named leave (leave counts subtracted). Zeros omitted.
pub fn compose_delta_formula(
    adds: Option<&str>,
    removes: Option<&str>,
    leave_formula: &BTreeMap<String, i32>,
) -> BTreeMap<String, i32> {
    let mut delta = bag_delta_formula(adds, removes);
    for (el, n) in leave_formula {
        *delta.entry(el.clone()).or_insert(0) -= n;
    }
    delta.retain(|_, n| *n != 0);
    delta
}

/// Named methyl leave (``cleave_side_group`` / ``leave_count=1`` methyl carbon).
pub fn leave_me() -> BTreeMap<String, i32> {
    BTreeMap::from([("C".into(), 1), ("H".into(), 3)])
}

/// Named methylene leave (benzodioxole CH2).
pub fn leave_ch2() -> BTreeMap<String, i32> {
    BTreeMap::from([("C".into(), 1), ("H".into(), 2)])
}

/// Named oxygen leave (nitroaromatic N–O cleavage / dehydration O leave).
pub fn leave_o() -> BTreeMap<String, i32> {
    BTreeMap::from([("O".into(), 1)])
}

/// Two-oxygen leave (nitro → amine style double O cleavage).
pub fn leave_oo() -> BTreeMap<String, i32> {
    BTreeMap::from([("O".into(), 2)])
}

/// Resolve a ``cleave_side_group`` leave label into a formula bag, if known.
pub fn named_leave_formula(leave: &str) -> Option<BTreeMap<String, i32>> {
    match leave {
        "Me" => Some(leave_me()),
        _ => None,
    }
}

/// Merge two delta maps (pair ends). Zero keys dropped.
pub fn merge_delta_formula(
    left: &BTreeMap<String, i32>,
    right: &BTreeMap<String, i32>,
) -> BTreeMap<String, i32> {
    let mut out = left.clone();
    for (el, n) in right {
        *out.entry(el.clone()).or_insert(0) += n;
    }
    out.retain(|_, n| *n != 0);
    out
}

/// One concrete outcome. Filters read these fields.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Effect {
    pub adds: Option<String>,
    pub removes: Option<String>,
    /// Declared net formula change (element → delta). Zeros omitted.
    ///
    /// Sealed from junction ``adds`` / ``removes`` bags minus
    /// [`Self::leave_formula`]. Cleavage: leave as negative, O/H at the cut
    /// from the bags. When OR arms disagree (halogen removal), each arm
    /// carries its own map under a [`When`].
    pub delta_formula: BTreeMap<String, i32>,
    /// Named leaving-piece formula (positive counts). Cleavage only.
    ///
    /// Sealed into [`Self::delta_formula`] as a negative contribution. Empty
    /// when the leave is open (`leave_count` None) or already encoded in
    /// ``removes`` (e.g. dehydration ``OH``).
    pub leave_formula: BTreeMap<String, i32>,
    pub cleaves: bool,
    /// Named leaving heavy-atom count (methyl dealkylation = 1). `None` = open.
    pub leave_count: Option<u16>,
    /// Effect bit, not a `pathways=("methide",)` switch.
    pub methide: bool,
    /// Capability: pair/path may dearomatize. Resolved against system aromaticity.
    pub dearomatizes: bool,
    /// Methide / alkyl partner element hint (`"C"`). Filters read this.
    pub partner: Option<String>,
    /// Branch constraint when this effect is one arm of a SMARTS OR.
    pub when: Option<When>,
}

impl Effect {
    /// Fill [`Self::delta_formula`] from bags − leave when still empty.
    pub fn sealed(mut self) -> Self {
        if self.delta_formula.is_empty() {
            self.delta_formula = compose_delta_formula(
                self.adds.as_deref(),
                self.removes.as_deref(),
                &self.leave_formula,
            );
        }
        self
    }

    /// Recompute [`Self::delta_formula`] from bags − leave (after mutating leave).
    pub fn reseal_delta(mut self) -> Self {
        self.delta_formula = compose_delta_formula(
            self.adds.as_deref(),
            self.removes.as_deref(),
            &self.leave_formula,
        );
        self
    }

    /// Declared delta, deriving from bags − leave if the map was never sealed.
    pub fn resolved_delta_formula(&self) -> BTreeMap<String, i32> {
        if self.delta_formula.is_empty()
            && (self.adds.as_ref().is_some_and(|s| !s.is_empty())
                || self.removes.as_ref().is_some_and(|s| !s.is_empty())
                || !self.leave_formula.is_empty())
        {
            compose_delta_formula(
                self.adds.as_deref(),
                self.removes.as_deref(),
                &self.leave_formula,
            )
        } else {
            self.delta_formula.clone()
        }
    }

    /// Attach a named leave formula and reseal delta.
    pub fn with_leave_formula(mut self, leave: BTreeMap<String, i32>) -> Self {
        self.leave_formula = leave;
        self.reseal_delta()
    }
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
    /// When OR arms need distinct effect data (especially ``delta_formula``).
    /// Empty ⇒ use [`Self::effect`] alone. Non-empty ⇒ resolve via [`When`].
    pub possibilities: Vec<Effect>,
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
    /// Heap preference for `find_path` (higher pops sooner on the max-heap).
    /// Default `0`. Negative demotes patterns that are real but less likely /
    /// counter-directional (e.g. Hydrogenation). Soft only — never drops.
    /// Schema trial; see HEURISTICS.
    pub search_bias: i8,
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
            effect: effect.sealed(),
            possibilities: Vec::new(),
            skip_same_rings: false,
            cleave_side_group: None,
            search_bias: 0,
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
                ..Effect::default()
            },
        )
    }

    /// Attach When-branched possibilities (each sealed). Span [`Self::effect`]
    /// stays the collapsed / representative effect for rule-level filters.
    pub fn with_possibilities(mut self, branches: impl IntoIterator<Item = Effect>) -> Self {
        self.possibilities = branches.into_iter().map(Effect::sealed).collect();
        self
    }

    /// Set cleavage side groups (leave, keep). Equal labels ⇒ swappable.
    ///
    /// Known leave labels (e.g. ``Me``) fill [`Effect::leave_formula`] and
    /// reseal [`Effect::delta_formula`] when leave was empty.
    pub fn with_cleave_side_group(
        mut self,
        leave: impl Into<String>,
        keep: impl Into<String>,
    ) -> Self {
        let leave = leave.into();
        let keep = keep.into();
        if self.effect.leave_formula.is_empty() {
            if let Some(formula) = named_leave_formula(&leave) {
                self.effect.leave_formula = formula;
                self.effect = self.effect.reseal_delta();
            }
        }
        for arm in &mut self.possibilities {
            if arm.leave_formula.is_empty() {
                if let Some(formula) = named_leave_formula(&leave) {
                    arm.leave_formula = formula.clone();
                    arm.delta_formula = compose_delta_formula(
                        arm.adds.as_deref(),
                        arm.removes.as_deref(),
                        &arm.leave_formula,
                    );
                }
            }
        }
        self.cleave_side_group = Some((leave, keep));
        self
    }

    /// Soft heap preference (higher first). Does not filter or drop.
    pub fn with_search_bias(mut self, bias: i8) -> Self {
        self.search_bias = bias;
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

    /// Effects to consider for this pattern (possibilities, or the single effect).
    pub fn effect_arms(&self) -> Vec<&Effect> {
        if self.possibilities.is_empty() {
            vec![&self.effect]
        } else {
            self.possibilities.iter().collect()
        }
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
    /// Site-local shell deltas for filters — **same shape** as
    /// [`crate::matched_atom::AlignedShells`] scoped to this site
    /// ([`crate::matched_atom::AlignedShells::at_sites`] or a forecast).
    /// Site selection reads these shells; product closeness still uses the
    /// full-molecule align cost. `None` until search fills it.
    pub shell_forecast: Option<crate::matched_atom::AlignedShells>,
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
    /// From [`PatternInfo::search_bias`] (pair: min of both ends).
    pub search_bias: i8,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bag_delta_hydroxyl() {
        let d = bag_delta_formula(Some("O"), Some("H"));
        assert_eq!(d.get("O"), Some(&1));
        assert_eq!(d.get("H"), Some(&-1));
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn bag_delta_omits_zeros() {
        let d = bag_delta_formula(Some("O"), Some("O"));
        assert!(d.is_empty());
    }

    #[test]
    fn bag_counts_two_letter_halogen() {
        assert_eq!(bag_counts("Cl").get("Cl"), Some(&1));
        assert_eq!(bag_counts("Br").get("Br"), Some(&1));
        assert_eq!(bag_counts("OH").get("O"), Some(&1));
        assert_eq!(bag_counts("OH").get("H"), Some(&1));
    }

    #[test]
    fn seal_fills_delta_formula() {
        let e = Effect {
            adds: Some("OO".into()),
            removes: Some("HH".into()),
            ..Effect::default()
        }
        .sealed();
        assert_eq!(e.delta_formula.get("O"), Some(&2));
        assert_eq!(e.delta_formula.get("H"), Some(&-2));
    }

    #[test]
    fn hydroxyl_pattern_carries_delta() {
        let p = PatternInfo::hydroxyl("h", "[#6h1:1]");
        assert_eq!(p.effect.delta_formula.get("O"), Some(&1));
        assert_eq!(p.effect.delta_formula.get("H"), Some(&-1));
    }

    #[test]
    fn cleavage_me_leave_is_negative_plus_junction_o() {
        let d = compose_delta_formula(Some("OO"), None, &leave_me());
        assert_eq!(d.get("C"), Some(&-1));
        assert_eq!(d.get("H"), Some(&-3));
        assert_eq!(d.get("O"), Some(&2));
    }

    #[test]
    fn cleave_side_group_me_fills_leave_formula() {
        let p = PatternInfo::new(
            "methyl_alcohol",
            "[#6H3:1][#8H0:2]",
            Edit::Smirks("[C:1][O:2]>>[O:2].[C:1]O".into()),
            Effect {
                adds: Some("O".into()),
                cleaves: true,
                leave_count: Some(1),
                ..Effect::default()
            },
        )
        .with_cleave_side_group("Me", "hetero");
        assert_eq!(p.effect.leave_formula, leave_me());
        assert_eq!(p.effect.delta_formula.get("C"), Some(&-1));
        assert_eq!(p.effect.delta_formula.get("H"), Some(&-3));
        assert_eq!(p.effect.delta_formula.get("O"), Some(&1));
    }
}
