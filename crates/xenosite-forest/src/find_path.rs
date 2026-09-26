//! First-run plan-guided search over a [`crate::ruleset::RuleSet`].
//!
//! Walks carry tagged [`crate::forest_mol::ForestMol`] (structure/`csmi` cache,
//! atom tags through edits). Expands via [`RuleSet::candidates`]: site–pattern–
//! parent triples. A search reads [`PatternInfo`] / effect fields to decide,
//! then materializes only for survivors. Nested sets stay namespaces on each
//! step's leaf-first `rule_path`.
//!
//! Outcomes carry [`crate::canonical_plan::Deps`] plans (elementary steps +
//! precedes + [`crate::canonical_plan::Maybe`] cleavage bags). Composite leaves
//! own a `canonical_plan` hook (Python) that returns steps named after existing
//! rules. Closer uses atom-diff cost. Child diffs: tag-lift + extend; keep only
//! cost 0 (skip MCS). Otherwise fresh MCS (`mcs_lift_rematch`). Tag-lift
//! impossible → `mcs_lift_fallback` (expect zero).

use std::cmp::Ordering;
use std::collections::{BTreeSet, BinaryHeap, HashMap, HashSet};
use std::rc::Rc;

use crate::ForestError;
use crate::candidate::Candidate;
use crate::canonical_plan::{CanonicalStep, CleavageSide, Deps, Maybe, as_deps};
use crate::forest_mol::ForestMol;
use crate::labels::Tag;
use crate::mol::{canon_of, parse_mol};
use crate::pattern::{CleaveFoldKey, CleaveSideSig, PatternInfo, SiteInfo};
use crate::rules::default_ruleset;
use crate::ruleset::RuleSet;

/// Billed work for one search (Python `PathCounters` subset).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PathCounters {
    pub nodes: usize,
    pub mol_edits: usize,
    pub expansions: usize,
    /// Target hit not emitted: same linearizations **or** same rule/Maybe
    /// skeleton as an already-yielded plan. Total wasted duplicate finds
    /// ([`Self::dropped_exact_plan`] + [`Self::dropped_skeleton_twin`]).
    pub dropped_duplicate_plan: usize,
    /// Subset of [`Self::dropped_duplicate_plan`]: exact
    /// [`Deps::same_linearizations`] match (same step identities).
    pub dropped_exact_plan: usize,
    /// Subset of [`Self::dropped_duplicate_plan`]: remapped free-step twin via
    /// [`Deps::same_rule_maybe_skeleton`] but **not** exact linearizations.
    /// How often search reaches the target again only to learn indices moved.
    pub dropped_skeleton_twin: usize,
    /// Hit that *would* match [`Deps::dominates_extension_of`] against a yield.
    /// Counted only — not used to drop or abort (that lever over-collapsed
    /// multipath; HEURISTICS).
    pub signal_contained_plan: usize,
    /// Expand skipped: O-add then O-remove (or reverse) with **equal** site
    /// sets (orbit-aware for singletons). Allowed when sites differ.
    pub blocked_circular_oxygen: usize,
    /// Product lacked a Chematic `canonical_smiles_stable_key` — CSMI dedup
    /// skipped for that child (fail-closed; may re-explore).
    pub unstable_csmi_key: usize,
    /// Declared [`crate::pattern::Effect::delta_formula`] ≠ observed product
    /// Δformula (heavy). Soft — warn only; never drops.
    pub formula_delta_mismatch: usize,
    /// Structured mismatch details (Python ``formula_delta_mismatches``).
    pub formula_delta_mismatches: Vec<crate::formula_check::FormulaDeltaMismatch>,
    /// Child diff fell back to a fresh MCS because tag-lift was impossible
    /// (no shared tags). Expect zero.
    pub mcs_lift_fallback: usize,
    /// Extended lift was not cost 0 → used fresh MCS instead (no Aut chase).
    /// Counted for derisk; not a Drop assert yet.
    pub mcs_lift_rematch: usize,
    /// Match diversity: heap item re-pushed because `−log(n+1)` went stale.
    pub diversity_repush: usize,
    /// When true, [`Drop`] does not assert zero mismatches (intentional tests).
    #[cfg(test)]
    pub allow_formula_delta_mismatch: bool,
    /// When true, [`Drop`] does not assert zero MCS-lift fallbacks.
    #[cfg(test)]
    pub allow_mcs_lift_fallback: bool,
}

impl PathCounters {
    pub fn billed(&self) -> usize {
        self.mol_edits + self.nodes
    }

    /// Optimization signal: duplicate yield drops + contained-extension matches.
    /// Do not use to abort walks; more careful heuristics might reduce later.
    pub fn plan_drops(&self) -> usize {
        self.dropped_duplicate_plan + self.signal_contained_plan
    }

    /// Record a soft formula-delta mismatch (count + detail list).
    pub fn record_formula_delta_mismatch(
        &mut self,
        detail: crate::formula_check::FormulaDeltaMismatch,
    ) {
        self.formula_delta_mismatch += 1;
        self.formula_delta_mismatches.push(detail);
    }

    /// Fundamental check: no soft formula-delta mismatches on this search.
    pub fn assert_formula_delta_clean(&self) {
        assert_eq!(
            self.formula_delta_mismatch, 0,
            "formula_delta_mismatch must be zero; recorded {:?}",
            self.formula_delta_mismatches
        );
    }

    /// Fundamental check: child diffs must lift via tags + generators.
    pub fn assert_mcs_lift_clean(&self) {
        assert_eq!(
            self.mcs_lift_fallback, 0,
            "mcs_lift_fallback must be zero (tag-lift + extend or MCS)"
        );
    }
}

#[cfg(test)]
impl Drop for PathCounters {
    fn drop(&mut self) {
        if std::thread::panicking() {
            return;
        }
        if !self.allow_formula_delta_mismatch {
            self.assert_formula_delta_clean();
        }
        if !self.allow_mcs_lift_fallback {
            self.assert_mcs_lift_clean();
        }
    }
}

/// One accepted edit on a walk: RuleSet namespace path + pattern + kept product.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathStep {
    pub rule_path: Vec<Option<String>>,
    pub pattern_name: String,
    pub site: usize,
    /// Unique-edit primary-map orbit (includes `site`). Plan equivalence under
    /// automorphism compares sites via [`crate::same_site_orbit`].
    pub site_orbit: Vec<usize>,
    /// Kept fragment CSMI (the search node).
    pub product: String,
    /// Cleaved-off fragment CSMIs (not expanded).
    pub sides: Vec<String>,
}

impl PathStep {
    /// Named segments of [`Self::rule_path`] (unnamed sets omitted).
    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }

    /// Same rule/pattern and sites in one unique-edit orbit.
    pub fn same_site_class(&self, other: &Self) -> bool {
        self.pattern_name == other.pattern_name
            && self.leaf_rule() == other.leaf_rule()
            && crate::same_site_orbit(self.site, &self.site_orbit, other.site, &other.site_orbit)
    }

    fn from_emission(emission: &ForestEmission, product: String, sides: Vec<String>) -> Self {
        Self {
            rule_path: emission.rule_path.clone(),
            pattern_name: emission.pattern_name.clone(),
            site: emission.site,
            site_orbit: emission.site_orbit.clone(),
            product,
            sides,
        }
    }
}

/// One reactant→target hit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathOutcome {
    pub steps: Vec<PathStep>,
    /// Elementary plan with precedes and cleavage [`Maybe`] (on the plan).
    pub plan: Deps,
    pub smiles: String,
}

impl PathOutcome {
    /// Cleavage-side bags on the plan (Python `PathOutcome.maybe`).
    pub fn maybe(&self) -> &Maybe {
        self.plan.maybe()
    }

    /// [`Deps::allows`] / [`Maybe::allows`] — discarded side or overlapping site.
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        self.plan.allows(site, side)
    }
}

#[derive(Clone)]
struct Walk {
    mol: ForestMol,
    steps: Vec<PathStep>,
    plan: Vec<CanonicalStep>,
    /// Accumulated cleavage bags (attached to [`Deps`] at hit yield).
    maybe: Vec<CleavageSide>,
    /// Uncleared ring-open sites (Python `_Walk.opens`).
    opens: Vec<BTreeSet<usize>>,
    /// O-adding hops on this walk (hydroxylation / hydrate / …).
    o_added: Vec<OxygenSite>,
    /// O-removing hops on this walk (dehydration).
    o_removed: Vec<OxygenSite>,
    /// Dedup keys from root through [`Self::mol`] (inclusive). A child whose
    /// stable key is already here is a cycle — refuse enqueue.
    ancestors: HashSet<String>,
    /// Parent's [`crate::atom_diff::AtomDiff::cost`] when this walk was
    /// enqueued. `None` = root (always expand).
    parent_cost: Option<usize>,
    /// Diff of **this** mol vs target. Filled at enqueue by
    /// [`try_atom_diff_for_child`] / cleavage lift when safe; else `None`
    /// and pop runs full MCS.
    diff: Option<crate::atom_diff::AtomDiff>,
}

/// Site of an oxygen add/remove hop. Multi-atom `atoms` (from `site_map`)
/// compare by set equality; singletons use unique-edit / top orbits.
#[derive(Clone, Debug)]
struct OxygenSite {
    site: usize,
    orbit: Vec<usize>,
    atoms: BTreeSet<usize>,
}

/// How closeness / improvement combine into one soft score.
///
/// Live [`Close`] / [`Improve`] / [`Add`] are **log-space** of the same
/// positive factors [`Product`] multiplies: per metric
/// `imp = max(0, parent−child)+1` and `close = SCALE/(1+child)`.
/// `Add = log(∏ factors)` ≡ ranking of `Product` (fixed-point ln).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MatchCombine {
    /// Σ `ln(SCALE/(1+child))` over active metrics.
    Close,
    /// Σ `ln(max(0, parent−child)+1)` over active metrics.
    Improve,
    /// Close + Improve = `ln` of the Product factors (ablation / former default).
    Add,
    /// Direct multiply of the same factors (ablation / `--score product-*`).
    Product,
    /// Linear costs → score: Σ `−child` (close only).
    LinNegC,
    /// Linear: Σ `−parent − child`.
    LinNegPC,
    /// Linear: Σ `parent − 2·child` (close + signed Δ).
    LinPNeg2C,
    /// Linear: Σ `−parent − 2·child`.
    LinNegPNeg2C,
    /// Log costs, no improve clamp: Σ `−ln(1+child) − ln(1+parent)`.
    /// Separates bad children; does not use `max(0, p−c)`. **Default.**
    #[default]
    LogNegPC,
}

/// Which distance(s) feed the match score.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MatchMetric {
    /// Atom-diff cost only.
    Atom,
    /// Formula L1 vs target (includes H) only.
    Formula,
    /// Atom + formula terms (default; **added** under Add/Close/Improve).
    #[default]
    Both,
}

/// Data recipe for the match-family heap score (not SoftStack).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MatchScoreSpec {
    pub combine: MatchCombine,
    pub metric: MatchMetric,
}

impl MatchScoreSpec {
    /// Default live recipe: `−ln(1+child) − ln(1+parent)` (atom+formula).
    pub const fn log_neg_pc() -> Self {
        Self {
            combine: MatchCombine::LogNegPC,
            metric: MatchMetric::Both,
        }
    }

    /// Former default: log-space product of close×imp factors (atom+formula).
    pub const fn add_both() -> Self {
        Self {
            combine: MatchCombine::Add,
            metric: MatchMetric::Both,
        }
    }

    /// Legacy product×both (SCALE/(1+c) mul). Ablation / CLI only.
    pub const fn product_both() -> Self {
        Self {
            combine: MatchCombine::Product,
            metric: MatchMetric::Both,
        }
    }

    /// Bench / CLI label: `add-both`, `close-atom`, `product-both`, …
    pub fn label(self) -> &'static str {
        match (self.combine, self.metric) {
            (MatchCombine::Close, MatchMetric::Atom) => "close-atom",
            (MatchCombine::Close, MatchMetric::Formula) => "close-formula",
            (MatchCombine::Close, MatchMetric::Both) => "close-both",
            (MatchCombine::Improve, MatchMetric::Atom) => "improve-atom",
            (MatchCombine::Improve, MatchMetric::Formula) => "improve-formula",
            (MatchCombine::Improve, MatchMetric::Both) => "improve-both",
            (MatchCombine::Add, MatchMetric::Atom) => "add-atom",
            (MatchCombine::Add, MatchMetric::Formula) => "add-formula",
            (MatchCombine::Add, MatchMetric::Both) => "add-both",
            (MatchCombine::Product, MatchMetric::Atom) => "product-atom",
            (MatchCombine::Product, MatchMetric::Formula) => "product-formula",
            (MatchCombine::Product, MatchMetric::Both) => "product-both",
            (MatchCombine::LinNegC, MatchMetric::Atom) => "lin-neg-c-atom",
            (MatchCombine::LinNegC, MatchMetric::Formula) => "lin-neg-c-formula",
            (MatchCombine::LinNegC, MatchMetric::Both) => "lin-neg-c",
            (MatchCombine::LinNegPC, MatchMetric::Atom) => "lin-neg-pc-atom",
            (MatchCombine::LinNegPC, MatchMetric::Formula) => "lin-neg-pc-formula",
            (MatchCombine::LinNegPC, MatchMetric::Both) => "lin-neg-pc",
            (MatchCombine::LinPNeg2C, MatchMetric::Atom) => "lin-p-neg2c-atom",
            (MatchCombine::LinPNeg2C, MatchMetric::Formula) => "lin-p-neg2c-formula",
            (MatchCombine::LinPNeg2C, MatchMetric::Both) => "lin-p-neg2c",
            (MatchCombine::LinNegPNeg2C, MatchMetric::Atom) => "lin-neg-p-neg2c-atom",
            (MatchCombine::LinNegPNeg2C, MatchMetric::Formula) => "lin-neg-p-neg2c-formula",
            (MatchCombine::LinNegPNeg2C, MatchMetric::Both) => "lin-neg-p-neg2c",
            (MatchCombine::LogNegPC, MatchMetric::Atom) => "log-neg-pc-atom",
            (MatchCombine::LogNegPC, MatchMetric::Formula) => "log-neg-pc-formula",
            (MatchCombine::LogNegPC, MatchMetric::Both) => "log-neg-pc",
        }
    }

    /// Linear cost→score weight pair `(parent_w, child_w)` for [`MatchCombine::Lin*`].
    pub fn linear_weights(self) -> Option<(i64, i64)> {
        match self.combine {
            MatchCombine::LinNegC => Some((0, -1)),
            MatchCombine::LinNegPC => Some((-1, -1)),
            MatchCombine::LinPNeg2C => Some((1, -2)),
            MatchCombine::LinNegPNeg2C => Some((-1, -2)),
            _ => None,
        }
    }

    /// Nine live combines (close|improve|add × atom|formula|both). Legacy
    /// `Product` is opt-in via CLI, not in this matrix.
    pub fn matrix() -> [Self; 9] {
        let mut out = [Self::default(); 9];
        let mut i = 0;
        for combine in [
            MatchCombine::Close,
            MatchCombine::Improve,
            MatchCombine::Add,
        ] {
            for metric in [MatchMetric::Atom, MatchMetric::Formula, MatchMetric::Both] {
                out[i] = Self { combine, metric };
                i += 1;
            }
        }
        out
    }
}

/// Site-application identity for the diversity term: rule + parent site tags +
/// tags minted for atoms the effect added.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct DiversityKey {
    rule: String,
    site_tags: BTreeSet<Tag>,
    added_tags: BTreeSet<Tag>,
}

fn diversity_key_for(
    parent: &ForestMol,
    product: &ForestMol,
    pattern_name: &str,
    site_atoms: &BTreeSet<usize>,
) -> DiversityKey {
    let site_tags: BTreeSet<Tag> = site_atoms
        .iter()
        .filter_map(|&i| parent.tag_of(i))
        .collect();
    let added_tags: BTreeSet<Tag> = (0..product.mol().atom_count())
        .filter_map(|i| product.tag_of(i))
        .filter(|tag| parent.index_of(*tag).is_none())
        .collect();
    DiversityKey {
        rule: pattern_name.to_string(),
        site_tags,
        added_tags,
    }
}

/// Fixed-point `−n` diversity penalty (same scale as ln score terms).
/// Stronger than `−ln(n+1)`: each prior accept costs a full unit.
pub fn diversity_penalty(times_applied: usize) -> i64 {
    -(times_applied as i64) * (MATCH_LOG_SCALE as i64)
}

/// How the find_path frontier ranks walks (score, then `seq`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeapScoreMode {
    /// Soft stack: `search_bias`, site H-progress, `cost_gain`, then `seq`.
    /// Opt-in via `FindPathConfig` / `--score soft`.
    SoftStack,
    /// Match-family score from [`MatchScoreSpec`] (default: log-neg-pc).
    /// Ord key: `(score + diversity, score)`, then `seq`.
    Match(MatchScoreSpec),
}

impl Default for HeapScoreMode {
    fn default() -> Self {
        Self::Match(MatchScoreSpec::log_neg_pc())
    }
}

impl HeapScoreMode {
    /// Default match recipe (`−ln(1+c) − ln(1+p)`, formula+atom).
    pub const fn match_log_neg_pc() -> Self {
        Self::Match(MatchScoreSpec::log_neg_pc())
    }

    /// Former default: log of close×imp product.
    pub const fn match_add() -> Self {
        Self::Match(MatchScoreSpec::add_both())
    }

    /// Alias of [`Self::match_log_neg_pc`].
    pub const fn match_product() -> Self {
        Self::match_log_neg_pc()
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::SoftStack => "soft",
            Self::Match(spec) => spec.label(),
        }
    }

    /// Parse a bench / CLI score label (`log-neg-pc`, `soft`, `add-both`, …).
    pub fn from_label(label: &str) -> Option<Self> {
        if label == "soft" {
            return Some(Self::SoftStack);
        }
        for combine in [
            MatchCombine::Close,
            MatchCombine::Improve,
            MatchCombine::Add,
            MatchCombine::Product,
            MatchCombine::LinNegC,
            MatchCombine::LinNegPC,
            MatchCombine::LinPNeg2C,
            MatchCombine::LinNegPNeg2C,
            MatchCombine::LogNegPC,
        ] {
            for metric in [MatchMetric::Atom, MatchMetric::Formula, MatchMetric::Both] {
                let spec = MatchScoreSpec { combine, metric };
                if spec.label() == label {
                    return Some(Self::Match(spec));
                }
            }
        }
        None
    }
}

/// Closeness scale shared by Product multiply and log-space Add/Close.
const MATCH_CLOSE_SCALE: f64 = 1_000_000.0;
/// Fixed-point for `ln` terms (`round(ln(x) * this)`).
const MATCH_LOG_SCALE: f64 = 1_000_000.0;

/// Fixed-point `ln(x)` for positive factors. Higher is better.
fn ln_fixed(x: f64) -> i64 {
    debug_assert!(x > 0.0);
    (x.ln() * MATCH_LOG_SCALE).round() as i64
}

/// Log-space closeness factor: `ln(SCALE / (1+cost))`.
/// Equals [`neg_log1p_score`] plus a metric-count constant.
pub fn log_close_term(cost: usize) -> i64 {
    ln_fixed(MATCH_CLOSE_SCALE / (1.0 + cost as f64))
}

/// Fixed-point `−ln(1 + cost)` (= `log_close_term` − `ln(SCALE)`).
/// Kept as a named closeness shape; prefer [`log_close_term`] in score.
pub fn neg_log1p_score(cost: usize) -> i64 {
    (-(1.0 + cost as f64).ln() * MATCH_LOG_SCALE).round() as i64
}

/// Log-space improvement factor: `ln(max(0, parent−child)+1)`.
pub fn log_improve_term(parent: usize, child: usize) -> i64 {
    let imp = (parent as i64 - child as i64).max(0) + 1;
    ln_fixed(imp as f64)
}

fn match_factors(
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
    metric: MatchMetric,
) -> (Vec<f64>, Vec<f64>) {
    let use_formula = matches!(metric, MatchMetric::Formula | MatchMetric::Both);
    let use_atom = matches!(metric, MatchMetric::Atom | MatchMetric::Both);
    let mut imps = Vec::new();
    let mut closes = Vec::new();
    if use_formula {
        imps.push(((parent_formula_dist as i64 - child_formula_dist as i64).max(0) + 1) as f64);
        closes.push(MATCH_CLOSE_SCALE / (1.0 + child_formula_dist as f64));
    }
    if use_atom {
        match (parent_atom_cost, child_atom_cost) {
            (Some(p), Some(c)) => {
                imps.push(((p as i64 - c as i64).max(0) + 1) as f64);
                closes.push(MATCH_CLOSE_SCALE / (1.0 + c as f64));
            }
            // Same as Product: missing side → factor 1 (log 0 / mul 1).
            _ => {
                imps.push(1.0);
                closes.push(1.0);
            }
        }
    }
    (imps, closes)
}

fn legacy_product_score(
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
    metric: MatchMetric,
) -> i64 {
    let (imps, closes) = match_factors(
        parent_formula_dist,
        child_formula_dist,
        parent_atom_cost,
        child_atom_cost,
        metric,
    );
    let mut prod = 1i64;
    for x in imps.into_iter().chain(closes) {
        prod = prod.saturating_mul(x as i64);
    }
    prod
}

/// Parent→child match-family heap score from [`MatchScoreSpec`].
///
/// [`MatchCombine::Add`] / Close / Improve use fixed-point `ln` of the **same**
/// positive factors [`MatchCombine::Product`] multiplies — ranking-identical
/// to Product (up to float rounding). Metrics contribute one imp factor and
/// one close factor each when active. A target hit is cost/formula 0.
fn linear_cost_score(
    parent_w: i64,
    child_w: i64,
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
    metric: MatchMetric,
) -> i64 {
    let use_formula = matches!(metric, MatchMetric::Formula | MatchMetric::Both);
    let use_atom = matches!(metric, MatchMetric::Atom | MatchMetric::Both);
    let mut score = 0i64;
    if use_formula {
        score += parent_w * parent_formula_dist as i64 + child_w * child_formula_dist as i64;
    }
    if use_atom {
        match (parent_atom_cost, child_atom_cost) {
            (Some(p), Some(c)) => {
                score += parent_w * p as i64 + child_w * c as i64;
            }
            // Child only: parent term 0 (unknown).
            (None, Some(c)) => {
                score += child_w * c as i64;
            }
            _ => {}
        }
    }
    score
}

/// Σ `−ln(1+child) − ln(1+parent)` over active metrics (no `max(0,p−c)`).
fn log_neg_pc_score(
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
    metric: MatchMetric,
) -> i64 {
    let use_formula = matches!(metric, MatchMetric::Formula | MatchMetric::Both);
    let use_atom = matches!(metric, MatchMetric::Atom | MatchMetric::Both);
    let mut score = 0i64;
    if use_formula {
        score += neg_log1p_score(child_formula_dist) + neg_log1p_score(parent_formula_dist);
    }
    if use_atom {
        match (parent_atom_cost, child_atom_cost) {
            (Some(p), Some(c)) => {
                score += neg_log1p_score(c) + neg_log1p_score(p);
            }
            (None, Some(c)) => {
                score += neg_log1p_score(c);
            }
            _ => {}
        }
    }
    score
}

pub fn hop_match_score(
    spec: MatchScoreSpec,
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
) -> i64 {
    if matches!(spec.combine, MatchCombine::Product) {
        return legacy_product_score(
            parent_formula_dist,
            child_formula_dist,
            parent_atom_cost,
            child_atom_cost,
            spec.metric,
        );
    }
    if let Some((pw, cw)) = spec.linear_weights() {
        return linear_cost_score(
            pw,
            cw,
            parent_formula_dist,
            child_formula_dist,
            parent_atom_cost,
            child_atom_cost,
            spec.metric,
        );
    }
    if matches!(spec.combine, MatchCombine::LogNegPC) {
        return log_neg_pc_score(
            parent_formula_dist,
            child_formula_dist,
            parent_atom_cost,
            child_atom_cost,
            spec.metric,
        );
    }

    let (imps, closes) = match_factors(
        parent_formula_dist,
        child_formula_dist,
        parent_atom_cost,
        child_atom_cost,
        spec.metric,
    );
    let close_sum: i64 = closes.iter().map(|&c| ln_fixed(c)).sum();
    let improve_sum: i64 = imps.iter().map(|&i| ln_fixed(i)).sum();

    match spec.combine {
        MatchCombine::Close => close_sum,
        MatchCombine::Improve => improve_sum,
        MatchCombine::Add => close_sum + improve_sum,
        MatchCombine::Product
        | MatchCombine::LinNegC
        | MatchCombine::LinNegPC
        | MatchCombine::LinPNeg2C
        | MatchCombine::LinNegPNeg2C
        | MatchCombine::LogNegPC => unreachable!("handled above"),
    }
}

/// Convenience: [`MatchScoreSpec::add_both`] (former default; log-Product).
pub fn hop_match_add_score(
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
) -> i64 {
    hop_match_score(
        MatchScoreSpec::add_both(),
        parent_formula_dist,
        child_formula_dist,
        parent_atom_cost,
        child_atom_cost,
    )
}

/// Alias of [`hop_match_add_score`] (old product-both name).
pub fn hop_match_product_score(
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
) -> i64 {
    hop_match_add_score(
        parent_formula_dist,
        child_formula_dist,
        parent_atom_cost,
        child_atom_cost,
    )
}

fn match_score_for(
    mode: HeapScoreMode,
    parent_formula_dist: usize,
    child_formula_dist: usize,
    parent_atom_cost: Option<usize>,
    child_atom_cost: Option<usize>,
) -> i64 {
    match mode {
        HeapScoreMode::SoftStack => 0,
        HeapScoreMode::Match(spec) => hop_match_score(
            spec,
            parent_formula_dist,
            child_formula_dist,
            parent_atom_cost,
            child_atom_cost,
        ),
    }
}

/// Heap entry: mode score(s), then `seq` as a pure tiebreak. Target hits rank
/// via atom_diff / formula closeness (cost 0), not a bool tier or score
/// sentinel. Pop is plain [`BinaryHeap::pop`] — best Ord value only (no
/// DFS/BFS alternation). BinaryHeap is max-heap.
///
/// Match Ord: `(match_priority, match_score)` where
/// `match_priority = match_score + diversity_penalty(n)` and `n` is how often
/// that [`DiversityKey`] has already been accepted. SoftStack ignores both.
#[derive(Clone)]
struct HeapItem {
    mode: HeapScoreMode,
    /// SoftStack: from [`PatternInfo::search_bias`] (higher preferred).
    search_bias: i8,
    /// SoftStack: site H-progress vs parent diff (higher preferred).
    site_progress: i32,
    /// SoftStack: `parent_cost - child_cost` (higher preferred).
    cost_gain: i32,
    /// Match-family base score from [`MatchScoreSpec`] (no diversity).
    match_score: i64,
    /// Match Ord primary: `match_score + −n` (fixed-point) at enqueue / last refresh.
    match_priority: i64,
    /// Hop that created this walk (`None` = root). Match diversity only.
    diversity_key: Option<DiversityKey>,
    seq: usize,
    walk: Walk,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.mode == other.mode
            && self.search_bias == other.search_bias
            && self.site_progress == other.site_progress
            && self.cost_gain == other.cost_gain
            && self.match_score == other.match_score
            && self.match_priority == other.match_priority
            && self.seq == other.seq
    }
}

impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.mode {
            HeapScoreMode::SoftStack => self
                .search_bias
                .cmp(&other.search_bias)
                .then_with(|| self.site_progress.cmp(&other.site_progress))
                .then_with(|| self.cost_gain.cmp(&other.cost_gain)),
            HeapScoreMode::Match(_) => self
                .match_priority
                .cmp(&other.match_priority)
                .then_with(|| self.match_score.cmp(&other.match_score)),
        }
        // Tiebreak only — higher seq preferred among equal scores.
        .then_with(|| self.seq.cmp(&other.seq))
    }
}

fn ha_distance(ha: usize, target_ha: usize) -> usize {
    ha.abs_diff(target_ha)
}

/// Insert into `seen` only when Chematic proves a stable key (fail-closed).
fn remember_seen(seen: &mut HashSet<String>, mol: &ForestMol) -> bool {
    match mol.stable_csmi_key() {
        Some(k) => seen.insert(k.as_ref().to_string()),
        None => false,
    }
}

/// True when `mol` has a stable key already in `seen`. Unstable → never "seen".
fn already_seen(seen: &HashSet<String>, mol: &ForestMol) -> bool {
    mol.stable_csmi_key()
        .is_some_and(|k| seen.contains(k.as_ref()))
}

/// True when ``child``'s stable key already appears on the walk (cycle).
/// Unstable keys do not trigger — fail-closed skip, not a cycle refuse.
fn repeats_ancestor(ancestors: &HashSet<String>, child: &ForestMol) -> bool {
    child
        .stable_csmi_key()
        .is_some_and(|k| ancestors.contains(k.as_ref()))
}

fn with_child_ancestor(ancestors: &HashSet<String>, child: &ForestMol) -> HashSet<String> {
    let mut next = ancestors.clone();
    if let Some(k) = child.stable_csmi_key() {
        next.insert(k.as_ref().to_string());
    }
    next
}

fn root_ancestors(start: &ForestMol) -> HashSet<String> {
    let mut ancestors = HashSet::new();
    if let Some(k) = start.stable_csmi_key() {
        ancestors.insert(k.as_ref().to_string());
    }
    ancestors
}

/// Product index of the atom that carried ``parent_idx``'s forest tag.
fn parent_to_product_idx(
    parent: &ForestMol,
    product: &ForestMol,
    parent_idx: usize,
) -> Option<usize> {
    let tag = parent.tag_of(parent_idx)?;
    product.index_of(tag)
}

/// After DH: each end's heavy neighbors on the **product** match the target.
fn dh_product_ends_match_target(
    parent: &ForestMol,
    product: &ForestMol,
    ends: (usize, usize),
    target: &crate::Molecule,
) -> bool {
    let (a, b) = ends;
    let mut product_ends = Vec::with_capacity(2);
    for end in [a, b] {
        let Some(idx) = parent_to_product_idx(parent, product, end) else {
            return false;
        };
        product_ends.push(idx);
    }
    crate::atom_diff::dh_product_ends_match(product.mol(), &product_ends, target)
}

/// Tagged emission: ForestMol products + elementary plan + cleavage site data.
#[derive(Clone)]
struct ForestEmission {
    site: usize,
    site_orbit: Vec<usize>,
    /// Discovery site atoms (Python frozenset site for CleavageSide).
    site_atoms: BTreeSet<usize>,
    cleaves: bool,
    /// Cross-rule Or fold signature (from [`PatternInfo::cleave_side_group`]).
    cleave_side_sig: CleaveSideSig,
    adds_oxygen: bool,
    removes_oxygen: bool,
    oxygen_site: OxygenSite,
    pattern_name: String,
    /// From [`PatternInfo::search_bias`] (pair: min of ends).
    search_bias: i8,
    /// Site H-progress vs parent diff at emit (0 if no diff). Soft heap score.
    site_progress: i32,
    /// Parent end atoms for post-application DH neighbor match (None if not DH).
    dh_ends: Option<(usize, usize)>,
    rule_path: Vec<Option<String>>,
    products: Vec<ForestMol>,
    plan: Vec<CanonicalStep>,
}

impl ForestEmission {
    fn product_csmis(&self) -> Vec<String> {
        self.products
            .iter()
            .map(|p| p.csmi().as_ref().to_string())
            .collect()
    }

    fn cleave_fold_key(&self) -> CleaveFoldKey {
        self.cleave_side_sig.fold_key(&self.product_csmis())
    }
}

/// Keep fragments that are worth continuing toward ``target``.
///
/// On a bifurcation with a parent diff, every product that passes the cleavage
/// expand gate (strict cost drop via tag-lift when possible, ha ≥ target) is
/// kept — **both** sides if both match. Uses tagged [`ForestMol`] products;
/// never re-parses SMILES.
fn keep_fragments(
    parent: &ForestMol,
    products: &[ForestMol],
    target_csmi: &str,
    target_ha: usize,
    parent_diff: Option<&crate::atom_diff::AtomDiff>,
    target_mol: Option<&crate::Molecule>,
    mcs: &mut LiftMcsCounters<'_>,
) -> Vec<(ForestMol, Vec<String>, Option<crate::atom_diff::AtomDiff>)> {
    if products.is_empty() {
        return Vec::new();
    }
    if products.len() >= 2 {
        if let (Some(pdiff), Some(tmol)) = (parent_diff, target_mol) {
            let mut kept = Vec::new();
            for (i, child) in products.iter().enumerate() {
                let csmi = child.csmi();
                let is_hit = csmi.as_ref() == target_csmi;
                if !is_hit && child.heavy_atom_count() < target_ha {
                    continue;
                }
                let child_diff = crate::atom_diff::atom_diff_after_cleavage_tracked(
                    parent,
                    pdiff,
                    child,
                    tmol,
                    Some(&mut mcs.fallback),
                    Some(&mut mcs.rematch),
                );
                if is_hit || child_diff.cost() < pdiff.cost() {
                    let mut sides = Vec::new();
                    for (j, other) in products.iter().enumerate() {
                        if j != i && other.csmi().as_ref() != csmi.as_ref() {
                            sides.push(other.csmi().as_ref().to_string());
                        }
                    }
                    kept.push((child.clone(), sides, Some(child_diff)));
                }
            }
            if !kept.is_empty() {
                return kept;
            }
        }
    }
    // Closest fragment by HA/CSMI; still attach a lifted diff when possible so
    // enqueue does not rematch MCS.
    keep_fragment(products, target_csmi, target_ha)
        .into_iter()
        .map(|(mol, sides)| {
            let child_diff = match (parent_diff, target_mol) {
                (Some(pdiff), Some(tmol)) => {
                    Some(crate::atom_diff::atom_diff_after_cleavage_tracked(
                        parent,
                        pdiff,
                        &mol,
                        tmol,
                        Some(&mut mcs.fallback),
                        Some(&mut mcs.rematch),
                    ))
                }
                _ => None,
            };
            (mol, sides, child_diff)
        })
        .collect()
}

/// Counters for tag-lift vs MCS on cleavage / child diffs.
struct LiftMcsCounters<'a> {
    fallback: &'a mut usize,
    rematch: &'a mut usize,
}

/// Keep the fragment closest to ``target``. Prefer exact CSMI hit.
fn keep_fragment(
    products: &[ForestMol],
    target_csmi: &str,
    target_ha: usize,
) -> Option<(ForestMol, Vec<String>)> {
    if products.is_empty() {
        return None;
    }
    let mut best: Option<(usize, Rc<str>, usize)> = None;
    for (i, mol) in products.iter().enumerate() {
        let csmi = mol.csmi();
        let cost = if csmi.as_ref() == target_csmi {
            0
        } else {
            1 + ha_distance(mol.heavy_atom_count(), target_ha)
        };
        match &best {
            None => best = Some((i, Rc::clone(&csmi), cost)),
            Some((_, _, best_cost)) if cost < *best_cost => {
                best = Some((i, Rc::clone(&csmi), cost));
            }
            Some((_, kept, best_cost)) if cost == *best_cost && csmi.as_ref() < kept.as_ref() => {
                best = Some((i, Rc::clone(&csmi), cost));
            }
            _ => {}
        }
    }
    let (best_i, kept_csmi, _) = best?;
    let kept = products[best_i].clone();
    let mut sides = Vec::new();
    for (i, mol) in products.iter().enumerate() {
        if i == best_i {
            continue;
        }
        let csmi = mol.csmi();
        if csmi.as_ref() != kept_csmi.as_ref() {
            sides.push(csmi.as_ref().to_string());
        }
    }
    Some((kept, sides))
}

fn closer(parent_ha: usize, child_ha: usize, target_ha: usize, target_hit: bool) -> bool {
    // Refuse only walks that grow more distant in heavy-atom count. Equal
    // distance (same-size DH / oxidation hops) must stay open — a strict `<`
    // drops propane→propene→epoxide. When atom_diff is on, this is the cheap
    // enqueue gate; full cost closer runs on pop (lazy).
    target_hit || ha_distance(child_ha, target_ha) <= ha_distance(parent_ha, target_ha)
}

fn cost_closer(parent_cost: usize, child_cost: usize, target_hit: bool) -> bool {
    target_hit || child_cost < parent_cost
}

/// Parent-relative atom_diff cost drop for the heap. Higher preferred.
/// ≤0 counters DFS vs positive peers. A target hit is just child_cost 0.
/// Requires known costs — callers compute child diff when parent cost is known.
fn hop_cost_gain(parent_cost: Option<usize>, child_cost: Option<usize>) -> i32 {
    match (parent_cost, child_cost) {
        (Some(p), Some(c)) => p as i32 - c as i32,
        _ => 0,
    }
}

/// Search bounds. Defaults match Python `find_path` knobs.
#[derive(Clone, Copy, Debug)]
pub struct FindPathConfig {
    pub max_paths: usize,
    pub max_nodes: usize,
    /// When true, build an [`crate::atom_diff::AtomDiff`] once per expansion
    /// and refuse candidates that cannot help (no filter closures).
    pub use_atom_diff: bool,
    /// When true with `use_atom_diff`, defer full cost closer to pop: enqueue
    /// every kept fragment, verify `cost() < parent_cost` before expanding.
    /// Default **false**: child cost is already resolved at enqueue for heap
    /// `cost_gain`, so refuse non-closer there (DFS among improvers only).
    pub lazy_closer: bool,
    /// Frontier ranking (score, then seq).
    pub heap_score: HeapScoreMode,
    /// Yield-drop remapped free-step twins via [`Deps::same_rule_maybe_skeleton`]
    /// beside exact [`Deps::same_linearizations`]. Ablation: set false for
    /// exact-only.
    pub drop_skeleton_twins: bool,
    /// Match diversity: count accepted applications per site key, add
    /// fixed-point `−n` into the heap primary (`score + diversity`). On pop,
    /// refresh or accept. SoftStack ignores. Default **false** (opt-in).
    pub diversity: bool,
}

impl Default for FindPathConfig {
    fn default() -> Self {
        Self {
            max_paths: 1,
            max_nodes: 800,
            // Match Python live `use_filters=True`.
            use_atom_diff: true,
            lazy_closer: false,
            heap_score: HeapScoreMode::match_log_neg_pc(),
            drop_skeleton_twins: true,
            diversity: false,
        }
    }
}

fn accept_all_candidates(_c: &Candidate) -> bool {
    true
}

/// [`FindPath`] with the built-in accept-all keep predicate.
pub type OpenFindPath<'a, 'b> = FindPath<'a, 'b, fn(&Candidate) -> bool>;

/// HEURISTICS: a later walk that is only a reordering of an already-yielded
/// [`Deps`] is not a new path. Optionally also drop remapped-index free-step
/// twins ([`Deps::same_rule_maybe_skeleton`]) when `drop_skeleton_twins`.
/// Dominated-extension is **not** a yield drop (over-collapsed multipath); it
/// only bumps [`PathCounters::signal_contained_plan`].
fn plan_already_yielded(found: &[PathOutcome], plan: &Deps, drop_skeleton_twins: bool) -> bool {
    found.iter().any(|h| {
        h.plan.same_linearizations(plan)
            || (drop_skeleton_twins && h.plan.same_rule_maybe_skeleton(plan))
    })
}

/// Classify why a hit would be a duplicate (exact first, else skeleton).
fn duplicate_plan_kind(
    found: &[PathOutcome],
    plan: &Deps,
    drop_skeleton_twins: bool,
) -> Option<DuplicatePlanKind> {
    if found.iter().any(|h| h.plan.same_linearizations(plan)) {
        return Some(DuplicatePlanKind::Exact);
    }
    if drop_skeleton_twins && found.iter().any(|h| h.plan.same_rule_maybe_skeleton(plan)) {
        return Some(DuplicatePlanKind::Skeleton);
    }
    None
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DuplicatePlanKind {
    Exact,
    Skeleton,
}

fn record_yield_plan_signals(
    counters: &mut PathCounters,
    found: &[PathOutcome],
    plan: &Deps,
    drop_skeleton_twins: bool,
) {
    match duplicate_plan_kind(found, plan, drop_skeleton_twins) {
        Some(DuplicatePlanKind::Exact) => {
            counters.dropped_duplicate_plan += 1;
            counters.dropped_exact_plan += 1;
            return;
        }
        Some(DuplicatePlanKind::Skeleton) => {
            counters.dropped_duplicate_plan += 1;
            counters.dropped_skeleton_twin += 1;
            return;
        }
        None => {}
    }
    if found.iter().any(|h| h.plan.dominates_extension_of(plan)) {
        counters.signal_contained_plan += 1;
    }
}

fn oxygen_site_from_parts(site: usize, orbit: &[usize], atoms: BTreeSet<usize>) -> OxygenSite {
    OxygenSite {
        site,
        orbit: orbit.to_vec(),
        atoms,
    }
}

/// Hydration↔dehydration is circular only when site sets match. Multi-atom
/// sites (e.g. beta-elim with adjacent C) compare by set equality — any atom
/// difference is allowed. Singletons use top / unique-edit orbits
/// ([`same_site_orbit`]).
fn oxygen_sites_equal(a: &OxygenSite, b: &OxygenSite) -> bool {
    if a.atoms.len() > 1 || b.atoms.len() > 1 {
        return a.atoms == b.atoms;
    }
    crate::same_site_orbit(a.site, &a.orbit, b.site, &b.orbit)
}

fn undoes_oxygen_edit(
    effect: &crate::pattern::Effect,
    site: &OxygenSite,
    o_added: &[OxygenSite],
    o_removed: &[OxygenSite],
) -> bool {
    if crate::atom_diff::effect_removes_oxygen(effect) {
        return o_added.iter().any(|prior| oxygen_sites_equal(prior, site));
    }
    if crate::atom_diff::effect_adds_oxygen(effect) {
        return o_removed
            .iter()
            .any(|prior| oxygen_sites_equal(prior, site));
    }
    false
}

fn child_oxygen_lists(
    parent: &Walk,
    site: &OxygenSite,
    adds_oxygen: bool,
    removes_oxygen: bool,
) -> (Vec<OxygenSite>, Vec<OxygenSite>) {
    let mut o_added = parent.o_added.clone();
    let mut o_removed = parent.o_removed.clone();
    if adds_oxygen {
        o_added.push(site.clone());
    }
    if removes_oxygen {
        o_removed.push(site.clone());
    }
    (o_added, o_removed)
}

/// Yield walks that turn ``reactant`` into ``target``.
///
/// Pull iterator (Python generator parity): each [`Iterator::next`] resumes the
/// search until the next [`PathOutcome`]. Nested sets keep step namespaces on
/// each `rule_path`. Callers that want a list use `.collect()`.
pub fn find_path<'a, 'b>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'a, 'b>, ForestError> {
    find_path_with(
        reactant,
        target,
        ruleset,
        counters,
        FindPathConfig::default(),
        accept_all_candidates,
    )
}

/// [`find_path`] with the Phase-I default catalog.
pub fn find_path_default<'b>(
    reactant: &str,
    target: &str,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'static, 'b>, ForestError> {
    // Leak-free: borrow the process-wide default via once / static ruleset.
    find_path(reactant, target, default_ruleset_ref(), counters)
}

fn default_ruleset_ref() -> &'static RuleSet {
    use std::sync::OnceLock;
    static DEFAULT: OnceLock<RuleSet> = OnceLock::new();
    DEFAULT.get_or_init(default_ruleset)
}

/// [`find_path`] with atom-diff candidate gating (no filter closures).
pub fn find_path_diff<'a, 'b>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
) -> Result<OpenFindPath<'a, 'b>, ForestError> {
    find_path_with(
        reactant,
        target,
        ruleset,
        counters,
        FindPathConfig {
            use_atom_diff: true,
            ..FindPathConfig::default()
        },
        accept_all_candidates,
    )
}

/// Same as [`find_path`], with caller bounds and a candidate predicate.
///
/// `keep` reads the triple (site, pattern, parent) **before** materialize.
/// Dropped candidates do not count as `mol_edits`.
pub fn find_path_with<'a, 'b, K>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    config: FindPathConfig,
    keep: K,
) -> Result<FindPath<'a, 'b, K>, ForestError>
where
    K: Fn(&Candidate) -> bool,
{
    let start = ForestMol::parse(reactant)?;
    let target_csmi = canon_of(target)?;
    let target_mol = parse_mol(&target_csmi)?;
    let target_ha = target_mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    let mut seen = HashSet::new();
    remember_seen(&mut seen, &start);
    let ancestors = root_ancestors(&start);
    let target_formula = crate::forest::molecule_formula(&target_mol);
    let start_formula_dist = crate::forest::formula_l1(&start.formula(), &target_formula);
    let match_score = match_score_for(
        config.heap_score,
        start_formula_dist,
        start_formula_dist,
        None,
        None,
    );
    heap.push(HeapItem {
        mode: config.heap_score,
        search_bias: 0,
        site_progress: 0,
        cost_gain: 0,
        match_score,
        match_priority: match_score,
        diversity_key: None,
        seq,
        walk: Walk {
            mol: start,
            steps: Vec::new(),
            plan: Vec::new(),
            maybe: Vec::new(),
            opens: Vec::new(),
            o_added: Vec::new(),
            o_removed: Vec::new(),
            ancestors,
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    Ok(FindPath {
        ruleset,
        counters,
        keep,
        config,
        target_csmi,
        target_mol,
        target_formula,
        target_ha,
        heap,
        seq,
        seen,
        diversity_counts: HashMap::new(),
        yielded: Vec::new(),
        done: false,
    })
}

/// Pull search over a [`RuleSet`]: yields one [`PathOutcome`] per `next`.
///
/// Expand pulls candidates, buffers survivors only for `order_key` sort, then
/// materializes one emission at a time (never a full emission list).
pub struct FindPath<'a, 'b, K> {
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    keep: K,
    config: FindPathConfig,
    target_csmi: String,
    target_mol: crate::Molecule,
    target_formula: crate::forest::Formula,
    target_ha: usize,
    heap: BinaryHeap<HeapItem>,
    seq: usize,
    seen: HashSet<String>,
    /// Accepted applications per diversity key (Match + `config.diversity`).
    diversity_counts: HashMap<DiversityKey, usize>,
    /// Already-yielded hits (for [`plan_already_yielded`] only).
    yielded: Vec<PathOutcome>,
    done: bool,
}

impl<K> FindPath<'_, '_, K>
where
    K: Fn(&Candidate) -> bool,
{
    /// Drain the search into a list (Python `list(find_path(...))`).
    pub fn collect_all(self) -> Result<Vec<PathOutcome>, ForestError> {
        self.collect()
    }
}

impl<K> Iterator for FindPath<'_, '_, K>
where
    K: Fn(&Candidate) -> bool,
{
    type Item = Result<PathOutcome, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let FindPathConfig {
            max_paths,
            max_nodes,
            use_atom_diff,
            lazy_closer,
            heap_score,
            drop_skeleton_twins,
            diversity,
        } = self.config;

        while let Some(mut item) = self.heap.pop() {
            if self.yielded.len() >= max_paths || self.counters.nodes >= max_nodes {
                break;
            }
            // Lazy diversity: priority must still equal score + −n (fixed-point).
            if diversity {
                if let (HeapScoreMode::Match(_), Some(key)) =
                    (item.mode, item.diversity_key.clone())
                {
                    let n = self.diversity_counts.get(&key).copied().unwrap_or(0);
                    let fresh = item.match_score.saturating_add(diversity_penalty(n));
                    if fresh != item.match_priority {
                        item.match_priority = fresh;
                        self.heap.push(item);
                        self.counters.diversity_repush += 1;
                        continue;
                    }
                    *self.diversity_counts.entry(key).or_insert(0) += 1;
                }
            }
            let walk = item.walk;
            let here = walk.mol.csmi();
            if here.as_ref() == self.target_csmi.as_str() {
                self.counters.nodes += 1;
                let plan = as_deps(walk.plan).with_maybe(Maybe::new(walk.maybe));
                record_yield_plan_signals(self.counters, &self.yielded, &plan, drop_skeleton_twins);
                if plan_already_yielded(&self.yielded, &plan, drop_skeleton_twins) {
                    continue;
                }
                let outcome = PathOutcome {
                    steps: walk.steps,
                    plan,
                    smiles: here.as_ref().to_string(),
                };
                self.yielded.push(outcome.clone());
                return Some(Ok(outcome));
            }

            let diff = if use_atom_diff {
                let d = match walk.diff {
                    Some(ref d) => d.clone(),
                    None => crate::atom_diff::atom_diff(walk.mol.mol(), &self.target_mol),
                };
                if lazy_closer {
                    if let Some(pc) = walk.parent_cost {
                        if !cost_closer(pc, d.cost(), false) {
                            continue;
                        }
                    }
                }
                Some(d)
            } else {
                None
            };
            self.counters.nodes += 1;
            self.counters.expansions += 1;
            let parent_cost = diff.as_ref().map(|d| d.cost());
            let parent_ha = walk.mol.heavy_atom_count();
            let mut hits_from_here = 0usize;
            // Cross-rule cleave Or: enqueue each (fold_key, continue_csmi) once.
            let mut seen_cleave_continues: HashSet<(CleaveFoldKey, String)> = HashSet::new();
            let mut unstable_csmi = 0usize;
            // Local: Expand already borrows `self.counters` for the loop.
            let mut mcs_lift_fb = 0usize;
            let mut mcs_lift_rm = 0usize;

            let expand = match Expand::new(
                self.counters,
                ExpandInput {
                    ruleset: self.ruleset,
                    parent: &walk.mol,
                    target: &self.target_mol,
                    keep: &self.keep,
                    diff: diff.as_ref(),
                    o_added: &walk.o_added,
                    o_removed: &walk.o_removed,
                },
            ) {
                Ok(e) => e,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            };

            for emission in expand {
                let emission = match emission {
                    Ok(e) => e,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                let keeps = keep_fragments(
                    &walk.mol,
                    &emission.products,
                    &self.target_csmi,
                    self.target_ha,
                    diff.as_ref(),
                    Some(&self.target_mol),
                    &mut LiftMcsCounters {
                        fallback: &mut mcs_lift_fb,
                        rematch: &mut mcs_lift_rm,
                    },
                );
                let cleave_key = if emission.cleaves && emission.products.len() >= 2 {
                    Some(emission.cleave_fold_key())
                } else {
                    None
                };
                for (kept, sides, lifted_diff) in keeps {
                    let kept_csmi = kept.csmi().as_ref().to_string();
                    if let Some(key) = &cleave_key {
                        if !seen_cleave_continues.insert((key.clone(), kept_csmi.clone())) {
                            // Another arm (possibly another rule) already queued
                            // this continuation for the same Or bag.
                            continue;
                        }
                    }
                    let child_ha = kept.heavy_atom_count();
                    let target_hit = kept_csmi == self.target_csmi;

                    let mut child_diff = if use_atom_diff {
                        lifted_diff.or_else(|| {
                            diff.as_ref().and_then(|parent_d| {
                                crate::atom_diff::try_atom_diff_for_child_tracked(
                                    &walk.mol,
                                    parent_d,
                                    &kept,
                                    &self.target_mol,
                                    Some(&mut mcs_lift_rm),
                                )
                            })
                        })
                    } else {
                        None
                    };

                    // Known child cost for heap cost_gain (lack of improvement
                    // counters DFS). Reused on walk so pop does not re-MCS.
                    // Prefer cost-0 tag-lift+extend; else MCS (counted rematch).
                    if use_atom_diff && child_diff.is_none() && parent_cost.is_some() {
                        mcs_lift_fb += 1;
                        child_diff =
                            Some(crate::atom_diff::atom_diff(kept.mol(), &self.target_mol));
                    }

                    let allow = if use_atom_diff {
                        if lazy_closer {
                            true
                        } else if let Some(pc) = parent_cost {
                            cost_closer(pc, child_diff.as_ref().unwrap().cost(), target_hit)
                        } else {
                            true
                        }
                    } else {
                        closer(parent_ha, child_ha, self.target_ha, target_hit)
                    };
                    if !allow {
                        continue;
                    }
                    // DH / QF path ends: each end must match target neighbors
                    // **after** the edit (product), not on the reactant before.
                    if let Some(ends) = emission.dh_ends {
                        if !dh_product_ends_match_target(&walk.mol, &kept, ends, &self.target_mol) {
                            continue;
                        }
                    }
                    if repeats_ancestor(&walk.ancestors, &kept) {
                        continue;
                    }
                    if already_seen(&self.seen, &kept) && !target_hit {
                        continue;
                    }
                    if kept.stable_csmi_key().is_none() {
                        unstable_csmi += 1;
                    }
                    remember_seen(&mut self.seen, &kept);

                    let (child_maybe, child_opens) = accumulate_maybe(
                        &walk.maybe,
                        &walk.opens,
                        &emission.site_atoms,
                        emission.cleaves,
                        emission.products.len(),
                        &sides,
                    );
                    let (o_added, o_removed) = child_oxygen_lists(
                        &walk,
                        &emission.oxygen_site,
                        emission.adds_oxygen,
                        emission.removes_oxygen,
                    );
                    let mut steps = walk.steps.clone();
                    steps.push(PathStep::from_emission(&emission, kept_csmi.clone(), sides));
                    let mut plan = walk.plan.clone();
                    plan.extend(emission.plan.iter().cloned());
                    let cost_gain =
                        hop_cost_gain(parent_cost, child_diff.as_ref().map(|d| d.cost()));
                    let parent_f =
                        crate::forest::formula_l1(&walk.mol.formula(), &self.target_formula);
                    let child_f = crate::forest::formula_l1(&kept.formula(), &self.target_formula);
                    let match_score = match_score_for(
                        heap_score,
                        parent_f,
                        child_f,
                        parent_cost,
                        child_diff.as_ref().map(|d| d.cost()),
                    );
                    let diversity_key =
                        if diversity && matches!(heap_score, HeapScoreMode::Match(_)) {
                            Some(diversity_key_for(
                                &walk.mol,
                                &kept,
                                &emission.pattern_name,
                                &emission.site_atoms,
                            ))
                        } else {
                            None
                        };
                    let match_priority = match &diversity_key {
                        Some(key) => {
                            let n = self.diversity_counts.get(key).copied().unwrap_or(0);
                            match_score.saturating_add(diversity_penalty(n))
                        }
                        None => match_score,
                    };
                    let ancestors = with_child_ancestor(&walk.ancestors, &kept);
                    self.heap.push(HeapItem {
                        mode: heap_score,
                        search_bias: emission.search_bias,
                        site_progress: emission.site_progress,
                        cost_gain,
                        match_score,
                        match_priority,
                        diversity_key,
                        seq: self.seq,
                        walk: Walk {
                            mol: kept,
                            steps,
                            plan,
                            maybe: child_maybe,
                            opens: child_opens,
                            o_added,
                            o_removed,
                            ancestors,
                            parent_cost,
                            diff: child_diff,
                        },
                    });
                    self.seq += 1;
                    if target_hit {
                        hits_from_here += 1;
                        if self.yielded.len() + hits_from_here >= max_paths {
                            break;
                        }
                    }
                }
                if self.yielded.len() + hits_from_here >= max_paths {
                    break;
                }
            }
            self.counters.unstable_csmi_key += unstable_csmi;
            self.counters.mcs_lift_fallback += mcs_lift_fb;
            self.counters.mcs_lift_rematch += mcs_lift_rm;
        }

        self.done = true;
        None
    }
}

/// Python walk update: ring-open accumulates `opens`; bifurcation appends bags.
fn accumulate_maybe(
    parent_maybe: &[CleavageSide],
    parent_opens: &[BTreeSet<usize>],
    site_atoms: &BTreeSet<usize>,
    cleaves: bool,
    n_products: usize,
    discarded: &[String],
) -> (Vec<CleavageSide>, Vec<BTreeSet<usize>>) {
    if n_products == 1 && cleaves {
        let mut opens = parent_opens.to_vec();
        opens.push(site_atoms.clone());
        (parent_maybe.to_vec(), opens)
    } else if n_products > 1 {
        let mut maybe = parent_maybe.to_vec();
        for side in discarded {
            maybe.push(CleavageSide::new(
                site_atoms.iter().copied(),
                side.clone(),
                parent_opens.iter().cloned(),
            ));
        }
        (maybe, parent_opens.to_vec())
    } else {
        (parent_maybe.to_vec(), parent_opens.to_vec())
    }
}

fn candidate_site_atoms(candidate: &Candidate) -> BTreeSet<usize> {
    let mut atoms: BTreeSet<usize> = candidate
        .pattern
        .site_map
        .iter()
        .filter_map(|m| candidate.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.insert(candidate.site);
    }
    atoms
}

/// One DFS frame while streaming ResonancePair emissions.
struct PairFrame<'a> {
    set: &'a RuleSet,
    member_i: usize,
}

/// Pending pair after filter + order_key sort (discovery only).
struct PendingPair<'a> {
    pair: crate::pair_edit::PairCandidate,
    set: &'a RuleSet,
}

/// Inputs for one expand (bundled so `Expand::new` stays under clippy's arity cap).
struct ExpandInput<'a, K> {
    ruleset: &'a RuleSet,
    parent: &'a ForestMol,
    target: &'a crate::Molecule,
    keep: &'a K,
    diff: Option<&'a crate::atom_diff::AtomDiff>,
    o_added: &'a [OxygenSite],
    o_removed: &'a [OxygenSite],
}

/// Pull tagged emissions one at a time.
///
/// Candidate survivors may be buffered for `order_key` sort (discovery, no
/// edit). Materialize runs per `next`. Pair leaves are walked depth-first the
/// same way — never a full emission `Vec`.
struct Expand<'a, K> {
    parent: &'a ForestMol,
    target: &'a crate::Molecule,
    counters: &'a mut PathCounters,
    keep: &'a K,
    diff: Option<&'a crate::atom_diff::AtomDiff>,
    o_added: &'a [OxygenSite],
    o_removed: &'a [OxygenSite],
    deferred: std::vec::IntoIter<Candidate>,
    pair_stack: Vec<PairFrame<'a>>,
    pair_pending: std::vec::IntoIter<PendingPair<'a>>,
    done: bool,
}

impl<'a, K> Expand<'a, K>
where
    K: Fn(&Candidate) -> bool,
{
    fn new(counters: &'a mut PathCounters, input: ExpandInput<'a, K>) -> Result<Self, ForestError> {
        let ExpandInput {
            ruleset,
            parent,
            target,
            keep,
            diff,
            o_added,
            o_removed,
        } = input;
        let mol = parent.mol();
        // Pull candidates; buffer only survivors for order_key sort.
        let mut deferred = Vec::new();
        for c in ruleset.candidates(mol) {
            let c = c?;
            let atoms = candidate_site_atoms(&c);
            let oxy = oxygen_site_from_parts(c.site, &c.orbit, atoms);
            if undoes_oxygen_edit(&c.pattern.effect, &oxy, o_added, o_removed) {
                counters.blocked_circular_oxygen += 1;
                continue;
            }
            let ok = if let Some(d) = diff {
                keep(&c)
                    && crate::atom_diff::candidate_could_help_on(&c, d, Some(mol), Some(target))
            } else {
                keep(&c)
            };
            if ok {
                deferred.push(c);
            }
        }
        if let Some(d) = diff {
            deferred.sort_by_key(|cand| {
                let (a, b, cname, h_prog, pname) =
                    crate::atom_diff::candidate_order_key_on(cand, d, Some(mol), Some(target));
                // Higher search_bias first (negated so sort ascending prefers high).
                // Then site H-progress (already negated in key: applying helps).
                (a, b, cname, -cand.pattern.search_bias, h_prog, pname)
            });
        }
        Ok(Self {
            parent,
            target,
            counters,
            keep,
            diff,
            o_added,
            o_removed,
            deferred: deferred.into_iter(),
            pair_stack: vec![PairFrame {
                set: ruleset,
                member_i: 0,
            }],
            pair_pending: Vec::new().into_iter(),
            done: false,
        })
    }

    fn emit_candidate(
        &mut self,
        candidate: &Candidate,
    ) -> Result<Option<ForestEmission>, ForestError> {
        let pieces = candidate.materialize_mols(self.parent.mol())?;
        if pieces.is_empty() {
            return Ok(None);
        }
        if let Some(detail) = crate::formula_check::check_effect_delta_formula(
            self.parent.mol(),
            &candidate.pattern.effect,
            &pieces,
            &candidate.pattern.name,
        ) {
            self.counters.record_formula_delta_mismatch(detail);
        }
        let products: Vec<ForestMol> = pieces
            .into_iter()
            .map(|piece| self.parent.adopt_product(piece))
            .collect();
        let site_atoms = candidate_site_atoms(candidate);
        let atoms: Vec<usize> = site_atoms.iter().copied().collect();
        let site_progress = self
            .diff
            .map(|d| {
                crate::atom_diff::site_h_progress_best_placement(
                    &candidate.pattern.effect,
                    &atoms,
                    &[],
                    d,
                    self.parent.mol(),
                    self.target,
                )
            })
            .unwrap_or(0);
        let plan = match candidate
            .leaf_rule()
            .and_then(crate::rules::leaf_rule)
            .filter(|leaf| leaf.has_plan_hook())
        {
            Some(leaf) => leaf.canonical_plan(self.parent.mol(), &atoms, None),
            None => candidate.identity_plan_with_gens(
                &self.parent.atom_bond_generators(),
                self.parent.mol().atom_count(),
                self.parent.mol(),
            ),
        };
        Ok(Some(ForestEmission {
            site: candidate.site,
            site_orbit: candidate.orbit.clone(),
            site_atoms: site_atoms.clone(),
            cleaves: candidate.pattern.effect.cleaves,
            cleave_side_sig: candidate.pattern.cleave_side_sig(),
            adds_oxygen: crate::atom_diff::effect_adds_oxygen(&candidate.pattern.effect),
            removes_oxygen: crate::atom_diff::effect_removes_oxygen(&candidate.pattern.effect),
            oxygen_site: oxygen_site_from_parts(candidate.site, &candidate.orbit, site_atoms),
            pattern_name: candidate.pattern.name.clone(),
            search_bias: candidate.pattern.search_bias,
            site_progress,
            dh_ends: None,
            rule_path: candidate.rule_path.clone(),
            products,
            plan,
        }))
    }

    fn emit_pair(
        &mut self,
        pending: &PendingPair<'a>,
    ) -> Result<Option<ForestEmission>, ForestError> {
        let mol = self.parent.mol();
        let pair = &pending.pair;
        // Python ResonancePair bumps mol_edits only after alternating paths exist.
        let pieces = pair.materialize_mols(mol)?;
        if pieces.is_empty() {
            return Ok(None);
        }
        self.counters.mol_edits += 1;
        if let Some(detail) = crate::formula_check::check_effect_delta_formula(
            mol,
            &pair.effect,
            &pieces,
            &pair.pattern_name,
        ) {
            self.counters.record_formula_delta_mismatch(detail);
        }
        let products: Vec<ForestMol> = pieces
            .into_iter()
            .map(|piece| self.parent.adopt_product(piece))
            .collect();
        let site_atoms = pair.plan_site_atoms();
        let site_atoms_set: BTreeSet<usize> = site_atoms.iter().copied().collect();
        let ends = [&pair.left.effect, &pair.right.effect];
        let plan = pending.set.canonical_plan(mol, &site_atoms, Some(&ends));
        let site_progress = self
            .diff
            .map(|d| crate::atom_diff::pair_site_h_progress(pair, d, mol, self.target))
            .unwrap_or(0);
        let dh_ends = if crate::atom_diff::is_dehydrogenation_effect(&pair.effect) {
            pair.end_atoms()
        } else {
            None
        };
        Ok(Some(ForestEmission {
            site: pair.site,
            site_orbit: vec![pair.site],
            site_atoms: site_atoms_set.clone(),
            cleaves: pair.effect.cleaves,
            cleave_side_sig: {
                let left = pair.left.cleave_side_sig();
                let right = pair.right.cleave_side_sig();
                if left == right {
                    left
                } else {
                    CleaveSideSig::Ungrouped
                }
            },
            adds_oxygen: crate::atom_diff::effect_adds_oxygen(&pair.effect),
            removes_oxygen: crate::atom_diff::effect_removes_oxygen(&pair.effect),
            oxygen_site: oxygen_site_from_parts(pair.site, &[pair.site], site_atoms_set),
            pattern_name: pair.pattern_name.clone(),
            search_bias: pair.left.search_bias.min(pair.right.search_bias),
            site_progress,
            dh_ends,
            rule_path: pending.pair.rule_path.clone(),
            products,
            plan,
        }))
    }

    fn load_leaf_pairs(&mut self, set: &'a RuleSet) -> Result<(), ForestError> {
        let mol = self.parent.mol();
        let mut pairs = set.pair_candidates_leaf(mol)?;
        if let Some(d) = self.diff {
            let mut blocked = 0usize;
            pairs.retain(|pair| {
                if !keep_pair(pair, self.keep) {
                    return false;
                }
                let atoms: BTreeSet<usize> = pair.plan_site_atoms().into_iter().collect();
                let oxy = oxygen_site_from_parts(pair.site, &[pair.site], atoms);
                if undoes_oxygen_edit(&pair.effect, &oxy, self.o_added, self.o_removed) {
                    blocked += 1;
                    return false;
                }
                crate::atom_diff::pair_could_help(pair, d, mol, self.target)
            });
            self.counters.blocked_circular_oxygen += blocked;
            pairs.sort_by_key(|p| {
                let cleave = if p.effect.cleaves { 0u8 } else { 1 };
                let dear = if p.effect.dearomatizes { 0u8 } else { 1 };
                let oxy = if p.effect.adds.as_deref().is_some_and(|a| a.contains('O')) {
                    0u8
                } else {
                    1
                };
                let want_cleave = d.target_smaller() || d.has_cleavage();
                let want_dear = !d.loses_aromaticity.is_empty();
                let want_oxy = crate::atom_diff::any_needs_oxygen(self.target, d);
                (
                    if want_cleave { cleave } else { 0 },
                    if want_dear { dear } else { 0 },
                    if want_oxy { oxy } else { 0 },
                    p.pattern_name.clone(),
                )
            });
        } else {
            let mut blocked = 0usize;
            pairs.retain(|p| {
                if !keep_pair(p, self.keep) {
                    return false;
                }
                let atoms: BTreeSet<usize> = p.plan_site_atoms().into_iter().collect();
                let oxy = oxygen_site_from_parts(p.site, &[p.site], atoms);
                if undoes_oxygen_edit(&p.effect, &oxy, self.o_added, self.o_removed) {
                    blocked += 1;
                    return false;
                }
                true
            });
            self.counters.blocked_circular_oxygen += blocked;
        }
        // Leaf-first path: this set, then ancestors still on the stack (root last).
        let mut rule_path = vec![set.name.clone()];
        for frame in self.pair_stack.iter().rev() {
            rule_path.push(frame.set.name.clone());
        }
        let pending: Vec<PendingPair<'a>> = pairs
            .into_iter()
            .map(|pair| {
                let pair = pair.with_rule_path(rule_path.clone());
                PendingPair { pair, set }
            })
            .collect();
        self.pair_pending = pending.into_iter();
        Ok(())
    }

    fn advance_pairs(&mut self) -> Result<(), ForestError> {
        loop {
            if !self.pair_pending.as_slice().is_empty() {
                return Ok(());
            }
            let action = {
                let Some(frame) = self.pair_stack.last_mut() else {
                    return Ok(());
                };
                let members = frame.set.members();
                if frame.member_i < members.len() {
                    let member = &members[frame.member_i];
                    frame.member_i += 1;
                    match member {
                        crate::ruleset::RuleMember::Set(child) => {
                            // `'a` outlives the frame; child lives in the ruleset tree.
                            let child: &'a RuleSet = child;
                            Some(PairAction::Push(child))
                        }
                        crate::ruleset::RuleMember::Pattern(_) => Some(PairAction::Continue),
                    }
                } else {
                    let set = frame.set;
                    Some(PairAction::LoadLeaf(set))
                }
            };
            match action {
                None => return Ok(()),
                Some(PairAction::Continue) => {}
                Some(PairAction::Push(child)) => {
                    self.pair_stack.push(PairFrame {
                        set: child,
                        member_i: 0,
                    });
                }
                Some(PairAction::LoadLeaf(set)) => {
                    self.pair_stack.pop();
                    self.load_leaf_pairs(set)?;
                }
            }
        }
    }
}

enum PairAction<'a> {
    Continue,
    Push(&'a RuleSet),
    LoadLeaf(&'a RuleSet),
}

impl<K> Iterator for Expand<'_, K>
where
    K: Fn(&Candidate) -> bool,
{
    type Item = Result<ForestEmission, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            while let Some(candidate) = self.deferred.next() {
                self.counters.mol_edits += 1;
                match self.emit_candidate(&candidate) {
                    Ok(Some(emission)) => return Some(Ok(emission)),
                    Ok(None) => {}
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            if let Some(pending) = self.pair_pending.next() {
                match self.emit_pair(&pending) {
                    Ok(Some(emission)) => return Some(Ok(emission)),
                    Ok(None) => continue,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            match self.advance_pairs() {
                Ok(()) => {
                    if self.pair_pending.as_slice().is_empty() && self.pair_stack.is_empty() {
                        self.done = true;
                        return None;
                    }
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

fn keep_pair<K>(pair: &crate::pair_edit::PairCandidate, keep: &K) -> bool
where
    K: Fn(&Candidate) -> bool,
{
    let mut stand_in = Candidate {
        site: pair.site,
        orbit: vec![pair.site],
        pattern: pair.left.clone(),
        rule_path: Vec::new(),
        mapped: Default::default(),
        parent: crate::candidate::ParentRef::Context,
    };
    stand_in.pattern.effect = pair.effect.clone();
    stand_in.pattern.name = pair.pattern_name.clone();
    keep(&stand_in)
}

/// Python-style filter closures (optional). Prefer [`find_path_with`] + reading
/// [`Candidate`] fields when the predicate does not need the leaf `RuleSet`.
///
/// Pull iterator: metabolize emissions are consumed one at a time (no full
/// emission list). Callers that want a list use `.collect()`.
pub fn find_path_with_filters<'a, 'b, R, S>(
    reactant: &str,
    target: &str,
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    config: FindPathConfig,
    filter_rules: R,
    filter_sites: S,
) -> Result<FindPathFilters<'a, 'b, R, S>, ForestError>
where
    R: Fn(&crate::Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&crate::Molecule, usize, &SiteInfo) -> bool,
{
    let start = ForestMol::parse(reactant)?;
    let target_csmi = canon_of(target)?;
    let target_ha = ForestMol::parse(&target_csmi)?.heavy_atom_count();

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    let mut seen = HashSet::new();
    remember_seen(&mut seen, &start);
    let ancestors = root_ancestors(&start);
    heap.push(HeapItem {
        mode: config.heap_score,
        search_bias: 0,
        site_progress: 0,
        cost_gain: 0,
        match_score: 0,
        match_priority: 0,
        diversity_key: None,
        seq,
        walk: Walk {
            mol: start,
            steps: Vec::new(),
            plan: Vec::new(),
            maybe: Vec::new(),
            opens: Vec::new(),
            o_added: Vec::new(),
            o_removed: Vec::new(),
            ancestors,
            parent_cost: None,
            diff: None,
        },
    });
    seq += 1;

    Ok(FindPathFilters {
        ruleset,
        counters,
        filter_rules,
        filter_sites,
        config,
        target_csmi,
        target_ha,
        heap,
        seq,
        seen,
        yielded: Vec::new(),
        done: false,
    })
}

/// Filter-closure search: pulls [`RuleSet::metabolize`] one emission at a time.
pub struct FindPathFilters<'a, 'b, R, S> {
    ruleset: &'a RuleSet,
    counters: &'b mut PathCounters,
    filter_rules: R,
    filter_sites: S,
    config: FindPathConfig,
    target_csmi: String,
    target_ha: usize,
    heap: BinaryHeap<HeapItem>,
    seq: usize,
    seen: HashSet<String>,
    yielded: Vec<PathOutcome>,
    done: bool,
}

impl<R, S> Iterator for FindPathFilters<'_, '_, R, S>
where
    R: Fn(&crate::Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&crate::Molecule, usize, &SiteInfo) -> bool,
{
    type Item = Result<PathOutcome, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let FindPathConfig {
            max_paths,
            max_nodes,
            drop_skeleton_twins,
            ..
        } = self.config;

        while let Some(item) = self.heap.pop() {
            if self.yielded.len() >= max_paths || self.counters.nodes >= max_nodes {
                break;
            }
            self.counters.nodes += 1;
            let walk = item.walk;
            let here = walk.mol.csmi();
            if here.as_ref() == self.target_csmi.as_str() {
                let plan = as_deps(walk.plan).with_maybe(Maybe::new(walk.maybe));
                record_yield_plan_signals(self.counters, &self.yielded, &plan, drop_skeleton_twins);
                if plan_already_yielded(&self.yielded, &plan, drop_skeleton_twins) {
                    continue;
                }
                let outcome = PathOutcome {
                    steps: walk.steps,
                    plan,
                    smiles: here.as_ref().to_string(),
                };
                self.yielded.push(outcome.clone());
                return Some(Ok(outcome));
            }

            let mol = walk.mol.mol();
            self.counters.expansions += 1;
            let mut hits_from_here = 0usize;

            // Pull metabolize one emission at a time — no full list.
            let emissions =
                self.ruleset
                    .metabolize(mol, &self.filter_rules, &self.filter_sites, true);
            for emission in emissions {
                let emission = match emission {
                    Ok(e) => e,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                self.counters.mol_edits += 1;
                // Filter-path emissions are CSMI strings (no tag continuity). Prefer
                // find_path_with + Candidate materialize for tagged walks.
                let products: Result<Vec<_>, _> = emission
                    .products
                    .iter()
                    .map(|s| ForestMol::parse(s))
                    .collect();
                let products = match products {
                    Ok(p) => p,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                let keeps = keep_fragments(
                    &walk.mol,
                    &products,
                    &self.target_csmi,
                    self.target_ha,
                    None,
                    None,
                    &mut LiftMcsCounters {
                        fallback: &mut self.counters.mcs_lift_fallback,
                        rematch: &mut self.counters.mcs_lift_rematch,
                    },
                );
                for (kept, sides, _) in keeps {
                    let kept_csmi = kept.csmi().as_ref().to_string();
                    let child_ha = kept.heavy_atom_count();
                    let target_hit = kept_csmi == self.target_csmi;
                    if !closer(
                        walk.mol.heavy_atom_count(),
                        child_ha,
                        self.target_ha,
                        target_hit,
                    ) {
                        continue;
                    }
                    if repeats_ancestor(&walk.ancestors, &kept) {
                        continue;
                    }
                    if already_seen(&self.seen, &kept) && !target_hit {
                        continue;
                    }
                    if kept.stable_csmi_key().is_none() {
                        self.counters.unstable_csmi_key += 1;
                    }
                    remember_seen(&mut self.seen, &kept);

                    let site_atoms: BTreeSet<usize> = emission.site_atoms.iter().copied().collect();
                    let (child_maybe, child_opens) = accumulate_maybe(
                        &walk.maybe,
                        &walk.opens,
                        &site_atoms,
                        emission.cleaves,
                        products.len(),
                        &sides,
                    );
                    let mut steps = walk.steps.clone();
                    steps.push(PathStep {
                        rule_path: emission.rule_path.clone(),
                        pattern_name: emission.pattern_name.clone(),
                        site: emission.site,
                        site_orbit: emission.site_orbit.clone(),
                        product: kept_csmi.clone(),
                        sides,
                    });
                    let mut plan = walk.plan.clone();
                    plan.extend(emission.plan.iter().cloned());
                    let ancestors = with_child_ancestor(&walk.ancestors, &kept);
                    self.heap.push(HeapItem {
                        mode: self.config.heap_score,
                        search_bias: emission.search_bias,
                        site_progress: 0,
                        // No atom_diff on this path — HA closer already gated.
                        cost_gain: 0,
                        match_score: 0,
                        match_priority: 0,
                        diversity_key: None,
                        seq: self.seq,
                        walk: Walk {
                            mol: kept,
                            steps,
                            plan,
                            maybe: child_maybe,
                            opens: child_opens,
                            o_added: walk.o_added.clone(),
                            o_removed: walk.o_removed.clone(),
                            ancestors,
                            parent_cost: None,
                            diff: None,
                        },
                    });
                    self.seq += 1;
                    if target_hit {
                        hits_from_here += 1;
                        if self.yielded.len() + hits_from_here >= max_paths {
                            break;
                        }
                    }
                }
                if self.yielded.len() + hits_from_here >= max_paths {
                    break;
                }
            }
        }

        self.done = true;
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydroxylation::hydroxylation;
    use crate::mol::canon_of;
    use crate::ruleset::o_dealkylation;

    #[test]
    fn heap_prefers_most_recently_queued_among_peers() {
        // Among equal soft scores, higher seq pops first (Ord tiebreak).
        let older = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 1,
            walk: Walk {
                mol: ForestMol::parse("CC").unwrap(),
                steps: vec![],
                plan: vec![],
                maybe: vec![],
                opens: vec![],
                o_added: vec![],
                o_removed: vec![],
                ancestors: HashSet::new(),
                parent_cost: None,
                diff: None,
            },
        };
        let newer = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 2,
            walk: older.walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(older);
        heap.push(newer);
        assert_eq!(heap.pop().unwrap().seq, 2);
        assert_eq!(heap.pop().unwrap().seq, 1);
    }

    #[test]
    fn heap_pops_best_ord_value_only() {
        // Plain BinaryHeap::pop — among equal soft scores, higher seq wins.
        let walk = Walk {
            mol: ForestMol::parse("CC").unwrap(),
            steps: vec![],
            plan: vec![],
            maybe: vec![],
            opens: vec![],
            o_added: vec![],
            o_removed: vec![],
            ancestors: HashSet::new(),
            parent_cost: None,
            diff: None,
        };
        let mk = |seq| HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,
            match_priority: 0,
            diversity_key: None,
            seq,
            walk: walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(mk(1));
        heap.push(mk(2));
        heap.push(mk(3));
        assert_eq!(heap.pop().unwrap().seq, 3);
        assert_eq!(heap.pop().unwrap().seq, 2);
        assert_eq!(heap.pop().unwrap().seq, 1);
        assert!(heap.pop().is_none());
    }

    #[test]
    fn heap_prefers_higher_search_bias_over_seq() {
        // Good scores override DFS: demoted bias loses even if enqueued later.
        let demoted = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: -1,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 99,
            walk: Walk {
                mol: ForestMol::parse("CC").unwrap(),
                steps: vec![],
                plan: vec![],
                maybe: vec![],
                opens: vec![],
                o_added: vec![],
                o_removed: vec![],
                ancestors: HashSet::new(),
                parent_cost: None,
                diff: None,
            },
        };
        let preferred = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 1,
            walk: demoted.walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(demoted);
        heap.push(preferred);
        assert_eq!(heap.pop().unwrap().search_bias, 0);
        assert_eq!(heap.pop().unwrap().search_bias, -1);
    }

    #[test]
    fn heap_prefers_higher_site_progress_over_seq() {
        // Site H-progress overrides DFS among equal search_bias.
        let low = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 99,
            walk: Walk {
                mol: ForestMol::parse("CC").unwrap(),
                steps: vec![],
                plan: vec![],
                maybe: vec![],
                opens: vec![],
                o_added: vec![],
                o_removed: vec![],
                ancestors: HashSet::new(),
                parent_cost: None,
                diff: None,
            },
        };
        let high = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 2,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 1,
            walk: low.walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(low);
        heap.push(high);
        assert_eq!(heap.pop().unwrap().site_progress, 2);
        assert_eq!(heap.pop().unwrap().site_progress, 0);
    }

    #[test]
    fn heap_lack_of_improvement_counters_dfs() {
        // Non-positive cost_gain loses to an older positive peer despite LIFO.
        let flat_newer = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 0,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 99,
            walk: Walk {
                mol: ForestMol::parse("CC").unwrap(),
                steps: vec![],
                plan: vec![],
                maybe: vec![],
                opens: vec![],
                o_added: vec![],
                o_removed: vec![],
                ancestors: HashSet::new(),
                parent_cost: None,
                diff: None,
            },
        };
        let gain_older = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 1,
            walk: flat_newer.walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(flat_newer);
        heap.push(gain_older);
        let first = heap.pop().unwrap();
        assert_eq!(first.cost_gain, 1);
        assert_eq!(first.seq, 1);
        assert_eq!(heap.pop().unwrap().cost_gain, 0);
    }

    #[test]
    fn heap_prefers_larger_cost_gain_over_seq() {
        let small_newer = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 1,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 99,
            walk: Walk {
                mol: ForestMol::parse("CC").unwrap(),
                steps: vec![],
                plan: vec![],
                maybe: vec![],
                opens: vec![],
                o_added: vec![],
                o_removed: vec![],
                ancestors: HashSet::new(),
                parent_cost: None,
                diff: None,
            },
        };
        let big_older = HeapItem {
            mode: HeapScoreMode::SoftStack,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 3,
            match_score: 0,

            match_priority: 0,
            diversity_key: None,
            seq: 1,
            walk: small_newer.walk.clone(),
        };
        let mut heap = BinaryHeap::new();
        heap.push(small_newer);
        heap.push(big_older);
        assert_eq!(heap.pop().unwrap().cost_gain, 3);
        assert_eq!(heap.pop().unwrap().cost_gain, 1);
    }

    #[test]
    fn diversity_penalty_is_neg_n_fixed() {
        assert_eq!(diversity_penalty(0), 0);
        assert_eq!(diversity_penalty(1), -(MATCH_LOG_SCALE as i64));
        assert_eq!(diversity_penalty(4), -4 * (MATCH_LOG_SCALE as i64));
        // Stronger than −ln(n+1) for n ≥ 1.
        assert!(diversity_penalty(1) < neg_log1p_score(1));
        assert!(diversity_penalty(3) < neg_log1p_score(3));
    }

    #[test]
    fn match_heap_key_is_priority_then_score() {
        // Ord: (score + diversity, score). Higher priority wins; equal priority
        // prefers higher base score.
        let walk = Walk {
            mol: ForestMol::parse("CC").unwrap(),
            steps: vec![],
            plan: vec![],
            maybe: vec![],
            opens: vec![],
            o_added: vec![],
            o_removed: vec![],
            ancestors: HashSet::new(),
            parent_cost: None,
            diff: None,
        };
        let mode = HeapScoreMode::match_log_neg_pc();
        let high_pri = HeapItem {
            mode,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 0,
            match_score: 10,
            match_priority: 10 + diversity_penalty(0),
            diversity_key: None,
            seq: 1,
            walk: walk.clone(),
        };
        let low_pri_better_score = HeapItem {
            mode,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 0,
            match_score: 100,
            match_priority: 100 + diversity_penalty(5),
            diversity_key: None,
            seq: 2,
            walk: walk.clone(),
        };
        // diversity_penalty(5) = −5·SCALE → priority << base score.
        assert!(high_pri > low_pri_better_score);

        let same_pri_low_score = HeapItem {
            mode,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 0,
            match_score: 5,
            match_priority: 10,
            diversity_key: None,
            seq: 3,
            walk: walk.clone(),
        };
        let same_pri_high_score = HeapItem {
            mode,
            search_bias: 0,
            site_progress: 0,
            cost_gain: 0,
            match_score: 10,
            match_priority: 10,
            diversity_key: None,
            seq: 1,
            walk,
        };
        assert!(same_pri_high_score > same_pri_low_score);
    }

    #[test]
    fn diversity_key_uses_site_and_added_tags() {
        let parent = ForestMol::parse("C").unwrap();
        let site = parent.tag_of(0).unwrap();
        // Index-stable apply: parent C → product atom 0; O at 1 is born.
        let co = parse_mol("CO").unwrap();
        assert_eq!(co.atom_count(), 2);
        let product = parent.from_apply(co, &[Some(0)]);
        let site_atoms: BTreeSet<usize> = [0].into_iter().collect();
        let key = diversity_key_for(&parent, &product, "hydroxylation", &site_atoms);
        assert_eq!(key.rule, "hydroxylation");
        assert_eq!(key.site_tags, BTreeSet::from([site]));
        assert_eq!(key.added_tags.len(), 1);
        let added = *key.added_tags.iter().next().unwrap();
        assert!(parent.index_of(added).is_none());
        assert_eq!(product.tag_of(1), Some(added));
    }

    #[test]
    fn diversity_opt_in_finds_simple_path() {
        let rules = o_dealkylation();
        let mut with = PathCounters::default();
        let mut without = PathCounters::default();
        let a = find_path_with(
            "COc1ccccc1",
            "Oc1ccccc1",
            &rules,
            &mut with,
            FindPathConfig {
                max_paths: 1,
                diversity: true,
                ..FindPathConfig::default()
            },
            accept_all_candidates,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        let b = find_path_with(
            "COc1ccccc1",
            "Oc1ccccc1",
            &rules,
            &mut without,
            FindPathConfig {
                max_paths: 1,
                ..FindPathConfig::default()
            },
            accept_all_candidates,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert_eq!(a.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(a[0].smiles, b[0].smiles);
        assert_eq!(without.diversity_repush, 0);
        assert!(!FindPathConfig::default().diversity);
    }

    #[test]
    fn hop_cost_gain_is_parent_minus_child() {
        assert_eq!(hop_cost_gain(Some(5), Some(0)), 5);
        assert_eq!(hop_cost_gain(Some(5), Some(4)), 1);
        assert_eq!(hop_cost_gain(Some(5), Some(5)), 0);
        assert_eq!(hop_cost_gain(Some(5), Some(6)), -1);
        assert_eq!(hop_cost_gain(Some(5), None), 0);
    }

    #[test]
    fn neg_log1p_score_decreases_with_cost() {
        assert_eq!(neg_log1p_score(0), 0);
        assert!(neg_log1p_score(0) > neg_log1p_score(1));
        assert!(neg_log1p_score(1) > neg_log1p_score(10));
        // log_close = ln(SCALE) + neg_log1p (fixed-point rounding may differ by 1).
        let ln_scale = ln_fixed(MATCH_CLOSE_SCALE);
        assert!((log_close_term(3) - (ln_scale + neg_log1p_score(3))).abs() <= 1);
    }

    #[test]
    fn match_add_prefers_joint_improvement_and_closeness() {
        // Better atom+formula improvement at same closeness → higher score.
        let flat = hop_match_add_score(4, 4, Some(10), Some(10));
        let better = hop_match_add_score(4, 2, Some(10), Some(5));
        assert!(better > flat, "better={better} flat={flat}");
        // Closer child beats farther at equal hop gain.
        let close = hop_match_add_score(6, 4, Some(12), Some(10));
        let far = hop_match_add_score(8, 6, Some(14), Some(12));
        assert!(close > far, "close={close} far={far}");
        // Exact match (cost/formula 0) outranks a partial improvement via closeness.
        assert!(hop_match_add_score(9, 0, Some(9), Some(0)) > better);
    }

    #[test]
    fn match_add_ranks_identically_to_product() {
        // Add = ln(Product factors); pairwise order must match Product.
        let both_add = MatchScoreSpec::add_both();
        let both_prod = MatchScoreSpec::product_both();
        let cases = [
            (4, 4, Some(10usize), Some(10usize)),
            (4, 2, Some(10), Some(5)),
            (6, 4, Some(12), Some(10)),
            (8, 6, Some(14), Some(12)),
            (9, 0, Some(9), Some(0)),
            (8, 4, Some(20), Some(10)),
            (0, 0, Some(0), Some(0)),
            (5, 5, None, None),
            (3, 1, Some(7), None),
        ];
        for (i, a) in cases.iter().enumerate() {
            for (j, b) in cases.iter().enumerate() {
                if i >= j {
                    continue;
                }
                let pa = hop_match_score(both_prod, a.0, a.1, a.2, a.3);
                let pb = hop_match_score(both_prod, b.0, b.1, b.2, b.3);
                let aa = hop_match_score(both_add, a.0, a.1, a.2, a.3);
                let ab = hop_match_score(both_add, b.0, b.1, b.2, b.3);
                assert_eq!(
                    pa.cmp(&pb),
                    aa.cmp(&ab),
                    "product {pa:?} vs {pb:?} but add {aa:?} vs {ab:?} at {a:?} / {b:?}"
                );
            }
        }
    }

    #[test]
    fn heap_score_from_label_roundtrips() {
        assert_eq!(
            HeapScoreMode::from_label("log-neg-pc"),
            Some(HeapScoreMode::match_log_neg_pc())
        );
        assert_eq!(
            HeapScoreMode::from_label("soft"),
            Some(HeapScoreMode::SoftStack)
        );
        assert_eq!(
            HeapScoreMode::from_label("add-both"),
            Some(HeapScoreMode::match_add())
        );
        assert!(HeapScoreMode::from_label("nope").is_none());
    }

    #[test]
    fn match_combine_and_metric_axes() {
        let both_add = MatchScoreSpec::add_both();
        let close_both = MatchScoreSpec {
            combine: MatchCombine::Close,
            metric: MatchMetric::Both,
        };
        let improve_both = MatchScoreSpec {
            combine: MatchCombine::Improve,
            metric: MatchMetric::Both,
        };
        let add_atom = MatchScoreSpec {
            combine: MatchCombine::Add,
            metric: MatchMetric::Atom,
        };
        let add_formula = MatchScoreSpec {
            combine: MatchCombine::Add,
            metric: MatchMetric::Formula,
        };
        // Same residual, different hop gain: improve ranks the gain; close ties.
        let gain = hop_match_score(improve_both, 8, 4, Some(20), Some(10));
        let flat = hop_match_score(improve_both, 4, 4, Some(10), Some(10));
        assert!(gain > flat);
        let close_a = hop_match_score(close_both, 8, 4, Some(20), Some(10));
        let close_b = hop_match_score(close_both, 4, 4, Some(10), Some(10));
        assert_eq!(close_a, close_b, "close ignores hop gain at equal residual");
        // Atom-only ignores formula; formula-only ignores atom.
        let atom = hop_match_score(add_atom, 0, 9, Some(10), Some(2));
        let atom_same_a = hop_match_score(add_atom, 9, 0, Some(10), Some(2));
        assert_eq!(atom, atom_same_a);
        let formula = hop_match_score(add_formula, 10, 2, Some(0), Some(9));
        let formula_same_a = hop_match_score(add_formula, 10, 2, Some(9), Some(0));
        assert_eq!(formula, formula_same_a);
        // Default add-both beats flat when both axes improve.
        let joint = hop_match_score(both_add, 4, 2, Some(10), Some(5));
        let joint_flat = hop_match_score(both_add, 4, 4, Some(10), Some(10));
        assert!(joint > joint_flat);
        assert_eq!(MatchScoreSpec::matrix().len(), 9);
        // Add = close + improve (log factors).
        assert_eq!(
            hop_match_score(both_add, 8, 4, Some(20), Some(10)),
            hop_match_score(close_both, 8, 4, Some(20), Some(10))
                + hop_match_score(improve_both, 8, 4, Some(20), Some(10))
        );
    }

    #[test]
    fn ethane_to_ethanol_is_one_hydroxylation() {
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let outcome = &hits[0];
        assert_eq!(outcome.smiles, canon_of("CCO").unwrap());
        assert_eq!(outcome.steps.len(), 1);
        assert_eq!(outcome.steps[0].leaf_rule(), Some("Hydroxylation"));
        assert_eq!(outcome.steps[0].namespace(), vec!["Hydroxylation"]);
        assert!(outcome.steps[0].sides.is_empty());
        assert_eq!(outcome.plan.len(), 1);
        assert_eq!(outcome.plan[0].rule, "Hydroxylation");
        assert_eq!(counters.mol_edits, 1);
    }

    #[test]
    fn anisole_to_phenol_cleaves_and_records_side() {
        let mut counters = PathCounters::default();
        let hits = find_path("COc1ccccc1", "Oc1ccccc1", &o_dealkylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let outcome = &hits[0];
        assert_eq!(outcome.smiles, canon_of("Oc1ccccc1").unwrap());
        assert_eq!(outcome.steps.len(), 1);
        assert_eq!(outcome.steps[0].leaf_rule(), Some("Dealkylation"));
        assert_eq!(outcome.steps[0].pattern_name, "O-Me");
        assert!(
            !outcome.steps[0].sides.is_empty(),
            "cleavage should leave a side fragment"
        );
        assert!(!outcome.maybe().is_empty(), "maybe bags on the plan");
        assert!(
            outcome.allows(None, Some(outcome.maybe().entries[0].side.as_str())),
            "allows discarded side"
        );
        assert_eq!(outcome.plan[0].rule, "Dealkylation");
    }

    #[test]
    fn quinone_formation_cleavage_records_side_and_maybe() {
        // Anisole → methane via QF dealkylate: bifurcation (n_products≥2),
        // not ring-open. Leaf is QuinoneFormation; quinone sits in maybe/sides.
        use crate::rules::quinone_formation;

        let mut counters = PathCounters::default();
        let hits = find_path("COc1ccccc1", "C", &quinone_formation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let outcome = &hits[0];
        assert_eq!(outcome.smiles, canon_of("C").unwrap());
        assert_eq!(outcome.steps[0].leaf_rule(), Some("QuinoneFormation"));
        assert!(
            !outcome.steps[0].sides.is_empty(),
            "QF cleavage must record the quinone side: {:?}",
            outcome.steps[0].sides
        );
        assert!(
            !outcome.maybe().is_empty(),
            "maybe bags on the plan for discarded quinone"
        );
        let want_para = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        let want_ortho = canon_of("O=C1C=CC=CC1=O").unwrap();
        let is_quinone = |s: &str| {
            let c = canon_of(s).unwrap();
            c == want_para || c == want_ortho
        };
        let in_sides = outcome.steps[0].sides.iter().any(|s| is_quinone(s));
        let in_maybe = outcome.maybe().sides().iter().any(|s| is_quinone(s));
        assert!(
            in_sides || in_maybe,
            "quinone in sides/maybe; sides={:?} maybe={:?}",
            outcome.steps[0].sides,
            outcome.maybe().sides()
        );
    }

    #[test]
    fn composed_ruleset_path_keeps_leaf_and_outer_namespace() {
        let set = RuleSet::compose(Some("Forest".into()), [hydroxylation(), o_dealkylation()]);
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &set, &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        let step = &hits[0].steps[0];
        assert_eq!(step.namespace(), vec!["Hydroxylation", "Forest"]);
        assert_eq!(step.leaf_rule(), Some("Hydroxylation"));
    }

    #[test]
    fn nested_ruleset_path_appends_each_set() {
        let inner = RuleSet::compose(Some("Inner".into()), [hydroxylation()]);
        let outer = RuleSet::compose(Some("Outer".into()), [inner]);
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &outer, &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        assert_eq!(
            hits[0].steps[0].namespace(),
            vec!["Hydroxylation", "Inner", "Outer"]
        );
    }

    #[test]
    fn already_at_target_yields_empty_plan() {
        let mut counters = PathCounters::default();
        let hits = find_path("CCO", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert_eq!(hits.len(), 1);
        assert!(hits[0].steps.is_empty());
        assert!(hits[0].plan.is_empty());
        assert!(hits[0].maybe().is_empty());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        assert_eq!(counters.nodes, 1);
        assert_eq!(counters.mol_edits, 0);
    }

    #[test]
    fn keep_predicate_skips_materialize() {
        let mut counters = PathCounters::default();
        let refuse_h2 = |c: &Candidate| c.pattern.name != "h2";
        let hits = find_path_with(
            "CC",
            "CCO",
            &hydroxylation(),
            &mut counters,
            FindPathConfig::default(),
            refuse_h2,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        // Ethane only matches h2; refusing it yields no path and no edits.
        assert!(hits.is_empty());
        assert_eq!(counters.mol_edits, 0);
    }

    #[test]
    fn default_ruleset_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("CC", "CCO", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Hydroxylation"));
        assert!(!hits[0].plan.is_empty());
    }

    #[test]
    fn default_ruleset_anisole_to_phenol() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("COc1ccccc1", "Oc1ccccc1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("Oc1ccccc1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Dealkylation"));
        assert!(!hits[0].steps[0].sides.is_empty());
    }

    #[test]
    fn phase_one_ethene_to_epoxide() {
        let mut counters = PathCounters::default();
        let hits = find_path("C=C", "C1CO1", &crate::rules::phase_one(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("C1CO1").unwrap());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Epoxidation"));
    }

    #[test]
    fn default_ruleset_hydroquinone_to_quinone() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("O=C1C=CC(=O)C=C1").unwrap());
        let leaf = hits[0].steps[0].leaf_rule();
        assert!(
            leaf == Some("Dehydrogenation") || leaf == Some("QuinoneFormation"),
            "leaf={leaf:?} steps={:?}",
            hits[0].steps
        );
        assert!(
            !hits[0].plan.is_empty(),
            "expected elementary plan, got {:?}",
            hits[0].plan
        );
    }

    #[test]
    fn benzene_to_quinone_plan_precedes_both_oh_before_dh() {
        // Python test_find_path_uses_atom_diff_filters: precedes (0,2),(1,2).
        let mut counters = PathCounters::default();
        let hits = find_path_default("c1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let names: Vec<_> = hits[0].plan.iter().map(|s| s.rule.as_str()).collect();
        assert_eq!(
            names.iter().filter(|&&n| n == "Hydroxylation").count(),
            2,
            "plan={names:?}"
        );
        assert_eq!(
            names.iter().filter(|&&n| n == "Dehydrogenation").count(),
            1,
            "plan={names:?}"
        );
        assert!(!names.contains(&"QuinoneFormation"));
        let hydroxyl: Vec<_> = names
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n == "Hydroxylation").then_some(i))
            .collect();
        let dh = names.iter().position(|&n| n == "Dehydrogenation").unwrap();
        let edges: HashSet<_> = hits[0].plan.precedes().iter().copied().collect();
        assert!(!edges.contains(&(hydroxyl[0], hydroxyl[1])));
        assert!(!edges.contains(&(hydroxyl[1], hydroxyl[0])));
        assert!(edges.contains(&(hydroxyl[0], dh)));
        assert!(edges.contains(&(hydroxyl[1], dh)));
    }

    #[test]
    fn phenol_to_quinone_plan_oh_precedes_dh() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("Oc1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let names: Vec<_> = hits[0].plan.iter().map(|s| s.rule.as_str()).collect();
        assert!(!names.contains(&"QuinoneFormation"));
        assert_eq!(names.iter().filter(|&&n| n == "Hydroxylation").count(), 1);
        assert_eq!(names.iter().filter(|&&n| n == "Dehydrogenation").count(), 1);
        let oh = names.iter().position(|&n| n == "Hydroxylation").unwrap();
        let dh = names.iter().position(|&n| n == "Dehydrogenation").unwrap();
        let edges: HashSet<_> = hits[0].plan.precedes().iter().copied().collect();
        assert!(edges.contains(&(oh, dh)));
    }

    #[test]
    fn path_hit_plan_replays_to_hit() {
        let mut counters = PathCounters::default();
        let hits = find_path_default("c1ccccc1", "O=C1C=CC(=O)C=C1", &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert!(
            hits[0].plan.reaches("c1ccccc1", &hits[0].smiles).unwrap(),
            "plan={:?} smiles={}",
            hits[0].plan,
            hits[0].smiles
        );
    }

    #[test]
    fn ethane_plan_replays_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        assert!(hits[0].plan.reaches("CC", "CCO").unwrap());
    }

    #[test]
    fn atom_diff_gates_ethane_to_ethanol() {
        let mut counters = PathCounters::default();
        let hits = find_path_diff("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("CCO").unwrap());
        // Diff gating should not inflate edits beyond the one helpful site.
        assert_eq!(counters.mol_edits, 1);
    }

    #[test]
    fn atom_diff_default_ruleset_hydroquinone() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "Oc1ccc(O)cc1",
            "O=C1C=CC(=O)C=C1",
            &default_ruleset(),
            &mut counters,
            FindPathConfig {
                use_atom_diff: true,
                ..FindPathConfig::default()
            },
            accept_all_candidates,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert_eq!(hits[0].smiles, canon_of("O=C1C=CC(=O)C=C1").unwrap());
    }

    #[test]
    fn walk_preserves_tags_through_hydroxylation() {
        let parent = ForestMol::parse("CC").unwrap();
        let t0 = parent.tag_of(0).expect("stamped");
        let t1 = parent.tag_of(1).expect("stamped");
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert!(!hits.is_empty());
        // Product mol is not on the outcome; check adopt via a fresh emit.
        let set = hydroxylation();
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let pieces = cands[0].materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        // Surviving carbon tags are stable; born O is remapped via labels.
        // adopt_product keeps the shared tag_gen with the parent.
        assert!(child.shares_tag_gen(&parent));
        assert_eq!(child.tag_of(child.index_of(t0).unwrap()), Some(t0));
        assert_eq!(child.tag_of(child.index_of(t1).unwrap()), Some(t1));
        assert_eq!(child.mol().atom_count(), 3);
    }

    #[test]
    fn cleavage_side_and_maybe_allows_on_deps() {
        use crate::canonical_plan::{CleavageSide, Maybe};
        use std::collections::BTreeSet;

        let side = CleavageSide::new([1, 2], canon_of("CCO").unwrap(), []);
        let empty = Maybe::default();
        assert!(empty.is_empty());

        let filled = Maybe::new([side]);
        assert!(!filled.is_empty());
        let demethyl: BTreeSet<_> = [0usize, 1].into_iter().collect();
        let formation: BTreeSet<_> = [1usize, 2].into_iter().collect();
        assert!(filled.allows(Some(&demethyl), None));
        assert!(!filled.allows(Some(&formation), None));
        assert!(filled.allows(None, Some("CCO")));
        let far: BTreeSet<_> = [9usize, 10].into_iter().collect();
        assert!(!filled.allows(Some(&far), None));

        let plan = Deps::bind([crate::canonical_plan::Step::new(
            "Dealkylation",
            [crate::canonical_plan::PlanAtom::index(1)],
        )])
        .with_maybe(filled);
        assert!(plan.allows(Some(&demethyl), None));
        assert!(!plan.allows(Some(&formation), None));
    }

    #[test]
    fn hydrolysis_emits_maybe_bag_on_plan() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            &crate::rules::hydrolysis(),
            &mut counters,
            FindPathConfig {
                max_paths: 3,
                max_nodes: 40,
                use_atom_diff: false,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let with_bag: Vec<_> = hits.iter().filter(|h| !h.maybe().is_empty()).collect();
        assert!(!with_bag.is_empty(), "expected CleavageSide on hydrolysis");
        let outcome = with_bag[0];
        assert_eq!(
            outcome
                .plan
                .iter()
                .map(|s| s.rule.as_str())
                .collect::<Vec<_>>(),
            vec!["Hydrolysis"]
        );
        let want_side = canon_of("CC(C)(C)O").unwrap();
        assert!(
            outcome.maybe().entries.iter().any(|e| e.side == want_side),
            "sides={:?}",
            outcome.maybe().sides()
        );
        assert!(outcome.allows(None, Some(&want_side)));
    }

    #[test]
    fn dealkylation_emits_maybe_and_multipath() {
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::dealkylation(),
            &mut counters,
            FindPathConfig {
                max_paths: 4,
                max_nodes: 80,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        assert!(
            hits.iter().any(|h| !h.maybe().is_empty()),
            "expected maybe bag"
        );
        let parent_ha = ForestMol::parse("CN(C)Cc1ccccc1")
            .unwrap()
            .heavy_atom_count();
        for h in &hits {
            for e in &h.maybe().entries {
                let side_ha = ForestMol::parse(&e.side).unwrap().heavy_atom_count();
                assert!(side_ha < parent_ha, "side={} ha={side_ha}", e.side);
            }
        }
    }

    #[test]
    fn multipath_emits_up_to_max_paths() {
        let mut counters = PathCounters::default();
        let one = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 1,
                max_nodes: 100,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert_eq!(one.len(), 1);

        let mut counters = PathCounters::default();
        let many = find_path_with(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 4,
                max_nodes: 200,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(
            many.len() > 1 && many.len() <= 4,
            "max_paths=1 → {} hits; max_paths=4 → {} (want >1)",
            one.len(),
            many.len()
        );
        assert!(many.iter().all(|h| h.smiles == one[0].smiles));
        // HEURISTICS: a reordering of the same Deps is not a new path.
        // Stronger observational check: multipath hits share no total orders.
        for (i, a) in many.iter().enumerate() {
            for b in many.iter().skip(i + 1) {
                assert!(
                    !a.plan.same_linearizations(&b.plan),
                    "hit{i} same_linearizations as later hit"
                );
                assert_eq!(
                    a.plan.linearization_overlap(&b.plan),
                    0,
                    "hit{i} linearization_overlap={} with later hit",
                    a.plan.linearization_overlap(&b.plan)
                );
            }
        }
    }

    #[test]
    fn plan_drop_counters_signal_only_no_contained_prune() {
        // Duplicate → count + drop. Contained-extension → count only, keep yield.
        use crate::canonical_plan::{PlanAtom, Step};
        let short = as_deps([Step::new("Dealkylation", [PlanAtom::index(0)])]);
        let longer = as_deps([
            Step::new("Dealkylation", [PlanAtom::index(0)]),
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
        ]);
        let twin = short.clone();
        let found = vec![PathOutcome {
            steps: vec![],
            plan: short,
            smiles: "C".into(),
        }];
        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &twin, true);
        assert_eq!(counters.dropped_duplicate_plan, 1);
        assert_eq!(counters.dropped_exact_plan, 1);
        assert_eq!(counters.dropped_skeleton_twin, 0);
        assert_eq!(counters.signal_contained_plan, 0);
        assert!(plan_already_yielded(&found, &twin, true));

        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &longer, true);
        assert_eq!(counters.dropped_duplicate_plan, 0);
        assert_eq!(counters.dropped_exact_plan, 0);
        assert_eq!(counters.dropped_skeleton_twin, 0);
        assert_eq!(counters.signal_contained_plan, 1);
        assert!(!plan_already_yielded(&found, &longer, true));
        assert_eq!(counters.plan_drops(), 1);
    }

    #[test]
    fn skeleton_twin_drop_counts_separately_from_exact() {
        use crate::canonical_plan::{CleavageSide, Maybe, PlanAtom, Step};
        let maybe = Maybe::new([
            CleavageSide::new([0], "OC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new([9], "OC", std::iter::empty::<BTreeSet<usize>>()),
        ]);
        let first = as_deps([
            Step::new("Dealkylation", [PlanAtom::index(0)]),
            Step::new("Dealkylation", [PlanAtom::index(9)]),
        ])
        .with_maybe(maybe.clone());
        // Remapped indices, same rule multiset + Maybe sides → skeleton only.
        let twin = as_deps([
            Step::new("Dealkylation", [PlanAtom::index(12)]),
            Step::new("Dealkylation", [PlanAtom::index(0)]),
        ])
        .with_maybe(maybe);
        assert!(!first.same_linearizations(&twin));
        assert!(first.same_rule_maybe_skeleton(&twin));

        let found = vec![PathOutcome {
            steps: vec![],
            plan: first,
            smiles: "C".into(),
        }];
        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &twin, true);
        assert_eq!(counters.dropped_duplicate_plan, 1);
        assert_eq!(counters.dropped_exact_plan, 0);
        assert_eq!(counters.dropped_skeleton_twin, 1);

        let mut counters = PathCounters::default();
        record_yield_plan_signals(&mut counters, &found, &twin, false);
        assert_eq!(counters.dropped_duplicate_plan, 0);
        assert_eq!(counters.dropped_skeleton_twin, 0);
        assert!(!plan_already_yielded(&found, &twin, false));
    }

    #[test]
    fn tba_cleavage_maybe_allows_overlapping_site() {
        // Terbinafine → TBA: formation Dealkylation leaves naphthyl-amine bag.
        const TERB: &str = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12";
        const TBA: &str = r"C(#C/C=C/C=O)C(C)(C)C";
        let mut counters = PathCounters::default();
        let hits = find_path_with(
            TERB,
            TBA,
            &crate::rules::phase_one(),
            &mut counters,
            FindPathConfig {
                max_paths: 3,
                max_nodes: 80,
                ..FindPathConfig::default()
            },
            |_| true,
        )
        .unwrap()
        .collect_all()
        .unwrap();
        assert!(!hits.is_empty(), "billed={}", counters.billed());
        let direct: Vec<_> = hits
            .iter()
            .filter(|h| h.plan.len() == 1 && h.plan[0].rule == "Dealkylation")
            .collect();
        assert!(!direct.is_empty(), "expected Dealkylation → TBA");
        let outcome = direct[0];
        assert!(!outcome.maybe().is_empty());
        assert!(
            outcome
                .maybe()
                .entries
                .iter()
                .any(|e| e.side.contains("cccc")),
            "sides={:?}",
            outcome.maybe().sides()
        );
        let formation = outcome.maybe().entries[0].site.clone();
        assert!(!outcome.allows(Some(&formation), None));
        // Overlapping different site (share one atom) passes Maybe.
        let mut demethyl = BTreeSet::new();
        demethyl.insert(*formation.iter().next().unwrap());
        demethyl.insert(formation.iter().next().unwrap() + 1000);
        assert!(outcome.allows(Some(&demethyl), None));
        assert!(outcome.allows(None, Some(outcome.maybe().entries[0].side.as_str())));
    }

    /// Ethane → ethanol → ethene with **only** Hydroxylation + Dehydration.
    /// Beta-elim site includes the adjacent carbon (`site_map` [1,3]), so the
    /// site set differs from hydroxylation — allowed. Alcohol dehydration at
    /// the lone OH carbon would be equal-site circular (see
    /// `oxygen_site_equal_uses_orbits_and_allows_site_diff`).
    #[test]
    fn hydroxylation_then_dehydration_same_site_yields_ethene() {
        use crate::rules::dehydration;
        use crate::ruleset::{accept_all_rules, accept_all_sites};

        let set = RuleSet::compose(
            Some("HydrateDehydrate".into()),
            [hydroxylation(), dehydration()],
        );
        let ethane = ForestMol::parse("CC").unwrap();
        let ethane_c0 = ethane.tag_of(0).expect("tag");

        let oh_cands: Vec<_> = hydroxylation()
            .candidates(ethane.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!oh_cands.is_empty());
        let oh = &oh_cands[0];
        let oh_pieces = oh.materialize_mols(ethane.mol()).unwrap();
        assert_eq!(oh_pieces.len(), 1);
        let ethanol = ethane.adopt_product(oh_pieces[0].clone());
        assert_eq!(ethanol.csmi().as_ref(), canon_of("CCO").unwrap());
        let oh_carbon = ethanol.index_of(ethane_c0).expect("tagged carbon");
        assert!(
            oh.site == oh_carbon || oh.orbit.contains(&oh_carbon),
            "OH site {} orbit {:?} should name tagged carbon {}",
            oh.site,
            oh.orbit,
            oh_carbon
        );

        let ethene_csmi = canon_of("C=C").unwrap();
        let mut saw_ethene = false;
        let mut shared_site = false;
        for emission in set.metabolize(ethanol.mol(), accept_all_rules, accept_all_sites, true) {
            let emission = emission.unwrap();
            if emission.leaf_rule() != Some("Dehydration") {
                continue;
            }
            let site_atoms: BTreeSet<usize> = emission.site_atoms.iter().copied().collect();
            if site_atoms.contains(&oh_carbon) || emission.site == oh_carbon {
                shared_site = true;
            }
            for p in &emission.products {
                if p.as_str() == ethene_csmi.as_str() {
                    saw_ethene = true;
                }
            }
        }
        assert!(
            shared_site,
            "dehydration should act on the hydroxylated carbon {oh_carbon}"
        );
        assert!(
            saw_ethene,
            "hydroxylation then dehydration at that carbon yields ethene (productive)"
        );
    }

    #[test]
    fn oxygen_site_equal_uses_orbits_and_allows_site_diff() {
        // Singletons: top orbit membership.
        let oh = oxygen_site_from_parts(0, &[0, 1], BTreeSet::from([0]));
        let alcohol_same = oxygen_site_from_parts(0, &[0], BTreeSet::from([0]));
        let alcohol_peer = oxygen_site_from_parts(1, &[1], BTreeSet::from([1]));
        assert!(oxygen_sites_equal(&oh, &alcohol_same));
        assert!(oxygen_sites_equal(&oh, &alcohol_peer)); // peer in OH orbit
        // Multi-atom beta-elim differs by adjacent C → not circular.
        let beta = oxygen_site_from_parts(0, &[0], BTreeSet::from([0, 1]));
        assert!(!oxygen_sites_equal(&oh, &beta));
        let rem = crate::pattern::Effect {
            removes: Some("OH".into()),
            ..Default::default()
        };
        assert!(undoes_oxygen_edit(&rem, &alcohol_same, &[oh.clone()], &[]));
        assert!(!undoes_oxygen_edit(&rem, &beta, &[oh], &[]));
    }
}
