//! [`BoundPattern`]: a leaf [`RuleSet`] narrowed to one [`PatternInfo`].
//!
//! Behaves like the owning rule for metabolize / member access, with a built-in
//! filter that keeps this pattern and (for two-site candidates) requires at
//! least one site atom to be among this pattern's defined site map.

use crate::ForestMol;
use crate::mol::Molecule;
use crate::pattern::{PatternInfo, SiteInfo};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};

/// Pattern handle bound to its owning leaf rule (inseparable).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundPattern {
    rule: RuleSet,
    pattern: PatternInfo,
}

impl BoundPattern {
    /// Construct from a leaf rule and a pattern that must be among its patterns.
    pub fn new(rule: RuleSet, pattern: PatternInfo) -> Self {
        Self { rule, pattern }
    }

    pub fn rule(&self) -> &RuleSet {
        &self.rule
    }

    pub fn pattern(&self) -> &PatternInfo {
        &self.pattern
    }

    pub fn name(&self) -> &str {
        &self.pattern.name
    }

    pub fn rule_name(&self) -> Option<&str> {
        self.rule.name.as_deref()
    }

    /// `xf:Rule/pattern` CURIE.
    pub fn curie(&self) -> String {
        match self.rule_name() {
            Some(rule) => format!("xf:{rule}/{}", self.pattern.name),
            None => format!("xf:?/{}", self.pattern.name),
        }
    }

    /// Absolute IRI under `https://w3id.org/xenosite/forest/`.
    pub fn iri(&self) -> String {
        crate::mapping::expand_iri(&self.curie())
    }

    /// Length of the pattern view (always 1).
    pub fn len(&self) -> usize {
        1
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    pub fn get(&self, index: usize) -> Option<BoundPattern> {
        (index == 0).then(|| self.clone())
    }

    pub fn get_str(&self, name: &str) -> Option<BoundPattern> {
        (name == self.pattern.name).then(|| self.clone())
    }

    pub fn contains_name(&self, name: &str) -> bool {
        name == self.pattern.name
    }

    /// One-pattern view for compose / callers that need a bare [`RuleSet`].
    pub fn as_ruleset(&self) -> RuleSet {
        RuleSet::new(self.rule.name.clone(), [self.pattern.clone()])
    }

    /// `filter_rules` that accepts only this pattern name, AND-composed with `extra`.
    pub fn compose_rule_filter<'a, R>(
        &'a self,
        extra: R,
    ) -> impl Fn(&Molecule, &RuleSet, &PatternInfo) -> bool + 'a
    where
        R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool + 'a,
    {
        let pattern_name = self.pattern.name.as_str();
        move |mol, rule, pattern| pattern.name == pattern_name && extra(mol, rule, pattern)
    }

    /// Site filter AND-composed with `extra`.
    ///
    /// Multi-map patterns (two sites): do **not** require every site atom to be
    /// owned by this pattern — pair partners may come from elsewhere on the leaf.
    /// The pattern-name rule filter already restricts which pattern emits; site
    /// filtering only forwards to `extra`.
    pub fn compose_site_filter<'a, S>(
        &'a self,
        extra: S,
    ) -> impl Fn(&Molecule, usize, &SiteInfo) -> bool + 'a
    where
        S: Fn(&Molecule, usize, &SiteInfo) -> bool + 'a,
    {
        move |mol, site, info| extra(mol, site, info)
    }

    /// Metabolize through the owning leaf with the BoundPattern filters applied.
    #[allow(clippy::type_complexity)]
    pub fn metabolize<'a, R, S>(
        &'a self,
        mol: &'a ForestMol,
        filter_rules: R,
        filter_sites: S,
        unique_csmi: bool,
    ) -> crate::stream::Metabolize<
        'a,
        impl Fn(&Molecule, &RuleSet, &PatternInfo) -> bool + 'a,
        impl Fn(&Molecule, usize, &SiteInfo) -> bool + 'a,
    >
    where
        R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool + 'a,
        S: Fn(&Molecule, usize, &SiteInfo) -> bool + 'a,
    {
        let rules = self.compose_rule_filter(filter_rules);
        let sites = self.compose_site_filter(filter_sites);
        self.rule.metabolize(mol, rules, sites, unique_csmi)
    }

    /// Metabolize with accept-all caller filters (BoundPattern filter only).
    #[allow(clippy::type_complexity)]
    pub fn metabolize_default<'a>(
        &'a self,
        mol: &'a ForestMol,
        unique_csmi: bool,
    ) -> crate::stream::Metabolize<
        'a,
        impl Fn(&Molecule, &RuleSet, &PatternInfo) -> bool + 'a,
        impl Fn(&Molecule, usize, &SiteInfo) -> bool + 'a,
    > {
        self.metabolize(mol, accept_all_rules, accept_all_sites, unique_csmi)
    }
}
