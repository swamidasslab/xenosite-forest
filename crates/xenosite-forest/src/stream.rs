//! Pull iterators for [`RuleSet`] candidates and metabolize emissions.
//!
//! [`Candidates`] discovers every SOM (atom/bond + composed pairs). Unique-edit
//! may buffer SMARTS hits **per pattern**; pair ends compose via
//! [`crate::pair_edit::compose_pair_sites`]. [`Metabolize`] is a thin filter →
//! [`DeferredSite::apply`] walk over the same candidate stream. Products are
//! tagged [`crate::ForestMol`] adopted from the discovery parent.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::ForestError;
use crate::ForestMol;
use crate::candidate::{DeferredSite, ParentRef};
use crate::kekule::kekule_forms;
use crate::mol::Molecule;
use crate::pair_edit::compose_candidates_from_endpoints;
use crate::pattern::{Edit, Emission, PatternInfo, SiteInfo};
use crate::ruleset::{RuleMember, RuleSet, site_bond_is_exclusive_double};
use crate::unique_edit::{unique_sites, unique_sites_on_forms};

/// Discover site–pattern triples one at a time (no edit).
///
/// Yields every SOM on this set and nested children: unique-edit atom/bond
/// sites, then ResonancePair sites for each leaf's endpoint patterns.
pub struct Candidates<'a> {
    set: &'a RuleSet,
    mol: Rc<ForestMol>,
    /// When nested under a parent set, append that name on each yield (leaf-first).
    parent_link: Option<Option<String>>,
    member_i: usize,
    pending: std::vec::IntoIter<DeferredSite>,
    child: Option<Box<Candidates<'a>>>,
    pairs_loaded: bool,
    done: bool,
}

impl<'a> Candidates<'a> {
    pub(crate) fn new(set: &'a RuleSet, mol: &'a ForestMol) -> Self {
        Self {
            set,
            mol: Rc::new(mol.copy_mol()),
            parent_link: None,
            member_i: 0,
            pending: Vec::new().into_iter(),
            child: None,
            pairs_loaded: false,
            done: false,
        }
    }

    fn nested(set: &'a RuleSet, mol: Rc<ForestMol>, parent_name: Option<String>) -> Self {
        Self {
            set,
            mol,
            parent_link: Some(parent_name),
            member_i: 0,
            pending: Vec::new().into_iter(),
            child: None,
            pairs_loaded: false,
            done: false,
        }
    }

    fn finish_candidate(&self, mut c: DeferredSite) -> DeferredSite {
        if let Some(parent_name) = &self.parent_link {
            c.rule_path.push(parent_name.clone());
        }
        c
    }

    fn load_pattern(&mut self, pattern: &PatternInfo) -> Result<(), ForestError> {
        let batch = pattern_candidate_batch(self.set, &self.mol, pattern)?;
        self.pending = batch.into_iter();
        Ok(())
    }

    fn load_pairs(&mut self) -> Result<(), ForestError> {
        let endpoints = self.set.leaf_pair_endpoints();
        if endpoints.is_empty() {
            return Ok(());
        }
        let mut pairs = compose_candidates_from_endpoints(Rc::clone(&self.mol), &endpoints)?;
        for p in &mut pairs {
            if p.rule_path.is_empty() {
                p.rule_path.push(self.set.name.clone());
            }
        }
        self.pending = pairs.into_iter();
        Ok(())
    }
}

impl Iterator for Candidates<'_> {
    type Item = Result<DeferredSite, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            if let Some(c) = self.pending.next() {
                return Some(Ok(self.finish_candidate(c)));
            }
            if let Some(child) = self.child.as_mut() {
                match child.next() {
                    Some(Ok(c)) => return Some(Ok(self.finish_candidate(c))),
                    Some(Err(e)) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                    None => {
                        self.child = None;
                        continue;
                    }
                }
            }
            let members = self.set.members();
            if self.member_i < members.len() {
                let member = &members[self.member_i];
                self.member_i += 1;
                match member {
                    RuleMember::Pattern(pattern) => {
                        if matches!(pattern.edit, Edit::PairEndpoint(_)) {
                            continue;
                        }
                        if let Err(e) = self.load_pattern(pattern) {
                            self.done = true;
                            return Some(Err(e));
                        }
                    }
                    RuleMember::Set(child) => {
                        self.child = Some(Box::new(Candidates::nested(
                            child,
                            Rc::clone(&self.mol),
                            self.set.name.clone(),
                        )));
                    }
                }
                continue;
            }
            if !self.pairs_loaded {
                self.pairs_loaded = true;
                if let Err(e) = self.load_pairs() {
                    self.done = true;
                    return Some(Err(e));
                }
                continue;
            }
            self.done = true;
            return None;
        }
    }
}

/// Per-pattern unique-edit buffer only (orbits need all hits for that SMARTS).
pub(crate) fn pattern_candidate_batch(
    set: &RuleSet,
    forest: &Rc<ForestMol>,
    pattern: &PatternInfo,
) -> Result<Vec<DeferredSite>, ForestError> {
    let mol = forest.mol();
    let mut out = Vec::new();
    let use_forms = site_bond_is_exclusive_double(&pattern.smarts, &pattern.site_map)
        && mol.atoms().any(|(_, atom)| atom.aromatic);
    if use_forms {
        let forms = kekule_forms(mol)?;
        let hits = unique_sites_on_forms(
            mol,
            &forms,
            &pattern.smarts,
            pattern.site_kind,
            &pattern.site_map,
        )?;
        for (hit, form_i) in hits {
            let Some(&site) = hit.mapped.get(&pattern.primary_map()) else {
                continue;
            };
            // Resolve on aromatic context `mol`, not the Kekulé form.
            out.push(DeferredSite::atom(
                Rc::clone(forest),
                site,
                hit.orbit,
                pattern.resolve_for_match(mol, &hit.mapped),
                vec![set.name.clone()],
                hit.mapped,
                ParentRef::Form(Box::new(forms[form_i].clone())),
            ));
        }
    } else {
        for hit in unique_sites(mol, &pattern.smarts, pattern.site_kind, &pattern.site_map)? {
            let Some(&site) = hit.mapped.get(&pattern.primary_map()) else {
                continue;
            };
            out.push(DeferredSite::atom(
                Rc::clone(forest),
                site,
                hit.orbit,
                pattern.resolve_for_match(mol, &hit.mapped),
                vec![set.name.clone()],
                hit.mapped,
                ParentRef::Context,
            ));
        }
    }
    Ok(out)
}

/// Materialize emissions: iterate [`Candidates`], filter, [`DeferredSite::apply`].
pub struct Metabolize<'a, R, S> {
    sites: Candidates<'a>,
    filter_rules: R,
    filter_sites: S,
    unique_csmi: bool,
    seen_csmi: BTreeMap<BTreeSet<String>, String>,
    seen_leaf: BTreeSet<(String, BTreeSet<String>)>,
    done: bool,
}

/// [`Metabolize`] with accept-all filters (no closures).
pub type OpenMetabolize<'a> = Metabolize<
    'a,
    fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
    fn(&Molecule, usize, &SiteInfo) -> bool,
>;

impl<'a, R, S> Metabolize<'a, R, S>
where
    R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&Molecule, usize, &SiteInfo) -> bool,
{
    pub(crate) fn new(
        set: &'a RuleSet,
        mol: &'a ForestMol,
        filter_rules: R,
        filter_sites: S,
        unique_csmi: bool,
    ) -> Self {
        Self {
            sites: Candidates::new(set, mol),
            filter_rules,
            filter_sites,
            unique_csmi,
            seen_csmi: BTreeMap::new(),
            seen_leaf: BTreeSet::new(),
            done: false,
        }
    }

    fn leaf_for_filter(&self, site: &DeferredSite) -> RuleSet {
        site.leaf_rule()
            .and_then(crate::rules::leaf_rule)
            .unwrap_or_else(|| RuleSet::new(site.leaf_rule().map(str::to_string), []))
    }

    fn keep_site(&self, site: &DeferredSite) -> bool {
        let leaf = self.leaf_for_filter(site);
        let mol = site.mol();
        if site.is_pair() {
            if !(self.filter_rules)(mol, &leaf, &site.left) {
                return false;
            }
            if !(self.filter_rules)(mol, &leaf, &site.right) {
                return false;
            }
        } else if !(self.filter_rules)(mol, &leaf, &site.pattern) {
            return false;
        }
        (self.filter_sites)(mol, site.info.site, &site.info)
    }

    fn take_emission(&mut self, emission: Emission) -> Option<Emission> {
        if !self.unique_csmi {
            return Some(emission);
        }
        // Fail-closed: only dedup when every product has a stable Chematic key.
        let Some(emission_key) = emission_stable_product_key(&emission.products) else {
            return Some(emission);
        };
        // Within-leaf dedup (same pattern + product multiset).
        let leaf_key = (emission.pattern_name.clone(), emission_key.clone());
        if !self.seen_leaf.insert(leaf_key) {
            return None;
        }
        // Cross-leaf: first leaf to produce this product multiset wins.
        let leaf_name = emission.leaf_rule().unwrap_or("").to_string();
        if let Some(kept) = self.seen_csmi.get(&emission_key) {
            if kept != &leaf_name {
                return None;
            }
        } else {
            self.seen_csmi.insert(emission_key, leaf_name);
        }
        Some(emission)
    }
}

/// Product multiset for yield dedup — `None` if any fragment lacks a stable key.
fn emission_stable_product_key(products: &[ForestMol]) -> Option<BTreeSet<String>> {
    let mut key = BTreeSet::new();
    for p in products {
        key.insert(p.stable_csmi_key()?.as_ref().to_string());
    }
    Some(key)
}

impl<'a, R, S> Iterator for Metabolize<'a, R, S>
where
    R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
    S: Fn(&Molecule, usize, &SiteInfo) -> bool,
{
    type Item = Result<Emission, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            match self.sites.next() {
                None => {
                    self.done = true;
                    return None;
                }
                Some(Err(e)) => {
                    self.done = true;
                    return Some(Err(e));
                }
                Some(Ok(site)) => {
                    if !self.keep_site(&site) {
                        continue;
                    }
                    match site.apply() {
                        Ok(None) => continue,
                        Ok(Some(emission)) => {
                            if let Some(e) = self.take_emission(emission) {
                                return Some(Ok(e));
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
    }
}
