//! Pull iterators for [`RuleSet`] candidates and metabolize emissions.
//!
//! Real generators: each [`Iterator::next`] discovers or materializes the next
//! item. Unique-edit may buffer SMARTS hits **per pattern** (orbits need the
//! full match set for that pattern); the ruleset walk does not buffer every
//! pattern's candidates or every emission before yielding. Pair discovery may
//! buffer SMARTS hits for one leaf's endpoints before yielding pairs.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::ForestError;
use crate::candidate::{Candidate, ParentRef};
use crate::kekule::kekule_forms;
use crate::mol::{Molecule, stable_csmi_key_of};
use crate::pair_edit::{PairCandidate, pair_candidates as discover_pairs};
use crate::pattern::{Edit, Emission, PatternInfo, SiteInfo};
use crate::ruleset::{RuleMember, RuleSet, site_bond_is_exclusive_double};
use crate::unique_edit::{unique_sites, unique_sites_on_forms};

/// Discover site–pattern triples one at a time (no edit).
pub struct Candidates<'a> {
    set: &'a RuleSet,
    mol: &'a Molecule,
    /// When nested under a parent set, append that name on each yield (leaf-first).
    parent_link: Option<Option<String>>,
    member_i: usize,
    pending: std::vec::IntoIter<Candidate>,
    child: Option<Box<Candidates<'a>>>,
    done: bool,
}

impl<'a> Candidates<'a> {
    pub(crate) fn new(set: &'a RuleSet, mol: &'a Molecule) -> Self {
        Self {
            set,
            mol,
            parent_link: None,
            member_i: 0,
            pending: Vec::new().into_iter(),
            child: None,
            done: false,
        }
    }

    fn nested(set: &'a RuleSet, mol: &'a Molecule, parent_name: Option<String>) -> Self {
        Self {
            set,
            mol,
            parent_link: Some(parent_name),
            member_i: 0,
            pending: Vec::new().into_iter(),
            child: None,
            done: false,
        }
    }

    fn finish_candidate(&self, mut c: Candidate) -> Candidate {
        if let Some(parent_name) = &self.parent_link {
            c.rule_path.push(parent_name.clone());
        }
        c
    }

    fn load_pattern(&mut self, pattern: &PatternInfo) -> Result<(), ForestError> {
        let batch = pattern_candidate_batch(self.set, self.mol, pattern)?;
        self.pending = batch.into_iter();
        Ok(())
    }
}

impl Iterator for Candidates<'_> {
    type Item = Result<Candidate, ForestError>;

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
            if self.member_i >= members.len() {
                self.done = true;
                return None;
            }
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
                        self.mol,
                        self.set.name.clone(),
                    )));
                }
            }
        }
    }
}

/// Per-pattern unique-edit buffer only (orbits need all hits for that SMARTS).
pub(crate) fn pattern_candidate_batch(
    set: &RuleSet,
    mol: &Molecule,
    pattern: &PatternInfo,
) -> Result<Vec<Candidate>, ForestError> {
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
            out.push(Candidate {
                site,
                orbit: hit.orbit,
                pattern: pattern.resolve_for_match(mol, &hit.mapped),
                rule_path: vec![set.name.clone()],
                mapped: hit.mapped,
                parent: ParentRef::Form(Box::new(forms[form_i].clone())),
            });
        }
    } else {
        for hit in unique_sites(mol, &pattern.smarts, pattern.site_kind, &pattern.site_map)? {
            let Some(&site) = hit.mapped.get(&pattern.primary_map()) else {
                continue;
            };
            out.push(Candidate {
                site,
                orbit: hit.orbit,
                pattern: pattern.resolve_for_match(mol, &hit.mapped),
                rule_path: vec![set.name.clone()],
                mapped: hit.mapped,
                parent: ParentRef::Context,
            });
        }
    }
    Ok(out)
}

/// Materialize emissions one at a time (filter → edit → yield).
pub struct Metabolize<'a, R, S> {
    set: &'a RuleSet,
    mol: &'a Molecule,
    filter_rules: Rc<R>,
    filter_sites: Rc<S>,
    unique_csmi: bool,
    member_i: usize,
    pattern_pending: std::vec::IntoIter<Candidate>,
    child: Option<Box<Metabolize<'a, R, S>>>,
    pairs: std::vec::IntoIter<PairCandidate>,
    pairs_loaded: bool,
    pair_patterns: Vec<PatternInfo>,
    seen_csmi: BTreeMap<BTreeSet<String>, String>,
    seen_leaf: BTreeSet<(String, BTreeSet<String>)>,
    /// `Some(parent_name)` when this walk is nested under a parent set.
    parent_link: Option<Option<String>>,
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
        mol: &'a Molecule,
        filter_rules: R,
        filter_sites: S,
        unique_csmi: bool,
    ) -> Self {
        Self {
            set,
            mol,
            filter_rules: Rc::new(filter_rules),
            filter_sites: Rc::new(filter_sites),
            unique_csmi,
            member_i: 0,
            pattern_pending: Vec::new().into_iter(),
            child: None,
            pairs: Vec::new().into_iter(),
            pairs_loaded: false,
            pair_patterns: Vec::new(),
            seen_csmi: BTreeMap::new(),
            seen_leaf: BTreeSet::new(),
            parent_link: None,
            done: false,
        }
    }

    fn nested(
        set: &'a RuleSet,
        mol: &'a Molecule,
        filter_rules: Rc<R>,
        filter_sites: Rc<S>,
        parent_name: Option<String>,
    ) -> Self {
        Self {
            set,
            mol,
            filter_rules,
            filter_sites,
            unique_csmi: false,
            member_i: 0,
            pattern_pending: Vec::new().into_iter(),
            child: None,
            pairs: Vec::new().into_iter(),
            pairs_loaded: false,
            pair_patterns: Vec::new(),
            seen_csmi: BTreeMap::new(),
            seen_leaf: BTreeSet::new(),
            parent_link: Some(parent_name),
            done: false,
        }
    }

    fn take_emission(
        &mut self,
        mut emission: Emission,
        from_nested_child: bool,
    ) -> Option<Emission> {
        if let Some(parent_name) = &self.parent_link {
            emission.rule_path.push(parent_name.clone());
        }
        // Fail-closed: only dedup when every product has a stable Chematic key.
        let emission_key = emission_stable_product_key(&emission.products);
        if from_nested_child {
            // Parent unique_csmi is the cross-child CSMI layer.
            if self.unique_csmi {
                if let Some(key) = &emission_key {
                    let leaf_name = emission.leaf_rule().unwrap_or("").to_string();
                    if let Some(kept) = self.seen_csmi.get(key) {
                        if kept != &leaf_name {
                            return None;
                        }
                    } else {
                        self.seen_csmi.insert(key.clone(), leaf_name);
                    }
                }
            }
            return Some(emission);
        }
        // Leaf pattern / pair on this set: within-leaf unique_csmi yield.
        if self.unique_csmi {
            if let Some(key) = emission_key {
                let leaf_key = (emission.pattern_name.clone(), key);
                if !self.seen_leaf.insert(leaf_key) {
                    return None;
                }
            }
        }
        Some(emission)
    }
}

/// Product multiset for yield dedup — `None` if any fragment lacks a stable key.
fn emission_stable_product_key(products: &[String]) -> Option<BTreeSet<String>> {
    let mut key = BTreeSet::new();
    for p in products {
        key.insert(stable_csmi_key_of(p)?);
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
            while let Some(c) = self.pattern_pending.next() {
                let info = SiteInfo {
                    site: c.site,
                    orbit: c.orbit.clone(),
                    pattern: c.pattern.clone(),
                    shell_forecast: None,
                };
                if !(self.filter_sites)(self.mol, c.site, &info) {
                    continue;
                }
                match c.emit(self.mol) {
                    Ok(Some(mut emission)) => {
                        if self.set.has_plan_hook() {
                            emission.plan =
                                self.set
                                    .canonical_plan(self.mol, &emission.site_atoms, None);
                        }
                        if let Some(e) = self.take_emission(emission, false) {
                            return Some(Ok(e));
                        }
                    }
                    Ok(None) => {}
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            if let Some(child) = self.child.as_mut() {
                match child.next() {
                    Some(Ok(emission)) => {
                        if let Some(e) = self.take_emission(emission, true) {
                            return Some(Ok(e));
                        }
                        continue;
                    }
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
                        if !(self.filter_rules)(self.mol, self.set, pattern) {
                            continue;
                        }
                        match pattern_candidate_batch(self.set, self.mol, pattern) {
                            Ok(batch) => self.pattern_pending = batch.into_iter(),
                            Err(e) => {
                                self.done = true;
                                return Some(Err(e));
                            }
                        }
                    }
                    RuleMember::Set(child) => {
                        self.child = Some(Box::new(Metabolize::nested(
                            child,
                            self.mol,
                            Rc::clone(&self.filter_rules),
                            Rc::clone(&self.filter_sites),
                            self.set.name.clone(),
                        )));
                    }
                }
                continue;
            }

            if !self.pairs_loaded {
                self.pairs_loaded = true;
                let endpoints: Vec<PatternInfo> = self
                    .set
                    .leaf_pair_endpoints()
                    .into_iter()
                    .filter(|p| (self.filter_rules)(self.mol, self.set, p))
                    .collect();
                self.pair_patterns = endpoints.clone();
                if !endpoints.is_empty() {
                    match discover_pairs(self.mol, &endpoints) {
                        Ok(pairs) => self.pairs = pairs.into_iter(),
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
                        }
                    }
                }
                continue;
            }

            while let Some(pair) = self.pairs.next() {
                let pattern = self.pair_patterns.first().cloned().unwrap_or_else(|| {
                    PatternInfo::new(
                        "pair",
                        "[#6:1]",
                        Edit::PairEndpoint("x".into()),
                        Default::default(),
                    )
                });
                let info = SiteInfo {
                    site: pair.site,
                    orbit: vec![pair.site],
                    pattern,
                    shell_forecast: None,
                };
                if !(self.filter_sites)(self.mol, pair.site, &info) {
                    continue;
                }
                match pair.emit(self.mol) {
                    Ok(Some(emission)) => {
                        let site_atoms = pair.plan_site_atoms();
                        let ends = [&pair.left.effect, &pair.right.effect];
                        let plan = self.set.canonical_plan(self.mol, &site_atoms, Some(&ends));
                        let emission = Emission {
                            site: emission.site,
                            site_orbit: vec![emission.site],
                            site_atoms: site_atoms.clone(),
                            cleaves: pair.effect.cleaves,
                            pattern_name: emission.pattern_name,
                            search_bias: pair.left.search_bias.min(pair.right.search_bias),
                            rule_path: vec![self.set.name.clone()],
                            products: emission.products,
                            plan,
                        };
                        if let Some(e) = self.take_emission(emission, false) {
                            return Some(Ok(e));
                        }
                    }
                    Ok(None) => {}
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }

            self.done = true;
            return None;
        }
    }
}

/// Streaming pair-candidate discovery across this set and nested children.
pub struct PairCandidates<'a> {
    sets: Vec<&'a RuleSet>,
    set_i: usize,
    mol: &'a Molecule,
    pending: std::vec::IntoIter<PairCandidate>,
    done: bool,
}

impl<'a> PairCandidates<'a> {
    pub(crate) fn new(set: &'a RuleSet, mol: &'a Molecule) -> Self {
        let mut sets = Vec::new();
        fn walk<'b>(s: &'b RuleSet, out: &mut Vec<&'b RuleSet>) {
            out.push(s);
            for m in s.members() {
                if let RuleMember::Set(child) = m {
                    walk(child, out);
                }
            }
        }
        walk(set, &mut sets);
        Self {
            sets,
            set_i: 0,
            mol,
            pending: Vec::new().into_iter(),
            done: false,
        }
    }
}

impl Iterator for PairCandidates<'_> {
    type Item = Result<PairCandidate, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            if let Some(p) = self.pending.next() {
                return Some(Ok(p));
            }
            if self.set_i >= self.sets.len() {
                self.done = true;
                return None;
            }
            let set = self.sets[self.set_i];
            self.set_i += 1;
            let endpoints = set.leaf_pair_endpoints();
            if endpoints.is_empty() {
                continue;
            }
            match discover_pairs(self.mol, &endpoints) {
                Ok(pairs) => self.pending = pairs.into_iter(),
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

/// Stream ResonancePair emissions (materialize one pair at a time).
pub struct PairEmissions<'a> {
    stack: Vec<PairEmissionFrame<'a>>,
    mol: &'a Molecule,
    pending: std::vec::IntoIter<PendingPairEmission<'a>>,
    done: bool,
}

struct PairEmissionFrame<'a> {
    set: &'a RuleSet,
    member_i: usize,
}

struct PendingPairEmission<'a> {
    pair: PairCandidate,
    set: &'a RuleSet,
    rule_path: Vec<Option<String>>,
}

enum PairEmissionAction<'a> {
    Continue,
    Push(&'a RuleSet),
    LoadLeaf(&'a RuleSet),
}

impl<'a> PairEmissions<'a> {
    pub(crate) fn new(set: &'a RuleSet, mol: &'a Molecule) -> Self {
        Self {
            stack: vec![PairEmissionFrame { set, member_i: 0 }],
            mol,
            pending: Vec::new().into_iter(),
            done: false,
        }
    }

    fn load_leaf(&mut self, set: &'a RuleSet) -> Result<(), ForestError> {
        let pairs = set.pair_candidates_leaf(self.mol)?;
        let mut rule_path = vec![set.name.clone()];
        for frame in self.stack.iter().rev() {
            rule_path.push(frame.set.name.clone());
        }
        let pending: Vec<_> = pairs
            .into_iter()
            .map(|pair| PendingPairEmission {
                pair,
                set,
                rule_path: rule_path.clone(),
            })
            .collect();
        self.pending = pending.into_iter();
        Ok(())
    }

    fn advance(&mut self) -> Result<(), ForestError> {
        loop {
            if !self.pending.as_slice().is_empty() {
                return Ok(());
            }
            let action = {
                let Some(frame) = self.stack.last_mut() else {
                    return Ok(());
                };
                let members = frame.set.members();
                if frame.member_i < members.len() {
                    let member = &members[frame.member_i];
                    frame.member_i += 1;
                    match member {
                        RuleMember::Set(child) => PairEmissionAction::Push(child),
                        RuleMember::Pattern(_) => PairEmissionAction::Continue,
                    }
                } else {
                    PairEmissionAction::LoadLeaf(frame.set)
                }
            };
            match action {
                PairEmissionAction::Continue => {}
                PairEmissionAction::Push(child) => {
                    self.stack.push(PairEmissionFrame {
                        set: child,
                        member_i: 0,
                    });
                }
                PairEmissionAction::LoadLeaf(set) => {
                    self.stack.pop();
                    self.load_leaf(set)?;
                }
            }
        }
    }
}

impl Iterator for PairEmissions<'_> {
    type Item = Result<Emission, ForestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            if let Some(pending) = self.pending.next() {
                match pending.pair.emit(self.mol) {
                    Ok(Some(emission)) => {
                        let site_atoms = pending.pair.plan_site_atoms();
                        let ends = [&pending.pair.left.effect, &pending.pair.right.effect];
                        let plan = pending
                            .set
                            .canonical_plan(self.mol, &site_atoms, Some(&ends));
                        return Some(Ok(Emission {
                            site: emission.site,
                            site_orbit: vec![emission.site],
                            site_atoms: site_atoms.clone(),
                            cleaves: pending.pair.effect.cleaves,
                            pattern_name: emission.pattern_name,
                            search_bias: pending
                                .pair
                                .left
                                .search_bias
                                .min(pending.pair.right.search_bias),
                            rule_path: pending.rule_path,
                            products: emission.products,
                            plan,
                        }));
                    }
                    Ok(None) => continue,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                }
            }
            match self.advance() {
                Ok(()) => {
                    if self.pending.as_slice().is_empty() && self.stack.is_empty() {
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
