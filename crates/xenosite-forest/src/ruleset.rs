//! A set of [`PatternInfo`] records and nested sets, as namespaces.
//!
//! Nested sets stay nested: [`RuleSet::compose`] does not flatten. Each set
//! appends itself onto the emission [`crate::pattern::Emission::rule_path`]
//! after the child that emitted (leaf first, outer last), matching Python
//! `info["rule"]` / addition chain order.
//!
//! Primary walk: [`RuleSet::candidates`] yields site–pattern–[`ParentRef`]
//! triples without applying edits. A search reads [`PatternInfo`] / [`Effect`]
//! to filter, then [`Candidate::materialize`] only for survivors. Filter
//! closures are optional convenience on [`RuleSet::metabolize`], not required.

use std::collections::{BTreeMap, BTreeSet};

use chematic::core::{Atom, BondOrder, Element};
use chematic::smarts::{BondPrimitive, BondQuery, parse_smarts};

use crate::ForestError;
use crate::candidate::{Candidate, ParentRef};
use crate::canonical_plan::{CanonicalPlanFn, Step, steps_for_leaf};
use crate::kekule::kekule_forms;
use crate::mol::{Molecule, atom_idx, canon_smiles};
use crate::pair_edit::pair_candidates;
use crate::pattern::{Edit, Emission, PatternInfo, SiteInfo};
use crate::smirks::apply_smirks_at;
use crate::unique_edit::{unique_sites, unique_sites_on_forms};
use crate::valence::accept_product;

/// True when the SMARTS bond between the first two `site_map` atoms is an
/// exclusive double (`=`), not `=,:`. Epoxidation matches that on Kekulé forms.
fn site_bond_is_exclusive_double(smarts: &str, site_map: &[u16]) -> bool {
    if site_map.len() < 2 {
        return false;
    }
    let Ok(query) = parse_smarts(smarts) else {
        return false;
    };
    let mut idxs = BTreeMap::new();
    for (i, atom) in query.atoms.iter().enumerate() {
        if let Some(mapno) = atom.atom_map {
            idxs.insert(mapno, i);
        }
    }
    let (Some(&left), Some(&right)) = (idxs.get(&site_map[0]), idxs.get(&site_map[1])) else {
        return false;
    };
    let Some(bond) = query.bonds.iter().find(|bond| {
        (bond.atom1 == left && bond.atom2 == right) || (bond.atom1 == right && bond.atom2 == left)
    }) else {
        return false;
    };
    matches!(bond.query, BondQuery::Primitive(BondPrimitive::Double))
}

/// `filter_rules(mol, rule, pattern) -> bool` before SMARTS runs.
pub fn accept_all_rules(_mol: &Molecule, _rule: &RuleSet, _pattern: &PatternInfo) -> bool {
    true
}

/// `filter_sites(mol, site, info) -> bool` after unique-edit, before the edit.
pub fn accept_all_sites(_mol: &Molecule, _site: usize, _info: &SiteInfo) -> bool {
    true
}

/// Python `FilterRules`. Prefer filtering [`Candidate`]s by reading pattern data.
pub type FilterRules = dyn Fn(&Molecule, &RuleSet, &PatternInfo) -> bool;
/// Python `FilterSites`. Prefer filtering [`Candidate`]s by reading site + effect.
pub type FilterSites = dyn Fn(&Molecule, usize, &SiteInfo) -> bool;

/// Stored filters (`Box<dyn Fn>`). Capture-by-move; `'static`.
pub struct BoxedFilters {
    pub rules: Box<FilterRules>,
    pub sites: Box<FilterSites>,
}

impl Default for BoxedFilters {
    fn default() -> Self {
        Self {
            rules: Box::new(accept_all_rules),
            sites: Box::new(accept_all_sites),
        }
    }
}

/// One child of a [`RuleSet`]: a leaf pattern or a nested set (a namespace).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RuleMember {
    Pattern(PatternInfo),
    Set(RuleSet),
}

/// Container of patterns and nested sets. Nested sets stay nested.
#[derive(Clone, Debug)]
pub struct RuleSet {
    pub name: Option<String>,
    /// Leaf-owned expander for [`Self::canonical_plan`] (Python method).
    plan_fn: Option<CanonicalPlanFn>,
    members: Vec<RuleMember>,
}

impl PartialEq for RuleSet {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.members == other.members
            && match (self.plan_fn, other.plan_fn) {
                (None, None) => true,
                (Some(a), Some(b)) => std::ptr::fn_addr_eq(a, b),
                _ => false,
            }
    }
}

impl Eq for RuleSet {}

impl RuleSet {
    pub fn new(name: Option<String>, patterns: impl IntoIterator<Item = PatternInfo>) -> Self {
        Self {
            name,
            plan_fn: None,
            members: patterns.into_iter().map(RuleMember::Pattern).collect(),
        }
    }

    /// Nest the given sets as members. Does not flatten their patterns.
    pub fn compose(name: Option<String>, sets: impl IntoIterator<Item = RuleSet>) -> Self {
        Self {
            name,
            plan_fn: None,
            members: sets.into_iter().map(RuleMember::Set).collect(),
        }
    }

    /// Attach a plan expander (composite leaves: return elementary rule steps).
    pub fn with_canonical_plan(mut self, f: CanonicalPlanFn) -> Self {
        self.plan_fn = Some(f);
        self
    }

    /// Elementary steps for one accepted hop (Python `canonical_plan`).
    ///
    /// Default: identity — this rule at the discovery site. Composites return
    /// steps named after existing catalog rules (`Hydroxylation`, …).
    pub fn canonical_plan(
        &self,
        mol: &Molecule,
        site_atoms: &[usize],
        end_effects: Option<&[&crate::pattern::Effect]>,
    ) -> Vec<Step> {
        let leaf = self.name.as_deref().unwrap_or("");
        steps_for_leaf(self.plan_fn, leaf, mol, site_atoms, end_effects)
    }

    pub fn members(&self) -> &[RuleMember] {
        &self.members
    }

    /// Flat walk of every leaf [`PatternInfo`] under this set (including nested).
    pub fn patterns(&self) -> Vec<&PatternInfo> {
        let mut out = Vec::new();
        self.collect_patterns(&mut out);
        out
    }

    fn collect_patterns<'a>(&'a self, out: &mut Vec<&'a PatternInfo>) {
        for member in &self.members {
            match member {
                RuleMember::Pattern(pattern) => out.push(pattern),
                RuleMember::Set(set) => set.collect_patterns(out),
            }
        }
    }

    pub fn push(&mut self, pattern: PatternInfo) {
        self.members.push(RuleMember::Pattern(pattern));
    }

    pub fn push_set(&mut self, set: RuleSet) {
        self.members.push(RuleMember::Set(set));
    }

    /// Discover site–pattern–parent triples without applying edits.
    ///
    /// Nested sets append themselves onto each candidate's `rule_path`.
    /// Pair-endpoint patterns are collected and resolved into pair emissions
    /// only at [`Self::metabolites`] / [`Self::metabolize`] (path flip still
    /// needs both ends together).
    pub fn candidates(&self, mol: &Molecule) -> Result<Vec<Candidate>, ForestError> {
        self.candidates_inner(mol)
    }

    fn candidates_inner(&self, mol: &Molecule) -> Result<Vec<Candidate>, ForestError> {
        let mut out = Vec::new();
        for member in &self.members {
            match member {
                RuleMember::Pattern(pattern) => {
                    if matches!(pattern.edit, Edit::PairEndpoint(_)) {
                        continue;
                    }
                    out.extend(pattern_candidates(self, mol, pattern)?);
                }
                RuleMember::Set(child) => {
                    for mut c in child.candidates_inner(mol)? {
                        c.rule_path.push(self.name.clone());
                        out.push(c);
                    }
                }
            }
        }
        Ok(out)
    }

    /// Pair-endpoint patterns on this leaf set (not nested).
    pub(crate) fn leaf_pair_endpoints(&self) -> Vec<PatternInfo> {
        self.members
            .iter()
            .filter_map(|m| match m {
                RuleMember::Pattern(p) if matches!(p.edit, Edit::PairEndpoint(_)) => {
                    Some(p.clone())
                }
                _ => None,
            })
            .collect()
    }

    /// Materialize every candidate (and leaf pair paths). No filter closures.
    pub fn metabolites(
        &self,
        mol: &Molecule,
        unique_csmi: bool,
    ) -> Result<Vec<Emission>, ForestError> {
        self.metabolize(mol, accept_all_rules, accept_all_sites, unique_csmi)
    }

    /// ResonancePair path candidates for this set and nested children.
    ///
    /// Discovery only — no path flip. Each candidate carries a merged
    /// [`crate::pattern::Effect`] for filtering before
    /// [`crate::pair_edit::PairCandidate::materialize`].
    pub fn pair_candidates(
        &self,
        mol: &Molecule,
    ) -> Result<Vec<crate::pair_edit::PairCandidate>, ForestError> {
        let mut out = self.pair_candidates_leaf(mol)?;
        for member in &self.members {
            if let RuleMember::Set(child) = member {
                out.extend(child.pair_candidates(mol)?);
            }
        }
        Ok(out)
    }

    /// Pair candidates from this leaf's own endpoint patterns only.
    pub fn pair_candidates_leaf(
        &self,
        mol: &Molecule,
    ) -> Result<Vec<crate::pair_edit::PairCandidate>, ForestError> {
        let endpoints = self.leaf_pair_endpoints();
        if endpoints.is_empty() {
            return Ok(Vec::new());
        }
        pair_candidates(mol, &endpoints)
    }

    /// ResonancePair path emissions for this set and nested children.
    ///
    /// Materializes [`Self::pair_candidates`]. Each emission's `rule_path` is
    /// leaf-first with this set appended when nested.
    pub fn pair_emissions(&self, mol: &Molecule) -> Result<Vec<Emission>, ForestError> {
        let mut out = Vec::new();
        for member in &self.members {
            match member {
                RuleMember::Pattern(_) => {}
                RuleMember::Set(child) => {
                    for mut emission in child.pair_emissions(mol)? {
                        emission.rule_path.push(self.name.clone());
                        out.push(emission);
                    }
                }
            }
        }
        for pair in self.pair_candidates_leaf(mol)? {
            if let Some(emission) = pair.emit(mol)? {
                let site_atoms = pair.plan_site_atoms();
                let ends = [&pair.left.effect, &pair.right.effect];
                let plan = self.canonical_plan(mol, &site_atoms, Some(&ends));
                out.push(Emission {
                    site: emission.site,
                    site_atoms: site_atoms.clone(),
                    cleaves: pair.effect.cleaves,
                    pattern_name: emission.pattern_name,
                    rule_path: vec![self.name.clone()],
                    products: emission.products,
                    plan,
                });
            }
        }
        Ok(out)
    }

    pub fn metabolize_boxed(
        &self,
        mol: &Molecule,
        filters: &BoxedFilters,
        unique_csmi: bool,
    ) -> Result<Vec<Emission>, ForestError> {
        self.metabolize(
            mol,
            |m, r, p| (filters.rules)(m, r, p),
            |m, s, i| (filters.sites)(m, s, i),
            unique_csmi,
        )
    }

    /// Filter candidates by closures (Python parity), then materialize.
    ///
    /// Prefer [`Self::candidates`] + reading [`PatternInfo`] when the search
    /// can filter without callbacks. Nested sets receive `unique_csmi=false`
    /// so alternate children bubble; this set's caller `unique_csmi` is the
    /// cross-child CSMI layer.
    pub fn metabolize<R, S>(
        &self,
        mol: &Molecule,
        filter_rules: R,
        filter_sites: S,
        unique_csmi: bool,
    ) -> Result<Vec<Emission>, ForestError>
    where
        R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
        S: Fn(&Molecule, usize, &SiteInfo) -> bool,
    {
        self.metabolize_inner(mol, &filter_rules, &filter_sites, unique_csmi)
    }

    fn metabolize_inner<R, S>(
        &self,
        mol: &Molecule,
        filter_rules: &R,
        filter_sites: &S,
        unique_csmi: bool,
    ) -> Result<Vec<Emission>, ForestError>
    where
        R: Fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
        S: Fn(&Molecule, usize, &SiteInfo) -> bool,
    {
        let mut emissions = Vec::new();
        let mut seen_csmi: BTreeMap<BTreeSet<String>, String> = BTreeMap::new();
        let mut seen_leaf: BTreeSet<(String, BTreeSet<String>)> = BTreeSet::new();

        // Nested children first (same order as members), then this leaf's
        // patterns — walk members so namespace order matches compose order.
        for member in &self.members {
            match member {
                RuleMember::Pattern(pattern) => {
                    if matches!(pattern.edit, Edit::PairEndpoint(_)) {
                        continue;
                    }
                    if !filter_rules(mol, self, pattern) {
                        continue;
                    }
                    for c in pattern_candidates(self, mol, pattern)? {
                        let info = SiteInfo {
                            site: c.site,
                            pattern: c.pattern.clone(),
                        };
                        if !filter_sites(mol, c.site, &info) {
                            continue;
                        }
                        push_emission(c.emit(mol)?, unique_csmi, &mut seen_leaf, &mut emissions);
                    }
                }
                RuleMember::Set(child) => {
                    let child_emissions =
                        child.metabolize_inner(mol, filter_rules, filter_sites, false)?;
                    for mut emission in child_emissions {
                        emission.rule_path.push(self.name.clone());
                        if unique_csmi {
                            let emission_key: BTreeSet<String> =
                                emission.products.iter().cloned().collect();
                            let leaf_name = emission.leaf_rule().unwrap_or("").to_string();
                            if let Some(kept) = seen_csmi.get(&emission_key) {
                                if kept != &leaf_name {
                                    continue;
                                }
                            } else {
                                seen_csmi.insert(emission_key, leaf_name);
                            }
                        }
                        emissions.push(emission);
                    }
                }
            }
        }

        let pair_endpoints: Vec<PatternInfo> = self
            .leaf_pair_endpoints()
            .into_iter()
            .filter(|p| filter_rules(mol, self, p))
            .collect();
        if !pair_endpoints.is_empty() {
            for pair in crate::pair_edit::pair_candidates(mol, &pair_endpoints)? {
                let info = SiteInfo {
                    site: pair.site,
                    pattern: pair_endpoints[0].clone(),
                };
                if !filter_sites(mol, pair.site, &info) {
                    continue;
                }
                let Some(emission) = pair.emit(mol)? else {
                    continue;
                };
                let site_atoms = pair.plan_site_atoms();
                let ends = [&pair.left.effect, &pair.right.effect];
                let plan = self.canonical_plan(mol, &site_atoms, Some(&ends));
                let emission = Emission {
                    site: emission.site,
                    site_atoms: site_atoms.clone(),
                    cleaves: pair.effect.cleaves,
                    pattern_name: emission.pattern_name,
                    rule_path: vec![self.name.clone()],
                    products: emission.products,
                    plan,
                };
                push_emission(Some(emission), unique_csmi, &mut seen_leaf, &mut emissions);
            }
        }

        Ok(emissions)
    }
}

fn pattern_candidates(
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
        for (mapped, form_i) in hits {
            let Some(&site) = mapped.get(&pattern.primary_map()) else {
                continue;
            };
            out.push(Candidate {
                site,
                pattern: pattern.clone(),
                rule_path: vec![set.name.clone()],
                mapped,
                parent: ParentRef::Form(Box::new(forms[form_i].clone())),
            });
        }
    } else {
        for mapped in unique_sites(mol, &pattern.smarts, pattern.site_kind, &pattern.site_map)? {
            let Some(&site) = mapped.get(&pattern.primary_map()) else {
                continue;
            };
            out.push(Candidate {
                site,
                pattern: pattern.clone(),
                rule_path: vec![set.name.clone()],
                mapped,
                parent: ParentRef::Context,
            });
        }
    }
    Ok(out)
}

fn push_emission(
    emission: Option<Emission>,
    unique_csmi: bool,
    seen_leaf: &mut BTreeSet<(String, BTreeSet<String>)>,
    emissions: &mut Vec<Emission>,
) {
    let Some(emission) = emission else {
        return;
    };
    if unique_csmi {
        let emission_key: BTreeSet<String> = emission.products.iter().cloned().collect();
        let leaf_key = (emission.pattern_name.clone(), emission_key);
        if !seen_leaf.insert(leaf_key) {
            return;
        }
    }
    emissions.push(emission);
}

fn add_hydroxyl(mol: &Molecule, carbon: usize) -> Result<Molecule, ForestError> {
    let (mut product, oxygen) = mol.with_atom_added(Atom::organic(Element::O));
    product
        .add_bond(atom_idx(carbon), oxygen, BondOrder::Single)
        .map_err(|err| ForestError::Smirks(err.to_string()))?;
    Ok(product)
}

/// Same as [`apply_edit_mols`], returning product CSMIs.
pub(crate) fn apply_edit_for_candidate(
    mol: &Molecule,
    pattern: &PatternInfo,
    mapped: &BTreeMap<u16, usize>,
) -> Result<Vec<String>, ForestError> {
    Ok(apply_edit_mols(mol, pattern, mapped)?
        .iter()
        .map(canon_smiles)
        .collect())
}

/// Apply a pattern edit at a mapped site, keeping chematic products (tags intact).
pub(crate) fn apply_edit_mols(
    mol: &Molecule,
    pattern: &PatternInfo,
    mapped: &BTreeMap<u16, usize>,
) -> Result<Vec<Molecule>, ForestError> {
    match &pattern.edit {
        Edit::Hydroxyl => {
            let Some(&carbon) = mapped.get(&1) else {
                return Ok(Vec::new());
            };
            let product = add_hydroxyl(mol, carbon)?;
            if accept_product(&product) {
                Ok(vec![product])
            } else {
                Ok(Vec::new())
            }
        }
        Edit::Smirks(smirks) => {
            let mut cache = crate::kekule::KekuleCache::default();
            let work = crate::kekule::reactant_parent(mol, mapped, smirks, &mut cache)?;
            apply_smirks_at(smirks, &work, mapped)
        }
        Edit::PairEndpoint(_) => Ok(Vec::new()),
    }
}

/// O-dealkylation of a methyl ether (anisole-shaped SMARTS / SMIRKS).
pub fn o_dealkylation() -> RuleSet {
    RuleSet::new(
        Some("Dealkylation".into()),
        [PatternInfo::new(
            "O-Me",
            "[#6H3:1][#8H0:2]",
            Edit::Smirks("[C:1][O:2]>>[O:2].[C:1](=O)O".into()),
            crate::pattern::Effect {
                adds: None,
                removes: None,
                cleaves: true,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
            },
        )],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydroxylation::hydroxylation;
    use crate::mol::{canon_of, parse_mol};
    use crate::pattern::Effect;
    use std::collections::BTreeSet;

    fn canon_set(smiles: impl IntoIterator<Item = impl AsRef<str>>) -> BTreeSet<String> {
        smiles
            .into_iter()
            .map(|s| canon_of(s.as_ref()).unwrap())
            .collect()
    }

    fn products_of(
        set: &RuleSet,
        smiles: &str,
        filter_rules: impl Fn(&Molecule, &RuleSet, &PatternInfo) -> bool,
        filter_sites: impl Fn(&Molecule, usize, &SiteInfo) -> bool,
    ) -> BTreeSet<String> {
        let mol = parse_mol(smiles).unwrap();
        set.metabolize(&mol, filter_rules, filter_sites, true)
            .unwrap()
            .into_iter()
            .flat_map(|e| e.products)
            .map(|s| canon_of(&s).unwrap())
            .collect()
    }

    #[test]
    fn hydroxylation_set_matches_door_on_propane() {
        let got = products_of(&hydroxylation(), "CCC", accept_all_rules, accept_all_sites);
        assert_eq!(got, canon_set(["CCCO", "CC(C)O"]));
    }

    #[test]
    fn filter_rules_closure_drops_a_pattern_by_name() {
        let only_h = |_m: &Molecule, _r: &RuleSet, p: &PatternInfo| p.name == "h";
        let benzene = products_of(&hydroxylation(), "c1ccccc1", only_h, accept_all_sites);
        assert_eq!(benzene, canon_set(["Oc1ccccc1"]));
        let propane = products_of(&hydroxylation(), "CCC", only_h, accept_all_sites);
        assert!(propane.is_empty());
    }

    #[test]
    fn filter_sites_closure_keeps_primary_carbons() {
        let primary = |mol: &Molecule, site: usize, _info: &SiteInfo| {
            mol.neighbors(crate::mol::atom_idx(site))
                .filter(|(nbr, _)| mol.atom(*nbr).element.atomic_number() == 6)
                .count()
                <= 1
        };
        let got = products_of(&hydroxylation(), "CCC", accept_all_rules, primary);
        assert_eq!(got, canon_set(["CCCO"]));
    }

    #[test]
    fn boxed_filters_capture_a_name() {
        let keep = String::from("h2");
        let filters = BoxedFilters {
            rules: Box::new(move |_m, _r, p| p.name == keep),
            sites: Box::new(accept_all_sites),
        };
        let mol = parse_mol("CCC").unwrap();
        let got: BTreeSet<String> = hydroxylation()
            .metabolize_boxed(&mol, &filters, true)
            .unwrap()
            .into_iter()
            .flat_map(|e| e.products)
            .map(|s| canon_of(&s).unwrap())
            .collect();
        assert_eq!(got, canon_set(["CCCO", "CC(C)O"]));
        let benzene = parse_mol("c1ccccc1").unwrap();
        assert!(
            hydroxylation()
                .metabolize_boxed(&benzene, &filters, true)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn filter_rules_reads_methide_as_an_effect_field() {
        let mut pattern = PatternInfo::hydroxyl("methide-ish", "[#6h3:1]");
        pattern.effect = Effect {
            adds: None,
            removes: Some("H".into()),
            cleaves: false,
            methide: true,
            dearomatizes: false,
            leave_count: None,
            partner: None,
        };
        let set = RuleSet::new(Some("probe".into()), [pattern]);
        let refuse_methide = |_m: &Molecule, _r: &RuleSet, p: &PatternInfo| !p.effect.methide;
        assert!(products_of(&set, "CC", refuse_methide, accept_all_sites).is_empty());
        assert!(!products_of(&set, "CC", accept_all_rules, accept_all_sites).is_empty());
    }

    #[test]
    fn compose_nests_sets_as_namespaces() {
        let set = RuleSet::compose(
            Some("PhaseI-probe".into()),
            [hydroxylation(), o_dealkylation()],
        );
        assert_eq!(set.members().len(), 2);
        assert_eq!(set.patterns().len(), 3);
        let anisole = parse_mol("COc1ccccc1").unwrap();
        let emissions = set
            .metabolize(&anisole, accept_all_rules, accept_all_sites, true)
            .unwrap();
        let names: BTreeSet<String> = emissions.iter().map(|e| e.pattern_name.clone()).collect();
        assert!(names.contains("O-Me"));
        assert!(names.contains("h"));
        let smiles: BTreeSet<String> = emissions
            .iter()
            .flat_map(|e| e.products.iter().cloned())
            .map(|s| canon_of(&s).unwrap())
            .collect();
        assert!(smiles.contains(&canon_of("Oc1ccccc1").unwrap()));
        for emission in &emissions {
            assert_eq!(
                emission.rule_path.last().unwrap().as_deref(),
                Some("PhaseI-probe")
            );
            assert!(emission.rule_path.len() >= 2);
        }
        let dealk = emissions.iter().find(|e| e.pattern_name == "O-Me").unwrap();
        assert_eq!(dealk.namespace(), vec!["Dealkylation", "PhaseI-probe"]);
        let hyd = emissions.iter().find(|e| e.pattern_name == "h").unwrap();
        assert_eq!(hyd.namespace(), vec!["Hydroxylation", "PhaseI-probe"]);
    }

    #[test]
    fn nested_ruleset_appends_each_set_on_the_path() {
        let inner = RuleSet::compose(Some("Inner".into()), [hydroxylation()]);
        let outer = RuleSet::compose(Some("Outer".into()), [inner]);
        let mol = parse_mol("CC").unwrap();
        let emissions = outer
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
        assert!(!emissions.is_empty());
        for emission in &emissions {
            assert_eq!(
                emission.namespace(),
                vec!["Hydroxylation", "Inner", "Outer"]
            );
        }
    }

    #[test]
    fn unnamed_outer_stays_on_path_but_drops_from_namespace() {
        let outer = RuleSet::compose(None, [hydroxylation()]);
        let mol = parse_mol("CC").unwrap();
        let emissions = outer
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
        let emission = &emissions[0];
        assert_eq!(emission.rule_path, vec![Some("Hydroxylation".into()), None]);
        assert_eq!(emission.namespace(), vec!["Hydroxylation"]);
    }

    /// Two overlapping leaf rules, same SMARTS / same product — outer unique_csmi
    /// keeps the first leaf (Python RuleSet cross-rule yield).
    #[test]
    fn compose_unique_csmi_keeps_first_leaf_on_overlap() {
        let a = RuleSet::new(
            Some("OverlapOhA".into()),
            [PatternInfo::hydroxyl("oh", "[#6h3:1]")],
        );
        let b = RuleSet::new(
            Some("OverlapOhB".into()),
            [PatternInfo::hydroxyl("oh", "[#6h3:1]")],
        );
        let set = RuleSet::compose(Some("OverlapSet".into()), [a, b]);
        let mol = parse_mol("CC").unwrap();
        let with_dedup = set
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
        assert_eq!(with_dedup.len(), 1);
        assert_eq!(with_dedup[0].leaf_rule(), Some("OverlapOhA"));
        assert_eq!(with_dedup[0].namespace(), vec!["OverlapOhA", "OverlapSet"]);
        assert_eq!(
            canon_set(with_dedup[0].products.iter().cloned()),
            canon_set(["CCO"])
        );

        let without = set
            .metabolize(&mol, accept_all_rules, accept_all_sites, false)
            .unwrap();
        assert_eq!(without.len(), 2);
        assert_eq!(without[0].leaf_rule(), Some("OverlapOhA"));
        assert_eq!(without[1].leaf_rule(), Some("OverlapOhB"));
    }

    #[test]
    fn filter_rules_sees_leaf_set_not_outer_compose() {
        let set = RuleSet::compose(Some("Forest".into()), [hydroxylation(), o_dealkylation()]);
        let seen = std::cell::RefCell::new(BTreeSet::new());
        let mol = parse_mol("CC").unwrap();
        let _ = set
            .metabolize(
                &mol,
                |_m, rule, _p| {
                    if let Some(name) = &rule.name {
                        seen.borrow_mut().insert(name.clone());
                    }
                    true
                },
                accept_all_sites,
                true,
            )
            .unwrap();
        let seen = seen.into_inner();
        assert!(seen.contains("Hydroxylation"));
        assert!(seen.contains("Dealkylation"));
        assert!(!seen.contains("Forest"));
    }

    #[test]
    fn nested_outermost_yield_only_dedups_across_inner_children() {
        let a = RuleSet::new(
            Some("OverlapOhA".into()),
            [PatternInfo::hydroxyl("oh", "[#6h3:1]")],
        );
        let b = RuleSet::new(
            Some("OverlapOhB".into()),
            [PatternInfo::hydroxyl("oh", "[#6h3:1]")],
        );
        let inner = RuleSet::compose(Some("Inner".into()), [a, b]);
        let outer = RuleSet::compose(Some("Outer".into()), [inner]);
        let mol = parse_mol("CC").unwrap();
        let products = outer
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
        assert_eq!(products.len(), 1);
        assert_eq!(
            products[0].namespace(),
            vec!["OverlapOhA", "Inner", "Outer"]
        );
    }

    #[test]
    fn dehydrogenation_metabolize_emits_quinone_via_pair_door() {
        use crate::rules::dehydrogenation;
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let emissions = dehydrogenation()
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .unwrap();
        let want = canon_of("O=C1C=CC(=O)C=C1").unwrap();
        assert!(
            emissions.iter().any(|e| {
                e.products.iter().any(|p| canon_of(p).unwrap() == want)
                    && e.leaf_rule() == Some("Dehydrogenation")
            }),
            "{emissions:?}"
        );
    }

    #[test]
    fn epoxidation_matches_kekule_forms_on_benzene() {
        use crate::rules::epoxidation;
        let mol = parse_mol("c1ccccc1").unwrap();
        let candidates = epoxidation().candidates(&mol).unwrap();
        assert_eq!(candidates.len(), 1, "{candidates:?}");
        assert!(matches!(
            candidates[0].parent,
            crate::candidate::ParentRef::Form(_)
        ));
        let emissions = epoxidation().metabolites(&mol, true).unwrap();
        assert_eq!(emissions.len(), 1, "{emissions:?}");
        assert_eq!(emissions[0].pattern_name, "epoxide");
        assert_eq!(emissions[0].leaf_rule(), Some("Epoxidation"));
        let want = canon_of("C1=CC2OC2C=C1").unwrap();
        assert!(
            emissions[0]
                .products
                .iter()
                .any(|p| canon_of(p).unwrap() == want),
            "want {want}, got {:?}",
            emissions[0].products
        );
    }

    #[test]
    fn epoxidation_on_alkene_does_not_need_forms() {
        use crate::rules::epoxidation;
        let got = products_of(&epoxidation(), "C=C", accept_all_rules, accept_all_sites);
        assert_eq!(got, canon_set(["C1CO1"]));
    }

    #[test]
    fn candidates_defer_materialize() {
        let mol = parse_mol("CC").unwrap();
        let set = hydroxylation();
        let cands = set.candidates(&mol).unwrap();
        assert_eq!(cands.len(), 1);
        assert_eq!(cands[0].pattern.name, "h2");
        // Filtering by pattern data needs no closure into the rule.
        let refuse: Vec<_> = cands
            .iter()
            .filter(|c| c.pattern.effect.adds.as_deref() != Some("O"))
            .collect();
        assert!(refuse.is_empty());
        let products = cands[0].materialize(&mol).unwrap();
        assert_eq!(canon_of(&products[0]).unwrap(), canon_of("CCO").unwrap());
    }
}
