//! A set of [`PatternInfo`] records and nested sets, as namespaces.
//!
//! Nested sets stay nested: [`RuleSet::compose`] does not flatten. Each set
//! appends itself onto the emission [`crate::pattern::Emission::rule_path`]
//! after the child that emitted (leaf first, outer last), matching Python
//! `info["rule"]` / addition chain order.
//!
//! Python `FilterRules` / `FilterSites` are `Callable`. Here they are
//! `impl Fn` on [`RuleSet::metabolize`], or [`BoxedFilters`] when a search
//! needs to store them. Built-in filters should read [`PatternInfo`] fields.
//! Filters still see each leaf set (not the outer container). A Python lambda
//! still crosses the GIL.

use std::collections::{BTreeMap, BTreeSet};

use chematic::core::{Atom, BondOrder, Element};

use crate::ForestError;
use crate::mol::{Molecule, atom_idx, canon_smiles};
use crate::pattern::{Edit, Emission, PatternInfo, SiteInfo};
use crate::smirks::apply_smirks_at;
use crate::unique_edit::unique_atom_sites;
use crate::valence::accept_product;

/// `filter_rules(mol, rule, pattern) -> bool` before SMARTS runs.
pub fn accept_all_rules(_mol: &Molecule, _rule: &RuleSet, _pattern: &PatternInfo) -> bool {
    true
}

/// `filter_sites(mol, site, info) -> bool` after unique-edit, before the edit.
pub fn accept_all_sites(_mol: &Molecule, _site: usize, _info: &SiteInfo) -> bool {
    true
}

/// Python `FilterRules`. Use `impl Fn` on [`RuleSet::metabolize`] when possible.
pub type FilterRules = dyn Fn(&Molecule, &RuleSet, &PatternInfo) -> bool;
/// Python `FilterSites`.
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RuleSet {
    pub name: Option<String>,
    members: Vec<RuleMember>,
}

impl RuleSet {
    pub fn new(name: Option<String>, patterns: impl IntoIterator<Item = PatternInfo>) -> Self {
        Self {
            name,
            members: patterns
                .into_iter()
                .map(RuleMember::Pattern)
                .collect(),
        }
    }

    /// Nest the given sets as members. Does not flatten their patterns.
    pub fn compose(name: Option<String>, sets: impl IntoIterator<Item = RuleSet>) -> Self {
        Self {
            name,
            members: sets.into_iter().map(RuleMember::Set).collect(),
        }
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

    /// Run each member. Nested sets receive `unique_csmi=false` so alternate
    /// children bubble up; this set's caller `unique_csmi` is the cross-child
    /// CSMI layer. Each emission's [`Emission::rule_path`] ends with this set.
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
        // Cross-child yield: emission CSMI frozenset → kept leaf rule name.
        let mut seen_csmi: BTreeMap<BTreeSet<String>, String> = BTreeMap::new();
        // Within this leaf's own patterns (when unique_csmi): (pattern, emission).
        let mut seen_leaf: BTreeSet<(String, BTreeSet<String>)> = BTreeSet::new();

        for member in &self.members {
            match member {
                RuleMember::Pattern(pattern) => {
                    // Pair-endpoint paths are data until the pair door is wired.
                    if matches!(pattern.edit, Edit::PairEndpoint(_)) {
                        continue;
                    }
                    if !filter_rules(mol, self, pattern) {
                        continue;
                    }
                    for mapped in unique_atom_sites(mol, &pattern.smarts)? {
                        let Some(&site) = mapped.get(&pattern.primary_map()) else {
                            continue;
                        };
                        let info = SiteInfo {
                            site,
                            pattern: pattern.clone(),
                        };
                        if !filter_sites(mol, site, &info) {
                            continue;
                        }
                        let products = apply_edit(mol, pattern, &mapped)?;
                        if products.is_empty() {
                            continue;
                        }
                        let emission_key: BTreeSet<String> = products.iter().cloned().collect();
                        if unique_csmi {
                            let leaf_key = (pattern.name.clone(), emission_key.clone());
                            if !seen_leaf.insert(leaf_key) {
                                continue;
                            }
                        }
                        emissions.push(Emission {
                            site,
                            pattern_name: pattern.name.clone(),
                            rule_path: vec![self.name.clone()],
                            products,
                        });
                    }
                }
                RuleMember::Set(child) => {
                    // Children: yield off so alternate rules / nested sets bubble.
                    let child_emissions =
                        child.metabolize_inner(mol, filter_rules, filter_sites, false)?;
                    for mut emission in child_emissions {
                        emission.rule_path.push(self.name.clone());
                        if unique_csmi {
                            let emission_key: BTreeSet<String> =
                                emission.products.iter().cloned().collect();
                            let leaf_name = emission
                                .leaf_rule()
                                .unwrap_or("")
                                .to_string();
                            if let Some(kept) = seen_csmi.get(&emission_key) {
                                if kept != &leaf_name {
                                    // Overlapping coverage: keep first leaf only.
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
        Ok(emissions)
    }
}

fn add_hydroxyl(mol: &Molecule, carbon: usize) -> Result<Molecule, ForestError> {
    let (mut product, oxygen) = mol.with_atom_added(Atom::organic(Element::O));
    product
        .add_bond(atom_idx(carbon), oxygen, BondOrder::Single)
        .map_err(|err| ForestError::Smirks(err.to_string()))?;
    Ok(product)
}

fn apply_edit(
    mol: &Molecule,
    pattern: &PatternInfo,
    mapped: &std::collections::BTreeMap<u16, usize>,
) -> Result<Vec<String>, ForestError> {
    match &pattern.edit {
        Edit::Hydroxyl => {
            let Some(&carbon) = mapped.get(&1) else {
                return Ok(Vec::new());
            };
            let product = add_hydroxyl(mol, carbon)?;
            if accept_product(&product) {
                Ok(vec![canon_smiles(&product)])
            } else {
                Ok(Vec::new())
            }
        }
        Edit::Smirks(smirks) => {
            let pieces = apply_smirks_at(smirks, mol, mapped)?;
            Ok(pieces.iter().map(canon_smiles).collect())
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
        assert_eq!(
            dealk.namespace(),
            vec!["Dealkylation", "PhaseI-probe"]
        );
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
        assert_eq!(
            emission.rule_path,
            vec![Some("Hydroxylation".into()), None]
        );
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
        assert_eq!(
            with_dedup[0].namespace(),
            vec!["OverlapOhA", "OverlapSet"]
        );
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
        let set = RuleSet::compose(
            Some("Forest".into()),
            [hydroxylation(), o_dealkylation()],
        );
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
}
