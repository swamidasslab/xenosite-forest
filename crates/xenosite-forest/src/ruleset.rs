//! A set of [`PatternInfo`] records, and the closures that filter them.
//!
//! Python `FilterRules` / `FilterSites` are `Callable`. Here they are
//! `impl Fn` on [`RuleSet::metabolize`], or [`BoxedFilters`] when a search
//! needs to store them. Built-in filters should read [`PatternInfo`] fields.
//! A Python lambda still crosses the GIL.

use std::collections::BTreeSet;

use chematic::core::{Atom, BondOrder, Element};

use crate::ForestError;
use crate::mol::{Molecule, atom_idx, canon_smiles};
use crate::pattern::{Edit, Emission, PatternInfo, SiteInfo, SiteKind};
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

/// Container of patterns. Nested sets flatten: compose concatenates data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RuleSet {
    pub name: Option<String>,
    patterns: Vec<PatternInfo>,
}

impl RuleSet {
    pub fn new(name: Option<String>, patterns: impl IntoIterator<Item = PatternInfo>) -> Self {
        Self {
            name,
            patterns: patterns.into_iter().collect(),
        }
    }

    /// Concatenate sets. The algorithm still walks a flat `Vec<PatternInfo>`.
    pub fn compose(name: Option<String>, sets: impl IntoIterator<Item = RuleSet>) -> Self {
        Self {
            name,
            patterns: sets.into_iter().flat_map(|set| set.patterns).collect(),
        }
    }

    pub fn patterns(&self) -> &[PatternInfo] {
        &self.patterns
    }

    pub fn push(&mut self, pattern: PatternInfo) {
        self.patterns.push(pattern);
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

    /// Run each pattern: `filter_rules` → unique-edit → `filter_sites` → edit.
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
        let mut emissions = Vec::new();
        let mut seen: BTreeSet<BTreeSet<String>> = BTreeSet::new();
        for pattern in &self.patterns {
            if !filter_rules(mol, self, pattern) {
                continue;
            }
            let SiteKind::Atom = pattern.site_kind;
            for mapped in unique_atom_sites(mol, &pattern.smarts)? {
                let Some(&site) = mapped.get(&1) else {
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
                if unique_csmi {
                    let key: BTreeSet<String> = products.iter().cloned().collect();
                    if !seen.insert(key) {
                        continue;
                    }
                }
                emissions.push(Emission {
                    site,
                    pattern_name: pattern.name.clone(),
                    products,
                });
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
    fn compose_concatenates_pattern_data() {
        let set = RuleSet::compose(
            Some("PhaseI-probe".into()),
            [hydroxylation(), o_dealkylation()],
        );
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
    }
}
