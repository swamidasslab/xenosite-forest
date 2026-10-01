//! Deferred site hits: any [`SiteKind`] SOM, edit on demand.
//!
//! A search filters these (via [`Self::info`] or the soft field aliases), then
//! [`DeferredSite::apply`] — run the edit and get an [`Emission`] whose
//! products are [`ForestMol`] adopted from the discovery parent. Tags and
//! kekulé caches propagate transparently; CSMI is never the default product
//! type. Pair and atom/bond sites share this type.

use std::collections::{BTreeMap, HashSet};
use std::rc::Rc;

use crate::ForestError;
use crate::ForestMol;
use crate::canonical_plan::{Step, identity_plan_on_forest};
use crate::mol::Molecule;
use crate::pattern::{Edit, Effect, Emission, PatternInfo, SiteInfo, SiteKind};
use crate::ruleset::apply_edit_mols;

/// Soft rename — prefer [`DeferredSite`].
pub type Candidate = DeferredSite;

/// Which mol to apply the edit on.
///
/// [`ParentRef::Context`] is the aromatic (or input) parent — ResonanceRule
/// SMIRKS still go through Kekulé `reactant_parent` at materialize time.
/// [`ParentRef::Form`] is an exclusive-double match already on a Kekulé writing
/// (Epoxidation); apply runs on that form directly.
#[derive(Clone)]
pub enum ParentRef {
    Context,
    Form(Box<Molecule>),
}

impl std::fmt::Debug for ParentRef {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Context => f.write_str("Context"),
            Self::Form(mol) => write!(f, "Form(atoms={})", mol.atom_count()),
        }
    }
}

/// Resonance-pair payload (absent for atom/bond sites).
#[derive(Clone, Debug)]
pub(crate) struct PairPayload {
    pub left: PatternInfo,
    pub right: PatternInfo,
    /// Merged end effects (dearomatizes resolved against system aromaticity).
    pub effect: Effect,
    pub map2: BTreeMap<u16, usize>,
    pub start: usize,
    pub end: usize,
    pub system: HashSet<usize>,
}

/// One discovered site that has not yet been edited.
///
/// Holds any SOM. Soft field aliases (`site`, `orbit`, `pattern`, `mapped`,
/// `parent`) match the old [`Candidate`] shape so call sites migrate gradually
/// onto [`Self::info`] / methods. The discovery [`ForestMol`] is held behind
/// [`Rc`]; apply adopts products from it so atom tags and caches continue.
#[derive(Clone)]
pub struct DeferredSite {
    mol: Rc<ForestMol>,
    /// Filter bag — same content as the soft aliases below.
    pub info: SiteInfo,
    /// Soft alias of [`SiteInfo::site`].
    pub site: usize,
    /// Soft alias of [`SiteInfo::orbit`].
    pub orbit: Vec<usize>,
    /// Soft alias of [`SiteInfo::pattern`].
    pub pattern: PatternInfo,
    /// Leaf-first rule namespace (emitting set, then containers).
    pub rule_path: Vec<Option<String>>,
    /// Atom/bond SMARTS map; for pairs, the left-end map.
    pub mapped: BTreeMap<u16, usize>,
    pub parent: ParentRef,
    /// Soft alias — combined name for pairs (`left+right`), else pattern name.
    pub pattern_name: String,
    /// Soft alias — left endpoint (atom sites: same as [`Self::pattern`]).
    pub left: PatternInfo,
    /// Soft alias — right endpoint (atom sites: same as [`Self::pattern`]).
    pub right: PatternInfo,
    /// Soft alias — merged pair effect, else pattern effect.
    pub effect: Effect,
    pub(crate) pair: Option<PairPayload>,
}

impl std::fmt::Debug for DeferredSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DeferredSite")
            .field("site", &self.site)
            .field("orbit", &self.orbit)
            .field("pattern", &self.pattern.name)
            .field("rule_path", &self.rule_path)
            .field("is_pair", &self.pair.is_some())
            .field("mol_atoms", &self.mol.mol().atom_count())
            .finish()
    }
}

impl DeferredSite {
    pub(crate) fn atom(
        mol: Rc<ForestMol>,
        site: usize,
        orbit: Vec<usize>,
        pattern: PatternInfo,
        rule_path: Vec<Option<String>>,
        mapped: BTreeMap<u16, usize>,
        parent: ParentRef,
    ) -> Self {
        let info = SiteInfo {
            site,
            orbit: orbit.clone(),
            pattern: pattern.clone(),
            shell_forecast: None,
        };
        Self {
            mol,
            info,
            site,
            orbit,
            pattern: pattern.clone(),
            rule_path,
            mapped,
            parent,
            pattern_name: pattern.name.clone(),
            left: pattern.clone(),
            right: pattern.clone(),
            effect: pattern.effect.clone(),
            pair: None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pair(
        mol: Rc<ForestMol>,
        site: usize,
        pattern_name: String,
        left: PatternInfo,
        right: PatternInfo,
        effect: Effect,
        map1: BTreeMap<u16, usize>,
        map2: BTreeMap<u16, usize>,
        start: usize,
        end: usize,
        system: HashSet<usize>,
        rule_path: Vec<Option<String>>,
    ) -> Self {
        let search_bias = left.search_bias.min(right.search_bias);
        let mut pattern = left.clone();
        pattern.name = pattern_name.clone();
        pattern.effect = effect.clone();
        pattern.search_bias = search_bias;
        pattern.site_kind = SiteKind::AtomPair;
        let orbit = vec![site];
        let info = SiteInfo {
            site,
            orbit: orbit.clone(),
            pattern: pattern.clone(),
            shell_forecast: None,
        };
        Self {
            mol,
            info,
            site,
            orbit,
            pattern,
            rule_path,
            mapped: map1,
            parent: ParentRef::Context,
            pattern_name,
            left: left.clone(),
            right: right.clone(),
            effect: effect.clone(),
            pair: Some(PairPayload {
                left,
                right,
                effect,
                map2,
                start,
                end,
                system,
            }),
        }
    }

    /// Chemistry graph (discovery parent).
    pub fn mol(&self) -> &Molecule {
        self.mol.mol()
    }

    /// Discovery [`ForestMol`] (tags + caches). Apply adopts from this.
    pub fn forest(&self) -> &ForestMol {
        &self.mol
    }

    pub fn forest_rc(&self) -> &Rc<ForestMol> {
        &self.mol
    }

    pub fn pattern_name(&self) -> &str {
        &self.pattern_name
    }

    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path.first().and_then(|n| n.as_deref())
    }

    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    pub fn is_pair(&self) -> bool {
        self.pair.is_some()
    }

    pub fn is_pair_endpoint(&self) -> bool {
        matches!(self.pattern.edit, Edit::PairEndpoint(_))
    }

    /// Discovery site atoms for each pair end (Python `end_atoms`).
    pub fn end_atoms(&self) -> Option<(usize, usize)> {
        let p = self.pair.as_ref()?;
        let a = self.mapped.get(&p.left.primary_map()).copied()?;
        let b = p.map2.get(&p.right.primary_map()).copied()?;
        Some((a, b))
    }

    /// Conjugated-system anchors for the alternating path (Python `path_ends`).
    pub fn path_ends(&self) -> (usize, usize) {
        match &self.pair {
            Some(p) => (p.start, p.end),
            None => (self.site, self.site),
        }
    }

    /// Elementary plan site atoms (pair ends or unique-edit site maps).
    pub fn site_atoms(&self) -> Vec<usize> {
        if self.pair.is_some() {
            match self.end_atoms() {
                Some((a, b)) => vec![a, b],
                None => vec![self.site],
            }
        } else {
            let mut atoms: Vec<usize> = self
                .pattern
                .site_map
                .iter()
                .filter_map(|m| self.mapped.get(m).copied())
                .collect();
            if atoms.is_empty() {
                atoms.push(self.site);
            }
            atoms
        }
    }

    /// Alias for [`Self::site_atoms`] (pair call-site parity).
    pub fn plan_site_atoms(&self) -> Vec<usize> {
        self.site_atoms()
    }

    /// Atoms recorded on a hop / step: unique-edit orbit for atom/bond sites,
    /// sorted pair end atoms for AtomPair.
    pub fn discovery_atoms(&self) -> Vec<usize> {
        if self.pair.is_some() {
            let mut v = self.site_atoms();
            v.sort_unstable();
            v.dedup();
            v
        } else {
            self.orbit.clone()
        }
    }

    /// Pattern bag for pathway bookkeeping — left endpoint for pairs (combined
    /// name lives on [`Self::pattern_name`]), else the site pattern.
    pub fn bookkeeping_pattern(&self) -> &PatternInfo {
        &self.left
    }

    /// Low-level graph edit only (chematic products, no ForestMol adopt).
    ///
    /// Prefer [`Self::apply`] — it returns [`Emission`] with tagged
    /// [`ForestMol`] products. Use this only when a caller already has a
    /// different adopt parent (legacy find_path emit helpers).
    ///
    /// Non-cleaving Effect-formula filter is **off** while chasing native↔Rust
    /// product parity (C18). Catalog Effect accuracy stays in
    /// `pattern_info_catalog`; re-enable the Keep-H drop after parity is green.
    pub fn materialize_mols(&self) -> Result<Vec<Molecule>, ForestError> {
        let mols = match &self.pair {
            None => {
                let work = match &self.parent {
                    ParentRef::Context => self.mol(),
                    ParentRef::Form(form) => form.as_ref(),
                };
                apply_edit_mols(work, &self.pattern, &self.mapped)?
            }
            Some(p) => crate::pair_edit::materialize_pair_mols(
                self.mol(),
                &p.left,
                &p.right,
                &p.effect,
                &self.mapped,
                &p.map2,
                p.start,
                p.end,
                &p.system,
            )?,
        };
        Ok(mols)
    }

    /// Apply the edit; return a metabolize [`Emission`] with tagged products.
    ///
    /// Products are [`ForestMol`] adopted from the discovery parent — atom
    /// tracking and kekulé caches continue without an extra call. Returns
    /// `None` when the edit yields no acceptable products.
    ///
    /// CSMI is not the product type. Use [`Emission::product_csmis`] only when
    /// a string identity is an explicit choice.
    pub fn apply(&self) -> Result<Option<Emission>, ForestError> {
        let mols = self.materialize_mols()?;
        if mols.is_empty() {
            return Ok(None);
        }
        // Effect-formula soft check deferred while chasing product parity.
        let star_label = self
            .leaf_rule()
            .and_then(crate::star_conjugate::conjugate_star_label);
        let products = mols
            .into_iter()
            .map(|piece| {
                let product = self.mol.from_edit_product(piece);
                match star_label {
                    Some(label) => {
                        crate::star_conjugate::collapse_conjugate_to_star(&product, label)
                    }
                    None => product,
                }
            })
            .collect::<Vec<_>>();
        Ok(Some(Emission {
            site: self.site,
            site_orbit: self.orbit.clone(),
            site_atoms: self.site_atoms(),
            cleaves: self.effect.cleaves,
            pattern_name: self.pattern_name.clone(),
            search_bias: self.pattern.search_bias,
            rule_path: self.rule_path.clone(),
            products,
            reactant: (*self.mol).clone(),
            plan: self.elementary_plan(),
        }))
    }

    /// Soft alias of [`Self::apply`].
    pub fn emit(&self) -> Result<Option<Emission>, ForestError> {
        self.apply()
    }

    /// Elementary plan: leaf `canonical_plan` hook when present, else identity.
    ///
    /// Pair sites pass both end effects into the hook (quinone-shaped expansion).
    /// Atom/bond sites pass `None` for ends. One door for Expand / metabolize.
    pub fn elementary_plan(&self) -> Vec<Step> {
        let atoms = self.site_atoms();
        match self
            .leaf_rule()
            .and_then(crate::rules::leaf_rule)
            .filter(|leaf| leaf.has_plan_hook())
        {
            Some(leaf) => {
                if self.is_pair() {
                    let ends = [&self.left.effect, &self.right.effect];
                    leaf.canonical_plan(self.mol(), &atoms, Some(&ends))
                } else {
                    leaf.canonical_plan(self.mol(), &atoms, None)
                }
            }
            None => self.identity_plan(),
        }
    }

    /// Cleavage side signature: both ends must agree for pairs.
    pub fn cleave_side_sig(&self) -> crate::pattern::CleaveSideSig {
        if self.is_pair() {
            let left = self.left.cleave_side_sig();
            let right = self.right.cleave_side_sig();
            if left == right {
                left
            } else {
                crate::pattern::CleaveSideSig::Ungrouped
            }
        } else {
            self.pattern.cleave_side_sig()
        }
    }

    /// Elementary plan for this site (cached generators on the discovery ForestMol).
    pub fn identity_plan(&self) -> Vec<Step> {
        self.identity_plan_with_gens(&self.mol.atom_bond_generators(), self.mol().atom_count())
    }

    /// Like [`Self::identity_plan`], reusing supplied generators.
    pub fn identity_plan_with_gens(
        &self,
        generators: &[crate::orbits::AtomBondGenerator],
        n_atoms: usize,
    ) -> Vec<Step> {
        let rule = self
            .leaf_rule()
            .unwrap_or(self.pattern.name.as_str())
            .to_string();
        if self.pair.is_some() {
            let atoms = self.site_atoms();
            identity_plan_on_forest(rule, &self.mol, atoms.clone(), atoms)
        } else {
            let orbit = crate::orbits::atom_orbit_with_gens(generators, n_atoms, self.site);
            identity_plan_on_forest(rule, &self.mol, [self.site], orbit)
        }
    }
}
