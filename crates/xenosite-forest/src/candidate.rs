//! Deferred site hits: edit or ResonancePair, materialize on demand.
//!
//! Prefer [`crate::ruleset::RuleSet::candidates`] /
//! [`crate::ruleset::RuleSet::metabolites`]. Those doors yield SMIRKS edits and
//! ResonancePair hits under [`Candidate`] — callers do not call pair-specific
//! discovery. Leaf-first [`Candidate::rule_path`] is stamped by the RuleSet.

use std::collections::BTreeMap;

use crate::ForestError;
use crate::canonical_plan::{Step, identity_plan_at_indexes};
use crate::mol::Molecule;
use crate::pair_edit::PairCandidate;
use crate::pattern::{Edit, Effect, Emission, PatternInfo};
use crate::ruleset::apply_edit_mols;

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

/// One SMIRKS / atom-edit discovery (not a ResonancePair).
#[derive(Clone, Debug)]
pub struct EditCandidate {
    pub site: usize,
    /// Primary-map atoms collapsed into this unique-edit class (sorted).
    pub orbit: Vec<usize>,
    pub pattern: PatternInfo,
    /// Leaf-first rule namespace (emitting set, then containers).
    pub rule_path: Vec<Option<String>>,
    pub mapped: BTreeMap<u16, usize>,
    pub parent: ParentRef,
}

/// Deferred hit from [`crate::ruleset::RuleSet::candidates`]: SMIRKS edit or
/// ResonancePair. Shared accessors are the polymorphic surface.
#[derive(Clone, Debug)]
pub enum Candidate {
    Edit(EditCandidate),
    Pair(PairCandidate),
}

impl Candidate {
    pub fn edit(
        site: usize,
        orbit: Vec<usize>,
        pattern: PatternInfo,
        rule_path: Vec<Option<String>>,
        mapped: BTreeMap<u16, usize>,
        parent: ParentRef,
    ) -> Self {
        Self::Edit(EditCandidate {
            site,
            orbit,
            pattern,
            rule_path,
            mapped,
            parent,
        })
    }

    pub fn from_pair(pair: PairCandidate) -> Self {
        Self::Pair(pair)
    }

    pub fn as_edit(&self) -> Option<&EditCandidate> {
        match self {
            Self::Edit(e) => Some(e),
            Self::Pair(_) => None,
        }
    }

    pub fn as_pair(&self) -> Option<&PairCandidate> {
        match self {
            Self::Edit(_) => None,
            Self::Pair(p) => Some(p),
        }
    }

    pub fn is_pair(&self) -> bool {
        matches!(self, Self::Pair(_))
    }

    pub fn site(&self) -> usize {
        match self {
            Self::Edit(e) => e.site,
            Self::Pair(p) => p.site,
        }
    }

    pub fn orbit(&self) -> &[usize] {
        match self {
            Self::Edit(e) => &e.orbit,
            Self::Pair(p) => std::slice::from_ref(&p.site),
        }
    }

    pub fn pattern_name(&self) -> &str {
        match self {
            Self::Edit(e) => e.pattern.name.as_str(),
            Self::Pair(p) => p.pattern_name.as_str(),
        }
    }

    /// SMIRKS pattern when this is an edit hit (`None` for ResonancePair).
    pub fn pattern(&self) -> Option<&PatternInfo> {
        self.as_edit().map(|e| &e.pattern)
    }

    pub fn effect(&self) -> &Effect {
        match self {
            Self::Edit(e) => &e.pattern.effect,
            Self::Pair(p) => &p.effect,
        }
    }

    pub fn search_bias(&self) -> i8 {
        match self {
            Self::Edit(e) => e.pattern.search_bias,
            Self::Pair(p) => p.left.search_bias.min(p.right.search_bias),
        }
    }

    pub fn rule_path(&self) -> &[Option<String>] {
        match self {
            Self::Edit(e) => &e.rule_path,
            Self::Pair(p) => &p.rule_path,
        }
    }

    pub fn rule_path_mut(&mut self) -> &mut Vec<Option<String>> {
        match self {
            Self::Edit(e) => &mut e.rule_path,
            Self::Pair(p) => &mut p.rule_path,
        }
    }

    pub fn leaf_rule(&self) -> Option<&str> {
        self.rule_path().first().and_then(|n| n.as_deref())
    }

    /// Hop / emission rule label: leaf set name, else pattern / pair name.
    pub fn rule_name(&self) -> &str {
        self.leaf_rule().unwrap_or_else(|| self.pattern_name())
    }

    pub fn namespace(&self) -> Vec<&str> {
        self.rule_path()
            .iter()
            .filter_map(|name| name.as_deref())
            .collect()
    }

    /// Discovery site atoms (pattern maps, or pair ends).
    pub fn site_atoms(&self) -> Vec<usize> {
        match self {
            Self::Edit(e) => {
                let mut atoms: Vec<usize> = e
                    .pattern
                    .site_map
                    .iter()
                    .filter_map(|m| e.mapped.get(m).copied())
                    .collect();
                if atoms.is_empty() {
                    atoms.push(e.site);
                }
                atoms
            }
            Self::Pair(p) => p.plan_site_atoms(),
        }
    }

    /// Run the edit / path flip and return product molecules.
    pub fn materialize_mols(&self, context: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        match self {
            Self::Edit(e) => {
                let work = match &e.parent {
                    ParentRef::Context => context,
                    ParentRef::Form(form) => form,
                };
                apply_edit_mols(work, &e.pattern, &e.mapped)
            }
            Self::Pair(p) => p.materialize_mols(context),
        }
    }

    /// Run the edit and return product CSMIs. Call only after filtering.
    pub fn materialize(&self, context: &Molecule) -> Result<Vec<String>, ForestError> {
        match self {
            Self::Edit(e) => {
                use crate::ruleset::apply_edit_for_candidate;
                let work = match &e.parent {
                    ParentRef::Context => context,
                    ParentRef::Form(form) => form,
                };
                apply_edit_for_candidate(work, &e.pattern, &e.mapped)
            }
            Self::Pair(p) => p.materialize(context),
        }
    }

    /// Cleavage Or-fold signature (edit pattern, or merged pair ends).
    pub fn cleave_side_sig(&self) -> crate::pattern::CleaveSideSig {
        match self {
            Self::Edit(e) => e.pattern.cleave_side_sig(),
            Self::Pair(p) => {
                let left = p.left.cleave_side_sig();
                let right = p.right.cleave_side_sig();
                if left == right {
                    left
                } else {
                    crate::pattern::CleaveSideSig::Ungrouped
                }
            }
        }
    }

    /// Path ends for dehydrogenation pair progress (`None` for edits / non-DH).
    pub fn dh_ends(&self) -> Option<(usize, usize)> {
        match self {
            Self::Edit(_) => None,
            Self::Pair(p) => {
                if crate::atom_diff::is_dehydrogenation_effect(&p.effect) {
                    p.end_atoms()
                } else {
                    None
                }
            }
        }
    }

    /// Path anchors for H-progress (pair conjugated ends; empty for edits).
    pub fn path_ends(&self) -> Vec<usize> {
        match self {
            Self::Edit(_) => Vec::new(),
            Self::Pair(p) => {
                let (a, b) = p.path_ends();
                vec![a, b]
            }
        }
    }

    /// Plan for this hit: leaf `canonical_plan` when hooked, else identity.
    pub fn plan(&self, mol: &Molecule) -> Vec<Step> {
        match self {
            Self::Edit(e) => {
                match e
                    .rule_path
                    .first()
                    .and_then(|n| n.as_deref())
                    .and_then(crate::rules::leaf_rule)
                    .filter(|leaf| leaf.has_plan_hook())
                {
                    Some(leaf) => {
                        let atoms = self.site_atoms();
                        leaf.canonical_plan(mol, &atoms, None)
                    }
                    None => e.identity_plan(mol),
                }
            }
            Self::Pair(p) => pair_plan(p, mol),
        }
    }

    /// Like [`Self::plan`], preferring ForestMol-cached generators for identity.
    pub fn plan_with_gens(
        &self,
        generators: &[crate::orbits::AtomBondGenerator],
        n_atoms: usize,
        mol: &Molecule,
    ) -> Vec<Step> {
        match self {
            Self::Edit(e) => {
                match e
                    .rule_path
                    .first()
                    .and_then(|n| n.as_deref())
                    .and_then(crate::rules::leaf_rule)
                    .filter(|leaf| leaf.has_plan_hook())
                {
                    Some(leaf) => {
                        let atoms = self.site_atoms();
                        leaf.canonical_plan(mol, &atoms, None)
                    }
                    None => e.identity_plan_with_gens(generators, n_atoms, mol),
                }
            }
            Self::Pair(p) => pair_plan(p, mol),
        }
    }

    /// Materialize into an [`Emission`] (site + namespace + products).
    pub fn emit(&self, context: &Molecule) -> Result<Option<Emission>, ForestError> {
        let mols = self.materialize_mols(context)?;
        if mols.is_empty() {
            return Ok(None);
        }
        crate::formula_check::check_effect_delta_formula(
            context,
            self.effect(),
            &mols,
            self.pattern_name(),
        );
        let products = mols
            .iter()
            .map(crate::mol::canon_smiles)
            .collect::<Vec<_>>();
        Ok(Some(Emission {
            site: self.site(),
            site_orbit: self.orbit().to_vec(),
            site_atoms: self.site_atoms(),
            cleaves: self.effect().cleaves,
            pattern_name: self.pattern_name().to_string(),
            search_bias: self.search_bias(),
            rule_path: self.rule_path().to_vec(),
            mols,
            products,
            cleave_side_sig: self.cleave_side_sig(),
            plan: self.plan(context),
        }))
    }

    pub fn is_pair_endpoint(&self) -> bool {
        match self {
            Self::Edit(e) => matches!(e.pattern.edit, Edit::PairEndpoint(_)),
            Self::Pair(_) => true,
        }
    }

    /// Elementary plan (edit identity, or pair leaf `canonical_plan`).
    pub fn identity_plan(&self, mol: &Molecule) -> Vec<Step> {
        match self {
            Self::Edit(e) => e.identity_plan(mol),
            Self::Pair(p) => pair_plan(p, mol),
        }
    }

    pub fn identity_plan_with_gens(
        &self,
        generators: &[crate::orbits::AtomBondGenerator],
        n_atoms: usize,
        mol: &Molecule,
    ) -> Vec<Step> {
        match self {
            Self::Edit(e) => e.identity_plan_with_gens(generators, n_atoms, mol),
            Self::Pair(p) => pair_plan(p, mol),
        }
    }
}

fn pair_plan(p: &PairCandidate, mol: &Molecule) -> Vec<Step> {
    let site_atoms = p.plan_site_atoms();
    let ends = [&p.left.effect, &p.right.effect];
    match p.leaf_rule().and_then(crate::rules::leaf_rule) {
        Some(leaf) => leaf.canonical_plan(mol, &site_atoms, Some(&ends)),
        None => {
            let rule = p.rule_name().to_string();
            crate::canonical_plan::identity_plan(rule, site_atoms)
        }
    }
}

impl EditCandidate {
    pub fn identity_plan(&self, mol: &Molecule) -> Vec<Step> {
        let gens = crate::orbits::atom_bond_generators(mol);
        self.identity_plan_with_gens(&gens, mol.atom_count(), mol)
    }

    pub fn identity_plan_with_gens(
        &self,
        generators: &[crate::orbits::AtomBondGenerator],
        n_atoms: usize,
        mol: &Molecule,
    ) -> Vec<Step> {
        let rule = self
            .rule_path
            .first()
            .and_then(|n| n.as_deref())
            .unwrap_or(self.pattern.name.as_str())
            .to_string();
        let mut site_atoms: Vec<usize> = self
            .pattern
            .site_map
            .iter()
            .filter_map(|m| self.mapped.get(m).copied())
            .collect();
        if site_atoms.is_empty() {
            site_atoms.push(self.site);
        }
        let orbit = crate::orbits::atom_orbit_with_gens(generators, n_atoms, self.site);
        identity_plan_at_indexes(
            rule.clone(),
            mol,
            site_atoms.iter().copied(),
            orbit.iter().copied(),
        )
        .unwrap_or_else(|_| {
            crate::canonical_plan::identity_plan_with_orbit(rule, site_atoms, orbit)
        })
    }
}
