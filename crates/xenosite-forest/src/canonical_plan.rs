//! Canonical elementary plans for composite reaction hops.
//!
//! Mirrors Python `xenosite.forest.canonical_plan`. A search records these
//! steps on each walk; ordinary rules are identity (one step at the site).
//! Quinone-shaped hops expand via [`PlanKind`] data on the emitting
//! [`crate::ruleset::RuleSet`] — not a rule-name branch in `find_path`.

use crate::pattern::Effect;

/// One atom note in a [`CanonicalStep`] site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PlanAtom {
    /// Known index on the mol at emit time.
    Index(usize),
    /// Atom a prior elementary step will add (origin index + element).
    Ref {
        idx: usize,
        element: String,
        depth: u8,
    },
}

impl PlanAtom {
    pub fn index(idx: usize) -> Self {
        Self::Index(idx)
    }

    pub fn oxygen_ref(idx: usize) -> Self {
        Self::Ref {
            idx,
            element: "O".into(),
            depth: 0,
        }
    }

    /// Origin / index when this note points at a known or future atom.
    pub fn anchor(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::Ref { idx, .. } => Some(*idx),
        }
    }
}

/// One elementary step in a canonical plan.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalStep {
    pub rule: String,
    pub site: Vec<PlanAtom>,
}

impl CanonicalStep {
    pub fn new(rule: impl Into<String>, site: impl IntoIterator<Item = PlanAtom>) -> Self {
        Self {
            rule: rule.into(),
            site: site.into_iter().collect(),
        }
    }
}

/// How a leaf [`crate::ruleset::RuleSet`] expands an accepted hop into
/// elementary [`CanonicalStep`]s. Read by emit; not a search `if` on the name.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PlanKind {
    /// The rule is already elementary: one step at the discovery site.
    #[default]
    Identity,
    /// Prep missing oxygens (hydroxylation / oxidative dehalogenation), then
    /// one dehydrogenation — quinone-formation shape.
    HydroxylationThenDehydrogenation,
}

/// Identity plan: `rule` at the given site atoms.
pub fn identity_canonical_plan(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
) -> Vec<CanonicalStep> {
    vec![CanonicalStep::new(
        rule,
        site.into_iter().map(PlanAtom::index),
    )]
}

fn end_needs_oxygen(effect: &Effect) -> bool {
    let partner = effect.partner.as_deref().unwrap_or("");
    let needs_o = effect.adds.as_deref().is_some_and(|a| a.contains('O'));
    // Python also reads end["needs"]; Rust Effect has adds/partner only.
    needs_o && partner != "O"
}

fn bonded(mol: &crate::Molecule, idx: usize, atomic_num: u8) -> Option<usize> {
    use crate::mol::atom_idx;
    mol.neighbors(atom_idx(idx)).find_map(|(nbr, _)| {
        let n = crate::mol::atom_usize(nbr);
        (mol.atom(nbr).element.atomic_number() == atomic_num).then_some(n)
    })
}

const HALOGEN: &[&str] = &["F", "Cl", "Br", "I", "At"];

fn halogen_z(partner: &str) -> Option<u8> {
    match partner {
        "F" => Some(9),
        "Cl" => Some(17),
        "Br" => Some(35),
        "I" => Some(53),
        "At" => Some(85),
        _ => None,
    }
}

/// Preps that supply missing oxygens, then one dehydrogenation.
pub fn hydroxylation_then_dehydrogenation(
    mol: &crate::Molecule,
    ends: &[&Effect],
    end_atoms: &[usize],
) -> Vec<CanonicalStep> {
    let mut preps = Vec::new();
    let mut dh_refs = Vec::new();
    for (end, &atom) in ends.iter().zip(end_atoms.iter()) {
        let partner = end.partner.as_deref().unwrap_or("");
        if end_needs_oxygen(end) {
            let anchor = PlanAtom::index(atom);
            if HALOGEN.contains(&partner) {
                let Some(z) = halogen_z(partner) else {
                    continue;
                };
                let Some(halo) = bonded(mol, atom, z) else {
                    continue;
                };
                preps.push(CanonicalStep::new(
                    "OxidativeDehalogenation",
                    [anchor.clone(), PlanAtom::index(halo)],
                ));
            } else {
                preps.push(CanonicalStep::new("Hydroxylation", [anchor.clone()]));
            }
            dh_refs.push(PlanAtom::oxygen_ref(atom));
            continue;
        }
        let atomic_num = match partner {
            "O" => Some(8),
            "N" => Some(7),
            "C" => Some(6),
            "S" => Some(16),
            _ => None,
        };
        if let Some(z) = atomic_num {
            if let Some(hetero) = bonded(mol, atom, z) {
                dh_refs.push(PlanAtom::index(hetero));
            }
        }
    }
    if dh_refs.is_empty() {
        return Vec::new();
    }
    preps.push(CanonicalStep::new("Dehydrogenation", dh_refs));
    preps
}

/// Build elementary steps from [`PlanKind`] and emit-time site data.
pub fn steps_for_kind(
    kind: PlanKind,
    rule_name: &str,
    mol: &crate::Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<CanonicalStep> {
    match kind {
        PlanKind::Identity => identity_canonical_plan(rule_name, site_atoms.iter().copied()),
        PlanKind::HydroxylationThenDehydrogenation => {
            if let (Some(ends), true) = (end_effects, site_atoms.len() >= 2) {
                let plan = hydroxylation_then_dehydrogenation(mol, ends, site_atoms);
                if !plan.is_empty() {
                    return plan;
                }
            }
            identity_canonical_plan("Dehydrogenation", site_atoms.iter().copied())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::pattern::Effect;

    #[test]
    fn identity_is_one_step() {
        let plan = identity_canonical_plan("Hydroxylation", [0]);
        assert_eq!(plan.len(), 1);
        assert_eq!(plan[0].rule, "Hydroxylation");
        assert_eq!(plan[0].site, vec![PlanAtom::index(0)]);
    }

    #[test]
    fn hydroquinone_pair_needs_no_prep() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        // End atoms are the ring carbons (map 1); partner O is found by bond.
        let carbons: Vec<usize> = (0..mol.atom_count())
            .filter(|&i| {
                let a = crate::mol::atom_idx(i);
                mol.atom(a).element.atomic_number() == 6
                    && mol
                        .neighbors(a)
                        .any(|(n, _)| mol.atom(n).element.atomic_number() == 8)
            })
            .collect();
        assert_eq!(carbons.len(), 2);
        let phenol = Effect {
            adds: None,
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: Some("O".into()),
        };
        let ends = [&phenol, &phenol];
        let plan = hydroxylation_then_dehydrogenation(&mol, &ends, &carbons);
        // Partner O already bonded → DH only, pointing at the oxygens.
        assert_eq!(plan.len(), 1);
        assert_eq!(plan[0].rule, "Dehydrogenation");
        assert_eq!(plan[0].site.len(), 2);
    }

    #[test]
    fn bare_carbon_end_preps_hydroxylation() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let end = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            cleaves: false,
            leave_count: None,
            methide: false,
            dearomatizes: true,
            partner: None,
        };
        let ends = [&end, &end];
        let plan = hydroxylation_then_dehydrogenation(&mol, &ends, &[0, 3]);
        assert_eq!(plan.len(), 3);
        assert_eq!(plan[0].rule, "Hydroxylation");
        assert_eq!(plan[1].rule, "Hydroxylation");
        assert_eq!(plan[2].rule, "Dehydrogenation");
        assert!(matches!(plan[2].site[0], PlanAtom::Ref { .. }));
    }
}
