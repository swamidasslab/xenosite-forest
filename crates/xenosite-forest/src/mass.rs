//! Monoisotopic masses for formula → m/z (MS1).
//!
//! Exact-mass tables for common bioorganic elements. Used by
//! [`crate::find_path_ms1`] and shared tests that pin PatternInfo
//! `delta_formula` to a mass delta (drift guard vs structure `find_path`).
//!
//! Unlabeled atoms use the most abundant isotope ([`chematic::core::Element::atomic_mass`]).
//! Explicit `atom.isotope` labels use that nuclide's exact mass when known.

use std::collections::BTreeMap;

use chematic::core::Atom;

use crate::forest::Formula;
use crate::mol::Molecule;

/// Most abundant isotope mass (Da). Sources: IUPAC / common proteomics tables.
/// Kept in sync with [`chematic::core::Element::atomic_mass`] for the organic subset.
pub fn element_mono_mass(symbol: &str) -> Option<f64> {
    Some(match symbol {
        "H" => 1.007_825_032_24,
        "C" => 12.0,
        "N" => 14.003_074_004_43,
        "O" => 15.994_914_619_57,
        "F" => 18.998_403_162_73,
        "Na" => 22.989_769_282_0,
        "P" => 30.973_761_998_42,
        "S" => 31.972_071_174_4,
        "Cl" => 34.968_852_682, // ³⁵Cl
        "K" => 38.963_706_486_4,
        "Br" => 78.918_337_6, // ⁷⁹Br
        "I" => 126.904_471_9,
        _ => return None,
    })
}

/// Exact mass of a labeled nuclide (mass number), when known.
///
/// Unlabeled atoms should call [`element_mono_mass`] / `Element::atomic_mass`
/// instead — those are the common (most abundant) isotope.
pub fn isotope_exact_mass(symbol: &str, mass_number: u16) -> Option<f64> {
    Some(match (symbol, mass_number) {
        ("H", 1) => 1.007_825_032_24,
        ("H", 2) => 2.014_101_778_11,
        ("H", 3) => 3.016_049_277_9,
        ("C", 12) => 12.0,
        ("C", 13) => 13.003_354_835_07,
        ("C", 14) => 14.003_241_988_4,
        ("N", 14) => 14.003_074_004_43,
        ("N", 15) => 15.000_108_898_88,
        ("O", 16) => 15.994_914_619_57,
        ("O", 17) => 16.999_131_756_50,
        ("O", 18) => 17.999_159_612_86,
        ("S", 32) => 31.972_071_174_4,
        ("S", 33) => 32.971_458_909_8,
        ("S", 34) => 33.967_867_004,
        ("Cl", 35) => 34.968_852_682,
        ("Cl", 37) => 36.965_902_602,
        ("Br", 79) => 78.918_337_6,
        ("Br", 81) => 80.916_289_7,
        _ => return None,
    })
}

/// Proton mass for [M+H]⁺ / [M−H]⁻ adduct arithmetic (Da).
pub const PROTON_MASS: f64 = 1.007_276_466_578_9;

fn atom_mono_mass(atom: &Atom) -> Option<f64> {
    let sym = atom.element.symbol();
    match atom.isotope {
        Some(n) => isotope_exact_mass(sym, n),
        None => Some(atom.element.atomic_mass()),
    }
}

/// Monoisotopic neutral mass of a molecule.
///
/// Walks atoms: unlabeled → most abundant isotope; labeled → that nuclide when
/// tabulated. Implicit hydrogens are unlabeled common ¹H. Returns `None` when
/// any atom (or its label) lacks a mass entry. Charge is ignored (adducts via
/// [`mz_of`] / [`Ms1Adduct`]).
pub fn molecule_mono_mass(mol: &Molecule) -> Option<f64> {
    let mut mass = 0.0;
    for (idx, atom) in mol.atoms() {
        mass += atom_mono_mass(atom)?;
        let h = mol.implicit_hydrogen_count(idx);
        if h > 0 {
            mass += element_mono_mass("H")? * f64::from(h);
        }
    }
    Some(mass)
}

/// Monoisotopic neutral mass of a formula (sum of element mono masses × counts).
///
/// Element symbols only — no per-atom isotope labels. Prefer
/// [`molecule_mono_mass`] when the mol may carry labels. Returns `None` when
/// any element lacks a table entry. Charge is ignored (adducts by [`mz_of`]).
pub fn formula_mono_mass(formula: &Formula) -> Option<f64> {
    let mut mass = 0.0;
    for (el, &n) in &formula.counts {
        if n == 0 {
            continue;
        }
        if n < 0 {
            return None;
        }
        let m = element_mono_mass(el)?;
        mass += m * f64::from(n);
    }
    Some(mass)
}

/// Apply an element delta map onto a formula (zeros dropped).
pub fn formula_apply_delta(base: &Formula, delta: &BTreeMap<String, i32>) -> Formula {
    let mut counts = base.counts.clone();
    for (el, &d) in delta {
        let n = counts.get(el).copied().unwrap_or(0) + d;
        if n == 0 {
            counts.remove(el);
        } else {
            counts.insert(el.clone(), n);
        }
    }
    Formula {
        counts,
        charge: base.charge,
    }
}

/// How the observed m/z relates to the neutral monoisotopic mass.
///
/// MS1 search defaults to [`Self::MPlusH`] (positive mode).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Ms1Adduct {
    /// Neutral molecule mass (tests / debugging).
    Neutral,
    /// [M+H]⁺ — default MS1 positive mode.
    #[default]
    MPlusH,
    /// [M−H]⁻.
    MMinusH,
}

impl Ms1Adduct {
    /// Convert neutral mono mass → expected m/z for this adduct.
    pub fn mz_from_neutral(self, neutral: f64) -> f64 {
        match self {
            Self::Neutral => neutral,
            Self::MPlusH => neutral + PROTON_MASS,
            Self::MMinusH => neutral - PROTON_MASS,
        }
    }
}

/// Expected m/z for a formula under `adduct`.
pub fn mz_of(formula: &Formula, adduct: Ms1Adduct) -> Option<f64> {
    Some(adduct.mz_from_neutral(formula_mono_mass(formula)?))
}

/// Expected m/z for a molecule under `adduct` (isotope-aware).
pub fn mz_of_mol(mol: &Molecule, adduct: Ms1Adduct) -> Option<f64> {
    Some(adduct.mz_from_neutral(molecule_mono_mass(mol)?))
}

/// Absolute |observed − target| in Da.
pub fn mz_abs_error(observed: f64, target: f64) -> f64 {
    (observed - target).abs()
}

/// True when |observed − target| ≤ `tol_da`.
pub fn mz_within(observed: f64, target: f64, tol_da: f64) -> bool {
    mz_abs_error(observed, target) <= tol_da
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    use chematic::core::AtomIdx;

    use crate::forest::{formula_delta, molecule_formula};
    use crate::mol::parse_mol;
    use crate::pattern::{bag_delta_formula, compose_delta_formula};
    use crate::rules::{epoxide_hydration, epoxide_opening, hydroxylation, sulfur_oxidation};

    #[test]
    fn water_mono_mass() {
        let mol = parse_mol("O").unwrap();
        let f = molecule_formula(&mol);
        let from_formula = formula_mono_mass(&f).unwrap();
        let from_mol = molecule_mono_mass(&mol).unwrap();
        // H2O ≈ 18.0106
        assert!(
            (from_formula - 18.010_564_684).abs() < 1e-6,
            "{from_formula}"
        );
        assert!(
            (from_mol - from_formula).abs() < 1e-6,
            "{from_mol} vs {from_formula}"
        );
    }

    #[test]
    fn ethane_ethanol_mass_delta_matches_observed_plus_oxygen() {
        // Sanitized mols: ethane C2H6 → ethanol C2H6O is +O only.
        // Catalog hydroxyl bags still say +O −H (junction story); mass must
        // follow the observed formula (HEURISTICS / formula_check heavy-only).
        let ethane = molecule_formula(&parse_mol("CC").unwrap());
        let ethanol = molecule_formula(&parse_mol("CCO").unwrap());
        let observed = formula_delta(&ethane, &ethanol);
        assert_eq!(observed.counts.get("O"), Some(&1));
        assert_eq!(observed.counts.get("H"), None);
        let m0 = formula_mono_mass(&ethane).unwrap();
        let m1 = formula_mono_mass(&ethanol).unwrap();
        let from_obs = formula_mono_mass(&formula_apply_delta(&ethane, &observed.counts)).unwrap();
        assert!(
            (from_obs - m1).abs() < 1e-6,
            "pred={from_obs} obs={m1} from {m0}"
        );
    }

    #[test]
    fn hydroxylation_catalog_heavy_delta_matches_observed_oxygen() {
        let set = hydroxylation();
        let info = &set.patterns()[0];
        let declared = info.effect.resolved_delta_formula();
        assert_eq!(declared.get("O"), Some(&1));
        // H bag may be −1 while live mols keep H count (valence).
        let ethane = molecule_formula(&parse_mol("CC").unwrap());
        let ethanol = molecule_formula(&parse_mol("CCO").unwrap());
        let observed = formula_delta(&ethane, &ethanol);
        assert_eq!(declared.get("O"), observed.counts.get("O"));
        let mz_target = mz_of(&ethanol, Ms1Adduct::MPlusH).unwrap();
        // Mass from observed delta (not raw bag H) hits ethanol mz.
        let pred = formula_apply_delta(&ethane, &observed.counts);
        let mz_pred = mz_of(&pred, Ms1Adduct::MPlusH).unwrap();
        assert!(mz_within(mz_pred, mz_target, 1e-6));
    }

    #[test]
    fn epoxide_hydration_declared_delta_matches_mol_mass() {
        // Fixed bags: OOHH ≡ ethene → glycol (+O2 +H2).
        let set = epoxide_hydration();
        let info = set.patterns()[0];
        let declared = info.effect.resolved_delta_formula();
        assert_eq!(declared.get("O"), Some(&2));
        assert_eq!(declared.get("H"), Some(&2));
        let ethene = parse_mol("C=C").unwrap();
        let parent_f = molecule_formula(&ethene);
        let cand = set.candidates(&ethene).next().unwrap().unwrap();
        let product = &cand.materialize_mols(&ethene).unwrap()[0];
        let observed = formula_delta(&parent_f, &molecule_formula(product)).counts;
        assert_eq!(declared, observed);
        let m0 = molecule_mono_mass(&ethene).unwrap();
        let m1 = molecule_mono_mass(product).unwrap();
        let pred = formula_mono_mass(&formula_apply_delta(&parent_f, &declared)).unwrap();
        assert!((pred - m1).abs() < 1e-4, "pred={pred} obs={m1} from {m0}");
    }

    #[test]
    fn epoxide_opening_hydrate_declared_matches_mol() {
        let set = epoxide_opening();
        let hydrate = set
            .patterns()
            .into_iter()
            .find(|p| p.name == "hydrate")
            .unwrap();
        let declared = hydrate.effect.resolved_delta_formula();
        assert_eq!(declared.get("O"), Some(&1));
        assert_eq!(declared.get("H"), Some(&2));
        let epox = parse_mol("C1OC1").unwrap();
        let parent_f = molecule_formula(&epox);
        let cand = set
            .candidates(&epox)
            .map(|c| c.unwrap())
            .find(|c| c.pattern_name() == "hydrate")
            .unwrap();
        let product = &cand.materialize_mols(&epox).unwrap()[0];
        let observed = formula_delta(&parent_f, &molecule_formula(product)).counts;
        assert_eq!(declared, observed);
    }

    #[test]
    fn sulfur_hydroxy_h_is_substrate_dependent_keep_o_only() {
        // CSC (thioether) may show +H; CCS (thiol) may not. Bags stay O-only.
        let set = sulfur_oxidation();
        let hydroxy = set
            .patterns()
            .into_iter()
            .find(|p| p.name == "hydroxy")
            .unwrap();
        let declared = hydroxy.effect.resolved_delta_formula();
        assert_eq!(declared.get("O"), Some(&1));
        assert_eq!(declared.get("H"), None);
        let mol = parse_mol("CSC").unwrap();
        let parent_f = molecule_formula(&mol);
        let cand = set
            .candidates(&mol)
            .map(|c| c.unwrap())
            .find(|c| c.pattern_name() == "hydroxy")
            .unwrap();
        let product = &cand.materialize_mols(&mol).unwrap()[0];
        let observed = formula_delta(&parent_f, &molecule_formula(product)).counts;
        assert_eq!(observed.get("O"), Some(&1));
        // H may disagree; heavy O is what formula_check seals on.
    }

    #[test]
    fn compose_delta_matches_bag_for_simple_oh() {
        let a = bag_delta_formula(Some("O"), Some("H"));
        let b = compose_delta_formula(Some("O"), Some("H"), &BTreeMap::new());
        assert_eq!(a, b);
    }

    #[test]
    fn m_plus_h_adduct_shifts_by_proton() {
        let f = molecule_formula(&parse_mol("CC").unwrap());
        let neutral = formula_mono_mass(&f).unwrap();
        let mz = mz_of(&f, Ms1Adduct::MPlusH).unwrap();
        assert!((mz - neutral - PROTON_MASS).abs() < 1e-12);
    }

    #[test]
    fn labeled_carbon_13_shifts_mono_mass() {
        let mut mol = parse_mol("C").unwrap(); // methane
        let unlabeled = molecule_mono_mass(&mol).unwrap();
        mol.set_isotope(AtomIdx(0), Some(13));
        let labeled = molecule_mono_mass(&mol).unwrap();
        let delta = isotope_exact_mass("C", 13).unwrap() - element_mono_mass("C").unwrap();
        assert!(
            (labeled - unlabeled - delta).abs() < 1e-9,
            "unlabeled={unlabeled} labeled={labeled} want_delta={delta}"
        );
        // Formula path cannot see the label — stays at common mass.
        let from_formula = formula_mono_mass(&molecule_formula(&mol)).unwrap();
        assert!((from_formula - unlabeled).abs() < 1e-9);
    }

    #[test]
    fn mz_of_mol_matches_formula_when_unlabeled() {
        let mol = parse_mol("CCO").unwrap();
        let a = mz_of_mol(&mol, Ms1Adduct::MPlusH).unwrap();
        let b = mz_of(&molecule_formula(&mol), Ms1Adduct::MPlusH).unwrap();
        assert!((a - b).abs() < 1e-4, "{a} vs {b}");
    }
}
