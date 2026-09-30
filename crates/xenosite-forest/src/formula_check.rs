//! Check declared [`Effect::delta_formula`] against observed product formulas.
//!
//! Mismatches emit [`log::warn!`] (Python analogue: ``warnings.warn`` /
//! ``FormulaDeltaMismatchWarning``) and return a structured
//! [`FormulaDeltaMismatch`] so callers can append to
//! [`crate::find_path::PathCounters::formula_delta_mismatches`] and bump
//! [`crate::find_path::PathCounters::formula_delta_mismatch`]. Soft only —
//! never drops chemistry.
//!
//! Comparisons include hydrogen. Do not strip H to hide Effect/mol disagreement.

use std::collections::BTreeMap;

use crate::forest::{Formula, formula_delta, molecule_formula};
use crate::mol::Molecule;
use crate::pattern::{Effect, bag_delta_formula};

/// One soft formula-delta disagreement (declared vs observed, including H).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FormulaDeltaMismatch {
    pub pattern_name: String,
    pub declared: BTreeMap<String, i32>,
    pub observed: BTreeMap<String, i32>,
    pub cleaves: bool,
    pub adds: Option<String>,
    pub removes: Option<String>,
    pub leave: BTreeMap<String, i32>,
}

/// Zero-free element→delta map (includes H when nonzero).
fn nonzero(counts: &BTreeMap<String, i32>) -> BTreeMap<String, i32> {
    counts
        .iter()
        .filter(|(_, n)| **n != 0)
        .map(|(el, n)| (el.clone(), *n))
        .collect()
}

fn sum_formulas(mols: &[Molecule]) -> Formula {
    let mut counts = BTreeMap::new();
    let mut positive = 0i32;
    let mut negative = 0i32;
    for mol in mols {
        let f = molecule_formula(mol);
        for (el, n) in f.counts {
            if el == crate::forest::CHARGE_PLUS {
                positive += n;
            } else if el == crate::forest::CHARGE_MINUS {
                negative += n;
            } else {
                *counts.entry(el).or_insert(0) += n;
            }
        }
    }
    Formula { counts, charge: 0 }.with_charge_units(positive, negative)
}

/// Expected net for multi-fragment cleavage: junction bags only.
///
/// Named ``leave_formula`` atoms stay in a product fragment (cancel in
/// ``sum(products) − parent``). ``removes`` that are eliminated (halide, not
/// kept as a fragment) appear in the net — use ``adds − removes``, not adds
/// alone.
fn cleavage_net_expected(effect: &Effect) -> BTreeMap<String, i32> {
    nonzero(&bag_delta_formula(
        effect.adds.as_deref(),
        effect.removes.as_deref(),
    ))
}

fn cleavage_net_actual(parent: &Formula, products: &[Molecule]) -> BTreeMap<String, i32> {
    let after = sum_formulas(products);
    nonzero(&formula_delta(parent, &after).counts)
}

fn single_expected(effect: &Effect) -> BTreeMap<String, i32> {
    nonzero(&effect.resolved_delta_formula())
}

fn single_actual(parent: &Formula, product: &Molecule) -> BTreeMap<String, i32> {
    let after = molecule_formula(product);
    nonzero(&formula_delta(parent, &after).counts)
}

fn format_map(map: &BTreeMap<String, i32>) -> String {
    if map.is_empty() {
        return "{}".into();
    }
    let parts: Vec<String> = map.iter().map(|(k, v)| format!("{k}:{v}")).collect();
    format!("{{{}}}", parts.join(", "))
}

/// Compare observed product formula change to declared delta.
///
/// Returns ``None`` when they match (or the check was skipped). On mismatch
/// logs a warning and returns the structured record.
///
/// Check declared [`Effect::delta_formula`] against observed product formulas
/// (including H).
///
/// - One product, or several **non-cleaving** alternatives (e.g. ResonancePair
///   kekulé / path variants): each product is a full molecule; compare
///   ``product − parent`` to ``delta_formula`` per piece.
/// - Cleavage (`effect.cleaves` with 2+ products): ``sum(products) − parent``
///   vs junction ``adds − removes``.
pub fn check_effect_delta_formula(
    parent: &Molecule,
    effect: &Effect,
    products: &[Molecule],
    pattern_name: &str,
) -> Option<FormulaDeltaMismatch> {
    if products.is_empty() {
        return None;
    }
    let parent_f = molecule_formula(parent);

    // Multiple non-cleaving products are alternative full outcomes, not fragments.
    if !effect.cleaves {
        for product in products {
            if let Some(detail) = check_single_delta(&parent_f, effect, product, pattern_name) {
                return Some(detail);
            }
        }
        return None;
    }

    let (expected, actual) = if products.len() == 1 {
        (
            single_expected(effect),
            single_actual(&parent_f, &products[0]),
        )
    } else {
        (
            cleavage_net_expected(effect),
            cleavage_net_actual(&parent_f, products),
        )
    };

    finish_mismatch(effect, pattern_name, expected, actual, products.len())
}

fn check_single_delta(
    parent_f: &crate::forest::Formula,
    effect: &Effect,
    product: &Molecule,
    pattern_name: &str,
) -> Option<FormulaDeltaMismatch> {
    let expected = single_expected(effect);
    let actual = single_actual(parent_f, product);
    finish_mismatch(effect, pattern_name, expected, actual, 1)
}

fn finish_mismatch(
    effect: &Effect,
    pattern_name: &str,
    expected: BTreeMap<String, i32>,
    actual: BTreeMap<String, i32>,
    n_products: usize,
) -> Option<FormulaDeltaMismatch> {
    // Star conjugates use dummy ``*`` atoms — bag stoichiometry ≠ mol formula.
    if expected.contains_key("*") || actual.contains_key("*") {
        return None;
    }
    // Open leave (cleaves, empty leave_formula) with unexplained heavy loss:
    // annotation incomplete, not a sealed-delta bug.
    if effect.cleaves
        && effect.leave_formula.is_empty()
        && expected.is_empty()
        && actual.values().any(|&n| n < 0)
    {
        return None;
    }
    // Ring-retained leave: named leave still on the single product (e.g.
    // isoxazole N–O open). Incomplete leave vs fragment split — skip.
    if effect.cleaves && !effect.leave_formula.is_empty() && n_products == 1 && actual.is_empty() {
        let leave_as_delta: BTreeMap<String, i32> = effect
            .leave_formula
            .iter()
            .filter(|(_, n)| **n != 0)
            .map(|(el, n)| (el.clone(), -n))
            .collect();
        if expected == nonzero(&leave_as_delta) {
            return None;
        }
    }

    if expected == actual {
        return None;
    }
    log::warn!(
        "Formula delta mismatch for pattern {pattern_name}: declared \
         delta {} ≠ observed {} (cleaves={}, adds={:?}, removes={:?}, leave={:?})",
        format_map(&expected),
        format_map(&actual),
        effect.cleaves,
        effect.adds,
        effect.removes,
        effect.leave_formula,
    );
    Some(FormulaDeltaMismatch {
        pattern_name: pattern_name.to_string(),
        declared: expected,
        observed: actual,
        cleaves: effect.cleaves,
        adds: effect.adds.clone(),
        removes: effect.removes.clone(),
        leave: effect.leave_formula.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::pattern::{Effect, leave_me};

    #[test]
    fn hydroxyl_delta_matches_including_h() {
        // Ethane → ethanol: net +O; H count unchanged (C2H6 → C2H6O).
        let parent = parse_mol("CC").unwrap();
        let product = parse_mol("CCO").unwrap();
        let effect = Effect {
            adds: Some("O".into()),
            ..Effect::default()
        }
        .sealed();
        assert!(check_effect_delta_formula(&parent, &effect, &[product], "h").is_none());
    }

    #[test]
    fn hydroxyl_wrong_h_remove_mismatches() {
        // Old edit-stoichiometry bags (O:+1,H:-1) disagree with mol formula.
        let parent = parse_mol("CC").unwrap();
        let product = parse_mol("CCO").unwrap();
        let effect = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            ..Effect::default()
        }
        .sealed();
        let detail =
            check_effect_delta_formula(&parent, &effect, &[product], "h").expect("mismatch");
        assert_eq!(detail.declared.get("H"), Some(&-1));
        assert!(!detail.observed.contains_key("H"));
        assert_eq!(detail.observed.get("O"), Some(&1));
    }

    #[test]
    fn wrong_declared_oxygen_mismatches() {
        let parent = parse_mol("CC").unwrap();
        let product = parse_mol("CCO").unwrap();
        let effect = Effect {
            adds: Some("OO".into()),
            ..Effect::default()
        }
        .sealed();
        let detail =
            check_effect_delta_formula(&parent, &effect, &[product], "bad").expect("mismatch");
        assert_eq!(detail.pattern_name, "bad");
        assert_eq!(detail.declared.get("O"), Some(&2));
        assert_eq!(detail.observed.get("O"), Some(&1));
    }

    #[test]
    fn cleavage_net_matches_external_adds() {
        let parent = parse_mol("COc1ccccc1").unwrap();
        let phenol = parse_mol("Oc1ccccc1").unwrap();
        let formic = parse_mol("O=CO").unwrap();
        let effect = Effect {
            adds: Some("OO".into()),
            cleaves: true,
            leave_count: Some(1),
            leave_formula: leave_me(),
            ..Effect::default()
        }
        .sealed();
        assert!(
            check_effect_delta_formula(&parent, &effect, &[phenol, formic], "methyl_carboxylic")
                .is_none()
        );
    }

    #[test]
    fn noncleaving_multi_product_alternatives_checked_each() {
        // ResonancePair-style: several full-molecule alternatives, not fragments.
        let parent = parse_mol("O=C1CCCCC1").unwrap();
        let enol_a = parse_mol("OC1=CCCCC1").unwrap();
        let enol_b = parse_mol("OC1=CCCCC1").unwrap();
        let effect = Effect::default().sealed();
        assert!(
            check_effect_delta_formula(&parent, &effect, &[enol_a, enol_b], "tautomer").is_none(),
            "summing alternatives as cleavage would false-positive"
        );
        // One bad alternative still reports mismatch.
        let bad = parse_mol("CCO").unwrap();
        let enol = parse_mol("OC1=CCCCC1").unwrap();
        assert!(check_effect_delta_formula(&parent, &effect, &[enol, bad], "tautomer").is_some());
    }

    #[test]
    fn hydrogenation_h_delta_matches() {
        let parent = parse_mol("C=C").unwrap();
        let product = parse_mol("CC").unwrap();
        let effect = Effect {
            adds: Some("HH".into()),
            ..Effect::default()
        }
        .sealed();
        assert!(check_effect_delta_formula(&parent, &effect, &[product], "h").is_none());
    }
}

#[cfg(test)]
mod alkene_probe {
    use super::*;
    use crate::ForestMol;
    use crate::mol::canon_smiles;
    use crate::rules::hydrogenation;
    use std::collections::BTreeMap;

    #[test]
    fn alkene_patterns_do_not_duplicate_site_product() {
        // Two patterns must not match the same site to yield the same product.
        let substrates = [
            "C=C",
            "c1ccccc1",
            "C=C=O",
            "C=Cc1ccccc1",
            "CC=C",
            "c1ccc:C(:c1)=C=O",
        ];
        let mut dups = Vec::new();
        for smi in substrates {
            let Ok(parent) = ForestMol::parse(smi) else {
                continue;
            };
            let mut by_key: BTreeMap<(Vec<usize>, String), Vec<String>> = BTreeMap::new();
            for site in hydrogenation().candidates(&parent).filter_map(Result::ok) {
                if !site.pattern_name().starts_with("alkene") {
                    continue;
                }
                let Ok(pieces) = site.materialize_mols() else {
                    continue;
                };
                let mut atoms = site.site_atoms();
                atoms.sort_unstable();
                for p in pieces {
                    let key = (atoms.clone(), canon_smiles(&p));
                    by_key
                        .entry(key)
                        .or_default()
                        .push(site.pattern_name().to_string());
                }
            }
            for ((atoms, prod), names) in by_key {
                let mut u = names.clone();
                u.sort();
                u.dedup();
                if u.len() > 1 {
                    dups.push(format!(
                        "{smi} atoms={atoms:?} product={prod} patterns={u:?}"
                    ));
                }
            }
        }
        assert!(
            dups.is_empty(),
            "overlapping alkene patterns on same site/product:\n  {}",
            dups.join("\n  ")
        );
    }

    #[test]
    fn alkene_effect_bags_match_on_cumulene_ethene_benzene() {
        for smi in ["C=C=O", "C=C", "c1ccccc1"] {
            let parent = ForestMol::parse(smi).unwrap();
            let mut saw = false;
            for site in hydrogenation().candidates(&parent).filter_map(Result::ok) {
                if !site.pattern_name().starts_with("alkene") {
                    continue;
                }
                saw = true;
                let Ok(pieces) = site.materialize_mols() else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                assert!(
                    check_effect_delta_formula(
                        parent.mol(),
                        &site.effect,
                        &pieces,
                        site.pattern_name(),
                    )
                    .is_none(),
                    "bag mismatch on {smi} / {}",
                    site.pattern_name()
                );
            }
            assert!(saw, "expected alkene* on {smi}");
        }
    }
}
