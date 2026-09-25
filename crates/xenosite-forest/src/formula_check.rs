//! Check declared [`Effect::delta_formula`] against observed product formulas.
//!
//! Mismatches emit [`log::warn!`] (Python analogue: ``warnings.warn`` /
//! ``FormulaDeltaMismatchWarning``). The host app chooses a logger
//! (``env_logger``, etc.). Soft only — never drops chemistry.

use std::collections::BTreeMap;

use crate::forest::{Formula, formula_delta, molecule_formula};
use crate::mol::Molecule;
use crate::pattern::{Effect, bag_delta_formula};

/// Drop hydrogens — edit bags and sanitized mols disagree on H for hydroxyl.
fn heavy(counts: &BTreeMap<String, i32>) -> BTreeMap<String, i32> {
    counts
        .iter()
        .filter(|(el, n)| *el != "H" && **n != 0)
        .map(|(el, n)| (el.clone(), *n))
        .collect()
}

fn sum_formulas(mols: &[Molecule]) -> Formula {
    let mut counts = BTreeMap::new();
    let mut charge = 0i32;
    for mol in mols {
        let f = molecule_formula(mol);
        for (el, n) in f.counts {
            *counts.entry(el).or_insert(0) += n;
        }
        charge += f.charge;
    }
    Formula { counts, charge }
}

/// Expected heavy-atom net for multi-fragment cleavage: junction bags only.
///
/// Named ``leave_formula`` atoms stay in a product fragment (cancel in
/// ``sum(products) − parent``). ``removes`` that are eliminated (halide, not
/// kept as a fragment) appear in the net — use ``adds − removes``, not adds
/// alone.
fn cleavage_net_expected(effect: &Effect) -> BTreeMap<String, i32> {
    heavy(&bag_delta_formula(
        effect.adds.as_deref(),
        effect.removes.as_deref(),
    ))
}

fn cleavage_net_actual(parent: &Formula, products: &[Molecule]) -> BTreeMap<String, i32> {
    let after = sum_formulas(products);
    heavy(&formula_delta(parent, &after).counts)
}

fn single_expected(effect: &Effect) -> BTreeMap<String, i32> {
    heavy(&effect.resolved_delta_formula())
}

fn single_actual(parent: &Formula, product: &Molecule) -> BTreeMap<String, i32> {
    let after = molecule_formula(product);
    heavy(&formula_delta(parent, &after).counts)
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
/// Returns ``true`` when they match (or the check was skipped). On mismatch
/// logs a warning and returns ``false``.
///
/// Single product: heavy ``product − parent`` vs heavy ``delta_formula``.
/// Cleavage (2+ products): heavy ``sum(products) − parent`` vs ``adds − removes``.
pub fn check_effect_delta_formula(
    parent: &Molecule,
    effect: &Effect,
    products: &[Molecule],
    pattern_name: &str,
) -> bool {
    if products.is_empty() {
        return true;
    }
    let parent_f = molecule_formula(parent);

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

    // Star conjugates use dummy ``*`` atoms — bag stoichiometry ≠ mol formula.
    if expected.contains_key("*") || actual.contains_key("*") {
        return true;
    }
    // Open leave (cleaves, empty leave_formula) with unexplained heavy loss:
    // annotation incomplete, not a sealed-delta bug.
    if effect.cleaves
        && effect.leave_formula.is_empty()
        && expected.is_empty()
        && actual.values().any(|&n| n < 0)
    {
        return true;
    }

    if expected == actual {
        return true;
    }
    log::warn!(
        "Formula delta mismatch for pattern {pattern_name}: declared heavy \
         delta {} ≠ observed {} (cleaves={}, adds={:?}, removes={:?}, leave={:?})",
        format_map(&expected),
        format_map(&actual),
        effect.cleaves,
        effect.adds,
        effect.removes,
        effect.leave_formula,
    );
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::pattern::{Effect, leave_me};

    #[test]
    fn hydroxyl_heavy_delta_matches() {
        let parent = parse_mol("CC").unwrap();
        let product = parse_mol("CCO").unwrap();
        let effect = Effect {
            adds: Some("O".into()),
            removes: Some("H".into()),
            ..Effect::default()
        }
        .sealed();
        assert!(check_effect_delta_formula(
            &parent,
            &effect,
            &[product],
            "h"
        ));
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
        assert!(!check_effect_delta_formula(
            &parent,
            &effect,
            &[product],
            "bad"
        ));
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
        assert!(check_effect_delta_formula(
            &parent,
            &effect,
            &[phenol, formic],
            "methyl_carboxylic"
        ));
    }
}
