//! Always-on PatternInfo catalog audit (part of `cargo test --lib`).
//!
//! Structural invariants + chemistry probe: if applying a pattern clears
//! aromaticity at `site_map` atoms, the catalog must declare
//! `Effect.dearomatizes` capability (resolved false on aliphatic matches).

use std::collections::{BTreeMap, BTreeSet};

use crate::as_forest_mol;
use crate::ForestMol;
use crate::mol::{atom_idx, parse_mol};
use crate::pattern::{Edit, Effect, PatternInfo, SiteKind, compose_delta_formula};
use crate::rules::{catalog_names, leaf_rule};

/// Probe mols that expose aromatic (and a few aliphatic) sites.
const PROBES: &[&str] = &[
    "c1ccccc1",
    "COc1ccccc1",
    "COc1ccc(O)cc1",
    "Oc1ccc(O)cc1",
    "c1ccc2ccccc2c1",
    "COc1ccc2ccccc2c1",
    "C=C",
    "CC=O",
    "C1OC1",
];

/// Conjugation adducts — not publicly published; Effect bags still unfinished.
/// Tracked as xfail: [`catalog_adduct_effect_and_atom_diff_match_materialized_products`].
const ADDUCT_LEAVES: &[&str] = &[
    "Acetylation",
    "Sulfation",
    "Glucuronidation",
    "Glutathionation",
];

fn expected_delta(effect: &Effect) -> BTreeMap<String, i32> {
    compose_delta_formula(
        effect.adds.as_deref(),
        effect.removes.as_deref(),
        &effect.leave_formula,
    )
}

fn assert_effect_sealed(leaf: &str, pattern: &str, effect: &Effect) {
    let want = expected_delta(effect);
    assert_eq!(
        effect.delta_formula, want,
        "{leaf}/{pattern}: delta_formula not sealed from adds/removes/leave; got {:?} want {:?}",
        effect.delta_formula, want
    );
}

fn assert_site_map_matches_kind(leaf: &str, info: &PatternInfo) {
    assert!(
        !info.site_map.is_empty(),
        "{leaf}/{}: empty site_map",
        info.name
    );
    let primary = info.primary_map();
    assert!(
        info.site_map.contains(&primary),
        "{leaf}/{}: primary_map {primary} not in site_map {:?}",
        info.name,
        info.site_map
    );
    let n = info.site_map.len();
    match info.site_kind {
        // Atom primary; extra maps allowed (e.g. Dehydration beta-elim adjacent C).
        SiteKind::Atom => assert!(n >= 1, "{leaf}/{}: Atom site_map empty", info.name),
        SiteKind::Bond | SiteKind::DirectedBond => assert_eq!(
            n, 2,
            "{leaf}/{}: {:?} site_map must be length 2, got {:?}",
            info.name, info.site_kind, info.site_map
        ),
        // ResonancePair ends are one atom; one-bond SMARTS may list both.
        SiteKind::AtomPair => assert!(
            n == 1 || n == 2,
            "{leaf}/{}: AtomPair site_map must be length 1 or 2, got {:?}",
            info.name,
            info.site_map
        ),
    }
}

fn assert_edit_present(leaf: &str, info: &PatternInfo) {
    match &info.edit {
        Edit::Smirks(s) => assert!(
            s.contains(">>"),
            "{leaf}/{}: Smirks missing >>: {s}",
            info.name
        ),
        Edit::Hydroxyl => {}
        Edit::PairEndpoint(s) => assert!(
            !s.is_empty(),
            "{leaf}/{}: empty PairEndpoint token",
            info.name
        ),
    }
    assert!(
        !info.smarts.is_empty(),
        "{leaf}/{}: empty smarts",
        info.name
    );
}

#[test]
fn catalog_pattern_info_structural() {
    let mut seen_leaves = BTreeSet::new();
    for &name in catalog_names() {
        assert!(seen_leaves.insert(name), "duplicate catalog name {name}");
        let set = leaf_rule(name).unwrap_or_else(|| panic!("leaf_rule({name})"));
        assert_eq!(set.name.as_deref(), Some(name));
        let patterns = set.patterns();
        assert!(
            !patterns.is_empty(),
            "{name}: leaf must declare at least one PatternInfo"
        );
        let mut names = BTreeSet::new();
        for info in &patterns {
            assert!(
                !info.name.is_empty(),
                "{name}: PatternInfo.name must be non-empty"
            );
            assert!(
                names.insert(info.name.as_str()),
                "{name}: duplicate PatternInfo.name {:?}",
                info.name
            );
            assert_site_map_matches_kind(name, info);
            assert_edit_present(name, info);
            assert_effect_sealed(name, &info.name, &info.effect);
            for (i, arm) in info.possibilities.iter().enumerate() {
                assert_effect_sealed(name, &format!("{}#poss{i}", info.name), arm);
            }
        }
        if name == "EpoxideHydration" {
            assert_eq!(patterns.len(), 1);
            assert_eq!(patterns[0].name, "diol");
            assert_eq!(patterns[0].site_kind, SiteKind::Bond);
            assert_eq!(patterns[0].effect.adds.as_deref(), Some("OOHH"));
            assert_eq!(patterns[0].effect.delta_formula.get("O"), Some(&2));
            assert_eq!(patterns[0].effect.delta_formula.get("H"), Some(&2));
            assert!(patterns[0].effect.dearomatizes);
        }
    }
    assert_eq!(
        seen_leaves.len(),
        catalog_names().len(),
        "catalog walk missed leaves"
    );
}

/// If a catalog pattern clears aromaticity at site_map atoms on a probe,
/// `Effect.dearomatizes` capability must be true on the catalog record.
#[test]
fn catalog_dearomatizes_capability_matches_chemistry() {
    let mut misses: Vec<String> = Vec::new();
    for &name in catalog_names() {
        let set = leaf_rule(name).expect(name);
        let catalog: BTreeMap<&str, &PatternInfo> = set
            .patterns()
            .into_iter()
            .map(|p| (p.name.as_str(), p))
            .collect();

        for &smi in PROBES {
            let Ok(parent) = ForestMol::parse(smi) else {
                continue;
            };
            let Ok(cands) = set.candidates(&parent).collect::<Result<Vec<_>, _>>() else {
                continue;
            };
            for c in cands {
                let Some(catalog_info) = catalog.get(c.pattern.name.as_str()) else {
                    continue;
                };
                if catalog_info.effect.dearomatizes {
                    continue;
                }
                // Cleavage can destroy aromatic sites without being a
                // dearomatizing pathway — capability is for kept-site ring loss.
                if catalog_info.effect.cleaves {
                    continue;
                }
                let site_atoms: Vec<usize> = catalog_info
                    .site_map
                    .iter()
                    .filter_map(|m| c.mapped.get(m).copied())
                    .collect();
                if site_atoms.is_empty() {
                    continue;
                }
                let aromatic_before: Vec<bool> = site_atoms
                    .iter()
                    .map(|&i| parent.mol().atom(atom_idx(i)).aromatic)
                    .collect();
                if !aromatic_before.iter().any(|&a| a) {
                    continue;
                }
                let Ok(pieces) = c.materialize_mols() else {
                    continue;
                };
                let Some(product) = pieces.first() else {
                    continue;
                };
                let child = parent.from_edit_product(product.clone());
                let lost = site_atoms.iter().enumerate().any(|(k, &r)| {
                    if !aromatic_before[k] {
                        return false;
                    }
                    let Some(tag) = parent.tag_of(r) else {
                        return false;
                    };
                    // Atom must still be present; gone/cleaved is not dearomatize.
                    let Some(j) = child.index_of(tag) else {
                        return false;
                    };
                    !child.mol().atom(atom_idx(j)).aromatic
                });
                if lost {
                    misses.push(format!(
                        "{name}/{} on {smi} site={:?}: clears aromaticity but catalog dearomatizes=false",
                        c.pattern.name, site_atoms
                    ));
                }
            }
        }
    }

    assert!(
        misses.is_empty(),
        "PatternInfo dearomatizes capability missing for chemistry that clears aromaticity:\n  {}",
        misses.join("\n  ")
    );
}

/// Hits + miss messages for Effect / atom_diff accuracy on `leaves` × PROBES.
fn effect_atom_diff_accuracy(leaves: &[&str]) -> (usize, Vec<String>) {
    use crate::atom_diff::{effect_adds_oxygen, effect_removes_h};
    use crate::forest::{formula_delta, molecule_formula};
    use crate::formula_check::check_effect_delta_formula;

    fn effect_adds_h(effect: &Effect) -> bool {
        effect.adds.as_deref().is_some_and(|a| a.contains('H'))
    }

    let mut misses: Vec<String> = Vec::new();
    let mut hits = 0usize;
    for &name in leaves {
        let set = leaf_rule(name).expect(name);
        for &smi in PROBES {
            let Ok(parent) = ForestMol::parse(smi) else {
                continue;
            };
            let Ok(cands) = set.candidates(&parent).collect::<Result<Vec<_>, _>>() else {
                continue;
            };
            let parent_f = molecule_formula(parent.mol());
            for c in cands {
                let Ok(pieces) = c.materialize_mols() else {
                    continue;
                };
                if pieces.is_empty() {
                    continue;
                }
                hits += 1;
                if let Some(detail) = check_effect_delta_formula(
                    parent.mol(),
                    &c.effect,
                    &pieces,
                    &c.pattern_name,
                ) {
                    misses.push(format!(
                        "{name}/{} on {smi}: formula mismatch declared {:?} observed {:?}",
                        c.pattern_name, detail.declared, detail.observed
                    ));
                    continue;
                }
                // Bag flags vs observed net (cleavage: sum products − parent).
                let observed = if c.effect.cleaves && pieces.len() > 1 {
                    let mut counts = BTreeMap::new();
                    for piece in &pieces {
                        for (el, n) in molecule_formula(piece).counts {
                            *counts.entry(el).or_insert(0) += n;
                        }
                    }
                    formula_delta(
                        &parent_f,
                        &crate::forest::Formula {
                            counts,
                            charge: 0,
                        },
                    )
                    .counts
                } else {
                    formula_delta(&parent_f, &molecule_formula(&pieces[0])).counts
                };
                let dh = observed.get("H").copied().unwrap_or(0);
                let d_o = observed.get("O").copied().unwrap_or(0);
                if dh > 0 && !effect_adds_h(&c.effect) && !c.effect.cleaves {
                    // Cleavage nets can gain H from caps without an adds bag when
                    // leave bookkeeping cancels; non-cleaving must declare HH.
                    misses.push(format!(
                        "{name}/{} on {smi}: observed H:+{dh} but Effect adds has no H ({:?})",
                        c.pattern_name, c.effect.adds
                    ));
                }
                if dh < 0 && !effect_removes_h(&c.effect) && !c.effect.cleaves {
                    misses.push(format!(
                        "{name}/{} on {smi}: observed H:{dh} but Effect removes has no H ({:?})",
                        c.pattern_name, c.effect.removes
                    ));
                }
                if d_o > 0 && !effect_adds_oxygen(&c.effect) && !c.effect.cleaves {
                    misses.push(format!(
                        "{name}/{} on {smi}: observed O:+{d_o} but Effect does not add oxygen",
                        c.pattern_name
                    ));
                }
            }
        }
    }
    (hits, misses)
}

#[test]
fn catalog_effect_and_atom_diff_match_materialized_products() {
    //! Every non-adduct LEAF_CTORS hit on PROBES: declared Effect (incl. H)
    //! matches the product formula delta, and atom_diff bag flags agree.
    //! When OR arms disagree on H/O, annotate each arm with [`When`] — do not
    //! hide disagreement by dropping H from the check.
    //! Adducts: see xfail [`catalog_adduct_effect_and_atom_diff_match_materialized_products`].
    let leaves: Vec<&str> = catalog_names()
        .iter()
        .copied()
        .filter(|n| !ADDUCT_LEAVES.contains(n))
        .collect();
    let (hits, misses) = effect_atom_diff_accuracy(&leaves);
    assert!(
        hits > 0,
        "expected at least one materialize hit across catalog probes"
    );
    assert!(
        misses.is_empty(),
        "PatternInfo Effect / atom_diff bag flags disagree with materialized products (incl. H).\n\
         When OR arms disagree, use When possibilities — do not drop H from the check.\n  {}",
        misses.join("\n  ")
    );
}

/// Same Effect/atom_diff accuracy gate as the catalog test, for conjugation
/// adducts only. `#[ignore]` = xfail until adduct Effect bags are refactored
/// (not published yet). Run: `cargo test -p xenosite-forest catalog_adduct -- --ignored`
#[test]
#[ignore = "xfail: adduct Effect bags unfinished / not published; remove ignore when refactor seals formula Effects"]
fn catalog_adduct_effect_and_atom_diff_match_materialized_products() {
    let (hits, misses) = effect_atom_diff_accuracy(ADDUCT_LEAVES);
    assert!(
        hits > 0,
        "expected at least one adduct materialize hit on probes"
    );
    assert!(
        misses.is_empty(),
        "adduct Effect / atom_diff bag flags disagree with products (incl. H).\n  {}",
        misses.join("\n  ")
    );
}

#[test]
fn catalog_resolve_dearomatizes_on_aromatic_probes() {
    // Capability true + aromatic site_map → resolved true; aliphatic → false.
    for &name in catalog_names() {
        let set = leaf_rule(name).expect(name);
        let capable: Vec<&str> = set
            .patterns()
            .into_iter()
            .filter(|p| p.effect.dearomatizes)
            .map(|p| p.name.as_str())
            .collect();
        if capable.is_empty() {
            continue;
        }
        for &smi in PROBES {
            let Ok(mol) = as_forest_mol(smi) else {
                continue;
            };
            let Ok(cands) = set.candidates(&mol).collect::<Result<Vec<_>, _>>() else {
                continue;
            };
            for c in cands {
                if !capable.contains(&c.pattern.name.as_str()) {
                    continue;
                }
                let site_aromatic = c
                    .pattern
                    .site_map
                    .iter()
                    .filter_map(|m| c.mapped.get(m).copied())
                    .any(|i| mol.mol().atom(atom_idx(i)).aromatic);
                // Candidate carries resolved effect (context mol).
                assert_eq!(
                    c.pattern.effect.dearomatizes, site_aromatic,
                    "{name}/{} on {smi}: resolved dearomatizes={} but site_aromatic={site_aromatic}",
                    c.pattern.name, c.pattern.effect.dearomatizes
                );
            }
        }
    }
}
