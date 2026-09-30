//! Always-on PatternInfo catalog audit (part of `cargo test --lib`).
//!
//! Structural invariants + chemistry probe: if applying a pattern clears
//! aromaticity at `site_map` atoms, the catalog must declare
//! `Effect.dearomatizes` capability (resolved false on aliphatic matches).
//!
//! Effect / When coverage scans the shared native-parity pool
//! [`crate::substrate_library::coverage_candidates`]
//! (`substrate_library` + `pattern_substrates` from
//! `tests/data/coverage_substrates.txt`) — same split as
//! historical native `SUBSTRATE_LIBRARY` + `_PATTERN_SUBSTRATES`. Per-leaf
//! [`RuleSet::example_substrates`] / [`crate::rules::LEAF_EXAMPLE_SUBSTRATES`]
//! stay the short site_kind list (native `_example_substrates`); do not grow
//! that table to chase mute patterns.

use std::collections::{BTreeMap, BTreeSet};

use crate::ForestMol;
use crate::as_forest_mol;
use crate::mol::atom_idx;
use crate::pattern::{Edit, Effect, PatternInfo, SiteKind, compose_delta_formula};
use crate::rules::{LEAF_EXAMPLE_SUBSTRATES, catalog_names, leaf_rule};
use crate::substrate_library::coverage_candidates;

/// Conjugation / reactivity adducts — products collapse to chematic `*`.
/// Effect bags declare star stoichiometry (see [`crate::star_conjugate`]).
const ADDUCT_LEAVES: &[&str] = &[
    "Acetylation",
    "Sulfation",
    "Glucuronidation",
    "Glutathionation",
    "GlutathionationNoThiol",
    "GSH",
    "Protein",
    "DNA",
    "Cyanide",
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
        SiteKind::Atom => assert!(n >= 1, "{leaf}/{}: Atom site_map empty", info.name),
        SiteKind::Bond | SiteKind::DirectedBond => assert_eq!(
            n, 2,
            "{leaf}/{}: {:?} site_map must be length 2, got {:?}",
            info.name, info.site_kind, info.site_map
        ),
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
    let example_names: BTreeSet<&str> = LEAF_EXAMPLE_SUBSTRATES.iter().map(|(n, _)| *n).collect();
    assert_eq!(
        example_names.len(),
        LEAF_EXAMPLE_SUBSTRATES.len(),
        "duplicate leaf in LEAF_EXAMPLE_SUBSTRATES"
    );
    for &name in catalog_names() {
        assert!(seen_leaves.insert(name), "duplicate catalog name {name}");
        assert!(
            example_names.contains(name),
            "{name}: missing from LEAF_EXAMPLE_SUBSTRATES (beside LEAF_CTORS in rules.rs)"
        );
        let set = leaf_rule(name).unwrap_or_else(|| panic!("leaf_rule({name})"));
        assert_eq!(set.name.as_deref(), Some(name));
        assert!(
            !set.example_substrates().is_empty(),
            "{name}: example_substrates empty after seal_leaf"
        );
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

/// If a catalog pattern clears aromaticity at site_map atoms on a coverage
/// candidate, `Effect.dearomatizes` capability must be true on the catalog record.
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

        for &smi in coverage_candidates() {
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
                if catalog_info.effect.dearomatizes || catalog_info.effect.cleaves {
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

fn expected_coverage_keys(leaves: &[&str]) -> (BTreeSet<String>, BTreeSet<String>) {
    let mut patterns = BTreeSet::new();
    let mut whens = BTreeSet::new();
    for &name in leaves {
        let set = leaf_rule(name).expect(name);
        for info in set.patterns() {
            patterns.insert(format!("{name}/{}", info.name));
            for (i, arm) in info.possibilities.iter().enumerate() {
                if arm.when.is_some() {
                    whens.insert(format!("{name}/{}#when{i}", info.name));
                }
            }
        }
    }
    (patterns, whens)
}

/// Effect / atom_diff accuracy + pattern/When coverage on
/// [`coverage_candidates`] (crate data file: library + pattern sections).
fn effect_atom_diff_accuracy(leaves: &[&str]) -> (usize, Vec<String>) {
    use crate::atom_diff::{effect_adds_oxygen, effect_removes_h};
    use crate::forest::{formula_delta, molecule_formula};
    use crate::formula_check::check_effect_delta_formula;

    fn effect_adds_h(effect: &Effect) -> bool {
        effect.adds.as_deref().is_some_and(|a| a.contains('H'))
    }

    let (want_patterns, want_whens) = expected_coverage_keys(leaves);
    let mut hit_patterns: BTreeSet<String> = BTreeSet::new();
    let mut hit_whens: BTreeSet<String> = BTreeSet::new();
    let mut misses: Vec<String> = Vec::new();
    let mut hits = 0usize;

    for &name in leaves {
        let set = leaf_rule(name).expect(name);
        let catalog: BTreeMap<String, PatternInfo> = set
            .patterns()
            .into_iter()
            .cloned()
            .map(|p| (p.name.clone(), p))
            .collect();

        for &smi in coverage_candidates() {
            let Ok(parent) = ForestMol::parse(smi) else {
                continue;
            };
            let Ok(cands) = set.candidates(&parent).collect::<Result<Vec<_>, _>>() else {
                continue;
            };
            let parent_f = molecule_formula(parent.mol());
            for c in cands {
                if let Some(pair) = c.pair.as_ref() {
                    hit_patterns.insert(format!("{name}/{}", pair.left.name));
                    hit_patterns.insert(format!("{name}/{}", pair.right.name));
                } else {
                    hit_patterns.insert(format!("{name}/{}", c.pattern.name));
                }

                let infos: Vec<&PatternInfo> = if let Some(pair) = c.pair.as_ref() {
                    vec![&pair.left, &pair.right]
                } else {
                    catalog.get(&c.pattern.name).into_iter().collect()
                };
                for info in infos {
                    if info.possibilities.is_empty() {
                        continue;
                    }
                    let selected = info.possibilities.iter().position(|a| {
                        a.when.as_ref().is_some_and(|w| {
                            w.matches(parent.mol(), &c.mapped)
                                || c.pair
                                    .as_ref()
                                    .is_some_and(|p| w.matches(parent.mol(), &p.map2))
                        })
                    });
                    if let Some(i) = selected {
                        hit_whens.insert(format!("{name}/{}#when{i}", info.name));
                    }
                }

                let Ok(pieces) = c.materialize_mols() else {
                    continue;
                };
                if pieces.is_empty() {
                    // Match without materialize is soft here — same as native
                    // coverage (resolution on SMARTS hit, not product emit).
                    continue;
                }
                // Conjugation adducts collapse to `*` on apply; formula gate
                // must see the same products.
                let pieces: Vec<_> = match crate::star_conjugate::conjugate_star_label(name) {
                    Some(label) => pieces
                        .into_iter()
                        .map(|piece| {
                            let product = parent.from_edit_product(piece);
                            crate::star_conjugate::collapse_conjugate_to_star(&product, label)
                                .mol()
                                .clone()
                        })
                        .collect(),
                    None => pieces,
                };
                hits += 1;
                if let Some(detail) =
                    check_effect_delta_formula(parent.mol(), &c.effect, &pieces, &c.pattern_name)
                {
                    misses.push(format!(
                        "{name}/{} on {smi}: formula mismatch declared {:?} observed {:?}",
                        c.pattern_name, detail.declared, detail.observed
                    ));
                    continue;
                }
                let observed = if c.effect.cleaves && pieces.len() > 1 {
                    let mut counts = BTreeMap::new();
                    for piece in &pieces {
                        for (el, n) in molecule_formula(piece).counts {
                            *counts.entry(el).or_insert(0) += n;
                        }
                    }
                    formula_delta(&parent_f, &crate::forest::Formula { counts, charge: 0 }).counts
                } else {
                    formula_delta(&parent_f, &molecule_formula(&pieces[0])).counts
                };
                let dh = observed.get("H").copied().unwrap_or(0);
                let d_o = observed.get("O").copied().unwrap_or(0);
                if dh > 0 && !effect_adds_h(&c.effect) && !c.effect.cleaves {
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

    let mute_patterns: Vec<_> = want_patterns.difference(&hit_patterns).cloned().collect();
    if !mute_patterns.is_empty() {
        misses.push(format!(
            "patterns never matched on coverage_candidates ({}):\n  {}\n  \
             (expand tests/data/coverage_substrates.txt [pattern])",
            mute_patterns.len(),
            mute_patterns.join("\n  ")
        ));
    }
    let mute_whens: Vec<_> = want_whens.difference(&hit_whens).cloned().collect();
    if !mute_whens.is_empty() {
        misses.push(format!(
            "When arms never selected on coverage_candidates ({}):\n  {}\n  \
             (expand tests/data/coverage_substrates.txt [pattern])",
            mute_whens.len(),
            mute_whens.join("\n  ")
        ));
    }

    (hits, misses)
}

#[test]
fn catalog_effect_and_atom_diff_match_materialized_products() {
    //! Every non-adduct PatternInfo + When arm on
    //! `coverage_candidates()` (crate data file): match coverage, then
    //! declared Effect (incl. H) matches product deltas when materialize succeeds.
    let leaves: Vec<&str> = catalog_names()
        .iter()
        .copied()
        .filter(|n| !ADDUCT_LEAVES.contains(n))
        .collect();
    let (hits, misses) = effect_atom_diff_accuracy(&leaves);
    assert!(
        hits > 0,
        "expected at least one materialize hit across coverage_candidates"
    );
    assert!(
        misses.is_empty(),
        "PatternInfo Effect / coverage disagree with coverage_substrates.txt.\n\
         When OR arms disagree, use When possibilities — do not drop H from the check.\n  {}",
        misses.join("\n  ")
    );
}

/// Same Effect/atom_diff + coverage gate for conjugation adducts only.
#[test]
fn catalog_adduct_effect_and_atom_diff_match_materialized_products() {
    let (hits, misses) = effect_atom_diff_accuracy(ADDUCT_LEAVES);
    assert!(
        hits > 0,
        "expected at least one adduct materialize hit on coverage_candidates"
    );
    assert!(
        misses.is_empty(),
        "adduct Effect / coverage disagree with coverage_substrates.txt.\n  {}",
        misses.join("\n  ")
    );
}

#[test]
fn catalog_resolve_dearomatizes_on_aromatic_probes() {
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
        for &smi in coverage_candidates() {
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
                assert_eq!(
                    c.effect.dearomatizes, site_aromatic,
                    "{name}/{} on {smi}: resolved dearomatizes={} site_aromatic={}",
                    c.pattern.name, c.effect.dearomatizes, site_aromatic
                );
            }
        }
    }
}
