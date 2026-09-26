//! MS1 m/z-targeted search.
//!
//! Sibling of [`crate::find_path`]: same RuleSet expand / PatternInfo filters
//! spirit, but the goal is a mass (m/z ± tol), not a structure CSMI. Plan
//! constraints use [`crate::canonical_plan::ApplyN`] (OR of transforms, apply
//! exactly N) composable with [`crate::canonical_plan::Deps`] / [`Maybe`].
//!
//! No MCS atom_diff toward a missing structure target. Closer / hit use
//! [`crate::mass`] + PatternInfo `delta_formula` (shared with structure
//! find_path — see `mass` / catalog tests for drift guards).

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use crate::ForestError;
use crate::canonical_plan::{ApplyN, Step, as_deps};
use crate::find_path::{PathCounters, PathOutcome, PathStep};
use crate::forest_mol::ForestMol;
use crate::mass::{Ms1Adduct, formula_apply_delta, mz_abs_error, mz_of, mz_of_mol, mz_within};
use crate::ruleset::RuleSet;

/// MS1 search bounds + mass target.
#[derive(Clone, Copy, Debug)]
pub struct Ms1Config {
    pub mz: f64,
    /// Absolute tolerance in Da (not ppm).
    pub tol_da: f64,
    pub adduct: Ms1Adduct,
    pub max_paths: usize,
    pub max_nodes: usize,
}

impl Default for Ms1Config {
    fn default() -> Self {
        Self {
            mz: 0.0,
            tol_da: 0.01,
            adduct: Ms1Adduct::MPlusH,
            max_paths: 1,
            max_nodes: 800,
        }
    }
}

impl Ms1Config {
    pub fn with_mz(mz: f64) -> Self {
        Self {
            mz,
            ..Self::default()
        }
    }
}

#[derive(Clone)]
struct Ms1Walk {
    mol: ForestMol,
    steps: Vec<PathStep>,
    plan: Vec<Step>,
    /// Applications counted per [`ApplyN`] pool index.
    pool_used: Vec<u16>,
    /// |mz − target| at this mol (Da).
    mz_err: f64,
    ancestors: HashSet<String>,
}

#[derive(Clone)]
struct Ms1HeapItem {
    /// Higher is better: −mz_err as fixed-point µDa, then seq.
    neg_err_ud: i64,
    seq: usize,
    walk: Ms1Walk,
}

impl PartialEq for Ms1HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.neg_err_ud == other.neg_err_ud && self.seq == other.seq
    }
}
impl Eq for Ms1HeapItem {}
impl PartialOrd for Ms1HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Ms1HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.neg_err_ud
            .cmp(&other.neg_err_ud)
            .then_with(|| other.seq.cmp(&self.seq))
    }
}

fn err_key(err: f64) -> i64 {
    -((err * 1_000_000.0).round() as i64)
}

fn remember_seen(seen: &mut HashSet<String>, mol: &ForestMol) -> bool {
    match mol.stable_csmi_key() {
        Some(k) => seen.insert(k.as_ref().to_string()),
        None => true,
    }
}

fn root_ancestors(start: &ForestMol) -> HashSet<String> {
    let mut a = HashSet::new();
    if let Some(k) = start.stable_csmi_key() {
        a.insert(k.as_ref().to_string());
    }
    a
}

fn pools_satisfied(pools: &[ApplyN], used: &[u16]) -> bool {
    pools.iter().zip(used.iter()).all(|(p, &u)| u == p.count)
}

fn pools_can_accept(pools: &[ApplyN], used: &[u16], rule: &str) -> Option<usize> {
    for (i, pool) in pools.iter().enumerate() {
        if pool.allows(rule) && used[i] < pool.count {
            return Some(i);
        }
    }
    None
}

fn any_pool_open(pools: &[ApplyN], used: &[u16]) -> bool {
    pools.iter().zip(used.iter()).any(|(p, &u)| u < p.count)
}

fn walk_mz(mol: &ForestMol, adduct: Ms1Adduct) -> Option<f64> {
    // Prefer atom walk so explicit isotope labels shift mass; formula path
    // is element-symbol only (common isotopes).
    mz_of_mol(mol.mol(), adduct).or_else(|| mz_of(mol.formula().as_ref(), adduct))
}

/// Predict child mz from PatternInfo `delta_formula` (no materialize).
pub fn predicted_mz_after_delta(
    parent: &ForestMol,
    delta: &std::collections::BTreeMap<String, i32>,
    adduct: Ms1Adduct,
) -> Option<f64> {
    let next = formula_apply_delta(parent.formula().as_ref(), delta);
    mz_of(&next, adduct)
}

/// Yield walks whose formula m/z matches `config` and whose [`ApplyN`] pools
/// are exactly filled.
///
/// `pools` empty ⇒ any ruleset chemistry may fire until mass hits
/// (unconstrained). Non-empty ⇒ each hop must feed some unsatisfied pool.
pub fn find_path_ms1(
    reactant: &str,
    ruleset: &RuleSet,
    pools: &[ApplyN],
    counters: &mut PathCounters,
    config: Ms1Config,
) -> Result<Vec<PathOutcome>, ForestError> {
    let start = ForestMol::parse(reactant)?;
    let start_mz = walk_mz(&start, config.adduct).ok_or_else(|| {
        ForestError::Plan("reactant formula has unsupported elements for mono mass".into())
    })?;
    let start_err = mz_abs_error(start_mz, config.mz);
    let ancestors = root_ancestors(&start);

    let mut heap = BinaryHeap::new();
    let mut seq = 0usize;
    let mut seen = HashSet::new();
    remember_seen(&mut seen, &start);
    heap.push(Ms1HeapItem {
        neg_err_ud: err_key(start_err),
        seq,
        walk: Ms1Walk {
            mol: start,
            steps: Vec::new(),
            plan: Vec::new(),
            pool_used: vec![0u16; pools.len()],
            mz_err: start_err,
            ancestors,
        },
    });
    seq += 1;

    let mut found = Vec::new();
    while let Some(item) = heap.pop() {
        if found.len() >= config.max_paths || counters.nodes >= config.max_nodes {
            break;
        }
        counters.nodes += 1;
        let walk = item.walk;

        let here_mz = match walk_mz(&walk.mol, config.adduct) {
            Some(m) => m,
            None => continue,
        };
        let hit_mass = mz_within(here_mz, config.mz, config.tol_da);
        let hit_pools = pools_satisfied(pools, &walk.pool_used);
        if hit_mass && hit_pools {
            let plan = as_deps(walk.plan.clone()).with_apply_n(pools.iter().cloned());
            found.push(PathOutcome {
                steps: walk.steps,
                plan,
                smiles: walk.mol.csmi().as_ref().to_string(),
            });
            continue;
        }
        if pools.is_empty() {
            if hit_mass {
                continue;
            }
        } else if !any_pool_open(pools, &walk.pool_used) {
            continue;
        }

        counters.expansions += 1;
        let mol = walk.mol.mol();
        let gens = walk.mol.atom_bond_generators();
        for cand in ruleset.candidates(mol) {
            let cand = cand?;
            let rule_name = cand.leaf_rule().unwrap_or(cand.pattern.name.as_str());

            let pool_idx = if pools.is_empty() {
                None
            } else {
                match pools_can_accept(pools, &walk.pool_used, rule_name) {
                    Some(i) => Some(i),
                    None => continue,
                }
            };

            let delta = cand.pattern.effect.resolved_delta_formula();
            // Soft pre-filter (non-cleavage only): try full declared delta and
            // heavy-only (strip H). Hydroxyl bags are +O −H while live mols are
            // +O — heavy wins. H-only nets (epoxide rearrange +2H) need full.
            // Cleavage leave bags do not predict either fragment's mono mass —
            // skip the hint and let materialize + closer decide.
            if !cand.pattern.effect.cleaves {
                let full_mz = predicted_mz_after_delta(&walk.mol, &delta, config.adduct);
                let mut heavy_delta = delta.clone();
                heavy_delta.remove("H");
                let heavy_mz = if heavy_delta != delta {
                    predicted_mz_after_delta(&walk.mol, &heavy_delta, config.adduct)
                } else {
                    None
                };
                let pred_mz = match (full_mz, heavy_mz) {
                    (Some(a), Some(b)) => {
                        let ea = mz_abs_error(a, config.mz);
                        let eb = mz_abs_error(b, config.mz);
                        Some(if ea <= eb { a } else { b })
                    }
                    (a, b) => a.or(b),
                };
                if let Some(pred_mz) = pred_mz {
                    let pred_err = mz_abs_error(pred_mz, config.mz);
                    if pred_err > walk.mz_err + 1e-9 && pred_err > config.tol_da {
                        continue;
                    }
                }
            }

            let pieces = cand.materialize_mols(mol)?;
            if pieces.is_empty() {
                continue;
            }
            counters.mol_edits += 1;
            let child_mol = walk.mol.adopt_product(pieces[0].clone());
            if let Some(k) = child_mol.stable_csmi_key() {
                if walk.ancestors.contains(k.as_ref()) {
                    continue;
                }
            }
            if !remember_seen(&mut seen, &child_mol) {
                continue;
            }

            let child_mz = match walk_mz(&child_mol, config.adduct) {
                Some(m) => m,
                None => continue,
            };
            let child_err = mz_abs_error(child_mz, config.mz);
            if child_err > walk.mz_err + 1e-9 {
                continue;
            }

            let mut pool_used = walk.pool_used.clone();
            if let Some(i) = pool_idx {
                pool_used[i] = pool_used[i].saturating_add(1);
            }

            let plan_steps = cand.identity_plan_with_gens(&gens, mol.atom_count());
            let path_step = PathStep {
                rule_path: cand.rule_path.clone(),
                pattern_name: cand.pattern.name.clone(),
                site: cand.site,
                site_orbit: cand.orbit.clone(),
                product: child_mol.csmi().as_ref().to_string(),
                sides: Vec::new(),
            };

            let mut ancestors = walk.ancestors.clone();
            if let Some(k) = child_mol.stable_csmi_key() {
                ancestors.insert(k.as_ref().to_string());
            }
            let mut steps = walk.steps.clone();
            steps.push(path_step);
            let mut plan = walk.plan.clone();
            plan.extend(plan_steps);

            heap.push(Ms1HeapItem {
                neg_err_ud: err_key(child_err),
                seq,
                walk: Ms1Walk {
                    mol: child_mol,
                    steps,
                    plan,
                    pool_used,
                    mz_err: child_err,
                    ancestors,
                },
            });
            seq += 1;
        }
    }

    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical_plan::ApplyN;
    use crate::forest::molecule_formula;
    use crate::mass::{Ms1Adduct, mz_of, mz_of_mol, mz_within};
    use crate::mol::parse_mol;
    use crate::rules::{
        dealkylation, dehydrogenation, epoxidation, epoxide_hydration, epoxide_opening,
        hydroxylation, nitrogen_oxidation, phase_one, sulfur_oxidation,
    };
    use crate::ruleset::RuleSet;

    /// Every emitted hit's final product (and last path-step product) must
    /// lie within `tol_da` of the target m/z. Steps must get monotonically
    /// closer. Intermediates need not equal the target mass.
    fn assert_hits_satisfy_mz(hits: &[PathOutcome], reactant: &str, mz: f64, tol_da: f64) {
        assert!(!hits.is_empty(), "expected at least one MS1 hit");
        let start = ForestMol::parse(reactant).unwrap();
        let start_mz = mz_of_mol(start.mol(), Ms1Adduct::MPlusH).unwrap();
        for h in hits {
            let mut prev_err = crate::mass::mz_abs_error(start_mz, mz);
            for (i, step) in h.steps.iter().enumerate() {
                let step_mol = ForestMol::parse(&step.product).unwrap();
                let step_mz = mz_of_mol(step_mol.mol(), Ms1Adduct::MPlusH).unwrap();
                let err = crate::mass::mz_abs_error(step_mz, mz);
                assert!(
                    err <= prev_err + 1e-6,
                    "step {i} mz err grew: {prev_err} → {err} ({})",
                    step.product
                );
                prev_err = err;
                if i + 1 == h.steps.len() {
                    assert!(
                        mz_within(step_mz, mz, tol_da),
                        "last step product {} mz={step_mz} outside tol of {mz}",
                        step.product
                    );
                    assert_eq!(
                        step.product, h.smiles,
                        "last path step product must equal hit smiles"
                    );
                }
            }
            let mol = ForestMol::parse(&h.smiles).unwrap();
            let got = mz_of_mol(mol.mol(), Ms1Adduct::MPlusH).unwrap();
            assert!(
                mz_within(got, mz, tol_da),
                "hit smiles {} mz={got} outside tol of target {mz}",
                h.smiles
            );
        }
    }

    /// Apply `leaf` once to `reactant`; return distinct (csmi, mz) products.
    fn products_from_apply(reactant: &str, set: &RuleSet) -> Vec<(String, f64)> {
        let parent = ForestMol::parse(reactant).unwrap();
        let mol = parent.mol();
        let mut out = Vec::new();
        for cand in set.candidates(mol) {
            let cand = cand.unwrap();
            let pieces = cand.materialize_mols(mol).unwrap();
            if pieces.is_empty() {
                continue;
            }
            let child = parent.adopt_product(pieces[0].clone());
            let mz = mz_of_mol(child.mol(), Ms1Adduct::MPlusH).unwrap();
            out.push((child.csmi().as_ref().to_string(), mz));
        }
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out.dedup_by(|a, b| a.0 == b.0);
        out
    }

    /// Materialize products from `set`, then assert each product CSMI appears
    /// in `find_path_ms1` hits when targeting that product's [M+H]⁺ m/z.
    fn assert_applied_products_in_ms1(
        reactant: &str,
        set: &RuleSet,
        pool_arms: &[&str],
        count: u16,
    ) {
        let products = products_from_apply(reactant, set);
        assert!(
            !products.is_empty(),
            "{reactant} + {:?} produced nothing",
            set.name
        );
        let pools = [ApplyN::new(pool_arms.iter().copied(), count)];
        for (csmi, mz) in &products {
            let mut counters = PathCounters::default();
            let hits = find_path_ms1(
                reactant,
                set,
                &pools,
                &mut counters,
                Ms1Config {
                    mz: *mz,
                    tol_da: 0.001,
                    adduct: Ms1Adduct::MPlusH,
                    max_paths: 16,
                    max_nodes: 400,
                },
            )
            .unwrap();
            assert!(
                hits.iter().any(|h| h.smiles == *csmi),
                "{reactant} → {csmi} (mz={mz}) missing from MS1 hits {:?}; billed={}",
                hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
                counters.billed()
            );
            assert_hits_satisfy_mz(&hits, reactant, *mz, 0.001);
        }
    }

    #[test]
    fn applied_hydroxylation_products_appear_in_ms1() {
        assert_applied_products_in_ms1("CC", &hydroxylation(), &["Hydroxylation"], 1);
        assert_applied_products_in_ms1("c1ccccc1", &hydroxylation(), &["Hydroxylation"], 1);
        // Propane: primary + secondary alcohols share mass; both must hit.
        assert_applied_products_in_ms1("CCC", &hydroxylation(), &["Hydroxylation"], 1);
    }

    #[test]
    fn applied_epoxidation_product_appears_in_ms1() {
        assert_applied_products_in_ms1("C=C", &epoxidation(), &["Epoxidation"], 1);
    }

    #[test]
    fn applied_epoxide_hydration_product_appears_in_ms1() {
        assert_applied_products_in_ms1("C=C", &epoxide_hydration(), &["EpoxideHydration"], 1);
    }

    #[test]
    fn applied_epoxide_opening_products_appear_in_ms1() {
        // rearrange (+2H → ethanol) and hydrate (+O+2H → glycol).
        assert_applied_products_in_ms1("C1OC1", &epoxide_opening(), &["EpoxideOpening"], 1);
    }

    #[test]
    fn applied_sulfur_oxidation_products_appear_in_ms1() {
        assert_applied_products_in_ms1("CCS", &sulfur_oxidation(), &["SulfurOxidation"], 1);
    }

    #[test]
    fn applied_nitrogen_oxidation_products_appear_in_ms1() {
        assert_applied_products_in_ms1("CN", &nitrogen_oxidation(), &["NitrogenOxidation"], 1);
    }

    #[test]
    fn applied_dehydrogenation_products_appear_in_ms1() {
        assert_applied_products_in_ms1("CCO", &dehydrogenation(), &["Dehydrogenation"], 1);
    }

    /// Apply `leaves` in order; return final (csmi, mz).
    fn chain_apply(reactant: &str, leaves: &[&str]) -> (String, f64) {
        let mut cur = ForestMol::parse(reactant).unwrap();
        for leaf in leaves {
            let set = crate::rules::leaf_rule(leaf).unwrap_or_else(|| panic!("missing {leaf}"));
            let mol = cur.mol();
            let mut next = None;
            for cand in set.candidates(mol) {
                let cand = cand.unwrap();
                let pieces = cand.materialize_mols(mol).unwrap();
                if pieces.is_empty() {
                    continue;
                }
                next = Some(cur.adopt_product(pieces[0].clone()));
                break;
            }
            cur = next.unwrap_or_else(|| panic!("{reactant} chain failed at {leaf}"));
        }
        let mz = mz_of_mol(cur.mol(), Ms1Adduct::MPlusH).unwrap();
        (cur.csmi().as_ref().to_string(), mz)
    }

    fn assert_chain_in_ms1(
        reactant: &str,
        set: &RuleSet,
        pool_arms: &[&str],
        count: u16,
        leaves: &[&str],
        max_nodes: usize,
    ) {
        let (csmi, mz) = chain_apply(reactant, leaves);
        let pools = [ApplyN::new(pool_arms.iter().copied(), count)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            reactant,
            set,
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 16,
                max_nodes,
            },
        )
        .unwrap();
        assert!(
            hits.iter().any(|h| h.smiles == csmi),
            "chain {leaves:?} → {csmi} (mz={mz}) missing from {:?}; billed={}",
            hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
            counters.billed()
        );
        assert_hits_satisfy_mz(&hits, reactant, mz, 0.001);
    }

    #[test]
    fn harder_two_hydroxylations_ethane_and_benzene() {
        assert_chain_in_ms1(
            "CC",
            &hydroxylation(),
            &["Hydroxylation"],
            2,
            &["Hydroxylation", "Hydroxylation"],
            200,
        );
        assert_chain_in_ms1(
            "c1ccccc1",
            &hydroxylation(),
            &["Hydroxylation"],
            2,
            &["Hydroxylation", "Hydroxylation"],
            200,
        );
        assert_chain_in_ms1(
            "CCC",
            &hydroxylation(),
            &["Hydroxylation"],
            2,
            &["Hydroxylation", "Hydroxylation"],
            400,
        );
    }

    #[test]
    fn harder_ethene_epoxidation_then_opening() {
        // Two elementary hops under PhaseOne OR pool.
        assert_chain_in_ms1(
            "C=C",
            &phase_one(),
            &[
                "Epoxidation",
                "EpoxideOpening",
                "EpoxideHydration",
                "Hydroxylation",
            ],
            2,
            &["Epoxidation", "EpoxideOpening"],
            800,
        );
    }

    #[test]
    fn harder_ethene_to_glycol_via_phase_one_or_pool() {
        // One-hop EpoxideHydration inside a broad PhaseOne OR pool (count=1).
        let glycol = ForestMol::parse("OCCO").unwrap();
        let mz = mz_of_mol(glycol.mol(), Ms1Adduct::MPlusH).unwrap();
        let pools = [ApplyN::new(
            [
                "Hydroxylation",
                "Epoxidation",
                "EpoxideOpening",
                "EpoxideHydration",
            ],
            1,
        )];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            "C=C",
            &phase_one(),
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 8,
                max_nodes: 800,
            },
        )
        .unwrap();
        assert!(
            hits.iter().any(|h| h.smiles == glycol.csmi().as_ref()),
            "glycol missing from {:?}; billed={}",
            hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
            counters.billed()
        );
        assert_hits_satisfy_mz(&hits, "C=C", mz, 0.001);
    }

    #[test]
    fn harder_anisole_dealkylation_phenol_in_ms1() {
        // Cleavage: soft delta hint skipped; keep-side phenol must still hit.
        let set = dealkylation();
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let mut phenol = None;
        for cand in set.candidates(parent.mol()) {
            let cand = cand.unwrap();
            if !cand.pattern.effect.cleaves {
                continue;
            }
            for p in cand.materialize_mols(parent.mol()).unwrap() {
                let child = parent.adopt_product(p);
                let f = molecule_formula(child.mol());
                if f.counts.get("C") == Some(&6)
                    && f.counts.get("O") == Some(&1)
                    && f.counts.get("H") == Some(&6)
                {
                    let mz = mz_of_mol(child.mol(), Ms1Adduct::MPlusH).unwrap();
                    phenol = Some((child.csmi().as_ref().to_string(), mz));
                    break;
                }
            }
            if phenol.is_some() {
                break;
            }
        }
        let (csmi, mz) = phenol.expect("anisole dealk should emit phenol");
        let pools = [ApplyN::new(["Dealkylation"], 1)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            "COc1ccccc1",
            &set,
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 16,
                max_nodes: 800,
            },
        )
        .unwrap();
        assert!(
            hits.iter().any(|h| h.smiles == csmi),
            "phenol {csmi} missing from {:?}; billed={}",
            hits.iter().map(|h| h.smiles.as_str()).collect::<Vec<_>>(),
            counters.billed()
        );
        assert_hits_satisfy_mz(&hits, "COc1ccccc1", mz, 0.001);
    }

    #[test]
    fn harder_three_hydroxylations_on_benzene() {
        assert_chain_in_ms1(
            "c1ccccc1",
            &hydroxylation(),
            &["Hydroxylation"],
            3,
            &["Hydroxylation", "Hydroxylation", "Hydroxylation"],
            2000,
        );
    }

    #[test]
    fn harder_toluene_two_oh_under_phase_one_or() {
        assert_chain_in_ms1(
            "Cc1ccccc1",
            &phase_one(),
            &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dehydrogenation",
            ],
            2,
            &["Hydroxylation", "Hydroxylation"],
            3000,
        );
    }

    #[test]
    fn harder_butylbenzene_two_hydroxylations() {
        assert_chain_in_ms1(
            "CCCCc1ccccc1",
            &hydroxylation(),
            &["Hydroxylation"],
            2,
            &["Hydroxylation", "Hydroxylation"],
            3000,
        );
    }

    #[test]
    fn harder_ethene_glycol_then_dehydrogenation() {
        assert_chain_in_ms1(
            "C=C",
            &phase_one(),
            &[
                "EpoxideHydration",
                "Dehydrogenation",
                "Hydroxylation",
                "Epoxidation",
                "EpoxideOpening",
            ],
            2,
            &["EpoxideHydration", "Dehydrogenation"],
            3000,
        );
    }

    #[test]
    fn harder_veratrole_dealkylation() {
        assert_chain_in_ms1(
            "COc1ccc(OC)cc1",
            &dealkylation(),
            &["Dealkylation"],
            1,
            &["Dealkylation"],
            2000,
        );
    }

    #[test]
    fn harder_n_dealkylation_dimethylaniline() {
        assert_chain_in_ms1(
            "CN(C)c1ccccc1",
            &crate::rules::n_dealkylation(),
            &["NDealkylation"],
            1,
            &["NDealkylation"],
            2000,
        );
    }

    /// Deterministic regression corpus for the apply→MS1 property.
    /// Full proptest draw lives in `tests/ms1_apply_fuzz.rs`.
    #[test]
    fn fuzz_applied_products_recoverable_by_ms1() {
        let cases: &[(&str, &str)] = &[
            ("CC", "Hydroxylation"),
            ("CCC", "Hydroxylation"),
            ("c1ccccc1", "Hydroxylation"),
            ("Cc1ccccc1", "Hydroxylation"),
            ("CCCCc1ccccc1", "Hydroxylation"),
            ("CCO", "Hydroxylation"),
            ("C=C", "Epoxidation"),
            ("C=C", "EpoxideHydration"),
            ("C/C=C/C", "Epoxidation"),
            ("c1ccccc1", "Epoxidation"),
            ("C1OC1", "EpoxideOpening"),
            ("CC1OC1C", "EpoxideOpening"),
            ("CCS", "SulfurOxidation"),
            ("CSC", "SulfurOxidation"),
            ("CN", "NitrogenOxidation"),
            ("CCN", "NitrogenOxidation"),
            ("CCO", "Dehydrogenation"),
            ("CC(O)C", "Dehydrogenation"),
            ("COc1ccccc1", "Dealkylation"),
            ("COc1ccc(OC)cc1", "Dealkylation"),
            ("CN(C)c1ccccc1", "NDealkylation"),
            ("C=C", "Hydrogenation"),
            ("C#C", "Hydrogenation"),
        ];
        let mut checked = 0usize;
        for &(reactant, leaf) in cases {
            let Some(set) = crate::rules::leaf_rule(leaf) else {
                panic!("missing leaf {leaf}");
            };
            let products = products_from_apply(reactant, &set);
            if products.is_empty() {
                continue;
            }
            assert_applied_products_in_ms1(reactant, &set, &[leaf], 1);
            checked += products.len();
        }
        assert!(
            checked >= 20,
            "fuzz property covered too few products: {checked}"
        );
    }

    #[test]
    fn ethene_to_glycol_ms1_one_epoxide_hydration() {
        let glycol = molecule_formula(&parse_mol("OCCO").unwrap());
        let mz = mz_of(&glycol, Ms1Adduct::MPlusH).unwrap();
        let pools = [ApplyN::new(["EpoxideHydration"], 1)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            "C=C",
            &epoxide_hydration(),
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 1,
                max_nodes: 50,
            },
        )
        .unwrap();
        assert_eq!(hits.len(), 1, "{hits:?} billed={}", counters.billed());
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("EpoxideHydration"));
        // Catalog bags are mass-faithful for this leaf — prediction ≈ hit.
        let parent = ForestMol::parse("C=C").unwrap();
        let delta = epoxide_hydration().patterns()[0]
            .effect
            .resolved_delta_formula();
        let pred = predicted_mz_after_delta(&parent, &delta, Ms1Adduct::MPlusH).unwrap();
        assert!(mz_within(pred, mz, 0.001), "pred={pred} target={mz}");
    }

    #[test]
    fn ethane_to_ethanol_ms1_one_hydroxylation() {
        let ethanol = molecule_formula(&parse_mol("CCO").unwrap());
        let mz = mz_of(&ethanol, Ms1Adduct::MPlusH).unwrap();
        let pools = [ApplyN::new(["Hydroxylation"], 1)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            "CC",
            &hydroxylation(),
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 1,
                max_nodes: 50,
            },
        )
        .unwrap();
        assert_eq!(hits.len(), 1, "{hits:?} billed={}", counters.billed());
        assert_eq!(hits[0].plan.steps().len(), 1);
        assert_eq!(hits[0].plan.steps()[0].rule, "Hydroxylation");
        assert_eq!(hits[0].plan.apply_n().len(), 1);
        assert_eq!(hits[0].plan.apply_n()[0].count, 1);
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Hydroxylation"));
    }

    #[test]
    fn predicted_catalog_delta_is_soft_hint_vs_materialized() {
        // Catalog +O −H does not equal sanitized +O; prediction is a hint.
        // Hit / closer use materialized formula mass (see ethane_to_ethanol_ms1).
        let parent = ForestMol::parse("CC").unwrap();
        let set = hydroxylation();
        let delta = set.patterns()[0].effect.resolved_delta_formula();
        let pred = predicted_mz_after_delta(&parent, &delta, Ms1Adduct::MPlusH).unwrap();
        let c = set.candidates(parent.mol()).next().unwrap().unwrap();
        let child = parent.adopt_product(c.materialize_mols(parent.mol()).unwrap()[0].clone());
        let obs = walk_mz(&child, Ms1Adduct::MPlusH).unwrap();
        assert!(
            (pred - obs).abs() > 0.5,
            "expected H-bag drift: pred={pred} obs={obs}"
        );
        // Heavy O agrees: both move toward ethanol-scale mass from ethane.
        let ethane_mz = walk_mz(&parent, Ms1Adduct::MPlusH).unwrap();
        assert!(pred > ethane_mz && obs > ethane_mz);
    }

    #[test]
    fn refuse_when_apply_n_count_exceeds_reachable_at_mass() {
        let ethanol = molecule_formula(&parse_mol("CCO").unwrap());
        let mz = mz_of(&ethanol, Ms1Adduct::MPlusH).unwrap();
        let pools = [ApplyN::new(["Hydroxylation"], 3)];
        let mut counters = PathCounters::default();
        let hits = find_path_ms1(
            "CC",
            &hydroxylation(),
            &pools,
            &mut counters,
            Ms1Config {
                mz,
                tol_da: 0.001,
                adduct: Ms1Adduct::MPlusH,
                max_paths: 1,
                max_nodes: 80,
            },
        )
        .unwrap();
        assert!(hits.is_empty(), "{hits:?}");
    }

    #[test]
    fn structure_find_path_ethane_ethanol_still_works() {
        // Drift guard: MS1 work must not break structure find_path.
        use crate::find_path::find_path;
        let mut counters = PathCounters::default();
        let hits = find_path("CC", "CCO", &hydroxylation(), &mut counters)
            .unwrap()
            .collect_all()
            .unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].steps[0].leaf_rule(), Some("Hydroxylation"));
    }
}
