//! Seeded random walks over a [`RuleSet`].
//!
//! Helper for tests and sampling: pick a candidate (or ResonancePair), apply
//! it, repeat. Deterministic for a given seed. Does not search toward a target.
//!
//! Product filtering (skip multi-component / skip already-seen CSMIs) is
//! controlled by [`crate::pathway::PathwayOptions`] — defaults are off.

use std::collections::HashSet;

use crate::ForestError;
use crate::candidate::Candidate;
use crate::mol::{Molecule, canon_smiles, parse_mol};
use crate::pair_edit::PairCandidate;
use crate::pathway::PathwayOptions;
use crate::pattern::PatternInfo;
use crate::ruleset::RuleSet;

/// One applied hop in a [`random_path`] walk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RandomPathStep {
    /// Leaf rule name when known (else pattern name).
    pub rule: String,
    pub pattern_name: String,
    /// Discovery site atoms (sorted unique-edit orbit, or pair end atoms).
    pub site: Vec<usize>,
    /// Product CSMIs from this materialize (before picking one to continue).
    pub products: Vec<String>,
    /// Index into [`Self::products`] that became the next reactant.
    pub chosen: usize,
}

/// Outcome of [`random_path`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RandomPathOutcome {
    /// Final molecule CSMI (same as [`Self::path`] last entry).
    pub smiles: String,
    /// Reactant CSMI followed by each chosen product along the walk.
    pub path: Vec<String>,
    pub steps: Vec<RandomPathStep>,
    /// PatternInfo for each applied hop (pair hops use the left endpoint pattern
    /// plus the right name in [`RandomPathStep::pattern_name`]).
    pub patterns: Vec<PatternInfo>,
}

/// Seeded XorShift — no `rand` dependency (WASM-clean).
#[derive(Clone, Debug)]
struct XorShift64(u64);

impl XorShift64 {
    fn new(seed: u64) -> Self {
        // Avoid the all-zero fixed point.
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        if state == 0 {
            state = 0xA5A5_A5A5_A5A5_A5A5;
        }
        Self(state)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn index(&mut self, len: usize) -> usize {
        debug_assert!(len > 0);
        (self.next_u64() as usize) % len
    }
}

enum Pending {
    Candidate(Box<Candidate>),
    Pair(Box<PairCandidate>),
}

/// Randomly apply up to `max_steps` rules from `ruleset` starting at `reactant`.
///
/// Same as [`random_path_with`] with [`PathwayOptions::default`] (no product
/// filters).
pub fn random_path(
    reactant: &str,
    seed: u64,
    ruleset: &RuleSet,
    max_steps: usize,
) -> Result<RandomPathOutcome, ForestError> {
    random_path_with(
        reactant,
        seed,
        ruleset,
        max_steps,
        PathwayOptions::default(),
    )
}

/// [`random_path`] with explicit [`PathwayOptions`].
///
/// Stops early when no pending edit yields an eligible product under `options`,
/// or when a materialize is empty.
pub fn random_path_with(
    reactant: &str,
    seed: u64,
    ruleset: &RuleSet,
    max_steps: usize,
    options: PathwayOptions,
) -> Result<RandomPathOutcome, ForestError> {
    let mut mol = parse_mol(reactant)?;
    let mut rng = XorShift64::new(seed);
    let mut steps = Vec::new();
    let mut patterns = Vec::new();
    let start_smi = canon_smiles(&mol);
    let mut path = vec![start_smi.clone()];
    let mut seen = HashSet::new();
    seen.insert(start_smi);

    for _ in 0..max_steps {
        let mut pending = collect_pending(ruleset, &mol)?;
        let mut advanced = false;
        while !pending.is_empty() {
            let pick = rng.index(pending.len());
            let item = pending.swap_remove(pick);
            let (step, pattern, next) = match item {
                Pending::Candidate(c) => {
                    apply_candidate(c.as_ref(), &mol, &mut rng, &seen, options)?
                }
                Pending::Pair(p) => apply_pair(p.as_ref(), &mol, &mut rng, &seen, options)?,
            };
            let Some(next) = next else {
                // Try another pending (empty materialize or all products filtered).
                continue;
            };
            let chosen_smi = step.products[step.chosen].clone();
            seen.insert(chosen_smi.clone());
            path.push(chosen_smi);
            steps.push(step);
            patterns.push(pattern);
            mol = next;
            advanced = true;
            break;
        }
        if !advanced {
            break;
        }
    }

    Ok(RandomPathOutcome {
        smiles: path.last().cloned().unwrap_or_default(),
        path,
        steps,
        patterns,
    })
}

fn collect_pending(ruleset: &RuleSet, mol: &Molecule) -> Result<Vec<Pending>, ForestError> {
    let mut out = Vec::new();
    for c in ruleset.candidates(mol) {
        out.push(Pending::Candidate(Box::new(c?)));
    }
    for p in ruleset.pair_candidates(mol) {
        out.push(Pending::Pair(Box::new(p?)));
    }
    // Stable order so the same seed picks the same hop across runs (HashMap
    // iteration order is otherwise process-randomized).
    out.sort_by_key(pending_sort_key);
    Ok(out)
}

fn pending_sort_key(p: &Pending) -> (u8, String, usize, Vec<usize>) {
    match p {
        Pending::Candidate(c) => (0, c.pattern.name.clone(), c.site, c.orbit.clone()),
        Pending::Pair(p) => (1, p.pattern_name.clone(), p.site, pair_site(p)),
    }
}

fn product_indices(
    products: &[String],
    seen: &HashSet<String>,
    options: PathwayOptions,
) -> Vec<usize> {
    products
        .iter()
        .enumerate()
        .filter(|(_, smi)| options.allows(smi, seen))
        .map(|(i, _)| i)
        .collect()
}

/// Sort pieces by CSMI so product choice is seed-stable.
fn sorted_pieces(pieces: Vec<Molecule>) -> (Vec<Molecule>, Vec<String>) {
    let mut pairs: Vec<(String, Molecule)> =
        pieces.into_iter().map(|m| (canon_smiles(&m), m)).collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    let products: Vec<String> = pairs.iter().map(|(s, _)| s.clone()).collect();
    let mols: Vec<Molecule> = pairs.into_iter().map(|(_, m)| m).collect();
    (mols, products)
}

fn apply_candidate(
    c: &Candidate,
    mol: &Molecule,
    rng: &mut XorShift64,
    seen: &HashSet<String>,
    options: PathwayOptions,
) -> Result<(RandomPathStep, PatternInfo, Option<Molecule>), ForestError> {
    let pieces = c.materialize_mols(mol)?;
    let rule = c.leaf_rule().unwrap_or(c.pattern_name()).to_string();
    if pieces.is_empty() {
        return Ok((
            RandomPathStep {
                rule,
                pattern_name: c.pattern_name().to_string(),
                site: c.orbit.clone(),
                products: Vec::new(),
                chosen: 0,
            },
            c.pattern.clone(),
            None,
        ));
    }
    let (pieces, products) = sorted_pieces(pieces);
    let eligible = product_indices(&products, seen, options);
    if eligible.is_empty() {
        return Ok((
            RandomPathStep {
                rule,
                pattern_name: c.pattern_name().to_string(),
                site: c.orbit.clone(),
                products,
                chosen: 0,
            },
            c.pattern.clone(),
            None,
        ));
    }
    let chosen = eligible[rng.index(eligible.len())];
    Ok((
        RandomPathStep {
            rule,
            pattern_name: c.pattern_name().to_string(),
            site: c.orbit.clone(),
            products,
            chosen,
        },
        c.pattern.clone(),
        Some(pieces[chosen].clone()),
    ))
}

fn apply_pair(
    p: &PairCandidate,
    mol: &Molecule,
    rng: &mut XorShift64,
    seen: &HashSet<String>,
    options: PathwayOptions,
) -> Result<(RandomPathStep, PatternInfo, Option<Molecule>), ForestError> {
    let pieces = p.materialize_mols(mol)?;
    if pieces.is_empty() {
        return Ok((
            RandomPathStep {
                rule: p.pattern_name.clone(),
                pattern_name: p.pattern_name.clone(),
                site: pair_site(p),
                products: Vec::new(),
                chosen: 0,
            },
            p.left.clone(),
            None,
        ));
    }
    let (pieces, products) = sorted_pieces(pieces);
    let eligible = product_indices(&products, seen, options);
    if eligible.is_empty() {
        return Ok((
            RandomPathStep {
                rule: p.pattern_name.clone(),
                pattern_name: p.pattern_name.clone(),
                site: pair_site(p),
                products,
                chosen: 0,
            },
            p.left.clone(),
            None,
        ));
    }
    let chosen = eligible[rng.index(eligible.len())];
    Ok((
        RandomPathStep {
            rule: p.pattern_name.clone(),
            pattern_name: p.pattern_name.clone(),
            site: pair_site(p),
            products,
            chosen,
        },
        p.left.clone(),
        Some(pieces[chosen].clone()),
    ))
}

fn pair_site(p: &PairCandidate) -> Vec<usize> {
    match p.end_atoms() {
        Some((a, b)) => {
            let mut v = vec![a, b];
            v.sort_unstable();
            v.dedup();
            v
        }
        None => {
            let (a, b) = p.path_ends();
            let mut v = vec![a, b];
            v.sort_unstable();
            v.dedup();
            v
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydroxylation::hydroxylation;
    use crate::rules::phase_one;
    use crate::ruleset::o_dealkylation;

    #[test]
    fn ethane_hydroxylation_is_deterministic() {
        let rules = hydroxylation();
        let a = random_path("CC", 1, &rules, 1).unwrap();
        let b = random_path("CC", 1, &rules, 1).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.steps.len(), 1);
        assert_eq!(a.steps[0].rule, "Hydroxylation");
        assert!(
            a.smiles == "CCO" || a.smiles == "C(C)O",
            "unexpected ethanol spelling {}",
            a.smiles
        );
        assert_eq!(a.path.len(), 2);
        assert_eq!(a.path[0], "CC");
        assert_eq!(a.path[1], a.smiles);
        assert_eq!(a.patterns.len(), 1);
    }

    #[test]
    fn same_seed_is_identical() {
        let rules = phase_one();
        let a = random_path("COc1ccccc1", 42, &rules, 3).unwrap();
        let b = random_path("COc1ccccc1", 42, &rules, 3).unwrap();
        assert_eq!(a, b);
        assert!(!a.steps.is_empty(), "expected a non-empty walk");
    }

    #[test]
    fn different_seeds_diverge() {
        let rules = phase_one();
        let a = random_path("COc1ccccc1", 1, &rules, 3).unwrap();
        let b = random_path("COc1ccccc1", 99, &rules, 3).unwrap();
        assert!(
            a != b,
            "seeds 1 and 99 produced identical walks: {:?}",
            a.path
        );
        assert_eq!(a.path.len(), a.steps.len() + 1);
        assert_eq!(b.path.len(), b.steps.len() + 1);
    }

    #[test]
    fn no_loops_preset_skips_multicomponent() {
        let rules = o_dealkylation();
        let out = random_path_with(
            "COc1ccccc1",
            7,
            &rules,
            1,
            PathwayOptions::no_loops_or_fragments(),
        )
        .unwrap();
        assert_eq!(out.steps.len(), 1);
        assert!(!out.smiles.contains('.'));
    }

    #[test]
    fn default_options_do_not_filter_products() {
        // Contrasts with no_loops_or_fragments: default walk keeps every
        // non-empty piece (including ones that would be filtered when opted in).
        let rules = o_dealkylation();
        let open = random_path("COc1ccccc1", 7, &rules, 1).unwrap();
        let filtered = random_path_with(
            "COc1ccccc1",
            7,
            &rules,
            1,
            PathwayOptions::no_loops_or_fragments(),
        )
        .unwrap();
        assert_eq!(open.steps.len(), 1);
        assert_eq!(filtered.steps.len(), 1);
        assert!(!filtered.smiles.contains('.'));
        // Open path may choose any piece; filters only constrain when set.
        assert!(PathwayOptions::default().allows("Oc1ccccc1.C=O", &HashSet::new()));
    }

    #[test]
    fn anisole_o_dealk_one_step() {
        let rules = o_dealkylation();
        let out = random_path("COc1ccccc1", 7, &rules, 1).unwrap();
        assert_eq!(out.steps.len(), 1);
        assert!(!out.steps[0].products.is_empty());
        let joined = out.steps[0].products.join("|");
        assert!(
            joined.contains("Oc1ccccc1") || out.smiles.contains("O"),
            "got smiles={} products={:?}",
            out.smiles,
            out.steps[0].products
        );
    }

    #[test]
    fn zero_steps_returns_reactant() {
        let rules = hydroxylation();
        let out = random_path("CC", 1, &rules, 0).unwrap();
        assert!(out.steps.is_empty());
        assert_eq!(out.smiles, "CC");
        assert_eq!(out.path, vec!["CC".to_string()]);
    }
}
