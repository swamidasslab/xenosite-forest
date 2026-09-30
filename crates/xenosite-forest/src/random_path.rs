//! Seeded random walks over a [`RuleSet`].
//!
//! Helper for tests and sampling: pick a site from [`RuleSet::sites`], apply
//! it, repeat. Deterministic for a given seed. Does not search toward a target.
//!
//! Product filtering (skip multi-component / skip already-seen CSMIs) is
//! controlled by [`crate::pathway::PathwayOptions`] — defaults are off.

use std::collections::HashSet;

use crate::ForestError;
use crate::candidate::Candidate;
use crate::forest_mol::{IntoForestMol, as_forest_mol};
use crate::pathway::PathwayOptions;
use crate::pattern::PatternInfo;
use crate::ruleset::RuleSet;

/// One applied hop in a [`random_path`] walk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RandomPathStep {
    /// Leaf rule name when known (else pattern name).
    pub rule: String,
    pub pattern_name: String,
    /// Discovery site atoms (see [`Candidate::discovery_atoms`]).
    pub site: Vec<usize>,
    /// Product CSMIs from this apply (explicit downgrade; before picking one).
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
    /// PatternInfo for each applied hop (see [`Candidate::bookkeeping_pattern`]).
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

/// Randomly apply up to `max_steps` rules from `ruleset` starting at `reactant`.
///
/// Same as [`random_path_with`] with [`PathwayOptions::default`] (no product
/// filters).
pub fn random_path(
    reactant: impl IntoForestMol,
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
/// Stops early when no site yields an eligible product under `options`, or when
/// apply yields nothing.
pub fn random_path_with(
    reactant: impl IntoForestMol,
    seed: u64,
    ruleset: &RuleSet,
    max_steps: usize,
    options: PathwayOptions,
) -> Result<RandomPathOutcome, ForestError> {
    // Tracing continues when reactant is already a ForestMol.
    let mut mol = as_forest_mol(reactant)?;
    let mut rng = XorShift64::new(seed);
    let mut steps = Vec::new();
    let mut patterns = Vec::new();
    let start_smi = mol.csmi().as_ref().to_string();
    let mut path = vec![start_smi.clone()];
    let mut seen = HashSet::new();
    seen.insert(start_smi);

    for _ in 0..max_steps {
        // Stable order so the same seed picks the same hop across runs.
        let mut sites: Vec<Candidate> = ruleset.sites(&mol).collect::<Result<_, _>>()?;
        sites.sort_by_key(|c| {
            (
                u8::from(c.is_pair()),
                c.pattern_name.clone(),
                c.site,
                c.discovery_atoms(),
            )
        });
        let mut advanced = false;
        while !sites.is_empty() {
            let site = sites.swap_remove(rng.index(sites.len()));
            let Some(emission) = site.apply()? else {
                continue;
            };
            let rule = site
                .leaf_rule()
                .unwrap_or(site.pattern_name.as_str())
                .to_string();
            let pattern_name = site.pattern_name.clone();
            let discovery = site.discovery_atoms();
            let pattern = site.bookkeeping_pattern().clone();
            let mut ranked = emission.products;
            ranked.sort_by(|a, b| a.csmi().cmp(&b.csmi()));
            // Explicit CSMI downgrade for pathway options / step record.
            let products: Vec<String> = ranked
                .iter()
                .map(|m| m.csmi().as_ref().to_string())
                .collect();
            let eligible: Vec<usize> = products
                .iter()
                .enumerate()
                .filter(|(_, smi)| options.allows(smi, &seen))
                .map(|(i, _)| i)
                .collect();
            if eligible.is_empty() {
                continue;
            }
            let chosen = eligible[rng.index(eligible.len())];
            let next = ranked[chosen].clone();
            let chosen_smi = products[chosen].clone();
            seen.insert(chosen_smi.clone());
            path.push(chosen_smi);
            steps.push(RandomPathStep {
                rule,
                pattern_name,
                site: discovery,
                products,
                chosen,
            });
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
