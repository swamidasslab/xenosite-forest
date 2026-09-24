//! Cleavage product graph: Or of arms that share the same fragment multiset.
//!
//! Both sides of a bifurcation are first-class products (nodes). Arms that
//! produce the same sorted fragment CSMI multiset fold into one [`CleavageOr`].
//! Choosing an arm **and** a continuation fragment builds the linked
//! [`Maybe`] from the other side(s).
//!
//! MCS / [`crate::atom_diff`] gates which fragments stay expandable — the graph
//! does not bypass diff.
//!
//! This is **not** archive And/Or plan trees: Deps stay the plan language.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::rc::Rc;

use crate::ForestError;
use crate::atom_diff::{AtomDiff, atom_diff, candidate_could_help, pair_could_help};
use crate::canonical_plan::{CleavageSide, Maybe};
use crate::forest_mol::ForestMol;
use crate::mol::{Molecule, canon_of};
use crate::pair_edit::PairCandidate;
use crate::ruleset::RuleSet;

/// One way to produce a fragment multiset: rule/site on the parent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageArm {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    pub site_atoms: Vec<usize>,
    /// Both (all) bifurcation fragments, sorted CSMIs.
    pub products: Vec<String>,
}

impl CleavageArm {
    /// Linked Maybe when continuing the walk on `continue_csmi`.
    ///
    /// The other product CSMIs become [`CleavageSide`] entries for this arm.
    pub fn maybe_for(&self, continue_csmi: &str) -> Option<Maybe> {
        if !self.products.iter().any(|p| p == continue_csmi) {
            return None;
        }
        let entries: Vec<_> = self
            .products
            .iter()
            .filter(|p| p.as_str() != continue_csmi)
            .map(|side| {
                CleavageSide::new(
                    self.site_atoms.iter().copied(),
                    side.clone(),
                    std::iter::empty::<BTreeSet<usize>>(),
                )
            })
            .collect();
        Some(Maybe::new(entries))
    }
}

/// Or over arms that share the same sorted fragment multiset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageOr {
    /// Canonical key: sorted fragment CSMIs (both sides kept).
    pub fragments: Vec<String>,
    pub arms: Vec<CleavageArm>,
}

impl CleavageOr {
    pub fn n_arms(&self) -> usize {
        self.arms.len()
    }

    /// Choose an arm and which fragment to continue on; Maybe = the other side(s).
    pub fn choose(&self, arm_index: usize, continue_csmi: &str) -> Option<(&CleavageArm, Maybe)> {
        let arm = self.arms.get(arm_index)?;
        let maybe = arm.maybe_for(continue_csmi)?;
        Some((arm, maybe))
    }
}

/// One-hop cleavage layer: folded Or products from a single mol.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CleavageLayer {
    pub products: Vec<CleavageOr>,
}

impl CleavageLayer {
    pub fn n_products(&self) -> usize {
        self.products.len()
    }

    pub fn n_arms(&self) -> usize {
        self.products.iter().map(|p| p.arms.len()).sum()
    }

    pub fn max_or_fanin(&self) -> usize {
        self.products
            .iter()
            .map(|p| p.arms.len())
            .max()
            .unwrap_or(0)
    }

    /// Distinct fragment CSMIs appearing on either side of any Or.
    pub fn distinct_fragments(&self) -> BTreeSet<&str> {
        self.products
            .iter()
            .flat_map(|o| o.fragments.iter().map(|s| s.as_str()))
            .collect()
    }
}

/// Node in the multi-hop cleavage product graph (one fragment CSMI).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageNode {
    pub csmi: String,
    /// Incoming Or folds that produced this fragment (among others).
    /// `(parent CSMI, Or, this fragment was one of Or.fragments)`.
    pub via: Vec<(String, CleavageOr)>,
}

/// Multi-hop cleavage product graph (BFS of cleaving edits only).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CleavageGraph {
    pub nodes: Vec<CleavageNode>,
}

impl CleavageGraph {
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn n_ors(&self) -> usize {
        self.nodes.iter().map(|n| n.via.len()).sum()
    }

    pub fn n_arms(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|n| n.via.iter())
            .map(|(_, o)| o.arms.len())
            .sum()
    }

    pub fn max_or_fanin(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|n| n.via.iter())
            .map(|(_, o)| o.arms.len())
            .max()
            .unwrap_or(0)
    }
}

/// Options for building a cleavage layer / graph.
#[derive(Clone, Debug)]
pub struct CleavageGraphConfig {
    /// When set, MCS-gate fragments toward this target.
    pub target: Option<String>,
    /// Max distinct fragment CSMIs in a BFS graph (including the root).
    pub max_nodes: usize,
    /// Max BFS depth (root = 0).
    pub max_depth: usize,
}

impl Default for CleavageGraphConfig {
    fn default() -> Self {
        Self {
            target: None,
            max_nodes: 64,
            // First-pass cleave net stays shallow (HEURISTICS).
            max_depth: 4,
        }
    }
}

impl CleavageGraphConfig {
    /// Shallow first-pass toward `target` (depth 4).
    pub fn first_pass(target: impl Into<String>) -> Self {
        Self {
            target: Some(target.into()),
            max_nodes: 64,
            max_depth: 4,
        }
    }
}

fn sorted_fragments(products: &[String]) -> Vec<String> {
    let mut v = products.to_vec();
    v.sort();
    v
}

fn fold_layer(arms: Vec<CleavageArm>) -> CleavageLayer {
    let mut buckets: BTreeMap<Vec<String>, Vec<CleavageArm>> = BTreeMap::new();
    for arm in arms {
        let key = arm.products.clone();
        buckets.entry(key).or_default().push(arm);
    }
    let products = buckets
        .into_iter()
        .map(|(fragments, arms)| CleavageOr { fragments, arms })
        .collect();
    CleavageLayer { products }
}

/// Whether this fragment should be BFS-expanded toward the target.
///
/// Both sides of a split stay on [`CleavageOr::fragments`]. Expansion is
/// per-fragment: if **both** match (strict MCS cost drop and large enough to
/// still reach the target), **both** are enqueued — we do not pick a winner.
pub(crate) fn fragment_worth_expanding(
    parent_diff: Option<&AtomDiff>,
    fragment_csmi: &str,
    target: &Molecule,
    target_csmi: &str,
    target_ha: usize,
) -> bool {
    if fragment_csmi == target_csmi {
        return true;
    }
    let Ok(frag) = ForestMol::parse(fragment_csmi) else {
        return false;
    };
    // Cleavage only shrinks; a piece smaller than the target cannot be a core.
    if frag.heavy_atom_count() < target_ha {
        return false;
    }
    let child_diff = atom_diff(frag.mol(), target);
    match parent_diff {
        Some(parent) => child_diff.cost() < parent.cost(),
        None => true,
    }
}

/// A split is recorded if it bifurcates and at least one fragment is worth
/// expanding (when a target is set). Both fragment CSMIs stay on the Or.
fn split_usable(
    products: &[String],
    parent_diff: Option<&AtomDiff>,
    target: Option<&Molecule>,
    target_csmi: Option<&str>,
    target_ha: Option<usize>,
) -> bool {
    if products.len() < 2 {
        return false;
    }
    match (parent_diff, target, target_csmi, target_ha) {
        (diff, Some(t), Some(tc), Some(tha)) => products
            .iter()
            .any(|p| fragment_worth_expanding(diff, p, t, tc, tha)),
        _ => true,
    }
}

/// One-hop cleavage layer from `mol`, folding arms by sorted fragment multiset.
///
/// Both sides of every bifurcation are retained on [`CleavageArm::products`] /
/// [`CleavageOr::fragments`]. When `config.target` is set, a split is kept if
/// at least one fragment is worth expanding; both CSMIs stay on the Or either way.
pub fn cleavage_layer(
    mol: &Molecule,
    ruleset: &RuleSet,
    config: &CleavageGraphConfig,
) -> Result<CleavageLayer, ForestError> {
    let target_mol = match &config.target {
        Some(t) => Some(crate::mol::parse_mol(t)?),
        None => None,
    };
    let target_csmi = match &config.target {
        Some(t) => Some(canon_of(t)?),
        None => None,
    };
    let target_ha = target_mol.as_ref().map(|t| {
        t.atoms()
            .filter(|(_, a)| a.element.atomic_number() > 1)
            .count()
    });
    let parent_diff = target_mol.as_ref().map(|t| atom_diff(mol, t));

    let mut raw: Vec<CleavageArm> = Vec::new();

    for c in ruleset.candidates(mol)? {
        if !c.pattern.effect.cleaves {
            continue;
        }
        if let Some(diff) = &parent_diff {
            if !candidate_could_help(&c, diff) {
                continue;
            }
        }
        let Some(emission) = c.emit(mol)? else {
            continue;
        };
        if !emission.cleaves {
            continue;
        }
        let products = sorted_fragments(&emission.products);
        if !split_usable(
            &products,
            parent_diff.as_ref(),
            target_mol.as_ref(),
            target_csmi.as_deref(),
            target_ha,
        ) {
            continue;
        }
        let rule = emission
            .leaf_rule()
            .unwrap_or(emission.pattern_name.as_str())
            .to_string();
        raw.push(CleavageArm {
            rule,
            pattern_name: emission.pattern_name,
            site: emission.site,
            site_orbit: emission.site_orbit,
            site_atoms: emission.site_atoms,
            products,
        });
    }

    for pair in ruleset.pair_candidates(mol)? {
        if !pair.effect.cleaves {
            continue;
        }
        if let (Some(diff), Some(t)) = (&parent_diff, &target_mol) {
            if !pair_could_help(&pair, diff, mol, t) {
                continue;
            }
        }
        push_pair_arm(
            mol,
            &pair,
            parent_diff.as_ref(),
            target_mol.as_ref(),
            target_csmi.as_deref(),
            target_ha,
            &mut raw,
        )?;
    }

    Ok(fold_layer(raw))
}

fn push_pair_arm(
    mol: &Molecule,
    pair: &PairCandidate,
    parent_diff: Option<&AtomDiff>,
    target_mol: Option<&Molecule>,
    target_csmi: Option<&str>,
    target_ha: Option<usize>,
    raw: &mut Vec<CleavageArm>,
) -> Result<(), ForestError> {
    let Some(emission) = pair.emit(mol)? else {
        return Ok(());
    };
    let products = sorted_fragments(&emission.products);
    if !split_usable(&products, parent_diff, target_mol, target_csmi, target_ha) {
        return Ok(());
    }
    raw.push(CleavageArm {
        rule: pair
            .pattern_name
            .split('+')
            .next()
            .unwrap_or("Pair")
            .to_string(),
        pattern_name: emission.pattern_name,
        site: emission.site,
        site_orbit: vec![emission.site],
        site_atoms: pair.plan_site_atoms(),
        products,
    });
    Ok(())
}

/// BFS cleavage product graph from `start`.
///
/// Both fragment CSMIs stay on each Or and as graph nodes. BFS **expands**
/// every fragment that matches the target progress gate (strict MCS cost drop,
/// ha ≥ target). If both sides match, both are enqueued — no single-winner
/// prune. Non-matching sides stay on the Or for Maybe / [`CleavageOr::choose`].
pub fn cleavage_product_graph(
    start: &str,
    ruleset: &RuleSet,
    config: &CleavageGraphConfig,
) -> Result<CleavageGraph, ForestError> {
    let root_csmi = canon_of(start)?;
    let mut nodes = vec![CleavageNode {
        csmi: root_csmi.clone(),
        via: Vec::new(),
    }];
    let mut index: BTreeMap<String, usize> = BTreeMap::new();
    index.insert(root_csmi.clone(), 0);

    let mut queue: VecDeque<(usize, usize)> = VecDeque::new();
    let mut scheduled: BTreeSet<usize> = BTreeSet::new();
    queue.push_back((0, 0));
    scheduled.insert(0);

    let target_mol = match &config.target {
        Some(t) => Some(crate::mol::parse_mol(t)?),
        None => None,
    };
    let target_csmi = match &config.target {
        Some(t) => Some(canon_of(t)?),
        None => None,
    };
    let target_ha = target_mol.as_ref().map(|t| {
        t.atoms()
            .filter(|(_, a)| a.element.atomic_number() > 1)
            .count()
    });

    while let Some((ni, depth)) = queue.pop_front() {
        if depth >= config.max_depth || nodes.len() >= config.max_nodes {
            continue;
        }
        let parent_csmi = nodes[ni].csmi.clone();
        let parent_mol = crate::mol::parse_mol(&parent_csmi)?;
        let parent_diff = target_mol.as_ref().map(|t| atom_diff(&parent_mol, t));
        let layer = cleavage_layer(&parent_mol, ruleset, config)?;

        for or in layer.products {
            // Expandable fragments first so leaf scraps cannot consume the node
            // budget ahead of matching cores. When both sides match, both stay
            // expandable and both are enqueued — no single-winner prune.
            let mut ranked: Vec<(&String, bool)> = or
                .fragments
                .iter()
                .map(|frag| {
                    let expand = match (
                        parent_diff.as_ref(),
                        target_mol.as_ref(),
                        target_csmi.as_deref(),
                        target_ha,
                    ) {
                        (diff, Some(t), Some(tc), Some(tha)) => {
                            fragment_worth_expanding(diff, frag, t, tc, tha)
                        }
                        _ => true,
                    };
                    (frag, expand)
                })
                .collect();
            ranked.sort_by_key(|(_, expand)| !expand);

            for (frag, expand) in ranked {
                let child_i = match index.entry(frag.clone()) {
                    std::collections::btree_map::Entry::Occupied(o) => *o.get(),
                    std::collections::btree_map::Entry::Vacant(v) => {
                        if nodes.len() >= config.max_nodes {
                            // Still want both matching sides when possible; stop
                            // allocating new leaf/new nodes under the cap.
                            continue;
                        }
                        let i = nodes.len();
                        nodes.push(CleavageNode {
                            csmi: frag.clone(),
                            via: Vec::new(),
                        });
                        v.insert(i);
                        i
                    }
                };

                let already = nodes[child_i]
                    .via
                    .iter()
                    .any(|(p, o)| p == &parent_csmi && o.fragments == or.fragments);
                if !already {
                    nodes[child_i].via.push((parent_csmi.clone(), or.clone()));
                }

                if expand && depth < config.max_depth && scheduled.insert(child_i) {
                    queue.push_back((child_i, depth + 1));
                }
            }
        }
    }

    Ok(CleavageGraph { nodes })
}

/// One core reached by a shallow cleavage-only BFS (measurement helper).
///
/// Search itself does **not** consume these: [`crate::find_path`] runs the
/// depth-capped cleave-first phase on its own heap walks (tagged
/// [`ForestMol`] stays on [`Walk`]). This builder exists to size the net
/// standalone via `cleavage_graph_bench` / unit tests.
#[derive(Clone)]
pub struct CleavageSeed {
    pub mol: ForestMol,
    /// Diff of `mol` vs the search target (lifted through cleavage when possible).
    pub diff: AtomDiff,
    pub plan: Vec<crate::canonical_plan::Step>,
    pub maybe: Vec<crate::canonical_plan::CleavageSide>,
    /// Path hops for measurement / PathStep-shaped dumps.
    pub hops: Vec<CleavageSeedHop>,
    pub depth: usize,
}

impl CleavageSeed {
    pub fn csmi(&self) -> Rc<str> {
        self.mol.csmi()
    }
}

/// One cleavage hop on a [`CleavageSeed`] walk.
#[derive(Clone, Debug)]
pub struct CleavageSeedHop {
    pub rule: String,
    pub pattern_name: String,
    pub site: usize,
    pub site_orbit: Vec<usize>,
    pub product: String,
    pub sides: Vec<String>,
}

/// BFS cleavage-only seeds toward `target`, depth-capped (default first-pass: 4).
///
/// Returns every expandable core visited (including the root at depth 0). Leaf
/// scraps that fail the expand gate are recorded on Maybe of the continuing
/// seed, not as separate seeds.
///
/// `start` must already be a tagged [`ForestMol`] (same object the search will
/// walk). Child products are [`ForestMol::adopt_product`] — no SMILES round-trip.
/// Child diffs **lift** the parent MCS after cleavage
/// ([`crate::atom_diff::atom_diff_after_cleavage`]) when that already shows a
/// cost drop; otherwise one MCS.
pub fn cleavage_first_seeds(
    start: &ForestMol,
    target_csmi: &str,
    target_mol: &Molecule,
    ruleset: &RuleSet,
    max_depth: usize,
) -> Result<Vec<CleavageSeed>, ForestError> {
    let max_nodes = 64usize;
    let target_ha = target_mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .count();
    let root_diff = atom_diff(start.mol(), target_mol);
    let root_csmi = start.csmi().as_ref().to_string();

    let mut seeds = Vec::new();
    let mut seen: BTreeSet<String> = BTreeSet::new();
    let mut queue: VecDeque<CleavageSeed> = VecDeque::new();
    queue.push_back(CleavageSeed {
        mol: start.copy_mol(),
        diff: root_diff,
        plan: Vec::new(),
        maybe: Vec::new(),
        hops: Vec::new(),
        depth: 0,
    });
    seen.insert(root_csmi);

    while let Some(parent) = queue.pop_front() {
        seeds.push(parent.clone());
        if parent.depth >= max_depth || seeds.len() >= max_nodes {
            continue;
        }

        let mut expandable: Vec<(ForestMol, AtomDiff, CleavageArm)> = Vec::new();

        for c in ruleset.candidates(parent.mol.mol())? {
            if !c.pattern.effect.cleaves {
                continue;
            }
            if !candidate_could_help(&c, &parent.diff) {
                continue;
            }
            let pieces = c.materialize_mols(parent.mol.mol())?;
            if pieces.len() < 2 {
                continue;
            }
            let mut adopted: Vec<ForestMol> = pieces
                .into_iter()
                .map(|p| parent.mol.adopt_product(p))
                .collect();
            adopted.sort_by_key(|m| m.csmi().as_ref().to_string());
            let products: Vec<String> = adopted
                .iter()
                .map(|m| m.csmi().as_ref().to_string())
                .collect();
            let mut child_rows = Vec::new();
            for child in adopted {
                let csmi = child.csmi().as_ref().to_string();
                let is_hit = csmi == target_csmi;
                if !is_hit && child.heavy_atom_count() < target_ha {
                    continue;
                }
                let child_diff = crate::atom_diff::atom_diff_after_cleavage(
                    &parent.mol,
                    &parent.diff,
                    &child,
                    target_mol,
                );
                if is_hit || child_diff.cost() < parent.diff.cost() {
                    child_rows.push((child, child_diff));
                }
            }
            if child_rows.is_empty() {
                continue;
            }
            let rule = c.leaf_rule().unwrap_or(c.pattern.name.as_str()).to_string();
            let arm = CleavageArm {
                rule,
                pattern_name: c.pattern.name.clone(),
                site: c.site,
                site_orbit: c.orbit.clone(),
                site_atoms: {
                    let mut atoms: Vec<usize> = c
                        .pattern
                        .site_map
                        .iter()
                        .filter_map(|m| c.mapped.get(m).copied())
                        .collect();
                    if atoms.is_empty() {
                        atoms.push(c.site);
                    }
                    atoms
                },
                products,
            };
            for (child, child_diff) in child_rows {
                expandable.push((child, child_diff, arm.clone()));
            }
        }

        for pair in ruleset.pair_candidates(parent.mol.mol())? {
            if !pair.effect.cleaves {
                continue;
            }
            if !pair_could_help(&pair, &parent.diff, parent.mol.mol(), target_mol) {
                continue;
            }
            let pieces = pair.materialize_mols(parent.mol.mol())?;
            if pieces.len() < 2 {
                continue;
            }
            let mut adopted: Vec<ForestMol> = pieces
                .into_iter()
                .map(|p| parent.mol.adopt_product(p))
                .collect();
            adopted.sort_by_key(|m| m.csmi().as_ref().to_string());
            let products: Vec<String> = adopted
                .iter()
                .map(|m| m.csmi().as_ref().to_string())
                .collect();
            let mut child_rows = Vec::new();
            for child in adopted {
                let csmi = child.csmi().as_ref().to_string();
                let is_hit = csmi == target_csmi;
                if !is_hit && child.heavy_atom_count() < target_ha {
                    continue;
                }
                let child_diff = crate::atom_diff::atom_diff_after_cleavage(
                    &parent.mol,
                    &parent.diff,
                    &child,
                    target_mol,
                );
                if is_hit || child_diff.cost() < parent.diff.cost() {
                    child_rows.push((child, child_diff));
                }
            }
            if child_rows.is_empty() {
                continue;
            }
            let arm = CleavageArm {
                rule: pair
                    .pattern_name
                    .split('+')
                    .next()
                    .unwrap_or("Pair")
                    .to_string(),
                pattern_name: pair.pattern_name.clone(),
                site: pair.site,
                site_orbit: vec![pair.site],
                site_atoms: pair.plan_site_atoms(),
                products,
            };
            for (child, child_diff) in child_rows {
                expandable.push((child, child_diff, arm.clone()));
            }
        }

        expandable.sort_by_key(|(m, _, _)| m.csmi().as_ref().to_string());

        for (child_mol, child_diff, arm) in expandable {
            let frag = child_mol.csmi().as_ref().to_string();
            if !seen.insert(frag.clone()) {
                continue;
            }
            let Some(maybe_bag) = arm.maybe_for(&frag) else {
                continue;
            };
            let mut plan = parent.plan.clone();
            plan.push(
                crate::canonical_plan::Step::new(
                    arm.rule.clone(),
                    arm.site_atoms
                        .iter()
                        .copied()
                        .map(crate::canonical_plan::PlanAtom::index),
                )
                .with_orbit(arm.site_orbit.iter().copied()),
            );
            let mut maybe = parent.maybe.clone();
            maybe.extend(maybe_bag.entries.iter().cloned());
            let sides: Vec<String> = maybe_bag.sides().into_iter().map(str::to_string).collect();
            let mut hops = parent.hops.clone();
            hops.push(CleavageSeedHop {
                rule: arm.rule.clone(),
                pattern_name: arm.pattern_name.clone(),
                site: arm.site,
                site_orbit: arm.site_orbit.clone(),
                product: frag,
                sides,
            });
            queue.push_back(CleavageSeed {
                mol: child_mol,
                diff: child_diff,
                plan,
                maybe,
                hops,
                depth: parent.depth + 1,
            });
            if seeds.len() + queue.len() >= max_nodes {
                break;
            }
        }
    }

    Ok(seeds)
}

/// Compat: parse `start` once, then [`cleavage_first_seeds`].
pub fn cleavage_first_seeds_smiles(
    start: &str,
    target: &str,
    ruleset: &RuleSet,
    max_depth: usize,
) -> Result<Vec<CleavageSeed>, ForestError> {
    let root = ForestMol::parse(start)?;
    let target_csmi = canon_of(target)?;
    let target_mol = crate::mol::parse_mol(target)?;
    cleavage_first_seeds(&root, &target_csmi, &target_mol, ruleset, max_depth)
}

/// Summary counts for benches / tests.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageGraphStats {
    pub nodes: usize,
    pub ors: usize,
    pub arms: usize,
    pub max_or_fanin: usize,
    pub one_hop_products: usize,
    pub one_hop_arms: usize,
    pub one_hop_fragments: usize,
}

pub fn cleavage_graph_stats(
    start: &str,
    target: Option<&str>,
    ruleset: &RuleSet,
) -> Result<(CleavageGraphStats, CleavageGraph), ForestError> {
    let config = CleavageGraphConfig {
        target: target.map(|t| t.to_string()),
        ..CleavageGraphConfig::default()
    };
    let start_mol = crate::mol::parse_mol(start)?;
    let layer = cleavage_layer(&start_mol, ruleset, &config)?;
    let graph = cleavage_product_graph(start, ruleset, &config)?;
    Ok((
        CleavageGraphStats {
            nodes: graph.n_nodes(),
            ors: graph.n_ors(),
            arms: graph.n_arms(),
            max_or_fanin: graph.max_or_fanin(),
            one_hop_products: layer.n_products(),
            one_hop_arms: layer.n_arms(),
            one_hop_fragments: layer.distinct_fragments().len(),
        },
        graph,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::{leaf_rule, phase_one};

    #[test]
    fn anisole_keeps_both_cleavage_fragments() {
        let mol = crate::mol::parse_mol("COc1ccccc1").unwrap();
        let dealk = leaf_rule("Dealkylation").expect("Dealkylation");
        let layer = cleavage_layer(
            &mol,
            &dealk,
            &CleavageGraphConfig {
                target: Some("Oc1ccccc1".into()),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(layer.n_products() >= 1, "products={:?}", layer.products);
        let or = layer
            .products
            .iter()
            .find(|o| o.fragments.len() >= 2)
            .expect("expected a bifurcation Or");
        assert_eq!(or.fragments.len(), 2, "fragments={:?}", or.fragments);
        // Both sides retained — not a single kept/side pick.
        assert_eq!(or.arms[0].products, or.fragments);
    }

    #[test]
    fn choose_continuation_puts_other_side_in_maybe() {
        let mol = crate::mol::parse_mol("COc1ccc(OC)cc1").unwrap();
        let layer = cleavage_layer(
            &mol,
            &phase_one(),
            &CleavageGraphConfig {
                target: Some("Oc1ccc(O)cc1".into()),
                ..Default::default()
            },
        )
        .unwrap();
        let or = layer
            .products
            .iter()
            .find(|o| o.fragments.len() >= 2)
            .unwrap();
        let continue_on = or.fragments[0].as_str();
        let other = or.fragments[1].as_str();
        let (arm, maybe) = or.choose(0, continue_on).unwrap();
        assert!(arm.products.contains(&continue_on.to_string()));
        assert_eq!(maybe.sides(), vec![other]);
    }

    #[test]
    fn tetramethoxy_biphenyl_or_fanin_and_both_sides() {
        let start = "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC";
        let target = "Oc1ccc(-c2ccc(O)c(O)c2)cc1O";
        let (stats, graph) = cleavage_graph_stats(start, Some(target), &phase_one()).unwrap();
        assert!(
            stats.one_hop_arms >= stats.one_hop_products,
            "stats={stats:?}"
        );
        assert!(
            stats.max_or_fanin >= 2 || stats.one_hop_arms >= 2,
            "stats={stats:?}"
        );
        assert!(
            stats.one_hop_fragments >= 2,
            "both sides should appear, stats={stats:?}"
        );
        assert!(stats.nodes >= 2, "graph should keep multiple fragments");
        assert!(stats.nodes < 40, "nodes={}", stats.nodes);
        // At least one via Or on a non-root node carries ≥2 fragments.
        let both_sides = graph
            .nodes
            .iter()
            .any(|n| n.via.iter().any(|(_, o)| o.fragments.len() >= 2));
        assert!(both_sides, "no Or retained both cleavage fragments");
    }

    #[test]
    fn cleavage_first_seeds_depth_capped_both_sides() {
        let seeds =
            cleavage_first_seeds_smiles("COc1ccc(OC)cc1", "Oc1ccc(O)cc1", &phase_one(), 4).unwrap();
        assert!(seeds.iter().any(|s| s.depth == 0));
        assert!(seeds.iter().any(|s| s.depth >= 1));
        assert!(seeds.iter().all(|s| s.depth <= 4));
        let cores: Vec<_> = seeds.iter().filter(|s| s.depth > 0).collect();
        assert!(!cores.is_empty(), "expected cleavage cores");
        assert!(
            cores.iter().all(|s| s.mol.heavy_atom_count() >= 8),
            "scraps should not be seeds"
        );
    }

    #[test]
    fn skip_non_cleaving_patterns() {
        let mol = crate::mol::parse_mol("c1ccccc1").unwrap();
        let layer = cleavage_layer(
            &mol,
            &crate::hydroxylation::hydroxylation(),
            &CleavageGraphConfig::default(),
        )
        .unwrap();
        assert!(layer.products.is_empty());
    }

    #[test]
    fn both_matching_sides_stay_as_nodes() {
        // Symmetric dimethoxy: cleaving either methoxy toward hydroquinone keeps
        // a core that still matches. Both fragment CSMIs stay as nodes — we do
        // not prune to a single winner when both sides match.
        let start = "COc1ccc(OC)cc1";
        let target = "Oc1ccc(O)cc1";
        let graph = cleavage_product_graph(
            start,
            &phase_one(),
            &CleavageGraphConfig {
                target: Some(target.into()),
                max_nodes: 32,
                max_depth: 4,
            },
        )
        .unwrap();
        let or = graph
            .nodes
            .iter()
            .flat_map(|n| n.via.iter())
            .find(|(_, o)| o.fragments.len() >= 2)
            .map(|(_, o)| o)
            .expect("expected a bifurcation Or with both sides");
        assert_eq!(or.fragments.len(), 2);
        let as_nodes = or
            .fragments
            .iter()
            .filter(|f| graph.nodes.iter().any(|n| n.csmi == **f))
            .count();
        assert_eq!(
            as_nodes,
            2,
            "both matching sides must be first-class nodes, fragments={:?} nodes={:?}",
            or.fragments,
            graph.nodes.iter().map(|n| &n.csmi).collect::<Vec<_>>()
        );
    }
}
