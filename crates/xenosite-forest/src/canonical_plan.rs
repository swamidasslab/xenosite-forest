//! Elementary plans: the plan **is** [`Deps`] (steps + precedes).
//!
//! No parallel `CanonicalStep` dialect. A hop emits elementary [`Step`]s;
//! [`Deps::bind`] rewrites [`PlanAtom::WillAdd`] → [`PlanAtom::AddedBy`] and
//! builds precedes. Composite hops expand via [`PlanKind`] data on the leaf
//! [`crate::ruleset::RuleSet`] — not a rule-name branch in search.
//!
//! Replay: [`Deps::linearizations`] → [`Linearization::apply`] through named
//! elementary rules at resolved sites.

use std::collections::{HashSet, VecDeque};
use std::ops::Deref;

use crate::ForestError;
use crate::mol::{Molecule, atom_idx, atom_usize, canon_of, canon_smiles, parse_mol};
use crate::pattern::Effect;

/// One atom note in a [`Step`] site.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum PlanAtom {
    /// Known index on the mol at emit / resolve time.
    Index(usize),
    /// Atom a prep step will add: element at this anchor index.
    /// Becomes [`PlanAtom::AddedBy`] when the plan is bound.
    WillAdd { element: String, at: usize },
    /// Atom created by an earlier elementary step (rule + that step's anchors).
    AddedBy { rule: String, anchors: Vec<usize> },
}

impl PlanAtom {
    pub fn index(idx: usize) -> Self {
        Self::Index(idx)
    }

    pub fn will_add(element: impl Into<String>, at: usize) -> Self {
        Self::WillAdd {
            element: element.into(),
            at,
        }
    }

    pub fn oxygen_at(at: usize) -> Self {
        Self::will_add("O", at)
    }

    /// Anchor index this note depends on (known atom or will-add site).
    pub fn anchor(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::WillAdd { at, .. } => Some(*at),
            Self::AddedBy { anchors, .. } => anchors.first().copied(),
        }
    }
}

/// One elementary reaction at a site. The plan language — not a search-only note.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Step {
    pub rule: String,
    pub site: Vec<PlanAtom>,
}

/// Compat alias while call sites migrate.
pub type CanonicalStep = Step;

impl Step {
    pub fn new(rule: impl Into<String>, site: impl IntoIterator<Item = PlanAtom>) -> Self {
        Self {
            rule: rule.into(),
            site: site.into_iter().collect(),
        }
    }

    /// Origin / will-add anchors named by this step's site notes.
    pub fn anchors(&self) -> HashSet<usize> {
        self.site.iter().filter_map(PlanAtom::anchor).collect()
    }

    /// Resolve site notes to current heavy-atom indices on `mol`.
    pub fn resolve_site(&self, mol: &Molecule) -> Result<HashSet<usize>, ForestError> {
        let mut out = HashSet::new();
        for note in &self.site {
            out.insert(resolve_atom(note, mol)?);
        }
        Ok(out)
    }

    /// Run this elementary rule at the resolved site; return product mols.
    ///
    /// Products keep atom indices from the edit (no SMILES round-trip) so
    /// later steps' anchors still resolve.
    pub fn apply(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        let wanted = self.resolve_site(mol)?;
        let Some(rule) = crate::rules::leaf_rule(&self.rule) else {
            return Err(ForestError::Plan(format!(
                "unknown plan rule {}",
                self.rule
            )));
        };
        let mut products = Vec::new();
        let mut seen = HashSet::new();
        for c in rule.candidates(mol)? {
            if !wanted.contains(&c.site) {
                continue;
            }
            for p in c.materialize_mols(mol)? {
                let smi = canon_smiles(&p);
                if seen.insert(smi) {
                    products.push(p);
                }
            }
        }
        let endpoints = rule.leaf_pair_endpoints();
        if !endpoints.is_empty() {
            for pair in crate::pair_edit::pair_candidates(mol, &endpoints)? {
                if !pair_matches_wanted(mol, &pair, &wanted) {
                    continue;
                }
                for p in pair.materialize_mols(mol)? {
                    let smi = canon_smiles(&p);
                    if seen.insert(smi) {
                        products.push(p);
                    }
                }
            }
        }
        Ok(products)
    }
}

fn pair_matches_wanted(
    mol: &Molecule,
    pair: &crate::pair_edit::PairCandidate,
    wanted: &HashSet<usize>,
) -> bool {
    let Some((a, b)) = pair.end_atoms() else {
        return false;
    };
    if wanted.len() == 2 && wanted.contains(&a) && wanted.contains(&b) {
        return true;
    }
    // Plan named partner heteroatoms (phenol O, etc.).
    let partners: HashSet<usize> = [a, b]
        .into_iter()
        .filter_map(|carbon| {
            mol.neighbors(atom_idx(carbon)).find_map(|(n, _)| {
                let z = mol.atom(n).element.atomic_number();
                (z == 8 || z == 7 || z == 16).then_some(atom_usize(n))
            })
        })
        .collect();
    partners.len() == wanted.len() && partners.iter().all(|p| wanted.contains(p))
}

fn resolve_atom(note: &PlanAtom, mol: &Molecule) -> Result<usize, ForestError> {
    match note {
        PlanAtom::Index(i) => {
            if *i >= mol.atom_count() {
                return Err(ForestError::Plan(format!("site index {i} out of range")));
            }
            Ok(*i)
        }
        PlanAtom::WillAdd { element, at } => resolve_added_element(mol, *at, element),
        PlanAtom::AddedBy { anchors, .. } => {
            let Some(&at) = anchors.first() else {
                return Err(ForestError::Plan("AddedBy with empty anchors".into()));
            };
            resolve_added_element(mol, at, "O")
                .or_else(|_| resolve_added_element(mol, at, "N"))
                .or_else(|_| resolve_added_element(mol, at, "S"))
        }
    }
}

fn resolve_added_element(mol: &Molecule, at: usize, element: &str) -> Result<usize, ForestError> {
    let z = match element {
        "O" => 8u8,
        "N" => 7,
        "S" => 16,
        "C" => 6,
        _ => {
            return Err(ForestError::Plan(format!(
                "cannot resolve added element {element}"
            )));
        }
    };
    bonded(mol, at, z).ok_or_else(|| {
        ForestError::Plan(format!(
            "no {element} bonded to atom {at} (will-add / added-by unresolved)"
        ))
    })
}

/// Flat elementary steps plus precedes (transitive reduction).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Deps {
    steps: Vec<Step>,
    precedes: Vec<(usize, usize)>,
}

impl Deps {
    /// Build from already-bound steps + raw precedes; stores the reduction.
    pub fn new(
        steps: impl IntoIterator<Item = Step>,
        precedes: impl IntoIterator<Item = (usize, usize)>,
    ) -> Self {
        let steps: Vec<_> = steps.into_iter().collect();
        let n = steps.len();
        let raw: Vec<_> = precedes.into_iter().collect();
        for &(a, b) in &raw {
            assert!(
                a < n && b < n && a != b,
                "invalid precedes ({a}, {b}) for n={n}"
            );
        }
        let precedes = if n == 0 {
            Vec::new()
        } else {
            canonical_dependency_edges(n, &raw).expect("cycle in precedes")
        };
        Self { steps, precedes }
    }

    /// Bind will-add notes → added-by, then build precedes (the plan identity).
    pub fn bind(steps: impl IntoIterator<Item = Step>) -> Self {
        bind_deps(steps.into_iter().collect())
    }

    pub fn steps(&self) -> &[Step] {
        &self.steps
    }

    pub fn precedes(&self) -> &[(usize, usize)] {
        &self.precedes
    }

    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Every topological sort under precedes (tiny graphs).
    pub fn linearizations(&self) -> Vec<Linearization> {
        let n = self.steps.len();
        if n == 0 {
            return vec![Linearization { steps: Vec::new() }];
        }
        let mut outgoing = vec![Vec::new(); n];
        let mut indeg = vec![0usize; n];
        for &(a, b) in &self.precedes {
            outgoing[a].push(b);
            indeg[b] += 1;
        }
        let mut out = Vec::new();
        let mut path = Vec::new();
        let mut indeg_work = indeg.clone();
        fn rec(
            steps: &[Step],
            outgoing: &[Vec<usize>],
            indeg: &mut [usize],
            path: &mut Vec<usize>,
            out: &mut Vec<Linearization>,
        ) {
            if path.len() == steps.len() {
                out.push(Linearization {
                    steps: path.iter().map(|&i| steps[i].clone()).collect(),
                });
                return;
            }
            let ready: Vec<_> = (0..steps.len()).filter(|&i| indeg[i] == 0).collect();
            for i in ready {
                indeg[i] = usize::MAX; // mark used
                path.push(i);
                for &b in &outgoing[i] {
                    indeg[b] -= 1;
                }
                rec(steps, outgoing, indeg, path, out);
                for &b in &outgoing[i] {
                    indeg[b] += 1;
                }
                path.pop();
                indeg[i] = 0;
            }
        }
        rec(&self.steps, &outgoing, &mut indeg_work, &mut path, &mut out);
        out
    }

    /// True if some linearization apply reaches `target` CSMI (or canon spelling).
    pub fn reaches(&self, reactant: &str, target: &str) -> Result<bool, ForestError> {
        let want = canon_of(target)?;
        let mol = parse_mol(reactant)?;
        for lin in self.linearizations() {
            for product in lin.apply(&mol)? {
                if canon_smiles(&product) == want {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

impl Deref for Deps {
    type Target = [Step];

    fn deref(&self) -> &Self::Target {
        &self.steps
    }
}

/// Ordered elementary steps; apply chains resolve → rule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Linearization {
    pub steps: Vec<Step>,
}

impl Linearization {
    /// Apply steps in order. Returns product molecules after the last step.
    pub fn apply(&self, mol: &Molecule) -> Result<Vec<Molecule>, ForestError> {
        if self.steps.is_empty() {
            return Ok(vec![mol.clone()]);
        }
        let mut currents = vec![mol.clone()];
        for step in &self.steps {
            let mut next = Vec::new();
            let mut seen = HashSet::new();
            for cur in &currents {
                for product in step.apply(cur)? {
                    let smi = canon_smiles(&product);
                    if seen.insert(smi) {
                        next.push(product);
                    }
                }
            }
            if next.is_empty() {
                return Ok(Vec::new());
            }
            currents = next;
        }
        Ok(currents)
    }
}

/// Unique minimal dependency edge set (transitive reduction).
pub fn canonical_dependency_edges(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<(usize, usize)>, &'static str> {
    let closure = transitive_closure(n, edges)?;
    let mut canonical = Vec::new();
    for a in 0..n {
        let descendants = &closure[a];
        let mut through: HashSet<usize> = HashSet::new();
        for &b in descendants {
            through.extend(&closure[b]);
        }
        for &b in descendants {
            if !through.contains(&b) {
                canonical.push((a, b));
            }
        }
    }
    canonical.sort_unstable();
    Ok(canonical)
}

fn transitive_closure(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<HashSet<usize>>, &'static str> {
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut indegree = vec![0usize; n];
    let mut seen = HashSet::new();
    for &(a, b) in edges {
        if a >= n || b >= n {
            return Err("invalid edge");
        }
        if a == b {
            return Err("cycle in precedes");
        }
        if seen.insert((a, b)) {
            outgoing[a].push(b);
            indegree[b] += 1;
        }
    }
    let mut queue: VecDeque<usize> = (0..n).filter(|&v| indegree[v] == 0).collect();
    let mut topo = Vec::with_capacity(n);
    while let Some(a) = queue.pop_front() {
        topo.push(a);
        for &b in &outgoing[a] {
            indegree[b] -= 1;
            if indegree[b] == 0 {
                queue.push_back(b);
            }
        }
    }
    if topo.len() != n {
        return Err("cycle in precedes");
    }
    let mut closure: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for &a in topo.iter().rev() {
        for &b in &outgoing[a] {
            closure[a].insert(b);
            let child = closure[b].clone();
            closure[a].extend(child);
        }
    }
    Ok(closure)
}

/// Bind will-add → added-by and collect precedes from those notes.
pub fn bind_deps(steps: Vec<Step>) -> Deps {
    let mut edges = Vec::new();
    let mut bound = Vec::with_capacity(steps.len());
    for (later, step) in steps.iter().enumerate() {
        let mut site = Vec::with_capacity(step.site.len());
        for item in &step.site {
            match item {
                PlanAtom::WillAdd { element: _, at } => {
                    let mut bound_note = PlanAtom::Index(*at);
                    for (earlier, previous) in steps[..later].iter().enumerate() {
                        if previous.anchors().contains(at) {
                            edges.push((earlier, later));
                            let mut anchors: Vec<_> = previous.anchors().into_iter().collect();
                            anchors.sort_unstable();
                            bound_note = PlanAtom::AddedBy {
                                rule: previous.rule.clone(),
                                anchors,
                            };
                        }
                    }
                    site.push(bound_note);
                }
                PlanAtom::AddedBy { rule, anchors } => {
                    let wanted: HashSet<_> = anchors.iter().copied().collect();
                    for (earlier, previous) in steps.iter().enumerate() {
                        if previous.rule == *rule && previous.anchors() == wanted {
                            edges.push((earlier, later));
                        }
                    }
                    site.push(item.clone());
                }
                PlanAtom::Index(_) => site.push(item.clone()),
            }
        }
        bound.push(Step {
            rule: step.rule.clone(),
            site,
        });
    }
    Deps::new(bound, edges)
}

/// [`Deps::bind`] alias (Python `as_deps`).
pub fn as_deps(steps: impl IntoIterator<Item = Step>) -> Deps {
    Deps::bind(steps)
}

/// How a leaf expands an accepted hop into elementary [`Step`]s.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PlanKind {
    /// Already elementary: one step at the discovery site.
    #[default]
    Identity,
    /// Prep missing oxygens, then one dehydrogenation (quinone-formation shape).
    HydroxylationThenDehydrogenation,
}

/// Identity plan: `rule` at the given site atoms.
pub fn identity_plan(rule: impl Into<String>, site: impl IntoIterator<Item = usize>) -> Vec<Step> {
    vec![Step::new(rule, site.into_iter().map(PlanAtom::index))]
}

/// Compat name.
pub fn identity_canonical_plan(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    identity_plan(rule, site)
}

fn end_needs_oxygen(effect: &Effect) -> bool {
    let partner = effect.partner.as_deref().unwrap_or("");
    let needs_o = effect.adds.as_deref().is_some_and(|a| a.contains('O'));
    needs_o && partner != "O"
}

fn bonded(mol: &Molecule, idx: usize, atomic_num: u8) -> Option<usize> {
    mol.neighbors(atom_idx(idx)).find_map(|(nbr, _)| {
        let n = atom_usize(nbr);
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
    mol: &Molecule,
    ends: &[&Effect],
    end_atoms: &[usize],
) -> Vec<Step> {
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
                preps.push(Step::new(
                    "OxidativeDehalogenation",
                    [anchor.clone(), PlanAtom::index(halo)],
                ));
            } else {
                preps.push(Step::new("Hydroxylation", [anchor.clone()]));
            }
            dh_refs.push(PlanAtom::oxygen_at(atom));
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
    preps.push(Step::new("Dehydrogenation", dh_refs));
    preps
}

/// Elementary steps from [`PlanKind`] (unbound; call [`Deps::bind`] for the plan).
pub fn steps_for_kind(
    kind: PlanKind,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    match kind {
        PlanKind::Identity => identity_plan(rule_name, site_atoms.iter().copied()),
        PlanKind::HydroxylationThenDehydrogenation => {
            if let (Some(ends), true) = (end_effects, site_atoms.len() >= 2) {
                let plan = hydroxylation_then_dehydrogenation(mol, ends, site_atoms);
                if !plan.is_empty() {
                    return plan;
                }
            }
            identity_plan("Dehydrogenation", site_atoms.iter().copied())
        }
    }
}

/// Bound [`Deps`] for one hop.
pub fn plan_for_kind(
    kind: PlanKind,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Deps {
    Deps::bind(steps_for_kind(
        kind,
        rule_name,
        mol,
        site_atoms,
        end_effects,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pattern::Effect;

    #[test]
    fn identity_is_one_step() {
        let plan = identity_plan("Hydroxylation", [0]);
        assert_eq!(plan.len(), 1);
        assert_eq!(plan[0].rule, "Hydroxylation");
        assert_eq!(plan[0].site, vec![PlanAtom::index(0)]);
    }

    #[test]
    fn hydroquinone_pair_needs_no_prep() {
        let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
        let carbons: Vec<usize> = (0..mol.atom_count())
            .filter(|&i| {
                let a = atom_idx(i);
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
        assert!(matches!(plan[2].site[0], PlanAtom::WillAdd { .. }));
    }

    #[test]
    fn bind_benzene_qf_precedes_both_oh_before_dh() {
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
        let deps = Deps::bind(plan);
        assert_eq!(deps.len(), 3);
        let edges: HashSet<_> = deps.precedes().iter().copied().collect();
        assert!(!edges.contains(&(0, 1)));
        assert!(!edges.contains(&(1, 0)));
        assert!(edges.contains(&(0, 2)));
        assert!(edges.contains(&(1, 2)));
        assert!(matches!(deps[2].site[0], PlanAtom::AddedBy { .. }));
        assert_eq!(deps.linearizations().len(), 2); // OH arms commute
    }

    #[test]
    fn bind_phenol_one_oh_precedes_dh() {
        let plan = vec![
            Step::new("Hydroxylation", [PlanAtom::index(1)]),
            Step::new(
                "Dehydrogenation",
                [PlanAtom::index(0), PlanAtom::oxygen_at(1)],
            ),
        ];
        let deps = Deps::bind(plan);
        assert_eq!(deps.precedes(), &[(0, 1)]);
    }

    #[test]
    fn canonical_edges_drop_transitive() {
        let reduced = canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert_eq!(reduced, vec![(0, 1), (1, 2)]);
    }

    #[test]
    fn replay_benzene_qf_plan_reaches_quinone() {
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
        let deps = Deps::bind(hydroxylation_then_dehydrogenation(&mol, &ends, &[0, 3]));
        assert!(
            deps.reaches("c1ccccc1", "O=C1C=CC(=O)C=C1").unwrap(),
            "plan={deps:?}"
        );
    }

    #[test]
    fn replay_hydroxylation_ethane() {
        let mol = parse_mol("CC").unwrap();
        let rule = crate::rules::hydroxylation();
        let em = rule
            .metabolize(&mol, |_, _, _| true, |_, _, _| true, true)
            .unwrap();
        assert!(!em.is_empty());
        let deps = Deps::bind(identity_plan("Hydroxylation", [em[0].site]));
        assert!(deps.reaches("CC", "CCO").unwrap());
    }
}
