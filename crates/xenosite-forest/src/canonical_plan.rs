//! Elementary plans: the plan **is** [`Deps`] (steps + precedes + maybe).
//!
//! No parallel `CanonicalStep` dialect. A hop emits elementary [`Step`]s;
//! [`Deps::bind`] rewrites [`PlanAtom::WillAdd`] → [`PlanAtom::AddedBy`] and
//! builds precedes. Composite leaves own a [`CanonicalPlanFn`] (Python
//! `canonical_plan`) that returns steps named after existing elementary
//! rules — not a `PlanKind` enum in search.
//!
//! Cleavage fragments discarded by the walk live on the plan as [`Maybe`]
//! (not a sibling on the path outcome, and not searched).
//!
//! Replay: [`Deps::linearizations`] → [`Linearization::apply`] through named
//! elementary rules at resolved sites.

use std::collections::{BTreeSet, HashSet, VecDeque};
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

fn plan_atom_sort_key(a: &PlanAtom) -> (u8, String, Vec<usize>) {
    match a {
        PlanAtom::Index(i) => (0, String::new(), vec![*i]),
        PlanAtom::WillAdd { element, at } => (1, element.clone(), vec![*at]),
        PlanAtom::AddedBy { rule, anchors } => (2, rule.clone(), anchors.clone()),
    }
}

/// One elementary reaction at a site. The plan language — not a search-only note.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Step {
    pub rule: String,
    pub site: Vec<PlanAtom>,
    /// Automorphism orbit of the site atoms (generator closure), sorted.
    /// Empty means unknown / not filled; treat as the resolved site indexes.
    /// Filled from ForestMol-cached gens when the plan is emitted.
    pub orbit: Vec<usize>,
}

/// Compat alias while call sites migrate.
pub type CanonicalStep = Step;

impl Step {
    pub fn new(rule: impl Into<String>, site: impl IntoIterator<Item = PlanAtom>) -> Self {
        let mut site: Vec<_> = site.into_iter().collect();
        site.sort_by_key(plan_atom_sort_key);
        site.dedup();
        Self {
            rule: rule.into(),
            site,
            orbit: Vec::new(),
        }
    }

    pub fn with_orbit(mut self, orbit: impl IntoIterator<Item = usize>) -> Self {
        let mut orbit: Vec<_> = orbit.into_iter().collect();
        orbit.sort_unstable();
        orbit.dedup();
        self.orbit = orbit;
        self
    }

    /// Origin / will-add anchors named by this step's site notes.
    pub fn anchors(&self) -> HashSet<usize> {
        self.site.iter().filter_map(PlanAtom::anchor).collect()
    }

    /// Orbit for equivalence checks: filled orbit, or resolved index anchors.
    pub fn site_orbit(&self) -> Vec<usize> {
        if !self.orbit.is_empty() {
            return self.orbit.clone();
        }
        let mut atoms: Vec<_> = self.anchors().into_iter().collect();
        atoms.sort_unstable();
        atoms
    }

    /// Same rule and site classes under passed / filled orbits.
    pub fn same_site_class(&self, other: &Self) -> bool {
        if self.rule != other.rule {
            return false;
        }
        let a = self.site_orbit();
        let b = other.site_orbit();
        match (a.first(), b.first()) {
            (Some(&ai), Some(&bi)) => crate::same_site_orbit(ai, &a, bi, &b),
            (None, None) => true,
            _ => false,
        }
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
        for c in rule.candidates(mol) {
            let c = c?;
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

/// Fragment discarded at a cleavage. Not searched.
///
/// `site` is the cleavage site on the parent. `opens` are earlier ring-open
/// sites on the same walk, oldest first. `side` is the discarded fragment CSMI.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CleavageSide {
    pub site: BTreeSet<usize>,
    pub side: String,
    pub opens: Vec<BTreeSet<usize>>,
}

impl CleavageSide {
    pub fn new(
        site: impl IntoIterator<Item = usize>,
        side: impl Into<String>,
        opens: impl IntoIterator<Item = BTreeSet<usize>>,
    ) -> Self {
        Self {
            site: site.into_iter().collect(),
            side: side.into(),
            opens: opens.into_iter().collect(),
        }
    }

    /// Formation site plus prior ring-opens (Python `span_sites`).
    pub fn span_sites(&self) -> Vec<&BTreeSet<usize>> {
        let mut out: Vec<&BTreeSet<usize>> = self.opens.iter().collect();
        out.push(&self.site);
        out
    }
}

/// Uncleared cleavage fragments carried on a [`Deps`] plan. Not a step.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Maybe {
    pub entries: Vec<CleavageSide>,
}

impl Maybe {
    pub fn new(entries: impl IntoIterator<Item = CleavageSide>) -> Self {
        Self {
            entries: entries.into_iter().collect(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn sides(&self) -> Vec<&str> {
        self.entries.iter().map(|e| e.side.as_str()).collect()
    }

    /// True when `side` is a discarded fragment, or `site` overlaps a bag span.
    ///
    /// The bifurcating cleavage site itself does not pass: that step is already
    /// in the required plan. A different reaction on overlapping atoms does.
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        if self.entries.is_empty() {
            return false;
        }
        if let Some(want) = side {
            let want = canon_of(want).unwrap_or_else(|_| want.to_string());
            return self.entries.iter().any(|e| e.side == want);
        }
        let Some(keys) = site else {
            return true;
        };
        if keys.is_empty() {
            return false;
        }
        for entry in &self.entries {
            if keys == &entry.site {
                continue;
            }
            for span in entry.span_sites() {
                if !keys.is_disjoint(span) {
                    return true;
                }
            }
        }
        false
    }
}

/// Flat elementary steps plus precedes (transitive reduction) and [`Maybe`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Deps {
    steps: Vec<Step>,
    precedes: Vec<(usize, usize)>,
    /// Discarded cleavage fragments (Python `PathOutcome.maybe`, on the plan).
    maybe: Maybe,
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
        Self {
            steps,
            precedes,
            maybe: Maybe::default(),
        }
    }

    /// Bind will-add notes → added-by, then build precedes (the plan identity).
    pub fn bind(steps: impl IntoIterator<Item = Step>) -> Self {
        bind_deps(steps.into_iter().collect())
    }

    /// Attach cleavage-side bags (Python `Maybe` on the outcome, here on Deps).
    pub fn with_maybe(mut self, maybe: Maybe) -> Self {
        self.maybe = maybe;
        self
    }

    pub fn steps(&self) -> &[Step] {
        &self.steps
    }

    pub fn precedes(&self) -> &[(usize, usize)] {
        &self.precedes
    }

    pub fn maybe(&self) -> &Maybe {
        &self.maybe
    }

    /// Delegate to [`Maybe::allows`] (site overlap or discarded side SMILES).
    pub fn allows(&self, site: Option<&BTreeSet<usize>>, side: Option<&str>) -> bool {
        self.maybe.allows(site, side)
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

    /// Number of topological sorts under precedes.
    ///
    /// Counts via bitmask DP (no materializing [`Linearization`]s). `n > 20`
    /// falls back to enumerating [`Self::linearizations`] (plans that large are
    /// not expected in PhaseOne multipath).
    pub fn n_linearizations(&self) -> usize {
        let n = self.steps.len();
        if n == 0 {
            return 1;
        }
        if n > 20 {
            return self.linearizations().len();
        }
        count_topological_sorts(n, &self.precedes)
    }

    /// True iff `other` admits exactly the same total orders.
    ///
    /// Requires the same multiset of leaf [`Step`]s (exact rule + site notes —
    /// unique-edit already emits one canonical site per class, so orbit-aware
    /// align is not needed), then compares canonical precedes after aligning
    /// indices. Prefer this over `==` when declaration order of steps may
    /// differ; construction already stores the transitive reduction, so `==`
    /// also sees reduced edges when step order matches.
    pub fn same_linearizations(&self, other: &Deps) -> bool {
        let Some(aligned) = align_deps_indices(&self.steps, &other.steps) else {
            return false;
        };
        let mut edges_b: Vec<_> = other
            .precedes
            .iter()
            .map(|&(a, b)| (aligned[a], aligned[b]))
            .collect();
        edges_b.sort_unstable();
        self.precedes == edges_b
    }

    /// Yield-key when site indices remapped across free-step reorderings: same
    /// rule multiset, same Maybe side CSMI multiset, and precedes isomorphic
    /// under rule-name alignment. Complements exact [`Self::same_linearizations`]
    /// when Index notes differ only because intermediates renumbered.
    pub fn same_rule_maybe_skeleton(&self, other: &Deps) -> bool {
        if self.steps.len() != other.steps.len() {
            return false;
        }
        let mut sides_a: Vec<_> = self.maybe.sides();
        let mut sides_b: Vec<_> = other.maybe.sides();
        sides_a.sort_unstable();
        sides_b.sort_unstable();
        if sides_a != sides_b {
            return false;
        }
        let Some(aligned) = align_deps_indices_by_rule(&self.steps, &other.steps) else {
            return false;
        };
        let mut edges_b: Vec<_> = other
            .precedes
            .iter()
            .map(|&(a, b)| (aligned[a], aligned[b]))
            .collect();
        edges_b.sort_unstable();
        self.precedes == edges_b
    }

    /// True when `other` is a longer walk that only adds steps beyond `self`
    /// (Maybe sides of `self` ⊆ `other`; each of `self`'s steps matches a
    /// distinct step of `other` by rule name). Drops dominated multipath hits.
    pub fn dominates_extension_of(&self, other: &Deps) -> bool {
        if self.steps.len() >= other.steps.len() {
            return false;
        }
        let mut sides_a: Vec<_> = self.maybe.sides();
        let mut sides_b: Vec<_> = other.maybe.sides();
        sides_a.sort_unstable();
        sides_b.sort_unstable();
        if !sides_a.iter().all(|s| sides_b.contains(s)) {
            return false;
        }
        let mut used = vec![false; other.steps.len()];
        for step in &self.steps {
            let found = other
                .steps
                .iter()
                .enumerate()
                .find_map(|(j, o)| (!used[j] && step.rule == o.rule).then_some(j));
            let Some(j) = found else {
                return false;
            };
            used[j] = true;
        }
        true
    }

    /// Count of shared total orders: `|L(self) ∩ L(other)|`.
    ///
    /// Same step multiset required (else `0`). After aligning indices, the
    /// intersection of topological sorts is the sorts of the **edge union**;
    /// a cycle in that union means empty intersection (`0`). Equal to
    /// [`Self::n_linearizations`] on both sides iff [`Self::same_linearizations`].
    pub fn linearization_overlap(&self, other: &Deps) -> usize {
        let Some(aligned) = align_deps_indices(&self.steps, &other.steps) else {
            return 0;
        };
        let n = self.steps.len();
        if n == 0 {
            return 1;
        }
        let mut edges = self.precedes.clone();
        for &(a, b) in &other.precedes {
            edges.push((aligned[a], aligned[b]));
        }
        let Ok(reduced) = canonical_dependency_edges(n, &edges) else {
            return 0;
        };
        // Same nodes + union edges; maybe does not affect required orders.
        Deps {
            steps: self.steps.clone(),
            precedes: reduced,
            maybe: Maybe::default(),
        }
        .n_linearizations()
    }
}

/// Count topological sorts of a DAG on `0..n` (bitmask DP).
fn count_topological_sorts(n: usize, precedes: &[(usize, usize)]) -> usize {
    debug_assert!(n <= 20);
    let mut preds = vec![0u32; n];
    for &(a, b) in precedes {
        if a < n && b < n {
            preds[b] |= 1u32 << a;
        }
    }
    let full = 1usize << n;
    let mut dp = vec![0u128; full];
    dp[0] = 1;
    for mask in 0..full {
        let ways = dp[mask];
        if ways == 0 {
            continue;
        }
        for (v, pred) in preds.iter().enumerate() {
            let bit = 1usize << v;
            if mask & bit != 0 {
                continue;
            }
            if (mask & *pred as usize) != *pred as usize {
                continue;
            }
            dp[mask | bit] = dp[mask | bit].saturating_add(ways);
        }
    }
    usize::try_from(dp[full - 1]).unwrap_or(usize::MAX)
}

/// Map indices in `steps_b` → indices in `steps_a` by [`Step`] equality.
///
/// `None` if the leaf multisets differ. Duplicate equal steps match greedily.
/// Unique-edit canonical sites make exact equality the right identity — do not
/// widen to [`Step::same_site_class`] here.
pub fn align_deps_indices(steps_a: &[Step], steps_b: &[Step]) -> Option<Vec<usize>> {
    if steps_a.len() != steps_b.len() {
        return None;
    }
    let mut used = vec![false; steps_b.len()];
    // remap[j_in_b] = i_in_a
    let mut remap = vec![0usize; steps_b.len()];
    for (i, step) in steps_a.iter().enumerate() {
        let found = steps_b
            .iter()
            .enumerate()
            .find_map(|(j, other)| (!used[j] && other == step).then_some(j));
        let j = found?;
        used[j] = true;
        remap[j] = i;
    }
    Some(remap)
}

/// Align by rule name only (greedy). Used for remapped-index yield keys.
fn align_deps_indices_by_rule(steps_a: &[Step], steps_b: &[Step]) -> Option<Vec<usize>> {
    if steps_a.len() != steps_b.len() {
        return None;
    }
    let mut used = vec![false; steps_b.len()];
    let mut remap = vec![0usize; steps_b.len()];
    for (i, step) in steps_a.iter().enumerate() {
        let found = steps_b
            .iter()
            .enumerate()
            .find_map(|(j, other)| (!used[j] && other.rule == step.rule).then_some(j));
        let j = found?;
        used[j] = true;
        remap[j] = i;
    }
    Some(remap)
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

/// One bitmask per node: bit `b` set iff `a` must precede `b`.
///
/// Port of Python `transitive_closure_masks`. Nodes are opaque index labels
/// `0..n-1` — this does **not** check that two plans' steps are the same
/// reactions/sites. Align node identity first (see [`Deps::same_linearizations`]).
///
/// `# Errors`
/// Invalid edge, self-loop, or cycle. `n > 128` (bitmask width).
pub fn transitive_closure_masks(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<u128>, &'static str> {
    if n > 128 {
        return Err("n > 128");
    }
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
    let mut closure = vec![0u128; n];
    for &a in topo.iter().rev() {
        for &b in &outgoing[a] {
            // closure[a] |= (1 << b) | closure[b]
            closure[a] |= (1u128 << b) | closure[b];
        }
    }
    Ok(closure)
}

/// Unique minimal dependency edge set (transitive reduction) for this DAG.
///
/// Port of Python `canonical_dependency_edges`. Two DAGs over the **same
/// labeled nodes** `0..n-1` have exactly the same valid orderings iff this
/// returns the same edge list. [`Deps`] stores this form on construction.
pub fn canonical_dependency_edges(
    n: usize,
    edges: &[(usize, usize)],
) -> Result<Vec<(usize, usize)>, &'static str> {
    let closure = transitive_closure_masks(n, edges)?;
    let mut canonical_edges = Vec::new();
    for a in 0..n {
        let descendants = closure[a];
        let mut through = 0u128;
        let mut remaining = descendants;
        while remaining != 0 {
            let bit = remaining & remaining.wrapping_neg();
            let b = bit.trailing_zeros() as usize;
            remaining -= bit;
            through |= closure[b];
        }
        let mut direct = descendants & !through;
        while direct != 0 {
            let bit = direct & direct.wrapping_neg();
            let b = bit.trailing_zeros() as usize;
            direct -= bit;
            canonical_edges.push((a, b));
        }
    }
    canonical_edges.sort_unstable();
    Ok(canonical_edges)
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
                    let mut anchors = anchors.clone();
                    anchors.sort_unstable();
                    site.push(PlanAtom::AddedBy {
                        rule: rule.clone(),
                        anchors,
                    });
                }
                PlanAtom::Index(_) => site.push(item.clone()),
            }
        }
        bound.push(Step::new(step.rule.clone(), site).with_orbit(step.orbit.iter().copied()));
    }
    Deps::new(bound, edges)
}

/// [`Deps::bind`] alias (Python `as_deps`).
pub fn as_deps(steps: impl IntoIterator<Item = Step>) -> Deps {
    Deps::bind(steps)
}

/// Leaf-owned plan expander (Python `ReactionRule.canonical_plan`).
///
/// Returns elementary [`Step`]s named after catalog rules (`Hydroxylation`,
/// `Dehydrogenation`, …). `None` on a [`crate::ruleset::RuleSet`] means
/// identity: one step at the discovery site.
pub type CanonicalPlanFn = fn(
    mol: &Molecule,
    rule_name: &str,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step>;

/// Identity plan: `rule` at the given site atoms.
pub fn identity_plan(rule: impl Into<String>, site: impl IntoIterator<Item = usize>) -> Vec<Step> {
    vec![Step::new(rule, site.into_iter().map(PlanAtom::index))]
}

/// Identity plan with automorphism orbit of the site (from generators).
pub fn identity_plan_with_orbit(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
    orbit: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    let site: Vec<_> = site.into_iter().collect();
    vec![Step::new(rule, site.into_iter().map(PlanAtom::index)).with_orbit(orbit)]
}

/// Compat name.
pub fn identity_canonical_plan(
    rule: impl Into<String>,
    site: impl IntoIterator<Item = usize>,
) -> Vec<Step> {
    identity_plan(rule, site)
}

/// Resolve a leaf's plan hook (or identity).
pub fn steps_for_leaf(
    plan: Option<CanonicalPlanFn>,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    match plan {
        Some(f) => f(mol, rule_name, site_atoms, end_effects),
        None => identity_plan(rule_name, site_atoms.iter().copied()),
    }
}

/// Bound [`Deps`] for one hop.
pub fn plan_for_leaf(
    plan: Option<CanonicalPlanFn>,
    rule_name: &str,
    mol: &Molecule,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Deps {
    Deps::bind(steps_for_leaf(
        plan,
        rule_name,
        mol,
        site_atoms,
        end_effects,
    ))
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

/// Python `QuinoneFormation.canonical_plan`: prep missing oxygens, then DH.
///
/// Steps name existing elementary rules (`Hydroxylation`,
/// `OxidativeDehalogenation`, `Dehydrogenation`). Wired on the QF leaf via
/// [`crate::ruleset::RuleSet::with_canonical_plan`].
pub fn quinone_canonical_plan(
    mol: &Molecule,
    _rule_name: &str,
    site_atoms: &[usize],
    end_effects: Option<&[&Effect]>,
) -> Vec<Step> {
    if let (Some(ends), true) = (end_effects, site_atoms.len() >= 2) {
        let plan = hydroxylation_then_dehydrogenation(mol, ends, site_atoms);
        if !plan.is_empty() {
            return plan;
        }
    }
    identity_plan("Dehydrogenation", site_atoms.iter().copied())
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
        assert!(
            plan[2]
                .site
                .iter()
                .any(|a| matches!(a, PlanAtom::WillAdd { .. }))
        );
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
        assert!(
            deps[2]
                .site
                .iter()
                .any(|a| matches!(a, PlanAtom::AddedBy { .. }))
        );
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
            .metabolize(&mol, |_, _, _| true, |_, _, _| true, true).collect::<Result<Vec<_>, _>>().unwrap();
        assert!(!em.is_empty());
        let deps = Deps::bind(identity_plan("Hydroxylation", [em[0].site]));
        assert!(deps.reaches("CC", "CCO").unwrap());
    }

    #[test]
    fn canonical_edges_drop_transitive() {
        let reduced = canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert_eq!(reduced, vec![(0, 1), (1, 2)]);
    }

    #[test]
    fn transitive_closure_and_canonical() {
        // Python test_transitive_closure_and_canonical
        let closure = transitive_closure_masks(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert!(closure[0] & (1 << 1) != 0 && closure[0] & (1 << 2) != 0);
        assert!(closure[1] & (1 << 2) != 0);
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap(),
            vec![(0, 1), (1, 2)]
        );
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2)]).unwrap(),
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap()
        );
    }

    #[test]
    fn transitive_closure_cycle_raises() {
        assert!(
            transitive_closure_masks(2, &[(0, 1), (1, 0)])
                .unwrap_err()
                .contains("cycle")
        );
        assert!(
            transitive_closure_masks(2, &[(0, 0)])
                .unwrap_err()
                .contains("cycle")
        );
    }

    #[test]
    fn deps_same_linearizations_via_transitive_closure() {
        // Python test_deps_same_linearizations_via_transitive_closure:
        // lin-set identity is canonical edges — not == / raw precedes.
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        let chain = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2)]);
        let with_transitive =
            Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2), (0, 2)]);
        assert_eq!(chain, with_transitive);
        assert_eq!(chain.precedes(), &[(0, 1), (1, 2)]);
        assert_eq!(with_transitive.precedes(), &[(0, 1), (1, 2)]);
        assert!(chain.same_linearizations(&with_transitive));
        assert_eq!(
            canonical_dependency_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap(),
            vec![(0, 1), (1, 2)]
        );

        let flipped = Deps::new([c.clone(), a.clone(), b.clone()], [(1, 2), (2, 0)]); // a≺b≺c
        assert!(chain.same_linearizations(&flipped));
        assert_ne!(chain, flipped); // == is order-sensitive

        let layered = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 2), (1, 2)]);
        assert!(layered.same_linearizations(&Deps::new(
            [a.clone(), b.clone(), c.clone()],
            [(0, 2), (1, 2)]
        )));

        let free_dealk = Deps::new([a.clone(), b.clone(), c.clone()], [(1, 2)]); // only b≺c
        assert!(!free_dealk.same_linearizations(&chain));

        let other_nodes = Deps::new(
            [
                Step::new("X", [PlanAtom::index(0)]),
                Step::new("Y", [PlanAtom::index(1)]),
                Step::new("Z", [PlanAtom::index(2)]),
            ],
            [(0, 1), (1, 2)],
        );
        assert!(!chain.same_linearizations(&other_nodes));
    }

    #[test]
    fn deps_stores_canonical_precedes() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        let d = Deps::new([a, b, c], [(0, 1), (1, 2), (0, 2)]);
        assert_eq!(d.precedes(), &[(0, 1), (1, 2)]);
    }

    #[test]
    fn same_linearizations_two_prep_then_final() {
        let h0 = Step::new("Hydroxylation", [PlanAtom::index(0)]);
        let h3 = Step::new("Hydroxylation", [PlanAtom::index(3)]);
        let dh = Step::new(
            "Dehydrogenation",
            [
                PlanAtom::AddedBy {
                    rule: "Hydroxylation".into(),
                    anchors: vec![0],
                },
                PlanAtom::AddedBy {
                    rule: "Hydroxylation".into(),
                    anchors: vec![3],
                },
            ],
        );
        let layered = Deps::new([h0.clone(), h3.clone(), dh.clone()], [(0, 2), (1, 2)]);
        let swapped = Deps::new([h3, h0, dh], [(0, 2), (1, 2)]);
        assert!(layered.same_linearizations(&swapped));
        assert_eq!(layered.n_linearizations(), 2);
        assert_eq!(layered.linearization_overlap(&swapped), 2);
    }

    #[test]
    fn linearization_overlap_partial_and_conflicting() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let c = Step::new("C", [PlanAtom::index(2)]);
        // Free A∥B ≺ C → 2 orders
        let free = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 2), (1, 2)]);
        // Chain A≺B≺C → 1 order (subset of free)
        let chain = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1), (1, 2)]);
        assert_eq!(free.linearization_overlap(&chain), 1);
        assert_eq!(chain.linearization_overlap(&free), 1);
        assert!(!free.same_linearizations(&chain));

        // Opposite A/B orders → union cycles → 0
        let ab = Deps::new([a.clone(), b.clone(), c.clone()], [(0, 1)]);
        let ba = Deps::new([a.clone(), b.clone(), c.clone()], [(1, 0)]);
        assert_eq!(ab.linearization_overlap(&ba), 0);

        // Different step multiset → 0
        let other = Deps::new(
            [a, b, Step::new("D", [PlanAtom::index(2)])],
            [(0, 2), (1, 2)],
        );
        assert_eq!(free.linearization_overlap(&other), 0);
        assert_eq!(
            Deps::new([], []).linearization_overlap(&Deps::new([], [])),
            1
        );
    }

    #[test]
    fn same_linearizations_duplicate_equal_steps() {
        // Python test_align_duplicate_steps_greedy
        let a1 = Step::new("A", [PlanAtom::index(0)]);
        let a2 = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        let d1 = Deps::new([a1.clone(), a2.clone(), b.clone()], [(0, 2), (1, 2)]);
        let d2 = Deps::new([a2, b, a1], [(0, 1), (2, 1)]);
        assert!(d1.same_linearizations(&d2));
    }

    #[test]
    fn free_nodes_n_linearizations_is_factorial() {
        let steps: Vec<_> = (0..4)
            .map(|i| Step::new(format!("S{i}"), [PlanAtom::index(i)]))
            .collect();
        assert_eq!(Deps::new(steps, []).n_linearizations(), 24);
        assert_eq!(Deps::new([], []).n_linearizations(), 1);
        assert_eq!(
            Deps::new([Step::new("A", [PlanAtom::index(0)])], []).n_linearizations(),
            1
        );
    }

    #[test]
    fn same_rule_maybe_skeleton_collapses_remapped_free_dealks() {
        let a = Step::new("Dealkylation", [PlanAtom::index(0)]);
        let b = Step::new("Dealkylation", [PlanAtom::index(9)]);
        let c = Step::new("Dealkylation", [PlanAtom::index(12)]);
        let d = Step::new("Dealkylation", [PlanAtom::index(0)]);
        let maybe = Maybe::new([
            CleavageSide::new([0], "OC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new([9], "OC", std::iter::empty::<BTreeSet<usize>>()),
        ]);
        let p0 = Deps::new([a, b], []).with_maybe(maybe.clone());
        let p1 = Deps::new([c, d], []).with_maybe(maybe);
        assert!(!p0.same_linearizations(&p1)); // exact sites differ
        assert!(p0.same_rule_maybe_skeleton(&p1));
    }

    #[test]
    fn dominates_extension_drops_longer_same_maybe_walk() {
        let short = Deps::new(
            [
                Step::new("Dealkylation", [PlanAtom::index(8)]),
                Step::new("Dealkylation", [PlanAtom::index(1)]),
            ],
            [],
        )
        .with_maybe(Maybe::new([
            CleavageSide::new([8], "CNC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new(
                [1],
                "c1ccc(C(C)(C)C)cc1",
                std::iter::empty::<BTreeSet<usize>>(),
            ),
        ]));
        let longer = Deps::new(
            [
                Step::new("Dealkylation", [PlanAtom::index(8)]),
                Step::new("Dealkylation", [PlanAtom::index(1)]),
                Step::new("Dehydrogenation", [PlanAtom::index(6)]),
            ],
            [],
        )
        .with_maybe(Maybe::new([
            CleavageSide::new([8], "CNC", std::iter::empty::<BTreeSet<usize>>()),
            CleavageSide::new(
                [1],
                "c1ccc(C(C)(C)C)cc1",
                std::iter::empty::<BTreeSet<usize>>(),
            ),
        ]));
        assert!(short.dominates_extension_of(&longer));
        assert!(!longer.dominates_extension_of(&short));
        assert!(!short.same_linearizations(&longer));
    }

    #[test]
    fn n_linearizations_matches_enumeration_on_layered() {
        let h0 = Step::new("Hydroxylation", [PlanAtom::index(0)]);
        let h3 = Step::new("Hydroxylation", [PlanAtom::index(3)]);
        let dh = Step::new("Dehydrogenation", [PlanAtom::index(0)]);
        let layered = Deps::new([h0, h3, dh], [(0, 2), (1, 2)]);
        assert_eq!(layered.n_linearizations(), layered.linearizations().len());
        assert_eq!(layered.n_linearizations(), 2);
    }

    #[test]
    fn deps_invalid_edge_panics() {
        let a = Step::new("A", [PlanAtom::index(0)]);
        let b = Step::new("B", [PlanAtom::index(1)]);
        assert!(
            std::panic::catch_unwind(|| {
                let _ = Deps::new([a.clone(), b.clone()], [(0, 5)]);
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                let _ = Deps::new([a, b], [(0, 1), (1, 0)]);
            })
            .is_err()
        );
    }
}
