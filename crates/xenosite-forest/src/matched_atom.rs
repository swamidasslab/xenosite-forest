//! Per-atom neighborhoods (n0 / n1 / n2) for a whole molecule, and the
//! alignment diff of two such records.
//!
//! Schema:
//! - [`MoleculeShells`] — every heavy atom: `aromatic` plus heavy+H element bags
//!   at graph distance 0 / 1 / 2. Hydrogens count as neighbors of their heavy
//!   atom (implicit H included). Example ethane carbon: `n0=C:1`, `n1=C:1 H:3`,
//!   `n2=H:3`.
//! - [`align_shells`] — two [`MoleculeShells`] + a reactant→target map → the
//!   **same atom shape** with **deltas** (target − reactant) on aligned atoms,
//!   plus how many heavy atoms sit outside the alignment on each side.
//! - **Cost** — residual between **projected** and **target**, not |δ| itself
//!   ([`site_shell_cost`]). Sign: δ is always **product − reactant**
//!   ([`align_shells`] / [`edit_shells`]); `projected = current + δ`. Cost is
//!   Σ normalized |projected − target|. Each shell contributes
//!   `L1 / Σ max(|a|,|b|)` ∈ [0,1] ([`shell_norm_l1`]); plus |Δaromatic| ∈ {0,1}
//!   as a dearomatization hint. Per atom ≤ 4 (≤ 3 when aligned n0 matches).
//!   Site atoms matched as a **multiset** (orbit / MCS swap safe). `δ = 0` is
//!   distance now; a complete edit scores 0. Close pairs use the joint site
//!   atom list. Cleaving sites should pass [`site_atoms_with_leave`].
//! - **Site bag** — [`SiteShellBag`]: order-invariant multiset of kept site
//!   deltas plus cleaved/added counts. Cleavage (methyl leave) is first-class.
//! - [`check_site_shell_bags`] — warn or error on mismatch after an edit.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use crate::forest_mol::ForestMol;
use crate::mol::{Molecule, atom_idx, atom_usize};

/// Element → count (absolute ≥ 0) or signed delta. Includes `"H"`.
pub type Shell = BTreeMap<String, i32>;

/// Σ |a[k] − b[k]| over the key union; missing treated as 0.
pub fn shell_l1(a: &Shell, b: &Shell) -> usize {
    let mut keys: BTreeSet<&str> = a.keys().map(String::as_str).collect();
    keys.extend(b.keys().map(String::as_str));
    keys.into_iter()
        .map(|k| {
            (a.get(k).copied().unwrap_or(0) - b.get(k).copied().unwrap_or(0)).unsigned_abs()
                as usize
        })
        .sum()
}

/// Normalized shell distance: `L1 / Σ_k max(|a[k]|, |b[k]|)` ∈ [0,1].
///
/// Empty vs empty is 0. Each of n0/n1/n2 can contribute at most 1 to an atom.
pub fn shell_norm_l1(a: &Shell, b: &Shell) -> f64 {
    let mut keys: BTreeSet<&str> = a.keys().map(String::as_str).collect();
    keys.extend(b.keys().map(String::as_str));
    let mut num = 0i64;
    let mut den = 0i64;
    for k in keys {
        let av = a.get(k).copied().unwrap_or(0);
        let bv = b.get(k).copied().unwrap_or(0);
        num += i64::from((av - bv).unsigned_abs());
        den += i64::from(av.unsigned_abs().max(bv.unsigned_abs()));
    }
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

/// Local environment of one heavy atom: aromatic + shells n0/n1/n2.
///
/// Absolute shells use `aromatic` ∈ {0,1} and non-negative bag counts.
/// Aligned deltas use `aromatic` = target−reactant ∈ {−1,0,1} and signed bags.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct AtomNeighborhood {
    /// 0/1 on a molecule; target−reactant (−1/0/1) after alignment.
    pub aromatic: i8,
    /// Distance 0 — the center heavy element only (not its H).
    pub n0: Shell,
    /// Distance 1 — heavy neighbors + H on the center.
    pub n1: Shell,
    /// Distance 2 — heavies at dist 2 + H on heavies at dist 1.
    pub n2: Shell,
}

impl AtomNeighborhood {
    /// True when aromatic and all shells are zero (aligned atom unchanged).
    pub fn is_unchanged(&self) -> bool {
        self.aromatic == 0 && self.n0.is_empty() && self.n1.is_empty() && self.n2.is_empty()
    }

    /// Walk n0, n1, n2 (in that order).
    pub fn shells(&self) -> [(&str, &Shell); 3] {
        [("n0", &self.n0), ("n1", &self.n1), ("n2", &self.n2)]
    }

    /// Σ |δ| vs `other`: aromatic + each of n0/n1/n2 via [`shell_l1`] (missing = 0).
    pub fn l1(&self, other: &Self) -> usize {
        (i32::from(self.aromatic) - i32::from(other.aromatic)).unsigned_abs() as usize
            + shell_l1(&self.n0, &other.n0)
            + shell_l1(&self.n1, &other.n1)
            + shell_l1(&self.n2, &other.n2)
    }

    /// [`Self::l1`] against an empty neighborhood (cost of a stored delta record).
    pub fn abs_delta(&self) -> usize {
        self.l1(&Self::default())
    }

    /// Normalized distance vs `other`.
    ///
    /// Shells: Σ [`shell_norm_l1`] over n0/n1/n2 (each ∈ [0,1]). When
    /// `dearomatic` is set, also |Δaromatic| ∈ {0,1}. Max 3 or 4 accordingly;
    /// aligned n0 match → max 2 or 3.
    pub fn norm_l1_opts(&self, other: &Self, dearomatic: bool) -> f64 {
        let shells = shell_norm_l1(&self.n0, &other.n0)
            + shell_norm_l1(&self.n1, &other.n1)
            + shell_norm_l1(&self.n2, &other.n2);
        if dearomatic {
            let ar = (i32::from(self.aromatic) - i32::from(other.aromatic)).unsigned_abs() as f64;
            ar + shells
        } else {
            shells
        }
    }

    /// [`Self::norm_l1_opts`] with dearomatization hint on.
    pub fn norm_l1(&self, other: &Self) -> f64 {
        self.norm_l1_opts(other, true)
    }

    /// [`Self::norm_l1`] against empty (full miss on all shells).
    pub fn abs_norm(&self) -> f64 {
        self.norm_l1(&Self::default())
    }

    pub fn abs_norm_opts(&self, dearomatic: bool) -> f64 {
        self.norm_l1_opts(&Self::default(), dearomatic)
    }
}

/// Neighborhoods for **all** heavy atoms in one molecule (keyed by atom index).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeShells {
    pub atoms: BTreeMap<usize, AtomNeighborhood>,
}

/// Same atom records as deltas under an alignment, plus unaligned heavy counts.
///
/// On a full-molecule align, `unaligned_*` are absolute unmatched heavy counts.
/// On a site view ([`AlignedShells::at_sites`]), `unaligned_*` count **site**
/// atoms that are cleaved / added.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlignedShells {
    /// Target − reactant neighborhood for each (kept) aligned reactant atom.
    pub atoms: BTreeMap<usize, AtomNeighborhood>,
    /// Reactant → target atom index (for atoms present in `atoms`).
    pub alignment: BTreeMap<usize, usize>,
    /// Unmatched reactant heavies (absolute), or site cleaved count on a site view.
    pub unaligned_reactant: usize,
    /// Unmatched target heavies (absolute), or site added count on a site view.
    pub unaligned_target: usize,
}

impl AlignedShells {
    /// Drop aligned atoms whose neighborhood delta is all zero.
    pub fn without_unchanged(&self) -> Self {
        let mut atoms = BTreeMap::new();
        let mut alignment = BTreeMap::new();
        for (&r, env) in &self.atoms {
            if env.is_unchanged() {
                continue;
            }
            atoms.insert(r, env.clone());
            if let Some(&t) = self.alignment.get(&r) {
                alignment.insert(r, t);
            }
        }
        Self {
            atoms,
            alignment,
            unaligned_reactant: self.unaligned_reactant,
            unaligned_target: self.unaligned_target,
        }
    }

    /// How many aligned atoms still carry a nonzero neighborhood delta.
    pub fn mismatched_matched(&self) -> usize {
        self.atoms.values().filter(|e| !e.is_unchanged()).count()
    }

    /// Site-scoped view: kept site atoms with nonzero delta, plus how many site
    /// atoms are cleaved (unaligned on reactant).
    ///
    /// Methyl leave: the methyl index is absent from `alignment` → cleaved += 1,
    /// and the heteroatom’s shell delta (C:−1 H:+1) stays in `atoms`.
    pub fn at_sites(&self, site_atoms: &[usize]) -> Self {
        let site: HashSet<usize> = site_atoms.iter().copied().collect();
        let mut atoms = BTreeMap::new();
        let mut alignment = BTreeMap::new();
        let mut cleaved = 0usize;
        for &r in &site {
            if let Some(env) = self.atoms.get(&r) {
                if !env.is_unchanged() {
                    atoms.insert(r, env.clone());
                }
                if let Some(&t) = self.alignment.get(&r) {
                    alignment.insert(r, t);
                }
            } else if !self.alignment.contains_key(&r) {
                cleaved += 1;
            }
        }
        Self {
            atoms,
            alignment,
            unaligned_reactant: cleaved,
            unaligned_target: 0,
        }
    }

    /// Magnitude of stored deltas: Σ |δ| on kept atoms + cleaved + added.
    ///
    /// Site / search cost is [`site_shell_cost`] (Σ |projected − target| with
    /// `projected = current + δ`), not this magnitude.
    pub fn delta_magnitude(&self) -> usize {
        self.atoms
            .values()
            .map(AtomNeighborhood::abs_delta)
            .sum::<usize>()
            + self.unaligned_reactant
            + self.unaligned_target
    }

    /// Deprecated name for [`Self::delta_magnitude`]. Prefer [`site_shell_cost`].
    pub fn cost(&self) -> usize {
        self.delta_magnitude()
    }

    /// Order-invariant site bag (forecast ↔ applied edit).
    pub fn site_bag(&self, site_atoms: &[usize]) -> SiteShellBag {
        SiteShellBag::from_align(self, site_atoms)
    }
}

/// Order-invariant site shell pattern: kept deltas as a multiset + cleaved/added.
///
/// Indices are dropped so unique-edit orbit mates / MCS orientation swaps do not
/// spuriously mismatch. Close pairs must pass the **joint** site atom list so
/// each end’s n1/n2 (which may include the other end) is in the same bag.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SiteShellBag {
    /// Sorted neighborhood deltas for kept site atoms (nonzero only).
    pub kept: Vec<AtomNeighborhood>,
    /// Site atoms cleaved (reactant heavy absent from the alignment).
    pub cleaved: usize,
    /// Site-side additions (target heavies counted into the site view).
    pub added: usize,
}

impl SiteShellBag {
    pub fn from_align(align: &AlignedShells, site_atoms: &[usize]) -> Self {
        let site = align.at_sites(site_atoms);
        let mut kept: Vec<AtomNeighborhood> = site.atoms.into_values().collect();
        kept.sort();
        Self {
            kept,
            cleaved: site.unaligned_reactant,
            added: site.unaligned_target,
        }
    }

    /// Magnitude of the bag’s deltas (not residual cost — use [`site_shell_cost`]).
    pub fn delta_magnitude(&self) -> usize {
        self.kept
            .iter()
            .map(AtomNeighborhood::abs_delta)
            .sum::<usize>()
            + self.cleaved
            + self.added
    }

    /// Deprecated name for [`Self::delta_magnitude`].
    pub fn cost(&self) -> usize {
        self.delta_magnitude()
    }

    /// True when bags match exactly (sorted kept + cleaved/added counts).
    pub fn matches(&self, other: &Self) -> bool {
        self == other
    }

    /// Multiset L1 after greedy matching of kept neighborhoods; plus |Δcleaved|, |Δadded|.
    pub fn l1(&self, other: &Self) -> usize {
        let mut unused: Vec<AtomNeighborhood> = other.kept.clone();
        let mut cost = 0usize;
        for a in &self.kept {
            if let Some((i, _)) = unused.iter().enumerate().min_by_key(|(_, b)| a.l1(b)) {
                let b = unused.swap_remove(i);
                cost += a.l1(&b);
            } else {
                cost += a.abs_delta();
            }
        }
        for b in &unused {
            cost += b.abs_delta();
        }
        cost + self.cleaved.abs_diff(other.cleaved) + self.added.abs_diff(other.added)
    }
}

/// How [`check_site_shell_bags`] reports a mismatch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SiteShellCheck {
    /// `log::warn!` only.
    Warn,
    /// Return [`Err`]([`SiteShellMismatch`]).
    Error,
}

/// Forecast site bag did not match the applied edit’s site bag.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteShellMismatch {
    pub label: String,
    pub forecast: SiteShellBag,
    pub actual: SiteShellBag,
    pub l1: usize,
}

impl std::fmt::Display for SiteShellMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "site shell mismatch [{}]: l1={} forecast_cost={} actual_cost={} cleaved {}→{} kept {}→{}",
            self.label,
            self.l1,
            self.forecast.delta_magnitude(),
            self.actual.delta_magnitude(),
            self.forecast.cleaved,
            self.actual.cleaved,
            self.forecast.kept.len(),
            self.actual.kept.len()
        )
    }
}

impl std::error::Error for SiteShellMismatch {}

/// Compare forecast vs applied site bags. Close pairs: pass **both** ends in
/// `site_atoms` so mutual n1/n2 effects stay in one bag.
pub fn check_site_shell_bags(
    label: &str,
    forecast: &SiteShellBag,
    actual: &SiteShellBag,
    mode: SiteShellCheck,
) -> Result<(), SiteShellMismatch> {
    if forecast.matches(actual) {
        return Ok(());
    }
    let l1 = forecast.l1(actual);
    let miss = SiteShellMismatch {
        label: label.into(),
        forecast: forecast.clone(),
        actual: actual.clone(),
        l1,
    };
    match mode {
        SiteShellCheck::Warn => {
            log::warn!("{miss}");
            Ok(())
        }
        SiteShellCheck::Error => Err(miss),
    }
}

/// Tag-aligned reactant→product shell deltas ([`ForestMol`] labels).
///
/// Cleaved heavies (parent tag absent on child) → `unaligned_reactant`.
/// Added heavies (child tag absent on parent) → `unaligned_target`.
pub fn edit_shells(parent: &ForestMol, child: &ForestMol) -> AlignedShells {
    let mut alignment = BTreeMap::new();
    for i in 0..parent.mol().atom_count() {
        if parent.mol().atom(atom_idx(i)).element.atomic_number() <= 1 {
            continue;
        }
        let Some(tag) = parent.tag_of(i) else {
            continue;
        };
        if let Some(j) = child.index_of(tag) {
            if child.mol().atom(atom_idx(j)).element.atomic_number() > 1 {
                alignment.insert(i, j);
            }
        }
    }
    let parent_shells = molecule_shells(parent.mol());
    let child_shells = molecule_shells(child.mol());
    align_shells(&parent_shells, &child_shells, &alignment)
}

/// Site bag of an applied edit at `site_atoms` (parent indices).
///
/// Use the joint atom list for ResonancePair ends (close pairs share shells).
pub fn edit_site_bag(parent: &ForestMol, child: &ForestMol, site_atoms: &[usize]) -> SiteShellBag {
    let edit = edit_shells(parent, child);
    SiteShellBag::from_align(&edit, site_atoms)
}

/// Forecast site bag from a reactant→target align at `site_atoms`.
pub fn forecast_site_bag(align: &AlignedShells, site_atoms: &[usize]) -> SiteShellBag {
    SiteShellBag::from_align(align, site_atoms)
}

fn shell_add(base: &Shell, delta: &Shell) -> Shell {
    let mut keys: BTreeSet<&str> = base.keys().map(String::as_str).collect();
    keys.extend(delta.keys().map(String::as_str));
    let mut out = Shell::new();
    for key in keys {
        let v = base.get(key).copied().unwrap_or(0) + delta.get(key).copied().unwrap_or(0);
        if v != 0 {
            out.insert(key.to_string(), v);
        }
    }
    out
}

/// Apply a signed neighborhood δ: `projected = current + δ`.
///
/// δ must be **product − reactant** (same sign as [`align_shells`]). Negating δ
/// walks away from the product; cost then compares that wrong projection to
/// target and stays high.
fn apply_neighborhood(current: &AtomNeighborhood, delta: &AtomNeighborhood) -> AtomNeighborhood {
    AtomNeighborhood {
        aromatic: current.aromatic + delta.aromatic,
        n0: shell_add(&current.n0, &delta.n0),
        n1: shell_add(&current.n1, &delta.n1),
        n2: shell_add(&current.n2, &delta.n2),
    }
}

/// Options for [`site_shell_cost_opts`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiteShellCostOpts {
    /// When set, each atom also pays |Δaromatic| ∈ {0,1} between projected and
    /// target (aromatic-bit mismatch). Prefer on when the candidate's resolved
    /// [`crate::pattern::Effect::dearomatizes`] is true (Epoxidation on
    /// aromatic sites, DH/QF ends). Ablation may force it for all sites.
    pub dearomatic: bool,
    /// On-demand view: move every H bag one shell closer and drop n2.
    /// `n1["H"]→n0`, `n2["H"]→n1`; n2 heavies discarded. PatternInfo unchanged.
    pub h_closer_no_n2: bool,
}

/// Shift H one shell inward and clear n2 (on-demand view of full n0/n1/n2).
///
/// Ethane carbon `n0=C:1, n1=C:1 H:3, n2=H:3` → `n0=C:1 H:3, n1=C:1 H:3, n2=∅`.
pub fn neighborhood_h_closer_no_n2(env: &AtomNeighborhood) -> AtomNeighborhood {
    let mut n0 = env.n0.clone();
    let mut n1 = env.n1.clone();
    if let Some(h) = n1.remove("H") {
        *n0.entry("H".into()).or_default() += h;
    }
    if let Some(&h) = env.n2.get("H") {
        *n1.entry("H".into()).or_default() += h;
    }
    n1.retain(|_, v| *v != 0);
    n0.retain(|_, v| *v != 0);
    AtomNeighborhood {
        aromatic: env.aromatic,
        n0,
        n1,
        n2: Shell::new(),
    }
}

/// [`neighborhood_h_closer_no_n2`] on every atom.
pub fn molecule_shells_h_closer_no_n2(shells: &MoleculeShells) -> MoleculeShells {
    MoleculeShells {
        atoms: shells
            .atoms
            .iter()
            .map(|(&i, env)| (i, neighborhood_h_closer_no_n2(env)))
            .collect(),
    }
}

/// [`neighborhood_h_closer_no_n2`] on aligned deltas (linear on bags).
pub fn aligned_shells_h_closer_no_n2(align: &AlignedShells) -> AlignedShells {
    AlignedShells {
        atoms: align
            .atoms
            .iter()
            .map(|(&i, env)| (i, neighborhood_h_closer_no_n2(env)))
            .collect(),
        alignment: align.alignment.clone(),
        unaligned_reactant: align.unaligned_reactant,
        unaligned_target: align.unaligned_target,
    }
}

/// Site cost: Σ normalized |projected − target| at `site_atoms`.
///
/// `projected = current + δ` with δ signed **product − reactant**. The quantity
/// compared is the residual between projected and target — not |δ| alone.
/// Defaults to shells only (no |Δaromatic|). See [`site_shell_cost_opts`].
pub fn site_shell_cost(
    current: &MoleculeShells,
    delta: Option<&AlignedShells>,
    target: &MoleculeShells,
    reactant_to_target: &BTreeMap<usize, usize>,
    site_atoms: &[usize],
) -> f64 {
    site_shell_cost_opts(
        current,
        delta,
        target,
        reactant_to_target,
        site_atoms,
        SiteShellCostOpts::default(),
    )
}

/// [`site_shell_cost`] with explicit [`SiteShellCostOpts`].
///
/// Builds the **projected** site bag (`current + δ`), then measures distance to
/// the **target** bag ([`AtomNeighborhood::norm_l1_opts`], greedy multiset).
/// Close pairs: pass both ends. Cleaving sites should pass site∪leave via
/// [`site_atoms_with_leave`] and list leave indices in `leave_atoms` so those
/// heavies count as unmatched debt until cleaved (no MCS target mate — same
/// view as the gate’s leave fragment). Prefer [`site_shell_cost_best_map`] when
/// multiple MCS maps are available.
pub fn site_shell_cost_opts(
    current: &MoleculeShells,
    delta: Option<&AlignedShells>,
    target: &MoleculeShells,
    reactant_to_target: &BTreeMap<usize, usize>,
    site_atoms: &[usize],
    opts: SiteShellCostOpts,
) -> f64 {
    site_shell_cost_leave(
        current,
        delta,
        target,
        reactant_to_target,
        site_atoms,
        &[],
        opts,
    )
}

/// [`site_shell_cost_opts`] with an explicit leave fragment.
///
/// `leave_atoms` (subset of `site_atoms`) are Effect leave heavies: while still
/// present they contribute projected shells with **no** target mate (debt).
/// After the edit cleaves them they drop from both bags — residual falls when
/// the debt exceeds any kept-site shell change.
pub fn site_shell_cost_leave(
    current: &MoleculeShells,
    delta: Option<&AlignedShells>,
    target: &MoleculeShells,
    reactant_to_target: &BTreeMap<usize, usize>,
    site_atoms: &[usize],
    leave_atoms: &[usize],
    opts: SiteShellCostOpts,
) -> f64 {
    let zero = AtomNeighborhood::default();
    let leave: HashSet<usize> = leave_atoms.iter().copied().collect();
    let mut projected: Vec<AtomNeighborhood> = Vec::new();
    let mut target_envs: Vec<AtomNeighborhood> = Vec::new();

    for &r in site_atoms {
        let Some(cur) = current.atoms.get(&r) else {
            continue;
        };
        // Cleaved in the applied edit: leave/debt gone.
        let cleaved = delta.is_some_and(|d| !d.alignment.contains_key(&r));
        if cleaved {
            continue;
        }

        let projected_env = match delta {
            Some(d) => apply_neighborhood(cur, d.atoms.get(&r).unwrap_or(&zero)),
            None => cur.clone(),
        };
        let projected_env = if opts.h_closer_no_n2 {
            neighborhood_h_closer_no_n2(&projected_env)
        } else {
            projected_env
        };
        projected.push(projected_env);

        // Leave fragment: no target mate (gate discards it). Kept site: pair.
        if leave.contains(&r) {
            continue;
        }
        if let Some(t) = reactant_to_target.get(&r).copied() {
            if let Some(tgt) = target.atoms.get(&t) {
                let tgt = if opts.h_closer_no_n2 {
                    neighborhood_h_closer_no_n2(tgt)
                } else {
                    tgt.clone()
                };
                target_envs.push(tgt);
            }
        }
    }

    neighborhood_bag_norm_l1(&projected, &target_envs, opts.dearomatic)
}

/// [`site_shell_cost_opts`] minimized over reactant→target maps (MCS placements).
///
/// Unique-edit orbit / MCS orientation can assign different target mates to the
/// same site atoms; take the **minimum** residual so a bad orientation does not
/// inflate cost.
pub fn site_shell_cost_best_map(
    current: &MoleculeShells,
    delta: Option<&AlignedShells>,
    target: &MoleculeShells,
    maps: &[BTreeMap<usize, usize>],
    site_atoms: &[usize],
    opts: SiteShellCostOpts,
) -> f64 {
    if maps.is_empty() {
        return 0.0;
    }
    maps.iter()
        .map(|m| site_shell_cost_opts(current, delta, target, m, site_atoms, opts))
        .fold(f64::INFINITY, f64::min)
}

/// Expand `site_atoms` with heavies on the leaving side of a cleavage bond.
///
/// Bond priority (first match wins per site):
/// 1. site ↔ other SMARTS-mapped atoms (`mapped` beyond the site) — the edit
///    bond the pattern actually cleaves (e.g. Dehydration C–O when site_map
///    only lists C);
/// 2. bonded pair already inside `site_atoms`;
/// 3. [`AtomDiff`] `cleavage_bonds` touching the site;
/// 4. any heavy neighbor of a singleton site (open leave).
///
/// Then the smaller side (or the side whose heavy count matches `leave_count`)
/// is unioned in. No-op when no cleavage bond touches the site.
pub fn site_atoms_with_leave(
    mol: &Molecule,
    site_atoms: &[usize],
    leave_count: Option<usize>,
    cleavage_bonds: &BTreeSet<(usize, usize)>,
    mapped: &[usize],
) -> Vec<usize> {
    let mut out: BTreeSet<usize> = site_atoms.iter().copied().collect();
    let site: HashSet<usize> = site_atoms.iter().copied().collect();

    let mut bonds: Vec<(usize, usize)> = Vec::new();
    // 1. Pattern-mapped partners outside the declared site (leave end).
    for &a in site_atoms {
        for &b in mapped {
            if site.contains(&b) || a == b {
                continue;
            }
            if mol.bond_between(atom_idx(a), atom_idx(b)).is_some() {
                bonds.push(bond_key_usize(a, b));
            }
        }
    }
    // 2. Explicit bonded pair inside the site.
    if bonds.is_empty() && site_atoms.len() >= 2 {
        for i in 0..site_atoms.len() {
            for j in (i + 1)..site_atoms.len() {
                let a = site_atoms[i];
                let b = site_atoms[j];
                if mol.bond_between(atom_idx(a), atom_idx(b)).is_some() {
                    bonds.push(bond_key_usize(a, b));
                }
            }
        }
    }
    // 3. MCS cleavage bonds touching the site.
    if bonds.is_empty() {
        for &(a, b) in cleavage_bonds {
            if site.contains(&a) || site.contains(&b) {
                bonds.push((a, b));
            }
        }
    }
    // 4. Fall back: any heavy neighbor of a singleton site (open leave).
    if bonds.is_empty() && site_atoms.len() == 1 {
        let a = site_atoms[0];
        for (nbr, _) in mol.neighbors(atom_idx(a)) {
            let b = atom_usize(nbr);
            if mol.atom(nbr).element.atomic_number() > 1 {
                bonds.push(bond_key_usize(a, b));
            }
        }
    }

    for (a, b) in bonds {
        if let Some(n) = leave_count {
            // Sized leave: flood each side, keep the side with `n` heavies
            // (or the smaller if neither matches exactly).
            let side_a = heavy_side(mol, a, b);
            let side_b = heavy_side(mol, b, a);
            let leave = pick_leave_side(&side_a, &side_b, Some(n));
            out.extend(leave.iter().copied());
        } else {
            // Open leave: fragment size unknown — do **not** flood a ring half.
            // Only the bond ends; heavies that actually leave are those absent
            // from the applied edit alignment (see site_shell_cost_leave callers).
            out.insert(a);
            out.insert(b);
        }
    }
    let mut v: Vec<usize> = out.into_iter().collect();
    v.sort_unstable();
    v
}

fn bond_key_usize(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

fn heavy_side(mol: &Molecule, start: usize, blocked: usize) -> Vec<usize> {
    let mut seen = HashSet::from([start]);
    let mut stack = vec![start];
    while let Some(idx) = stack.pop() {
        for (nbr, _) in mol.neighbors(atom_idx(idx)) {
            let n = atom_usize(nbr);
            if n == blocked || !seen.insert(n) {
                continue;
            }
            if mol.atom(nbr).element.atomic_number() <= 1 {
                continue;
            }
            stack.push(n);
        }
    }
    seen.into_iter()
        .filter(|&i| mol.atom(atom_idx(i)).element.atomic_number() > 1)
        .collect()
}

fn pick_leave_side<'a>(a: &'a [usize], b: &'a [usize], leave_count: Option<usize>) -> &'a [usize] {
    if let Some(n) = leave_count {
        if a.len() == n && b.len() != n {
            return a;
        }
        if b.len() == n && a.len() != n {
            return b;
        }
    }
    if a.len() <= b.len() { a } else { b }
}

/// Greedy multiset normalized L1 between two neighborhood bags.
fn neighborhood_bag_norm_l1(
    a: &[AtomNeighborhood],
    b: &[AtomNeighborhood],
    dearomatic: bool,
) -> f64 {
    let mut unused: Vec<AtomNeighborhood> = b.to_vec();
    let mut cost = 0.0;
    for env in a {
        if let Some((i, _)) = unused.iter().enumerate().min_by(|(_, x), (_, y)| {
            env.norm_l1_opts(x, dearomatic)
                .partial_cmp(&env.norm_l1_opts(y, dearomatic))
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
            let other = unused.swap_remove(i);
            cost += env.norm_l1_opts(&other, dearomatic);
        } else {
            cost += env.abs_norm_opts(dearomatic);
        }
    }
    for other in &unused {
        cost += other.abs_norm_opts(dearomatic);
    }
    cost
}

fn hydrogens(mol: &Molecule, idx: usize) -> i32 {
    mol.implicit_hydrogen_count(atom_idx(idx)) as i32
}

fn add_el(shell: &mut Shell, el: &str, n: i32) {
    if n == 0 {
        return;
    }
    *shell.entry(el.to_string()).or_insert(0) += n;
}

fn add_heavy(shell: &mut Shell, mol: &Molecule, idx: usize) {
    let z = mol.atom(atom_idx(idx)).element.atomic_number();
    if z <= 1 {
        return;
    }
    add_el(shell, mol.atom(atom_idx(idx)).element.symbol(), 1);
}

/// Heavy+H bags at distance 0 / 1 / 2 from heavy `center`, plus aromatic.
pub fn atom_neighborhood(mol: &Molecule, center: usize) -> AtomNeighborhood {
    let atom = mol.atom(atom_idx(center));
    debug_assert!(atom.element.atomic_number() > 1);

    let mut n0 = Shell::new();
    add_heavy(&mut n0, mol, center);

    let mut dist: HashMap<usize, u8> = HashMap::new();
    let mut q = VecDeque::new();
    dist.insert(center, 0);
    q.push_back(center);
    while let Some(u) = q.pop_front() {
        let d = dist[&u];
        if d >= 2 {
            continue;
        }
        for (nbr, _) in mol.neighbors(atom_idx(u)) {
            let v = atom_usize(nbr);
            if mol.atom(nbr).element.atomic_number() <= 1 {
                continue;
            }
            if dist.contains_key(&v) {
                continue;
            }
            dist.insert(v, d + 1);
            q.push_back(v);
        }
    }

    let mut n1 = Shell::new();
    let mut n2 = Shell::new();
    add_el(&mut n1, "H", hydrogens(mol, center));
    for (&idx, &d) in &dist {
        if idx == center {
            continue;
        }
        match d {
            1 => {
                add_heavy(&mut n1, mol, idx);
                add_el(&mut n2, "H", hydrogens(mol, idx));
            }
            2 => {
                add_heavy(&mut n2, mol, idx);
            }
            _ => {}
        }
    }

    AtomNeighborhood {
        aromatic: i8::from(atom.aromatic),
        n0,
        n1,
        n2,
    }
}

fn heavy_indices(mol: &Molecule) -> Vec<usize> {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .map(|(idx, _)| atom_usize(idx))
        .collect()
}

/// [`AtomNeighborhood`] for every heavy atom in `mol`.
pub fn molecule_shells(mol: &Molecule) -> MoleculeShells {
    let mut atoms = BTreeMap::new();
    for idx in heavy_indices(mol) {
        atoms.insert(idx, atom_neighborhood(mol, idx));
    }
    MoleculeShells { atoms }
}

fn shell_sub(to: &Shell, from: &Shell) -> Shell {
    let mut keys: BTreeSet<&str> = to.keys().map(String::as_str).collect();
    keys.extend(from.keys().map(String::as_str));
    let mut out = Shell::new();
    for key in keys {
        let d = to.get(key).copied().unwrap_or(0) - from.get(key).copied().unwrap_or(0);
        if d != 0 {
            out.insert(key.to_string(), d);
        }
    }
    out
}

fn neighborhood_delta(to: &AtomNeighborhood, from: &AtomNeighborhood) -> AtomNeighborhood {
    AtomNeighborhood {
        aromatic: to.aromatic - from.aromatic,
        n0: shell_sub(&to.n0, &from.n0),
        n1: shell_sub(&to.n1, &from.n1),
        n2: shell_sub(&to.n2, &from.n2),
    }
}

/// Align two molecule shells: same per-atom shape, values are deltas, plus
/// unaligned heavy-atom counts on each side.
pub fn align_shells(
    reactant: &MoleculeShells,
    target: &MoleculeShells,
    alignment: &BTreeMap<usize, usize>,
) -> AlignedShells {
    let mut atoms = BTreeMap::new();
    for (&r, &t) in alignment {
        let Some(from) = reactant.atoms.get(&r) else {
            continue;
        };
        let Some(to) = target.atoms.get(&t) else {
            continue;
        };
        atoms.insert(r, neighborhood_delta(to, from));
    }
    let mapped_r: HashSet<usize> = alignment.keys().copied().collect();
    let mapped_t: HashSet<usize> = alignment.values().copied().collect();
    let unaligned_reactant = reactant
        .atoms
        .keys()
        .filter(|i| !mapped_r.contains(i))
        .count();
    let unaligned_target = target
        .atoms
        .keys()
        .filter(|i| !mapped_t.contains(i))
        .count();
    AlignedShells {
        atoms,
        alignment: alignment.clone(),
        unaligned_reactant,
        unaligned_target,
    }
}

/// Site delta forecast with caller-projected unmatched reductions.
pub fn site_delta_forecast(
    current: &AlignedShells,
    site_atoms: &[usize],
    reduce_unaligned_reactant: usize,
    reduce_unaligned_target: usize,
) -> AlignedShells {
    let mut site = current.at_sites(site_atoms);
    site.unaligned_reactant = reduce_unaligned_reactant.min(current.unaligned_reactant);
    site.unaligned_target = reduce_unaligned_target.min(current.unaligned_target);
    site
}

/// Compact shell for display: `C:1`, `C:1 H:3`, or signed `H:-1`.
pub fn format_shell(shell: &Shell) -> String {
    if shell.is_empty() {
        return "∅".into();
    }
    shell
        .iter()
        .map(|(el, n)| format!("{el}:{n}"))
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_diff::atom_diff;
    use crate::forest_mol::ForestMol;
    use crate::mol::{atom_idx, parse_mol};
    use crate::rules::dehydrogenation;
    use crate::rules::{dealkylation, epoxidation, hydroxylation};

    fn site_atoms_cand(c: &crate::candidate::Candidate) -> Vec<usize> {
        let mut atoms: Vec<usize> = c
            .pattern
            .site_map
            .iter()
            .filter_map(|m| c.mapped.get(m).copied())
            .collect();
        if atoms.is_empty() {
            atoms.push(c.site);
        }
        atoms.sort_unstable();
        atoms.dedup();
        atoms
    }

    fn shell_negate(shell: &Shell) -> Shell {
        shell.iter().map(|(k, &v)| (k.clone(), -v)).collect()
    }

    fn negate_neighborhood(env: &AtomNeighborhood) -> AtomNeighborhood {
        AtomNeighborhood {
            aromatic: -env.aromatic,
            n0: shell_negate(&env.n0),
            n1: shell_negate(&env.n1),
            n2: shell_negate(&env.n2),
        }
    }

    /// Flip every atom δ to reactant − product (wrong sign for projection).
    fn negate_aligned(align: &AlignedShells) -> AlignedShells {
        AlignedShells {
            atoms: align
                .atoms
                .iter()
                .map(|(&k, v)| (k, negate_neighborhood(v)))
                .collect(),
            alignment: align.alignment.clone(),
            unaligned_reactant: align.unaligned_reactant,
            unaligned_target: align.unaligned_target,
        }
    }

    #[test]
    fn ethane_carbon_shells_include_h() {
        let ethane = parse_mol("CC").unwrap();
        let shells = molecule_shells(&ethane);
        assert_eq!(shells.atoms.len(), 2);
        for env in shells.atoms.values() {
            assert_eq!(format_shell(&env.n0), "C:1");
            assert_eq!(format_shell(&env.n1), "C:1 H:3");
            assert_eq!(format_shell(&env.n2), "H:3");
        }
    }

    #[test]
    fn shell_l1_missing_is_zero() {
        let mut a = Shell::new();
        a.insert("O".into(), 1);
        a.insert("H".into(), 2);
        let mut b = Shell::new();
        b.insert("H".into(), 3);
        assert_eq!(shell_l1(&a, &b), 2);
    }

    #[test]
    fn h_closer_no_n2_shifts_ethane_carbon() {
        // Ethane C: n0=C:1, n1=C:1 H:3, n2=H:3 → n0=C:1 H:3, n1=C:1 H:3, n2=∅.
        let mol = parse_mol("CC").unwrap();
        let env = atom_neighborhood(&mol, 0);
        assert_eq!(env.n0.get("C").copied(), Some(1));
        assert_eq!(env.n1.get("H").copied(), Some(3));
        assert_eq!(env.n2.get("H").copied(), Some(3));
        let v = neighborhood_h_closer_no_n2(&env);
        assert_eq!(v.n0.get("C").copied(), Some(1));
        assert_eq!(v.n0.get("H").copied(), Some(3));
        assert_eq!(v.n1.get("C").copied(), Some(1));
        assert_eq!(v.n1.get("H").copied(), Some(3));
        assert!(v.n2.is_empty(), "{:?}", v.n2);
    }

    #[test]
    fn shell_norm_l1_unit_interval() {
        let mut a = Shell::new();
        a.insert("C".into(), 1);
        a.insert("H".into(), 3);
        let mut b = Shell::new();
        b.insert("C".into(), 1);
        b.insert("H".into(), 2);
        b.insert("O".into(), 1);
        let n = shell_norm_l1(&a, &b);
        assert!((0.0..1.0).contains(&n) || (n - 1.0).abs() < 1e-12, "{n}");
        assert_eq!(shell_norm_l1(&a, &a), 0.0);
        assert_eq!(shell_norm_l1(&Shell::new(), &Shell::new()), 0.0);
        // Full miss vs empty: 1.0
        assert!((shell_norm_l1(&a, &Shell::new()) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn aligned_atom_norm_at_most_three_with_dearomatic() {
        // Same element → n0 matches; n1/n2 + aromatic disagree → ≤ 3.
        let mut left = AtomNeighborhood {
            aromatic: 1,
            ..Default::default()
        };
        left.n0.insert("C".into(), 1);
        left.n1.insert("C".into(), 1);
        left.n1.insert("H".into(), 3);
        let mut right = left.clone();
        right.aromatic = 0;
        right.n1.insert("O".into(), 1);
        *right.n1.get_mut("H").unwrap() = 2;
        let n = left.norm_l1(&right);
        assert!(n <= 3.0 + 1e-12, "{n}");
        assert!(n > 1.0, "dearomatic hint should fire: {n}");
        assert_eq!(shell_norm_l1(&left.n0, &right.n0), 0.0);
        // Pure dearomatic, shells equal → exactly 1.
        let mut same_shells = left.clone();
        same_shells.aromatic = 0;
        assert!((left.norm_l1(&same_shells) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn demethylation_site_bag_matches_applied() {
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let align = aligned_shells_mol(parent.mol(), &target);
        let map = atom_diff(parent.mol(), &target).mapping;
        let set = dealkylation();
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let c = cands
            .iter()
            .find(|c| c.pattern.name.contains("methyl_alcohol"))
            .expect("methyl_alcohol");
        let atoms = site_atoms_cand(c);
        assert!(
            atoms.len() >= 2,
            "demethylation site should be the O–Me bond: {atoms:?}"
        );
        let forecast = forecast_site_bag(&align, &atoms);
        assert!(
            forecast.cleaved >= 1,
            "methyl leave must count as cleaved: {forecast:?}"
        );
        assert!(
            forecast.delta_magnitude() > 0,
            "edit delta magnitude: {forecast:?}"
        );

        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        assert!(
            site_shell_cost(&cur, None, &tgt, &map, &atoms) > 0.0,
            "distance now at demethylation site"
        );
        // Residual align as δ lands on target → cost 0.
        assert!(
            site_shell_cost(&cur, Some(&align), &tgt, &map, &atoms) < 1e-12,
            "perfect residual δ"
        );

        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let actual = edit_site_bag(&parent, &child, &atoms);
        check_site_shell_bags(
            "anisole demethylation",
            &forecast,
            &actual,
            SiteShellCheck::Error,
        )
        .unwrap();
        let edit = edit_shells(&parent, &child);
        assert!(
            site_shell_cost(&cur, Some(&edit), &tgt, &map, &atoms) < 1e-12,
            "applied demethylation lands on target at site"
        );
    }

    #[test]
    fn delta_sign_is_product_minus_reactant() {
        // δ := product − reactant; projected = current + δ; cost = |projected − target|.
        // Negating δ (reactant − product) must not land on target.
        let parent = ForestMol::parse("COc1ccccc1").unwrap();
        let target = parse_mol("Oc1ccccc1").unwrap();
        let map = atom_diff(parent.mol(), &target).mapping;
        let set = dealkylation();
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let c = cands
            .iter()
            .find(|c| c.pattern.name.contains("methyl_alcohol"))
            .expect("methyl_alcohol");
        let atoms = site_atoms_cand(c);
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let edit = edit_shells(&parent, &child);
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);

        let right = site_shell_cost(&cur, Some(&edit), &tgt, &map, &atoms);
        assert!(
            right < 1e-12,
            "product−reactant δ: projected ≈ target ({right})"
        );

        let flipped = negate_aligned(&edit);
        let wrong = site_shell_cost(&cur, Some(&flipped), &tgt, &map, &atoms);
        assert!(
            wrong > 0.5,
            "reactant−product δ must miss target (got {wrong})"
        );

        // Full residual align as δ is also product−reactant → lands at 0.
        let align = aligned_shells_mol(parent.mol(), &target);
        assert!(site_shell_cost(&cur, Some(&align), &tgt, &map, &atoms) < 1e-12);
        let flipped_align = negate_aligned(&align);
        assert!(
            site_shell_cost(&cur, Some(&flipped_align), &tgt, &map, &atoms) > 0.5,
            "negated residual align must miss"
        );
    }

    #[test]
    fn hydroxylation_orbit_bag_matches_applied() {
        let parent = ForestMol::parse("CC").unwrap();
        let target = parse_mol("CCO").unwrap();
        let align = aligned_shells_mol(parent.mol(), &target);
        let map = atom_diff(parent.mol(), &target).mapping;
        let set = hydroxylation();
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let c = &cands[0];
        // Joint orbit: both carbons — edit hydroxylates one; bags still match.
        let mut atoms = site_atoms_cand(c);
        for &i in &c.orbit {
            atoms.push(i);
        }
        atoms.sort_unstable();
        atoms.dedup();
        let forecast = forecast_site_bag(&align, &atoms);
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let actual = edit_site_bag(&parent, &child, &atoms);
        check_site_shell_bags("ethane OH orbit", &forecast, &actual, SiteShellCheck::Error)
            .unwrap();
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        let edit = edit_shells(&parent, &child);
        assert!(
            site_shell_cost(&cur, Some(&edit), &tgt, &map, &atoms) < 1e-12,
            "OH lands on target at joint orbit site"
        );
    }

    #[test]
    fn hydroquinone_pair_joint_bag_matches() {
        let parent = ForestMol::parse("Oc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C=C1").unwrap();
        let align = aligned_shells_mol(parent.mol(), &target);
        let map = atom_diff(parent.mol(), &target).mapping;
        let pairs = dehydrogenation()
            .pair_candidates_leaf(parent.mol())
            .unwrap();
        let pair = &pairs[0];
        let (a, b) = pair.end_atoms().expect("ends");
        // Joint site — each end’s shell sees the other when close.
        let atoms = [a, b];
        let forecast = forecast_site_bag(&align, &atoms);
        let pieces = pair.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let actual = edit_site_bag(&parent, &child, &atoms);
        check_site_shell_bags("HQ DH joint", &forecast, &actual, SiteShellCheck::Error).unwrap();
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        let edit = edit_shells(&parent, &child);
        assert!(
            site_shell_cost(&cur, Some(&edit), &tgt, &map, &atoms) < 1e-12,
            "joint pair δ lands on quinone at both ends"
        );
        // Solo end under-counts mutual shell effects when ends are close.
        let solo = site_shell_cost(&cur, Some(&edit), &tgt, &map, &[a]);
        let joint = site_shell_cost(&cur, Some(&edit), &tgt, &map, &atoms);
        assert!(joint < 1e-12);
        let _ = solo;
    }

    #[test]
    fn epoxidation_aromatic_residual_drops_with_dearomatic() {
        let parent = ForestMol::parse("COc1ccc(O)cc1").unwrap();
        let target = parse_mol("O=C1C=C(O)C(=O)C(O)=C1").unwrap();
        let ad = atom_diff(parent.mol(), &target);
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        let set = epoxidation();
        let c = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
            .into_iter()
            .find(|c| c.site == 4)
            .expect("epoxide at site 4");
        assert!(
            c.pattern.effect.dearomatizes,
            "aromatic MeOPhOH site must resolve dearomatizes"
        );
        let atoms = site_atoms_cand(&c);
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let edit = edit_shells(&parent, &child);
        let opts = SiteShellCostOpts {
            dearomatic: c.pattern.effect.dearomatizes,
            ..SiteShellCostOpts::default()
        };
        let before = site_shell_cost_opts(&cur, None, &tgt, &ad.mapping, &atoms, opts);
        let after = site_shell_cost_opts(&cur, Some(&edit), &tgt, &ad.mapping, &atoms, opts);
        assert!(
            before > after + 1e-12,
            "dearomatic residual should drop: {before:.3} → {after:.3}"
        );
        let shells_only = SiteShellCostOpts::default();
        let before_s = site_shell_cost_opts(&cur, None, &tgt, &ad.mapping, &atoms, shells_only);
        let after_s =
            site_shell_cost_opts(&cur, Some(&edit), &tgt, &ad.mapping, &atoms, shells_only);
        assert!(
            (before_s - after_s).abs() < 1e-9,
            "shells-only stays flat without PatternInfo dearomatizes: {before_s:.3} → {after_s:.3}"
        );
    }

    #[test]
    fn site_atoms_with_leave_includes_methyl() {
        let mol = parse_mol("COc1ccccc1").unwrap();
        let diff = atom_diff(&mol, &parse_mol("Oc1ccccc1").unwrap());
        // O–Me bond: find methyl (C with 1 heavy neighbor) and oxygen.
        let mut me = None;
        let mut oxy = None;
        for i in 0..mol.atom_count() {
            if mol.atom(atom_idx(i)).element.atomic_number() <= 1 {
                continue;
            }
            let sym = mol.atom(atom_idx(i)).element.symbol();
            let heavy_n = mol
                .neighbors(atom_idx(i))
                .filter(|(n, _)| mol.atom(*n).element.atomic_number() > 1)
                .count();
            if sym == "C" && heavy_n == 1 {
                me = Some(i);
            }
            if sym == "O" {
                oxy = Some(i);
            }
        }
        let me = me.expect("methyl");
        let oxy = oxy.expect("oxygen");
        // Unique-edit style: only heteroatom in the seed → leave expands to Me
        // via MCS cleavage_bonds / neighbor fallback.
        let expanded = site_atoms_with_leave(&mol, &[oxy], Some(1), &diff.cleavage_bonds, &[oxy]);
        assert!(expanded.contains(&me), "leave Me missing: {expanded:?}");
        assert!(expanded.contains(&oxy), "{expanded:?}");
        // With both SMARTS maps, prefer the O–Me bond over MCS noise.
        let mapped = site_atoms_with_leave(&mol, &[oxy], Some(1), &diff.cleavage_bonds, &[me, oxy]);
        assert!(mapped.contains(&me) && mapped.contains(&oxy), "{mapped:?}");
        assert_eq!(mapped.len(), 2, "leave_count=1 Me side: {mapped:?}");
    }

    #[test]
    fn dehydration_alcohol_residual_drops_with_leave() {
        use crate::rules::dehydration;
        let parent = ForestMol::parse("COc1ccc(CC=C)cc1O").unwrap();
        let target = parse_mol("O=C1C=CC(=O)C(CC=C)=C1").unwrap();
        let ad = atom_diff(parent.mol(), &target);
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        let set = dehydration();
        let c = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
            .into_iter()
            .find(|c| c.pattern.name == "alcohol")
            .expect("Dehydration alcohol");
        assert!(c.pattern.effect.cleaves);
        assert_eq!(c.pattern.effect.leave_count, Some(1));
        let site: Vec<usize> = c
            .pattern
            .site_map
            .iter()
            .filter_map(|m| c.mapped.get(m).copied())
            .collect();
        let mapped: Vec<usize> = c.mapped.values().copied().collect();
        let atoms = site_atoms_with_leave(
            parent.mol(),
            &site,
            c.pattern.effect.leave_count.map(|n| n as usize),
            &ad.cleavage_bonds,
            &mapped,
        );
        assert!(
            atoms.len() <= 3,
            "leave_count=1 should not pull the whole ring: {atoms:?}"
        );
        let leave_only: Vec<usize> = atoms
            .iter()
            .copied()
            .filter(|a| !site.contains(a))
            .collect();
        assert!(
            !leave_only.is_empty(),
            "expected OH leave beyond site {site:?}, got {atoms:?}"
        );
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let edit = edit_shells(&parent, &child);
        let before = site_shell_cost_leave(
            &cur,
            None,
            &tgt,
            &ad.mapping,
            &atoms,
            &leave_only,
            SiteShellCostOpts::default(),
        );
        let after = site_shell_cost_leave(
            &cur,
            Some(&edit),
            &tgt,
            &ad.mapping,
            &atoms,
            &leave_only,
            SiteShellCostOpts::default(),
        );
        assert!(
            before > after + 1e-12,
            "cleaving OH leave should drop residual: {before:.3} → {after:.3} atoms={atoms:?} leave={leave_only:?}"
        );
    }

    #[test]
    fn open_leave_does_not_flood_ring() {
        use crate::rules::dealkylation;
        let parent = ForestMol::parse("COc1ccc2c(OC)cccc2c1").unwrap();
        let target = parse_mol("O=C1C(=O)c2ccccc2C=C1").unwrap();
        let ad = atom_diff(parent.mol(), &target);
        let cur = molecule_shells(parent.mol());
        let tgt = molecule_shells(&target);
        let set = dealkylation();
        let c = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
            .into_iter()
            .find(|c| c.pattern.name == "cc_quaternary_alcohol" && c.site == 2)
            .expect("cc_quaternary at site 2");
        assert!(c.pattern.effect.cleaves);
        assert!(c.pattern.effect.leave_count.is_none());
        let site: Vec<usize> = c
            .pattern
            .site_map
            .iter()
            .filter_map(|m| c.mapped.get(m).copied())
            .collect();
        let mapped: Vec<usize> = c.mapped.values().copied().collect();
        let expanded = site_atoms_with_leave(
            parent.mol(),
            &site,
            None,
            &ad.cleavage_bonds,
            &mapped,
        );
        assert_eq!(
            expanded.len(),
            2,
            "open leave must stay at bond ends, not flood the ring: {expanded:?}"
        );
        let pieces = c.materialize_mols(parent.mol()).unwrap();
        let child = parent.adopt_product(pieces[0].clone());
        let edit = edit_shells(&parent, &child);
        // Leave = heavies the edit actually drops from the alignment.
        let leave_only: Vec<usize> = expanded
            .iter()
            .copied()
            .filter(|a| !edit.alignment.contains_key(a))
            .collect();
        assert!(
            !leave_only.is_empty(),
            "edit should cleave at least one bond end: align missing from {expanded:?}"
        );
        let mut atoms = site.clone();
        atoms.extend(&leave_only);
        atoms.sort_unstable();
        atoms.dedup();
        let before = site_shell_cost_leave(
            &cur,
            None,
            &tgt,
            &ad.mapping,
            &atoms,
            &leave_only,
            SiteShellCostOpts::default(),
        );
        let after = site_shell_cost_leave(
            &cur,
            Some(&edit),
            &tgt,
            &ad.mapping,
            &atoms,
            &leave_only,
            SiteShellCostOpts::default(),
        );
        assert!(
            before > after + 1e-12,
            "open-leave cleavage should drop residual: {before:.3} → {after:.3} leave={leave_only:?}"
        );
    }

    fn aligned_shells_mol(a: &Molecule, b: &Molecule) -> AlignedShells {
        let diff = atom_diff(a, b);
        align_shells(&molecule_shells(a), &molecule_shells(b), &diff.mapping)
    }
}
