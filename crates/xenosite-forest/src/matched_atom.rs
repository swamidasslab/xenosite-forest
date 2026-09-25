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
//! - [`site_delta_forecast`] — same shape for one site: **unchanged atoms
//!   dropped**; `unaligned_*` are **projected reductions** (matched-side change
//!   is the remaining atom map).
//! - [`AlignedShells::cost`] — Σ |δ| across aromatic + n0/n1/n2 on kept atoms.
//!   Distinguishes alcohol vs carbonyl addition and cleavage via shell shape.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use crate::mol::{Molecule, atom_idx, atom_usize};

/// Element → count (absolute ≥ 0) or signed delta. Includes `"H"`.
pub type Shell = BTreeMap<String, i32>;

/// Local environment of one heavy atom: aromatic + shells n0/n1/n2.
///
/// Absolute shells use `aromatic` ∈ {0,1} and non-negative bag counts.
/// Aligned deltas use `aromatic` = target−reactant ∈ {−1,0,1} and signed bags.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
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

    /// Σ |δ| of aromatic and all shell bag entries.
    pub fn abs_delta(&self) -> usize {
        let mut c = self.aromatic.unsigned_abs() as usize;
        for shell in [&self.n0, &self.n1, &self.n2] {
            for &v in shell.values() {
                c += v.unsigned_abs() as usize;
            }
        }
        c
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
/// On a [`site_delta_forecast`], unchanged atoms are omitted and `unaligned_*`
/// are **projected reductions** in those unmatched counts (matched-side change
/// is `atoms.len()` nonzero site deltas).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlignedShells {
    /// Target − reactant neighborhood for each (kept) aligned reactant atom.
    pub atoms: BTreeMap<usize, AtomNeighborhood>,
    /// Reactant → target atom index (for atoms present in `atoms`).
    pub alignment: BTreeMap<usize, usize>,
    /// Unmatched reactant heavies (absolute), or projected reduction in a forecast.
    pub unaligned_reactant: usize,
    /// Unmatched target heavies (absolute), or projected reduction in a forecast.
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

    /// Σ |δ| over aromatic + n0/n1/n2 on every kept aligned atom.
    ///
    /// Encodes alcohol vs carbonyl addition (O / H placement in the shells)
    /// and cleavage (shell change on the kept side). Unaligned counts are
    /// separate projected reductions — not part of this sum.
    pub fn cost(&self) -> usize {
        self.atoms.values().map(AtomNeighborhood::abs_delta).sum()
    }
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
///
/// H is not a separate field: center H lands in `n1`; H on a dist-1 heavy
/// lands in `n2`. Explicit H atoms in the mol are ignored as centers (heavy
/// only); their contribution is via the owning heavy's implicit count.
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

/// Site delta forecast: **same shape** as [`AlignedShells`].
///
/// - `atoms` / `alignment`: only site atoms that currently have a nonzero delta
///   (unchanged excluded). That set is the projected matched-side reduction.
/// - `unaligned_reactant` / `unaligned_target`: **projected reductions** in
///   unmatched heavies (not absolute remaining counts).
pub fn site_delta_forecast(
    current: &AlignedShells,
    site_atoms: &[usize],
    reduce_unaligned_reactant: usize,
    reduce_unaligned_target: usize,
) -> AlignedShells {
    let site: HashSet<usize> = site_atoms.iter().copied().collect();
    let mut atoms = BTreeMap::new();
    let mut alignment = BTreeMap::new();
    for (&r, env) in &current.atoms {
        if !site.contains(&r) || env.is_unchanged() {
            continue;
        }
        atoms.insert(r, env.clone());
        if let Some(&t) = current.alignment.get(&r) {
            alignment.insert(r, t);
        }
    }
    AlignedShells {
        atoms,
        alignment,
        unaligned_reactant: reduce_unaligned_reactant.min(current.unaligned_reactant),
        unaligned_target: reduce_unaligned_target.min(current.unaligned_target),
    }
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
    use crate::mol::parse_mol;

    fn aligned(a: &str, b: &str) -> AlignedShells {
        let ra = parse_mol(a).unwrap();
        let rb = parse_mol(b).unwrap();
        let diff = atom_diff(&ra, &rb);
        align_shells(
            &molecule_shells(&ra),
            &molecule_shells(&rb),
            &diff.mapping,
        )
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
            assert_eq!(env.aromatic, 0);
        }
    }

    #[test]
    fn ethene_carbon_shells_include_h() {
        let ethene = parse_mol("C=C").unwrap();
        let shells = molecule_shells(&ethene);
        for env in shells.atoms.values() {
            assert_eq!(format_shell(&env.n0), "C:1");
            assert_eq!(format_shell(&env.n1), "C:1 H:2");
            assert_eq!(format_shell(&env.n2), "H:2");
        }
    }

    #[test]
    fn ethane_to_ethene_h_delta_in_n1_n2() {
        let d = aligned("CC", "C=C");
        assert_eq!(d.unaligned_reactant, 0);
        assert_eq!(d.unaligned_target, 0);
        for env in d.atoms.values() {
            assert_eq!(format_shell(&env.n0), "∅");
            assert_eq!(format_shell(&env.n1), "H:-1");
            assert_eq!(format_shell(&env.n2), "H:-1");
        }
        assert_eq!(d.without_unchanged().atoms.len(), 2);
    }

    #[test]
    fn ethane_to_ethanol_site_forecast_places_oxygen() {
        let d = aligned("CC", "CCO");
        assert_eq!(d.unaligned_target, 1);
        // Hydroxylation on the CH2 carbon (index 1): clears its shell delta and
        // projects unaligned_target −1 (the O).
        let site = site_delta_forecast(&d, &[1], 0, 1);
        assert_eq!(site.unaligned_reactant, 0);
        assert_eq!(site.unaligned_target, 1);
        assert!(site.atoms.contains_key(&1));
        assert!(!site.atoms.contains_key(&0));
        assert!(!site.atoms.values().any(|e| e.is_unchanged()));
    }

    #[test]
    fn anisole_to_phenol_unaligned_methyl() {
        let d = aligned("COc1ccccc1", "Oc1ccccc1");
        assert_eq!(d.unaligned_reactant, 1);
        assert_eq!(d.unaligned_target, 0);
        let o = d
            .atoms
            .values()
            .find(|e| e.n1.get("H") == Some(&1))
            .expect("O n1 H:+1");
        assert_eq!(o.n1.get("H"), Some(&1));
        let cleave = site_delta_forecast(&d, &[0, 1], 1, 0);
        assert_eq!(cleave.unaligned_reactant, 1);
        assert_eq!(cleave.unaligned_target, 0);
    }

    #[test]
    fn cost_encodes_alcohol_vs_carbonyl_vs_cleavage() {
        // Alcohol add: ethane → ethanol. O in n1, one H lost on the CH2.
        let alcohol = aligned("CC", "CCO");
        // Carbonyl add: ethane → acetaldehyde. O in n1, two H lost on that C.
        let carbonyl = aligned("CC", "CC=O");
        // Cleavage: anisole → phenol. Kept O gains H; methyl unaligned.
        let cleave = aligned("COc1ccccc1", "Oc1ccccc1");

        let alcohol_c = alcohol.without_unchanged().cost();
        let carbonyl_c = carbonyl.without_unchanged().cost();
        let cleave_c = cleave.without_unchanged().cost();

        assert!(alcohol_c > 0, "alcohol cost={alcohol_c}");
        assert!(carbonyl_c > 0, "carbonyl cost={carbonyl_c}");
        assert!(cleave_c > 0, "cleave cost={cleave_c}");
        // Carbonyl loses more H on the oxidized carbon than alcohol.
        assert!(
            carbonyl_c > alcohol_c,
            "carbonyl ({carbonyl_c}) should cost more than alcohol ({alcohol_c})"
        );
        // Alcohol site places O on CH2 (atom 1): O:+1 H:−1 in n1.
        let oh_site = alcohol
            .atoms
            .iter()
            .find(|(_, e)| e.n1.get("O") == Some(&1))
            .expect("O appears in n1 at hydroxylation site");
        assert_eq!(oh_site.1.n1.get("H"), Some(&-1));
        // Carbonyl site: O:+1 and H:−2 (aldehyde).
        let co_site = carbonyl
            .atoms
            .iter()
            .find(|(_, e)| e.n1.get("O") == Some(&1))
            .expect("O appears in n1 at carbonyl site");
        assert_eq!(co_site.1.n1.get("H"), Some(&-2));
        // Cleavage: no new O on target; unaligned reactant methyl.
        assert_eq!(cleave.unaligned_reactant, 1);
        assert_eq!(cleave.unaligned_target, 0);
        assert!(
            !cleave
                .atoms
                .values()
                .any(|e| e.n1.get("O") == Some(&1) || e.n0.get("O") == Some(&1)),
            "cleavage must not look like oxygen addition on aligned atoms"
        );
    }
}
