//! Per-atom neighborhoods (n0 / n1 / n2) for a whole molecule, and the
//! alignment diff of two such records.
//!
//! Schema:
//! - [`MoleculeShells`] — every heavy atom: `aromatic`, center `h`, heavy-element
//!   bags at distance 0 / 1 / 2 (`n0` / `n1` / `n2`). Counts are non-negative.
//! - [`align_shells`] — two [`MoleculeShells`] + a reactant→target map → the
//!   **same atom shape** with **deltas** (target − reactant) on aligned atoms,
//!   plus how many heavy atoms sit outside the alignment on each side.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use crate::atom_diff::atom_diff;
use crate::mol::{Molecule, atom_idx, atom_usize};

/// Element → count (absolute ≥ 0) or signed delta.
pub type Shell = BTreeMap<String, i32>;

/// Local environment of one heavy atom: aromatic + H + shells n0/n1/n2.
///
/// Absolute shells use `aromatic` ∈ {0,1} and non-negative bag counts.
/// Aligned deltas use `aromatic` = target−reactant ∈ {−1,0,1} and signed bags.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomNeighborhood {
    /// 0/1 on a molecule; target−reactant (−1/0/1) after alignment.
    pub aromatic: i8,
    /// Hydrogens on this atom (absolute) or target−reactant after alignment.
    pub h: i32,
    /// Distance 0 — the center element.
    pub n0: Shell,
    /// Distance 1 — heavy neighbors.
    pub n1: Shell,
    /// Distance 2 — heavy atoms two bonds away.
    pub n2: Shell,
}

/// Neighborhoods for **all** heavy atoms in one molecule (keyed by atom index).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeShells {
    pub atoms: BTreeMap<usize, AtomNeighborhood>,
}

/// Same atom records as deltas under an alignment, plus unaligned heavy counts.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlignedShells {
    /// Target − reactant neighborhood for each aligned reactant atom.
    pub atoms: BTreeMap<usize, AtomNeighborhood>,
    /// Reactant → target atom index.
    pub alignment: BTreeMap<usize, usize>,
    /// Heavy atoms in the reactant with no image.
    pub unaligned_reactant: usize,
    /// Heavy atoms in the target with no preimage.
    pub unaligned_target: usize,
}

fn hydrogens(mol: &Molecule, idx: usize) -> i32 {
    mol.implicit_hydrogen_count(atom_idx(idx)) as i32
}

fn shell_insert(shell: &mut Shell, mol: &Molecule, idx: usize) {
    let z = mol.atom(atom_idx(idx)).element.atomic_number();
    if z <= 1 {
        return;
    }
    let sym = mol.atom(atom_idx(idx)).element.symbol().to_string();
    *shell.entry(sym).or_insert(0) += 1;
}

/// Heavy-atom bags at distance 0 / 1 / 2 from `center`, plus aromatic + H.
pub fn atom_neighborhood(mol: &Molecule, center: usize) -> AtomNeighborhood {
    let atom = mol.atom(atom_idx(center));
    let mut n0 = Shell::new();
    shell_insert(&mut n0, mol, center);

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
    for (&idx, &d) in &dist {
        if idx == center {
            continue;
        }
        match d {
            1 => shell_insert(&mut n1, mol, idx),
            2 => shell_insert(&mut n2, mol, idx),
            _ => {}
        }
    }

    AtomNeighborhood {
        aromatic: i8::from(atom.aromatic),
        h: hydrogens(mol, center),
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
        h: to.h - from.h,
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

/// MCS primary map, then [`align_shells`] on full-molecule shells.
pub fn aligned_shells(reactant: &Molecule, target: &Molecule) -> AlignedShells {
    let diff = atom_diff(reactant, target);
    align_shells(
        &molecule_shells(reactant),
        &molecule_shells(target),
        &diff.mapping,
    )
}

/// Compact shell for display: `C`, `C2,O`, or signed `C-1,O+1`.
pub fn format_shell(shell: &Shell) -> String {
    if shell.is_empty() {
        return "∅".into();
    }
    shell
        .iter()
        .map(|(el, n)| match *n {
            1 => el.clone(),
            -1 => format!("{el}-1"),
            n if n > 1 => format!("{el}{n}"),
            n => format!("{el}{n:+}"),
        })
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn molecule_shells_covers_all_heavy_atoms() {
        let ethanol = parse_mol("CCO").unwrap();
        let shells = molecule_shells(&ethanol);
        assert_eq!(shells.atoms.len(), 3);
        assert!(shells.atoms.contains_key(&0));
        assert!(shells.atoms.contains_key(&1));
        assert!(shells.atoms.contains_key(&2));
    }

    #[test]
    fn ethane_to_ethene_align_is_h_delta_only() {
        let a = parse_mol("CC").unwrap();
        let b = parse_mol("C=C").unwrap();
        let d = aligned_shells(&a, &b);
        assert_eq!(d.unaligned_reactant, 0);
        assert_eq!(d.unaligned_target, 0);
        assert_eq!(d.atoms.len(), 2);
        for env in d.atoms.values() {
            assert_eq!(env.h, -1);
            assert_eq!(env.aromatic, 0);
            assert!(env.n0.is_empty());
            assert!(env.n1.is_empty());
            assert!(env.n2.is_empty());
        }
    }

    #[test]
    fn ethane_to_ethanol_one_unaligned_oxygen() {
        let a = parse_mol("CC").unwrap();
        let b = parse_mol("CCO").unwrap();
        let d = aligned_shells(&a, &b);
        assert_eq!(d.unaligned_reactant, 0);
        assert_eq!(d.unaligned_target, 1);
        assert_eq!(d.atoms.len(), 2);
    }

    #[test]
    fn anisole_to_phenol_unaligned_methyl() {
        let a = parse_mol("COc1ccccc1").unwrap();
        let b = parse_mol("Oc1ccccc1").unwrap();
        let d = aligned_shells(&a, &b);
        assert_eq!(d.unaligned_reactant, 1);
        assert_eq!(d.unaligned_target, 0);
        let o = d
            .atoms
            .values()
            .find(|e| e.h == 1)
            .expect("phenol O gains H");
        assert_eq!(o.h, 1);
    }
}
