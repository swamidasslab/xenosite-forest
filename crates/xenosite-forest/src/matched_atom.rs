//! Per-matched-atom neighborhoods (n0 / n1 / n2) for MCS pairs.
//!
//! Schema (not a filter branch): each MCS-matched reactant atom carries its
//! local heavy-atom bags at graph distance 0, 1, and 2, plus aromaticity and
//! the center's H count — on **both** the reactant atom and its target image.
//! No cached loss terms; readers compare `from` vs `to`.

use std::collections::{BTreeMap, HashMap, VecDeque};

use crate::atom_diff::atom_diff;
use crate::mol::{Molecule, atom_idx, atom_usize};

/// Heavy-atom element → count at one graph distance from a center.
pub type Shell = BTreeMap<String, usize>;

/// Local environment of one atom: aromatic + H + shells n0/n1/n2.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomNeighborhood {
    pub aromatic: bool,
    /// Hydrogens on this atom (implicit_hydrogen_count).
    pub h: i32,
    /// Distance 0 — the center element (count 1).
    pub n0: Shell,
    /// Distance 1 — heavy neighbors.
    pub n1: Shell,
    /// Distance 2 — heavy atoms two bonds away.
    pub n2: Shell,
}

/// One MCS-matched atom: reactant index, target image, both neighborhoods.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatchedAtom {
    pub reactant: usize,
    pub target: usize,
    pub from: AtomNeighborhood,
    pub to: AtomNeighborhood,
}

/// All matched atoms under one (or the primary) MCS placement.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MatchedAtoms {
    pub atoms: Vec<MatchedAtom>,
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

    // BFS distances among heavy atoms only (H not a graph node here).
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
        aromatic: atom.aromatic,
        h: hydrogens(mol, center),
        n0,
        n1,
        n2,
    }
}

/// Build [`MatchedAtoms`] from an explicit reactant→target mapping.
pub fn matched_atoms_from_mapping(
    reactant: &Molecule,
    target: &Molecule,
    mapping: &BTreeMap<usize, usize>,
) -> MatchedAtoms {
    let mut atoms: Vec<MatchedAtom> = mapping
        .iter()
        .map(|(&r, &t)| MatchedAtom {
            reactant: r,
            target: t,
            from: atom_neighborhood(reactant, r),
            to: atom_neighborhood(target, t),
        })
        .collect();
    atoms.sort_by_key(|a| a.reactant);
    MatchedAtoms { atoms }
}

/// MCS then [`matched_atoms_from_mapping`] on the primary placement.
pub fn matched_atoms(reactant: &Molecule, target: &Molecule) -> MatchedAtoms {
    let diff = atom_diff(reactant, target);
    matched_atoms_from_mapping(reactant, target, &diff.mapping)
}

/// Compact shell for display: `C` or `C2,O`.
pub fn format_shell(shell: &Shell) -> String {
    if shell.is_empty() {
        return "∅".into();
    }
    shell
        .iter()
        .map(|(el, n)| {
            if *n == 1 {
                el.clone()
            } else {
                format!("{el}{n}")
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn ethane_to_ethene_shells() {
        let a = parse_mol("CC").unwrap();
        let b = parse_mol("C=C").unwrap();
        let m = matched_atoms(&a, &b);
        assert_eq!(m.atoms.len(), 2);
        for atom in &m.atoms {
            assert_eq!(format_shell(&atom.from.n0), "C");
            assert_eq!(format_shell(&atom.from.n1), "C");
            assert_eq!(format_shell(&atom.from.n2), "∅");
            assert_eq!(atom.from.h, 3);
            assert_eq!(atom.to.h, 2);
            assert!(!atom.from.aromatic);
            assert!(!atom.to.aromatic);
        }
    }

    #[test]
    fn anisole_methyl_cleaved_not_in_matched() {
        let a = parse_mol("COc1ccccc1").unwrap();
        let b = parse_mol("Oc1ccccc1").unwrap();
        let m = matched_atoms(&a, &b);
        // Methyl carbon is unmapped (cleaved); matched set is O + ring.
        assert!(
            m.atoms
                .iter()
                .all(|x| x.from.n0.contains_key("C") || x.from.n0.contains_key("O"))
        );
        assert!(
            m.atoms
                .iter()
                .any(|x| x.from.n0.contains_key("O") && x.from.h == 0 && x.to.h == 1)
        );
    }
}
