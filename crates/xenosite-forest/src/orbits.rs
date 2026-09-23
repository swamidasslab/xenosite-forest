//! Colored atom+bond nauty graphs via canonaut. Pair orbits for unique-edit.

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};

use canonaut::structs::{CanonautManager, DenseGraph};
use chematic::core::BondOrder;

use crate::mol::{Molecule, atom_usize};

thread_local! {
    static GENERATORS: RefCell<Vec<Vec<u32>>> = const { RefCell::new(Vec::new()) };
}

fn on_automorphism(
    n: u32,
    perm: &mut [u32],
    _orbits: &mut [u32],
    _num_orbits: u32,
    _stabvert: u32,
    _index: u32,
) {
    GENERATORS.with(|slot| slot.borrow_mut().push(perm[..n as usize].to_vec()));
}

fn color_id(colors: &mut HashMap<String, u32>, label: String) -> u32 {
    let next = colors.len() as u32;
    *colors.entry(label).or_insert(next)
}

fn bond_label(order: BondOrder, aromatic: bool) -> &'static str {
    if aromatic {
        return "aromatic";
    }
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => "single",
        BondOrder::Double => "double",
        BondOrder::Triple => "triple",
        BondOrder::Aromatic => "aromatic",
        _ => "other",
    }
}

/// Automorphism generators of the atom+bond-as-vertex colored graph.
///
/// Each generator is `(atom_perm, bond_perm)` with `perm[i] = image of i`.
pub fn atom_bond_generators(mol: &Molecule) -> Vec<(Vec<usize>, Vec<usize>)> {
    let n_atoms = mol.atom_count();
    let n_bonds = mol.bond_count();
    let n_vertices = n_atoms + n_bonds;
    let mut graph = DenseGraph::new(n_vertices);
    let mut color_names = HashMap::new();
    let mut colors = vec![0u32; n_vertices];

    for (idx, atom) in mol.atoms() {
        let stereo = atom.cip_code.map(|code| format!("{code:?}"));
        let label = format!(
            "atom:{}:{}:{:?}:{}:{}:{}:{stereo:?}",
            atom.element.atomic_number(),
            atom.charge,
            atom.isotope,
            atom.aromatic,
            mol.implicit_hydrogen_count(idx),
            atom.hydrogen_count.unwrap_or(0),
        );
        colors[atom_usize(idx)] = color_id(&mut color_names, label);
    }

    for (bond_idx, bond) in mol.bonds() {
        let v = n_atoms + bond_idx.0 as usize;
        let aromatic = bond.order == BondOrder::Aromatic
            || mol.atom(bond.atom1).aromatic && mol.atom(bond.atom2).aromatic;
        let label = format!("bond:{}:{aromatic}", bond_label(bond.order, aromatic));
        colors[v] = color_id(&mut color_names, label);
        graph.add_edge(atom_usize(bond.atom1), v);
        graph.add_edge(atom_usize(bond.atom2), v);
    }

    graph.set_colors(colors);
    GENERATORS.with(|slot| slot.borrow_mut().clear());
    let mut manager = CanonautManager::new(n_vertices)
        .with_canonization()
        .with_automorphism_callback(on_automorphism);
    manager.canonize_graph(&graph);

    GENERATORS.with(|slot| {
        slot.borrow()
            .iter()
            .map(|perm| {
                let atom_map = (0..n_atoms).map(|i| perm[i] as usize).collect();
                let bond_map = (0..n_bonds)
                    .map(|b| perm[n_atoms + b] as usize - n_atoms)
                    .collect();
                (atom_map, bond_map)
            })
            .collect()
    })
}

fn sorted_pair(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

/// Unordered atom–atom orbit partition (benzene ortho/meta/para lock).
pub fn unordered_atom_pair_groups(mol: &Molecule) -> Vec<Vec<(usize, usize)>> {
    let n = mol.atom_count();
    let generators = atom_bond_generators(mol);
    let candidates: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();
    let mut parent: HashMap<(usize, usize), (usize, usize)> =
        candidates.iter().copied().map(|p| (p, p)).collect();

    fn find(
        parent: &mut HashMap<(usize, usize), (usize, usize)>,
        mut item: (usize, usize),
    ) -> (usize, usize) {
        while parent[&item] != item {
            let next = parent[&item];
            parent.insert(item, parent[&next]);
            item = next;
        }
        item
    }

    for pair in &candidates {
        for (atom_map, _) in &generators {
            let image = sorted_pair(atom_map[pair.0], atom_map[pair.1]);
            let a = find(&mut parent, *pair);
            let b = find(&mut parent, image);
            if a != b {
                parent.insert(b, a);
            }
        }
    }

    let mut buckets: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
    for pair in candidates {
        buckets
            .entry(find(&mut parent, pair))
            .or_default()
            .push(pair);
    }
    buckets.into_values().collect()
}

pub fn unordered_atom_pair_orbit_sizes(mol: &Molecule) -> Vec<usize> {
    let mut sizes: Vec<usize> = unordered_atom_pair_groups(mol)
        .into_iter()
        .map(|group| group.len())
        .collect();
    sizes.sort_unstable();
    sizes
}

/// Stable id for an unordered atom pair's nauty orbit.
pub fn atom_pair_orbit_id(mol: &Molecule, left: usize, right: usize) -> usize {
    let want = sorted_pair(left, right);
    for (id, group) in unordered_atom_pair_groups(mol).iter().enumerate() {
        if group.contains(&want) {
            return id;
        }
    }
    usize::MAX
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{parse_mol, ranks};
    use std::collections::HashSet;

    fn graph_distance(mol: &Molecule, start: usize, end: usize) -> usize {
        use std::collections::VecDeque;
        let mut dist = vec![usize::MAX; mol.atom_count()];
        let mut q = VecDeque::new();
        dist[start] = 0;
        q.push_back(start);
        while let Some(node) = q.pop_front() {
            if node == end {
                return dist[node];
            }
            for (nbr, _) in mol.neighbors(crate::mol::atom_idx(node)) {
                let j = atom_usize(nbr);
                if dist[j] == usize::MAX {
                    dist[j] = dist[node] + 1;
                    q.push_back(j);
                }
            }
        }
        usize::MAX
    }

    #[test]
    fn benzene_meta_and_para_share_ranks_but_not_orbits() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let classes = ranks(&mol);
        assert!(classes.iter().all(|&c| c == classes[0]));
        assert_eq!(unordered_atom_pair_orbit_sizes(&mol), vec![3, 6, 6]);

        let ortho = atom_pair_orbit_id(&mol, 0, 1);
        let meta = atom_pair_orbit_id(&mol, 0, 2);
        let para = atom_pair_orbit_id(&mol, 0, 3);
        assert_ne!(meta, para);
        assert_ne!(ortho, meta);
        assert_eq!(atom_pair_orbit_id(&mol, 1, 3), meta);
        assert_eq!(atom_pair_orbit_id(&mol, 3, 0), para);

        let mut by_dist: BTreeMap<usize, HashSet<usize>> = BTreeMap::new();
        for i in 0..6 {
            for j in (i + 1)..6 {
                by_dist
                    .entry(graph_distance(&mol, i, j))
                    .or_default()
                    .insert(atom_pair_orbit_id(&mol, i, j));
            }
        }
        assert_eq!(by_dist[&1].len(), 1);
        assert_eq!(by_dist[&2].len(), 1);
        assert_eq!(by_dist[&3].len(), 1);
        assert_ne!(by_dist[&2], by_dist[&3]);
    }

    #[test]
    fn naphthalene_has_more_than_one_carbon_class() {
        let mol = parse_mol("c1ccc2ccccc2c1").unwrap();
        let classes: HashSet<usize> = ranks(&mol).into_iter().collect();
        assert!(classes.len() >= 2);
        let sizes = unordered_atom_pair_orbit_sizes(&mol);
        assert!(sizes.len() >= 3);
    }
}
