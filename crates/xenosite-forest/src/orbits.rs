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
    _count: u32,
    perm: &mut [u32],
    _orbits: &mut [u32],
    _num_orbits: u32,
    _stabvert: u32,
    n: u32,
) {
    // Canonaut mirrors nauty userautomproc: last argument is n, first is generator count.
    let n = n as usize;
    if perm.len() < n {
        return;
    }
    GENERATORS.with(|slot| slot.borrow_mut().push(perm[..n].to_vec()));
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

/// One automorphism: `(atom_perm, bond_perm)` with `perm[i] = image of i`.
pub type AtomBondGenerator = (Vec<usize>, Vec<usize>);

/// Automorphism generators of the atom+bond-as-vertex colored graph.
pub fn atom_bond_generators(mol: &Molecule) -> Vec<AtomBondGenerator> {
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
    unordered_atom_pair_groups_with_gens(&atom_bond_generators(mol), mol.atom_count())
}

/// Like [`unordered_atom_pair_groups`], with generators supplied (e.g. ForestMol cache).
pub fn unordered_atom_pair_groups_with_gens(
    generators: &[AtomBondGenerator],
    n_atoms: usize,
) -> Vec<Vec<(usize, usize)>> {
    let candidates: Vec<(usize, usize)> = (0..n_atoms)
        .flat_map(|i| ((i + 1)..n_atoms).map(move |j| (i, j)))
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
        for (atom_map, _) in generators {
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
    atom_pair_orbit_id_with_gens(&atom_bond_generators(mol), mol.atom_count(), left, right)
}

/// Like [`atom_pair_orbit_id`], with generators supplied (e.g. ForestMol cache).
pub fn atom_pair_orbit_id_with_gens(
    generators: &[AtomBondGenerator],
    n_atoms: usize,
    left: usize,
    right: usize,
) -> usize {
    let want = sorted_pair(left, right);
    for (id, group) in unordered_atom_pair_groups_with_gens(generators, n_atoms)
        .iter()
        .enumerate()
    {
        if group.contains(&want) {
            return id;
        }
    }
    usize::MAX
}

/// Automorphism orbit of one atom (closure under atom+bond generators).
pub fn atom_orbit(mol: &Molecule, atom: usize) -> Vec<usize> {
    atom_orbit_with_gens(&atom_bond_generators(mol), mol.atom_count(), atom)
}

/// Like [`atom_orbit`], with generators supplied (e.g. ForestMol cache).
pub fn atom_orbit_with_gens(
    generators: &[AtomBondGenerator],
    n_atoms: usize,
    atom: usize,
) -> Vec<usize> {
    if atom >= n_atoms {
        return Vec::new();
    }
    let mut seen = vec![false; n_atoms];
    let mut stack = vec![atom];
    seen[atom] = true;
    while let Some(i) = stack.pop() {
        for (atom_map, _) in generators {
            let j = atom_map[i];
            if j < n_atoms && !seen[j] {
                seen[j] = true;
                stack.push(j);
            }
        }
    }
    (0..n_atoms).filter(|&i| seen[i]).collect()
}

/// Union of automorphism orbits of the given atoms (sorted, deduped).
pub fn atoms_orbit_with_gens(
    generators: &[AtomBondGenerator],
    n_atoms: usize,
    atoms: impl IntoIterator<Item = usize>,
) -> Vec<usize> {
    let mut seen = vec![false; n_atoms];
    for atom in atoms {
        if atom >= n_atoms || seen[atom] {
            continue;
        }
        for i in atom_orbit_with_gens(generators, n_atoms, atom) {
            seen[i] = true;
        }
    }
    (0..n_atoms).filter(|&i| seen[i]).collect()
}

/// Orbit-deduped unordered `k`-subsets of `eligible` site atoms.
///
/// Unique-edit collapses embeddings to one representative per class; ApplyN
/// with `count = k` still needs the distinct **combinations** of sites under
/// Aut(mol). Benzene carbons + `k = 2` → ortho / meta / para (3 reps), not
/// sequential unique-edit alone and not raw `C(6,2)`.
///
/// Each returned vector is a sorted representative `k`-tuple (atom indexes).
pub fn unordered_site_combinations_with_gens(
    generators: &[AtomBondGenerator],
    n_atoms: usize,
    eligible: &[usize],
    k: usize,
) -> Vec<Vec<usize>> {
    let mut eligible: Vec<usize> = eligible
        .iter()
        .copied()
        .filter(|&i| i < n_atoms)
        .collect();
    eligible.sort_unstable();
    eligible.dedup();
    if k == 0 {
        return vec![Vec::new()];
    }
    if k > eligible.len() {
        return Vec::new();
    }
    if k == 1 {
        // One rep per atom orbit intersected with eligible.
        let mut seen = vec![false; n_atoms];
        let mut reps = Vec::new();
        for &a in &eligible {
            if seen[a] {
                continue;
            }
            let orbit = atom_orbit_with_gens(generators, n_atoms, a);
            for &i in &orbit {
                if i < n_atoms {
                    seen[i] = true;
                }
            }
            reps.push(vec![a]);
        }
        return reps;
    }

    let subsets = k_subsets(&eligible, k);
    let mut parent: HashMap<Vec<usize>, Vec<usize>> =
        subsets.iter().cloned().map(|s| (s.clone(), s)).collect();

    fn find(parent: &mut HashMap<Vec<usize>, Vec<usize>>, mut item: Vec<usize>) -> Vec<usize> {
        while parent[&item] != item {
            let next = parent[&item].clone();
            let grand = parent[&next].clone();
            parent.insert(item.clone(), grand);
            item = next;
        }
        item
    }

    for subset in &subsets {
        for (atom_map, _) in generators {
            let mut image: Vec<usize> = subset.iter().map(|&i| atom_map[i]).collect();
            image.sort_unstable();
            // Skip images that left the eligible set.
            if !image.iter().all(|i| eligible.binary_search(i).is_ok()) {
                continue;
            }
            let a = find(&mut parent, subset.clone());
            let b = find(&mut parent, image);
            if a != b {
                parent.insert(b, a);
            }
        }
    }

    let mut buckets: BTreeMap<Vec<usize>, Vec<usize>> = BTreeMap::new();
    for subset in subsets {
        let root = find(&mut parent, subset.clone());
        buckets.entry(root).or_insert(subset);
    }
    buckets.into_values().collect()
}

/// All unordered `k`-subsets of `items` (items must be sorted for stable output).
fn k_subsets(items: &[usize], k: usize) -> Vec<Vec<usize>> {
    let n = items.len();
    if k > n {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        out.push(idx.iter().map(|&i| items[i]).collect());
        // Next combination in colex / std next_combination order.
        let mut i = k;
        while i > 0 && idx[i - 1] == n - k + i - 1 {
            i -= 1;
        }
        if i == 0 {
            break;
        }
        idx[i - 1] += 1;
        for j in i..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
    out
}

/// Convenience: generators from `mol`.
pub fn unordered_site_combinations(
    mol: &Molecule,
    eligible: &[usize],
    k: usize,
) -> Vec<Vec<usize>> {
    unordered_site_combinations_with_gens(
        &atom_bond_generators(mol),
        mol.atom_count(),
        eligible,
        k,
    )
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
    fn benzene_atom_orbit_is_all_six_carbons() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let gens = atom_bond_generators(&mol);
        let orbit = atom_orbit_with_gens(&gens, mol.atom_count(), 0);
        assert_eq!(orbit, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(atom_orbit(&mol, 3), orbit);
    }

    #[test]
    fn propane_methyl_orbit_is_two_ends() {
        let mol = parse_mol("CCC").unwrap();
        let ends = atom_orbit(&mol, 0);
        assert_eq!(ends.len(), 2);
        assert!(ends.contains(&0));
        assert!(!ends.contains(&1)); // middle carbon
    }

    #[test]
    fn naphthalene_has_more_than_one_carbon_class() {
        let mol = parse_mol("c1ccc2ccccc2c1").unwrap();
        let classes: HashSet<usize> = ranks(&mol).into_iter().collect();
        assert!(classes.len() >= 2);
        let sizes = unordered_atom_pair_orbit_sizes(&mol);
        assert!(sizes.len() >= 3);
    }

    #[test]
    fn benzene_site_combinations_k2_are_ortho_meta_para() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let eligible: Vec<_> = (0..6).collect();
        let combos = unordered_site_combinations(&mol, &eligible, 2);
        assert_eq!(combos.len(), 3, "{combos:?}");
        let mut dists: Vec<_> = combos
            .iter()
            .map(|c| graph_distance(&mol, c[0], c[1]))
            .collect();
        dists.sort_unstable();
        assert_eq!(dists, vec![1, 2, 3]);
    }

    #[test]
    fn benzene_site_combinations_k1_is_one_orbit() {
        let mol = parse_mol("c1ccccc1").unwrap();
        let eligible: Vec<_> = (0..6).collect();
        assert_eq!(unordered_site_combinations(&mol, &eligible, 1).len(), 1);
    }

    #[test]
    fn benzene_site_combinations_k3_count() {
        // Unordered triples of benzene carbons up to Aut(D6h): three classes.
        let mol = parse_mol("c1ccccc1").unwrap();
        let eligible: Vec<_> = (0..6).collect();
        let combos = unordered_site_combinations(&mol, &eligible, 3);
        assert_eq!(combos.len(), 3, "{combos:?}");
    }

    #[test]
    fn ethane_site_combinations_k1_and_k2() {
        let mol = parse_mol("CC").unwrap();
        let eligible = vec![0, 1];
        assert_eq!(unordered_site_combinations(&mol, &eligible, 1).len(), 1);
        assert_eq!(unordered_site_combinations(&mol, &eligible, 2).len(), 1);
    }
}
