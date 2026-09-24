//! Generator-reuse timings on larger scaffolds (not benzene toys).
//!
//! cargo run -p xenosite-forest --example gens_reuse_bench --release

use std::time::Instant;

use xenosite_forest::{
    ForestMol, atom_bond_generators, atom_orbit_with_gens, atom_pair_orbit_id,
    atom_pair_orbit_id_with_gens, parse_mol,
};

fn mean_ns(iters: u32, mut body: impl FnMut()) -> f64 {
    for _ in 0..iters.min(8) {
        body();
    }
    let t0 = Instant::now();
    for _ in 0..iters {
        body();
    }
    t0.elapsed().as_nanos() as f64 / f64::from(iters)
}

fn fmt_ns(ns: f64) -> String {
    if ns >= 1_000_000.0 {
        format!("{:.2} ms", ns / 1e6)
    } else if ns >= 1_000.0 {
        format!("{:.1} µs", ns / 1e3)
    } else {
        format!("{ns:.0} ns")
    }
}

fn main() {
    // Same spirit as find_path_bench LARGER / HARD — heavy enough that gens matter.
    let cases = [
        (
            "tBu-bis-ND",
            "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
            40usize,
        ),
        ("tribenzylamine", "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1", 30),
        ("tetraMeO-biphenyl", "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC", 50),
        ("tetraMeO-naph", "COc1cc(OC)c2c(OC)cc(OC)cc2c1", 45),
        ("macrocycle-ND", "C1CCCCCCNC2CCCC(CC2)NCCCC1", 35),
        (
            "TBA",
            "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
            40,
        ),
    ];

    println!("Generator reuse on larger mols (release)\n");
    println!(
        "{:<18} {:>8} {:>10} {:>12} {:>8} {:>12} {:>12} {:>7}",
        "mol", "atoms", "1× gens", "cache hit", "N", "old N×gens", "new 1×gens", "×"
    );

    for (name, smiles, n_queries) in cases {
        let mol = parse_mol(smiles).unwrap();
        let n_atoms = mol.atom_count();
        let held = ForestMol::parse(smiles).unwrap();
        let gens_once = mean_ns(40, || {
            let _ = atom_bond_generators(&mol);
        });
        let _ = held.atom_bond_generators();
        let hit = mean_ns(2000, || {
            let _ = held.atom_bond_generators();
        });
        println!(
            "{:<18} {:>8} {:>10} {:>12} {:>8} {:>12} {:>12} {:>6.0}×",
            name,
            n_atoms,
            fmt_ns(gens_once),
            fmt_ns(hit),
            n_queries,
            fmt_ns(gens_once * n_queries as f64),
            fmt_ns(gens_once),
            n_queries as f64
        );
    }

    println!("\nPair-orbit ids on tetraMeO-biphenyl (first 20 unordered pairs among heavy atoms 0..10)");
    let mol = parse_mol("COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC").unwrap();
    let gens = atom_bond_generators(&mol);
    let n = mol.atom_count();
    let pairs: Vec<(usize, usize)> = (0..10)
        .flat_map(|i| ((i + 1)..10).map(move |j| (i, j)))
        .take(20)
        .collect();
    let shared = mean_ns(80, || {
        for &(i, j) in &pairs {
            let _ = atom_pair_orbit_id_with_gens(&gens, n, i, j);
        }
    });
    let per = mean_ns(20, || {
        for &(i, j) in &pairs {
            let _ = atom_pair_orbit_id(&mol, i, j);
        }
    });
    println!("  shared gens: {}", fmt_ns(shared));
    println!("  gens/pair:   {}", fmt_ns(per));
    println!("  speedup:     {:.1}×", per / shared);

    println!("\nStep.orbit fills on TBA (8 atom orbits, shared vs each)");
    let smiles = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12";
    let mol = parse_mol(smiles).unwrap();
    let held = ForestMol::parse(smiles).unwrap();
    let gens = held.atom_bond_generators();
    let n = mol.atom_count();
    let cached = mean_ns(400, || {
        for a in 0..8 {
            let _ = atom_orbit_with_gens(&gens, n, a);
        }
    });
    let each = mean_ns(20, || {
        for a in 0..8 {
            let g = atom_bond_generators(&mol);
            let _ = atom_orbit_with_gens(&g, n, a);
        }
    });
    println!("  shared gens: {}", fmt_ns(cached));
    println!("  gens each:   {}", fmt_ns(each));
    println!("  speedup:     {:.0}×", each / cached);
}
