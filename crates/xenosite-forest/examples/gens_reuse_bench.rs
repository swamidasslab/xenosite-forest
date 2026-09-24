//! How much does generator reuse save?
//!
//! cargo run -p xenosite-forest --example gens_reuse_bench --release

use std::time::Instant;

use xenosite_forest::{
    ForestMol, atom_bond_generators, atom_orbit_with_gens, atom_pair_orbit_id,
    atom_pair_orbit_id_with_gens, dehydrogenate_hydroquinone, parse_mol,
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
    let cases = [
        ("benzene", "c1ccccc1", 15usize), // C(6,2) unordered pairs
        ("hydroquinone", "Oc1ccc(O)cc1", 6), // rough: few endpoint pairs
        ("phenol", "Oc1ccccc1", 10),
        ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 20),
    ];

    println!("Generator cost (release)\n");
    println!(
        "{:<14} {:>10} {:>12} {:>8} {:>12} {:>12} {:>8}",
        "mol", "1× gens", "cache hit", "N pairs", "old N×gens", "new 1×gens", "speedup"
    );

    for (name, smiles, n_pairs) in cases {
        let mol = parse_mol(smiles).unwrap();
        let held = ForestMol::parse(smiles).unwrap();
        let gens_once = mean_ns(100, || {
            let _ = atom_bond_generators(&mol);
        });
        let _ = held.atom_bond_generators();
        let hit = mean_ns(5000, || {
            let _ = held.atom_bond_generators();
        });
        let old = gens_once * n_pairs as f64;
        println!(
            "{:<14} {:>10} {:>12} {:>8} {:>12} {:>12} {:>7.0}×",
            name,
            fmt_ns(gens_once),
            fmt_ns(hit),
            n_pairs,
            fmt_ns(old),
            fmt_ns(gens_once),
            n_pairs as f64
        );
    }

    println!("\nBenzene: all 15 unordered pair-orbit ids");
    let mol = parse_mol("c1ccccc1").unwrap();
    let gens = atom_bond_generators(&mol);
    let shared = mean_ns(300, || {
        for i in 0..6 {
            for j in (i + 1)..6 {
                let _ = atom_pair_orbit_id_with_gens(&gens, 6, i, j);
            }
        }
    });
    let per = mean_ns(40, || {
        for i in 0..6 {
            for j in (i + 1)..6 {
                let _ = atom_pair_orbit_id(&mol, i, j);
            }
        }
    });
    println!("  shared gens (current style): {}", fmt_ns(shared));
    println!("  gens per pair (old style):   {}", fmt_ns(per));
    println!("  measured speedup:            {:.0}×", per / shared);

    println!("\nBenzene: 6 atom orbits for Step.orbit fills");
    let held = ForestMol::parse("c1ccccc1").unwrap();
    let gens = held.atom_bond_generators();
    let n = 6usize;
    let mol = parse_mol("c1ccccc1").unwrap();
    let cached = mean_ns(2000, || {
        for a in 0..6 {
            let _ = atom_orbit_with_gens(&gens, n, a);
        }
    });
    let each = mean_ns(60, || {
        for a in 0..6 {
            let g = atom_bond_generators(&mol);
            let _ = atom_orbit_with_gens(&g, n, a);
        }
    });
    println!("  shared gens: {}", fmt_ns(cached));
    println!("  gens each:   {}", fmt_ns(each));
    println!("  speedup:     {:.0}×", each / cached);

    let hq = parse_mol("Oc1ccc(O)cc1").unwrap();
    let dh = mean_ns(80, || {
        let _ = dehydrogenate_hydroquinone(&hq).unwrap();
    });
    let gens_hq = mean_ns(100, || {
        let _ = atom_bond_generators(&hq);
    });
    println!("\nhydroquinone dehydrogenate_hydroquinone: {}", fmt_ns(dh));
    println!("  of which 1× gens alone:                 {}", fmt_ns(gens_hq));
    println!(
        "  gens share of that door:                {:.0}%",
        100.0 * gens_hq / dh
    );
}
