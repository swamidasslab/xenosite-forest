//! Microbench of the chemistry door. Not find_path.
//!
//! cargo run -p xenosite-forest --example door_bench --release

use std::time::Instant;

use xenosite_forest::{
    ForestMol, apply_smirks_at, dehydrogenate_hydroquinone, hydroxylate, parse_mol, smarts_matches,
    unique_atom_sites, unordered_atom_pair_orbit_sizes,
};

fn mean_ns(iters: u32, mut body: impl FnMut()) -> f64 {
    // Warmup
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
        let ms = ns / 1e6;
        format!("{ms:.2} ms")
    } else if ns >= 1_000.0 {
        let us = ns / 1e3;
        format!("{us:.1} µs")
    } else {
        format!("{ns:.0} ns")
    }
}

fn main() {
    let cases = [
        ("benzene", "c1ccccc1"),
        ("phenol", "Oc1ccccc1"),
        ("hydroquinone", "Oc1ccc(O)cc1"),
        ("anisole", "COc1ccccc1"),
        ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ];

    println!("Rust chematic/canonaut door (release)\n");
    println!(
        "{:<14} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "mol", "parse", "csmi", "smarts", "unique", "hydrox", "orbits"
    );

    for (name, smiles) in cases {
        let parse = mean_ns(400, || {
            let _ = parse_mol(smiles).unwrap();
        });
        let mol = parse_mol(smiles).unwrap();
        let held = ForestMol::parse(smiles).unwrap();
        let csmi = mean_ns(800, || {
            held.xf().clear_structure();
            let _ = held.xf().csmi();
        });
        let smarts = mean_ns(400, || {
            let _ = smarts_matches(&mol, "[#6h1:1]").unwrap();
        });
        let unique = mean_ns(400, || {
            let _ = unique_atom_sites(&mol, "[#6h1:1]").unwrap();
        });
        let hydrox = mean_ns(200, || {
            let _ = hydroxylate(&mol).unwrap();
        });
        let orbits = mean_ns(80, || {
            let _ = unordered_atom_pair_orbit_sizes(&mol);
        });
        println!(
            "{:<14} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
            name,
            fmt_ns(parse),
            fmt_ns(csmi),
            fmt_ns(smarts),
            fmt_ns(unique),
            fmt_ns(hydrox),
            fmt_ns(orbits)
        );
    }

    println!("\nHydroxylate breakdown (sites × canon, not a faster graph edit)\n");
    println!(
        "{:<14} {:>6} {:>6} {:>10} {:>10} {:>10}",
        "mol", "sites", "prods", "unique", "hydrox", "per prod"
    );
    for (name, smiles) in cases {
        let mol = parse_mol(smiles).unwrap();
        let n_h1 = unique_atom_sites(&mol, "[#6h1:1]").unwrap().len();
        let n_h2 = unique_atom_sites(&mol, "[#6h2,#6h3:1]").unwrap().len();
        let sites = n_h1 + n_h2;
        let products = hydroxylate(&mol).unwrap();
        let unique_both = mean_ns(200, || {
            let _ = unique_atom_sites(&mol, "[#6h1:1]").unwrap();
            let _ = unique_atom_sites(&mol, "[#6h2,#6h3:1]").unwrap();
        });
        let hydrox = mean_ns(200, || {
            let _ = hydroxylate(&mol).unwrap();
        });
        let n = products.len().max(1) as f64;
        println!(
            "{:<14} {:>6} {:>6} {:>10} {:>10} {:>10}",
            name,
            sites,
            products.len(),
            fmt_ns(unique_both),
            fmt_ns(hydrox),
            fmt_ns(hydrox / n)
        );
    }

    let hq = parse_mol("Oc1ccc(O)cc1").unwrap();
    let dh = mean_ns(80, || {
        let _ = dehydrogenate_hydroquinone(&hq).unwrap();
    });
    println!("\nhydroquinone pair-edit: {}", fmt_ns(dh));

    let anisole = parse_mol("COc1ccccc1").unwrap();
    let hits = smarts_matches(&anisole, "[#6H3:1][#7,#8H0,#16:2]").unwrap();
    let smirks = mean_ns(200, || {
        let _ = apply_smirks_at("[C:1][O:2]>>[O:2].[C:1](=O)O", &anisole, &hits[0]).unwrap();
    });
    println!("anisole dealkylation apply: {}", fmt_ns(smirks));
}
