//! Compare Rust kekule parents to Python ResonanceMolSupplier / reactant_parent.
use std::collections::BTreeMap;
use xenosite_forest::kekule::{KekuleCache, ensure_kekule_parents, kekule_forms, reactant_parent};
use xenosite_forest::mol::{atom_idx, canon_smiles, parse_mol};
use xenosite_forest::smarts::smarts_matches;

fn dump_forms(label: &str, smi: &str) {
    let mol = parse_mol(smi).unwrap();
    let forms = kekule_forms(&mol).unwrap();
    println!("\n=== Rust kekule_forms | {label} | n={} ===", forms.len());
    for (i, f) in forms.iter().enumerate() {
        let charges: Vec<_> = f
            .atoms()
            .filter(|(_, a)| a.charge != 0)
            .map(|(i, a)| format!("{}:{}{:+}", i.0, a.element, a.charge))
            .collect();
        println!("  [{i}] {}  {:?}", canon_smiles(f), charges);
    }
}

fn dump_reactant_parent(
    label: &str,
    smi: &str,
    smarts: &str,
    smirks: &str,
    map: [(u16, usize); 2],
) {
    let mol = parse_mol(smi).unwrap();
    let mapped: BTreeMap<u16, usize> = map.into_iter().collect();
    let mut cache = KekuleCache::default();
    let parent = reactant_parent(&mol, &mapped, smirks, &mut cache).unwrap();
    let (l, r) = (mapped[&1], mapped[&2]);
    let (_, bond) = parent.bond_between(atom_idx(l), atom_idx(r)).unwrap();
    let charges: Vec<_> = parent
        .atoms()
        .filter(|(i, a)| a.charge != 0 || i.0 <= 4)
        .map(|(i, a)| {
            format!(
                "{}:{} ch={} h={:?}",
                i.0, a.element, a.charge, a.hydrogen_count
            )
        })
        .collect();
    println!(
        "\n=== Rust reactant_parent | {label} map={:?} ===\n  {} bond={:?}\n  {:?}",
        (l, r),
        canon_smiles(&parent),
        bond.order,
        charges
    );
    let _ = smarts;
    let slot = ensure_kekule_parents(&mol, l, r, &mut cache);
    println!("  assignments={}", slot.borrow().assignments.len());
}

fn main() {
    for (label, smi) in [
        ("benzene", "c1ccccc1"),
        ("aminophenol", "Nc1ccc(O)cc1"),
        ("pyridinium", "C[n+]1ccccc1"),
        ("indole", "c1ccc2[nH]ccc2c1"),
        ("APAP", "CC(=O)Nc1ccc(O)cc1"),
        ("nitrobenzene", "[O-][N+](=O)c1ccccc1"),
    ] {
        dump_forms(label, smi);
    }

    dump_reactant_parent(
        "APAP quat",
        "CC(=O)Nc1ccc(O)cc1",
        "[#6H0:1][#7:2]",
        "[#6H0:1][#7:2]>>([*:2].[*:1]-O)",
        [(1, 4), (2, 3)],
    );
    for (label, site) in [
        ("indole 5-6", (5, 6)),
        ("indole 6-7", (6, 7)),
        ("indole 0-1", (0, 1)),
    ] {
        dump_reactant_parent(
            label,
            "c1ccc2[nH]ccc2c1",
            "[#6h:1][#6:2]",
            "[#6h:1][#6:2]>>(O=[*:1].[*:2])",
            [(1, site.0), (2, site.1)],
        );
    }
}
