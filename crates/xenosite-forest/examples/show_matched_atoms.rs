//! Print molecule shells and alignment deltas with H inside n0/n1/n2.
use xenosite_forest::{
    AlignedShells, MoleculeShells, aligned_shells, format_shell, molecule_shells, parse_mol,
};

fn dump_mol(label: &str, smi: &str, shells: &MoleculeShells) {
    println!("-- {label} ({smi}) --");
    for (&i, env) in &shells.atoms {
        let ar = if env.aromatic != 0 { "ar" } else { "al" };
        println!(
            "  a{i} [{ar}]  n0={{{}}}  n1={{{}}}  n2={{{}}}",
            format_shell(&env.n0),
            format_shell(&env.n1),
            format_shell(&env.n2)
        );
    }
}

fn dump_aligned(label: &str, a_smi: &str, b_smi: &str, d: &AlignedShells) {
    println!(
        "=== {label}  ({a_smi} → {b_smi})  unaligned R={} T={} ===",
        d.unaligned_reactant, d.unaligned_target
    );
    for (&r, env) in &d.atoms {
        let t = d.alignment[&r];
        println!(
            "  r{r}→t{t}  arΔ={}  n0={{{}}}  n1={{{}}}  n2={{{}}}",
            env.aromatic,
            format_shell(&env.n0),
            format_shell(&env.n1),
            format_shell(&env.n2)
        );
    }
    println!();
}

fn case(label: &str, a_smi: &str, b_smi: &str) {
    let a = parse_mol(a_smi).unwrap();
    let b = parse_mol(b_smi).unwrap();
    dump_mol("reactant", a_smi, &molecule_shells(&a));
    dump_mol("target", b_smi, &molecule_shells(&b));
    dump_aligned(label, a_smi, b_smi, &aligned_shells(&a, &b));
}

fn main() {
    println!("=== ethane alone ===");
    dump_mol("ethane", "CC", &molecule_shells(&parse_mol("CC").unwrap()));
    println!();
    case("ethane→ethene", "CC", "C=C");
    case("ethene→ethane", "C=C", "CC");
    case("ethane→ethanol", "CC", "CCO");
    case("acetaldehyde→ethanol", "CC=O", "CCO");
}
