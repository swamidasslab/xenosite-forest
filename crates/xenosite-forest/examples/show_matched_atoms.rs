//! Focus: how center H sits beside n0/n1/n2 on molecules and after align.
use xenosite_forest::{
    aligned_shells, format_shell, molecule_shells, parse_mol, AlignedShells, MoleculeShells,
};

fn dump_mol(label: &str, smi: &str, shells: &MoleculeShells) {
    println!("-- {label} ({smi}) --");
    for (&i, env) in &shells.atoms {
        let ar = if env.aromatic != 0 { "ar" } else { "al" };
        println!(
            "  a{i} [{ar}]  H={:<2}  n0={:<4} n1={:<6} n2={}",
            env.h,
            format_shell(&env.n0),
            format_shell(&env.n1),
            format_shell(&env.n2)
        );
    }
}

fn dump_aligned(label: &str, a_smi: &str, b_smi: &str, a: &MoleculeShells, b: &MoleculeShells, d: &AlignedShells) {
    println!(
        "=== {label}  ({a_smi} → {b_smi})  unaligned R={} T={} ===",
        d.unaligned_reactant, d.unaligned_target
    );
    println!("  {:>8}  {:>10}  {:>10}  {:>6}  shells Δ", "pair", "H from→to", "HΔ", "arΔ");
    for (&r, env) in &d.atoms {
        let t = d.alignment[&r];
        let hf = a.atoms[&r].h;
        let ht = b.atoms[&t].h;
        println!(
            "  r{r:>2}→t{t:<2}   {hf}→{ht:<7}  {:+6}  {:+4}  n0={} n1={} n2={}",
            env.h,
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
    let sa = molecule_shells(&a);
    let sb = molecule_shells(&b);
    dump_mol("reactant", a_smi, &sa);
    dump_mol("target", b_smi, &sb);
    dump_aligned(label, a_smi, b_smi, &sa, &sb, &aligned_shells(&a, &b));
}

fn main() {
    case("ethane→ethene (lose H)", "CC", "C=C");
    case("ethene→ethane (gain H)", "C=C", "CC");
    case("ethane→ethanol (O + H shift)", "CC", "CCO");
    case("acetaldehyde→ethanol", "CC=O", "CCO");
    case("phenol→quinone-ish H on O", "Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1");
}
