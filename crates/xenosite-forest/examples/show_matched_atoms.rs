//! Print molecule shells and alignment deltas (n0/n1/n2) for small cases.
use xenosite_forest::{
    AlignedShells, MoleculeShells, aligned_shells, format_shell, molecule_shells, parse_mol,
};

fn dump_mol(label: &str, smi: &str, shells: &MoleculeShells) {
    println!("-- {label} ({smi}) all heavy atoms --");
    for (&i, env) in &shells.atoms {
        let ar = if env.aromatic != 0 { "ar" } else { "al" };
        println!(
            "  a{i} [{ar}] H={}  n0={}  n1={}  n2={}",
            env.h,
            format_shell(&env.n0),
            format_shell(&env.n1),
            format_shell(&env.n2)
        );
    }
}

fn dump_aligned(label: &str, reactant_smi: &str, target_smi: &str, d: &AlignedShells) {
    println!(
        "=== {label}  ({reactant_smi} → {target_smi})  unaligned R={} T={} ===",
        d.unaligned_reactant, d.unaligned_target
    );
    for (&r, env) in &d.atoms {
        let t = d.alignment[&r];
        println!(
            "  r{r}→t{t}  aromaticΔ={}  HΔ={:+}  n0={}  n1={}  n2={}",
            env.aromatic,
            env.h,
            format_shell(&env.n0),
            format_shell(&env.n1),
            format_shell(&env.n2)
        );
    }
    println!();
}

fn case(label: &str, reactant_smi: &str, target_smi: &str) {
    let reactant = parse_mol(reactant_smi).unwrap();
    let target = parse_mol(target_smi).unwrap();
    dump_mol("reactant", reactant_smi, &molecule_shells(&reactant));
    dump_mol("target", target_smi, &molecule_shells(&target));
    dump_aligned(
        label,
        reactant_smi,
        target_smi,
        &aligned_shells(&reactant, &target),
    );
}

fn main() {
    case("ethane→ethene", "CC", "C=C");
    case("ethane→ethanol", "CC", "CCO");
    case("anisole→phenol", "COc1ccccc1", "Oc1ccccc1");
}
