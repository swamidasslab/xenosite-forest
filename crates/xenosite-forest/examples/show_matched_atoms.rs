//! Print matched-atom n0/n1/n2 shells for small reactant→target cases.
use xenosite_forest::{format_shell, matched_atoms, parse_mol};

fn dump(label: &str, reactant_smi: &str, target_smi: &str) {
    let reactant = parse_mol(reactant_smi).unwrap();
    let target = parse_mol(target_smi).unwrap();
    let matched = matched_atoms(&reactant, &target);
    println!("=== {label}  ({reactant_smi} → {target_smi}) ===");
    for a in &matched.atoms {
        let ar = if a.from.aromatic { "ar" } else { "al" };
        let at = if a.to.aromatic { "ar" } else { "al" };
        println!(
            "  r{} → t{}  [{ar}→{at}]  H {}→{}",
            a.reactant, a.target, a.from.h, a.to.h
        );
        println!(
            "    from  n0={}  n1={}  n2={}",
            format_shell(&a.from.n0),
            format_shell(&a.from.n1),
            format_shell(&a.from.n2)
        );
        println!(
            "    to    n0={}  n1={}  n2={}",
            format_shell(&a.to.n0),
            format_shell(&a.to.n1),
            format_shell(&a.to.n2)
        );
    }
    println!();
}

fn main() {
    dump("ethane→ethene", "CC", "C=C");
    dump("ethene→ethane", "C=C", "CC");
    dump("ethane→ethanol", "CC", "CCO");
    dump("anisole→phenol", "COc1ccccc1", "Oc1ccccc1");
    dump("hydroquinone→quinone", "Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1");
}
