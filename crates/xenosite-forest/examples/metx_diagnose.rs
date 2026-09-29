//! Diagnose a few MetX misses: root residual, resolvable?, expand survivors.
use xenosite_forest::{
    atom_diff, canon_of, phase_one, residual_from_diff, residual_resolvable, as_forest_mol,
};

fn diagnose(label: &str, r: &str, p: &str, expect: &str) {
    println!("\n======== {label}");
    println!("expect: {expect}");
    let Ok(rc) = canon_of(r) else { println!("reactant canon FAIL"); return; };
    let Ok(pc) = canon_of(p) else { println!("product canon FAIL"); return; };
    println!("R_csmi={rc}");
    println!("P_csmi={pc}");
    let reactant = as_forest_mol(r).unwrap();
    let target = as_forest_mol(p).unwrap();
    let diff = atom_diff(reactant.mol(), target.mol());
    let res = residual_from_diff(&diff, Some(reactant.mol()), Some(target.mol()));
    println!(
        "root cost={} n_extra={} cleaved={} arom_loss={} bond_raises={} cats={:?}",
        res.cost, res.n_extra, res.cleavage_bonds.len(), res.loses_aromaticity.len(),
        res.bond_raises.len(), res.categories
    );
    let set = phase_one();
    let resolvable = residual_resolvable(&reactant, target.mol(), &diff, &set).unwrap();
    println!("residual_resolvable={resolvable}");
    // Count candidates that could_help
    let mut n = 0usize;
    let mut help = 0usize;
    let mut help_names = Vec::new();
    for c in set.candidates(&reactant) {
        let c = c.unwrap();
        n += 1;
        if c.could_help_on(&diff, Some(target.mol())) {
            help += 1;
            if help_names.len() < 12 {
                help_names.push(format!("{}:{}", c.pattern_name(), c.leaf_rule().unwrap_or("?")));
            }
        }
    }
    println!("candidates={n} could_help={help} sample={help_names:?}");
}

fn main() {
    // curated sample from MetX misses
    diagnose(
        "BIOTID00198 chromenone 'aliphatic OH' (H14→H24)",
        "CCCCOc1ccc2ccc(=O)oc2c1",
        "CC(O)CCOC1CCC2CCC(O)OC2C1",
        "DB product is fully saturated +O — not mono-OH; PhaseOne cannot dearomatize-saturate",
    );
    diagnose(
        "BIOTID00138 benzothiophene 'S-ox' (H8→H18)",
        "CC(=O)c1cc2ccccc2s1",
        "CC(O)C1CC2CCCCC2S1=O",
        "DB product saturated + S=O + carbonyl→CHOH — multi-change",
    );
    diagnose(
        "BIOTID00321 caffeine→'theobromine'",
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "CC1C(O)=NC(=O)c2c1ncn2C",
        "Formula loses N not Me; not real theobromine (C7H8N4O2)",
    );
    diagnose(
        "BIOTID00005 N-OH-IQ",
        "Cn1c(=N)[nH]c2c3cccnc3ccc21",
        "Cn1c(NO)nc2c3cccnc3ccc21",
        "N-hydroxylation: PhaseOne may lack N-OH rule",
    );
    diagnose(
        "BIOTID00966 tacrine→7-OH",
        "N=c1c2c([nH]c3ccccc13)CCCC2",
        "Nc1c2c(nc3c(O)cccc13)CCCC2",
        "Aromatic OH — should be in PhaseOne; diagnose gate/tautomer",
    );
    diagnose(
        "BIOTID00040 thalidomide arene oxide",
        "O=C1N=C(O)CC[C@@H]1N1C(=O)c2ccccc2C1=O",
        "O=C1N=C(O)CC[C@@H]1N1C(=O)C2=CC3OC3C=C2C1=O",
        "Arene epoxidation — check epoxidation rule + tautomer lactam",
    );
    diagnose(
        "BIOTID01270 leflunomide N-O cleavage",
        "Cc1oncc1C(=O)Nc1ccc(C(F)(F)F)cc1",
        r"C/C(O)=C(/C=N)C(=O)Nc1ccc(C(F)(F)F)cc1",
        "Isoxazole ring opening — PhaseOne may lack this cleavage",
    );
    // control: real theobromine spelling
    diagnose(
        "CTRL caffeine→theobromine (chematic-typical)",
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "Cn1cnc2c(=O)[nH]c(=O)n(C)c12",
        "True N-demethylation (−CH2)",
    );
}
