//! Diagnose remaining MetX misses (call-site nostereo).
use chematic::chem::remove_stereo;
use xenosite_forest::{
    as_forest_mol, atom_diff, phase_one, residual_from_diff, residual_resolvable,
};

fn nostereo(s: &str) -> String {
    let fm = as_forest_mol(s).unwrap();
    as_forest_mol(remove_stereo(fm.mol()))
        .unwrap()
        .csmi()
        .as_ref()
        .to_string()
}

fn diagnose(label: &str, r: &str, p: &str, note: &str) {
    println!("\n======== {label}");
    println!("note: {note}");
    let r = nostereo(r);
    let p = nostereo(p);
    println!("R={r}");
    println!("P={p}");
    let reactant = as_forest_mol(&r).unwrap();
    let target = as_forest_mol(&p).unwrap();
    println!(
        "ha R={} P={} Δ={}",
        reactant.heavy_atom_count(),
        target.heavy_atom_count(),
        target.heavy_atom_count() as i32 - reactant.heavy_atom_count() as i32
    );
    let diff = atom_diff(reactant.mol(), target.mol());
    let res = residual_from_diff(&diff, Some(reactant.mol()), Some(target.mol()));
    println!(
        "root cost={} n_extra={} cleaved={} arom_loss={} bond_raises={} cats={:?}",
        res.cost,
        res.n_extra,
        res.cleavage_bonds.len(),
        res.loses_aromaticity.len(),
        res.bond_raises.len(),
        res.categories
    );
    let set = phase_one();
    let resolvable = residual_resolvable(&reactant, target.mol(), &diff, &set).unwrap();
    println!("residual_resolvable={resolvable}");

    // leaf → product hits toward target CSMI
    let want = target.csmi().as_ref().to_string();
    let mut leaf_hits: Vec<String> = Vec::new();
    let mut help_leaves: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    for c in set.candidates(&reactant) {
        let c = c.unwrap();
        let leaf = c.leaf_rule().unwrap_or("?").to_string();
        if c.could_help_on(&diff, Some(target.mol())) {
            *help_leaves.entry(leaf.clone()).or_default() += 1;
        }
        if let Ok(Some(em)) = c.apply() {
            // Emission products
            for prod in em.products.iter() {
                let cs = prod.csmi().as_ref().to_string();
                if cs == want {
                    leaf_hits.push(format!("{}:{}", c.pattern_name(), leaf));
                }
            }
        }
    }
    println!("could_help leaves={help_leaves:?}");
    println!("one-step exact product hits={leaf_hits:?}");
}

fn main() {
    diagnose(
        "00138 S-ox",
        "CC(=O)c1cc2ccccc2s1",
        "CC(O)C1CC2CCCCC2S1=O",
        "labeled S-Oxidation",
    );
    diagnose(
        "00198 chromenone",
        "CCCCOc1ccc2ccc(=O)oc2c1",
        "CC(O)CCOC1CCC2CCC(O)OC2C1",
        "labeled Aliphatic Hydroxylation",
    );
    diagnose(
        "00055 amphetamine",
        "COc1ccccc1C[C@@H](C)N",
        "COC1CCC(O)CC1CC(C)N",
        "labeled Aromatic Hydroxylation",
    );
    diagnose(
        "00966 tacrine amine",
        "N=c1c2c([nH]c3ccccc13)CCCC2",
        "Nc1c2c(nc3c(O)cccc13)CCCC2",
        "ArOH + amine tautomer product",
    );
    diagnose(
        "01270 leflunomide",
        "Cc1oncc1C(=O)Nc1ccc(C(F)(F)F)cc1",
        r"C/C(O)=C(/C=N)C(=O)Nc1ccc(C(F)(F)F)cc1",
        "N-O bond cleavage / isoxazole open",
    );
    diagnose(
        "00154 dichlorobenzonitrile",
        "N#Cc1c(Cl)cccc1Cl",
        "N#CC1C(Cl)=CC=C2OC21Cl",
        "Epoxidation",
    );
    diagnose(
        "00059 TZD",
        "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(COc1ccc(C[C@@H]3SC(=O)N=C3O)cc1)O2",
        "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(COc1ccc(C[C@@H](C(=O)N=CO)S(=O)O)cc1)O2",
        "TZD ring opening",
    );
}
