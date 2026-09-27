//! Diagnostic: Chematic Hückel vs RdkitLike vs RDKit-parity experimental.
//!
//! Forest production already uses `apply_aromaticity_rdkit_parity_experimental`
//! (`mol::aromatize`). On remaining leaf-parity fail substrates the three
//! engines agree with each other and with RDKit — see RUST_PYTHON_PARITY C17.
use chematic::core::BondOrder;
use chematic::perception::{
    AromaticityAlgorithm, apply_aromaticity_ex, apply_aromaticity_rdkit_parity_experimental,
};
use chematic::smiles::parse;
use std::collections::BTreeSet;

fn arom_atoms(mol: &chematic::core::Molecule) -> BTreeSet<u32> {
    mol.atoms()
        .filter(|(_, a)| a.aromatic)
        .map(|(i, _)| i.0)
        .collect()
}

fn arom_bonds(mol: &chematic::core::Molecule) -> BTreeSet<(u32, u32)> {
    mol.bonds()
        .filter(|(_, b)| matches!(b.order, BondOrder::Aromatic))
        .map(|(_, b)| {
            let (x, y) = (b.atom1.0, b.atom2.0);
            if x < y {
                (x, y)
            } else {
                (y, x)
            }
        })
        .collect()
}

fn main() {
    let cases: &[(&str, &str)] = &[
        ("PhNCO", "O=C=Nc1ccccc1"),
        ("PhNCS", "S=C=Nc1ccccc1"),
        ("PhNCN", "N=C=Nc1ccccc1"),
        ("nitrobenzene", "[O-][N+](=O)c1ccccc1"),
        ("olsalazine", "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O"),
        ("cinnoline", "c1ccc2nnccc2c1"),
        ("sulfamethoxazole", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"),
        ("benzene-oxide", "C1=CC2OC2C=C1"),
        ("dimethylaniline", "CN(C)c1ccccc1"),
        ("aziridine-Ph", "c1ccccc1N1CC1"),
        ("Ph2NMe", "c1ccc(N(C)c2ccccc2)cc1"),
        ("dihydroacridine", "c1ccc2c(c1)Nc1ccccc1C2"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("chloramphenicol", "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl"),
        ("nitrosobenzene", "O=Nc1ccccc1"),
        ("dimethylaminophenol", "CN(C)c1ccc(O)cc1"),
    ];

    println!("{:<22} {:^8} {:^8}", "name", "Huckel=E", "RdkitL=E");
    for (name, smi) in cases {
        let raw = parse(smi).expect(smi);
        let huckel = apply_aromaticity_ex(&raw, AromaticityAlgorithm::Huckel);
        let rdkit_like = apply_aromaticity_ex(&raw, AromaticityAlgorithm::RdkitLike);
        let exp = apply_aromaticity_rdkit_parity_experimental(&raw)
            .unwrap_or_else(|e| panic!("{name}: experimental refused: {e:?}"));
        let ha = arom_atoms(&huckel);
        let da = arom_atoms(&rdkit_like);
        let ea = arom_atoms(&exp);
        let hb = arom_bonds(&huckel);
        let db = arom_bonds(&rdkit_like);
        let eb = arom_bonds(&exp);
        let he = if ha == ea && hb == eb { "SAME" } else { "DIFF" };
        let de = if da == ea && db == eb { "SAME" } else { "DIFF" };
        println!("{name:<22} {he:^8} {de:^8}  atoms={ea:?}");
        if he == "DIFF" {
            println!("  Huckel atoms={ha:?}");
        }
        if de == "DIFF" {
            println!("  RdkitLike atoms={da:?}");
        }
    }
}
