//! Cross-seam derisk tests: MCS, WASM-shaped public API, SMARTS used by live rules.

use xenosite_forest::{
    ForestMol, apply_smirks_at, canon_of, canon_smiles, dehydrogenate_hydroquinone, hydroxylate,
    parse_mol, smarts_matches, unique_atom_sites, unordered_atom_pair_orbit_sizes,
};

#[test]
fn mcs_ethane_in_propane_covers_two_carbons() {
    use chematic::smarts::find_mcs;
    let ethane = parse_mol("CC").unwrap();
    let propane = parse_mol("CCC").unwrap();
    let query = find_mcs(&[&ethane, &propane]);
    assert!(query.atom_count() >= 2);
}

#[test]
fn anisole_dealkylation_pattern_matches_methyl() {
    let mol = parse_mol("COc1ccccc1").unwrap();
    let hits = smarts_matches(&mol, "[#6H3:1][#7,#8H0,#16:2]").unwrap();
    assert_eq!(hits.len(), 1);
    let products = apply_smirks_at("[C:1][O:2]>>[O:2].[C:1](=O)O", &mol, &hits[0]).unwrap();
    let smiles: Vec<String> = products.iter().map(canon_smiles).collect();
    let phenol = canon_of("Oc1ccccc1").unwrap();
    assert!(
        smiles.iter().any(|s| canon_of(s).unwrap() == phenol),
        "{smiles:?}"
    );
}

#[test]
fn unique_edit_plus_hydroxylate_are_the_public_door() {
    let mol = parse_mol("c1ccccc1").unwrap();
    assert_eq!(unique_atom_sites(&mol, "[#6h1:1]").unwrap().len(), 1);
    assert_eq!(
        hydroxylate(&mol)
            .unwrap()
            .into_iter()
            .map(|s| canon_of(&s).unwrap())
            .collect::<Vec<_>>(),
        vec![canon_of("Oc1ccccc1").unwrap()]
    );
    assert_eq!(unordered_atom_pair_orbit_sizes(&mol), vec![3, 6, 6]);
}

#[test]
fn hydroquinone_pair_edit_is_nonempty() {
    let mol = parse_mol("Oc1ccc(O)cc1").unwrap();
    assert!(!dehydrogenate_hydroquinone(&mol).unwrap().is_empty());
}

#[test]
fn forest_mol_cache_is_per_held_mol() {
    let mol = ForestMol::parse("CCC").unwrap();
    assert!(!mol.xf().has_forest());
    let csmi = mol.xf().csmi();
    assert!(mol.xf().has_forest());
    assert!(std::rc::Rc::ptr_eq(&csmi, &mol.xf().csmi()));
    assert!(mol.copy_mol().shares_cache(&mol));
    assert!(!mol.rw_copy().has_forest());
}
