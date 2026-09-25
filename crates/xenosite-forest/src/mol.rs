//! Parse, aromatize (RDKit-parity), and write canonical SMILES.

use chematic::core::{AtomIdx, Molecule as CoreMolecule};
use chematic::perception::apply_aromaticity_rdkit_parity_experimental;
use chematic::smiles::{
    canonical_smiles, canonical_smiles_stable_key, parse, topological_equivalence_classes,
};

pub use chematic::core::Molecule;

/// Convert a chematic atom index to `usize`.
pub fn atom_usize(idx: AtomIdx) -> usize {
    idx.0 as usize
}

/// Convert a `usize` site index to a chematic atom index.
pub fn atom_idx(i: usize) -> AtomIdx {
    AtomIdx(i as u32)
}

#[derive(Debug)]
pub enum ForestError {
    Parse(String),
    Smarts(String),
    Smirks(String),
    Kekule(String),
    Plan(String),
}

impl std::fmt::Display for ForestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse(msg)
            | Self::Smarts(msg)
            | Self::Smirks(msg)
            | Self::Kekule(msg)
            | Self::Plan(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for ForestError {}

/// Apply RDKit-parity aromaticity. On refusal, return the input unchanged.
pub fn aromatize(mol: &CoreMolecule) -> CoreMolecule {
    match apply_aromaticity_rdkit_parity_experimental(mol) {
        Ok(aromatized) => aromatized,
        Err(_) => mol.clone(),
    }
}

/// Parse SMILES and apply RDKit-parity aromaticity.
///
/// Kekulé input is aromatized. Already-aromatic SMILES keep atom indexes.
/// If the parity engine refuses a structure, the parsed mol is kept.
pub fn parse_mol(smiles: &str) -> Result<CoreMolecule, ForestError> {
    let mol = parse(smiles).map_err(|err| ForestError::Parse(err.to_string()))?;
    Ok(aromatize(&mol))
}

pub fn canon_smiles(mol: &CoreMolecule) -> String {
    let raw = canonical_smiles(mol);
    match parse_mol(&raw) {
        Ok(reparsed) => canonical_smiles(&reparsed),
        Err(_) => raw,
    }
}

/// Fail-closed identity key for dedup / caches (Chematic docs).
///
/// [`canon_smiles`] is fine for display and for comparing to a known target
/// spelling. It is **not** always safe as a HashSet / yield key: coupled E/Z
/// systems can emit a non-idempotent canonical string. This wraps
/// [`canonical_smiles_stable_key`] — `None` means do not index that molecule
/// by SMILES identity (explore / yield without CSMI dedup).
pub fn stable_csmi_key(mol: &CoreMolecule) -> Option<String> {
    canonical_smiles_stable_key(mol)
}

/// Fail-closed key from a SMILES string (parse + aromatize, then stable key).
pub fn stable_csmi_key_of(smiles: &str) -> Option<String> {
    stable_csmi_key(&parse_mol(smiles).ok()?)
}

/// Canonical SMILES of any parseable writing of a structure.
///
/// Tests compare this, not a chematic spelling. `CCO` and `C(C)O` are the
/// same molecule once both sides go through here.
pub fn canon_of(smiles: &str) -> Result<String, ForestError> {
    Ok(canon_smiles(&parse_mol(smiles)?))
}

/// Topological equivalence classes (chematic's `CanonicalRankAtoms(breakTies=False)` analogue).
pub fn ranks(mol: &CoreMolecule) -> Vec<usize> {
    topological_equivalence_classes(mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benzene_aromatic_and_kekule_share_canonical_smiles() {
        let aromatic = parse_mol("c1ccccc1").unwrap();
        let kekule = parse_mol("C1=CC=CC=C1").unwrap();
        assert_eq!(canon_smiles(&aromatic), canon_smiles(&kekule));
        assert_eq!(
            canon_of("c1ccccc1").unwrap(),
            canon_of("C1=CC=CC=C1").unwrap()
        );
        assert!(aromatic.atoms().all(|(_, atom)| atom.aromatic));
        assert!(kekule.atoms().all(|(_, atom)| atom.aromatic));
    }

    #[test]
    fn writings_share_canonical_form() {
        assert_eq!(canon_of("CCO").unwrap(), canon_of("C(C)O").unwrap());
        assert_eq!(
            canon_of("Oc1ccccc1").unwrap(),
            canon_of("c1(O)ccccc1").unwrap()
        );
    }

    #[test]
    fn stable_csmi_key_accepts_simple_molecules() {
        let ethanol = parse_mol("CCO").unwrap();
        let key = stable_csmi_key(&ethanol).expect("ethanol has a stable key");
        assert_eq!(key, canon_smiles(&ethanol));
        assert_eq!(stable_csmi_key_of("C(C)O").as_deref(), Some(key.as_str()));
    }

    #[test]
    fn stable_csmi_key_fails_closed_on_coupled_ez() {
        // Chematic canonical_ez_residual fixture: coupled E/Z → None.
        let smiles = r"CC1CNC(/C=C\C=C/C=C\C2(C)C(=C(O)C(C2)=O)C(/C=C\C=C/C=C\C=C/C=C1)=O)=O";
        let mol = parse_mol(smiles).unwrap();
        assert_eq!(stable_csmi_key(&mol), None);
        // Display CSMI may still exist — just not a safe dedup key.
        assert!(!canon_smiles(&mol).is_empty());
    }

    #[test]
    fn ethane_and_propane_parse() {
        assert_eq!(canon_of("CC").unwrap(), canon_of("C-C").unwrap());
        let propane = parse_mol("CCC").unwrap();
        assert_eq!(propane.atom_count(), 3);
        assert_eq!(canon_of("CCC").unwrap(), canon_of("C(C)C").unwrap());
        let classes = ranks(&propane);
        assert_eq!(classes[0], classes[2]);
        assert_ne!(classes[0], classes[1]);
    }

    #[test]
    fn pyridone_is_aromatic_under_rdkit_parity() {
        let mol = parse_mol("O=c1cccc[nH]1").unwrap();
        let aromatic = mol.atoms().filter(|(_, atom)| atom.aromatic).count();
        assert!(
            aromatic >= 5,
            "RDKit-parity pyridone should keep the ring aromatic, got {aromatic}"
        );
    }
}
