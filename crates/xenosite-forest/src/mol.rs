//! Parse, aromatize (RDKit-parity), and write canonical SMILES.

use chematic::core::{AtomIdx, Molecule as CoreMolecule};
use chematic::perception::apply_aromaticity_rdkit_parity_experimental;
use chematic::smiles::{canonical_smiles, parse, topological_equivalence_classes};

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
}

impl std::fmt::Display for ForestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse(msg) | Self::Smarts(msg) | Self::Smirks(msg) | Self::Kekule(msg) => {
                f.write_str(msg)
            }
        }
    }
}

impl std::error::Error for ForestError {}

/// Parse SMILES and apply RDKit-parity aromaticity.
///
/// Kekulé input is aromatized. Already-aromatic SMILES keep atom indexes.
/// If the parity engine refuses a structure, the parsed mol is kept.
pub fn parse_mol(smiles: &str) -> Result<CoreMolecule, ForestError> {
    let mol = parse(smiles).map_err(|err| ForestError::Parse(err.to_string()))?;
    match apply_aromaticity_rdkit_parity_experimental(&mol) {
        Ok(aromatized) => Ok(aromatized),
        Err(_) => Ok(mol),
    }
}

pub fn canon_smiles(mol: &CoreMolecule) -> String {
    canonical_smiles(mol)
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
        assert!(aromatic.atoms().all(|(_, atom)| atom.aromatic));
        assert!(kekule.atoms().all(|(_, atom)| atom.aromatic));
    }

    #[test]
    fn ethane_and_propane_parse() {
        assert_eq!(canon_smiles(&parse_mol("CC").unwrap()), "CC");
        let propane = parse_mol("CCC").unwrap();
        assert_eq!(propane.atom_count(), 3);
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
