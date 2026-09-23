//! SMARTS matches as forest map-number dictionaries (`{mapno: atom_idx}`).

use std::collections::BTreeMap;

use chematic::smarts::{find_matches, parse_smarts};

use crate::mol::{ForestError, Molecule, atom_usize};

/// Forest `xf.smarts_matches`: each hit is map number → target atom index.
/// Hits without map 1 are dropped, matching the Python door.
pub fn smarts_matches(
    mol: &Molecule,
    smarts: &str,
) -> Result<Vec<BTreeMap<u16, usize>>, ForestError> {
    let query = parse_smarts(smarts).map_err(|err| ForestError::Smarts(err.to_string()))?;
    let mut hits = Vec::new();
    for embedding in find_matches(&query, mol) {
        let mut mapped = BTreeMap::new();
        for (query_idx, target) in embedding {
            if let Some(mapno) = query.atoms.get(query_idx).and_then(|atom| atom.atom_map) {
                if mapno != 0 {
                    mapped.insert(mapno, atom_usize(target));
                }
            }
        }
        if mapped.contains_key(&1) {
            hits.push(mapped);
        }
    }
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;

    #[test]
    fn hydroxylation_h_count_partition() {
        let propane = parse_mol("CCC").unwrap();
        let h1 = smarts_matches(&propane, "[#6h1:1]").unwrap();
        let h2 = smarts_matches(&propane, "[#6h2,#6h3:1]").unwrap();
        assert!(h1.is_empty(), "propane has no implicit-H=1 carbon");
        let sites: Vec<usize> = h2.into_iter().map(|m| m[&1]).collect();
        assert_eq!(sites.len(), 3);

        let benzene = parse_mol("c1ccccc1").unwrap();
        let aryl = smarts_matches(&benzene, "[#6h1:1]").unwrap();
        assert_eq!(aryl.len(), 6);
        let aliphatic = smarts_matches(&benzene, "[#6h2,#6h3:1]").unwrap();
        assert!(aliphatic.is_empty());
    }

    #[test]
    fn aromatic_or_double_bond_matches_benzene() {
        let benzene = parse_mol("c1ccccc1").unwrap();
        let hits = smarts_matches(&benzene, "[#6:1]=,:[#6:2]").unwrap();
        assert_eq!(hits.len(), 6, "six unique ring bonds (VF2 uniquify)");
    }

    #[test]
    fn total_h_vs_implicit_h() {
        let mol = parse_mol("CCO").unwrap();
        let total_h3 = smarts_matches(&mol, "[#6H3:1]").unwrap();
        let implicit_h3 = smarts_matches(&mol, "[#6h3:1]").unwrap();
        assert_eq!(total_h3.len(), 1);
        assert_eq!(implicit_h3.len(), 1);
        assert_eq!(total_h3[0][&1], implicit_h3[0][&1]);
    }
}
