//! Apply one SMIRKS match (chematic's replacement for isotope-pinned RunReactants).
//!
//! Python SMIRKS often use RDKit primitives (`[#6H3:1]`, `[#7,#8:2]`, product
//! grouping parens) that chematic's reaction parser rejects. Matching still
//! uses the full SMARTS door; apply rewrites the template to organic atoms
//! from the live match maps.

use std::collections::BTreeMap;

use chematic::core::Element;
use chematic::rxn::{apply_reaction_match, find_reaction_matches};

use crate::mol::{ForestError, Molecule, atom_idx, atom_usize};
use crate::valence::accept_product;

fn element_symbol(el: Element) -> Option<&'static str> {
    Some(match el {
        Element::C => "C",
        Element::N => "N",
        Element::O => "O",
        Element::S => "S",
        Element::P => "P",
        Element::F => "F",
        Element::CL => "Cl",
        Element::BR => "Br",
        Element::I => "I",
        Element::H => "H",
        _ => return None,
    })
}

/// Rewrite a Python/RDKit SMIRKS so chematic can parse it for this match.
///
/// - Drop product-side grouping parentheses after `>>`.
/// - Replace every reactant `[…:map]` with `[El:map]` from the matched atom.
/// - Leave product `[* :map]` / organic atoms alone when already chematic-clean.
pub fn specialize_smirks_for_maps(
    smirks: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<String, ForestError> {
    let (reactant, product) = smirks
        .split_once(">>")
        .ok_or_else(|| ForestError::Smirks(format!("SMIRKS missing >>: {smirks}")))?;
    let product = product.trim();
    let product = if product.starts_with('(') && product.ends_with(')') {
        &product[1..product.len() - 1]
    } else {
        product
    };

    let mut out = String::new();
    let bytes = reactant.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            let end = reactant[i..]
                .find(']')
                .map(|o| i + o)
                .ok_or_else(|| ForestError::Smirks(format!("unclosed [ in {smirks}")))?;
            let bracket = &reactant[i + 1..end];
            if let Some(colon) = bracket.rfind(':') {
                let mapno: u16 = bracket[colon + 1..]
                    .parse()
                    .map_err(|_| ForestError::Smirks(format!("bad map in [{bracket}]")))?;
                if let Some(&atom) = mapped.get(&mapno) {
                    let el = mol.atom(atom_idx(atom)).element;
                    let sym = element_symbol(el).ok_or_else(|| {
                        ForestError::Smirks(format!("unsupported element for map {mapno}"))
                    })?;
                    out.push('[');
                    out.push_str(sym);
                    out.push(':');
                    out.push_str(&mapno.to_string());
                    out.push(']');
                    i = end + 1;
                    continue;
                }
            }
            out.push_str(&reactant[i..=end]);
            i = end + 1;
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out.push_str(">>");
    out.push_str(product);
    Ok(out)
}

/// Run `smirks` on the match whose atom maps equal `mapped`.
///
/// Maps on `mapped` that the SMIRKS does not use are ignored. Aromatic SMARTS
/// may name a ring; the apply template names only the reacting atoms.
///
/// Products are split into connected fragments and dropped when the forest
/// valence gate refuses them.
pub fn apply_smirks_at(
    smirks: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<Vec<Molecule>, ForestError> {
    let specialized = match specialize_smirks_for_maps(smirks, mol, mapped) {
        Ok(s) => s,
        Err(_) => smirks.to_string(),
    };
    let try_forms = if specialized == smirks {
        vec![smirks.to_string()]
    } else {
        vec![specialized, smirks.to_string()]
    };
    for form in try_forms {
        match apply_smirks_raw(&form, mol, mapped) {
            Ok(pieces) if !pieces.is_empty() => return Ok(pieces),
            Ok(_) => continue,
            // Chematic rejects many RDKit SMIRKS primitives; skip that form.
            Err(_) => continue,
        }
    }
    Ok(Vec::new())
}

fn apply_smirks_raw(
    smirks: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<Vec<Molecule>, ForestError> {
    let matches = find_reaction_matches(smirks, &[mol])
        .map_err(|err| ForestError::Smirks(err.to_string()))?;
    for reaction_match in matches {
        let positions = reaction_match
            .atom_map_positions(smirks)
            .map_err(|err| ForestError::Smirks(err.to_string()))?;
        let same = mapped.iter().all(|(&mapno, &atom)| {
            positions
                .get(&mapno)
                .is_none_or(|(_slot, idx)| atom_usize(*idx) == atom)
        });
        if !same {
            continue;
        }
        let products = apply_reaction_match(smirks, &[mol], &reaction_match, true)
            .map_err(|err| ForestError::Smirks(err.to_string()))?;
        let Some(products) = products else {
            return Ok(Vec::new());
        };
        let mut pieces = Vec::new();
        for product in products {
            for frag in product.fragments() {
                if accept_product(&frag) {
                    pieces.push(frag);
                }
            }
        }
        return Ok(pieces);
    }
    Ok(Vec::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, canon_smiles, parse_mol};
    use crate::smarts::smarts_matches;
    use std::collections::BTreeSet;

    #[test]
    fn hydroxylation_smirks_on_named_atom() {
        let mol = parse_mol("CC").unwrap();
        let hits = smarts_matches(&mol, "[#6h3:1]").unwrap();
        assert_eq!(hits.len(), 2);
        let products = apply_smirks_at("[C:1]>>[C:1]O", &mol, &hits[0]).unwrap();
        assert_eq!(products.len(), 1);
        assert_eq!(
            canon_of(&canon_smiles(&products[0])).unwrap(),
            canon_of("CCO").unwrap()
        );
    }

    #[test]
    fn cleavage_smirks_yields_two_fragments() {
        let mol = parse_mol("CN").unwrap();
        let hits = smarts_matches(&mol, "[#6H3:1][#7:2]").unwrap();
        assert_eq!(hits.len(), 1);
        let products = apply_smirks_at("[C:1][N:2]>>[N:2].[C:1](=O)O", &mol, &hits[0]).unwrap();
        let got: BTreeSet<String> = products
            .iter()
            .map(|p| canon_of(&canon_smiles(p)).unwrap())
            .collect();
        let want = BTreeSet::from([canon_of("N").unwrap(), canon_of("O=CO").unwrap()]);
        assert_eq!(got, want);
    }

    #[test]
    fn specialize_rewrites_rdkit_dealkylation_for_anisole() {
        let mol = parse_mol("COc1ccccc1").unwrap();
        let hits = smarts_matches(&mol, "[#6H3:1][#7,#8H0,#16:2]").unwrap();
        assert_eq!(hits.len(), 1);
        let python = "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)";
        let form = specialize_smirks_for_maps(python, &mol, &hits[0]).unwrap();
        assert_eq!(form, "[C:1][O:2]>>[*:2].[*:1](=O)O");
        let products = apply_smirks_at(python, &mol, &hits[0]).unwrap();
        let phenol = canon_of("Oc1ccccc1").unwrap();
        assert!(
            products
                .iter()
                .any(|p| canon_of(&canon_smiles(p)).unwrap() == phenol),
            "form={form}"
        );
    }
}
