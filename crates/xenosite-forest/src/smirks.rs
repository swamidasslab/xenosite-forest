//! Apply one SMIRKS match (chematic's replacement for isotope-pinned RunReactants).
//!
//! Python SMIRKS often use RDKit primitives (`[#6H3:1]`, `[#7,#8:2]`, product
//! grouping parens) that chematic's reaction parser rejects. Matching still
//! uses the full SMARTS door; apply rewrites the template to organic atoms
//! from the live match maps.
//!
//! Product-side `[#N]` is expanded to **organic-subset** aliphatic and
//! aromatic spellings (`C`/`c`, …), not chematic's bracket `[C]`/`[c]` expand
//! (those miss or leave residual brackets — see CHEMATIC_ISSUE_acetyl…).
//! Aliphatic forms are tried first; aromatic forms are included when the
//! aliphatic spelling does not apply. `accept_product` keeps only chemically
//! valid results.

use std::collections::BTreeMap;

use chematic::core::Element;
use chematic::rxn::{apply_reaction_match, find_reaction_matches};

use crate::mol::{ForestError, Molecule, atom_idx, atom_usize};
use crate::valence::accept_product;

fn element_symbol(el: Element, aromatic: bool) -> Option<&'static str> {
    Some(match (el, aromatic) {
        (Element::C, true) => "c",
        (Element::N, true) => "n",
        (Element::O, true) => "o",
        (Element::S, true) => "s",
        (Element::P, true) => "p",
        (Element::C, false) => "C",
        (Element::N, false) => "N",
        (Element::O, false) => "O",
        (Element::S, false) => "S",
        (Element::P, false) => "P",
        (Element::F, _) => "F",
        (Element::CL, _) => "Cl",
        (Element::BR, _) => "Br",
        (Element::I, _) => "I",
        (Element::H, _) => "H",
        _ => return None,
    })
}

/// Elements chematic/SMARTS treat as aromaticable in `#N` expand.
fn aromatic_symbol(atomic_number: u8) -> Option<&'static str> {
    match atomic_number {
        6 => Some("c"),
        7 => Some("n"),
        8 => Some("o"),
        15 => Some("p"),
        16 => Some("s"),
        _ => None,
    }
}

fn aliphatic_symbol(atomic_number: u8) -> Option<&'static str> {
    Element::from_atomic_number(atomic_number).and_then(|el| element_symbol(el, false))
}

/// Rewrite a Python/RDKit SMIRKS so chematic can parse it for this match.
///
/// - Drop product-side grouping parentheses after `>>`.
/// - Replace every reactant `[…:map]` with `[El:map]` from the matched atom
///   (lowercase when that atom is aromatic — chematic `C` is aliphatic-only).
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
                    let a = mol.atom(atom_idx(atom));
                    let sym = element_symbol(a.element, a.aromatic).ok_or_else(|| {
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

/// One `[#N…]` site on the product side and its organic aliphatic/aromatic spellings.
struct HashSite {
    start: usize,
    end: usize, // exclusive, includes `]`
    aliphatic: String,
    aromatic: Option<String>,
}

/// Parse product-side `[#N]` / `[#N:map]` / `[#N;H1]` / `[#N;H1:map]` into organic spellings.
fn product_hash_sites(product: &str) -> Result<Vec<HashSite>, ForestError> {
    let mut sites = Vec::new();
    let bytes = product.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'[' {
            i += 1;
            continue;
        }
        let end = product[i..]
            .find(']')
            .map(|o| i + o)
            .ok_or_else(|| ForestError::Smirks(format!("unclosed [ in product {product}")))?;
        let inner = &product[i + 1..end];
        if !inner.starts_with('#') {
            i = end + 1;
            continue;
        }
        let number_end = inner[1..]
            .bytes()
            .position(|b| !b.is_ascii_digit())
            .map(|p| p + 1)
            .unwrap_or(inner.len());
        let number_text = &inner[1..number_end];
        let suffix = &inner[number_end..];
        let atomic_number: u8 = number_text
            .parse()
            .map_err(|_| ForestError::Smirks(format!("bad # atomic number in [{inner}]")))?;
        let Some(ali) = aliphatic_symbol(atomic_number) else {
            return Err(ForestError::Smirks(format!(
                "unsupported product # atomic number in [{inner}]"
            )));
        };
        // Supported suffixes: empty, :map, ;H1, ;H1:map (same boundary as chematic expand).
        let (h_count, map_suffix) = if suffix.is_empty() {
            (None, "")
        } else if suffix.starts_with(':')
            && suffix.len() > 1
            && suffix[1..].bytes().all(|b| b.is_ascii_digit())
        {
            (None, suffix)
        } else if let Some(map) = suffix.strip_prefix(";H1:")
            && !map.is_empty()
            && map.bytes().all(|b| b.is_ascii_digit())
        {
            (Some(1_u8), &suffix[3..])
        } else if suffix == ";H1" {
            (Some(1_u8), "")
        } else {
            // Complex product primitive — leave for chematic / other try forms.
            i = end + 1;
            continue;
        };
        let aliphatic = organic_atom_token(ali, h_count, map_suffix);
        let aromatic = aromatic_symbol(atomic_number)
            .map(|aro| organic_atom_token(aro, h_count, map_suffix));
        sites.push(HashSite {
            start: i,
            end: end + 1,
            aliphatic,
            aromatic,
        });
        i = end + 1;
    }
    Ok(sites)
}

fn organic_atom_token(sym: &str, h_count: Option<u8>, map_suffix: &str) -> String {
    // Unmapped bare organic (`C`, `O`) when no H/map constraints; otherwise bracket.
    if h_count.is_none() && map_suffix.is_empty() {
        return sym.to_string();
    }
    let mut out = String::from("[");
    out.push_str(sym);
    if let Some(h) = h_count {
        out.push('H');
        out.push_str(&h.to_string());
    }
    out.push_str(map_suffix);
    out.push(']');
    out
}

/// Product spellings: all-aliphatic organic first, then aromatic substitutions.
///
/// Chematic's `#` expand uses brackets (`[C]`/`[c]`) that miss on apply. Organic
/// `C`/`c` matches equivalent `#` coverage; aromatic variants are tried when
/// aliphatic does not apply. Callers still run [`accept_product`].
pub fn organic_product_variants(smirks: &str) -> Result<Vec<String>, ForestError> {
    let (reactant, product) = smirks
        .split_once(">>")
        .ok_or_else(|| ForestError::Smirks(format!("SMIRKS missing >>: {smirks}")))?;
    let product = product.trim();
    let product = if product.starts_with('(') && product.ends_with(')') {
        &product[1..product.len() - 1]
    } else {
        product
    };
    let sites = product_hash_sites(product)?;
    if sites.is_empty() {
        return Ok(vec![format!("{reactant}>>{product}")]);
    }

    // Variant 0: every site aliphatic.
    let mut variants = Vec::new();
    variants.push(rewrite_product(reactant, product, &sites, &vec![false; sites.len()]));

    // Further variants: flip aromaticable sites to aromatic (single flips, then
    // all-aromatic). Enough to cover chematic's aliphatic/aromatic branches
    // without a full 2^n blow-up for large products.
    let aromaticable: Vec<usize> = sites
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.aromatic.is_some().then_some(i))
        .collect();
    for &idx in &aromaticable {
        let mut flags = vec![false; sites.len()];
        flags[idx] = true;
        variants.push(rewrite_product(reactant, product, &sites, &flags));
    }
    if aromaticable.len() > 1 {
        let mut flags = vec![false; sites.len()];
        for &idx in &aromaticable {
            flags[idx] = true;
        }
        variants.push(rewrite_product(reactant, product, &sites, &flags));
    }

    // Stable unique order (aliphatic-all first).
    let mut seen = std::collections::BTreeSet::new();
    variants.retain(|v| seen.insert(v.clone()));
    Ok(variants)
}

fn rewrite_product(reactant: &str, product: &str, sites: &[HashSite], aromatic: &[bool]) -> String {
    let mut out = String::new();
    out.push_str(reactant);
    out.push_str(">>");
    let mut cursor = 0;
    for (i, site) in sites.iter().enumerate() {
        out.push_str(&product[cursor..site.start]);
        let token = if aromatic[i] {
            site.aromatic.as_ref().unwrap_or(&site.aliphatic)
        } else {
            &site.aliphatic
        };
        out.push_str(token);
        cursor = site.end;
    }
    out.push_str(&product[cursor..]);
    out
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
    let mut try_forms = Vec::new();
    // Organic product expand of `#` (aliphatic first, then aromatic) on the
    // specialized template, then on the original if different.
    for base in [&specialized, smirks] {
        match organic_product_variants(base) {
            Ok(vars) => {
                for v in vars {
                    if !try_forms.contains(&v) {
                        try_forms.push(v);
                    }
                }
            }
            Err(_) => {
                let owned = base.to_string();
                if !try_forms.contains(&owned) {
                    try_forms.push(owned);
                }
            }
        }
    }
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

    #[test]
    fn specialize_uses_aromatic_carbon_for_aryl_oxdehal() {
        let mol = parse_mol("Clc1ccccc1").unwrap();
        let hits = smarts_matches(&mol, "[#9,#17,#35,#53,#85:1]-[#6:2]").unwrap();
        assert_eq!(hits.len(), 1);
        let python = "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O";
        let form = specialize_smirks_for_maps(python, &mol, &hits[0]).unwrap();
        assert_eq!(form, "[Cl:1]-[c:2]>>[*:1].[*:2]O", "form={form}");
        let products = apply_smirks_at(python, &mol, &hits[0]).unwrap();
        let phenol = canon_of("Oc1ccccc1").unwrap();
        assert!(
            products
                .iter()
                .any(|p| canon_of(&canon_smiles(p)).unwrap() == phenol),
            "got {:?} form={form}",
            products.iter().map(canon_smiles).collect::<Vec<_>>()
        );
    }

    #[test]
    fn organic_product_variants_aliphatic_first_then_aromatic() {
        let vars = organic_product_variants("[O:1]>>[*:1][#6](=[#8])[#6]").unwrap();
        assert_eq!(vars[0], "[O:1]>>[*:1]C(=O)C", "{vars:?}");
        assert!(
            vars.iter().any(|v| v.contains('c') || v.contains('o')),
            "aromatic organic forms included: {vars:?}"
        );
        // Aliphatic-all is first so chemically correct acetyl wins before
        // aromatic carbonyl spellings.
        assert!(!vars[0].chars().any(|c| matches!(c, 'c' | 'n' | 'o' | 's' | 'p')));
    }

    /// Atomic product `#` applies via organic expand (aliphatic branch).
    #[test]
    fn acetylation_atomic_product_applies_via_organic_expand() {
        let mol = parse_mol("CCO").unwrap();
        let hits = smarts_matches(&mol, "[#8h1:1]").unwrap();
        assert_eq!(hits.len(), 1);
        let products =
            apply_smirks_at("[#8h1:1]>>[*:1][#6](=[#8])[#6]", &mol, &hits[0]).unwrap();
        assert_eq!(
            products
                .iter()
                .map(|p| canon_of(&canon_smiles(p)).unwrap())
                .collect::<Vec<_>>(),
            vec![canon_of("CC(=O)OCC").unwrap()]
        );
    }

    /// specialize rewrites reactant `#` to aliphatic `O`/`N`/`S` or aromatic
    /// `n`. Product `#` expands to organic aliphatic (and aromatic fallback).
    #[test]
    fn acetylation_covers_aliphatic_and_aromatic_heteroatom_branches() {
        let smirks = "[#7h1,#7h2,#8h1,#16h1:1]>>[*:1][#6](=[#8])[#6]";
        let cases = [
            ("CCO", "[O:1]>>[*:1][#6](=[#8])[#6]", "CC(=O)OCC"),
            ("CCN", "[N:1]>>[*:1][#6](=[#8])[#6]", "CCNC(C)=O"),
            ("CS", "[S:1]>>[*:1][#6](=[#8])[#6]", "CSC(C)=O"),
            ("Oc1ccccc1", "[O:1]>>[*:1][#6](=[#8])[#6]", "CC(=O)Oc1ccccc1"),
            ("Nc1ccccc1", "[N:1]>>[*:1][#6](=[#8])[#6]", "CC(=O)Nc1ccccc1"),
            ("Sc1ccccc1", "[S:1]>>[*:1][#6](=[#8])[#6]", "CC(=O)Sc1ccccc1"),
            ("[nH]1cccc1", "[n:1]>>[*:1][#6](=[#8])[#6]", "CC(=O)n1cccc1"),
        ];
        for (smiles, want_form, want_prod) in cases {
            let mol = parse_mol(smiles).unwrap();
            let hits = smarts_matches(&mol, "[#7h1,#7h2,#8h1,#16h1:1]").unwrap();
            assert_eq!(hits.len(), 1, "{smiles}");
            let form = specialize_smirks_for_maps(smirks, &mol, &hits[0]).unwrap();
            assert_eq!(form, want_form, "{smiles}");
            let products = apply_smirks_at(smirks, &mol, &hits[0]).unwrap();
            let got: BTreeSet<String> = products
                .iter()
                .map(|p| canon_of(&canon_smiles(p)).unwrap())
                .collect();
            let want = BTreeSet::from([canon_of(want_prod).unwrap()]);
            assert_eq!(got, want, "{smiles} form={form}");
        }
    }
}
