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
        (Element::AT, _) => "At",
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

/// Index of the `]` that closes the `[` at `open` (nested-bracket aware).
fn closing_bracket(s: &str, open: usize) -> Result<usize, ForestError> {
    let bytes = s.as_bytes();
    if open >= bytes.len() || bytes[open] != b'[' {
        return Err(ForestError::Smirks(format!("expected [ at {open} in {s}")));
    }
    let mut depth = 0_i32;
    for i in open..bytes.len() {
        match bytes[i] {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return Ok(i);
                }
            }
            _ => {}
        }
    }
    Err(ForestError::Smirks(format!("unclosed [ in {s}")))
}

/// Map number from a bracket body (`…:12`), if present.
fn map_from_bracket(inner: &str) -> Option<u16> {
    let colon = inner.rfind(':')?;
    let map_text = &inner[colon + 1..];
    if map_text.is_empty() || !map_text.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    map_text.parse().ok()
}

/// Chematic SMILES bond token for a live mol bond between two atoms.
fn bond_smiles_token(mol: &Molecule, a: usize, b: usize) -> Option<&'static str> {
    let (_idx, bond) = mol.bond_between(atom_idx(a), atom_idx(b))?;
    use chematic::core::BondOrder;
    // Prefer aromatic colon when either the bond or either atom is aromatic —
    // chematic apply templates reject `=,:` / `-,:` or-queries.
    if bond.order == BondOrder::Aromatic
        || mol.atom(atom_idx(a)).aromatic
        || mol.atom(atom_idx(b)).aromatic
    {
        return Some(":");
    }
    Some(match bond.order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => "-",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => ":",
        _ => return None,
    })
}

/// True when `s` (starting at `i`) is a SMARTS bond or-query (`-,:`, `=,:`, …).
fn bond_or_query_at(s: &str, i: usize) -> Option<usize> {
    let bytes = s.as_bytes();
    if i >= bytes.len() {
        return None;
    }
    // Bond chars then at least one comma-separated alternative.
    let mut j = i;
    let start = j;
    let is_bond = |b: u8| matches!(b, b'-' | b'=' | b'#' | b':' | b'/' | b'\\');
    if !is_bond(bytes[j]) {
        return None;
    }
    while j < bytes.len() && (is_bond(bytes[j]) || bytes[j] == b',') {
        j += 1;
    }
    let token = &s[start..j];
    if token.contains(',') && token.bytes().any(is_bond) {
        Some(j)
    } else {
        None
    }
}

/// Nearest map number in a specialized atom bracket left of `pos`.
fn map_left_of(s: &str, pos: usize) -> Option<u16> {
    let before = &s[..pos];
    let open = before.rfind('[')?;
    let close = closing_bracket(s, open).ok()?;
    if close >= pos {
        return None;
    }
    map_from_bracket(&s[open + 1..close])
}

/// Nearest map number in a specialized atom bracket right of `pos`.
fn map_right_of(s: &str, pos: usize) -> Option<u16> {
    let bytes = s.as_bytes();
    let mut i = pos;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            let close = closing_bracket(s, i).ok()?;
            return map_from_bracket(&s[i + 1..close]);
        }
        i += 1;
    }
    None
}

/// Resolve SMARTS bond or-queries using the live matched bond order.
fn resolve_bond_or_queries(
    reactant: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<String, ForestError> {
    let bytes = reactant.as_bytes();
    let mut out = String::new();
    let mut i = 0;
    while i < bytes.len() {
        if let Some(end) = bond_or_query_at(reactant, i) {
            let left = map_left_of(reactant, i);
            let right = map_right_of(reactant, end);
            let token = match (left, right) {
                (Some(a), Some(b)) => {
                    let Some(&ai) = mapped.get(&a) else {
                        return Err(ForestError::Smirks(format!(
                            "bond or-query missing map {a}"
                        )));
                    };
                    let Some(&bi) = mapped.get(&b) else {
                        return Err(ForestError::Smirks(format!(
                            "bond or-query missing map {b}"
                        )));
                    };
                    bond_smiles_token(mol, ai, bi).ok_or_else(|| {
                        ForestError::Smirks(format!(
                            "no bond between maps {a}={ai} and {b}={bi}"
                        ))
                    })?
                }
                _ => {
                    // No flanking maps — keep as-is (should be rare).
                    out.push_str(&reactant[i..end]);
                    i = end;
                    continue;
                }
            };
            out.push_str(token);
            i = end;
            continue;
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    Ok(out)
}

/// Product-side `[*…:map]` / charged wildcards → `[El±:map]` from the match.
fn specialize_product_wildcards(
    product: &str,
    mol: &Molecule,
    mapped: &BTreeMap<u16, usize>,
) -> Result<String, ForestError> {
    let bytes = product.as_bytes();
    let mut out = String::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'[' {
            out.push(bytes[i] as char);
            i += 1;
            continue;
        }
        let end = closing_bracket(product, i)?;
        let inner = &product[i + 1..end];
        let rewritten = if let Some(mapno) = map_from_bracket(inner) {
            if let Some(&atom) = mapped.get(&mapno) {
                // Wildcard or query-heavy product atom tied to a reactant map.
                // Plain `[*:map]` stays — chematic accepts it and existing
                // specialize tests rely on that. Charge / & / H0 queries need
                // an element (e.g. `[*&H0&+:1]` → `[NH0+:1]`).
                let plain_star = inner.starts_with('*')
                    && inner.len() > 1
                    && inner.as_bytes()[1] == b':'
                    && inner[2..].bytes().all(|b| b.is_ascii_digit());
                let needs = (inner.starts_with('*') && !plain_star)
                    || (inner.starts_with('#')
                        && inner.contains(|c: char| c == 'v' || c == ';'));
                if needs {                    let a = mol.atom(atom_idx(atom));
                    let sym = element_symbol(a.element, a.aromatic).ok_or_else(|| {
                        ForestError::Smirks(format!("unsupported element for map {mapno}"))
                    })?;
                    let charge = if inner.contains('+') {
                        "+"
                    } else if inner.contains('-') {
                        "-"
                    } else {
                        ""
                    };
                    // Prefer H0 when the query asked for it (n-oxide / S-oxide).
                    let h0 = inner.contains("H0") || inner.contains("h0");
                    let mut tok = String::from("[");
                    tok.push_str(sym);
                    if h0 {
                        tok.push_str("H0");
                    }
                    tok.push_str(charge);
                    tok.push(':');
                    tok.push_str(&mapno.to_string());
                    tok.push(']');
                    Some(tok)
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };
        if let Some(tok) = rewritten {
            out.push_str(&tok);
        } else {
            out.push_str(&product[i..=end]);
        }
        i = end + 1;
    }
    Ok(out)
}

/// Rewrite a Python/RDKit SMIRKS so chematic can parse it for this match.
///
/// - Drop product-side grouping parentheses after `>>`.
/// - Replace every reactant `[…:map]` with `[El:map]` from the matched atom
///   (lowercase when that atom is aromatic — chematic `C` is aliphatic-only).
///   Nested recursive `$()` brackets are handled (not truncated at the first `]`).
/// - Resolve SMARTS bond or-queries (`-,:`, `=,:`) to the live bond order.
/// - Rewrite product `[*&H0&+:map]`-style wildcards to `[ElH0+:map]`.
/// - Leave product organic atoms alone when already chematic-clean.
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
            let end = closing_bracket(reactant, i)?;
            let bracket = &reactant[i + 1..end];
            if let Some(mapno) = map_from_bracket(bracket) {
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
    let reactant = resolve_bond_or_queries(&out, mol, mapped)?;
    let product = specialize_product_wildcards(product, mol, mapped)?;
    Ok(format!("{reactant}>>{product}"))
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
        let end = closing_bracket(product, i)?;
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
    let product_side = smirks.split_once(">>").map(|(_, p)| p).unwrap_or("");
    let wants_h0_charge = product_side.contains("H0") && product_side.contains('+');
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
                let frag = if wants_h0_charge {
                    enforce_charged_h0(frag)
                } else {
                    frag
                };
                let frag = normalize_monatomic_astatine(frag);
                if accept_product(&frag) {
                    pieces.push(frag);
                }
            }
        }
        return Ok(pieces);
    }
    Ok(Vec::new())
}

/// Chematic often saturates `El+` with an implicit H (`[SH+]`) even when the
/// product template asked for H0 (`[*&H0&+:1]` → `[SH0+:1]`). Strip that H in
/// the SMILES round-trip so dialkyl S/N oxides match RDKit's closed-shell form.
fn enforce_charged_h0(mol: Molecule) -> Molecule {
    use crate::mol::{canon_smiles, parse_mol};
    let s = canon_smiles(&mol);
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            // Rewrite [XH+] / [xH+] → [X+] / [x+] (single-letter organic els).
            if i + 5 <= bytes.len()
                && bytes[i + 2] == b'H'
                && bytes[i + 3] == b'+'
                && bytes[i + 4] == b']'
                && bytes[i + 1].is_ascii_alphabetic()
            {
                out.push('[');
                out.push(bytes[i + 1] as char);
                out.push('+');
                out.push(']');
                i += 5;
                continue;
            }
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    if out == s {
        return mol;
    }
    parse_mol(&out).unwrap_or(mol)
}

/// RDKit writes monatomic At as ``[AtH]`` (not organic-subset); chematic yields
/// ``[At]``. Normalize so RDKit CSMI parity matches other halide leaves (Br/Cl
/// organic SMILES already imply the same valence-1 H accounting).
fn normalize_monatomic_astatine(mol: Molecule) -> Molecule {
    use crate::mol::{atom_idx, canon_smiles, parse_mol};
    use chematic::core::Element;
    if mol.atom_count() != 1 {
        return mol;
    }
    let atom = mol.atom(atom_idx(0));
    if atom.element != Element::AT {
        return mol;
    }
    if canon_smiles(&mol) == "[AtH]" {
        return mol;
    }
    parse_mol("[AtH]").unwrap_or(mol)
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

    #[test]
    fn specialize_nested_recursive_and_bond_or_query() {
        let mol = parse_mol("COP(=O)(O)O").unwrap();
        let smirks =
            "[#8;$([#8][#6]):1][#15:2](=[#8:3])([#8:4])[#8:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]";
        let hits = smarts_matches(&mol, smirks.split(">>").next().unwrap()).unwrap();
        assert_eq!(hits.len(), 1);
        let form = specialize_smirks_for_maps(smirks, &mol, &hits[0]).unwrap();
        assert_eq!(form, "[O:1][P:2](=[O:3])([O:4])[O:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]");
        let products = apply_smirks_at(smirks, &mol, &hits[0]).unwrap();
        let got: BTreeSet<_> = products.iter().map(|p| canon_of(&canon_smiles(p)).unwrap()).collect();
        assert_eq!(
            got,
            BTreeSet::from([canon_of("CO").unwrap(), canon_of("O=[PH](O)O").unwrap()])
        );
    }

    #[test]
    fn specialize_hydroxylamine_bond_or_and_n_oxide_charge() {
        let mol = parse_mol("CCNO").unwrap();
        let smirks = "[#7:1]-,:[#8:2]>>([*:1].[*:2])";
        let hits = smarts_matches(&mol, "[#7:1]-,:[#8:2]").unwrap();
        let form = specialize_smirks_for_maps(smirks, &mol, &hits[0]).unwrap();
        assert_eq!(form, "[N:1]-[O:2]>>[*:1].[*:2]");
        let products = apply_smirks_at(smirks, &mol, &hits[0]).unwrap();
        let got: BTreeSet<_> = products.iter().map(|p| canon_of(&canon_smiles(p)).unwrap()).collect();
        assert_eq!(
            got,
            BTreeSet::from([canon_of("CCN").unwrap(), canon_of("O").unwrap()])
        );

        let mol = parse_mol("CN(C)C").unwrap();
        let smirks = "[#7v3H0:1]>>[*&H0&+:1][O-]";
        let hits = smarts_matches(&mol, "[#7v3H0:1]").unwrap();
        let form = specialize_smirks_for_maps(smirks, &mol, &hits[0]).unwrap();
        assert_eq!(form, "[N:1]>>[NH0+:1][O-]");
        let products = apply_smirks_at(smirks, &mol, &hits[0]).unwrap();
        let got: BTreeSet<_> = products.iter().map(|p| canon_of(&canon_smiles(p)).unwrap()).collect();
        assert_eq!(got, BTreeSet::from([canon_of("C[N+](C)(C)[O-]").unwrap()]));
    }

    #[test]
    fn specialize_azo_bond_or_query() {
        let mol = parse_mol("c1ccc(/N=N/c2ccccc2)cc1").unwrap();
        let smirks = "[#7:1]=,:[#7:2]>>[*:1].[*:2]";
        let hits = smarts_matches(&mol, "[#7:1]=,:[#7:2]").unwrap();
        assert!(!hits.is_empty());
        let form = specialize_smirks_for_maps(smirks, &mol, &hits[0]).unwrap();
        assert!(
            form == "[N:1]=[N:2]>>[*:1].[*:2]" || form == "[N:1]:[N:2]>>[*:1].[*:2]",
            "got {form}"
        );
        let products = apply_smirks_at(smirks, &mol, &hits[0]).unwrap();
        assert!(!products.is_empty(), "azo apply empty form={form}");
    }

    #[test]
    fn sulfur_oxidation_ccs_matches_rdkit_forms() {
        let mol = parse_mol("CCS").unwrap();
        let hits = smarts_matches(&mol, "[#16;v2,v4:1]").unwrap();
        let zw = apply_smirks_at("[#16;v2,v4:1]>>[*&H0&+:1][O-]", &mol, &hits[0]).unwrap();
        let oh = apply_smirks_at("[#16;v2,v4:1]>>[*:1]O", &mol, &hits[0]).unwrap();
        assert_eq!(
            zw.iter().map(|p| canon_of(&canon_smiles(p)).unwrap()).collect::<Vec<_>>(),
            vec![canon_of("CC[S+][O-]").unwrap()]
        );
        assert_eq!(
            oh.iter().map(|p| canon_of(&canon_smiles(p)).unwrap()).collect::<Vec<_>>(),
            vec![canon_of("CCSO").unwrap()]
        );
    }
}
