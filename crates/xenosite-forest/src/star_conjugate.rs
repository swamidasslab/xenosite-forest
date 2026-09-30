//! Collapse conjugation adducts to chematic `*` and CX-label the new star.
//!
//! Two separate jobs:
//! - **Collapse** — multi-atom adducts from a conjugation edit become one `*`
//!   (drop born atoms past `stamp_end`; [`Molecule::set_wildcard`] on the
//!   attachment). Not used to invent a star just so a label can be stored.
//! - **Label** — free-text names (`GlcA`, `GSH`, …) live as **tag-keyed**
//!   entries on [`ForestMol`] via [`ForestMol::set_cx_label`]. Chematic
//!   `Molecule` has no string CX labels (only numeric `r_groups`, cleared by
//!   `set_wildcard`). Emit projects those tags into
//!   [`CxSmiles::atom_labels`] via [`ForestMol::write_cxsmiles`].

use std::collections::{HashMap, HashSet};
use std::sync::atomic::Ordering;

use chematic::core::AtomIdx;

use crate::forest_mol::ForestMol;
use crate::mol::{atom_idx, atom_usize};

/// Hardcoded CX atomLabel for each conjugation / reactivity leaf.
pub fn conjugate_star_label(leaf: &str) -> Option<&'static str> {
    Some(match leaf {
        "Acetylation" => "Ac",
        "Sulfation" => "SO3",
        "Glucuronidation" => "GlcA",
        "Glutathionation" | "GlutathionationNoThiol" | "GSH" => "GSH",
        "Protein" => "Protein",
        "DNA" => "DNA",
        "Cyanide" => "Cyanide",
        _ => return None,
    })
}

/// Collapse newly minted conjugate atoms to `*` and keep the product terminal.
///
/// Records `label` on each new star's [`Tag`] (prior stars keep their labels).
pub fn collapse_conjugate_to_star(product: &ForestMol, label: &str) -> ForestMol {
    let stamp = product.stamp_end();
    let mut mol = product.mol().clone();
    let n = mol.atom_count();

    let new_atoms: HashSet<usize> = (0..n)
        .filter(|&i| !mol.atom(atom_idx(i)).wildcard)
        .filter(|&i| match product.tag_of(i) {
            Some(tag) => tag.get() >= stamp,
            None => true,
        })
        .collect();

    if new_atoms.is_empty() {
        return label_terminal(product);
    }

    let mut star_of_attach: HashMap<usize, usize> = HashMap::new();
    for &ni in &new_atoms {
        for (nbr, _) in mol.neighbors(atom_idx(ni)) {
            let parent = atom_usize(nbr);
            if !new_atoms.contains(&parent) {
                star_of_attach.entry(parent).or_insert(ni);
            }
        }
    }

    if star_of_attach.is_empty() {
        return label_terminal(product);
    }

    let keep_stars: HashSet<usize> = star_of_attach.values().copied().collect();
    // Tags of stars before remove/remap (same Tag values survive from_apply).
    let new_star_tags: Vec<_> = keep_stars
        .iter()
        .filter_map(|&i| product.tag_of(i))
        .collect();

    for &star in &keep_stars {
        mol.set_wildcard(atom_idx(star));
    }

    let mut remove: Vec<usize> = new_atoms
        .into_iter()
        .filter(|i| !keep_stars.contains(i))
        .collect();
    remove.sort_unstable_by(|a, b| b.cmp(a));

    let mut src_to_new: Vec<Option<usize>> = (0..n).map(Some).collect();
    for i in remove {
        let remap = mol.remove_atom(AtomIdx(i as u32));
        src_to_new = src_to_new
            .into_iter()
            .map(|cur| {
                cur.and_then(|j| remap.get(j).copied().flatten().map(|a| a.0 as usize))
            })
            .collect();
    }

    let mut out = product.from_apply(mol, &src_to_new);
    for tag in new_star_tags {
        out.set_cx_label(tag, label);
    }
    out.is_terminal_product.store(true, Ordering::Relaxed);
    out
}

fn label_terminal(product: &ForestMol) -> ForestMol {
    let out = product.copy_mol();
    out.is_terminal_product.store(true, Ordering::Relaxed);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::parse_mol;
    use crate::rules::{
        acetylation, cyanide, dna, glucuronidation, glutathionation, gsh, protein, sulfation,
    };
    use crate::ruleset::RuleSet;
    use chematic::core::{Atom, BondOrder};
    use chematic::smiles::parse_cxsmiles;

    fn leaf_products(rule: RuleSet, smiles: &str) -> Vec<ForestMol> {
        use crate::ruleset::{accept_all_rules, accept_all_sites};
        let mol = ForestMol::new(parse_mol(smiles).unwrap());
        rule.metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .map(|e| e.expect("metabolize"))
            .flat_map(|e| e.products)
            .collect()
    }

    fn assert_star_labeled(products: &[ForestMol], label: &str) {
        assert!(!products.is_empty(), "expected products for label {label}");
        for p in products {
            assert!(
                p.mol().atoms().any(|(_, a)| a.wildcard),
                "expected * in product {}",
                p.csmi()
            );
            let cx = p.write_cxsmiles();
            assert!(cx.contains('*'), "CXSMILES missing *: {cx}");
            assert!(
                cx.contains(label),
                "CXSMILES missing label {label:?}: {cx}"
            );
            let stars: Vec<_> = p
                .mol()
                .atoms()
                .filter(|(_, a)| a.wildcard)
                .filter_map(|(i, _)| p.tag_of(atom_usize(i)))
                .collect();
            assert!(
                stars.iter().any(|t| p.cx_label(*t) == Some(label)),
                "tag-keyed label missing for {label} on {}",
                p.csmi()
            );
        }
    }

    #[test]
    fn conjugation_products_are_stars_with_hardcoded_labels() {
        let cases: &[(&str, RuleSet, &str)] = &[
            ("Ac", acetylation(), "CCO"),
            ("SO3", sulfation(), "CCO"),
            ("GlcA", glucuronidation(), "c1ccccc1O"),
            ("GSH", glutathionation(), "C1OC1c1ccccc1"),
            ("GSH", gsh(), "C1OC1c1ccccc1"),
            ("Protein", protein(), "C1OC1c1ccccc1"),
            ("DNA", dna(), "C1OC1c1ccccc1"),
            ("Cyanide", cyanide(), "C1OC1c1ccccc1"),
        ];
        for &(label, ref rule, smi) in cases {
            let products = leaf_products(rule.clone(), smi);
            assert_star_labeled(&products, label);
        }
    }

    #[test]
    fn product_csmis_emit_cx_labels() {
        use crate::ruleset::{accept_all_rules, accept_all_sites};
        let mol = ForestMol::new(parse_mol("c1ccccc1O").unwrap());
        let emissions: Vec<_> = glucuronidation()
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .map(|e| e.expect("metabolize"))
            .collect();
        assert!(!emissions.is_empty());
        for e in emissions {
            for cx in e.product_csmis() {
                assert!(cx.contains('*'), "{cx}");
                assert!(cx.contains("GlcA"), "{cx}");
                assert!(cx.contains('|'), "expected CX block: {cx}");
            }
        }
    }

    /// Labeling is tag-keyed: pick the star tag, set the string. No collapse.
    #[test]
    fn set_cx_label_assigns_to_chosen_star_tag() {
        let first = leaf_products(sulfation(), "Oc1ccc(O)cc1")
            .into_iter()
            .next()
            .expect("sulfation product");
        let first_star_tag = first
            .mol()
            .atoms()
            .find(|(_, a)| a.wildcard)
            .and_then(|(i, _)| first.tag_of(atom_usize(i)))
            .expect("SO3 star tag");
        assert_eq!(first.cx_label(first_star_tag), Some("SO3"));

        let free_o = first
            .mol()
            .atoms()
            .find_map(|(idx, a)| {
                if a.element.symbol() != "O" || a.wildcard {
                    return None;
                }
                let touches_star = first
                    .mol()
                    .neighbors(idx)
                    .any(|(n, _)| first.mol().atom(n).wildcard);
                if touches_star {
                    None
                } else {
                    Some(atom_usize(idx))
                }
            })
            .expect("free OH");

        // Attach a real `*` and label that tag — not a C that we then collapse.
        let mut mol = first.mol().clone();
        let new_idx = mol.add_atom(Atom::wildcard());
        mol.add_bond(atom_idx(free_o), new_idx, BondOrder::Single)
            .expect("attach star");
        let mut p = ForestMol::product(mol, &first);
        let new_star_tag = p
            .tag_of(atom_usize(new_idx))
            .expect("born star tagged");
        p.set_cx_label(new_star_tag, "GlcA");

        assert_eq!(p.cx_label(first_star_tag), Some("SO3"));
        assert_eq!(p.cx_label(new_star_tag), Some("GlcA"));

        let cx = p.write_cxsmiles();
        let parsed = parse_cxsmiles(&cx).expect(&cx);
        let old_i = p
            .mol()
            .atoms()
            .find(|(i, a)| a.wildcard && p.tag_of(atom_usize(*i)) == Some(first_star_tag))
            .map(|(i, _)| atom_usize(i))
            .expect("old star idx");
        let new_i = p
            .mol()
            .atoms()
            .find(|(i, a)| a.wildcard && p.tag_of(atom_usize(*i)) == Some(new_star_tag))
            .map(|(i, _)| atom_usize(i))
            .expect("new star idx");
        assert_eq!(parsed.atom_labels[old_i].as_deref(), Some("SO3"), "{cx}");
        assert_eq!(parsed.atom_labels[new_i].as_deref(), Some("GlcA"), "{cx}");
    }
}
