//! Tautomer (and related) normalization via chematic pick + Forest adopt.
//!
//! Chematic chooses the form (`normalize_zwitterion` → `remove_hydrogens` →
//! `canonical_tautomer`). Forest births the result with existing tracing
//! ([`ForestMol::product`] when atom indexes are stable, else
//! [`ForestMol::from_edit_product`]).
//!
//! Chematic-chem rebuilds drop caller-tag sidecars (no `copy_atom_tags_from`).
//! Until upstream carries tags, snapshot survivor tags before the pick and
//! restamp in order afterward (heavy + isotopic H keep relative order across
//! the whole stream). Probe: [`crate::chematic_features`] (ignored until fixed).
//! Prefer [`ForestMol::normalize_tautomer`] at call sites / in tests.

use chematic::chem::{canonical_tautomer, normalize_zwitterion, remove_hydrogens};
use chematic::core::Element;

use crate::chematic_tags::{get_label, set_label};
use crate::forest_mol::{ForestMol, IntoForestMol, as_forest_mol};
use crate::labels::Tag;
use crate::mol::{ForestError, Molecule, atom_idx, canon_smiles};

/// Outcome of [`normalize_tautomer`].
#[derive(Clone, Debug)]
pub struct NormalizedTautomer {
    /// Tagged molecule in the chematic-preferred form.
    pub mol: ForestMol,
    /// True when chemistry changed relative to the input.
    pub changed: bool,
}

/// Removable explicit H: element H with no isotope (matches chematic-chem).
fn is_removable_explicit_h(mol: &Molecule, i: usize) -> bool {
    let a = mol.atom(atom_idx(i));
    a.element == Element::H && a.isotope.is_none()
}

/// Tags of atoms that survive [`remove_hydrogens`], in survivor order.
fn survivor_tags(mol: &Molecule) -> Vec<Option<Tag>> {
    (0..mol.atom_count())
        .filter(|&i| !is_removable_explicit_h(mol, i))
        .map(|i| get_label(mol, atom_idx(i)))
        .collect()
}

/// Write `tags` onto `mol` by index (same length expected).
fn restamp_tags(mol: &mut Molecule, tags: &[Option<Tag>]) {
    debug_assert_eq!(
        mol.atom_count(),
        tags.len(),
        "pick must keep survivor atom order"
    );
    let n = mol.atom_count().min(tags.len());
    for (i, tag) in tags.iter().copied().take(n).enumerate() {
        set_label(mol, atom_idx(i), tag);
    }
}

/// Chematic preferred form: zwitterion → strip explicit H → canonical tautomer.
///
/// Tags: snapshot survivors → run chematic stream → restamp (temporary until
/// chematic-chem preserves sidecars on rebuild).
pub fn chematic_tautomer_pick(mol: &Molecule) -> Molecule {
    let tags = survivor_tags(mol);
    let mol = normalize_zwitterion(mol);
    let mol = remove_hydrogens(&mol);
    let mut mol = canonical_tautomer(&mol);
    restamp_tags(&mut mol, &tags);
    mol
}

/// Normalize tautomer / charge form; keep Forest atom tracing.
///
/// Intake is polymorphic ([`as_forest_mol`]). The chematic pick is adopted with
/// the same doors as other edit products — not `ForestMol::new` (that restamps).
pub fn normalize_tautomer<M: IntoForestMol>(mol: M) -> Result<NormalizedTautomer, ForestError> {
    let parent = as_forest_mol(mol)?;
    let picked = chematic_tautomer_pick(parent.mol());
    if canon_smiles(parent.mol()) == canon_smiles(&picked) {
        return Ok(NormalizedTautomer {
            mol: parent,
            changed: false,
        });
    }
    // Identity-preserving rebuilds keep atom order → index-stable product.
    // Explicit-H removal can drop atoms → tag-follow via from_edit_product.
    let mol = if picked.atom_count() == parent.mol().atom_count() {
        ForestMol::product(picked, &parent)
    } else {
        parent.from_edit_product(picked)
    };
    Ok(NormalizedTautomer { mol, changed: true })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::as_forest_mol;
    use crate::chematic_tags::get_label;
    use crate::mol::atom_idx;

    fn parent_tags(mol: &ForestMol) -> Vec<Option<u16>> {
        (0..mol.mol().atom_count())
            .map(|i| mol.tag_of(i).map(|t| t.get()))
            .collect()
    }

    fn heavy_tags(mol: &ForestMol) -> Vec<u16> {
        (0..mol.mol().atom_count())
            .filter(|&i| mol.mol().atom(atom_idx(i)).element.atomic_number() > 1)
            .filter_map(|i| mol.tag_of(i).map(|t| t.get()))
            .collect()
    }

    #[test]
    fn pass_through_when_already_canonical() {
        let out = normalize_tautomer("CCO").unwrap();
        assert!(!out.changed);
        assert_eq!(
            out.mol.csmi().as_ref(),
            as_forest_mol("CCO").unwrap().csmi().as_ref()
        );
    }

    #[test]
    fn keto_enol_agree() {
        let enol = normalize_tautomer("OC=C").unwrap();
        let keto = normalize_tautomer("CC=O").unwrap();
        assert_eq!(enol.mol.csmi().as_ref(), keto.mol.csmi().as_ref());
    }

    #[test]
    fn acetylacetone_enol_and_keto_agree() {
        let enol = normalize_tautomer("CC(O)=CC(=O)C").unwrap();
        let keto = normalize_tautomer("CC(=O)CC(=O)C").unwrap();
        assert_eq!(enol.mol.csmi().as_ref(), keto.mol.csmi().as_ref());
    }

    #[test]
    fn pyrazole_n1h_n2h_agree() {
        // Chematic tautomer corpus: c1cc[nH]n1
        let a = normalize_tautomer("c1cc[nH]n1").unwrap();
        let forms = chematic::chem::enumerate_tautomers(a.mol.mol());
        assert!(forms.len() >= 2, "pyrazole should enumerate ≥2 forms");
        let b = normalize_tautomer(forms[1].clone()).unwrap();
        assert_eq!(a.mol.csmi().as_ref(), b.mol.csmi().as_ref());
    }

    #[test]
    fn tetrazole_1h_2h_agree() {
        let a = normalize_tautomer("c1n[nH]nn1").unwrap();
        let b = normalize_tautomer("c1nnn[nH]1").unwrap();
        assert_eq!(a.mol.csmi().as_ref(), b.mol.csmi().as_ref());
    }

    #[test]
    fn pyridone_hydroxypyridine_agree() {
        let a = normalize_tautomer("O=c1cccc[nH]1").unwrap();
        let b = normalize_tautomer("Oc1ccccn1").unwrap();
        assert_eq!(a.mol.csmi().as_ref(), b.mol.csmi().as_ref());
    }

    #[test]
    fn bis_enol_multi_site_agree() {
        let a = normalize_tautomer("OC=CCC=C(O)C").unwrap();
        let b = normalize_tautomer("CC(O)=CCC=CO").unwrap();
        assert_eq!(a.mol.csmi().as_ref(), b.mol.csmi().as_ref());
    }

    #[test]
    fn amide_unchanged() {
        let out = normalize_tautomer("CC(=O)N").unwrap();
        assert!(!out.changed);
    }

    #[test]
    fn idempotent() {
        let once = normalize_tautomer("OC=C").unwrap();
        let twice = normalize_tautomer(&once.mol).unwrap();
        assert!(!twice.changed);
        assert_eq!(once.mol.csmi().as_ref(), twice.mol.csmi().as_ref());
    }

    #[test]
    fn index_stable_pick_keeps_parent_tags() {
        let parent = as_forest_mol("OC=C").unwrap();
        let before = parent_tags(&parent);
        let out = normalize_tautomer(&parent).unwrap();
        assert!(out.changed, "enol should move toward keto");
        assert_eq!(out.mol.mol().atom_count(), parent.mol().atom_count());
        assert!(out.mol.shares_tag_gen(&parent));
        let after = parent_tags(&out.mol);
        assert_eq!(
            before, after,
            "index-stable adopt must keep Forest tags: before={before:?} after={after:?}"
        );
        for (i, want) in after.iter().enumerate() {
            assert_eq!(
                get_label(out.mol.mol(), atom_idx(i)).map(|t| t.get()),
                *want
            );
        }
    }

    #[test]
    fn zwitterion_amino_acid_normalizes() {
        // Chematic standardize corpus: [NH3+]CC(=O)[O-]
        let out = normalize_tautomer("[NH3+]CC(=O)[O-]").unwrap();
        let smi = out.mol.csmi();
        assert!(
            !smi.as_ref().contains("[NH3+]") && !smi.as_ref().contains("[O-]"),
            "expected neutralized form, got {smi}"
        );
    }

    #[test]
    fn nitro_charge_sep_left_alone_or_stable() {
        // Permanent charge-sep; chematic leave alone or keep nitro form.
        let out = normalize_tautomer("C[N+](=O)[O-]").unwrap();
        let again = normalize_tautomer(&out.mol).unwrap();
        assert!(!again.changed);
        assert_eq!(out.mol.csmi().as_ref(), again.mol.csmi().as_ref());
    }

    #[test]
    fn polymorphic_molecule_intake() {
        let mol = as_forest_mol("CC=O").unwrap();
        let out = normalize_tautomer(mol.mol().clone()).unwrap();
        assert!(!out.changed);
    }

    #[test]
    fn explicit_h_removal_keeps_heavy_tags() {
        let parent = as_forest_mol("[H]OC([H])=C([H])[H]").unwrap();
        let before = heavy_tags(&parent);
        assert!(
            parent.mol().atom_count() > before.len(),
            "fixture must carry explicit H"
        );
        let out = normalize_tautomer(&parent).unwrap();
        assert!(out.changed);
        assert!(out.mol.shares_tag_gen(&parent));
        let mut after = heavy_tags(&out.mol);
        let mut expect = before.clone();
        after.sort_unstable();
        expect.sort_unstable();
        assert_eq!(
            after, expect,
            "heavy tag set must survive begin/end restamp + adopt (order may change)"
        );
    }

    #[test]
    fn pick_restamps_survivor_tags_end_to_end() {
        let parent = as_forest_mol("[H]OC([H])=C([H])[H]").unwrap();
        let want = survivor_tags(parent.mol());
        let picked = chematic_tautomer_pick(parent.mol());
        assert_eq!(picked.atom_count(), want.len());
        for (i, tag) in want.iter().enumerate() {
            assert_eq!(get_label(&picked, atom_idx(i)), *tag, "index {i}");
        }
    }

    #[test]
    fn chematic_pick_stages_compose() {
        let mol = as_forest_mol("[NH3+]CC(=O)[O-]").unwrap();
        let picked = chematic_tautomer_pick(mol.mol());
        let via_api = normalize_tautomer(&mol).unwrap();
        assert_eq!(canon_smiles(&picked).as_str(), via_api.mol.csmi().as_ref());
    }
}
