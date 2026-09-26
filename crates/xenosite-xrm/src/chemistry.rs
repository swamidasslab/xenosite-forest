//! Parse reactant/product; SMARTS match; formula delta; site extraction.
//! No forest types.

use crate::error::{Error, Result};
use crate::term::SiteRef;
use chematic::core::{AtomIdx, Molecule};
use chematic::smarts::{find_matches, parse_smarts};
use chematic::smiles::parse;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AtomMap {
    pub pairs: Vec<(usize, usize)>,
}

impl AtomMap {
    pub fn from_map_numbers(reactant: &Molecule, product: &Molecule) -> Self {
        let mut r_by_map: BTreeMap<u16, usize> = BTreeMap::new();
        for (idx, atom) in reactant.atoms() {
            if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                r_by_map.insert(m, idx.0 as usize);
            }
        }
        let mut pairs = Vec::new();
        for (idx, atom) in product.atoms() {
            if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                if let Some(&r) = r_by_map.get(&m) {
                    pairs.push((r, idx.0 as usize));
                }
            }
        }
        Self { pairs }
    }
}

/// One caller tag, optionally localized with `@map` / `@map1,map2`.
///
/// Examples: `chem:hydroxylation`, `chem:hydroxylation@1`, `forest.rule:NDealkylation@2,3`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalizedTag {
    pub base: String,
    pub site: SiteRef,
}

impl LocalizedTag {
    pub fn parse(raw: &str) -> Self {
        if let Some((base, rest)) = raw.split_once('@') {
            let mut map_nums = Vec::new();
            let mut reactant_atoms = Vec::new();
            for part in rest.split(',') {
                let p = part.trim();
                if p.is_empty() {
                    continue;
                }
                if let Some(a) = p.strip_prefix('a').and_then(|s| s.parse::<usize>().ok()) {
                    reactant_atoms.push(a);
                } else if let Ok(n) = p.parse::<u16>() {
                    if n != 0 {
                        map_nums.push(n);
                    }
                }
            }
            map_nums.sort_unstable();
            map_nums.dedup();
            reactant_atoms.sort_unstable();
            reactant_atoms.dedup();
            Self {
                base: base.to_string(),
                site: SiteRef {
                    map_nums,
                    reactant_atoms,
                },
            }
        } else {
            Self {
                base: raw.to_string(),
                site: SiteRef::default(),
            }
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MappedReaction {
    pub reactant_smiles: String,
    pub product_smiles: String,
    #[serde(default)]
    pub atom_map: AtomMap,
    /// Opaque CURIEs only (e.g. `forest.rule:Hydroxylation` or `chem:hydroxylation@1`).
    #[serde(default)]
    pub tags: Vec<String>,
}

impl MappedReaction {
    pub fn new(reactant: impl Into<String>, product: impl Into<String>) -> Self {
        Self {
            reactant_smiles: reactant.into(),
            product_smiles: product.into(),
            atom_map: AtomMap::default(),
            tags: Vec::new(),
        }
    }

    pub fn with_tags(mut self, tags: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.tags = tags.into_iter().map(Into::into).collect();
        self
    }
}

pub struct ReactionChemistry {
    pub reactant: Molecule,
    pub product: Molecule,
    pub atom_map: AtomMap,
    pub tags: BTreeSet<String>,
    pub localized_tags: Vec<LocalizedTag>,
    pub delta: BTreeMap<String, i32>,
}

impl ReactionChemistry {
    pub fn prepare(query: &MappedReaction) -> Result<Self> {
        let reactant =
            parse(&query.reactant_smiles).map_err(|e| Error::Chemistry(format!("reactant: {e}")))?;
        let product =
            parse(&query.product_smiles).map_err(|e| Error::Chemistry(format!("product: {e}")))?;
        let atom_map = if query.atom_map.pairs.is_empty() {
            AtomMap::from_map_numbers(&reactant, &product)
        } else {
            query.atom_map.clone()
        };
        let mut delta = BTreeMap::new();
        let rf = heavy_formula(&reactant);
        let pf = heavy_formula(&product);
        for k in rf.keys().chain(pf.keys()) {
            let d = pf.get(k).copied().unwrap_or(0) - rf.get(k).copied().unwrap_or(0);
            if d != 0 {
                delta.insert(k.clone(), d);
            }
        }
        let localized_tags: Vec<LocalizedTag> =
            query.tags.iter().map(|t| LocalizedTag::parse(t)).collect();
        Ok(Self {
            reactant,
            product,
            atom_map,
            tags: query.tags.iter().cloned().collect(),
            localized_tags,
            delta,
        })
    }

    /// True if SMARTS matches; when `require_map` non-empty, those map nums must hit.
    pub fn smarts_hits(
        &self,
        on_product: bool,
        smarts: &str,
        require_map: &[u16],
    ) -> Result<bool> {
        Ok(!self.smarts_sites(on_product, smarts, require_map)?.is_empty())
    }

    /// Site-localized SMARTS hits on reactant or product.
    pub fn smarts_sites(
        &self,
        on_product: bool,
        smarts: &str,
        require_map: &[u16],
    ) -> Result<Vec<SiteRef>> {
        let mol = if on_product {
            &self.product
        } else {
            &self.reactant
        };
        let query = parse_smarts(smarts).map_err(|e| Error::Chemistry(e.to_string()))?;
        let mut sites = Vec::new();
        for embedding in find_matches(&query, mol) {
            let mut query_maps = BTreeSet::new();
            let mut map_nums = BTreeSet::new();
            let mut reactant_atoms = BTreeSet::new();
            for (qi, target) in &embedding {
                if let Some(m) = query.atoms.get(*qi).and_then(|a| a.atom_map) {
                    if m != 0 {
                        query_maps.insert(m);
                    }
                }
                let atom = mol.atom(AtomIdx(target.0));
                if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                    map_nums.insert(m);
                }
                if !on_product {
                    reactant_atoms.insert(target.0 as usize);
                } else if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                    for (idx, a) in self.reactant.atoms() {
                        if a.atom_map == Some(m) {
                            reactant_atoms.insert(idx.0 as usize);
                        }
                    }
                }
            }
            if !require_map.is_empty() && !require_map.iter().all(|m| query_maps.contains(m)) {
                continue;
            }
            if !require_map.is_empty() {
                let mut focused_maps = BTreeSet::new();
                let mut focused_atoms = BTreeSet::new();
                for (qi, target) in &embedding {
                    let Some(qm) = query.atoms.get(*qi).and_then(|a| a.atom_map) else {
                        continue;
                    };
                    if !require_map.contains(&qm) {
                        continue;
                    }
                    let atom = mol.atom(AtomIdx(target.0));
                    if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                        focused_maps.insert(m);
                    }
                    if !on_product {
                        focused_atoms.insert(target.0 as usize);
                    } else if let Some(m) = atom.atom_map.filter(|&m| m != 0) {
                        for (idx, a) in self.reactant.atoms() {
                            if a.atom_map == Some(m) {
                                focused_atoms.insert(idx.0 as usize);
                            }
                        }
                    }
                }
                if !focused_maps.is_empty() || !focused_atoms.is_empty() {
                    map_nums = focused_maps;
                    reactant_atoms = focused_atoms;
                }
            }
            let site = SiteRef {
                map_nums: map_nums.into_iter().collect(),
                reactant_atoms: reactant_atoms.into_iter().collect(),
            };
            if !sites.iter().any(|s| s == &site) {
                sites.push(site);
            }
        }
        Ok(sites)
    }

    pub fn any_mapped_reactant_aromatic(&self, smarts: &str) -> Result<Option<bool>> {
        let query = parse_smarts(smarts).map_err(|e| Error::Chemistry(e.to_string()))?;
        for embedding in find_matches(&query, &self.reactant) {
            for (_, target) in embedding {
                return Ok(Some(self.reactant.atom(AtomIdx(target.0)).aromatic));
            }
        }
        Ok(None)
    }

    /// Match opaque tags by base name; return localized sites from `@…` suffixes.
    pub fn tag_sites(&self, wanted_bases: &[String]) -> Vec<SiteRef> {
        let mut out = Vec::new();
        for tag in &self.localized_tags {
            if wanted_bases.iter().any(|b| b == &tag.base) {
                if !out.iter().any(|s| s == &tag.site) {
                    out.push(tag.site.clone());
                }
            }
        }
        out
    }
}

fn heavy_formula(mol: &Molecule) -> BTreeMap<String, i32> {
    let mut counts = BTreeMap::new();
    for (_, atom) in mol.atoms() {
        if atom.element.atomic_number() <= 1 {
            continue;
        }
        *counts.entry(atom.element.symbol().to_string()).or_insert(0) += 1;
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_localized_tag() {
        let t = LocalizedTag::parse("chem:hydroxylation@1,2");
        assert_eq!(t.base, "chem:hydroxylation");
        assert_eq!(t.site.map_nums, vec![1, 2]);
    }
}
