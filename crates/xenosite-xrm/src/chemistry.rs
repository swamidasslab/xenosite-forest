//! Parse reactant/product; SMARTS match; formula delta. No forest types.

use crate::error::{Error, Result};
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MappedReaction {
    pub reactant_smiles: String,
    pub product_smiles: String,
    #[serde(default)]
    pub atom_map: AtomMap,
    /// Opaque CURIEs only (e.g. `forest.rule:Hydroxylation`).
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
        Ok(Self {
            reactant,
            product,
            atom_map,
            tags: query.tags.iter().cloned().collect(),
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
        let mol = if on_product {
            &self.product
        } else {
            &self.reactant
        };
        let query = parse_smarts(smarts).map_err(|e| Error::Chemistry(e.to_string()))?;
        for embedding in find_matches(&query, mol) {
            if require_map.is_empty() {
                return Ok(true);
            }
            let mut mapped = BTreeSet::new();
            for (qi, _) in &embedding {
                if let Some(m) = query.atoms.get(*qi).and_then(|a| a.atom_map) {
                    if m != 0 {
                        mapped.insert(m);
                    }
                }
            }
            if require_map.iter().all(|m| mapped.contains(m)) {
                return Ok(true);
            }
        }
        Ok(false)
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
