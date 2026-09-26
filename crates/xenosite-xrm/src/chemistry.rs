//! Thin chematic wrappers for mapped reactant→product queries.
//!
//! No forest types. Atom maps are caller-supplied or read from SMILES map numbers.

use crate::error::{Error, Result};
use chematic::core::{AtomIdx, Molecule};
use chematic::smarts::{find_matches, parse_smarts};
use chematic::smiles::parse;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Reactant atom index → product atom index.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AtomMap {
    pub pairs: Vec<(usize, usize)>,
}

impl AtomMap {
    pub fn from_pairs(pairs: impl IntoIterator<Item = (usize, usize)>) -> Self {
        Self {
            pairs: pairs.into_iter().collect(),
        }
    }

    /// Build from shared atom-map numbers on reactant and product molecules.
    pub fn from_map_numbers(reactant: &Molecule, product: &Molecule) -> Self {
        let mut r_by_map: BTreeMap<u16, usize> = BTreeMap::new();
        for (idx, atom) in reactant.atoms() {
            if let Some(m) = atom.atom_map {
                if m != 0 {
                    r_by_map.insert(m, idx.0 as usize);
                }
            }
        }
        let mut pairs = Vec::new();
        for (idx, atom) in product.atoms() {
            if let Some(m) = atom.atom_map {
                if m != 0 {
                    if let Some(&r) = r_by_map.get(&m) {
                        pairs.push((r, idx.0 as usize));
                    }
                }
            }
        }
        Self { pairs }
    }

    pub fn reactant_to_product(&self) -> BTreeMap<usize, usize> {
        self.pairs.iter().copied().collect()
    }
}

/// Input to the namer: a mapped transform plus optional opaque tags.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MappedReaction {
    pub reactant_smiles: String,
    pub product_smiles: String,
    /// Explicit atom correspondence. Empty → try SMILES atom-map numbers.
    #[serde(default)]
    pub atom_map: AtomMap,
    /// Opaque CURIEs (e.g. `forest.rule:Hydroxylation`). Never resolved in-process.
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

    pub fn with_atom_map(mut self, map: AtomMap) -> Self {
        self.atom_map = map;
        self
    }
}

/// Parsed chemistry used during assignment matching.
pub struct ReactionChemistry {
    pub reactant: Molecule,
    pub product: Molecule,
    pub atom_map: AtomMap,
    pub tags: BTreeSet<String>,
    pub reactant_formula: BTreeMap<String, i32>,
    pub product_formula: BTreeMap<String, i32>,
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
        let reactant_formula = heavy_formula(&reactant);
        let product_formula = heavy_formula(&product);
        let mut delta = BTreeMap::new();
        let keys: BTreeSet<_> = reactant_formula
            .keys()
            .chain(product_formula.keys())
            .cloned()
            .collect();
        for k in keys {
            let d = product_formula.get(&k).copied().unwrap_or(0)
                - reactant_formula.get(&k).copied().unwrap_or(0);
            if d != 0 {
                delta.insert(k, d);
            }
        }
        Ok(Self {
            reactant,
            product,
            atom_map,
            tags: query.tags.iter().cloned().collect(),
            reactant_formula,
            product_formula,
            delta,
        })
    }

    /// True if SMARTS matches; optional map numbers must be present on the hit.
    /// When an atom_map exists, mapped atoms should participate in it (soft when empty map).
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
        let r2p = self.atom_map.reactant_to_product();
        let p2r: BTreeMap<usize, usize> = r2p.iter().map(|(&r, &p)| (p, r)).collect();
        for embedding in find_matches(&query, mol) {
            let mut mapped: BTreeMap<u16, usize> = BTreeMap::new();
            for (query_idx, target) in embedding {
                if let Some(mapno) = query.atoms.get(query_idx).and_then(|a| a.atom_map) {
                    if mapno != 0 {
                        mapped.insert(mapno, target.0 as usize);
                    }
                }
            }
            if !require_map.is_empty() && !require_map.iter().all(|m| mapped.contains_key(m)) {
                continue;
            }
            if self.atom_map.pairs.is_empty() || require_map.is_empty() {
                return Ok(true);
            }
            let ok = require_map.iter().all(|m| {
                let idx = mapped[m];
                if on_product {
                    p2r.contains_key(&idx)
                } else {
                    r2p.contains_key(&idx)
                }
            });
            if ok {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub fn any_mapped_reactant_aromatic(&self, smarts: &str) -> Result<Option<bool>> {
        let query = parse_smarts(smarts).map_err(|e| Error::Chemistry(e.to_string()))?;
        for embedding in find_matches(&query, &self.reactant) {
            for (_q, target) in embedding {
                let atom = self.reactant.atom(AtomIdx(target.0));
                return Ok(Some(atom.aromatic));
            }
        }
        Ok(None)
    }
}

fn heavy_formula(mol: &Molecule) -> BTreeMap<String, i32> {
    let mut counts = BTreeMap::new();
    for (_, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        if z <= 1 {
            continue;
        }
        let sym = atom.element.symbol().to_string();
        *counts.entry(sym).or_insert(0) += 1;
    }
    let charge: i32 = mol.atoms().map(|(_, a)| i32::from(a.charge)).sum();
    if charge != 0 {
        counts.insert("charge".into(), charge);
    }
    counts
}
