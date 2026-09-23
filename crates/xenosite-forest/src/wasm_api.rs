//! Browser WASM exports. Enabled with `--features wasm`.
//!
//! `#[wasm_bindgen]` wraps [`crate::xf::ForestMol`] as a JS class the same
//! way PyO3 wraps it as a Python class: one JS object holds one Rust payload.

use wasm_bindgen::prelude::wasm_bindgen;

use crate::pattern::PatternInfo;
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};
use crate::xf::ForestMol as Held;
use crate::{canon_smiles, hydroxylate, parse_mol};

#[wasm_bindgen]
pub fn forest_canon_smiles(smiles: &str) -> Result<String, String> {
    let mol = parse_mol(smiles).map_err(|err| err.to_string())?;
    Ok(canon_smiles(&mol))
}

#[wasm_bindgen]
pub fn forest_hydroxylate(smiles: &str) -> Result<String, String> {
    let mol = parse_mol(smiles).map_err(|err| err.to_string())?;
    let products = hydroxylate(&mol).map_err(|err| err.to_string())?;
    Ok(products.join("\n"))
}

/// JS class wrapping [`ForestMol`]. Cache lives on the Rust payload.
#[wasm_bindgen(js_name = ForestMol)]
pub struct JsForestMol {
    inner: Held,
}

#[wasm_bindgen(js_class = ForestMol)]
impl JsForestMol {
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str) -> Result<JsForestMol, String> {
        Ok(Self {
            inner: Held::parse(smiles).map_err(|err| err.to_string())?,
        })
    }

    #[wasm_bindgen(getter)]
    pub fn has_forest(&self) -> bool {
        self.inner.xf().has_forest()
    }

    #[wasm_bindgen(getter)]
    pub fn csmi(&self) -> String {
        self.inner.xf().csmi().to_string()
    }

    #[wasm_bindgen]
    pub fn clear_structure(&self) {
        self.inner.xf().clear_structure();
    }

    #[wasm_bindgen]
    pub fn copy(&self) -> JsForestMol {
        Self {
            inner: self.inner.copy_mol(),
        }
    }

    #[wasm_bindgen(js_name = rwCopy)]
    pub fn rw_copy(&self) -> JsForestMol {
        Self {
            inner: self.inner.rw_copy(),
        }
    }

    #[wasm_bindgen(js_name = wipeForest)]
    pub fn wipe_forest(&self) {
        self.inner.wipe_forest();
    }
}

/// JS wrap of [`RuleSet`]. Patterns live on the Rust payload after `add_*`.
#[wasm_bindgen(js_name = RuleSet)]
pub struct JsRuleSet {
    inner: RuleSet,
}

#[wasm_bindgen(js_class = RuleSet)]
impl JsRuleSet {
    #[wasm_bindgen(constructor)]
    pub fn new() -> JsRuleSet {
        Self {
            inner: RuleSet::new(None, []),
        }
    }

    #[wasm_bindgen(js_name = hydroxylation)]
    pub fn hydroxylation() -> JsRuleSet {
        Self {
            inner: crate::hydroxylation::hydroxylation(),
        }
    }

    #[wasm_bindgen(js_name = addHydroxyl)]
    pub fn add_hydroxyl(&mut self, name: &str, smarts: &str) {
        self.inner.push(PatternInfo::hydroxyl(name, smarts));
    }

    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.inner.patterns().len()
    }

    pub fn metabolize(&self, smiles: &str) -> Result<String, String> {
        let mol = parse_mol(smiles).map_err(|err| err.to_string())?;
        let emissions = self
            .inner
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .map_err(|err| err.to_string())?;
        Ok(emissions
            .into_iter()
            .flat_map(|emission| emission.products)
            .collect::<Vec<_>>()
            .join("\n"))
    }
}

impl Default for JsRuleSet {
    fn default() -> Self {
        Self::new()
    }
}
