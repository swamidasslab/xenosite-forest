//! Browser WASM exports. Enabled with `--features wasm`.
//!
//! `#[wasm_bindgen]` wraps [`crate::xf::ForestMol`] as a JS class the same
//! way PyO3 wraps it as a Python class: one JS object holds one Rust payload.

use wasm_bindgen::prelude::wasm_bindgen;

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
