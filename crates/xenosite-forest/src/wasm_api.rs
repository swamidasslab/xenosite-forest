//! Browser WASM exports. Enabled with `--features wasm`.

use wasm_bindgen::prelude::wasm_bindgen;

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
