# Metabolic Forest notebooks

Requires the Rust extension (`uv sync` with maturin build backend, or
`maturin develop -m crates/xenosite-forest/Cargo.toml --features python,extension-module`).

| Notebook | Topic |
|----------|--------|
| [forest_product_door.ipynb](forest_product_door.ipynb) | Handle-first `xenosite.forest` API: `find_path`, BFS graph, shared `MetabolicNetwork`, metabolize |

Optional: `pip install xenosite-forest[rdkit]` for `ForestMol.to_rdkit()` drawing cells.

Regression smoke: `pytest tests/forest/rust/test_tutorial_snippets.py`.
