# Vendored chematic patch

`chematic-v1.0.21-atom-tag-visit-order.patch` applies on top of the
`vendor/chematic` submodule (tag `v1.0.21`) after a sparse checkout.

## What the patch adds

API (the real delta):

- `Atom.tag: Option<u32>` — non-chemical; copied by clone / apply / fragments;
  ignored by SMILES write and canon; not emitted as `:n`
- `Molecule::set_tag`
- `write_with_order` / `canonical_smiles_with_order`

Sparse-vendor glue (small):

- Neutralize nested `Cargo.toml` workspace (crates join this repo’s workspace)
- Facade features limited to `smiles` / `perception` / `smarts` / `rxn`
- Drop path-only dev-deps that point at non-sparse crates

## Apply

```bash
./scripts/vendor-chematic.sh
```

Do not commit patched files inside the submodule. `.gitmodules` sets
`ignore = dirty`.
