# Vendored chematic patch

`chematic-v1.0.21-atom-tag-visit-order.patch` is applied on top of the
`vendor/chematic` submodule (pinned to tag `v1.0.21`) after a sparse
checkout of the crates we use.

## What the patch adds

- `Atom.tag: Option<u32>` — non-chemical id copied by clone / apply /
  fragments; ignored by SMILES write and canon; not emitted as `:n`
- `Molecule::set_tag`
- `write_with_order` / `canonical_smiles_with_order` — SMILES string plus
  DFS visit order (`order[k]` = mol index of the k-th atom in the string)
- Workspace/Cargo.toml edits so the sparse tree builds as a path
  dependency inside this repo (inline package metadata; facade features
  limited to `smiles` / `perception` / `smarts` / `rxn`)

## Apply

```bash
./scripts/vendor-chematic.sh
```

Do not commit patched files inside the submodule. `.gitmodules` sets
`ignore = dirty` so a patched working tree does not dirty the parent.
