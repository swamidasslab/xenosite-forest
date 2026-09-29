# Forest ↔ XMET mapping (`xf:`)

Living SSSOM: [`mappings/xmet-forest.sssom.tsv`](../../mappings/xmet-forest.sssom.tsv).
Do **not** commit a `.gz` — `crates/xenosite-forest/build.rs` gzips into `OUT_DIR`
and `include_bytes!` embeds it. Runtime API decompresses once (`OnceLock`).

## CURIEs

| Prefix | IRI |
| --- | --- |
| `xf:` | `https://w3id.org/xenosite/forest/` |
| `xmet:` | `https://xenosite.org/ontology/xmet#` |

Object IDs are path segments:

- Catalog / leaf: `xf:Tautomerization`, `xf:PhaseOne`
- Pattern: `xf:Tautomerization/tautomer_h` → [`BoundPattern`](../../crates/xenosite-forest/src/bound_pattern.rs)

Locked mapping policy: prefer long `xf:` names that `resolve()`. Rainbow
short codes (`CJ`/`SO`/…) do not ship. Phase-color catalogs that exist in
Rust (`StableOxygenation`, `UnstableOxygenation`, `Reduction`) may appear;
`Conjugation` is deferred until that catalog is correct.

## Resolve

Rust / Python / WASM:

```python
from xenosite.forest import resolve, forest_xmet_sssom, BoundPattern

resolve("xf:Tautomerization")              # RuleSet
resolve("xf:Tautomerization/tautomer_h")   # BoundPattern
text = forest_xmet_sssom()                 # decompressed TSV
```

`resolve` walks public path segments via `LEAF_CTORS` + `ROOT_CATALOGS`.
Always returns product-door objects (never native/legacy).

`BoundPattern` is RuleSet-compatible for metabolize: pass-through to the owning
leaf with a pattern-name filter. Index / name access: `rule["tautomer_h"]`,
`rule[0]` (catalog → child set; leaf → bound pattern).

## Coverage ownership

| Concern | Owner |
| --- | --- |
| Every shipped `xf:` leaf/pattern resolves; inventory vs SSSOM; embed round-trip | **this repo** (`tests/sssom_coverage.rs`, `tests/forest/rust/test_sssom_resolve.py`) |
| Chemist home nesting / SPARQL cover / `always_with` vs ontology | **xenosite-xmet** |
| Namer / SMARTS emit | **xenosite-tagger** |

Chemist definition Forest-mention lint (when sibling `xenosite-xmet` is
checked out): `tests/forest/rust/test_xmet_definition_lint.py`. Semantic
patches land in `xmet.yaml` upstream — not by rewriting locked SSSOM rows.

## Release snapshot → xmet

On forest `vX.Y.Z`:

1. Copy `mappings/xmet-forest.sssom.tsv` into xenosite-xmet as
   `data/mappings/forest/xmet-forest.sssom.vX.Y.Z.tsv`.
2. Refresh living `data/mappings/xmet-forest.sssom.tsv` there to match.
3. Re-run xmet `make forest-pattern-coverage` (or equivalent) on the snapshot.
