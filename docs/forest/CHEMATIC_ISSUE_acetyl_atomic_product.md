# Chematic bug report draft — acetyl product atomic-number expand

Upstream filing draft (Chematic bug_report template). Vendored chematic
`v1.0.21`. Forest workaround: acetylation product side uses organic-subset
`C(=O)C` instead of `[#6](=[#8])[#6]` (RUST_PYTHON_PARITY work #5).

---

name: Bug report
about: Something isn't working correctly
labels: bug

## Description

SMIRKS apply with an **unmapped product-side** acetyl built from atomic-number
primitives (`[#6](=[#8])[#6]`) yields no usable products. `expand_atomic_number_primitives`
rewrites those atoms to **bracket** aliphatic/aromatic forms (`[C](=[O])[C]`,
`[c](=O)[c]`, …). Applying those bracket forms either returns empty after the
valence gate or leaves residual bracket carbons (e.g. `C(=O)(OCC)[C]`).

The same reaction with organic-subset product atoms (`C(=O)C`) applies
correctly on both aliphatic and aromatic attachment sites.

Matching / reactant-side `#` expansion is fine; the gap is product-side
`#` → bracket expand for newly created atoms that should be ordinary
organic-subset SMILES atoms.

## Reproduction

```rust
use chematic::rxn::{expand_atomic_number_primitives, run_reactants};
use chematic::smiles::parse;

let mol = parse("CCO").unwrap();

// Expand shows bracket C/c × O/o product variants — not organic `C(=O)C`:
let variants = expand_atomic_number_primitives("[O:1]>>[*:1][#6](=[#8])[#6]").unwrap();
// => [O:1]>>[*:1][C](=[O])[C], … [c](=[o])[c]

// Atomic / expanded bracket product: no good acetyl ester
let atomic = run_reactants("[O:1]>>[*:1][#6](=[#8])[#6]", &[&mol]);
// empty or products with leftover bracket [C] / [c]

// Organic-subset product works:
let organic = run_reactants("[O:1]>>[*:1]C(=O)C", &[&mol]).unwrap();
// => CC(=O)OCC (canonical)
```

Forest mirror (after specialize of reactant `#8h1` → `O`):

```text
apply "[#8h1:1]>>[*:1][#6](=[#8])[#6]" on CCO  → []
apply "[#8h1:1]>>[*:1]C(=O)C"           on CCO  → CC(=O)OCC
```

Same organic product also works for phenol (`Oc1ccccc1`) and pyrrole
(`[nH]1cccc1`) after reactant specialize to `O` / `n`.

## Expected behavior

Either:

1. `expand_atomic_number_primitives` should emit organic-subset product
   atoms (`C(=O)C` / `O`) for unmapped newly created atoms when bracket
   notation is unnecessary; or
2. `apply` / `run_reactants` on expanded bracket product atoms
   (`[C](=[O])[C]`) should produce the same molecule as organic
   `C(=O)C` (correct implicit H, no residual bracket atom in the product
   SMILES).

Aliphatic and aromatic expand branches that are chemically valid for the
product should both apply; invalid aromatic carbonyl product variants may
still be skipped, but the aliphatic bracket form must not silently fail
relative to `C(=O)C`.

## Actual behavior

- `[#6](=[#8])[#6]` / `[#6](=O)[#6]` product sides expand only to bracket
  `[C]`/`[c]`/`[O]`/`[o]` spellings.
- Those forms do not produce a clean acetyl adduct under apply (empty set
  or leftover `[C]`/`[c]` in the product).
- Hand-written organic `C(=O)C` applies and yields `CC(=O)OCC` from ethanol.

## Forest workaround

`organic_product_variants` rewrites product `#` to organic-subset aliphatic
then aromatic (`C`/`c`, …) and tries those forms in `apply_smirks_at` —
aliphatic first so chemically correct acetyl wins; aromatic forms included
when aliphatic does not apply. Catalog SMARTS keep `[*:1][#6](=[#8])[#6]`.

## Environment

- chematic version: v1.0.21 (vendored sparse checkout in xenosite-forest)
- Rust edition / caller: `chematic::rxn::{expand_atomic_number_primitives, apply_reaction_match, run_reactants}` via xenosite-forest `apply_smirks_at`
- OS: Linux
