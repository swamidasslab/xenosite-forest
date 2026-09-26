# xenosite-xrm

Config-driven **xenobiotic reaction naming** for Metabolic Forest.

## Boundary (hard)

This crate **must not** depend on `xenosite-forest` or any forest rule implementation.

- Naming reads only **SKOS** (JSON-LD), **SSSOM** (TSV), and **JSONL** assignment files.
- Opaque CURIEs such as `forest.rule:Hydroxylation` or `forest.pattern:Dealkylation/hemiaminal` may appear in SSSOM and in caller-supplied tags. They are **strings**, never resolved by importing forest code.
- A future optional adapter may call Rust `find_path` and pass those strings as tags (see [ISSUE_find_path_naming.md](ISSUE_find_path_naming.md)).

## Formats

| Role | Standard | Path |
| --- | --- | --- |
| Thesaurus | [SKOS](https://www.w3.org/TR/skos-reference/) JSON-LD | `data/ontology/xrm.skos.jsonld` |
| Crosswalks | [SSSOM](https://mapping-commons.github.io/sssom/) TSV | `data/mappings/*.sssom.tsv` |
| Assignment | JSON Lines | `data/assignments/xenobiotic.jsonl` |
| Bundle | JSON manifest | `data/manifest.json` |

XRM is a **separate** ontology (not an extension of MeSH or KEGG). MeSH Phase I/II and MOP process terms are linked via SSSOM `broadMatch` / `exactMatch`. Enzyme types are orthogonal and never primary labels.

## Usage

```rust
use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

let namer = Namer::from_manifest(DEFAULT_MANIFEST)?;
let terms = namer.name_smiles("CC", "CCO", &[])?;
for t in terms {
    println!("{} ({}) path={:?}", t.pref_label, t.id, t.path_labels);
}
```

With opaque forest tags (from an external caller or find_path adapter):

```rust
let terms = namer.name_smiles(
    "CCO",
    "CCO",
    &["forest.rule:Glucuronidation"],
)?;
```

Each [`Term`](src/term.rs) includes ontology identity, pref/alt labels, broader/narrower links, intra- and inter-ontology matches, specificity (depth / leaf), path labels, and assignment evidence ids.

## Design decisions

1. **Separate XRM + declared synonyms** — MeSH is too coarse (Phase I/II only); KEGG RCLASS is RDM-centric; MOP has good leaves but no xenobiotic Phase spine.
2. **SKOS + SSSOM + JSONL** — established standards; no OWL reasoner required on the hot path.
3. **Assignment without forest imports** — structural SMARTS / formula delta / opaque tags only.
4. **Phase I / Phase II expected**; enzyme names excluded from primary output.

## Feedback samples

Regenerate the proposed-term panel after ontology/assignment edits:

```bash
cargo run -p xenosite-xrm --example sample_terms -- --write
```

Output: [`data/samples/round-001.md`](data/samples/round-001.md). Comment on that file (or open an issue citing the case label) to steer the next round.

## Tests

```bash
cargo test -p xenosite-xrm
```

Includes schema validation, golden naming, enzyme-orthogonal guards, and dependency boundary checks.
