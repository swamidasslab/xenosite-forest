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
| Thesaurus (authoring) | YAML | `data/ontology/xrm.yaml` |
| Thesaurus (runtime) | [SKOS](https://www.w3.org/TR/skos-reference/) JSON-LD | `data/ontology/xrm.skos.jsonld` |
| Crosswalks | [SSSOM](https://mapping-commons.github.io/sssom/) TSV | `data/mappings/*.sssom.tsv` |
| Assignment | JSON Lines | `data/assignments/xenobiotic.jsonl` |
| Bundle | JSON manifest | `data/manifest.json` |

Edit `xrm.yaml`, then export SKOS:

```bash
python3 crates/xenosite-xrm/tools/yaml_to_skos.py
```

XRM is a **separate** med-chem SKOS thesaurus (not an extension of RXNO, MOP,
MeSH, GO, or KEGG). Those resources — plus ChEBI, Rhea, ECO, CHMO — are linked
via SSSOM. See [`data/ontology/RELATED_ONTOLOGIES.md`](data/ontology/RELATED_ONTOLOGIES.md).
Enzyme types are orthogonal (biological context) and never primary reaction labels.

## Usage

```rust
use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

let namer = Namer::from_manifest(DEFAULT_MANIFEST)?;
let terms = namer.name_smiles("CC", "CCO", &[])?;
for t in terms {
    println!("{} ({}) path={:?}", t.pref_label, t.id, t.path_labels);
}

// Faceted annotation bundles (combinations are not ontology concepts):
let bundles = namer.annotate_smiles("CC", "CCO", &[])?;
for b in bundles {
    println!("{:?}  site_label={:?}", b.transformation, b.site_label);
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
3. **Assignment without forest imports** — prefer structural SMARTS / formula delta; opaque tags only when structure cannot decide (and for Forest-map correspondence).
4. **Cross-cutting multi-spine tagging** — med-chem facets (metabolism phase, chemical transformation, Rainbow phase I family, phase II conjugation, medchem liability, reactive metabolite / product class, site type, structural delta, leaving group, product status, rule provenance, evidence, biological context) plus ambiguity and Metabolic Forest map as an alias spine. Combinations use annotation bundles, not combinatorial concepts.
5. **Site-localized terms** — caller tags may use `@map` (`chem:hydroxylation@1`); SMARTS hits attach `SiteRef` so multi-change cases disambiguate.
6. **Phase I / Phase II expected**; enzyme names excluded from primary output.

## Candidate harvest (offline Python)

```bash
python3 crates/xenosite-xrm/tools/harvest/harvest_candidates.py \
  --out crates/xenosite-xrm/data/candidates/round-001.jsonl
```

Sources: seed, ChEBI, PubChem, KEGG, Rhea, GO, Reactome. See
[`tools/harvest/GUIDANCE_SOURCES.md`](tools/harvest/GUIDANCE_SOURCES.md).

## Feedback samples

Regenerate the proposed-term panel after ontology/assignment edits:

```bash
cargo run -p xenosite-xrm --example sample_terms -- --write
```

Output: [`data/samples/round-005.md`](data/samples/round-005.md) (and prior rounds under `data/samples/`). Comment on that file (or open an issue citing the case label) to steer the next round.

Reactions are tagged with **many cross-cutting terms** from the spines in
[`data/ontology/SOURCES.md`](data/ontology/SOURCES.md). Med-chem style stacks
look like `phase I + unstable oxygenation + N-dealkylation + aldehyde forming +
bioactivation risk + tertiary amine site`.

Forest-map ruleset parents use **unabbreviated** labels ending in `ruleset`.
Abbreviations appear only as opaque `forest.ruleset:*` CURIE object ids in SSSOM.

## Tests

```bash
cargo test -p xenosite-xrm

# Refresh Forest SMIRKS harvest + coverage (validation only; namer does not import forest):
python3 crates/xenosite-xrm/tools/harvest_forest_smarts.py
```

Includes schema validation, golden naming, Forest opaque-tag coverage, competency
fixtures, enzyme-orthogonal guards, and dependency boundary checks.

```bash
pip install -r crates/xenosite-xrm/tools/requirements-competency.txt
# Optional: regenerate ~590 CQs + 150 gold rows from ontology/assignments
python3 crates/xenosite-xrm/tools/expand_competency.py
python3 crates/xenosite-xrm/tools/run_competency_tests.py --refresh-ttl
cargo test -p xenosite-xrm --test competency
cargo test -p xenosite-xrm --test gold_score
python3 crates/xenosite-xrm/tools/score_gold_set.py
```

See `data/competency/PLAN.md` (executable CQ design), `data/ontology/VALIDATION.md`,
and `data/ontology/SCOPE.md` (metabolism-only).
