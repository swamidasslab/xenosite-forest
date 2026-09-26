# Competency evaluation plan (executable questions)

Status: **executed** (executable CQs in CI; ~10× expanded; gold set scored; SHACL definitions complete).

Store of the design note (metabolism-specific XRM adaptation). Every competency
question is a regression test: if the ontology cannot answer it, add terms/links
or mark the question out of scope.

## Layers

1. **Ontology competency** — Can the vocabulary express the concept?
2. **Mapping competency** — Can it map cleanly to external ontologies/models?
3. **Tagging competency** — Given reactant/product/rule, does the system emit
   the right tags and localized names?

## Design note → XRM paths

```
# Design note                          # XRM adaptation
ontology/
  xmet.ttl / xmet.skos.ttl      →  data/ontology/xrm.skos.ttl (+ xrm.skos.jsonld)
  xmet.shacl.ttl                →  data/ontology/xrm.shacl.ttl
  mappings/mesh.tsv             →  data/mappings/views/mesh.tsv
  mappings/kegg.tsv             →  data/mappings/views/kegg.tsv
  mappings/rxno_mop.tsv         →  data/mappings/views/rxno_mop.tsv
tests/competency/
  cq.yml                        →  data/competency/cq.yml
  sparql/*.rq                   →  data/competency/sparql/
  fixtures/reactions.jsonl      →  data/competency/fixtures/reactions.jsonl
  fixtures/expected_tags.jsonl  →  data/competency/fixtures/expected_tags.jsonl
  gold/curated_reactions.tsv    →  data/competency/gold/curated_reactions.tsv
scripts/
  run_competency_tests.py       →  tools/run_competency_tests.py
  score_gold_set.py             →  tools/score_gold_set.py
```

Authoring remains YAML (`xrm.yaml`) → JSON-LD → Turtle. SSSOM TSVs stay the
canonical mapping store; `views/` are compact projections for CQ readability.

## CQ types

| Type | Layer | Runner |
| --- | --- | --- |
| `ontology_labels` / `ontology_broader` | ontology | Python over JSON-LD |
| `sparql` | ontology | rdflib over `xrm.skos.ttl` |
| `shacl` | ontology | pyshacl over `xrm.shacl.ttl` |
| `definition_coverage` | ontology | Python fraction check |
| `mapping_exists` | mapping | SSSOM TSV scan |
| `tagging` | tagging | exported fixtures → Rust `tests/competency.rs` |

## Gold-set scoring (per spine)

`transformation_exact`, `transformation_ancestor_ok`, `phase_exact`,
`product_class_exact`, `site_exact`, `site_equivalent_under_symmetry`,
`site_motif_ok`, `liability_exact`, `external_mapping_ok`

```bash
cargo test -p xenosite-xrm --test gold_score
python3 crates/xenosite-xrm/tools/score_gold_set.py
# → data/competency/gold/last_score.json + F1 report
```

## CQ budget (~10× initial ~50)

Regenerate with `python3 tools/expand_competency.py`.

| Area | Target | Notes |
| --- | --- | --- |
| Phase I transformations | 100 | SPARQL + broader + tagging mix |
| Phase II conjugations | 100 | attachment-atom SPARQL seeds retained |
| Reactive metabolite classes | 80 | quinone-like, GSH-trappable |
| Site localization | 80 | bond-centered epoxidation site_label |
| Structural deltas | 50 | mass shifts / oxygenation |
| External mappings | 50 | MeSH / MOP / GO / Forest / KEGG |
| Evidence | 20 | definition coverage + MS/literature |
| Provenance | 20 | includes SHACL |
| Leaving group | 25 | methyl/halide/… |
| Pharmacological role | 15 | active/inactive/prodrug |
| Annotation about | 10 | parent/product/reaction |
| Med-chem liability | 40 | soft spots / bioactivation |

Gold set target: **150** curated rows (`tools/expand_competency.py`).

## Run

```bash
pip install -r crates/xenosite-xrm/tools/requirements-competency.txt
python3 crates/xenosite-xrm/tools/run_competency_tests.py --refresh-ttl
cargo test -p xenosite-xrm --test competency
cargo test -p xenosite-xrm --test gold_score
python3 crates/xenosite-xrm/tools/score_gold_set.py
```

## SHACL

Core shape requires `skos:prefLabel` + `skos:inScheme` on every concept.
Definition completeness is a warning shape; `CQ-EV-002` enforces coverage ≥ 80%.
