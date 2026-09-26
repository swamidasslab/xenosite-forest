# Competency evaluation plan (executable questions)

Status: **implementing** against current XRM layouts (YAML/SKOS JSON-LD + JSONL
assignments + SSSOM). TTL/OWL/SHACL/SPARQL paths below are the target shape;
Rust/JSON equivalents run in CI today where SPARQL is not yet wired.

## Layers

1. **Ontology competency** — Can the vocabulary express the concept?
2. **Mapping competency** — Can it map cleanly to external ontologies/models?
3. **Tagging competency** — Given reactant/product/rule, does the system emit
   the right tags and localized names?

## Target repo shape (from design note)

```
ontology/
  xmet.ttl / xmet.skos.ttl / xmet.shacl.ttl
  mappings/
tests/
  competency/cq.yml + sparql/
  fixtures/reactions.jsonl + expected_tags.jsonl
  gold/curated_reactions.tsv
scripts/
  run_competency_tests.py
  score_gold_set.py
```

## Current XRM adaptation

| Target | Current path |
| --- | --- |
| SKOS thesaurus | `data/ontology/xrm.yaml` → `xrm.skos.jsonld` |
| Mappings | `data/mappings/*.sssom.tsv` |
| Competency questions | `data/competency/cq.yml` |
| Tagging fixtures | `data/competency/fixtures/reactions.jsonl` |
| Runner | `tools/run_competency_tests.py` + Rust `tests/competency.rs` |
| Forest validation | `tools/harvest_forest_smarts.py` + `validate_forest_coverage.py` |

## Gold-set scoring (per spine)

`transformation_exact`, `transformation_ancestor_ok`, `phase_exact`,
`product_class_exact`, `site_exact`, `site_equivalent_under_symmetry`,
`site_motif_ok`, `liability_exact`, `external_mapping_ok`

## Initial CQ budget (~50)

| Area | Count |
| --- | --- |
| Phase I transformations | 10 |
| Phase II conjugations | 10 |
| Reactive metabolite classes | 8 |
| Site localization | 8 |
| Structural deltas | 5 |
| External mappings | 5 |
| Evidence/provenance | 4 |

Every CQ is a regression test: if the ontology cannot answer it, add terms/links
or mark the question out of scope.
