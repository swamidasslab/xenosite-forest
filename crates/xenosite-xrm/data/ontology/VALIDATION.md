# Forest SMARTS validation (no namer coupling)

Forest rule SMIRKS are a **validation and candidate** source for XRM. The namer
must never import `xenosite-forest` or RDKit Forest code.

## Pipeline

```bash
# 1. Static text harvest from src/xenosite/forest/rules.py
python3 crates/xenosite-xrm/tools/harvest_forest_smarts.py
# → data/candidates/forest-smarts.jsonl

# 2. Coverage vs SSSOM + assignment tags
python3 crates/xenosite-xrm/tools/validate_forest_coverage.py
# → data/candidates/forest-coverage.json
```

## What validation asserts

1. Harvest is checked in and well-formed (SMIRKS, opaque `forest.*` tags).
2. Coverage report exists; most Forest rule classes have SSSOM and/or
   assignment-tag pathways into XRM.
3. For each SSSOM-mapped `forest.rule:*` / `forest.pattern:*` from the harvest,
   `Namer::name_smiles(..., &[tag])` succeeds (opaque tag path) — Forest code
   is not executed.
4. Crate `Cargo.toml` dependencies do not include forest crates.

## Competency evaluation (executable CQs)

```bash
pip install -r crates/xenosite-xrm/tools/requirements-competency.txt
python3 crates/xenosite-xrm/tools/run_competency_tests.py --refresh-ttl
cargo test -p xenosite-xrm --test competency
cargo test -p xenosite-xrm --test gold_score
python3 crates/xenosite-xrm/tools/score_gold_set.py
```

See `data/competency/PLAN.md`. SPARQL runs against `xrm.skos.ttl`; SHACL against
`xrm.shacl.ttl`. Tagging CQs export to `fixtures/reactions.jsonl` for Rust.

## Promoting SMARTS into assignments

Reviewed rows from `forest-smarts.jsonl` may become structural
`reactant_smarts` / `product_smarts` in `xenobiotic.jsonl` **after** chemist
simplification. Prefer simpler XRM SMARTS for naming; keep dense Forest SMIRKS
as the oracle for “does this Forest pattern family have an XRM handle?”
