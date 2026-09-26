# XRM candidate harvest (offline)

Python tooling to collect **candidate terms, synonyms, and example
reactant→product pairs** from public ontologies/databases. Output is JSONL for
human review before promotion into `data/ontology/` / `data/assignments/` /
SSSOM.

This does **not** import Metabolic Forest. Prefer structure-backed examples
(SMILES pairs, ChEBI/KEGG/Rhea ids) over reaction-tool tags.

## Sources

| Source | What we pull | API |
| --- | --- | --- |
| **ChEBI** | Synonyms, IUPAC names, parent links, SMILES | REST `chebi/backend/api` |
| **PubChem** | Example compound SMILES / titles by CID | PUG REST |
| **KEGG** | Reaction names containing xenobiotic keywords; RCLASS ids | REST `rest.kegg.jp` |
| **Rhea** | Reaction equations, ChEBI participants, EC, cross-refs | `rhea-db.org` TSV |
| **MeSH** | Phase I/II and biotransformation descriptors | NCBI E-utilities (optional) |
| **GO** | Xenobiotic metabolic process terms | QuickGO |
| **Reactome** | Pathway display names (glucuronidation, …) | Content service |
| **Seed lexicon** | Hand list of XRM-relevant reaction-type strings | local |

Also useful later (stubs documented, not all wired): MOP/OLS, NCIT, HMDB,
MetaCyc, UniProt enzyme names (enzyme labels stay orthogonal to primary XRM
prefLabels).

## Usage

Stdlib only (`urllib`, `json`, `csv`, `argparse`). Network required for live pulls.

```bash
# From repo root
python3 crates/xenosite-xrm/tools/harvest/harvest_candidates.py \
  --out crates/xenosite-xrm/data/candidates/round-001.jsonl \
  --sources seed,chebi,pubchem,kegg,rhea,go,reactome

# Forest SMIRKS/SMARTS (static text parse of rules.py — no RDKit, no namer import):
python3 crates/xenosite-xrm/tools/harvest_forest_smarts.py
# → data/candidates/forest-smarts.jsonl  (validation + promote-after-review)

# Faster / offline-ish: seed + cached fixtures only
python3 crates/xenosite-xrm/tools/harvest/harvest_candidates.py \
  --out crates/xenosite-xrm/data/candidates/seed-only.jsonl \
  --sources seed
```

Each JSONL record:

```json
{
  "candidate_id": "cand:chebi:hydroxylation:…",
  "pref_label": "aromatic hydroxylation",
  "synonyms": ["aryl hydroxylation", "…"],
  "definition": "…",
  "sources": [{"system": "chebi", "id": "CHEBI:…", "url": "…"}],
  "examples": [
    {
      "reactant_smiles": "c1ccccc1",
      "product_smiles": "Oc1ccccc1",
      "label": "benzene → phenol",
      "refs": [{"system": "pubchem", "id": "CID:241"}]
    }
  ],
  "suggested_spines": ["chemist reaction type", "redox polarity"],
  "notes": "…"
}
```

## Review workflow

1. Run harvest → `data/candidates/*.jsonl`
2. Diff / review; promote synonyms into SKOS `altLabel`, examples into
   `tests/naming_gold.rs` or `examples/sample_terms.rs`, external ids into SSSOM
3. Prefer SMARTS assignment rules for promoted examples; keep opaque tags for
   Forest-map correspondence only
