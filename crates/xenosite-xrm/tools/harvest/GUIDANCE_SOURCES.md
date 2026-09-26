# Guidance sources for XRM terms, synonyms, and examples

Checked-in harvester: `harvest_candidates.py` (stdlib Python). Promote reviewed
rows from `data/candidates/*.jsonl` into SKOS / SSSOM / goldens.

## Automated now

| Source | Guidance offered | Automation |
| --- | --- | --- |
| Seed lexicon | Classic xenobiotic pairs + preferred labels | local |
| ChEBI | Synonyms, IUPAC names, metabolite SMILES, parent links | REST |
| PubChem | Example compound titles + ConnectivitySMILES by CID | PUG REST |
| KEGG REACTION | Named hydroxylation / demethylation / glucuronidation / … reactions | `rest.kegg.jp` |
| Rhea | Expert reaction equations, EC, ChEBI participants, KEGG/Reactome xrefs | TSV REST |
| GO (QuickGO) | Process-level xenobiotic / drug metabolism terms | QuickGO search |
| Reactome | Phase I / Phase II / glucuronidation / GSH pathway names | ContentService |

## Closest external ontologies (map, do not extend)

See [`data/ontology/RELATED_ONTOLOGIES.md`](../../data/ontology/RELATED_ONTOLOGIES.md).
Order of relevance: **RXNO** (naming architecture) → **MOP** (broad
transformation parents) → **GO:0006805** (xenobiotic process anchor) →
**ChEBI / Rhea / KEGG RCLASS** (participants & structure deltas) → **MeSH**
(literature) → **ECO / CHMO** (evidence & methods). XRM stays separate SKOS
with SSSOM crosswalks.

## Strong next targets (stub / manual)

| Source | Why | How to automate later |
| --- | --- | --- |
| **RXNO** | Named-reaction architecture; low domain fit | OLS / OBO |
| **MOP** via OLS | Broad transformation parents; partial SSSOM already | EMBL-EBI OLS API |
| **MeSH** | Phase I/II descriptors; biotransformation | NCBI E-utilities / MeSH RDF |
| **ECO / CHMO** | Evidence and assay/method provenance | OLS |
| **NCIt** | Cancer/drug metabolism wording | NCI EVS API |
| **HMDB** | Human metabolite structures + biotransformations | HMDB XML/API dumps |
| **MetaCyc / BioCyc** | Curated reaction frames with compounds | Pathway Tools web services / flat files |
| **UniProt** (Rhea-linked) | Enzyme→reaction, not preferred as XRM prefLabels | UniProt REST |
| **SBO** | Systems biology reaction-type terms | OLS |
| **IUPAC xenobiotic glossary** (Pure Appl. Chem. 2021) | Authoritative chemist vocabulary | Manual PDF/table extract → seed |
| **DrugBank** metabolites | Clinical xenobiotic examples | Licensed dump / scrape policy check |
| **PubChem Pathway / BioAssay** | Additional example pairs | PUG + annotations |
| **Swamidass papers** (Rainbow, Forest, quinone) | 21 types, quinone subtypes, Forest rulesets | Already mirrored in SKOS spines |
| **Forest `rules.py` SMIRKS** (static text harvest) | Reactant/product SMARTS + pattern names for validation / candidate assignments | `tools/harvest_forest_smarts.py` → `data/candidates/forest-smarts.jsonl` (namer never imports forest) |

## Design rules when promoting

1. **Cross-cutting spines** — map each candidate to one or more ChatGPT spines
   (phase, transformation, Rainbow family, conjugation, medchem liability,
   reactive metabolite, site, structural delta, evidence, …), not a single parent.
2. **SMARTS over reaction tools** — examples with SMILES should become
   assignment SMARTS/delta rules; Forest tags only for Forest-map correspondence.
3. **Site localization** — multi-site examples should use mapped SMILES and/or
   `tag@map` so terms attach to sites; do not mint per-atom ontology concepts.
4. **Enzyme names stay orthogonal** — EC/CYP strings may appear under biological
   context or SSSOM, not as primary reaction `prefLabel`.
5. **Ambiguity** — when sources disagree, emit typed ambiguity candidates rather
   than forcing one label.
6. **Inventory size** — useful band ~350–500 chemist concepts; extra clear,
   well-motivated terms are fine (no hard cap to hit a target).
