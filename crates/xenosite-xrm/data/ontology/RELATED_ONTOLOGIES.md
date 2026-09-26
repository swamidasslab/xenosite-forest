# Related ontologies and XRM positioning

## Positioning

XRM is a **med-chem-oriented SKOS thesaurus** for xenobiotic metabolism reaction
naming, site-localized reaction labels, structural transformation tags,
product/liability classes, and mappings to existing biochemical and chemical
ontologies.

It fills a real gap: **not** an enzyme ontology, pathway ontology, compound
ontology, or named-synthesis ontology. It is the missing layer between
reaction-generating SMARTS/rules and human-useful med-chem metabolism names.

**RXNO + MOP** are the closest actual reaction ontologies, but they are
synthesis / name-reaction oriented, not xenobiotic-metabolism oriented. XRM
stays **separate and SKOS-first**, with explicit SSSOM mappings to
RXNO / MOP / GO / ChEBI / Rhea / KEGG / MeSH (and ECO / CHMO for evidence),
rather than extending any one of them. Do **not** grow a full chemical
ontology here — metabolism-specific terms only ([`SCOPE.md`](SCOPE.md)).

## Closest resources (ordered)

| Resource | Why it matters | Fit |
| --- | --- | --- |
| **RXNO** (Name Reaction Ontology) | Closest to “reaction naming.” Connects named organic reactions to roles in synthesis and to molecular processes in MOP. Useful for naming architecture and external mappings, but not xenobiotic/med-chem specific. | High structural fit, low domain fit |
| **MOP** (Molecular Process Ontology) | Closest for generic molecular processes: cyclization, methylation, demethylation, addition, etc. Best source for broad transformation parents. | High chemistry fit |
| **GO** xenobiotic metabolic process (`GO:0006805`) | Best existing biological anchor for xenobiotic metabolism (includes synonyms like drug metabolism). Process biology, not reaction naming. | High domain fit, low reaction-detail fit |
| **ChEBI** | Essential for product classes, conjugate groups, chemical roles, and small-molecule grounding. | Essential companion |
| **Rhea** | Best curated biochemical reaction knowledgebase; ChEBI participants; UniProtKB enzyme/transporter annotation. | Strong reaction grounding, not med-chem naming |
| **KEGG REACTION / RCLASS / Reaction Modules** | RCLASS / modules are based on chemical structure transformation patterns rather than enzyme identity. Good transformation crosswalks. | High operational relevance |
| **MeSH** | Literature indexing and synonyms (metabolism, biotransformation, hydroxylation, dealkylation, toxicokinetics). Too coarse for structural reaction tagging alone. | Literature/search fit |
| **ECO** (Evidence & Conclusion Ontology) | Evidence provenance: experimentally supported assertion, computational assertion, curator inference, etc. | Evidence spine |
| **CHMO** (Chemical Methods Ontology) | Assay/instrument/method evidence: MS, NMR, chromatography, spectroscopy. | Evidence/method spine |
| **PROV-O / Dublin Core / SIO-like patterns** | Provenance, versioning, assertions, sources, generated-by. | Infrastructure |

## Mapping policy

| Predicate | Use |
| --- | --- |
| `skos:exactMatch` | True identity only (rare) |
| `skos:closeMatch` | Similar but not identical MeSH/KEGG/Rhea/MOP/RXNO terms |
| `skos:broadMatch` / `skos:narrowMatch` | Coarser/finer external terms |
| `skos:relatedMatch` | Loose association (including opaque `forest.*` CURIEs) |

Checked-in SSSOM:

| File | Targets |
| --- | --- |
| [`../mappings/xrm-mop.sssom.tsv`](../mappings/xrm-mop.sssom.tsv) | MOP process leaves |
| [`../mappings/xrm-mesh.sssom.tsv`](../mappings/xrm-mesh.sssom.tsv) | MeSH literature descriptors |
| [`../mappings/xrm-forest.sssom.tsv`](../mappings/xrm-forest.sssom.tsv) | Opaque Metabolic Forest CURIEs (alias spine) |
| [`../mappings/xrm-external.sssom.tsv`](../mappings/xrm-external.sssom.tsv) | GO, RXNO, ChEBI, Rhea, KEGG RCLASS, ECO, CHMO |

Harvest tooling pulls ChEBI / KEGG / Rhea / GO / Reactome candidates; promote
reviewed rows into SSSOM rather than inventing parallel hierarchies.
