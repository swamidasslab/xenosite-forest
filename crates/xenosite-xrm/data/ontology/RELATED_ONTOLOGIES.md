# Related ontologies and XRM positioning

XRM is a **med-chem-oriented SKOS thesaurus** for xenobiotic metabolism reaction
naming, site-localized labels, structural transformation tags, product/liability
classes, and mappings to existing chemical/biochemical ontologies.

It fills the gap between reaction-generating SMARTS/rules and human-useful
med-chem metabolism names. It is **not** an enzyme ontology, pathway ontology,
compound ontology, or named-synthesis ontology — and it is **not an extension**
of any single external resource. Crosswalks use SSSOM.

## Closest resources (ordered)

| Resource | Why it matters | Fit |
| --- | --- | --- |
| **RXNO** (Name Reaction Ontology) | Closest to “reaction naming”; links named organic reactions to roles and to MOP | High structural fit, low xenobiotic/med-chem domain fit |
| **MOP** (Molecular Process Ontology) | Best source for broad transformation parents (methylation, demethylation, addition, …) | High chemistry fit |
| **GO** xenobiotic metabolic process (`GO:0006805`) | Best biological anchor for xenobiotic / drug metabolism | High domain fit, low reaction-detail fit |
| **ChEBI** | Product classes, conjugate groups, chemical roles, small-molecule grounding | Essential companion |
| **Rhea** | Curated biochemical reactions with ChEBI participants | Strong reaction grounding, not med-chem naming |
| **KEGG REACTION / RCLASS / modules** | Structure-transformation patterns (enzyme-independent) | High operational relevance |
| **MeSH** | Literature synonyms (biotransformation, hydroxylation, …) | Literature/search fit; too coarse alone |
| **ECO** | Evidence / assertion methods for biocuration | Evidence spine |
| **CHMO** | Assay/instrument methods (MS, NMR, chromatography) | Evidence/method spine |
| **PROV-O / Dublin Core** | Provenance, versioning, generated-by | Infrastructure |

## Mapping policy

| Predicate | Use |
| --- | --- |
| `skos:exactMatch` | True identity only (rare) |
| `skos:closeMatch` | Similar but not identical MeSH/KEGG/Rhea/MOP/RXNO terms |
| `skos:broadMatch` / `skos:narrowMatch` | Coarser/finer external terms |
| `skos:relatedMatch` | Loose association (including opaque `forest.*` CURIEs) |

Checked-in SSSOM files:

- [`../mappings/xrm-mop.sssom.tsv`](../mappings/xrm-mop.sssom.tsv)
- [`../mappings/xrm-mesh.sssom.tsv`](../mappings/xrm-mesh.sssom.tsv)
- [`../mappings/xrm-forest.sssom.tsv`](../mappings/xrm-forest.sssom.tsv)
- [`../mappings/xrm-external.sssom.tsv`](../mappings/xrm-external.sssom.tsv) (GO, RXNO, ECO, ChEBI seeds)

Harvest tooling already pulls ChEBI / KEGG / Rhea / GO / Reactome candidates;
promote reviewed rows into SSSOM rather than inventing parallel hierarchies.
