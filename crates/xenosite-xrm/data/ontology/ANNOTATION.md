# Annotation bundles and link types

SKOS concepts **name** things. Reaction annotations **combine** them.
Do not mint combinatorial concepts such as
`benzylic_phase_1_hydroxylation_clearance_liability`.

Site-localized display names (`C4 aromatic hydroxylation`, `para-hydroxylation
of anisole`) are generated from templates + `SiteRef`, not stored as ontology
terms.

## Inventory guidance (v0.1)

A useful first band is **~350–500** chemist-facing canonical concepts with
**1,000–2,000** labels/synonyms/mappings — enough to be useful, still curatable.
That band is a **floor for usefulness**, not a hard cap. Extra terms are welcome
when they are well motivated and clear (distinct chemist handle, SMARTS/delta
assignable, or a needed med-chem / product / site / evidence facet). Prefer
clear leaves over combinatorial compounds; site-localized display names stay as
templates. Forest-map rules/patterns are an alias spine counted separately.

## SKOS link types

| Link | Use |
| --- | --- |
| `skos:broader` / `skos:narrower` | Within a spine |
| `skos:relatedMatch` | Loose cross-spine association |
| `skos:exactMatch` | True identity to an external concept (rare) |
| `skos:closeMatch` | Similar MeSH/KEGG/Rhea/model term |
| `skos:broadMatch` / `skos:narrowMatch` | Legacy/model terms broader/narrower than XRM |

## Operational annotation properties (`xmet:`)

These are **assertion / event** links, not thesaurus identity. Emit them on
named reactions (assignment hits), not as SKOS concept parents.

| Property | Use |
| --- | --- |
| `xmet:hasPhase` | → metabolism phase |
| `xmet:hasTransformation` | → chemical transformation |
| `xmet:hasProductClass` | → reactive / product metabolite class |
| `xmet:hasSiteEnvironment` | → site type / environment |
| `xmet:hasStructuralDelta` | → structural delta |
| `xmet:hasMedChemInterpretation` | → medchem liability |
| `xmet:hasEvidenceType` | → evidence |
| `xmet:hasBiologicalContext` | → biological context (orthogonal) |
| `xmet:generatedByRule` | → rule provenance / SMARTS / Forest rule |
| `xmet:mapsModelOutput` | → legacy / Rainbow / model output |
| `xmet:localizesToSite` | → `SiteRef` / map nums |
| `xmet:mayPrecede` / `xmet:mayFollow` | pathway logic |
| `xmet:bioactivatesTo` / `xmet:detoxifiesTo` | liability product framing |
| `xmet:hasConjugateGroup` | glucuronide / sulfate / GSH / … |
| `xmet:hasAttachmentAtomType` | O / N / S / C / acyl |

## Bundle shape (illustrative)

```json
{
  "transformation": "xrm:0000106",
  "phase": "xrm:0000001",
  "site_environment": ["xrm:1600014"],
  "structural_delta": ["xrm:1700012"],
  "medchem_interpretation": ["xrm:1400010"],
  "site_label": "C7 benzylic hydroxylation"
}
```

Today the Rust `Term` list is a flat multi-tag emission with ancestor
expansion; bundles can be layered later without expanding the concept count.
