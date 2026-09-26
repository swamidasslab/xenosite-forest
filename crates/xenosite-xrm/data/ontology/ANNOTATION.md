# Annotation bundles and link types

SKOS concepts **name** things. Reaction annotations **combine** them.
Do not mint combinatorial concepts such as
`benzylic_phase_1_hydroxylation_clearance_liability`.

Site-localized display names (`C4 aromatic hydroxylation`, `benzylic hydroxylation @1`)
are generated from templates + `SiteRef` into
[`AnnotationBundle::site_label`](../../src/bundle.rs), not stored as ontology terms.

```rust
let bundles = namer.annotate_smiles("CC", "CCO", &[])?;
// bundles[i].transformation / phase / site_environment / medchem_interpretation / site_label
```

## Inventory guidance (v0.1)

A useful first band is **~350–500** chemist-facing canonical concepts with
**1,000–2,000** labels/synonyms/mappings — enough to be useful, still curatable.
That band is a **floor for usefulness**, not a hard cap. Extra terms are welcome
when they are well motivated and clear. Prefer clear leaves over combinatorial
compounds. Forest-map rules/patterns are an alias spine counted separately.

### Category bands (primary-spine counting)

Count a concept toward the spine of its **nearest spine-root ancestor** (not the
full polyhierarchy descendant closure — that double-counts leaves under both
chemical transformation and phase I family).

| Spine | Target | Notes |
| --- | --- | --- |
| Phase / pathway role | 10–20 | phase I/II, sequential metabolism, intermediate formation |
| Chemical transformation | 90–140 | Largest core |
| Phase I detailed classes | 25–40 | Rainbow-compatible + finer types |
| Phase II conjugation classes | 40–70 | By attachment atom where useful |
| Product / metabolite class | 70–110 | Quinone, epoxide, GSH adduct, … |
| Site environment | 60–100 | Benzylic, phenol, tertiary amine, heteroarene, … |
| Structural delta | 35–60 | +O, −2H, aromaticity loss, conjugate mass shift |
| Med-chem interpretation | 40–70 | Soft spot, bioactivation, blocking, … |
| Evidence / assertion | 30–50 | MS, NMR, predicted/observed/curated/conflict |
| Biological context | 40–80 | Enzyme/tissue/species/matrix (orthogonal) |
| Rule/model provenance | 20–40 | SMARTS, Rainbow class, XenoNet, legacy model |
| Localization templates | 20–40 | Templates only — not per-atom enumerations |
| Leaving group | 20–40 | Departing methyl/ethyl/halide/carboxylate/… fragments |
| Pharmacological role | 10–25 | Active/inactive metabolite, prodrug (metabolism-specific) |
| Annotation about | 5–15 | about parent / product / reaction (not combinatorial) |

## SKOS link types

| Link | Use |
| --- | --- |
| `skos:broader` / `skos:narrower` | Within a spine |
| `skos:related` / `skos:relatedMatch` | Loose cross-spine association |
| `skos:exactMatch` | True identity to an external concept (rare) |
| `skos:closeMatch` | Similar MeSH/KEGG/Rhea/model term |
| `skos:broadMatch` / `skos:narrowMatch` | Legacy/model terms broader/narrower than XRM |

## Operational annotation properties (`xmet:`)

These are **assertion / event** links on [`AnnotationBundle`], not thesaurus
identity parents.

| Property | Bundle field |
| --- | --- |
| `xmet:hasPhase` | `phase` |
| `xmet:hasTransformation` | `transformation` |
| `xmet:hasProductClass` | `product_class` |
| `xmet:hasSiteEnvironment` | `site_environment` |
| `xmet:hasStructuralDelta` | `structural_delta` |
| `xmet:hasMedChemInterpretation` | `medchem_interpretation` |
| `xmet:hasEvidenceType` | `evidence_type` |
| `xmet:hasBiologicalContext` | `biological_context` |
| `xmet:generatedByRule` | `generated_by_rule` |
| `xmet:mapsModelOutput` | (via SSSOM / provenance leaves) |
| `xmet:localizesToSite` | `site` + `site_label` |
| `xmet:mayPrecede` / `xmet:mayFollow` | SKOS `relatedMatch` pathway pairs |
| `xmet:bioactivatesTo` / `xmet:detoxifiesTo` | medchem ↔ product_class related pairs |
| `xmet:hasConjugateGroup` | phase2 / product_class conjugates |
| `xmet:hasAttachmentAtomType` | site environment attachment-atom leaves |
| `xmet:hasLeavingGroup` | `leaving_group` |
| `xmet:hasPharmacologicalRole` | `pharmacological_role` |
| `xmet:aboutParent` / `aboutProduct` / `aboutReaction` | distinguish tag target |

## Bundle shape

```json
{
  "transformation": "xrm:0000106",
  "phase": "xrm:0000001",
  "site_environment": ["xrm:1600014"],
  "structural_delta": ["xrm:1700012"],
  "medchem_interpretation": ["xrm:1400010"],
  "site_label": "benzylic hydroxylation @1",
  "site": {"map_nums": [1]},
  "xmet_properties": ["xmet:hasPhase", "xmet:hasTransformation", "xmet:localizesToSite"]
}
```
