# Sources for XRM term inventory

XRM **tags each reaction with many terms** from parallel, cross-cutting spines.
Authoring source: [`xrm.yaml`](xrm.yaml). Runtime export: [`xrm.skos.jsonld`](xrm.skos.jsonld).
Annotation bundles / `xmet:` link types: [`ANNOTATION.md`](ANNOTATION.md).

v0.1 useful band: ~350–500 chemist-facing concepts with 1,000–2,000
labels/synonyms/mappings — a floor, not a cap. Extra well-motivated terms are
fine. Site-localized names are templates, not concepts. Combinations emit via
[`AnnotationBundle`](ANNOTATION.md) (`namer.annotate`), not new concepts.
Positioning vs RXNO/MOP/GO/ChEBI/…: [`RELATED_ONTOLOGIES.md`](RELATED_ONTOLOGIES.md).

Primary-spine inventory (nearest spine-root ancestor; Forest map separate) is
tracked against the category bands in [`ANNOTATION.md`](ANNOTATION.md).

Spine vocabulary follows the ChatGPT design share (metabolism phase, chemical
transformation, Rainbow phase I family, phase II conjugation family, medchem
liability, reactive metabolite family, site type, structural delta, product
status, rule provenance, evidence, biological context), plus ambiguity and the
Metabolic Forest map as an **alias / mapping** spine.

Forest abbreviations (`SO`, `UO`, `DH`, `HD`, `RD`, …) are **never** XRM
prefLabels. They appear only as opaque `forest.*` CURIE object ids in SSSOM.
Forest ruleset concepts use unabbreviated labels ending in `ruleset`.

## Core spines (under `xenobiotic biotransformation`)

| Spine | Role | Auto-tag cues |
| --- | --- | --- |
| **metabolism phase** | Phase I / II / III framing | ancestors of typed reactions; `chem:phase-*` |
| **chemical transformation** | Enzyme-independent edit class | SMARTS / delta / typed chem tags |
| **phase I reaction family** | Rainbow five classes + detailed Phase I types | hydroxylation→stable oxygenation; dealkylation→unstable oxygenation |
| **phase II conjugation family** | Glucuronidation, sulfation, GSH, … | conjugation SMARTS / tags |
| **medchem liability** | Soft spot, clearance, bioactivation, blocking | typed liability emits |
| **reactive metabolite family** | Quinone, epoxide, aldehyde, acyl glucuronide, … | quinone / epoxide / GSH tags |
| **site type** | Atom / bond / ring / aromatic / benzylic / … | site aromaticity + typed sites |
| **structural delta** | Formula, mass, bond order, aromaticity, … | elemental delta; nested legacy delta spines |
| **product status** | Observed / predicted / intermediate / … | caller / harvest tags |
| **rule provenance** | SMARTS vs delta vs caller vs Forest map | assignment path |
| **evidence** | Literature, MS, NMR, incubation system, … | harvest / curation tags |
| **biological context** | Enzyme / tissue / species / matrix (orthogonal) | never primary reaction label |
| **ambiguity and underspecification** | Typed incomplete/conflicting evidence | `chem:*-ambiguity` |
| **Metabolic Forest map** | Alias spine: ruleset → rule → PatternInfo | `forest.rule:*`, `forest.pattern:*` |

Legacy facets (redox polarity, bond-edit topology, formula-delta class, ring
fate, oxygenation outcome, metabolite cardinality, aromatic impact,
pathway-step role, site atom class, site aromaticity, process facet) hang
**under** `structural delta`, `site type`, or `chemical transformation` rather
than as peer root spines.

## Metabolic Forest map (full names)

| XRM prefLabel | Opaque Forest CURIE |
| --- | --- |
| stable oxygenation ruleset | `forest.ruleset:SO` |
| unstable oxygenation ruleset | `forest.ruleset:UO` |
| dehydrogenation ruleset | `forest.ruleset:DH` |
| hydrolysis ruleset | `forest.ruleset:HD` |
| reduction ruleset | `forest.ruleset:RD` |
| quinone formation ruleset | `forest.ruleset:QF` |
| conjugation ruleset | `forest.ruleset:CJ` |
| tautomerization ruleset | `forest.ruleset:TT` |
| phase I ruleset | `forest.ruleset:PhaseOne` |
| bioactivation ruleset | `forest.ruleset:BA` |

Chemist Rainbow classes `stable oxygenation` / `unstable oxygenation` are under
**phase I reaction family**, with `skos:relatedMatch` to the Forest rulesets.

## Primary literature

1. **Rainbow** — Dang et al., JCIM 2020
   ([10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836)).
2. **Metabolic Forest** — metabolite enumeration; rulesets for Phase I classes,
   conjugation, quinone formation, tautomerization.
3. **Quinone formation** — Hughes & Swamidass, Chem. Res. Toxicol. 2017
   ([10.1021/acs.chemrestox.6b00385](https://doi.org/10.1021/acs.chemrestox.6b00385)).
4. **IUPAC** xenobiotic metabolism glossary (Pure Appl. Chem. 2021,
   [10.1515/pac-2018-0208](https://doi.org/10.1515/pac-2018-0208)).
5. DMPK / Phase II teaching literature for conjugation and process facets.

## Automated external guidance

Offline harvesters (Python, stdlib) collect candidate terms, synonyms, and
example pairs from **ChEBI, PubChem, KEGG, Rhea, GO, Reactome**, plus a local
seed lexicon. See:

- [`tools/harvest/README.md`](../tools/harvest/README.md)
- [`tools/harvest/GUIDANCE_SOURCES.md`](../tools/harvest/GUIDANCE_SOURCES.md)
- Output: [`data/candidates/`](../candidates/)

Promote reviewed candidates into YAML `synonyms`, SSSOM, and SMARTS-backed
goldens. Prefer structure (SMARTS/delta) over Forest reaction-tool tags.
