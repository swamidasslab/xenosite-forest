# Sources for XRM term inventory

XRM **tags each reaction with many terms** from parallel, cross-cutting spines.
That improves chemist readability and makes soft alignment to external
ontologies (MeSH, MOP, Forest opaque CURIEs) easier: each spine can map
independently.

Forest abbreviations (`SO`, `UO`, `DH`, `HD`, `RD`, …) are **never** XRM
prefLabels or altLabels. They appear only as opaque `forest.*` CURIE object
ids in SSSOM.

## Parallel spines (under `xenobiotic biotransformation`)

| Spine | Role | Auto-tag cues |
| --- | --- | --- |
| Chemist reaction type (`phase I` / `phase II`, named transformations) | Primary chemist vocabulary | SMARTS, delta, `chem:*` tags |
| **Metabolic Forest map** | Intentionally mirrors repo rulesets → rules → PatternInfo | `forest.rule:*`, `forest.pattern:*`, `forest.ruleset:*` |
| **Aromatic and conjugated-system impact** | π-system fate (dearomatize, quinoid, arene oxide, …) | quinone / arene-oxide / GSH-Michael / dearomatization tags |
| **Redox polarity** | Net oxidation / reduction / redox-neutral | delta O, typed oxidations/reductions/conjugations |
| **Site atom class** | C / N / O / S / halogen / multi-element site | typed dealkylation / oxidation tags |
| **Bond-edit topology** | Addition, cleavage, bond-order change, ring closure, … | SMARTS class of edit |
| **Metabolite cardinality** | Single metabolite vs fragmenting | hydroxylation vs dealkylation/hydrolysis |
| **Oxygenation outcome** | Stable oxygen addition vs oxygen-triggered cleavage | chemist parallel to Forest stable/unstable classes |
| **Electrophile role** | Generate / consume electrophile; nucleophile exposure | quinone, epoxide, GSH, bioactivation tags |
| **Ring fate** | Preserved / opened / formed / resized | epoxide, ring-opening facets |
| **Formula-delta class** | `+O`, −halogen, −C, … | elemental delta only |
| **Site aromaticity** | Aromatic vs aliphatic site | `site_aromatic` on typed rules |
| **Pathway-step role** | One-step / multi-step / preparatory / terminal | quinone one-/two-step and path tags |
| **Ambiguity and underspecification** | Typed incomplete/conflicting evidence | `chem:*-ambiguity`, `chem:*-underspecified` |
| Process facet | NIH shift, ipso, carbinolamine cleavage, acyl migration, … | mechanism tags |
| Bioactivation / detoxication | Toxicity-oriented outcome framing | reactive-metabolite tags |

Ambiguity subtypes include site/regio/stereo SOM ambiguity, reaction-type and
competing-type ambiguity, mechanism ambiguity, metabolite-structure and
atom-mapping underspecification, formula-only evidence, pathway-depth and
intermediate gaps, phase ambiguity, Forest-map correspondence ambiguity,
external-ontology alignment ambiguity, aromatic-impact and electrophile-role
ambiguity, and provenance underspecification — plus `fully specified` when
callers assert completeness.

## Metabolic Forest map (full names)

| XRM prefLabel | Opaque Forest CURIE |
| --- | --- |
| stable oxygenation | `forest.ruleset:SO` |
| unstable oxygenation | `forest.ruleset:UO` |
| dehydrogenation ruleset | `forest.ruleset:DH` |
| hydrolysis ruleset | `forest.ruleset:HD` |
| reduction ruleset | `forest.ruleset:RD` |
| quinone formation ruleset | `forest.ruleset:QF` |
| conjugation ruleset | `forest.ruleset:CJ` |
| tautomerization ruleset | `forest.ruleset:TT` |
| phase I ruleset | `forest.ruleset:PhaseOne` |
| bioactivation ruleset | `forest.ruleset:BA` |

Rules and patterns hang under those rulesets (`Hydroxylation rule`,
`Hydroxylation/h`, …) with `skos:exactMatch` to `forest.rule:*` /
`forest.pattern:*`.

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
