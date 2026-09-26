# Sources for XRM term inventory

Chemist-facing subclasses and cross-cutting process facets were drawn from
expert xenobiotic-metabolism vocabulary and Swamidass-lab papers, **not** from
this repository’s rule-implementation abbreviations.

## Metabolic Forest hierarchy spine

The Phase I portion intended to map to Metabolic Forest keeps the five
Rainbow / Forest **class** parents under `phase I`, using **unabbreviated**
prefLabels only:

| XRM prefLabel | Forest opaque CURIE (SSSOM) | Paper class |
| --- | --- | --- |
| stable oxygenation | `forest.ruleset:SO` | Rainbow red |
| unstable oxygenation | `forest.ruleset:UO` | Rainbow orange |
| dehydrogenation | `forest.ruleset:DH` | Rainbow green |
| hydrolysis | `forest.ruleset:HD` | Rainbow blue |
| reduction | `forest.ruleset:RD` | Rainbow purple |

Localized abbreviations (`SO`, `UO`, `DH`, `HD`, `RD`) and CamelCase
identifiers (`StableOxygenation`, …) are **not** XRM prefLabels or altLabels.
They may appear only as opaque `forest.*` CURIE object ids in SSSOM.

## Primary literature

1. **Rainbow** — Dang, Matlock, Hughes, Swamidass. *The Metabolic Rainbow:
   Deep Learning Phase I Metabolism in Five Colors*. JCIM 2020.
   DOI [10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836) /
   PMC [PMC8716320](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8716320/).
   Five classes + **21 reaction types** (aromatic/aliphatic hydroxylation;
   aromatic/aliphatic epoxidation; N-/S-oxidation; N-/O-/S-/C-dealkylation;
   oxidative deamination & dehalogenation; alcohol→carbonyl; single→double /
   double→triple bond; quinone & iminium formation; ester/amide/ether/cyanide
   hydrolysis; carbonyl/nitro/sulfo reduction; reductive dehalogenation;
   hydrogenation). Also cues the excluded remainder (tautomerization,
   isomerization, rearrangement, radical formation, hydration, deacylation,
   denitrogenation, decarbonylation).
2. **Metabolic Forest** — Hughes, Dang, Swamidass (and related). Metabolite
   structure enumeration with rulesets for the five Phase I classes,
   **conjugation** (acetylation, glucuronidation, glutathionation, sulfation),
   **quinone formation**, and **tautomerization**.
3. **Quinone formation** — Hughes, Miller, Swamidass. Computational prediction
   of quinone formation. Chem. Res. Toxicol. (quinone-formation model).
   Species: quinone, quinone-imine, quinone-methide, imine-methide; routes:
   **one-step** vs **two-step** quinone formation; bioactivation framing.
4. **IUPAC** — *Glossary and tutorial of xenobiotic metabolism terms* (Pure
   Appl. Chem. 2021, DOI [10.1515/pac-2018-0208](https://doi.org/10.1515/pac-2018-0208)).
5. **Medicinal chemistry / DMPK teaching notes** — arene oxide / NIH shift;
   N-/O-/S-dealkylation and oxidative deamination via α-carbon hydroxylation;
   ω / ω−1, allylic, benzylic hydroxylation.
6. **Phase II literature** — O-/N-/S-/C-glucuronidation; acyl glucuronides;
   phenolic vs alcoholic sulfation; GSH Michael / epoxide / halide paths;
   mercapturic acid; amino-acid conjugation; methylation.
7. **Cross-cutting process facets** — oxidative **dearomatization**,
   **rearomatization**, **NIH shift**, **ipso substitution**, conjugate
   addition, acyl migration.
