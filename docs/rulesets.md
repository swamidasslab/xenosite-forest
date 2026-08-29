# Built-in rulesets

Pass any **name** below to `bfs(..., ruleset=...)` or `load_ruleset(...)`. Aliases in parentheses are also registered.

These rules *enumerate structures*. The cited papers are the XenoSite models and labeling schemes the rulesets are written to match. Site-of-metabolism scores from those models are at [xenosite.org](https://xenosite.org).

Always cite Metabolic Forest when you use this package. Cite the matching paper as well when a specific ruleset or labeling scheme is central to the work.

## How they fit together

| Ruleset | Name | Rules | Matched publication |
| --- | --- | --- | --- |
| `Full` | `Full` | All groups below except `Bioactivation` / `TSO` as standalone sets | Metabolic Forest |
| `PhaseOneRS` | `PhaseOneRS` | Rainbow five colors: SO, DH, HD, RD, UO | Metabolic Rainbow |
| `QuinoneFormationRS` | `QF` (`QuinoneFormation`) | `QuinoneFormation` | Quinone formation |
| `Bioactivation` | `BA` (`BioactivationPathways`) | Quinone, epoxidation, nitroaromatic reduction, thiophene S-oxidation | Bioactivation |
| `StableOxygenationRS` | `SO` (`StableOxygenation`) | Epoxidation, hydroxylation, N-oxidation, S-oxidation | Rainbow (SO) |
| `DehydrogenationRS` | `DH` (`Dehydrogenation`) | Dehydrogenation | Rainbow (DH) |
| `HydrolysisRS` | `HD` (`Hydrolysis`) | Hydrolysis, epoxide opening, dephosphorylation, azo splitting | Rainbow (HD) |
| `ReductionRS` | `RD` (`Reduction`) | Hydrogenation, dehydration, N/S/O reduction, reductive dehalogenation, benzodioxole reduction | Rainbow (RD) |
| `UnstableOxygenationRS` | `UO` (`UnstableOxygenation`) | Dealkylation, oxidative dehalogenation | Rainbow (UO) |
| `ConjugationRS` | `CJ` (`Conjugation`) | Glucuronidation, sulfation, acetylation, glutathionation | UGT (glucuronidation) |
| `TautomerizationRS` | `TT` (`Tautomerization`) | Tautomerization | Metabolic Forest |
| `ThiopheneSulfurOxidationRS` | `TSO` (`ThiopheneSulfurOxidation`) | Thiophene S-oxidation | Bioactivation |

`PhaseOneRS` is the Rainbow-aligned Phase I set (no conjugations, quinone, or tautomerization). `Full` is the complete Metabolic Forest generator.

`PhaseOneRS` hydrolysis/reduction lists are slightly narrower than `HD`/`RD` in `Full` (no azo splitting or benzodioxole reduction). Use `PhaseOneRS` when you want the five-color Phase I scheme from the Rainbow paper.

## Phase I (Metabolic Rainbow)

The Rainbow paper labels Phase I metabolism in five reaction classes. `PhaseOneRS` implements that scheme:

| Color | Ruleset | Reaction types in this package |
| --- | --- | --- |
| SO | `SO` | Hydroxylation, epoxidation, nitrogen oxidation, sulfur oxidation |
| DH | `DH` | Dehydrogenation |
| HD | `HD` | Hydrolysis, epoxide opening, dephosphorylation |
| RD | `RD` | Hydrogenation, dehydration, nitrogen/sulfur/oxygen reduction, reductive dehalogenation |
| UO | `UO` | Dealkylation, oxidative dehalogenation |

Dang, N. L.; Matlock, M. K.; Hughes, T. B.; Swamidass, S. J.
The Metabolic Rainbow: Deep Learning Phase I Metabolism in Five Colors.
*J. Chem. Inf. Model.* **2020**, *60* (3), 1146–1164.
**DOI:** [10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836)

```bibtex
@article{Dang2020MetabolicRainbow,
  title = {The Metabolic Rainbow: Deep Learning Phase I Metabolism in Five Colors},
  author = {Dang, Na Le and Matlock, Matthew K. and Hughes, Tyler B. and Swamidass, S. Joshua},
  journal = {Journal of Chemical Information and Modeling},
  volume = {60},
  number = {3},
  pages = {1146--1164},
  year = {2020},
  doi = {10.1021/acs.jcim.9b00836},
  url = {https://doi.org/10.1021/acs.jcim.9b00836}
}
```

On the web: [xenosite.org/phase1](https://xenosite.org/phase1).

## Quinone formation

`QuinoneFormationRS` (`QF`) implements the specialized quinone / quinone-imine / quinone-methide structure algorithm described in Metabolic Forest and matched to the quinone XenoSite model.

Hughes, T. B.; Swamidass, S. J.
Deep Learning to Predict the Formation of Quinone Species in Drug Metabolism.
*Chem. Res. Toxicol.* **2017**, *30* (2), 642–656.
**DOI:** [10.1021/acs.chemrestox.6b00385](https://doi.org/10.1021/acs.chemrestox.6b00385)

```bibtex
@article{Hughes2017Quinone,
  title = {Deep Learning to Predict the Formation of Quinone Species in Drug Metabolism},
  author = {Hughes, Tyler B. and Swamidass, S. Joshua},
  journal = {Chemical Research in Toxicology},
  volume = {30},
  number = {2},
  pages = {642--656},
  year = {2017},
  doi = {10.1021/acs.chemrestox.6b00385},
  url = {https://doi.org/10.1021/acs.chemrestox.6b00385}
}
```

On the web: [xenosite.org/quinone](https://xenosite.org/quinone).

## Bioactivation

`Bioactivation` (`BA`) groups the four common bioactivation routes modeled jointly in the 2021 paper: quinone formation, epoxidation, nitroaromatic reduction, and thiophene sulfur oxidation.

Hughes, T. B.; Flynn, N.; Dang, N. L.; Swamidass, S. J.
Modeling the Bioactivation and Subsequent Reactivity of Drugs.
*Chem. Res. Toxicol.* **2021**, *34* (2), 584–600.
**DOI:** [10.1021/acs.chemrestox.0c00417](https://doi.org/10.1021/acs.chemrestox.0c00417)

```bibtex
@article{Hughes2021Bioactivation,
  title = {Modeling the Bioactivation and Subsequent Reactivity of Drugs},
  author = {Hughes, Tyler B. and Flynn, Noah and Dang, Na Le and Swamidass, S. Joshua},
  journal = {Chemical Research in Toxicology},
  volume = {34},
  number = {2},
  pages = {584--600},
  year = {2021},
  doi = {10.1021/acs.chemrestox.0c00417},
  url = {https://doi.org/10.1021/acs.chemrestox.0c00417}
}
```

## Other rules that match a XenoSite paper

| Rule / ruleset | Paper | DOI |
| --- | --- | --- |
| Epoxidation (`SO`) | Hughes, Miller, Swamidass. Modeling Epoxidation of Drug-like Molecules with a Deep Machine Learning Network. *ACS Cent. Sci.* **2015**, *1*, 168–180. | [10.1021/acscentsci.5b00131](https://doi.org/10.1021/acscentsci.5b00131) |
| Dealkylation (`UO`) | Dang, Hughes, Miller, Swamidass. Computationally Assessing the Bioactivation of Drugs by N-Dealkylation. *Chem. Res. Toxicol.* **2018**, *31*, 68–80. | [10.1021/acs.chemrestox.7b00191](https://doi.org/10.1021/acs.chemrestox.7b00191) |
| Glucuronidation (`CJ`) | Dang, Hughes, Krishnamurthy, Swamidass. A Simple Model Predicts UGT-Mediated Metabolism. *Bioinformatics* **2016**, *32*, 3183–3189. | [10.1093/bioinformatics/btw350](https://doi.org/10.1093/bioinformatics/btw350) |

Epoxidation on the web: [xenosite.org/epoxidation](https://xenosite.org/epoxidation). N-dealkylation: [xenosite.org/ndealk](https://xenosite.org/ndealk). UGT: [xenosite.org/ugt](https://xenosite.org/ugt).

## Metabolic Forest (this package)

The structure-enumeration system itself:

Hughes, T. B.; Dang, N. L.; Kumar, A.; Flynn, N. R.; Swamidass, S. J.
Metabolic Forest: Predicting the Diverse Structures of Drug Metabolites.
*J. Chem. Inf. Model.* **2020**, *60* (10), 4702–4716.
**DOI:** [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360)

BibTeX is in the [README](../README.md#citation).

Sequential metabolite networks (optional `network` extra) are related to XenoNet: Flynn, Dang, Ward, Swamidass. *J. Chem. Inf. Model.* **2020**, *60*, 3431–3449. **DOI:** [10.1021/acs.jcim.0c00361](https://doi.org/10.1021/acs.jcim.0c00361).

## Example

```python
from xenosite.metabolite import bfs, load_ruleset

# Rainbow Phase I
load_ruleset("PhaseOneRS")

# Quinone structures only
next(bfs(["Oc1ccc(O)cc1"], ruleset="QuinoneFormationRS"))

# Bioactivation routes
load_ruleset("Bioactivation")
```
