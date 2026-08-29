# Built-in rulesets

Pass any **name** below to `bfs(..., ruleset=...)` or `load_ruleset(...)`. Aliases in parentheses are also registered.

These rules *enumerate structures*. The cited papers are the XenoSite models and labeling schemes the rulesets are written to match. Site-of-metabolism scores from those models are at [xenosite.org](https://xenosite.org).

Always cite Metabolic Forest when you use this package. Cite the matching paper as well when a specific ruleset or labeling scheme is central to the work.

## How they fit together

Metabolic Forest groups 24 reaction rules into eight rulesets: the five Rainbow Phase I classes plus conjugation, quinone formation, and tautomerization. `Full` is that complete generator. `PhaseOneRS` is only the five Phase I classes (as implemented for structure search in this package).

| Ruleset | Name | Role |
| --- | --- | --- |
| `Full` | `Full` | All Metabolic Forest rules (Phase I + conjugation + quinone + tautomerization) |
| `PhaseOneRS` | `PhaseOneRS` | Phase I only: `SO` + `UO` + `DH` + `HD` + `RD` |
| `StableOxygenationRS` | `SO` (`StableOxygenation`) | Stable oxygenation |
| `UnstableOxygenationRS` | `UO` (`UnstableOxygenation`) | Unstable oxygenation |
| `DehydrogenationRS` | `DH` (`Dehydrogenation`) | Dehydrogenation |
| `HydrolysisRS` | `HD` (`Hydrolysis`) | Hydrolysis (Forest definition; see below) |
| `ReductionRS` | `RD` (`Reduction`) | Reduction (Forest definition; see below) |
| `QuinoneFormationRS` | `QF` (`QuinoneFormation`) | Quinone / quinone-imine / quinone-methide |
| `ConjugationRS` | `CJ` (`Conjugation`) | Phase II conjugations |
| `TautomerizationRS` | `TT` (`Tautomerization`) | Tautomerization |
| `Bioactivation` | `BA` (`BioactivationPathways`) | Four common bioactivation routes |
| `ThiopheneSulfurOxidationRS` | `TSO` (`ThiopheneSulfurOxidation`) | Thiophene S-oxidation alone |

`PhaseOneRS` hydrolysis/reduction lists are slightly narrower than the `HD` / `RD` groups inside `Full` (no azo splitting or benzodioxole reduction). See each class section for the exact rules and how they relate to the Rainbow paper’s five-color labeling.

## Metabolic Rainbow (five colors)

Dang, N. L.; Matlock, M. K.; Hughes, T. B.; Swamidass, S. J.
The Metabolic Rainbow: Deep Learning Phase I Metabolism in Five Colors.
*J. Chem. Inf. Model.* **2020**, *60* (3), 1146–1164.
**DOI:** [10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836)

The Rainbow paper labels Phase I sites of metabolism in five reaction classes. Each class has a mnemonic color. Those colors were mapped to colorblind-safe hex codes from Bang Wong’s Nature Methods palette (Okabe–Ito; Wong, *Nat. Methods* **2010**, *7*, 573).

| Abbr | Class | Color | Hex |
| --- | --- | --- | --- |
| **SO** | Stable Oxygenation | red (vermillion) | `#D55E00` |
| **UO** | Unstable Oxygenation | orange | `#E69F00` |
| **DH** | Dehydrogenation | green (bluish green) | `#009E73` |
| **HD** | Hydrolysis | blue (sky blue) | `#56B4E9` |
| **RD** | Reduction | purple (reddish purple) | `#CC79A7` |

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

Metabolic Forest (and this package) reuse the same five abbreviations for Phase I structure rules, but the **rules in each class are not identical** to Rainbow Table 1. Rainbow groups reaction *types for SOM labeling*; Metabolic Forest groups *structure-enumeration rules*. Differences are called out under each class.

### SO — Stable Oxygenation (red, `#D55E00`)

**Rainbow reaction types:** aromatic/aliphatic hydroxylation, aromatic/aliphatic epoxidation, N-oxidation, S-oxidation.

**Rules in this package** (`SO` / `StableOxygenationRS`, and in `PhaseOneRS` / `Full`):

- `Hydroxylation`
- `Epoxidation`
- `NitrogenOxidation`
- `SulfurOxidation`

This matches Rainbow SO and Metabolic Forest’s stable-oxygenation ruleset.

### UO — Unstable Oxygenation (orange, `#E69F00`)

**Rainbow reaction types:** N-/O-/S-/C-dealkylation, oxidative deamination, oxidative dehalogenation.

**Rules in this package** (`UO` / `UnstableOxygenationRS`, and in `PhaseOneRS` / `Full`):

- `Dealkylation` (covers N-/O-/S-/C-dealkylation and oxidative deamination patterns in the SMARTS)
- `OxidativeDehalogenation`

This matches Rainbow UO and Metabolic Forest’s unstable-oxygenation ruleset.

### DH — Dehydrogenation (green, `#009E73`)

**Rainbow reaction types:** dehydrogenation, alcohol → aldehyde/ketone, single → double bond, double → triple bond, **quinone formation**, iminium formation.

**Rules in this package** (`DH` / `DehydrogenationRS`, and in `PhaseOneRS` / `Full`):

- `Dehydrogenation` only

Quinone / imine / methide formation is **not** in `DH` here. It is a separate ruleset, `QuinoneFormationRS` (`QF`), because Metabolic Forest treats quinone structure generation as its own algorithm. Rainbow still counts quinone formation under the green DH label for SOM models.

### HD — Hydrolysis (blue, `#56B4E9`)

**Rainbow reaction types:** ester, amide, ether, and cyanide hydrolysis.

**Rules in this package:**

| Context | Rules |
| --- | --- |
| `PhaseOneRS` | `Hydrolysis`, `EpoxideOpening`, `Dephosphorylation` |
| `HD` / `HydrolysisRS` inside `Full` | those three, plus `AzoSplitting` |

Metabolic Forest places dephosphorylation, epoxide opening, carbonyl cleavage (`Hydrolysis`), and azo splitting in the hydrolysis ruleset. Rainbow’s blue HD label is narrower (hydrolysis of esters/amides/ethers/cyanide only). Use `PhaseOneRS` or `HD` for Forest-style enumeration; do not treat package `HD` as a 1:1 copy of Rainbow Table 1.

### RD — Reduction (purple, `#CC79A7`)

**Rainbow reaction types:** carbonyl, nitro, and sulfo reduction; reductive dehalogenation; hydrogenation.

**Rules in this package:**

| Context | Rules |
| --- | --- |
| `PhaseOneRS` | `Hydrogenation`, `Dehydration`, `NitrogenReduction`, `SulfurReduction`, `OxygenReduction`, `ReductiveDehalogenation` |
| `RD` / `ReductionRS` inside `Full` | those six, plus `BenzodioxoleReduction` |

Metabolic Forest’s reduction ruleset is broader than Rainbow’s purple RD label (it adds dehydration, oxygen reduction, and—in `Full`—benzodioxole reduction). Carbonyl reduction in Rainbow often corresponds to `Hydrogenation` / `OxygenReduction` patterns in the rules.

## Quinone formation

`QuinoneFormationRS` (`QF`) implements the specialized quinone / quinone-imine / quinone-methide structure algorithm from Metabolic Forest. It is matched to the quinone XenoSite model. In Rainbow labeling, quinone formation is a green (DH) reaction type; in this package it is its own ruleset, not part of `DH` or `PhaseOneRS`.

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

`Bioactivation` (`BA`) groups four common bioactivation routes modeled jointly in the 2021 paper:

- `QuinoneFormation`
- `Epoxidation`
- `NitroaromaticReduction`
- `ThiopheneSulfurOxidation`

This is not a Rainbow color and is not part of `Full`. It is a convenience ruleset for bioactivation pathway enumeration.

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

## Conjugation

`ConjugationRS` (`CJ`) is Phase II structure enumeration in Metabolic Forest:

- `Acetylation`
- `Glucuronidation`
- `Glutathionation`
- `Sulfation`

Rainbow Phase I does not include conjugations.

## Tautomerization

`TautomerizationRS` (`TT`) has a single `Tautomerization` rule. Metabolic Forest includes it because tautomers appear in known pathways even though tautomerization is not usually labeled as metabolism. Rainbow excluded tautomerization from the five colors (it was part of the ~7.7% of Phase I records left out of the Rainbow dataset).

```python
from xenosite.metabolite import rules, RuleSet

tt = RuleSet([rules.Tautomerization()], name="tautomerization")
```

## Thiophene sulfur oxidation

Rule: `ThiopheneSulfurOxidation`. Also registered as `TSO` / `ThiopheneSulfurOxidationRS`, and included in `Bioactivation`. Matched to the bioactivation paper (above).

```python
from xenosite.metabolite import rules, RuleSet, load_ruleset

tso = RuleSet([rules.ThiopheneSulfurOxidation()], name="thiophene_s_oxidation")
# or: load_ruleset("TSO")
```

## Epoxidation

Rule: `Epoxidation`. Matched to the epoxidation XenoSite model.

Hughes, T. B.; Miller, G. P.; Swamidass, S. J.
Modeling Epoxidation of Drug-like Molecules with a Deep Machine Learning Network.
*ACS Cent. Sci.* **2015**, *1* (4), 168–180.
**DOI:** [10.1021/acscentsci.5b00131](https://doi.org/10.1021/acscentsci.5b00131)

```python
from xenosite.metabolite import rules, RuleSet

epoxidation = RuleSet([rules.Epoxidation()], name="epoxidation")
```

```bibtex
@article{Hughes2015Epoxidation,
  title = {Modeling Epoxidation of Drug-like Molecules with a Deep Machine Learning Network},
  author = {Hughes, Tyler B. and Miller, Grover P. and Swamidass, S. Joshua},
  journal = {ACS Central Science},
  volume = {1},
  number = {4},
  pages = {168--180},
  year = {2015},
  doi = {10.1021/acscentsci.5b00131},
  url = {https://doi.org/10.1021/acscentsci.5b00131}
}
```

On the web: [xenosite.org/epoxidation](https://xenosite.org/epoxidation).

## N-dealkylation

Rule: `Dealkylation` (SMARTS also cover O-/S-/C-dealkylation and oxidative deamination patterns). Matched to the N-dealkylation XenoSite model.

Dang, N. L.; Hughes, T. B.; Miller, G. P.; Swamidass, S. J.
Computationally Assessing the Bioactivation of Drugs by N-Dealkylation.
*Chem. Res. Toxicol.* **2018**, *31* (1), 68–80.
**DOI:** [10.1021/acs.chemrestox.7b00191](https://doi.org/10.1021/acs.chemrestox.7b00191)

```python
from xenosite.metabolite import rules, RuleSet

dealkylation = RuleSet([rules.Dealkylation()], name="dealkylation")
```

```bibtex
@article{Dang2018NDealkylation,
  title = {Computationally Assessing the Bioactivation of Drugs by N-Dealkylation},
  author = {Dang, Na Le and Hughes, Tyler B. and Miller, Grover P. and Swamidass, S. Joshua},
  journal = {Chemical Research in Toxicology},
  volume = {31},
  number = {1},
  pages = {68--80},
  year = {2018},
  doi = {10.1021/acs.chemrestox.7b00191},
  url = {https://doi.org/10.1021/acs.chemrestox.7b00191}
}
```

On the web: [xenosite.org/ndealk](https://xenosite.org/ndealk).

## UGT (glucuronidation)

Rule: `Glucuronidation`. Matched to the UGT XenoSite model.

Dang, N. L.; Hughes, T. B.; Krishnamurthy, V.; Swamidass, S. J.
A Simple Model Predicts UGT-Mediated Metabolism.
*Bioinformatics* **2016**, *32* (20), 3183–3189.
**DOI:** [10.1093/bioinformatics/btw350](https://doi.org/10.1093/bioinformatics/btw350)

```python
from xenosite.metabolite import rules, RuleSet

ugt = RuleSet([rules.Glucuronidation()], name="ugt")
```

```bibtex
@article{Dang2016UGT,
  title = {A Simple Model Predicts UGT-Mediated Metabolism},
  author = {Dang, Na Le and Hughes, Tyler B. and Krishnamurthy, Varun and Swamidass, S. Joshua},
  journal = {Bioinformatics},
  volume = {32},
  number = {20},
  pages = {3183--3189},
  year = {2016},
  doi = {10.1093/bioinformatics/btw350},
  url = {https://doi.org/10.1093/bioinformatics/btw350}
}
```

On the web: [xenosite.org/ugt](https://xenosite.org/ugt).

## Glutathionation

Rule: `Glutathionation`. Enumerates glutathione conjugates (epoxide, halogen, and sulfur motifs). Matched to the glutathione site-of-reactivity XenoSite model.

Hughes, T. B.; Miller, G. P.; Swamidass, S. J.
Site of Reactivity Models Predict Molecular Reactivity of Diverse Chemicals with Glutathione.
*Chem. Res. Toxicol.* **2015**, *28* (4), 797–809.
**DOI:** [10.1021/acs.chemrestox.5b00017](https://doi.org/10.1021/acs.chemrestox.5b00017)

```python
from xenosite.metabolite import rules, RuleSet

gsh = RuleSet([rules.Glutathionation()], name="glutathionation")
```

```bibtex
@article{Hughes2015Glutathione,
  title = {Site of Reactivity Models Predict Molecular Reactivity of Diverse Chemicals with Glutathione},
  author = {Hughes, Tyler B. and Miller, Grover P. and Swamidass, S. Joshua},
  journal = {Chemical Research in Toxicology},
  volume = {28},
  number = {4},
  pages = {797--809},
  year = {2015},
  doi = {10.1021/acs.chemrestox.5b00017},
  url = {https://doi.org/10.1021/acs.chemrestox.5b00017}
}
```

The broader multitask reactivity model (GSH, cyanide, protein, DNA) is Hughes et al., *ACS Cent. Sci.* **2016**, *2*, 529–537. **DOI:** [10.1021/acscentsci.6b00162](https://doi.org/10.1021/acscentsci.6b00162). On the web: [xenosite.org](https://xenosite.org) (Reactivity).

## Metabolic Forest (this package)

The structure-enumeration system itself:

Hughes, T. B.; Dang, N. L.; Kumar, A.; Flynn, N. R.; Swamidass, S. J.
Metabolic Forest: Predicting the Diverse Structures of Drug Metabolites.
*J. Chem. Inf. Model.* **2020**, *60* (10), 4702–4716.
**DOI:** [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360)

BibTeX is in the [README](../README.md#citation).

Sequential metabolite networks (optional `network` extra) are related to XenoNet:

Flynn, N. R.; Dang, N. L.; Ward, M. D.; Swamidass, S. J.
XenoNet: Inference and Likelihood of Intermediate Metabolite Formation.
*J. Chem. Inf. Model.* **2020**, *60* (7), 3431–3449.
**DOI:** [10.1021/acs.jcim.0c00361](https://doi.org/10.1021/acs.jcim.0c00361)

## Example

```python
from xenosite.metabolite import bfs, load_ruleset

# Rainbow-aligned Phase I structure search
load_ruleset("PhaseOneRS")

# One color at a time
load_ruleset("SO")   # Stable Oxygenation
load_ruleset("UO")   # Unstable Oxygenation

# Quinone structures only
next(bfs(["Oc1ccc(O)cc1"], ruleset="QuinoneFormationRS"))

# Bioactivation routes
load_ruleset("Bioactivation")
```
