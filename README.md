# xenosite.forest

Python implementation of **Metabolic Forest**: enumerate explicit metabolite structures from reaction rules, and search pathways that connect a reactant to a putative product.

Site-of-metabolism *prediction* models (epoxidation, quinonation, Phase I, reactivity, and others) are on the web at **[xenosite.org](https://xenosite.org)**. This package is the structure-enumeration engine, not those neural-network models.

If you use this software, please cite the Metabolic Forest paper (DOI and BibTeX below).

## Install

```bash
uv add xenosite-forest
```

or

```bash
pip install xenosite-forest
```

Optional NetworkX helpers for building a metabolite graph:

```bash
uv add "xenosite-forest[network]"
```

Requires **Python 3.10+** and RDKit.

## Quick start

```python
from rdkit import Chem
from xenosite.forest import bfs, rules, PhaseOneRS

# Pathway from ethanol to acetaldehyde
smiles, steps, mols = next(bfs(["CCO", "CC=O"], ruleset="PhaseOneRS"))
print(smiles)
print(steps)

# Enumerate hydroxylation products of propane
for site, products in rules.Hydroxylation().metabolites(Chem.MolFromSmiles("CCC")):
    print(site, [Chem.MolToSmiles(p) for p in products])

# Named Phase I ruleset
print(sorted({rule.name for rule in PhaseOneRS}))
```

A longer walkthrough is in [`examples/tutorial.ipynb`](https://github.com/swamidasslab/xenosite-forest/blob/main/examples/tutorial.ipynb). API notes are in [`docs/usage.md`](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/usage.md). Rulesets and their papers are in [`docs/rulesets.md`](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/rulesets.md).

### Command line

```bash
xenosite-forest CCO CC=O --ruleset PhaseOneRS --depth 1
```

Useful flags: `--all-paths`, `--depth N`, `--phase1` (Phase I site strings), `--max N`, `--ruleset NAME`.

## Rulesets and papers

Pass these names to `bfs(..., ruleset=...)` or `load_ruleset(...)`. Full rule lists, Rainbow colors/hex codes, aliases, and BibTeX are in **[`docs/rulesets.md`](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/rulesets.md)**.

| Ruleset | What it enumerates | Matched paper |
| --- | --- | --- |
| `PhaseOneRS` | Phase I: SO, UO, DH, HD, RD | Dang et al., Metabolic Rainbow, *JCIM* 2020. DOI [10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836) |
| `SO` / `UO` / `DH` / `HD` / `RD` | One Rainbow color each (see hex table in docs) | Same Rainbow paper; Forest rule lists differ slightly for HD/RD/DH |
| `QuinoneFormationRS` (`QF`) | Quinone, quinone-imine, and quinone-methide structures | Hughes & Swamidass, *Chem. Res. Toxicol.* 2017. DOI [10.1021/acs.chemrestox.6b00385](https://doi.org/10.1021/acs.chemrestox.6b00385) |
| `Bioactivation` (`BA`) | Quinone, epoxidation, nitroaromatic reduction, thiophene S-oxidation | Hughes et al., *Chem. Res. Toxicol.* 2021. DOI [10.1021/acs.chemrestox.0c00417](https://doi.org/10.1021/acs.chemrestox.0c00417) |
| `Full` | Complete Metabolic Forest generator (Phase I, conjugations, quinone, tautomerization) | Hughes et al., Metabolic Forest, *JCIM* 2020. DOI [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360) |

Rainbow colorblind-safe hex (Wong / Okabe–Ito, closest to the paper figures): **SO** Stable Oxygenation red `#D55E00`, **UO** Unstable Oxygenation orange `#E69F00`, **DH** Dehydrogenation green `#009E73`, **HD** Hydrolysis blue `#56B4E9`, **RD** Reduction purple `#CC79A7`.

Related single-rule papers: epoxidation ([10.1021/acscentsci.5b00131](https://doi.org/10.1021/acscentsci.5b00131)), N-dealkylation ([10.1021/acs.chemrestox.7b00191](https://doi.org/10.1021/acs.chemrestox.7b00191)), UGT glucuronidation ([10.1093/bioinformatics/btw350](https://doi.org/10.1093/bioinformatics/btw350)), glutathione reactivity ([10.1021/acs.chemrestox.5b00017](https://doi.org/10.1021/acs.chemrestox.5b00017)).

## Documentation

- **[Rulesets and papers](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/rulesets.md)** — every built-in ruleset, Rainbow colors, and publication BibTeX
- **[Usage](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/usage.md)** — public API and pathway search
- **[Tutorial notebook](https://github.com/swamidasslab/xenosite-forest/blob/main/examples/tutorial.ipynb)** — interactive walkthrough
- **[xenosite.org](https://xenosite.org)** — XenoSite models for sites of metabolism and reactivity
- **[Source repository](https://github.com/swamidasslab/xenosite-forest)** — code, issues, and releases

Import the package as `xenosite.forest`. `xenosite` is a PEP 420 namespace, so other `xenosite.*` packages can be installed alongside this one.

## Citation

Please cite Metabolic Forest if you use this package:

Hughes, T. B.; Dang, N. L.; Kumar, A.; Flynn, N. R.; Swamidass, S. J.
Metabolic Forest: Predicting the Diverse Structures of Drug Metabolites.
*J. Chem. Inf. Model.* **2020**, *60* (10), 4702–4716.
**DOI:** [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360)

BibTeX (copy and paste):

```bibtex
@article{Hughes2020MetabolicForest,
  title = {Metabolic Forest: Predicting the Diverse Structures of Drug Metabolites},
  author = {Hughes, Tyler B. and Dang, Na Le and Kumar, Ayush and Flynn, Noah R. and Swamidass, S. Joshua},
  journal = {Journal of Chemical Information and Modeling},
  volume = {60},
  number = {10},
  pages = {4702--4716},
  year = {2020},
  doi = {10.1021/acs.jcim.0c00360},
  url = {https://doi.org/10.1021/acs.jcim.0c00360},
  publisher = {American Chemical Society}
}
```

A machine-readable citation is also in [`CITATION.cff`](https://github.com/swamidasslab/xenosite-forest/blob/main/CITATION.cff).

## Development

```bash
git clone https://github.com/swamidasslab/xenosite-forest.git
cd xenosite-forest
uv sync --extra network --group dev
uv run pytest -n auto
```

Versioning comes from git tags via [hatch-vcs](https://github.com/ofek/hatch-vcs) (setuptools-scm). Tag a release as `vX.Y.Z` (for example `v0.1.0`). On that commit the version is `X.Y.Z`. On later untagged commits it becomes the next patch with a dev suffix and short commit, for example `0.1.1.dev3+gabc1234`. Read it at runtime as `xenosite.forest.__version__`.

Pushing a `v*` tag runs [`.github/workflows/release.yml`](https://github.com/swamidasslab/xenosite-forest/blob/main/.github/workflows/release.yml): tests must pass and the resolved version must be a clean `X.Y.Z` before a GitHub Release and PyPI upload. A red tag workflow means do not treat that tag as released. To *block* creating tags unless checks pass, add a GitHub Ruleset on `refs/tags/v*` that requires the `release` / `test` status checks.

### One-time PyPI Trusted Publishing setup

No API token is stored in the repo. CI authenticates with [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) (OIDC).

1. Create the GitHub repo `swamidasslab/xenosite-forest` and push `main`.
2. GitHub → **Settings → Environments** → create an environment named exactly `pypi` (optional: require reviewers before deploy).
3. On [PyPI](https://pypi.org/manage/account/publishing/): add a **pending** trusted publisher (project does not need to exist yet):

   | Field | Value |
   | --- | --- |
   | PyPI project name | `xenosite-forest` |
   | Owner | `swamidasslab` |
   | Repository name | `xenosite-forest` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

4. Tag a release (`git tag -a v0.1.0 -m "0.1.0" && git push origin v0.1.0`). The `pypi-publish` job uploads only after tests and build succeed.

For a dry run, register the same publisher on [TestPyPI](https://test.pypi.org/manage/account/publishing/) first and temporarily point the publish action at TestPyPI.

## License

MIT. See [LICENSE](https://github.com/swamidasslab/xenosite-forest/blob/main/LICENSE).
