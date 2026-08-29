# xenosite.metabolite

Python implementation of **Metabolic Forest**: enumerate explicit metabolite structures from reaction rules, and search pathways that connect a reactant to a putative product.

Site-of-metabolism *prediction* models (epoxidation, quinonation, Phase I, reactivity, and others) are on the web at **[xenosite.org](https://xenosite.org)**. This package is the structure-enumeration engine, not those neural-network models.

If you use this software, please cite the Metabolic Forest paper (DOI and BibTeX below).

## Install

```bash
uv add xenosite-metabolite
```

or

```bash
pip install xenosite-metabolite
```

Optional NetworkX helpers for building a metabolite graph:

```bash
uv add "xenosite-metabolite[network]"
```

Requires **Python 3.10+** and RDKit.

## Quick start

```python
from rdkit import Chem
from xenosite.metabolite import bfs, rules, PhaseOneRS

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

A longer walkthrough is in [`examples/tutorial.ipynb`](examples/tutorial.ipynb). API notes are in [`docs/usage.md`](docs/usage.md).

### Command line

```bash
xenosite-metabolite CCO CC=O --ruleset PhaseOneRS --depth 1
```

Useful flags: `--all-paths`, `--depth N`, `--phase1` (Phase I site strings), `--max N`, `--ruleset NAME`.

## Documentation

- **[`docs/usage.md`](docs/usage.md)** — public API, rulesets, and pathway search
- **[`examples/tutorial.ipynb`](examples/tutorial.ipynb)** — interactive tutorial
- **[xenosite.org](https://xenosite.org)** — XenoSite models for sites of metabolism and reactivity

Import the package as `xenosite.metabolite`. `xenosite` is a PEP 420 namespace, so other `xenosite.*` packages can be installed alongside this one.

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

A machine-readable citation is also in [`CITATION.cff`](CITATION.cff).

## Development

```bash
git clone https://github.com/swamidasslab/xenosite-metabolite.git
cd xenosite-metabolite
uv sync --extra network --group dev
uv run pytest -n auto
```

## License

MIT. See [LICENSE](LICENSE).
