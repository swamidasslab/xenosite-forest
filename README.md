# xenosite.metabolite

Enumerate metabolite structures from reaction rules — the Metabolic Forest.

This package generates explicit metabolite structures (Phase I and common conjugations) from a reactant, and can search for pathways that connect a reactant to a putative product.

Site-of-metabolism *prediction* models live on the web at **[xenosite.org](https://xenosite.org)**. This library is the structure-enumeration engine behind Metabolic Forest, not those neural-network models.

## Install

```bash
uv add xenosite-metabolite
# or
pip install xenosite-metabolite
```

Optional NetworkX helpers:

```bash
uv add "xenosite-metabolite[network]"
```

Requires Python 3.10+ and RDKit.

## Quick start

```python
from rdkit import Chem
from xenosite.metabolite import bfs, rules

# Find a pathway from ethanol to acetaldehyde
smiles, steps, mols = next(bfs(["CCO", "CC=O"], ruleset="PhaseOneRS"))
print(smiles, steps)

# Enumerate hydroxylation metabolites of propane
for site, products in rules.Hydroxylation().metabolites(Chem.MolFromSmiles("CCC")):
    print(site, [Chem.MolToSmiles(p) for p in products])
```

A longer walkthrough is in [`examples/tutorial.ipynb`](examples/tutorial.ipynb).

Command line:

```bash
xenosite-metabolite CCO CC=O --ruleset PhaseOneRS --depth 1
```

## Citation

If you use this software, please cite the Metabolic Forest paper:

Hughes, T. B.; Dang, N. L.; Kumar, A.; Flynn, N. R.; Swamidass, S. J.
**Metabolic Forest: Predicting the Diverse Structures of Drug Metabolites.**
*J. Chem. Inf. Model.* **2020**, *60* (10), 4702–4716.
https://doi.org/10.1021/acs.jcim.0c00360

A machine-readable citation is in [`CITATION.cff`](CITATION.cff).

## Related

- [XenoSite](https://xenosite.org) — site-of-metabolism and reactivity models
- [Metabolic Forest (paper)](https://doi.org/10.1021/acs.jcim.0c00360)

## Development

```bash
git clone https://github.com/swamidasslab/xenosite-metabolite.git
cd xenosite-metabolite
uv sync --extra network --group dev
uv run pytest -n auto
```

The import path is `xenosite.metabolite`. `xenosite` is a PEP 420 namespace so other `xenosite.*` packages can be installed alongside this one.

## License

MIT. See [LICENSE](LICENSE).
