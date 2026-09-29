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

The default install is the **Rust** chematic door (no RDKit). For the RDKit
reference engine and the frozen 0.6.x archive:

```bash
uv add "xenosite-forest[rdkit]"
```

Optional NetworkX helpers for building a metabolite graph:

```bash
uv add "xenosite-forest[network]"
```

Optional [xenosite-predict](https://github.com/swamidasslab/xenosite-predict) (site-of-metabolism models; Python ≥3.11). Not installed in this repo’s CI:

```bash
uv add "xenosite-forest[predict]"
```

Requires **Python 3.11–3.14**. RDKit 2022.03+ only with the ``[rdkit]`` extra.

A WASM-clean Rust crate (chematic + canonaut) lives in [`crates/xenosite-forest`](crates/xenosite-forest). See [`docs/forest/RUST.md`](docs/forest/RUST.md) and [`docs/forest/MIGRATING_0.8.md`](docs/forest/MIGRATING_0.8.md).

## Quick start

```python
from xenosite.forest import find_path, PhaseOne, Epoxidation

hits, counters = find_path("CC", "CCO", max_paths=1)
print(hits[0]["steps"], counters["billed"])
```

RDKit reference engine (requires ``[rdkit]``):

```python
from rdkit import Chem
from xenosite.forest.native import bfs, rules, PhaseOneRS

smiles, steps, mols = next(bfs(["CCO", "CC=O"], ruleset="PhaseOneRS"))
print(smiles, steps)

for products, info in rules.Hydroxylation().metabolize(Chem.MolFromSmiles("CCC")):
    print(info["site"], [p.xf.csmi for p in products])
```

A longer walkthrough is in [`examples/tutorial.ipynb`](https://github.com/swamidasslab/xenosite-forest/blob/main/examples/tutorial.ipynb). API notes are in [`docs/usage.md`](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/usage.md). Rulesets and their papers are in [`docs/rulesets.md`](https://github.com/swamidasslab/xenosite-forest/blob/main/docs/rulesets.md).

### Command line

Requires the ``[rdkit]`` extra (temporary native bridge until a Rust CLI ships):

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
- **[Migrating to 0.7.0](docs/forest/MIGRATING_0.7.md)** — `metabolize` yield shape and other breaking changes
- **[Migrating to 0.8.0](docs/forest/MIGRATING_0.8.md)** — rust / native / legacy layout, optional RDKit
- **[Legacy freeze](docs/forest/LEGACY.md)** — frozen 0.6.x archive notes
- **[xenosite-predict notes](docs/xenosite-predict.md)** — RDKit 2026 valence caches and DNA/CN conjugation vs GSH (downstream adapters)
- **[Path-search performance](docs/forest/PERFORMANCE.md)** — archive BFS/DFS vs live `find_path` head-to-head
- **[Heuristics / divergences / dropped](docs/forest/)** — design notes (`HEURISTICS`, `DIVERGENCES`, `DROPPED`, `PAIR_ORBITS`)
- **[Tutorial notebook](https://github.com/swamidasslab/xenosite-forest/blob/main/examples/tutorial.ipynb)** — interactive walkthrough
- **[xenosite.org](https://xenosite.org)** — XenoSite models for sites of metabolism and reactivity
- **[Source repository](https://github.com/swamidasslab/xenosite-forest)** — code, issues, and releases

Import the package as `xenosite.forest` (Rust public stub). The RDKit reference
engine is `xenosite.forest.native`; the frozen 0.6.x archive is
`xenosite.forest.legacy` (see [`docs/forest/LEGACY.md`](docs/forest/LEGACY.md)).
`xenosite` is a PEP 420 namespace, so other `xenosite.*` packages can be
installed alongside this one.

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
uv sync --extra network --group dev   # maturin builds the Rust extension
uv run pytest -n auto
uv run pytest --cov=xenosite.forest
```

Versioning is the ``package.version`` in [`crates/xenosite-forest/Cargo.toml`](crates/xenosite-forest/Cargo.toml) (maturin). Tag a release as `vX.Y.Z` matching that Cargo version. Read it at runtime as `xenosite.forest.__version__`. User-facing notes go in [`changelog.d/`](changelog.d/); the tag workflow compiles them into `CHANGELOG.md`. See [docs/release.md](docs/release.md).

Pushing a `v*` tag runs [`.github/workflows/release.yml`](https://github.com/swamidasslab/xenosite-forest/blob/main/.github/workflows/release.yml): tests must pass, then platform wheels (manylinux / musllinux / macOS / Windows, x86_64 + aarch64) and an sdist are published to GitHub Releases and PyPI. A red tag workflow means do not treat that tag as released.

## License

MIT. See [LICENSE](https://github.com/swamidasslab/xenosite-forest/blob/main/LICENSE).
