# Using xenosite.metabolite

This library enumerates metabolite *structures* with the Metabolic Forest reaction rules. For site-of-metabolism scores, use [xenosite.org](https://xenosite.org).

Please cite Hughes et al., *Metabolic Forest*, *J. Chem. Inf. Model.* 2020, DOI [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360). BibTeX is in the [README](../README.md#citation).

## Public API

```python
from xenosite.metabolite import bfs, rules, RuleSet, PhaseOneRS, load_ruleset, RULESETS
```

| Symbol | Role |
| --- | --- |
| `bfs(mols, ruleset=..., **kwargs)` | Search pathways from a reactant, or between reactant and product |
| `rules` | Individual reaction classes (`Hydroxylation`, `Epoxidation`, `QuinoneFormation`, …) |
| `RuleSet` | Group of rules used together |
| `load_ruleset(name)` | Look up a named ruleset (`"Full"`, `"PhaseOneRS"`, `"QuinoneFormationRS"`, …) |
| `PhaseOneRS` | Phase I rules used in Metabolic Forest |
| `RULESETS` | Registry of built-in rulesets |

`xenosite.metabolite.net.MetaboliteNetwork` is optional and needs `pip install 'xenosite-metabolite[network]'`.

## Enumerate metabolites

```python
from rdkit import Chem
from xenosite.metabolite import rules

mol = Chem.MolFromSmiles("c1ccccc1O")
for site, products in rules.QuinoneFormation().metabolites(mol):
    print(site, [Chem.MolToSmiles(p) for p in products])
```

Each `site` is `(rule_name, atom_or_bond_indices)`. Products are RDKit molecules.

## Search a pathway

```python
from xenosite.metabolite import bfs

smiles, steps, mols = next(
    bfs(["CCO", "C=CO"], ruleset="PhaseOneRS", depth=1, phase1=True)
)
print(smiles)  # product SMILES
print(steps)   # [(rule_name, site_strings), ...]
```

With one molecule, `bfs` enumerates metabolites of that reactant. With two, it searches for a path from the first to the second. `depth` is the maximum number of sequential reactions.

`phase1=True` formats sites as Phase I strings (for example `1.h` or `2.3`) instead of frozensets of atom indices.

## Custom rulesets

```python
from xenosite.metabolite import rules, RuleSet

rs = RuleSet([rules.Epoxidation(), rules.EpoxideOpening()], name="epoxide")
path = next(rs.find_path(reactant, product, depth=2))
```

Built-in names include `Full`, `PhaseOneRS`, `Bioactivation`, `QuinoneFormationRS`, and the Phase I groups `SO`, `DH`, `HD`, `RD`, and `UO`. Pass any of those strings to `bfs(..., ruleset=...)` or `load_ruleset(...)`.

## Command line

```bash
xenosite-metabolite CCO CC=O --ruleset PhaseOneRS --depth 1 --phase1
```

See `xenosite-metabolite --help`.
