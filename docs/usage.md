# Using xenosite.forest

This library enumerates metabolite *structures* with the Metabolic Forest reaction rules. For site-of-metabolism scores, use [xenosite.org](https://xenosite.org).

Please cite Hughes et al., *Metabolic Forest*, *J. Chem. Inf. Model.* 2020, DOI [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360). BibTeX is in the [README](../README.md#citation).

## Public API

```python
from xenosite.forest import bfs, rules, RuleSet, PhaseOneRS, load_ruleset, RULESETS, AtomTrace
```

| Symbol | Role |
| --- | --- |
| `bfs(mols, ruleset=..., **kwargs)` | Search pathways from a reactant, or between reactant and product |
| `rules` | Individual reaction classes (`Hydroxylation`, `Epoxidation`, `QuinoneFormation`, …) |
| `RuleSet` | Group of rules used together |
| `load_ruleset(name)` | Look up a named ruleset (`"Full"`, `"PhaseOneRS"`, `"QuinoneFormationRS"`, …) |
| `PhaseOneRS` | Phase I rules used in Metabolic Forest |
| `RULESETS` | Registry of built-in rulesets |
| `AtomTrace(mol)` | 1-based atom-mapping history on a tagged metabolite |

`xenosite.forest.net.MetaboliteNetwork` is optional and needs `pip install 'xenosite-forest[network]'`. In a notebook, `net`, `net.draw()`, `net.grid()`, and `net.draw_path(path)` render the graph with structure drawings and Metabolic Rainbow edge colors. Defaults are publication-friendly (white figure background). Useful `draw(...)` options: `nodes=[...]`, `reaction_types`, `max_generation`, `background`, `mol_background`, `mol_border`; `prune(...)` returns a filtered copy; `Drawing.save("figure.pdf")` / `net.save_draw(...)` write SVG (always) or PDF/PS (via `rsvg-convert` or ImageMagick).

## Enumerate metabolites

```python
from rdkit import Chem
from xenosite.forest import rules

mol = Chem.MolFromSmiles("c1ccccc1O")
for site, products in rules.QuinoneFormation().metabolites(mol):
    print(site, [Chem.MolToSmiles(p) for p in products])
```

Each `site` is `(rule_name, atom_or_bond_indices)`. Products are RDKit molecules.

## Indexing

Two scales. Mixing them is the confusing part.

| Scale | What | Values |
| --- | --- | --- |
| **Atom number** (public) | `AtomTrace`, SMILES `:N`, `GetAtomMapNum()`, Phase I `1.h` | **1-based.** `0` on a map number means unmapped / new, not atom zero. |
| **Depth / step** | `t.depths`, `map(start_depth=0)`, reaction count | **0-based.** Depth 0 is the original reactant. These are not atom numbers. |
| **RDKit index** (internal only) | `GetIdx()`, `react_atom_idx`, `current_idx`, `ATOM_INDEX_PATHS` | **0-based.** Never export. |

Convert with `xenosite.forest.trace.atom_no` (GetIdx → atom number) and `rdkit_idx` (atom number → GetIdx).

## Atom tracing

After a reaction (unless `do_not_tag_atoms=True`), each product carries origin map numbers and an `AtomTrace`. Atom numbers are 1-based; depths are 0-based.

```python
from rdkit import Chem
from xenosite.forest import AtomTrace, rules

mol = Chem.MolFromSmiles("CCO")
_, products = next(rules.Hydroxylation().metabolize(mol))
p = products[0]
t = AtomTrace(p)

Chem.MolToSmiles(p, canonical=False)
# e.g. [CH2:1](O)[CH2:2][OH:3]  — :N is origin at depth 0; the new O has no map

t.map()       # {1: 1, 2: 3, 3: 4}  depth-0 atom number -> current atom number
t.follow(1)   # (1, 1)
t.origin(2)   # None if current atom 2 is the new O
t.added()     # frozenset of new 1-based atom numbers
t.removed()   # frozenset of depth-0 atom numbers that disappeared
t.depths      # (0, 1)
```

Canonical `MolToSmiles` may scramble atom order; reactant-aligned order uses `canonical=False`. Map numbers still appear in canonical SMILES, so `:N` alignment works either way.

`do_not_tag_atoms=True` skips tags, maps, and reorder.

## Search a pathway

```python
from xenosite.forest import bfs

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
from xenosite.forest import rules, RuleSet

rs = RuleSet([rules.Epoxidation(), rules.EpoxideOpening()], name="epoxide")
path = next(rs.find_path(reactant, product, depth=2))
```

Built-in names include `Full`, `PhaseOneRS`, `Bioactivation`, `QuinoneFormationRS`, and the Phase I groups `SO`, `DH`, `HD`, `RD`, and `UO`. Pass any of those strings to `bfs(..., ruleset=...)` or `load_ruleset(...)`.

Which ruleset matches which paper (Rainbow, quinone, bioactivation, and others) is in **[rulesets.md](rulesets.md)**.

## Command line

```bash
xenosite-forest CCO CC=O --ruleset PhaseOneRS --depth 1 --phase1
```

See `xenosite-forest --help`.
