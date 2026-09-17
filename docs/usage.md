# Using xenosite.forest

This library enumerates metabolite *structures* with the Metabolic Forest reaction rules. For site-of-metabolism scores, use [xenosite.org](https://xenosite.org).

Please cite Hughes et al., *Metabolic Forest*, *J. Chem. Inf. Model.* 2020, DOI [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360). BibTeX is in the [README](../README.md#citation).

## Public API

```python
from xenosite.forest import bfs, rules, RuleSet, PhaseOneRS, load_ruleset, RULESETS, AtomTrace, Step, StepPlan
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
| `Step` / `StepPlan` | Named reaction at a site; partial order with lazy linearizations |

`xenosite.forest.net.MetaboliteNetwork` is optional and needs `pip install 'xenosite-forest[network]'`.

## Phase1-equivalent steps

Some Forest rules correspond to one or more Phase I transformations. Ask a rule for plans with `phase1_steps(mol, site)` (0-based RDKit atom indices). Each plan is a `StepPlan` partial order; `iter_linearizations()` yields every consistent total order.

| Rule family | `phase1_steps` |
| --- | --- |
| Phase I (`Hydroxylation`, `Epoxidation`, `Dehydrogenation`, …) and `NDealkylation` | Degenerate: one singleton plan for a valid site |
| `QuinoneFormation` | Multi-step: prep layers (e.g. hydroxylation) then final dehydrogenation |
| Other rules (conjugates, …) | `NotImplementedError` |

```python
from rdkit import Chem
from xenosite.forest import StepPlan, rules

mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
plans = rules.QuinoneFormation().phase1_steps(mol, frozenset({4, 7}))
for plan in plans:
    for order in plan.iter_linearizations():
        print([(s.rule, sorted(s.site)) for s in order])

# Opt-in stamp on products (mol prop "phase1_steps")
_, products = next(
    rules.Epoxidation().metabolize(
        Chem.MolFromSmiles("C=C"), attach_phase1_steps=True, tag_atoms=False
    )
)
StepPlan.from_mol(products[0])
```

Sites in `Step` / `StepPlan` use **0-based** RDKit indices (same as reaction `site` frozensets), not 1-based atom numbers.

## Enumerate metabolites

```python
from rdkit import Chem
from xenosite.forest import rules

mol = Chem.MolFromSmiles("c1ccccc1O")
for site, products in rules.QuinoneFormation().metabolites(mol):
    print(site, [Chem.MolToSmiles(p) for p in products])
```

Each `site` is `(rule_name, atom_or_bond_indices)`. Products are RDKit molecules.

## Conjugation (Phase II)

`Acetylation`, `Glucuronidation`, `Glutathionation`, and `Sulfation` inherit `ConjugationRule`. **By default** products are bare `*` adducts (the conjugate group is collapsed). Opt in to the full chemical adduct or a CXSMILES label:

```python
from rdkit import Chem
from xenosite.forest import rules, load_ruleset
from xenosite.forest.utils import mol_to_cxsmiles

mol = Chem.MolFromSmiles("c1ccccc1O")

# Default: bare star
_, products = next(rules.Glucuronidation().metabolites(mol))
Chem.MolToSmiles(products[0])  # '*Oc1ccccc1'

# CXSMILES label on the dummy
_, products = next(rules.Glucuronidation(star_label="GlcA").metabolites(mol))
mol_to_cxsmiles(products[0])  # '*Oc1ccccc1 |$GlcA;;;;;;;$|'

# Full glucuronide
_, products = next(rules.Glucuronidation(as_star=False).metabolites(mol))
```

| Option | Meaning |
| --- | --- |
| `as_star=True` (default) | Collapse the conjugate group to a dummy `*` |
| `star_label=None` (default) | Bare `*` / plain SMILES |
| `star_label="GlcA"\|"GSH"\|…` | Set CX `atomLabel` on dummies; prefer `mol_to_cxsmiles` |
| `as_star=False` | Emit the full chemical adduct (GlcA, GSH, acetyl, sulfate) |

`Glutathionation` also takes `include_thiol=True` (default). Set `include_thiol=False` (or `load_ruleset("GlutathionationNoThiol")`) to drop the substrate-thiol disulfide SMARTS for DNA / cyanide-style reactivity.

| Label | Typical setup | Full structure? |
| --- | --- | --- |
| (none) | `Glutathionation()` | yes, with `as_star=False` |
| `GSH` | `Glutathionation(star_label="GSH")` | yes, with `as_star=False` |
| `Protein` | `Glutathionation(star_label="Protein")` | **no** (star-only) |
| `DNA` | `Glutathionation(include_thiol=False, star_label="DNA")` | **no** |
| `Cyanide` | `Glutathionation(include_thiol=False, star_label="Cyanide")` | **no** |
| `GlcA` | `Glucuronidation(star_label="GlcA")` | yes, with `as_star=False` |

The SMILES token before `|` is always valid on its own and depicts a bare `*`. See [rulesets](rulesets.md#conjugation) and [xenosite-predict notes](xenosite-predict.md).

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

With one molecule, `bfs` enumerates metabolites of that reactant. With two, it searches for a path from the first to the second. `depth` is the maximum number of sequential reactions. `dfs(...)` is the same API with depth-first order (better for sampling a few deep pathways). CLI: `--search dfs`.

Optional shared knobs: `max_paths=N` stops after N yields (default unlimited); `shuffle_rng=random.Random(seed)` randomizes frontier / reaction / product order on both BFS and DFS.

`phase1=True` formats sites as Phase I strings (for example `1.h` or `2.3`) instead of frozensets of atom indices.

Star conjugate adducts (`*` dummies from Phase II) are **not** metabolized further by default. Pass `expand_star_conjugates=True` (CLI: `--expand-star-conjugates`) to allow depth>1 expansion of those products.

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
