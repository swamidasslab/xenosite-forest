# Using xenosite.forest

This library enumerates metabolite *structures* with the Metabolic Forest reaction rules. For site-of-metabolism scores, use [xenosite.org](https://xenosite.org).

Please cite Hughes et al., *Metabolic Forest*, *J. Chem. Inf. Model.* 2020, DOI [10.1021/acs.jcim.0c00360](https://doi.org/10.1021/acs.jcim.0c00360). BibTeX is in the [README](../README.md#citation).

## Public API

```python
from xenosite.forest import (
    bfs, dfs, find_path, rules, RuleSet, PhaseOneRS, PhaseOneQF,
    load_ruleset, RULESETS, AtomTrace, AtomRef, Linearization, Step,
    StepPlan, And, Or, Deps, PathOutcome, PathSearchCounters, FormulaHint,
)
```

| Symbol | Role |
| --- | --- |
| `bfs` / `dfs` | Enumerate pathways from a reactant (or toward a product) |
| `find_path(R, T, …)` | MCS-guided search; yields `PathOutcome` (unstable API) |
| `PathOutcome` | Required `plan` + cleavage-side `maybe`; unpacks as `(smiles, steps, mols, plan, maybe)` |
| `PathSearchCounters` | Shared guided/classic instrumentation (`billed`, `sanitize_dropped`, …) |
| `rules` / `RuleSet` / `load_ruleset` / `RULESETS` | Reaction rules and named rulesets |
| `PhaseOneRS` / `PhaseOneQF` | Phase I; Phase I + `QuinoneFormation` (guided default) |
| `AtomTrace` | 1-based atom-mapping history on a tagged metabolite |
| `AtomRef` / `Step` / `StepPlan` / `And` / `Or` / `Deps` / `Linearization` | Pathway plans and self-apply linearizations |
| `FormulaHint` / `ADD_O` / `CLEAVE` / … | Per-SMARTS formula effects toward a target |

`xenosite.forest.net.MetaboliteNetwork` is optional (`pip install 'xenosite-forest[network]'`).
Site-of-metabolism models are optional via `pip install 'xenosite-forest[predict]'` (not used in Forest CI).

## Phase1-equivalent steps

Some Forest rules correspond to one or more Phase I transformations. Each rule
returns **one** `StepPlan` (nested `Step` / Seq / `And` / `Or`; guided emission
may use `Deps` for flat steps + precedes).

Sites use reactant-stable `AtomRef` values: a reactant **origin** index
(remapped through `atom_trace` on resolve), or an atom **created by** an
earlier step (`added_by=(rule, site)`). Both forms carry frame `depth` —
site idxs are GetIdx values in that tagged frame (created atoms often have
no depth-0 entry). Ints coerce to `AtomRef(origin=…)`.

| Rule family | `phase1_steps` |
| --- | --- |
| Phase I (`Hydroxylation`, `Epoxidation`, `Dehydrogenation`, …) and `NDealkylation` | Degenerate singleton `StepPlan` for a valid site |
| `QuinoneFormation` | One plan: `Or` of routes (e.g. DH-only \| dealk→DH \| `And(OH,…)`→DH) |
| Other rules (conjugates, …) | `NotImplementedError` |

### Ask a rule (no enumeration)

```python
from rdkit import Chem
from xenosite.forest import rules

mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
plan = rules.QuinoneFormation().phase1_steps(mol, frozenset({4, 7}))
# e.g. Dehydrogenation[3, 8] | Dealkylation[1, 3] → Dehydrogenation[3, 8]
for branch in plan.branches():
    print(branch)
for lin in plan.linearizations():  # every total order across Or/And/Seq
    ...

rules.Epoxidation().phase1_steps(Chem.MolFromSmiles("C=C"), frozenset({0, 1}))
# → Epoxidation[0, 1]
```

### Stamp while metabolizing

Opt in with `attach_phase1_steps=True`. Read stamps via `StepPlan` — do not inspect
mol props by name (`from_mol` raises if missing; `try_from_mol` returns `None`).

```python
from rdkit import Chem
from xenosite.forest import StepPlan, load_ruleset, rules

mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")

# Single rule
for site, products in rules.QuinoneFormation().metabolize(
    mol, attach_phase1_steps=True, tag_atoms=False
):
    for p in products:
        plan = StepPlan.try_from_mol(p)  # None on unstamped fragments
        if plan is not None:
            print(site, plan)

# Same flag on a ruleset (e.g. "QF", "PhaseOneRS", "Full")
for site, products in load_ruleset("QF").metabolize(
    mol, attach_phase1_steps=True, tag_atoms=False
):
    plan = StepPlan.try_from_mol(products[0])
    ...
```

### Resolve and replay

`AtomRef.resolve(mol)` / `Step.resolve_site(mol)` map refs to the current
`GetIdx` on a given mol (read-only; the plan object is unchanged). Origin refs
carry a frame `depth` — not every atom exists at depth 0 (created atoms like an
OH oxygen only appear mid-path). `Linearization.apply` copies the input mol,
runs Forest rules in order, and records created atoms on
`mol._forest["atom_refs"]` so later created-by refs resolve.

```python
from rdkit import Chem
from xenosite.forest import AtomRef, Linearization, StepPlan, rules

mol = Chem.MolFromSmiles("c1ccccc1")
plan = rules.QuinoneFormation().phase1_steps(mol, frozenset({0, 3}))
branch = next(b for b in plan.branches() if len(b) == 3)

for lin in branch.linearizations():
    products = lin.apply(mol)  # full plan
    mid = lin.apply(mol, drop_last=1)  # prep only; keep frags for final DH sites
    if mid:
        print(lin.steps[-1].resolve_site(mid[0]))

# Same AtomRef works on different mols; only the mol's maps / atom_refs change
ref = AtomRef(origin=3)
ref.resolve(mol)
```

`toward=` on `apply` retains fragments needed for AtomRefs that are not yet applied
(e.g. sites of later steps). `drop_last=N` omits the last N steps but passes their
sites as `toward`.

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

### Path search (`find_path`)

Package-level `find_path` is an MCS-guided search (distinct from classic
`RuleSet.find_path`). The API is **unstable**. Default ruleset is `PhaseOneQF`.

```python
from xenosite.forest import find_path, PathSearchCounters

counters = PathSearchCounters()
for outcome in find_path(
    "CN(C)Cc1ccccc1",
    "O=Cc1ccccc1",
    ruleset="ND",
    counters=counters,
):
    print(outcome.plan, outcome.maybe)
print(counters.as_dict())
```

Hits are `PathOutcome` values (Required `plan` + cleavage-side `maybe`); they
also unpack as `(smiles, steps, mols, plan, maybe)`. Prefer
`PathSearchCounters` when comparing search cost.

| Knob | Default | Role |
| --- | --- | --- |
| `ruleset` | `PhaseOneQF` | Phase I + QuinoneFormation; pass `PhaseOneRS` / `QF` / `ND` / … |
| `depth` | `None` | Unbounded hops; optional hard cap |
| `max_expansions` | `200` | Caps billed work; sets `budget_exhausted` |
| `expand_phase1_plans` | `True` | Expand `phase1_steps` instead of opaque composite hops |

Use classic `bfs` / `RuleSet.find_path` for exhaustive mol-BFS; use package
`find_path` when MCS-guided search is enough. Optimization claims should cite
`PathSearchCounters` shifts.

## Custom rulesets

```python
from xenosite.forest import rules, RuleSet

rs = RuleSet([rules.Epoxidation(), rules.EpoxideOpening()], name="epoxide")
path = next(rs.find_path(reactant, product, depth=2))
```

Built-in names include `Full`, `PhaseOneRS`, `PhaseOneQF` (Phase I + quinone; guided default), `Bioactivation`, `QuinoneFormationRS`, and the Phase I groups `SO`, `DH`, `HD`, `RD`, and `UO`. Pass any of those strings to `bfs(..., ruleset=...)` or `load_ruleset(...)`.

Which ruleset matches which paper (Rainbow, quinone, bioactivation, and others) is in **[rulesets.md](rulesets.md)**.

## Command line

```bash
xenosite-forest CCO CC=O --ruleset PhaseOneRS --depth 1 --phase1
```

See `xenosite-forest --help`.
