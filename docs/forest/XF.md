# `Mol.xf` — molecule accessor API

`mol.xf` is the **public convenience API** for questions about a molecule in
Metabolic Forest: canonical SMILES, SMARTS matches, rings, formula, pair-orbit
keys, finishing products, and atom-trace queries. It is **not** the forest
schema.

Implementation: `Xf` / `XfTracing` in
[`src/xenosite/forest/rdkitutil.py`](../../src/xenosite/forest/rdkitutil.py).
Typing brands (`ForestMol`, `TracingMol`, `NoForestMol`) and the `Xf` protocol
live in
[`src/xenosite/forest/rdkit_api.py`](../../src/xenosite/forest/rdkit_api.py).

Prefer `mol.xf` over deprecated `AtomTracker`. Migration notes for 0.6 → 0.7:
[`MIGRATING_0.7.md`](MIGRATING_0.7.md).

---

## `xf` vs `_forest`

Each access to `mol.xf` **mints a fresh** `Xf` with a strong reference to that
mol. Nothing is stored on the mol under the name `xf`, so facades cannot be
copied between molecules. Temporary chains like
`MolFromSmiles(...).xf.csmi` are safe.

Answers and tracing state are **cached on** `mol._forest` — a private bag the
facade installs lazily when a method needs it. Callers should use `xf`, not
poke `_forest`.

The TypedDict layout of that bag (`Forest`, nested `Structure` / `AtomTrace`,
…) is declared in
[`src/xenosite/forest/records.py`](../../src/xenosite/forest/records.py). Link
there when you need the shape of the cache. **Exact `_forest` field names and
nesting are unstable** — they move as the library evolves. Stable contracts
are the `xf` methods and properties below.

| Layer | Role |
| ----- | ---- |
| `mol.xf` | Public accessor; mint-on-read |
| `mol._forest` | Private cache / trace bag (do not rely on layout) |
| `records.Forest` / `Structure` / `AtomTrace` | TypedDict map of that bag (for maintainers) |

`has_forest` reports whether `_forest` is already present **without**
installing one. Almost every other property installs an empty forest if
missing. Wipe / constructor / reaction pieces return plain mols without a
forest; re-enter with `mol.xf.forestmol` (or any query that attaches).

---

## Quick start

```python
from rdkit import Chem
from xenosite.forest.rules import Hydroxylation
from xenosite.forest.records import Smarts

mol = Chem.MolFromSmiles("Oc1ccccc1")

# Structure answers (cached on first read)
print(mol.xf.csmi)
print(mol.xf.formula)
print(mol.xf.topol_equiv)

# Match-only SMARTS (not SMIRKS)
matches = mol.xf.smarts_matches(Smarts("[#6:1]"))
print(len(matches), matches[0] if matches else None)

# Metabolize products are finished through xf
for products, info in Hydroxylation().metabolize(mol):
    for p in products:
        print(info["site"], p.xf.csmi, p.xf.tracing.depth)
```

---

## Public surface (`Xf`)

### Identity and forest presence

| Member | Kind | Notes |
| ------ | ---- | ----- |
| `mol` | property | Parent RDKit mol (same object). |
| `has_forest` | property | `True` if `_forest` exists; does **not** install. |
| `forestmol` | property | Ensures forest; returns parent typed as `ForestMol`. |
| `forest` | property | The `_forest` dict (installs empty if missing). Prefer named accessors over reading keys. |
| `is_terminal` | property | Terminal product flag (search must not expand further). |
| `clear_structure` | method | Drops structure-dependent cache answers; labels / trace stay. |

```python
assert not Chem.MolFromSmiles("CCO").xf.has_forest  # fresh parse
_ = Chem.MolFromSmiles("CCO").xf.csmi               # attaches + caches
```

### Structure answers (cached)

These fill keys under `_forest["cache"]` (the `Structure` TypedDict) on first
use. After an in-place edit that changes bonding, call `clear_structure` (or
finish via `of_products`, which clears for you).

| Member | Returns | Notes |
| ------ | ------- | ----- |
| `csmi` | `str` | Canonical SMILES (`isomericSmiles=False`). Identity for dedup / search. |
| `formula` | `Formula` | Heavy-atom counts, total H, formal charge. |
| `topol_equiv` | `dict[int, int]` | Atom index → topological equivalence class. |
| `rings` | `dict[int, tuple[tuple[int, ...], ...]]` | Per-atom ring membership. |
| `aromatic_systems` | `tuple[frozenset[int], ...]` | Aromatic connected components. |
| `conjugated_systems` | `tuple[frozenset[int], ...]` | Conjugated atom sets (may load resonance). |
| `sanitize` | `int` | Sanitize a **copy**; caches RDKit status. Does not edit the parent. |
| `smarts_matches` | `tuple[dict[int, int], ...]` | Cached substructure matches for a `Smarts`, keyed by atom map. |

```python
from xenosite.forest.records import Smarts

mol = Chem.MolFromSmiles("c1ccccc1O")
assert mol.xf.csmi == "Oc1ccccc1"
# Map numbers in the SMARTS become keys in each match dict
hits = mol.xf.smarts_matches(Smarts("[#8:1]-[#6:2]"))
# e.g. ({1: oxygen_idx, 2: carbon_idx}, ...)
```

Emission identity in `metabolize` uses a frozenset of fragment `xf.csmi`
values (`info["csmi"]`). Reading `product.xf.csmi` after finishing is cheap
because the value is already cached.

### Pair-orbit keys

Used by unique-edit. Details: [`PAIR_ORBITS.md`](PAIR_ORBITS.md).

| Member | Notes |
| ------ | ----- |
| `pair_orbit_backend` | Resolved backend: `"nauty"` / `"smiles"` / `"none"`. Process-wide. |
| `set_pair_orbit_backend(backend)` | Classmethod; pass `None` to clear override. |
| `atom_pair_orbit_key(site)` | Atom–atom unique-edit signature, or `None`. |
| `bond_pair_orbit_key(bonds)` | Bond–bond signature, or `None`. |
| `site_pair_orbits(backend=None)` | Nested pair-orbit tables (or `None` for `"none"`). |

```python
site = frozenset({0, 3})
key = mol.xf.atom_pair_orbit_key(site)
```

### Finishing products

```python
finished = reactant.xf.of_products(products, site_info, executed=rule)
```

Stamps the reactant, copies each product, records the transform on the child
trace, marks terminal when the rule is a terminal rule, and clears structure
cache on each child. Prefer this over calling underscored tracing plumbing
directly. `metabolize` / `find_path` already finish this way.

---

## Atom tracing (`mol.xf.tracing`)

Nested facade (`XfTracing`). Trace state lives on `_forest["atom_trace"]`;
read it through these methods.

### Read-only

| Member | Notes |
| ------ | ----- |
| `active` | `True` when an atom trace is initialized. |
| `depth` | Transform depth from the root, or `None` if untraced. |
| `atom_indices(idx)` | Index frames from origin to now, or `None`. |
| `atom_depths(idx)` | Depth frames parallel to `atom_indices`. |
| `atom_origin(idx)` | Earliest recorded index (including birth index of an added atom). |
| `atom_root(idx)` | Depth-0 index, or `None` if the atom was created later. |
| `atom_added_by(idx)` | `(rule_name, site)` for a created atom, or `None`. |
| `removed_roots()` | Depth-0 indexes that left the molecule. |

```python
product = products[0]
t = product.xf.tracing
if t.active:
    print(t.depth)
    for atom in product.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        i = atom.GetIdx()
        added = t.atom_added_by(i)
        if added:
            print(i, "added by", added)
```

### Plumbing (underscored)

`_ensure` / `_stamp` / `_install` / `_trace` initialize labels and record
transforms. Application code should use `of_products` instead. Search entry
still calls `tracing._stamp()` on the root reactant.

---

## What not to do

- Do not treat `mol.xf` as synonymous with the `Forest` TypedDict or as a
  stable schema document.
- Do not read or write `mol._forest[...]` in application code. Field layout
  changes; `xf` is the contract.
- Do not assign to `mol.xf` or try to copy an `Xf` instance onto another mol
  (read-only property; each access mints a new facade).
- Do not keep a long-lived `xf` handle across edits that invalidate structure
  answers without `clear_structure` (or a fresh product from `of_products`).

---

## Related docs

| Doc | Topic |
| --- | ----- |
| [`MIGRATING_0.7.md`](MIGRATING_0.7.md) | 0.6 → 0.7 caller changes |
| [`PAIR_ORBITS.md`](PAIR_ORBITS.md) | Pair-orbit unique-edit |
| [`HEURISTICS.md`](HEURISTICS.md) | Find-path / filter policy |
| [`records.py`](../../src/xenosite/forest/records.py) | `_forest` TypedDicts (`Forest`, `Structure`, `AtomTrace`) — unstable layout |
