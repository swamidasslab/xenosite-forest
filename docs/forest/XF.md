# `Mol.xf` — molecule accessor API

`mol.xf` answers common questions about a molecule in Metabolic Forest:
canonical SMILES, SMARTS matches, rings, formula, and atom-trace queries. Each
access returns a short-lived facade over that mol; the methods and properties
below are the stable contract for callers.

Implementation: `Xf` / `XfTracing` in
[`src/xenosite/forest/rdkitutil.py`](../../src/xenosite/forest/rdkitutil.py).
Typing brands (`ForestMol`, `TracingMol`, `NoForestMol`) and the `Xf` protocol
live in
[`src/xenosite/forest/rdkit_api.py`](../../src/xenosite/forest/rdkit_api.py).

Prefer `mol.xf` over deprecated `AtomTracker`. Tutorial for callers migrating
off the tracker: [Replacing AtomTracker](#tutorial-replacing-atomtracker-with-molxf).
Migration notes for 0.6 → 0.7: [`MIGRATING_0.7.md`](MIGRATING_0.7.md).

---

## Caching on `_forest`

Each access to `mol.xf` mints a fresh `Xf` with a strong reference to that mol.
Nothing is stored on the mol under the name `xf`, so facades cannot be copied
between molecules. Temporary chains like `MolFromSmiles(...).xf.csmi` are safe.

Answers and tracing state are cached on `mol._forest` — private storage the
facade installs lazily when a method needs it. Callers should use `xf`, not
read or write `_forest` keys directly.

`has_forest` reports whether `_forest` is already present without installing
one. Almost every other property installs an empty bag if missing. Wipe /
constructor / reaction pieces return plain mols without a forest; re-enter with
`mol.xf.forestmol` (or any query that attaches).

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

These fill keys under `_forest["cache"]` on first use. After an in-place edit
that changes bonding, call `clear_structure` (or finish via `of_products`,
which clears for you).

| Member | Returns | Notes |
| ------ | ------- | ----- |
| `csmi` | `str` | Canonical SMILES (display / target compare). **Not** always safe for dedup — see Chematic `canonical_smiles_stable_key`. Rust: `ForestMol::stable_csmi_key` for search / `unique_csmi`. |
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

Emission identity in `metabolize` is the frozenset of **stable** product keys
when Chematic admits them (`canonical_smiles_stable_key`); if any fragment is
unstable, that emission is not CSMI-deduped (fail-closed). Display spelling
remains `product.xf.csmi` / Rust `csmi`. There is no `info["csmi"]`.

## Atom tracing (`mol.xf.tracing`)

Nested facade (`XfTracing`). Trace state lives on `_forest["atom_trace"]`;
read it through these methods.

### Read-only

| Member | Notes |
| ------ | ----- |
| `active` | `True` when an atom trace is initialized. |
| `depth` | Transform depth from the root, or `None` if untraced. |
| `depths()` | Sorted unique depth frames across heavy atoms (AtomTracker.depths). |
| `atom_indices(idx)` | Index frames from origin to now, or `None`. |
| `atom_depths(idx)` | Depth frames parallel to `atom_indices`. |
| `index_at(idx, depth)` | Index of atom `idx` at `depth`, or `None`. |
| `atom_origin(idx)` | Earliest recorded index (including birth index of an added atom). |
| `atom_root(idx)` | Depth-0 index, or `None` if the atom was created later. |
| `atom_added_by(idx)` | `(rule_name, site)` for a created atom, or `None`. |
| `added_indices()` | Current heavy-atom indexes with no depth-0 root. |
| `root_map()` | Current index → depth-0 root for atoms that have one. |
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
transforms. Application code should call `metabolize` / `find_path` (they
finish products internally) rather than these mutators. Search entry still
calls `tracing._stamp()` on the root reactant.

---

## Tutorial: replacing AtomTracker with `mol.xf`

`AtomTracker` was the 0.6 / archive way to stamp atom labels, read
label→`{idx, depth, …}` tag dicts, and ask topological ranks. In 0.7 the live
class at `xenosite.forest.AtomTracker` is a **deprecated thin facade**: most
helpers already call `mol.xf` / `mol.xf.tracing` and emit
`DeprecationWarning`. New code should use `xf` directly. The patterns below
match what callers used to write against the tracker.

Indexes on `xf.tracing` are **0-based RDKit `GetIdx()`** values (same scale as
`info["site"]`). The archived `AtomTrace` helper converted to 1-based SMILES
map numbers; if you still need that convention, add one yourself (`idx + 1`).

### Capability map

| AtomTracker (archive / facade) | `mol.xf` equivalent |
| ------------------------------ | ------------------- |
| `AtomTracker.topol_equiv(mol)` | `mol.xf.topol_equiv` |
| `AtomTracker.site_to_topol_site((name, idxs), te)` | Build ranks from `mol.xf.topol_equiv` (no dedicated helper; site shape also changed in 0.7) |
| `AtomTracker().initialize_tags(mol)` | `mol.xf.tracing._stamp()` — usually unnecessary; `metabolize` already stamps |
| `AtomTracker.tags(mol)` label dicts | Per-atom: `tracing.atom_indices(i)` / `atom_depths(i)` over heavy atoms |
| `tags(mol, depth=d)` / `tags(mol, idx=i)` | Filter those per-atom histories yourself |
| `tags(mol, compact=True)` / `compact_tags` | Zip depths↔indices and add `1` if you want 1-based map numbers |
| `AtomTracker.depths(mol)` | Current step: `tracing.depth`. Full set: `tracing.depths()` |
| “Atom has no depth 0” (added O, …) | `tracing.added_indices()` / `atom_added_by(i)` |
| Removed depth-0 atoms | `tracing.removed_roots()` |
| `metabolite_index_to_reversed_index_record` | `tracing.root_map()` / `tracing.index_at(i, depth)` |
| Tagging kwargs on `metabolize` (`tag_atoms`, …) | Gone; `metabolize` finishes products for you |

The facade still exposes `tags` / `depths` / `compact_tags` for call-site
compatibility. Prefer the per-atom tracing queries: they do not expose
internal `forestLabel` strings, and they stay stable as the cache layout
moves.

### Setup: a traced product

```python
from rdkit import Chem
from xenosite.forest.rules import Hydroxylation

reactant = Chem.MolFromSmiles("CCO")
products, info = next(Hydroxylation().metabolize(reactant))
product = products[0]
t = product.xf.tracing

assert t.active
assert t.depth == 1
assert info["site"] == frozenset({0})  # 0-based indexes on the reactant
print(product.xf.csmi)  # e.g. OCCO
```

Manual stamping (search entry / rare custom pipelines):

```python
mol = Chem.MolFromSmiles("CCO")
mol.xf.tracing._stamp()  # depth-0 records for each heavy atom
assert mol.xf.tracing.depth == 0
```

Application code almost never needs `_stamp`: `metabolize` installs the
reactant trace and records each child transform.

### Topological ranks and sites

**Before (AtomTracker):**

```python
from xenosite.forest import AtomTracker  # deprecated

te = AtomTracker.topol_equiv(mol)
# archive sites often looked like ("Hydroxylation_…", (0, 2))
bare, ranks = AtomTracker.site_to_topol_site(
    ("Hydroxylation_SmartsReactionRuleRxn0", (0, 2)), te
)
```

**After (`xf`):**

```python
te = mol.xf.topol_equiv  # dict[int, int]: GetIdx → class
site = info["site"]      # frozenset[int] on live metabolize
ranks = tuple(sorted(te[i] for i in site))
# Rule name lives on info["rule"][0].name — not packed into the site tuple
```

### Tag dicts → per-atom history

Archive callers walked a label→record map:

```python
# Before
tags = AtomTracker.tags(product)  # {label: {"idx": [...], "depth": [...]}, ...}
for label, rec in tags.items():
    if 0 not in rec["depth"]:
        print("added atom path", rec["idx"])
```

With `xf`, ask each current heavy atom:

```python
# After
t = product.xf.tracing
for atom in product.GetAtoms():
    if atom.GetAtomicNum() == 1:
        continue
    i = atom.GetIdx()
    idxs = t.atom_indices(i)    # e.g. (1,) for new O, or (0, 0) for a carried C
    depths = t.atom_depths(i)   # parallel frames, e.g. (1,) or (0, 1)
    root = t.atom_root(i)       # depth-0 index, or None if created later
    origin = t.atom_origin(i)   # earliest recorded index (birth index if added)
    added = t.atom_added_by(i)  # ("Hydroxylation", frozenset({0})) or None
```

On ethanol hydroxylation (`CCO` → `OCCO`), the new oxygen has
`atom_root is None`, `atom_depths == (1,)`, and
`atom_added_by == ("Hydroxylation", frozenset({0}))`. Carried carbons keep a
depth-0 root and a two-frame `atom_indices` / `atom_depths` history.

### Depths on a path product

**Before:** `AtomTracker.depths(product)` → sorted unique depths across all
tag records (e.g. `[0, 1]` after one step).

**After:**

```python
print(product.xf.tracing.depth)    # current step (None if untraced)
print(product.xf.tracing.depths())  # e.g. (0, 1) — same as AtomTracker.depths
```

### Parent mapping (reverse index)

`metabolite_index_to_reversed_index_record` built a compact 1-based
depth→index table and reversed the last few depths. The usual need is
“where was this product atom on the reactant?”:

```python
t = product.xf.tracing
parent_of = t.root_map()  # current GetIdx → depth-0 GetIdx on the reactant
# Specific earlier depth (not only 0):
prev = t.index_at(i, depth=0)
```

Atoms that left the molecule appear in `t.removed_roots()` as depth-0 indexes
(not current product indexes).

### Compact / 1-based maps

Old `compact=True` tags were `{label: {depth: idx+1, …}}` for SMILES map
convention. Rebuild only if something downstream still expects that shape:

```python
compact = {}
for atom in product.GetAtoms():
    if atom.GetAtomicNum() == 1:
        continue
    i = atom.GetIdx()
    depths = product.xf.tracing.atom_depths(i)
    idxs = product.xf.tracing.atom_indices(i)
    if not depths or not idxs:
        continue
    compact[i] = {d: j + 1 for d, j in zip(depths, idxs)}
```

### Worked example side by side

Ethanol → primary alcohol hydroxylation, asking “which atoms are new?”:

```python
from rdkit import Chem
from xenosite.forest.rules import Hydroxylation

# --- Before (deprecated facade; still works, warns) ---
from xenosite.forest import AtomTracker
import warnings

parent = Chem.MolFromSmiles("CCO")
products, _info = next(Hydroxylation().metabolize(parent))
child = products[0]
with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    tags = AtomTracker.tags(child)
    new_via_tags = [rec for rec in tags.values() if 0 not in rec["depth"]]
assert new_via_tags  # added oxygen has no depth-0 frame

# --- After (preferred) ---
t = child.xf.tracing
new_idxs = sorted(t.added_indices())
assert new_idxs
print(new_idxs[0], t.atom_added_by(new_idxs[0]))
# e.g. 1 ('Hydroxylation', frozenset({0}))
```

Same product, topological site identity for unique-edit style callers:

```python
te = child.xf.topol_equiv
# Live info["site"] is on the *reactant* frame for that emission; for ranks
# on the product molecule, use current indexes you care about:
product_site = frozenset(new_idxs)
ranks = tuple(sorted(te[i] for i in product_site))
```

### When the facade is still useful

Keep `AtomTracker` only while porting call sites that need the old
label-dict shape or `metabolite_index_to_reversed_index_record` unchanged.
It will be removed in a future release. Do not subclass it (live rules never
did; see [`DROPPED.md`](DROPPED.md)). For 0.6 → 0.7 yield / filter changes
around tagging kwargs, see [`MIGRATING_0.7.md`](MIGRATING_0.7.md).

---

## What not to do

- Do not read or write `mol._forest[...]` in application code; field layout
  changes as the library evolves.
- Do not assign to `mol.xf` or try to copy an `Xf` instance onto another mol
  (read-only property; each access mints a new facade).
- Do not keep a long-lived `xf` handle across edits that invalidate structure
  answers without `clear_structure` (or a fresh `metabolize` product).

---

## Internals

Maintainers: the TypedDict layout of `_forest` (`Forest`, nested `Structure` /
`AtomTrace`, …) is declared in
[`src/xenosite/forest/records.py`](../../src/xenosite/forest/records.py). Exact
field names and nesting are an unstable cache layout; prefer the `xf` surface
above when writing application code.

| Layer | Role |
| ----- | ---- |
| `mol.xf` | Public accessor; mint-on-read |
| `mol._forest` | Private cache / trace bag |
| `records.Forest` / `Structure` / `AtomTrace` | TypedDict map of that bag |

Product finishing (`xf._of_products`) and pair-orbit helpers
(`xf._atom_pair_orbit_key`, …) stay **forest-internal** — leading underscore,
not a stable public contract. Unique-edit still uses them from `metabolize`;
callers should not. Pair orbits require **pynauty** (always nauty):
[`PAIR_ORBITS.md`](PAIR_ORBITS.md).

---

## Related docs

| Doc | Topic |
| --- | ----- |
| [`MIGRATING_0.7.md`](MIGRATING_0.7.md) | 0.6 → 0.7 caller changes |
| [`PAIR_ORBITS.md`](PAIR_ORBITS.md) | Pair-orbit unique-edit |
| [`HEURISTICS.md`](HEURISTICS.md) | Find-path / filter policy |
| [`records.py`](../../src/xenosite/forest/records.py) | `_forest` TypedDicts — unstable cache layout |
