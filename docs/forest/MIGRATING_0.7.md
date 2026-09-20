# Migrating to 0.7.0

For callers of **0.6.x** moving to `xenosite.forest` **0.7**. The previous
forest implementation lives in the archive directory
`[src/xenosite/_archive_forest/](../../src/xenosite/_archive_forest/)` on
GitHub. Short examples only; design detail is in sibling notes under
`docs/forest/`.

## `metabolize` yield shape (breaking)

### Was (0.6 / archive)

Archive `ReactionRule.metabolize` yields `(site, metabolites)` where
`metabolites` is always a `list[Mol]` — never a bare mol
(`src/xenosite/_archive_forest/base.py` ~2384). Non-cleavage edits are still
a one-element list. Cleavage (e.g. azo) packs sibling fragments in the same
list. `__call__` matches (`base.py` ~2185–2187). Doctests unpack
`site, metabolites = next(...metabolize(mol))`
(`_archive_forest/rules.py`).

```python
for site, products in rule.metabolize(mol):
    # site: ("Hydroxylation", frozenset({2})) when format_output_site=True
    # products: list[Mol]  — len 1 for hydroxylation; len ≥ 2 for cleavage
    for p in products:
        ...
```

Common kwargs gone on live `metabolize`:
`only_emit_topologically_distinct_sites`, `tag_atoms` / `do_not_tag_atoms`,
`format_output_site`, `only_unique`, `only_largest_fragment`,
`attach_phase1_steps`, `include_sites` / `exclude_sites`.

### Is (0.7)

Live `ReactionRule.metabolize` /
`RuleSet.metabolize` yield `(products, info)` — always a `list[Mol]`
per emission (`src/xenosite/forest/rules.py` metabolize;
`rulesets.py`). Non-cleavage: one-element list. Cleavage: all sibling
fragments in that list. `info` is `ProductInfo` (= `SiteInfo`: `site`, `rule`,
`options`, …). Product SMILES are **not** on `info` — use `product.xf.csmi`,
or `frozenset(p.xf.csmi for p in products)` for emission identity (same `S`
as unique-edit check / `unique_csmi` yield). No `product_index` /
`product_count`. Default `unique_csmi=True` drops a later emission with the
same frozenset under the same `(rule, pattern)`; identical siblings inside
one cleavage list are kept.

```python
from rdkit import Chem
from xenosite.forest.rules import Hydroxylation

mol = Chem.MolFromSmiles("Oc1ccccc1")
for products, info in Hydroxylation().metabolize(mol):
    # products: list[Mol] — len 1 here
    print(info["site"], [p.xf.csmi for p in products])
```


|                | 0.6 / archive                      | 0.7                                            |
| -------------- | ---------------------------------- | ---------------------------------------------- |
| Yield          | `(site, list[Mol])`                | `(list[Mol], info)`                            |
| Site           | Often `(rule_name, frozenset[…])`  | Atom indexes only (`info["site"]`)             |
| Cleavage       | Sibling mols in one list           | Sibling mols in one list                       |
| Product SMILES | Caller computes                    | `product.xf.csmi` (emission: frozenset of those) |
| Dedup knobs    | Topo-site + optional `only_unique` | `unique_csmi` (default `True`, emission-level) |


`ReactionRule.__call__` matches `metabolize`.

### Adapting a 0.6 loop

```python
# Old — products is always a list (often [one_mol])
for site, products in rule.metabolize(mol, tag_atoms=False):
    for p in products:
        use(site, p)

# New — still a list per emission; site lives on info
for products, info in rule.metabolize(mol):
    site = info["site"]          # not (name, site)
    rule_name = info["rule"].name
    for product in products:
        use(site, product)
```

`bfs` / `dfs` still walk **one mol per node**; they unpack each emission
list into separate frontier children.

### `metabolites` (chemistry layer)

Yields `ProductsOfReaction(info=SiteInfo, products=list[Mol])`, not
`(site, products)`. Prefer `metabolize` for application code.

```python
for por in Hydroxylation().metabolites(mol):
    print(por.info["site"], por.products)
```

---



## Sites and `discovered_site`

- Public `info["site"]` is always unordered-shaped for bonds: a `frozenset`
of atom indexes (singleton atom sites are still a frozenset of one index
from SMIRKS rules).
- `directed_bond` rules (Dealkylation, NDealkylation, …): orientation is
on `info["discovered_site"]` as an ordered map-order `tuple`. Unique-edit
still uses directed ranks internally.
- Opt-in `canonical_emitted_sites=True`: emitted `site` may be the lex orbit
representative; `discovered_site` keeps discovery indexes when they differ.
`filter_sites` always sees **discovery**.

```python
for products, info in Dealkylation().metabolize(mol):
    site = info["site"]                      # frozenset
    oriented = info.get("discovered_site")   # tuple for directed_bond
```

Rules declare `site_kind`: `"atom"` | `"bond"` | `"directed_bond"` |
`"atom_pair"`. Callers usually read `info["site"]` only.

---



## Filters and tracing


| 0.6                               | 0.7                                                              |
| --------------------------------- | ---------------------------------------------------------------- |
| `include_sites` / `exclude_sites` | `filter_sites(mol, site, info) -> bool`                          |
| (no pattern gate)                 | `filter_rules(mol, rule, pattern_info) -> bool`                  |
| `AtomTracker` tagging kwargs      | Prefer `mol.xf` / `mol.xf.tracing` (`AtomTracker` is deprecated) |


`mol.xf` is the public molecule accessor (csmi, matches, tracing, …) — not
the forest schema. It caches into private `_forest`. Full surface:
[`XF.md`](XF.md). Cache layout TypedDicts (`Forest` / `Structure` /
`AtomTrace`) live in
[`records.py`](../../src/xenosite/forest/records.py); exact `_forest`
details are unstable — use `xf`, do not poke `_forest`.

Methide is always on Dehydrogenation / QuinoneFormation data. There is no
`pathways=("methide",)` flag — refuse with `filter_sites` that reads
`info["options"]["methide"]`.

`None` mol raises `ValueError` (no soft empty iterator).

---



## Product CSMI: check vs yield

- `unique_csmi=True` **(default):** drop duplicate
`(rule, pattern, emission frozenset)` emissions from the stream (same `S`
as check). Sibling fragments in one cleavage list stay together.
- **Check (always while site unique-edit is on):** same emission frozenset of
fragment CSMIs **and** matching site ranks under `(rule, pattern)` →
`SiteDeduplicationWarning` (unique-edit miss). Check does **not** drop.
- `RuleSet.metabolize`**:** children run with yield off so alternate rules
bubble up; the set’s own `unique_csmi` is the cross-child layer. Same
emission frozenset from two **different** child rules → **INFO** log
(`kept_rule` / `dropped_rule` / CSMI), not a Warning. Inventory:
`REDUNDANT_RULES.md`.

Do not treat cross-rule INFO as a unique-edit failure. Do not ignore
`SiteDeduplicationWarning` in tests that intend to catch signature bugs.

---



## Naming: SMARTS vs SMIRKS


| Was                         | Is                                              |
| --------------------------- | ----------------------------------------------- |
| `SmartsReactionRule`        | `SmirksReactionRule`                            |
| `.smarts` (reaction string) | `.smirks`                                       |
| Match-only queries          | Still `Smarts` / `smarts_matches` / `endpoints` |


---



## Plans and path search

- Live rules expose `canonical_plan(mol, info)`, not archive
`phase1_steps(mol, site)` on every rule. `StepPlan` / `AtomRef` /
`Linearization` are still re-exported from
`xenosite.forest` (implementation under the archive directory for now).
- `attach_phase1_steps=True` on `metabolize` is gone.
- `find_path(reactant, target, ruleset=..., max_nodes=..., counters=..., use_filters=..., **kwargs)` — plan-guided search. Older guided knobs
(`max_expansions` as the public mode, PathContext formula algebra, …) are
not the live API; see `DROPPED.md`.
- Package `bfs` / `dfs` / `find_path` live under `xenosite.forest`.

---



## Package layout

- Import `xenosite.forest`.
- Previous forest: `[src/xenosite/_archive_forest/](../../src/xenosite/_archive_forest/)`
on GitHub (not exercised by CI).
- Design notes: `docs/forest/` (`XF`, `HEURISTICS`, `DIVERGENCES`, `DROPPED`,
`PAIR_ORBITS`, `PERFORMANCE`, `REDUNDANT_RULES`).

---



## Quick checklist

1. Change loops from `site, products` → `products, info` (still a list per emission).
2. Read `info["site"]`; product SMILES via `product.xf.csmi` (or
  `frozenset(p.xf.csmi for p in products)` for emission identity). Drop
  `(rule_name, site)` unpacking — there is no `info["csmi"]`.
3. Replace `include_sites` / tagging kwargs with `filter_sites` / `mol.xf`.
4. Rename `SmartsReactionRule` / `.smarts` → `SmirksReactionRule` / `.smirks`.
5. Expect `SiteDeduplicationWarning` on true unique-edit misses; cross-rule overlaps
  are INFO only.
6. Prefer `canonical_plan` over `phase1_steps` on live rules.

