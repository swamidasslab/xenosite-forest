# Migrating to 0.7.0

For callers of **0.6.x** moving to **`xenosite.forest` 0.7**. The previous
forest implementation lives in the archive directory
[`src/xenosite/_archive_forest/`](../../src/xenosite/_archive_forest/) on
GitHub. Short examples only; design detail is in sibling notes under
`docs/forest/`.

## `metabolize` yield shape (breaking)

### Was (0.6 / archive)

Archive `ReactionRule.metabolize` yields **`(site, metabolites)`** where
`metabolites` is always a **`list[Mol]`** — never a bare mol
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
`RuleSet.metabolize` yield **`(product, info)`** — one mol per yield
(`src/xenosite/forest/rules.py` ~238–245, ~405; `rulesets.py` ~113–121).
`info` is `ProductInfo` (`site`, `rule`, `options`, `csmi`, …). Cleavage:
one yield per fragment (`product_index` / `product_count`); identical
fragment CSMIs may drop under default `unique_csmi=True`.

```python
from rdkit import Chem
from xenosite.forest.rules import Hydroxylation

mol = Chem.MolFromSmiles("Oc1ccccc1")
for product, info in Hydroxylation().metabolize(mol):
    print(info["site"], info["csmi"])
```

| | 0.6 / archive | 0.7 |
| --- | --- | --- |
| Yield | `(site, list[Mol])` | `(product, info)` |
| Site | Often `(rule_name, frozenset[…])` | Atom indexes only (`info["site"]`) |
| Cleavage | Sibling mols in one list | One yield per fragment |
| Product SMILES | Caller computes | `info["csmi"]` |
| Dedup knobs | Topo-site + optional `only_unique` | `unique_csmi` (default `True`) |

`ReactionRule.__call__` matches `metabolize`.

### Adapting a 0.6 loop

```python
# Old — products is always a list (often [one_mol])
for site, products in rule.metabolize(mol, tag_atoms=False):
    for p in products:
        use(site, p)

# New — one mol per iteration
for product, info in rule.metabolize(mol):
    site = info["site"]          # not (name, site)
    rule_name = info["rule"].name
    use(site, product)
```

To rebuild the old list-per-edit shape:

```python
from collections import defaultdict

by_edit = defaultdict(list)
for product, info in rule.metabolize(mol, unique_csmi=False):
    key = (
        info["site"],
        info["pattern"]["name"] if "pattern" in info else None,
    )
    by_edit[key].append(product)
```

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
- **`directed_bond`** rules (Dealkylation, NDealkylation, …): orientation is
  on `info["discovered_site"]` as an ordered map-order `tuple`. Unique-edit
  still uses directed ranks internally.
- Opt-in `canonical_emitted_sites=True`: emitted `site` may be the lex orbit
  representative; `discovered_site` keeps discovery indexes when they differ.
  `filter_sites` always sees **discovery**.

```python
for product, info in Dealkylation().metabolize(mol):
    site = info["site"]                      # frozenset
    oriented = info.get("discovered_site")   # tuple for directed_bond
```

Rules declare `site_kind`: `"atom"` | `"bond"` | `"directed_bond"` |
`"atom_pair"`. Callers usually read `info["site"]` only.

---

## Filters and tracing

| 0.6 | 0.7 |
| --- | --- |
| `include_sites` / `exclude_sites` | `filter_sites(mol, site, info) -> bool` |
| (no pattern gate) | `filter_rules(mol, rule, pattern_info) -> bool` |
| `AtomTracker` tagging kwargs | Prefer `mol.xf` / `mol.xf.tracing` (`AtomTracker` is deprecated) |

Methide is always on Dehydrogenation / QuinoneFormation data. There is no
`pathways=("methide",)` flag — refuse with `filter_sites` that reads
`info["options"]["methide"]`.

`None` mol raises `ValueError` (no soft empty iterator).

---

## Product CSMI: check vs yield

- **`unique_csmi=True` (default):** drop duplicate `(rule, pattern, csmi)`
  fragments from the stream (quiet for leaving groups / unequal-rank product
  iso).
- **Check (always while site unique-edit is on):** same emission frozenset of
  fragment CSMIs **and** matching site ranks under `(rule, pattern)` →
  `SiteDeduplicationWarning` (unique-edit miss). Check does **not** drop.
- **`RuleSet.metabolize`:** children run with yield off so alternate rules
  bubble up; the set’s own `unique_csmi` is the cross-child layer. Same CSMI
  from two **different** child rules → **INFO** log (`kept_rule` /
  `dropped_rule` / CSMI), not a Warning. Inventory: `REDUNDANT_RULES.md`.

Do not treat cross-rule INFO as a unique-edit failure. Do not ignore
`SiteDeduplicationWarning` in tests that intend to catch signature bugs.

---

## Naming: SMARTS vs SMIRKS

| Was | Is |
| --- | --- |
| `SmartsReactionRule` | `SmirksReactionRule` |
| `.smarts` (reaction string) | `.smirks` |
| Match-only queries | Still `Smarts` / `smarts_matches` / `endpoints` |

---

## Plans and path search

- Live rules expose `canonical_plan(mol, info)`, not archive
  `phase1_steps(mol, site)` on every rule. `StepPlan` / `AtomRef` /
  `Linearization` are still re-exported from
  `xenosite.forest` (implementation under the archive directory for now).
- `attach_phase1_steps=True` on `metabolize` is gone.
- `find_path(reactant, target, ruleset=..., max_nodes=..., counters=...,
  use_filters=..., **kwargs)` — plan-guided search. Older guided knobs
  (`max_expansions` as the public mode, PathContext formula algebra, …) are
  not the live API; see `DROPPED.md`.
- Package `bfs` / `dfs` / `find_path` live under `xenosite.forest`.

---

## Package layout

- Import **`xenosite.forest`**.
- Previous forest: [`src/xenosite/_archive_forest/`](../../src/xenosite/_archive_forest/)
  on GitHub (not exercised by CI).
- Design notes: `docs/forest/` (`HEURISTICS`, `DIVERGENCES`, `DROPPED`,
  `PAIR_ORBITS`, `PERFORMANCE`, `REDUNDANT_RULES`).

---

## Quick checklist

1. Change loops from `site, products` → `product, info`.
2. Read `info["site"]` / `info["csmi"]`; drop `(rule_name, site)` unpacking.
3. Replace `include_sites` / tagging kwargs with `filter_sites` / `mol.xf`.
4. Rename `SmartsReactionRule` / `.smarts` → `SmirksReactionRule` / `.smirks`.
5. Expect `SiteDeduplicationWarning` on true unique-edit misses; cross-rule overlaps
   are INFO only.
6. Prefer `canonical_plan` over `phase1_steps` on live rules.
