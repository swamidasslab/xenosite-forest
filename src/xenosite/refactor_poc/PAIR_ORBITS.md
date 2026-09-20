# Pair-orbit unique-edit (agent briefing)

Source of truth at HEAD: `graph_isomorphism.py` + signature types in `records.py`.
Rule chemistry stays in `rules.py` (thin call sites only). HEURISTICS Status: approved
for `swap_group` / nauty six-family / DH `unique_orbit="bond_atom"`.

This doc is for another agent. Facts and excerpts only.

---

## 1. Problem unique-edit solves

**Collapse isomorphic edits; do not collapse distinct chemistry.**

Two layers (do not conflate — HEURISTICS):

1. **Site unique-edit** — topological ranks + pair-orbit signature. Equivalent
   embeddings of the *same* edit share one key before `RunReactants`.
2. **Product csmi dedup** (`unique_csmi`) — drops site topology; keys
   `(rule name, PatternInfo.name | SMARTS, product csmi)`.

Pair orbits refine (1) when atom ranks alone are not enough. Classic negative:
benzene meta vs para share ranks but not atom–atom orbits.

```120:136:tests/refactor_poc/test_graph_isomorphism.py
def test_benzene_meta_and_para_share_ranks_but_not_atom_pair_orbits():
    """Negative: meta and para share ranks but not pair_group_id."""

    mol = _mol("c1ccccc1")
    meta, other_meta, para = (0, 2), (1, 3), (0, 3)
    assert _rank_key(mol, *meta) == _rank_key(mol, *para)
    assert atom_pair_orbit_isotope(mol, *meta) != atom_pair_orbit_isotope(mol, *para)
    ...
```

Core signature shape (not a `Site`): `((ga, gb), pair_group_id)`.

```431:456:src/xenosite/refactor_poc/records.py
# Unique-edit keys, **not** :data:`Site` definitions. Sites name concrete
# atom/bond indexes; these name isomorphism classes of an *edit*.
#
# Core shape: (groups=(ga, gb), pair_group_id).
# ``pair_group_id`` is a sequential int (``PairGroupId``) from CIP-sorted
# orbit membership tuples, or ``TRIVIAL_PAIR_GROUP`` when either end is a
# singleton topeqiv class. Never a bare atom/bond index.
#
# Ordered vs unordered is **tagged on the type**, not only which nauty family
# was queried. ``ordered`` is part of NamedTuple equality:
#   - ``ordered=False``: ends are swappable (shared resolved ``swap_group``);
#     ``end_ranks`` is ``()``.
#   - ``ordered=True``: ends have distinct roles; ``end_ranks`` holds site
#     topeqiv ranks in canonical ``PatternInfo.name`` order (atom/bond pairs).
```

---

## 2. Same-kind pairs (atom–atom / bond–bond)

Ordered vs unordered is decided by resolved **`swap_group`** (default = `name`):

```965:1003:src/xenosite/refactor_poc/graph_isomorphism.py
def resolved_swap_group(
    info: PatternInfo, effect: Effect | None = None
) -> str | None:
    """Resolved swap group: When override, else PatternInfo, else ``name``.

    Default is ``name`` so same-role pair ends are unordered without an
    explicit annotation. Set ``swap_group`` only when grouping differs from
    ``name`` (see HEURISTICS).
    """
    ...

def ends_swappable(...) -> bool:
    """True → unordered pair orbit; False → ordered.

    Reads resolved ``swap_group`` (When → PatternInfo → ``name``). Unordered
    when both ends share the same non-empty group. Unequal → ordered.
    Canonical pattern order for ordered keys uses ``PatternInfo.name``.
    """
    g1 = resolved_swap_group(info1, effect1)
    g2 = resolved_swap_group(info2, effect2)
    return g1 is not None and g1 == g2
```

Schema (approved):

```200:221:src/xenosite/refactor_poc/records.py
    ``swap_group`` (optional): pair unique-edit ordered vs unordered. Resolved
    as When override → this field → ``name`` (``None`` / absent means use
    ``name``). Same non-empty resolved group on both ends → unordered
    (swappable); unequal → ordered, ends sorted by ``name``. Set explicitly
    only when grouping should differ from ``name``. ...
    # Pair unique-edit group. Optional override when grouping differs from
    # ``name``; None / absent → use ``name``. Same resolved group on both
    # ends → unordered orbit; else ordered by ``name``. When may override.
    swap_group: str
```

Nauty families for same-kind (source of truth):

| Family | Candidates |
|--------|------------|
| `atom_atom_unordered` | combinations `(a1,a2)` with `a1 < a2` |
| `atom_atom_ordered` | permutations `(a1,a2)` with `a1 != a2` |
| `bond_bond_unordered` | combinations of bonds |
| `bond_bond_ordered` | permutations of bonds |

Asymmetric PatternInfo end roles → ordered family + `ordered=True` on the
signature; symmetric (shared `swap_group`) → unordered + `ordered=False`.

Examples: DH `phenol_end`×`phenol_end` unordered; `phenol_end`×`amine_end`
ordered (despite shared `edit=single_to_double`). QF same-name pairs
(`add_carbonyl_o`×2) unordered; cross-role ordered.

---

## 3. Bond–atom: directed, not a Site

Bond role ≠ atom role. A bond–atom unique-edit key is **directed** `(bond, atom)`.
It is **not** a `Site` definition. Sites name concrete indexes; this names an
isomorphism class of an *edit*.

**Surprising but intentional:** one directed unit captures **both**

- topological site identity (orbit `groups` + `pair_group`), and
- edit direction (bond vs atom role).

Direction is part of the edit's isomorphism class — not a separate Site field.
Same-kind atom–atom / bond–bond have no type-level direction; they use
`swap_group` ordered/unordered only.

```505:534:src/xenosite/refactor_poc/records.py
class BondAtomOrbitSignature(NamedTuple):
    """Directed (bond, atom) unique-edit key — **not** a :data:`Site`.

    ``groups`` is ``(bond_group, atom_group)``. Used when
    ``unique_orbit="bond_atom"`` (e.g. Dehydrogenation).

    Surprising but correct: this one directed unit captures **both**
    topological site identity (orbit groups) **and** edit direction (bond
    role vs atom role). Direction is part of the edit's isomorphism class,
    not a separate Site field. ...
    """

    groups: BondAtomGroupPair
    pair_group: PairGroupId


class BondAtomPairOrbitSignature(NamedTuple):
    """Two :class:`BondAtomOrbitSignature` ends (ResonancePair unique-edit).

    Unique-edit key, not a Site. Each end is already directed bond→atom;
    ``ordered`` only decides whether the *two ends* are swappable
    (``swap_group``), not whether bond and atom roles swap.
    """

    ends: tuple[BondAtomOrbitSignature, BondAtomOrbitSignature]
    ordered: Literal[True, False]
```

Dehydrogenation declares the mode as data:

```2057:2062:src/xenosite/refactor_poc/rules.py
    Unique-edit uses ``bond_atom`` orbits (the C–X bond plus the site atom).
    Pair-path ends with the same PatternInfo role are unordered; different
    roles (phenol vs amine) are ordered so ``(a,b)`` ≠ ``(b,a)``.
    """

    unique_orbit: UniqueOrbit = "bond_atom"
```

---

## 4. What nauty emits

`all_site_pair_orbits_nauty` is the source of truth (six families). Cross-kind
builds `product(atoms, bonds)` once; ordered ≡ unordered for that partition.

```417:430:src/xenosite/refactor_poc/graph_isomorphism.py
    # Cross-kind: one partition for both ordered and unordered keys.
    atom_bond_candidates = list(product(all_atoms, all_bonds))
    atom_bond_groups = _orbit_groups_from_generators(
        candidates,
        generators,
        kind1="atom",
        kind2="bond",
        ordered=True,  # product images stay (atom, bond); no sort
        ...
    )
    output["atom_bond_unordered"] = atom_bond_groups
    output["atom_bond_ordered"] = atom_bond_groups
```

Legacy three-mode tables flip to `(bond, atom)` for unique-edit compatibility:

```513:518:src/xenosite/refactor_poc/graph_isomorphism.py
        elif mode == "bond_atom":
            # Families use (atom, bond); PairMode tables use (bond, atom).
            raw = [
                [(b, a) for a, b in group]
                for group in families["atom_bond_unordered"]
            ]
```

Nested forest tables: `mode → (ga, gb) → {(end_a, end_b): PairGroupId}`.

- **`groups`**: topo class ids (`TopoGroupId`) — atom topeqiv or bond class.
  Separate from orbit membership. Same-kind: sorted `(ga, gb)`. Bond–atom:
  `(bond_group, atom_group)` directed.
- **`pair_group`**: sequential `PairGroupId` `0..n-1` from CIP-sorted orbit
  membership tuples within each mode (or `TRIVIAL_PAIR_GROUP` — see caveats).

```838:873:src/xenosite/refactor_poc/graph_isomorphism.py
def _nested_tables_from_groups(...) -> SitePairOrbitTables:
    """Build mode -> (ga, gb) -> {(a, b): PairGroupId}.

    Groups are numbered ``0..n-1`` by sorting CIP membership keys within each
    mode. Nauty and isotope agree on these ids when partitions agree.
    """
    ...
                by_groups[_group_pair(mode, g_left, g_right)][idx_pair] = pair_group
```

---

## 5. Types

| Type | Role |
|------|------|
| `AtomPairOrbitSignature` | atom–atom; `ordered` + `end_ranks` (`()` if unordered) |
| `BondPairOrbitSignature` | bond–bond; same shape |
| `BondAtomOrbitSignature` | one directed `(bond, atom)` unit |
| `BondAtomPairOrbitSignature` | two bond–atom units (ResonancePair); own `ordered` |

```487:545:src/xenosite/refactor_poc/records.py
class AtomPairOrbitSignature(NamedTuple):
    groups: AtomGroupPair
    pair_group: PairGroupId
    ordered: Literal[True, False]
    end_ranks: tuple[int, ...]  # () if unordered; (ra, rb) name-order if ordered


class BondPairOrbitSignature(NamedTuple):
    groups: BondGroupPair
    pair_group: PairGroupId
    ordered: Literal[True, False]
    end_ranks: tuple[int, ...]


class BondAtomOrbitSignature(NamedTuple):
    groups: BondAtomGroupPair
    pair_group: PairGroupId


class BondAtomPairOrbitSignature(NamedTuple):
    ends: tuple[BondAtomOrbitSignature, BondAtomOrbitSignature]
    ordered: Literal[True, False]


PairOrbitSignature: TypeAlias = (
    AtomPairOrbitSignature
    | BondPairOrbitSignature
    | BondAtomOrbitSignature
    | BondAtomPairOrbitSignature
)
```

Constructors in `graph_isomorphism`: `unordered_atom_pair_orbit` /
`ordered_atom_pair_orbit` / `unordered_bond_pair_orbit` /
`ordered_bond_pair_orbit` / `unordered_bond_atom_pair` /
`ordered_bond_atom_pair`. Lookup helpers: `atom_pair_orbit_key`,
`bond_pair_orbit_key`, `bond_atom_orbit_key`.

Assembly for emissions: `site_orbit` / `pair_orbit` → `site_signature` /
`pair_site_signature` (last field is the `PairOrbitSignature | None`).

---

## 6. How to assess (three concrete checks)

Run with `uv run` and pynauty available preferred (nauty default when importable).

### A. HQ collapse (isomorphic pair ends → one edit)

Hydroquinone `Oc1ccc(O)cc1` + Dehydrogenation: one pair unique-edit, one quinone.

```399:403:tests/refactor_poc/test_pair_signatures.py
def test_hydroquinone_dh_one_signature_one_quinone():
    mol = _mol("Oc1ccc(O)cc1")
    by_sig = _pair_product_csmi(Dehydrogenation(), mol)
    assert len(by_sig) == 1
    assert next(iter(by_sig.values())) == frozenset({"O=C1C=CC(=O)C=C1"})
```

Also: `test_hydroquinone_pair_emits_once_not_per_resonance_parent` in
`test_dh_bond_atom_unique_edit.py` (resonance parents must not double-emit).

### B. Phenol C–O+O ≠ C–O+C (direction is chemistry)

Same C–O bond; site atom O vs site atom C → distinct `BondAtomOrbitSignature`.
If these collapsed, direction would not be part of the edit class.

```python
from xenosite.refactor_poc.rdkitutil import MolFromSmiles
from xenosite.refactor_poc.graph_isomorphism import bond_atom_orbit_key

mol = MolFromSmiles("Oc1ccccc1")
o = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
c = next(a.GetIdx() for a in mol.GetAtomWithIdx(o).GetNeighbors())
b = mol.GetBondBetweenAtoms(o, c).GetIdx()
assert bond_atom_orbit_key(mol, b, o) != bond_atom_orbit_key(mol, b, c)
```

(At HEAD: groups `(0,0)` vs `(0,6)`; both may use `TRIVIAL_PAIR_GROUP` — see caveats.)

### C. Ethane either C same (isomorphic carbons)

Ethane `CC`: C–C bond + carbon 0 equals C–C bond + carbon 1 (automorphism).

```python
eth = MolFromSmiles("CC")
assert bond_atom_orbit_key(eth, 0, 0) == bond_atom_orbit_key(eth, 0, 1)
assert eth.xf.topol_equiv[0] == eth.xf.topol_equiv[1]
```

Related same-kind check: benzene meta share orbit, meta ≠ para
(`test_benzene_meta_and_para_share_ranks_but_not_atom_pair_orbits`).

Ordered/unordered wiring: `test_pair_signatures.py`,
`test_dh_bond_atom_unique_edit.py`, `test_nauty_six_family_orbits.py`.

---

## 7. Where code lives

| Location | Responsibility |
|----------|----------------|
| `src/xenosite/refactor_poc/graph_isomorphism.py` | Nauty/isotope orbits, tables, `TRIVIAL_PAIR_GROUP`, `resolved_swap_group` / `ends_swappable`, constructors, `site_signature` / `pair_site_signature` |
| `src/xenosite/refactor_poc/records.py` | `PatternInfo.swap_group`, `When.swap_group`, signature NamedTuples / NewTypes |
| `src/xenosite/refactor_poc/rules.py` | Thin call sites: pass `unique_orbit=getattr(self, "unique_orbit", "atom_atom")` into `site_signature` / `pair_site_signature`; DH sets `unique_orbit = "bond_atom"` |
| `tests/refactor_poc/test_graph_isomorphism.py` | Orbits, tables, TRIVIAL, oracle |
| `tests/refactor_poc/test_pair_signatures.py` | Signature ↔ product contract, order invariance |
| `tests/refactor_poc/test_dh_bond_atom_unique_edit.py` | DH bond_atom + ordered/unordered ends |
| `tests/refactor_poc/test_nauty_six_family_orbits.py` | Six families; atom_bond ordered≡unordered |

Rule call-site pattern:

```1187:1196:src/xenosite/refactor_poc/rules.py
                    signature = site_signature(
                        context,
                        work,
                        mapped,
                        ranks,
                        site,
                        rxn_num,
                        effect,
                        unique_orbit=getattr(self, "unique_orbit", "atom_atom"),
                    )
```

```1924:1935:src/xenosite/refactor_poc/rules.py
                    signature = pair_site_signature(
                        mol,
                        ranks,
                        map1,
                        map2,
                        site_a,
                        site_b,
                        info1,
                        info2,
                        preview,
                        unique_orbit=getattr(self, "unique_orbit", "atom_atom"),
                    )
```

Public xf: `mol.xf.atom_pair_orbit_key` / `bond_pair_orbit_key` /
`bond_atom_orbit_key` / `site_pair_orbits` / `pair_orbit_backend`.
Env: `XENOSITE_PAIR_ORBIT_BACKEND` ∈ `{nauty, smiles, none}`.

Do **not** invent a Site field for bond-vs-atom direction. Do **not** replace
`swap_group` with SMARTS narrowing for ResonancePair same-role couples
(HEURISTICS: not approved).

---

## 8. Open caveats

### `TRIVIAL_PAIR_GROUP` for singletons

```127:130:src/xenosite/refactor_poc/graph_isomorphism.py
# Documented trivial pair_group when either end is a singleton topeqiv group,
# or when backend is ``none`` for multi–multi (cheap path; no isotope tables).
# Reserved negative id so sequential orbit ids stay 0..n-1.
TRIVIAL_PAIR_GROUP: Final[PairGroupId] = PairGroupId(-1)
```

If either end's topo class has size 1, `pair_group = -1` and nauty/smiles
tables are not materialized for that query. Distinct chemistry can still
separate via `groups` (phenol C–O+O vs C–O+C: different atom groups). Backend
`none` also forces trivial for multi–multi (cheap path).

### `atom_bond` ordered ≡ unordered

Nauty emits one partition for both family names. Atoms and bonds are not
interchangeable same-kind sites; ordering does not split groups. Do not invent
a split. Legacy `PairMode "bond_atom"` stores flipped `(bond, atom)` pairs.

```28:37:tests/refactor_poc/test_nauty_six_family_orbits.py
def test_atom_bond_ordered_equals_unordered_partition():
    """Types already distinguish atom vs bond — ordering does not split groups."""
    ...
    assert families["atom_bond_unordered"] == families["atom_bond_ordered"]
```

### Other notes

- Isotope (`smiles` backend) remains callable for profiling; nauty is unique-edit
  authority when pynauty is importable.
- Narrowing SMARTS ≠ `swap_group` replacement for pair same-role couples.
- Map-rank embeddings stay in `pair_site_signature` (QF dealkylate methyl vs
  benzyl). Formula bag key makes `HCl` ≡ `ClH`.
