# Pair-orbit unique-edit (agent briefing)

Source of truth at HEAD: `graph_isomorphism.py` + signature types in `records.py`.
Rule chemistry stays in `rules.py` (thin call sites). HEURISTICS Status: approved
for `swap_group` / nauty six-family / DH `unique_orbit="bond_atom"`.

This doc is for another agent. Facts only. Do **not** claim novelty for
“automorphisms on pairs” — that construction is classical (see References).

---

## Methods (citation framing)

We computed automorphism generators of labeled molecular graphs using
nauty/Traces and derived equivalence classes as induced group orbits on typed
site sets, ordered tuples, and unordered subsets. Orbits on ordered pairs are
classically termed *orbitals*; we extended this construction recursively to
pairs of directed bond–atom sites.

Cite McKay–Piperno and Sharp after that sentence. Chemistry-facing novelty (if
any) is typed sites, directed bond–atom keys, recursive composite-site
orbiting, and stable canonical emission — not the group action on pairs.

---

## 1. Problem unique-edit solves

**Collapse isomorphic edits; do not collapse distinct chemistry.**

Two layers (do not conflate — HEURISTICS):

1. **Site unique-edit** — topological ranks + pair-orbit signature. Equivalent
   embeddings of the *same* edit share one key before `RunReactants`.
2. **Product csmi dedup** (`unique_csmi`) — drops site topology; keys
   `(rule name, PatternInfo.name | SMARTS, product csmi)`.

Pair orbits refine (1) when atom ranks alone are not enough. Classic negative:
benzene meta vs para share ranks but not atom–atom orbits
(`test_benzene_meta_and_para_share_ranks_but_not_atom_pair_orbits`).

Core signature shape (not a `Site`): `((ga, gb), pair_group_id)`.

---

## 2. Same-kind pairs (atom–atom / bond–bond)

Ordered vs unordered is decided by resolved **`swap_group`**
(When → PatternInfo → `name`):

- Same non-empty resolved group on both ends → unordered (`ordered=False`,
  `end_ranks=()`).
- Unequal → ordered (`ordered=True`, `end_ranks` in `PatternInfo.name` order).

Omit `swap_group` when it would equal `name`. Set it only when grouping
differs. Helpers: `resolved_swap_group`, `ends_swappable`.

Nauty families (source of truth: `all_site_pair_orbits_nauty`):

| Family | Candidates |
|--------|------------|
| `atom_atom_unordered` | combinations `(a1,a2)` with `a1 < a2` |
| `atom_atom_ordered` | permutations `(a1,a2)` with `a1 != a2` |
| `bond_bond_unordered` | combinations of bonds |
| `bond_bond_ordered` | permutations of bonds |
| `atom_bond_*` | one partition for both names (types already differ) |

Examples: DH `phenol_end`×`phenol_end` unordered; `phenol_end`×`amine_end`
ordered. QF same-name pairs unordered; cross-role ordered.

---

## 3. Bond–atom: directed, not a Site

Bond role ≠ atom role. A bond–atom unique-edit key is **directed** `(bond, atom)`.
It is **not** a `Site` definition.

One directed unit captures **both** topological site identity (orbit
`groups` + `pair_group`) **and** edit direction (bond vs atom role). Direction
is part of the edit's isomorphism class — not a separate Site field.

Dehydrogenation declares `unique_orbit = "bond_atom"` as data.

### 3b. Second-order orbit (ResonancePair)

A directed composite site is `s = (bond, atom)`. A ResonancePair selects two
such sites `((b1,a1), (b2,a2))`. Individual `BondAtomOrbitSignature` values
identify each site's *first-order* orbit; they do **not** identify the orbit
of the *pair*.

On benzene every endpoint `(bond, atom)` is equivalent, so an ends-only key
collapses adjacent vs opposite placements. Joint orbit lives on
`BondAtomPairOrbitSignature.pair_group`.

API: `endpoint_bond_atom_sites` (closed),
`bond_atom_orbits_from_nauty_generators`,
`bond_atom_pair_orbits_from_nauty_generators(..., ordered=...)`.
Wired through `pair_orbit` / `unordered_bond_atom_pair` /
`ordered_bond_atom_pair`. Regression:
`test_benzene_bond_atom_pair_needs_joint_orbit`.

---

## 4. Types

| Type | Role |
|------|------|
| `AtomPairOrbitSignature` | atom–atom; `ordered` + `end_ranks` |
| `BondPairOrbitSignature` | bond–bond; same shape |
| `BondAtomOrbitSignature` | one directed `(bond, atom)` unit |
| `BondAtomPairOrbitSignature` | two bond–atom units + joint `pair_group` |

Assembly: `site_orbit` / `pair_orbit` → `site_signature` /
`pair_site_signature`. Public xf: `atom_pair_orbit_key` /
`bond_pair_orbit_key` / `bond_atom_orbit_key` / `site_pair_orbits` /
`pair_orbit_backend`. Env: `XENOSITE_PAIR_ORBIT_BACKEND` ∈
`{nauty, smiles, none}`.

Legacy three-mode tables flip cross-kind to `(bond, atom)` for unique-edit.

---

## 5. Opt-in canonical emitted sites

Default **off**. Enable only via explicit kwargs
``metabolize(..., canonical_emitted_sites=True)`` (also
``find_path`` / ``bfs`` / ``dfs``, which splat down to metabolites).
No env or process-wide toggle.

Lex-smallest concrete member of each nauty orbit
(``LexicalOrbitRepresentatives`` on forest ``cache``); relative to current
RDKit indexing. Kinds: **singleton atom**, atom–atom, bond–bond, directed
bond–atom, pair of bond–atom sites (ordered vs unordered as for unique-edit).
One-atom unique-edit remaps onto the lex-smallest atom in its automorphism
orbit the same way pairs do.

### Contract (do not invert)

1. **`filter_sites` always sees the discovery site** — never the lex
   representative. Filters must **not** assume ``info["site"]`` after emit
   equals what they filtered. Hazard: search / allows / SiteInfo bags that
   assumed identity with emitted sites will be wrong when opt-in is on.
2. Unique-edit signatures are computed on the **discovery** embedding (orbit
   collapse unchanged).
3. After accept, chemistry remaps onto the lex representative
   (``canonicalize_smarts_match`` / ``canonicalize_pair_match``). Emitted
   **`site` = canonical**; optional **`discovered_site`** holds pre-canonical
   discovery indexes when they differ (SiteInfo, TraceAddition, TraceInfo).
4. Product-only: parent mol indexes are not mutated. After remap, product
   forest last ``idx`` layer is restamped (``restamp_product_forest_last_layer``).
   **`of_products` / ``clear_structure`` wipe the product's ``cache``** — lex-rep
   tables are read from the **parent** (``parent=`` on
   ``ensure_lexical_orbit_representatives`` / ``restamp_product_forest_last_layer``),
   never rebuilt against the cleared child.
5. Non-canonical embeddings are not separate products (unique-edit + chemistry
   at the rep).

Escape hatch for callers that need discovery: read ``discovered_site`` when
present. Tests: ``test_canonical_emitted_sites.py`` (includes
parent-cache-after-clear).
Profile (same find_path cases as forest_copy): OFF 3.862s / ON 4.007s ≈ **+3.7%** wall (commit 6f13ca9).

POC fuzz draws ``canonical_emitted_sites`` with Hypothesis ``st.booleans()`` /
``data.draw`` and passes the kwarg through ``find_path`` / ``bfs`` / ``dfs`` /
``metabolize``. Default stays off. Hardcoding ``True`` on individual flaky
tests (topo, guided gold, path hard, unique-edit) is a **future option** only
if stable lex sites prove helpful — not done yet.

Status: approved (opt-in emission policy; HEURISTICS).

---

## 6. How to assess (three checks)

Prefer `uv run` with pynauty (nauty default when importable).

### A. HQ collapse

Hydroquinone `Oc1ccc(O)cc1` + Dehydrogenation → one pair unique-edit, one
quinone (`test_hydroquinone_dh_one_signature_one_quinone`).

### B. Phenol C–O+O ≠ C–O+C

Same C–O bond; site atom O vs C → distinct `BondAtomOrbitSignature`
(`bond_atom_orbit_key`).

### C. Ethane either C same

`CC`: bond + carbon 0 equals bond + carbon 1 under automorphism.

Also: benzene meta share orbit, meta ≠ para; benzene bond–atom pairs need
joint orbit (§3b). Suites: `test_pair_signatures.py`,
`test_dh_bond_atom_unique_edit.py`, `test_nauty_six_family_orbits.py`,
`test_graph_isomorphism.py`.

---

## 7. Caveats (do not overclaim)

- `TRIVIAL_PAIR_GROUP` (`-1`) when either end is a singleton topeqiv class, or
  backend `none` for multi–multi, or second-order when nauty unavailable /
  all involved classes are singletons.
- `atom_bond` ordered ≡ unordered — do not invent a split.
- `topol_equiv` is RDKit `CanonicalRankAtoms(..., breakTies=False)` (CIP
  ranks), **not** proven automorphism orbits. `TRIVIAL_PAIR_GROUP` is safe
  only under the automorphism-orbit reading.
- Stereo: ordinary tetrahedral / E–Z via CIP / bond stereo labels on the
  colored graph; enhanced stereo groups and nonstandard stereo are not fully
  represented. Undirected nauty graph — directed dative bonds out of scope.
- Isotope (`smiles`) backend is for profiling; nauty is unique-edit authority
  when pynauty is importable.
- Narrowing SMARTS ≠ `swap_group` replacement for ResonancePair same-role
  couples (HEURISTICS: not approved).

---

## 8. Code map

| Location | Responsibility |
|----------|----------------|
| `graph_isomorphism.py` | Generators, six families, tables, TRIVIAL, swap helpers, signatures, lex-rep opt-in |
| `records.py` | `swap_group`, signature NamedTuples, `SiteInfo.discovered_site` |
| `rules.py` | Thin call sites; DH `unique_orbit="bond_atom"`; post-filter canonical remap (SMARTS + ResonancePair) |
| `rdkitutil.py` | `restamp_product_forest_last_layer` (product-only) |
| Tests §5 / §6 | Orbits, signatures, DH, six-family, `test_canonical_emitted_sites.py` |

Do **not** invent a Site field for bond-vs-atom direction. Do **not** replace
`swap_group` with SMARTS narrowing for same-role ResonancePair couples.

---

## References

1. **McKay BD, Piperno A (2014).** Practical graph isomorphism, II.
   *Journal of Symbolic Computation* 60:94–112.
   DOI: [10.1016/j.jsc.2013.09.003](https://doi.org/10.1016/j.jsc.2013.09.003).
   **Why:** nauty/Traces automorphism generators for (vertex-colored) graphs —
   our `atom_bond_generators_nauty` entrypoint.

2. **Sharp GR (1999).** Algorithmic recognition of group actions on orbitals.
   *LMS Journal of Computation and Mathematics* 2:1–33.
   DOI: [10.1112/S146115700000005X](https://doi.org/10.1112/s146115700000005x).
   **Why:** *orbital* = orbit of a group acting on ordered pairs of distinct
   points — the classical name for our ordered same-kind / directed constructions.

### Optional context (not load-bearing for unique-edit)

3. **McKay BD, Yirik MA, Steinbeck C (2022).** Surge: a fast open-source
   chemical graph generator. *Journal of Cheminformatics* 14:24.
   DOI: [10.1186/s13321-022-00604-9](https://doi.org/10.1186/s13321-022-00604-9).
   **Why (optional):** cheminformatics use of nauty automorphism groups for
   isomorphism-free chemical graph generation — adjacent, not site unique-edit.

4. **González Laffitte ME, Phan T-L, Stadler PF (2026).** Extension of partial
   atom-to-atom maps: uniqueness and algorithms. *Algorithms for Molecular
   Biology* 21.
   DOI: [10.1186/s13015-026-00300-5](https://doi.org/10.1186/s13015-026-00300-5).
   **Why (optional):** reaction AAM uniqueness / completion in chem
   informatics. Related domain, **not** a citation for molecular
   automorphism orbits on typed edit sites — do not overclaim.
