# Rust ↔ Python leaf-rule product parity

Tracking sheet for chematic (Rust) vs RDKit (Python) leaf parity. Product
identity is **RDKit CSMI** after round-trip. Site identity is **topological**
(engine-local ranks / pair orbits) — not raw atom indexes across engines.

Do not put exception name lists in test bodies. Unpaired or intentionally
divergent leaves use **rule attributes**:

- Python: ``ReactionRule.rust_parity_exception`` (non-empty reason string)
- Rust: ``RuleSet::parity_exception`` / ``with_parity_exception``

Corpus completeness (every PatternInfo possibility / every ``when``):
``tests/forest/test_rule_parity_corpus.py`` + ``PARITY_FUZZ_MOLS``.

Pairing completeness: ``tests/forest/test_rule_parity_pairs.py``.

Product/site parity: ``tests/forest/test_rule_parity_fuzz.py`` (parametric over
``PARITY_FUZZ_MOLS`` × paired leaves).

Native doors: ``xenosite_forest.RuleSet`` (``leaf`` / ``phase_one`` /
``default_ruleset`` / ``all_rules`` / ``catalog_names``), ``find_path``
(``ruleset=`` override), ``bfs`` / ``dfs`` / ``enumerate``. Python wrappers:
``xenosite.forest.find_path_rust``.

---

## Choices (append-only audit log)

Record every methodology or exception choice here **when it is made**. Do not
rewrite history: supersede with a new dated entry. Each entry has
``Status: approved``, ``not approved``, or ``not decided`` (same words as
[HEURISTICS](HEURISTICS.md)). A silent change in code or tests without an
entry here is drift.

### C1 — Product identity key is RDKit CSMI

Both engines' product SMILES are reparsed with RDKit
(``MolToSmiles(..., canonical=True, isomericSmiles=False)``, maps cleared).
Do not string-compare raw chematic display CSMI to RDKit CSMI.
Status: **approved**.

### C2 — Site identity is topological, engine-local

Compare unique-edit / orbit-collapsed site classes within each engine, then
align via product bags. Do not require chematic atom indexes to equal RDKit
indexes. Status: **approved**.

### C3 — Exceptions live on rule attributes, not test allowlists

Unpaired or intentionally divergent leaves set
``rust_parity_exception`` / ``parity_exception`` with a non-empty reason.
``test_rule_parity_pairs`` only discovers gaps; it does not hardcode names.
Status: **approved**.

### C4 — Corpus must cover every pattern possibility and every ``when``

``PARITY_FUZZ_MOLS`` is owned by ``test_rule_parity_corpus``. Missing cover is
a failure (no xfail). Grow the list; do not drop inventory rows.
Status: **approved**.

### C5 — Bounded corpus × paired leaves is parametric, not Hypothesis

Finite ``PARITY_FUZZ_MOLS`` × ``paired_rule_names()`` →
``pytest.mark.parametrize``. Hypothesis sampling over that grid is redundant.
Status: **approved**.

### C6 — Attribute exceptions recorded to date

| Leaf | Side | Reason (attribute text) | Status |
|------|------|-------------------------|--------|
| EpoxideHydration | Rust | Rust-only PhaseOne leaf (epoxide→diol one-hop); Python PhaseOne omits it | **approved** (unpaired) |
| TautomerRule | Python | design stub; not ported to Rust | **approved** (unpaired) |
| ConjugationRule | Python | base conjugation container; leaf ports are Acetylation / Sulfation / Glucuronidation / Glutathionation | **approved** (unpaired base) |

Product-mismatch excuses for paired leaves: **none approved**. Open gaps
stay in the work-order table until fixed or a new choice entry excuses them.

### C7 — Python pair topology dedup is incorrect for pairs

ResonancePair unique-edit (``pair_site_signature`` / pynauty) can disagree with
Rust ``atom_pair_orbit_id`` (canonaut). Observed: Hydrogenation ``C=C`` py 2 vs
rs 1. Options under consideration: (a) fix Python unique-edit, optionally by
calling a Rust orbit/dedup helper exposed to Python; (b) temporarily iterate
parity **without** unique-edit dedup on pair rules until (a) lands.
Status: **not decided** (blocks fair site-count parity for pair leaves).

### C8 — Adjudication when product sets disagree

Prefer **chemical correctness**, then **least ad hoc** (shared SMARTS/effect
data, not a one-off branch). Archive↔live lasting diffs stay in
[DIVERGENCES](DIVERGENCES.md). Rust↔Python lasting diffs stay in this file
(work order + Choices) until fixed or C6-style attribute excuse.
Status: **approved**.

### C9 — Default catalogs on native doors

``find_path`` default catalog: PhaseOne (override via ``ruleset=``).
``bfs`` / ``dfs`` / ``enumerate`` default catalog: ``default_ruleset``
(override via ``ruleset=``). Matches Rust ``*_default`` vs explicit-``&RuleSet``
split. Status: **approved**.

---

## Work order (tackle in this sequence)

| # | Issue | Status | Notes |
|---|--------|--------|-------|
| 0 | Corpus meta-test covers every pattern / every ``when`` | **done** | C4 |
| 1 | Discover all Rust + Python leaves; pair or attribute-exception | **done** | C3, C6 |
| 2 | Expose rulesets + ``find_path`` + enumerate to Python | **done** | C9; doors wired |
| 3 | Parametric parity (bounded corpus × paired leaves), not Hypothesis | **done** | C5 |
| 4 | **Python pair unique-edit / topology dedup is wrong for pairs** | **open** | C7 — decide (a) fix vs (b) no-dedup iteration |
| 5 | Acetylation: chematic SMIRKS apply miss | **open** | ``[*:1][#6](=[#8])[#6]`` no apply; ``[*:1]C(=O)C`` OK |
| 6 | Dephosphorylation: Rust emits nothing on ``COP(=O)(O)O`` | **open** | recursive SMARTS / P valence? |
| 7 | AzoSplitting / ThiopheneSulfurOxidation: Rust empty on aromatic examples | **open** | ResonanceRule ``=,:`` / Kekulé |
| 8 | NitrogenReduction: Rust empty on ``CCNO`` (hydroxylamine) | **open** | Resonance / pattern arm |
| 9 | BenzodioxoleReduction: wrong products | **open** | catechol vs ring-opened; C8 |
| 10 | SulfurOxidation: form / set mismatch on ``CCS`` | **open** | ``CCSO``+zwitterion vs ``CCS[O]`` |
| 11 | Dehydration: site-count mismatch with matching products | **open** | may clear after C7 / site_map align |
| 12 | Drive parametric suite green; no silent skips | **blocked on 4–11** | only C6-style excuses |

Depth-1 PhaseOne product diffs (separate probe, not leaf-only):
``tests/forest/probe_d1_diff.py`` / ``artifacts/d1_*`` — Dealkylation-heavy;
EpoxideHydration covered by C6.

---

## How to update this file

1. **New methodology choice** → append ``### C<n>`` with Status; never edit an
   old Status in place — add ``Superseded by C<n>`` on the old entry and a new
   entry.
2. **New attribute exception** → add a row under C6 (or a new C-entry if the
   policy for exceptions changes) **and** set the attribute on the rule.
3. **Gap fixed** → mark the work-order row **done** and, if a choice was
   involved, append a short C-entry noting the resolution.
4. **Gap excused without a fix** → only via rule attribute + C6 row (C3).
