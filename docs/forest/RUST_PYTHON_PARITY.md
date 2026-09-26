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
Status: **not decided**. Superseded by C12 (investigation + decision).

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

### C10 — Python sanitize cleanup of bad metabolites is a soft failure

Live Python often emits chemically bad fragments that only become acceptable
(or disappear) after RDKit ``SanitizeMol`` / ``sanitized_fragments`` /
``xf.sanitize``. That post-hoc cleanup is a **soft failure**: the edit path
produced junk that sanitize papered over. Ideally patterns, valence gates, and
apply would not need sanitize to rescue products (Rust already refuses some
shapes early, e.g. two-double nitrogen — see HEURISTICS).

Parity implications: comparing only post-sanitize RDKit CSMI can hide that
Python relied on sanitize while Rust never emitted (or emitted a different
raw graph). When investigating empty-Rust / extra-Python gaps, check whether
Python’s survivor required sanitize. Prefer fixing the emit path over
widening sanitize. Status: **approved** (goal); tightening emit so sanitize
is unnecessary remains open work, not an excuse for product mismatch.

### C11 — ``unique_csmi`` product collapse is a soft failure

Yield-layer CSMI dedup is a **soft failure**, same family as C10 sanitize:
sites / patterns / ``when`` should already have collapsed equivalent edits.
Product identity is not the design.

**Rule data:** ``ReactionRule.unique_csmi_compliant`` (default ``True``).
Yield CSMI drop runs only when the caller passes ``unique_csmi=True`` **and**
the rule is compliant. Non-compliant leaves emit every unique-edit survivor
(no yield papering). Goal: every leaf ``True``.

**Tests:** only ``tests/forest/test_csmi_dedup_soft_failure.py`` gates this —
compliant → hard-fail on duplicate ``(rule, pattern, CSMI)`` keys or
``SiteDeduplicationWarning``; non-compliant → ``pytest.xfail`` on hits (and
a meta check that the corpus still hits). Other suites do **not** xfail
non-compliant rules.

``SiteDeduplicationWarning`` remains the unique-edit-miss check signal.
Quiet unequal-rank same-CSMI pairs are C13 until ``product_equiv`` / identity
data lands. Status: **approved** (goal).

### C12 — C7 resolution: fix partition / unique-edit, not “no dedup”

Investigation of the C7 citation (Hydrogenation ``C=C`` py 2 vs rs 1):

- Pair ``pair_site_signature`` already keeps one ResonancePair survivor on
  ethene. The second Python emission is the **alkene SMIRKS** door on the
  same unordered atom-pair site (dual door, not pynauty↔canonaut orbit
  disagreement). Parametric site bags (C2) already collapse that to one
  topo class — site-count parity is green there.
- ``unique_csmi`` does **not** collapse those two (distinct pattern tokens:
  ``alkene`` vs pair ``None``). Relying on CSMI to paper that over would be
  C11; the real fix is **data partition** (or shared site unique-edit across
  SMIRKS + pair doors) so sites/patterns/whens alone yield one edit — same
  style as Hydroxylation ``h``/``h2`` and Dealkylation H0/h SMARTS splits.
- Option (b) from C7 (parity without unique-edit) is **not approved**: it
  hides the miss. Option (a) is **approved** as: fix unique-edit and/or
  partition overlapping SMIRKS↔pair coverage; do not use yield CSMI as the
  dedup. Longer conjugated-path Hydrogenation gaps vs Rust also track
  incomplete Rust ``conjugated_component`` (chemistry work order), not C7
  alone.

Status: **approved** (supersedes C7’s undecided fork).

### C13 — Product-identical distinct sites; Dealk one-side polymorphism

Two related gaps where ``unique_csmi`` quietly collapses (C11 soft failure)
or where pattern data over-emits. **Not decided** — need a schema/data
answer before code. Do not invent a silent branch.

#### A. Genuinely two sites → one product bag

Example (corpus): Dealkylation ``quaternary_alcohol`` on aspirin
``CC(=O)Oc1ccccc1C(=O)O``:

- directed site ``(1, 3)`` ranks ``(8, 9)`` — acetyl C–O
- directed site ``(4, 3)`` ranks ``(8, 11)`` — aryl C–O
- same emission CSMI frozenset ``{CC(=O)O, O=C(O)c1ccccc1O}``

Unequal ranks ⇒ not a unique-edit miss (no ``SiteDeduplicationWarning``);
yield CSMI still drops one. Same shape may appear for symmetric epoxide
hydration / opening when both carbons are chemically distinct embeddings
but products coincide (symmetric cases often already collapse via
``topol_equiv`` / unique-edit ``orbit`` — HEURISTICS pass-down orbits).

Automorphism-equivalent sites are **already** handled: one canonical emit +
``SiteInfo.orbit`` / ``Emission.site_orbit`` (HEURISTICS approved). This
gap is the **other** class: topologically distinct sites, identical
products.

#### B. Dealk map-1 polymorphism vs the other side

Dealkylation encodes carbon-side polymorphism as separate ``PatternInfo``
rows (H-count × alcohol/carbonyl/carboxylic) with ``when`` on the partner
(map 2 = N/O/S). That is good data. The failure mode: both carbons on a
bridging heteroatom can play map 1 (ester O bonded to two ``#6H0``), so the
**other** side of the same linker double-emits under the same pattern name
even when the chemistry outcome is one cleavage. Methyl vs aryl on
``Ph–O–Me`` correctly stay distinct (different patterns / products); the
bad case is both sides quaternary (ester) with one product bag.

``directed_bond`` must stay for regioisomers that **differ** in product
(DIVERGENCES / HEURISTICS). Undirected unique-edit for all Dealk is
**not** a candidate.

#### Options under consideration (pick later; do not ship a branch)

1. **Canonical site + ``product_equiv`` on SiteInfo** (user sketch): after
   apply, under ``(rule, pattern)``, group emissions with equal fragment
   CSMI frozenset; emit once at lex-smallest site key; record the other
   sites on a new field (name TBD: ``product_equiv`` / ``equiv_sites``).
   Distinct from unique-edit ``orbit`` (automorphism class). Discovery
   still uses product identity; yield is explicit one + data, not silent
   ``unique_csmi`` drop. Epoxide / ester both fit if they land in this
   class.

2. **SMARTS / ``when`` partition** so only one side of a Y-linker matches
   for product-symmetric cleavages (e.g. ester carbonyl-side only). Pure
   sites/patterns/whens; no product-layer fold. Risk: missing a real
   regioisomer when products would have differed; must be pattern-local.

3. **PatternInfo flag** that directed embeddings sharing map 2 (leave
   partner) are product-folded when products match — data on the rule,
   still needs apply or a proven identity. Close to (1) with a narrower
   key.

4. **Rely on ``unique_csmi``** for this class only. **Not approved** as the
   design (C11); at most a temporary xfail until (1)–(3) lands.

Related: ``cleave_side_group`` folds **across rules** at expand; it does
not solve within-leaf directed double-emit. ``canonical_emitted_sites`` is
automorphism lex remap, not product-class fold.

Status: **not decided**. Tracker: work order #15. Soft-fail xfail for
aspirin Dealkylation stays under C11 until this choice is approved and
implemented.

### C14 — ``exclusive_partner`` on Effect (bridging N/O)

ResonancePair ends that consume a heteroatom partner (iminium, O/N
``single_to_double``, dealkylate, replace_halogen) set
``Effect.exclusive_partner``. The couple loop refuses when an exclusive
partner atom appears on the other end's map. Chemistry decides the bit —
not a global shared-map disallow. Methide alkyl (``partner=="C"``) does
**not** set it; shared CH2-bridge methide couples remain open (work #19).

Same shape as ``methide`` / ``skip_same_rings``: data on the record, generic
gate reads it. Status: **approved** (schema + gate). Tracker: work #18.

---

## Work order (tackle in this sequence)

| # | Issue | Status | Notes |
|---|--------|--------|-------|
| 0 | Corpus meta-test covers every pattern / every ``when`` | **done** | C4 |
| 1 | Discover all Rust + Python leaves; pair or attribute-exception | **done** | C3, C6 |
| 2 | Expose rulesets + ``find_path`` + enumerate to Python | **done** | C9; doors wired |
| 3 | Parametric parity (bounded corpus × paired leaves), not Hypothesis | **done** | C5 |
| 4 | **Pair / dual-door unique-edit: partition SMIRKS↔pair overlap** | **open** | C12 — fix data/unique-edit (not no-dedup); H alkene+path_end |
| 5 | Acetylation: chematic SMIRKS apply miss | **done** | Product ``[*:1]C(=O)C`` (organic); ``[#6](=[#8])[#6]`` expands to bracket forms that miss. Reactant keeps ``#``; tests cover aliphatic + aromatic heteroatom specialize branches. Issue draft: ``CHEMATIC_ISSUE_acetyl_atomic_product.md`` |
| 6 | Dephosphorylation: Rust emits nothing on ``COP(=O)(O)O`` | **open** | recursive SMARTS / P valence? |
| 7 | AzoSplitting / ThiopheneSulfurOxidation: Rust empty on aromatic examples | **open** | ResonanceRule ``=,:`` / Kekulé |
| 8 | NitrogenReduction: Rust empty on ``CCNO`` (hydroxylamine) | **open** | Resonance / pattern arm |
| 9 | BenzodioxoleReduction: wrong products | **open** | catechol vs ring-opened; C8; watch C10 |
| 10 | SulfurOxidation: form / set mismatch on ``CCS`` | **open** | ``CCSO``+zwitterion vs ``CCS[O]``; watch C10 |
| 11 | Dehydration: site-count mismatch with matching products | **open** | site_map / topo key align |
| 12 | Drive parametric suite green; no silent skips | **blocked on 4–11** | only C6-style excuses |
| 13 | Reduce reliance on Python sanitize for product validity | **open** | C10 — soft failure; prefer emit-path fixes |
| 14 | ``unique_csmi_compliant`` + CSMI-dup parametric test | **done** | C11 — only that test xfails non-compliant; hard-fail if compliant |
| 15 | Product-identical distinct sites + Dealk other-side double-emit | **open** | C13 — not decided; options: ``product_equiv`` field vs SMARTS partition |
| 16 | Chematic atom-tracking migrate; drop vendored patch | **blocked on parity** | TODO.md § After parity — pass tracking tests, remove `vendor/chematic` |
| 17 | Centralize ForestMol cache + copy/edit behind ForestMol methods | **blocked on parity** | TODO.md § After parity — private-to-find-strays; keep user-visible cache |
| 18 | ``exclusive_partner`` bridging N/O: tests + corpus + Rust↔Py parity | **done** | C14 — Python ``test_exclusive_partner`` + Rust ``pair_edit`` refuse tests; corpus bridging N/O extras |
| 19 | Shared methide alkyl partner (same CH2 bridge) | **open** | C14 deliberately leaves ``partner=="C"`` off; decide if chemistry needs its own Effect bit |
| 20 | Pair-close / identical-partner corpus + meta-test green | **done** | ortho catechol/diamine/…; ``test_parity_fuzz_mols_cover_close_pair_ends`` green |

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
5. **Progress** → work-order Status is the gauge (``open`` / ``in progress`` /
   ``done`` / ``blocked on …``). Do not delete rows; flip Status.