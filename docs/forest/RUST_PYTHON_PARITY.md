# Rust ↔ Python leaf-rule product parity

Tracking sheet for chematic (Rust) vs RDKit (Python) leaf parity. Product
identity is **RDKit CSMI** after round-trip. Site identity is **topological**
(engine-local ranks / pair orbits) — not raw atom indexes across engines.

Do not put exception name lists in test bodies. Unpaired or intentionally
divergent leaves use **rule attributes**:

- Python: ``ReactionRule.rust_parity_exception`` (non-empty reason string)
- Rust: ``RuleSet::parity_exception`` / ``with_parity_exception``

Corpus completeness (every PatternInfo possibility / every ``when``):
``tests/forest/test_rule_parity_corpus.py`` + ``PARITY_CORPUS``
(``ParityEntry`` = mol + ``CoverIntent`` correspondences).

Pairing completeness: ``tests/forest/test_rule_parity_pairs.py``.

Product/site parity: ``tests/forest/test_rule_parity_fuzz.py`` (parametric
full rule×mol cartesian by default; focused CoverIntent-only with
``pytest --parity-focused`` / ``XENOSITE_PARITY_FULL=0`` — CI may disable).

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

``PARITY_CORPUS`` / ``PARITY_FUZZ_MOLS`` is owned by
``test_rule_parity_corpus``. Missing cover is a failure (no xfail). Grow the
list; do not drop inventory rows. Each inventory possibility has a designated
``CoverIntent`` on a corpus mol (verified by meta-tests).
Status: **approved**.

### C5 — Full rule×mol cartesian by default (CI keeps it when fast)

Default parametric parity / CSMI uses full ``paired/leaf × corpus`` cartesian
via ``parity_param_cases`` — ``pytest.mark.parametrize``, not Hypothesis.
Focused CoverIntent-only remains available via ``XENOSITE_PARITY_FULL=0`` /
``pytest --parity-focused`` for local speed experiments; CI should keep full
on while the suite stays fast enough (~seconds). Hypothesis sampling over
that grid is redundant. Status: **approved**.

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
a meta check that the corpus still hits). Cases use full rule×mol
cartesian by default; focused CoverIntent-only via ``--parity-focused`` /
``XENOSITE_PARITY_FULL=0``. Other suites do **not** xfail non-compliant rules.

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
or where pattern data over-emits.

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

#### Options under consideration (historical)

1. Canonical site + broad ``product_equiv`` on same product bag alone.
2. SMARTS / ``when`` partition (ester carbonyl-side only).
3. PatternInfo flag folding directed embeddings that share map 2.
4. Rely on ``unique_csmi`` — **not approved** (C11).

**Superseded by C18:** adopted criterion is **equivalent-embedding collapse**
on shared undirected scissile bond + same linker-cleavage family + same
normalized product bag — **not** broad same-product fold, **not** SMARTS
partition as the primary design, **not** ``unique_csmi``. Status of this
entry: **superseded**. Tracker remains work #15 until schema lands.

**Reproduce / live bags:** [PARITY_REMAINING.md](PARITY_REMAINING.md).

### C14 — ``exclusive_partner`` on Effect (bridging N/O)

ResonancePair ends that consume a heteroatom partner (iminium, O/N
``single_to_double``, dealkylate, replace_halogen) set
``Effect.exclusive_partner``. The couple loop refuses when an exclusive
partner atom appears on the other end's map. Chemistry decides the bit —
not a global shared-map disallow. Methide alkyl (``partner=="C"``) does
**not** set it; shared CH2-bridge methide couples remain open (work #19).

Same shape as ``methide`` / ``skip_same_rings``: data on the record, generic
gate reads it. Status: **approved** (schema + gate). Tracker: work #18.

### C16 — Rust-ahead QF products are the chemical guide (do not regress)

On phenyl isocyanate / isothiocyanate / carbodiimide
(``O=C=Nc1ccccc1``, ``S=C=Nc1ccccc1``, ``N=C=Nc1ccccc1``), Rust
``QuinoneFormation`` emits the shared ring-oxidized N=C=X products **and**
the dealkylate cleavage pair (``C=O``/``C=S``/``C=N`` leave + quinone-imine
``N=C1C=CC(=O)C=C1`` / ``N=C1C=CC=CC1=O``). Python emits only the
ring-oxidized set (parity ``rust_extra``).

Those Rust extras are chemically the QF dealkylate door on the exocyclic
``N–C(=X)`` bond — not junk sanitize survivors. **Closing product parity by
dropping Rust emissions is not approved** (C8: chemical correctness first).
Prefer lifting Python to the Rust bag, or a documented Python miss; pin Rust
with regression tests. Status: **approved**.

### C17 — Remaining leaf-parity fails are **not** RDKit↔Chematic arom-model

Checked (2026-09-27) the open parametric fail substrates (olsalazine, cinnoline,
nitro/nitroso, PhNCO/S/N, sulfamethoxazole, benzene-oxide, dialkylanilines,
aziridine-Ph, dihydroacridine, aspirin, chloramphenicol, …):

| Check | Result |
|-------|--------|
| Forest `parse_mol` / `aromatize` | already `apply_aromaticity_rdkit_parity_experimental` (`mol.rs`) |
| Chematic Hückel vs `RdkitLike` vs experimental parity | **identical** aromatic atom/bond sets on every fail substrate |
| Those sets vs live RDKit `GetIsAromatic` | **identical** (pinned in `mol.rs` test `remaining_parity_substrates_match_rdkit_aromatic_atoms`) |

**None** of the current leaf product-parity fails are marked arom-model.
Switching engines / “turning on” RDKit compat cannot close them — compat is
already the production path. Chematic SMARTS also exposes
``RdkitParityConfig::use_rdkit_parity_aromaticity`` (default ``false``); forest
does not need it because ``parse_mol`` / ``aromatize`` already stamp the
parity flags before ``find_matches``. Treat remaining gaps as edit / accept /
form / site-bag issues (H cumulated paths, QF extras, Dealk ring-open, etc.).
Status: **approved**.

### C18 — Chemistry-first parity plan (adopt / modify / defer)

Decision on the 2026-09-27 chemistry-first resolution plan. Full fail table:
[PARITY_REMAINING.md](PARITY_REMAINING.md). Does **not** target literal Python
parity (C8). Python leaf bags are historical evidence; approved chemistry is
the oracle. RDKit CSMI remains the **identity** key (C1) — that is spelling,
not chemistry authority.

#### Adopted posture — **approved**

1. Preserve chemically plausible Rust-ahead products.
2. Remove Rust products that violate transformation family or valence/redox.
3. Port Python behavior only when it is real intended metabolism.
4. Unresolved differences are explicit policy — not silent ``unique_csmi`` /
   CSMI churn (C11).

Conceptual split (implement with C13 schema, not all at once):
``raw_match_count`` / ``event_count`` / ``product_count``. Do **not** invent
the full product-metadata YAML surface until each field has a home on
``PatternInfo`` / ``Emission`` / SiteInfo (data-not-branches).

#### Phase 1 immediate fixes — **approved** (implement next)

| # | Case | Decision |
|---|------|----------|
| 1 | Nitro charge form (Dealk/NDealk) | Keep Rust chemistry; normalize to ``[N+](=O)[O-]``. Not a new reaction. |
| 2 | Cinnoline ring-open N=N | Keep Rust azo/imine; rematch must not add undeclared hetero H / charge. |
| 3 | Generic S=O hydrogenation | Refuse: generic H may not select hypervalent S/P π edges unless a dedicated rule opts in (data gate, not molecule ``if``). |
| 4 | SMx C–N cleavage | Rust-ahead accepted; pin or lift Python. |
| 5 | Cumulated ``X=CN=C1…`` H | Keep Rust refuse; quarantine Python extras. |
| 6 | Ph2NMe QF | Enumerate both carbon leaves (Me and Ph) when valence/accept pass. |
| 7 | Dihydroacridine → acridine | Keep under **Dehydrogenation**; not QF (see family contract). |

#### Policy layer — **approved** (write before patching)

**C13 event equivalence (supersedes C13 options):** collapse directed
embeddings only when they share the **same undirected scissile bond**, belong
to the **same declared linker-cleavage family**, and emit the **same
normalized product bag**. Emit one canonical event; retain equivalent
embeddings on the record (name TBD — not a broad ``product_equiv`` keyed
only on products). SMARTS partition is **not** the primary design.
``unique_csmi`` is **not** the mechanism.

Aspirin ``CC(=O)Oc1ccccc1C(=O)O``: one ester event, two directed embeddings
``(1,3)`` / ``(4,3)``.

Benzene-oxide DH ``C1=CC2OC2C=C1``: same architecture, but first compare
atom/bond-change sets. Product ``c1ccc2c(c1)O2`` is **aromatic benzene-oxide
framework writing**, not benzofuran (``c1ccc2occc2c1``). Fix naming in notes.

**Nitro hydrogenation:** remove from generic Hydrogenation; dedicated
NitrogenReduction (nitroso / hydroxylamine / amine). No radicals; no
preferred endpoint ``[O-][NH+](O)Ar``. Applies to nitrobenzene and
chloramphenicol.

**QF vs DH:** DH = net H loss / aromatization; QF = quinonoid carbonyl (or
explicitly permitted quinone-imine/methide patterns). Pure aromatization is
not QF. Dihydroacridine→acridine = DH only.

**Isoxazole N–O (SMx NitrogenReduction):** if atom-mapped tautomer
equivalents, canonicalize to keto–enamine
``CC(=O)C=C(N)NS(=O)(=O)c1ccc(N)cc1``; if connectivity differs, treat as
template difference — do not normalize away.

#### Deferred hard cases — **approved** (boundaries now; no local parity patch)

| Case | SMILES | Boundary now |
|------|--------|--------------|
| Olsalazine H / QF / DH | ``OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`` | No generic path across azo bridge; dedicated azo reduction later; quarantine Py multi-ring path extras |
| Cinnoline H / QF multi-ring | ``c1ccc2nnccc2c1`` | Explicit center + H delta + bounded component before enabling; keep ring-open N=N fix (Phase 1) separate |
| Benzene-oxide QF | ``C1=CC2OC2C=C1`` | Exclude mere re-aromatization / non-quinonoid mono-carbonyl until QF contract is code |

#### Not approved / modified

| Proposal | Status |
|----------|--------|
| Broad ``product_equiv`` = same product bag alone | **not approved** (too wide) |
| SMARTS orientation as primary C13 fix | **not approved** as primary (ok as later pattern-local aid) |
| ``unique_csmi`` as C13 design | **not approved** (C11) |
| Literal Python parity as target | **not approved** (C8) |
| Full product metadata schema in Phase 1 | **not decided** field-by-field; three counts are the goal |
| Retire Python tests immediately | **not approved**; demote as oracle via Rust-ahead pins + approved chemistry |

#### Parametric xfail + Rust standardization — **approved**

Open product-parity gaps are annotated on the corpus mol as
``ProductParityXfail(rule, reason)`` on ``ParityEntry.product_xfails``
(``tests/forest/rule_parity_corpus.py``). The parametric suite runs the
assertion and calls ``pytest.xfail(reason)`` when that pair fails — reasons
live with the example, not a side table. Meta-tests refuse stale PASS
annotations and unmarked live gaps. Remove the annotation when Phase 1–3
closes it against **approved chemistry**.

**Product / charge / tautomer standardization** (nitro ``[N+](=O)[O-]``,
isoxazole keto–enamine, etc.) must be implemented as **explicit Rust
transforms** with named invariants and unit tests. **Not approved:** RDKit
``SanitizeMol`` / charge-move black boxes, Chematic opaque reparse as the
normalizer, or silent ``unique_csmi`` as chemistry. Display identity may
still round-trip RDKit CSMI (C1); the chemical rewrite itself stays in
forest Rust.

**Milestone tag:** ``py-rust-parity-c18-baseline`` — parametric leaf product
parity green with ``C18_OPEN_PRODUCT_XFAIL_COUNT`` (29) known gaps annotated;
chemistry-first plan (this entry) recorded; arom-model ruled out (C17).
Decrement the constant when closing a row.

#### Work order impact

- Phase 1 → drive work #12 toward green on approved chemistry (not blind match).
- Work #15 → implement C18 C13 criterion (schema + aspirin/benzene-oxide tests).
- Hard cases stay open with instrumentation, not conjugated-path rewrite.

Status: **approved**.

---


### C15 — Specialize SMIRKS for chematic apply (nested, bonds, H0)

``specialize_smirks_for_maps`` is the shared apply door for Python/RDKit
SMIRKS primitives chematic rejects: nested recursive ``$()`` brackets,
bond or-queries (``-,:``, ``=,:``) resolved from the live match, and
product ``[*&H0&+:map]`` → ``[ElH0+:map]`` with post-apply H0 enforcement
when chematic saturates ``[SH+]``. Hydroxy S-ox product is organic ``O``
(not radical ``[O]``) on both catalogs. Status: **approved**.

## Work order (tackle in this sequence)

| # | Issue | Status | Notes |
|---|--------|--------|-------|
| 0 | Corpus meta-test covers every pattern / every ``when`` | **done** | C4 |
| 1 | Discover all Rust + Python leaves; pair or attribute-exception | **done** | C3, C6 |
| 2 | Expose rulesets + ``find_path`` + enumerate to Python | **done** | C9; doors wired |
| 3 | Parametric parity (full cartesian default; focused via ``--parity-focused``) | **done** | C5 |
| 4 | **Pair / dual-door unique-edit: partition SMIRKS↔pair overlap** | **open** | C12 — fix data/unique-edit (not no-dedup); H alkene+path_end |
| 5 | Acetylation: chematic SMIRKS apply miss | **done** | Product keeps ``[#6](=[#8])[#6]``; Rust ``organic_product_variants`` expands ``#`` to organic aliphatic then aromatic (bracket expand misses). Issue: ``CHEMATIC_ISSUE_acetyl_atomic_product.md`` |
| 6 | Dephosphorylation: Rust emits nothing on ``COP(=O)(O)O`` | **done** | specialize nested ``$()``; phosphate ester apply green |
| 7 | AzoSplitting / ThiopheneSulfurOxidation: Rust empty on aromatic examples | **done** | specialize bond or-queries ``=,:`` / ``-,:``; charge wildcards |
| 8 | NitrogenReduction: Rust empty on ``CCNO`` (hydroxylamine) | **done** | specialize ``-,:`` → live bond; hydroxylamine apply green |
| 9 | BenzodioxoleReduction: wrong products | **done** | CH2-leave atom-remove edit (chematic SMIRKS drops ring bonds) |
| 10 | SulfurOxidation: form / set mismatch on ``CCS`` | **done** | hydroxy ``[O]``→``O``; enforce H0 on charged S-oxide |
| 11 | Dehydration: site-count mismatch with matching products | **done** | beta_elim site_map=1 matches Python |
| 12 | Drive parametric suite green on **approved chemistry** (C18), not literal Py match | **in progress** | C8, C18 Phase 1 |
| 13 | Reduce reliance on Python sanitize for product validity | **open** | C10 — soft failure; prefer emit-path fixes |
| 14 | ``unique_csmi_compliant`` + CSMI-dup parametric test | **done** | C11 — only that test xfails non-compliant; hard-fail if compliant |
| 15 | C13 equivalent-embedding (scissile bond + family + bag); aspirin / benzene-oxide | **open** | C18 supersedes C13 options; schema TBD |
| 22 | C18 Phase 1: nitro form, cinnoline N=N, S=O H refuse, SMx C–N pin, cumulated H refuse, Ph2NMe both leaves, DH acridine | **open** | C18 |
| 23 | C18 Phase 2: nitro≠generic H; QF↔DH contract; isoxazole tautomer canon | **open** | C18 |
| 24 | C18 Phase 3: olsalazine/cinnoline/benzene-oxide QF — defer + instrument | **deferred** | C18 |
| 16 | Chematic atom-tracking migrate; drop vendored patch | **blocked on parity** | TODO.md § After parity — pass tracking tests, remove `vendor/chematic` |
| 17 | Centralize ForestMol cache + copy/edit behind ForestMol methods | **blocked on parity** | TODO.md § After parity — private-to-find-strays; keep user-visible cache |
| 18 | ``exclusive_partner`` bridging N/O: tests + corpus + Rust↔Py parity | **done** | C14 — Python ``test_exclusive_partner`` + Rust ``pair_edit`` refuse tests; corpus bridging N/O extras |
| 19 | Shared methide alkyl partner (same CH2 bridge) | **open** | C14 deliberately leaves ``partner=="C"`` off; decide if chemistry needs its own Effect bit |
| 20 | Pair-close / identical-partner corpus + meta-test green | **done** | ortho catechol/diamine/…; ``test_parity_fuzz_mols_cover_close_pair_ends`` green |
| 21 | Non-enzymatic rearrangements ruleset (find default, not enumerate) | **blocked on parity** | TODO.md § After parity #3 — TautomerRule + AzoSplitting + ring-closure; find variants include; enumerate excludes |

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