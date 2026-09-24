# Find-path heuristics and correctness

Consult this file before changing `find_path` behavior, and again when a search miss or a slow pair shows up. Add an entry when a new idea appears. Do not delete entries. Mark one done when the test that names it passes without `xfail`.

Follow `.cursor/rules/data-not-branches.mdc`. Prefer a generic decision that reads `PatternInfo` over a special case in the search.
`leave_count` and `breaks_ring` on `Effect` are that data: the same `cleaves` bit, with the leaving piece named or open, and ring opening on the bond.

## Exceptions

An exception to that rule is allowed only with a good reason. Record it here. Each entry has `Status: approved`, `not approved`, or `not decided`. A branch in the code without an `approved` entry is not an exception. None are approved.

A rule may later carry data that it is less likely, and a filter may read that. The schema is not decided. Do not invent the field, and do not treat the missing field as an exception.

## Migration

Status: done. Previous forest: `src/xenosite/_archive_forest/` (archive
directory on GitHub). Live public API is `xenosite.forest`. Prefer `mol.xf`
over `AtomTracker` — accessor API in `docs/forest/XF.md` (caches into
`_forest`; TypedDicts in `records.py` are unstable layout, not the public
contract). Callers: `docs/forest/MIGRATING_0.7.md`.

## Correctness, not yet done

- Do not drop a child only because `atom_diff` cost is not strictly lower. The `closer` check in `find_path` refuses a sideways step. Some real routes do not move the score down on every hop. Tried equal-cost sideways when the site touched cleavage / dearomatization / oxygen; PhCH2OH and multi-oxidation quinones blew up (`mol_edits` 21→68). Status: not decided — needs a tighter predicate than “touches the diff”.
- Do not walk a step plan that is the same steps in another order. Once a `Deps` has been yielded, a later walk that is only a reordering of those steps is not a new path.

## Speed, not yet done

- Put the isotope-bearing atom first in the reaction SMARTS. `SubstructMatch` starts at query atom 0 and does not reorder, so a common atom first is walked before the isotope is tested.
- **Rust door candidates (site–pattern–parent).** `RuleSet::candidates` yields deferred triples; `Candidate::materialize` runs the edit. `PairCandidate` is the ResonancePair analogue (merged `Effect` for filtering; path flip on materialize). `find_path` filters by reading pattern/effect data before materialize (`mol_edits` only counts survivors). Filter closures on `metabolize` stay for Python parity. Status: approved for the derisk door.
- Deferred metabolites. `metabolites` should be able to yield a light object instead of a sanitized mol. The object carries the anticipated formula (the same counts as `atom_trace["formula"]` and `delta_formula`) and says whether it is sure materialization can succeed. `metabolize` either returns that object or materializes it, so today's callers still get mols. A later hop materializes only when a filter or the plan needs the real mol. The same object is how a `StepPlan` can check that a later site was `added_by` an earlier step without replaying the edit.
- A richer walk-ranking key than target-hit + FIFO. Status: not decided. The frontier is already a lazy heap (see Already in). Sibling cost-sort before enqueue looked like a heap but preferred locally cheap dead ends on PhCH2OH→quinone; do not revive cost as the key without a better one. When a richer key lands, `_fresh_walk_priority` should recompute it from `atom_diff` / `order_key` data at pop.
- Pattern frequency. Count how often each pattern matches, and how often that match shows up in a real metabolite. Deprioritize or filter patterns that match almost everything and are rarely the one that was used. That number belongs on `PatternInfo`, not in a table inside `find_path`.
- Cleavage-first may not need its own branch once a heap prefers cleaving spans when the target is smaller. Today's `order_key` still reads `cleaves` on the span; do not add another cleavage predicate until a richer walk key has been tried.

## Already in the search

- Cleavage patterns run first when the target is smaller or a mapped bond is missing. Speed, not a new rule class. `order_key` also prefers dearomatizing spans when `loses_aromaticity`, and oxygen-adding spans when `needs_oxygen`.
- Oxygen is required only at atoms the diff says need it. Speed. Filters also see the live mol as their first argument (`FilterRules` / `FilterSites` take `ForestTracingMol`); a non-cleaving oxygen-only pattern is refused when the live oxygen count already meets the target. Mol is not stuffed into `SiteInfo`.
- Every full-size target placement stays open to filters and cost. `_mappings` / `_merge_views` keep each placement (not one union core; not only the smaller remainder). Tests: `test_multi_mcs_union_would_prune_both_benzyls`, `test_smaller_secondary_embedding_alone_misses_naphthaldehyde`.
- Ring-membership disagreement on a mapped atom adds the bridge into the ring to `cleavage_bonds`. MCS may keep an exocyclic atom on a ring image (PhCH2OH CH2 onto a quinone carbon); the cut is still that bridge. Reads the mapping, not a rule-name branch. Unblocks benzylic dealkylation on PhCH2OH→quinone.
- Target-hit early stop: a child that is already the target gets priority tier 0 on the heap, and further site edits on that node stop once `max_paths` hits are queued (TBA).
- Two-site unique-edit identity is a pair-orbit signature from `graph_isomorphism` / `records` (beside the rank key, not in it): same-kind → ``AtomPairOrbitSignature`` / ``BondPairOrbitSignature``; DH ResonancePair → ``AtomPairOrbitSignature`` on the two end atoms. One-site rules still collapse by atom rank. **Nauty six-family API** (`all_site_pair_orbits_nauty`) is the source of truth: `atom_atom` / `bond_bond` each have unordered and ordered families; `atom_bond` ordered ≡ unordered (types already distinguish — do not invent a split). Same-kind unordered = combinations; ordered = permutations. **Ordered vs unordered for pair unique-edit** reads resolved ``swap_group`` (When → PatternInfo → ``name``; absent/None means use ``name``): unordered iff both ends share the same non-empty group; else ordered with ends sorted by ``PatternInfo.name``. Map-rank embeddings are part of the signature (dealkylate methyl vs benzyl). Status: approved. **QF vs DH vs H:** distinct endpoint ``name``s are distinct groups by default — no cross-role unordered; same-name pairs (`add_carbonyl_o`×2, `phenol_end`×2, Hydrogenation `path_end`×2) unordered. DH `phenol_end` / `amine_end` / `methide_end` stay ordered across roles despite shared `edit=single_to_double`. Explicit ``swap_group`` only when grouping differs from ``name``. Dehydrogenation declares ``site_kind = "atom_pair"`` / ``sites_on = "atom_pairs"`` (one-bond emits both endpoints). Unique-edit is atom–atom only; ``swap_group`` orders asymmetric ends. Site shape is class data (``site_kind``), not ``Generic[SiteT]``. Narrowing SMARTS cannot replace ``swap_group`` for these couples (HEURISTICS: not approved). **Opt-in ``canonical_emitted_sites``** (default off; explicit kwargs on metabolize / find_path / bfs / dfs only — no env): ``filter_sites`` sees **discovery**; after accept, chemistry + emitted ``site`` remapped to lex orbit rep (singleton atom and pair kinds); ``discovered_site`` escape hatch when they differ; product forest last-layer restamp only; lex-rep cache read from **parent** after product ``clear_structure``. Filters must not assume ``site`` == discovery. Status: approved. Classical group orbits on pairs are not a novelty claim — see ``PAIR_ORBITS.md`` References. Legacy `site_pair_orbits_nauty` maps to unordered same-kind + atom_bond. Tests: `test_pair_signatures.py`, `test_dh_atom_pair_unique_edit.py`, `test_graph_isomorphism.py`, `test_canonical_emitted_sites.py`, `test_site_kind.py`.
- Two dedup layers (do not conflate): (1) **Site unique-edit** — topological ranks + pair orbit signature; keeps site topology so equivalent embeddings are one edit before `RunReactants`. (2) **Product csmi** — split into **check** vs **yield**. Check (always while site unique-edit is on): after unique-edit survivors, each emission's frozenset of fragment CSMIs + site ranks; a later emission under the same `(rule, PatternInfo.name | SMARTS)` with equal `S` and equal ranks → `SiteDeduplicationWarning` (unique-edit miss). Check does not drop. Yield (`unique_csmi`, default on): drop duplicate `(rule, pattern, emission frozenset)` emissions from the stream — same `S` as check; within one cleavage list, identical sibling CSMIs stay. `unique_csmi=False` yields every emission after site unique-edit; check may still warn. `RuleSet.metabolize` forces children `unique_csmi=False` (yield off, check on) so alternate rules bubble up; the RuleSet's own yield (caller setting) is the cross-child CSMI layer — nested sets do the same, so only each asked-to-uniquify level drops. Cross-rule same emission frozenset → INFO log (kept/dropped rule + shared CSMI), not a Warning; parent does not re-fire same-rule `SiteDeduplicationWarning`. Inventory: ``REDUNDANT_RULES.md``. Status: approved.
- Pattern SMARTS that implement the same edit must not nest applicability ranges (H-count nesting that double-emits under distinct names). Partition like Dealkylation `#6H0` vs `#6h`, Hydroxylation `#6h1` vs `#6h2,#6h3`, Glutathionation epoxide/aziridine `#6H0` for the substituted-carbon row. Residual multi-name same-csmi cases that are chemically distinct edits are listed under DIVERGENCES.md. Status: approved for partitions done; residuals documented.
- `find_path` and `ReactionRule.metabolize` / `RuleSet.metabolize` raise `ValueError` on `None` input (invalid parse must not soft-pass as any-path).
- Lazy heap / stale-priority handling on the find_path frontier. Status: approved. Minimal key: `(0 if target_hit else 1, seq)` — hits before non-hits, then FIFO. On pop, cheaply rescore; if the fresh key differs and is worse than peek, `heappush` back instead of expanding (compare-to-peek; no O(n) resort). Not a forest ranker. Test: `test_stale_priority_is_pushed_back`.
- `leave_count` on a cleaving effect is read by `filter_sites`: when it is an int, the smaller fragment across the cleaved bond must have that many heavy atoms. Test: `test_leave_count_one_refuses_a_larger_leaving_fragment`.
- A nitrogen with two double bonds is not a product. `C=[N+]=C` sanitizes (degree 2, valence 4) but it is not an iminium. An iminium is one double bond and two single bonds. `_sanitize_piece` drops that fragment, and the split fails with it. Carbamazepine site `{4, 17}` was iminium on one ring carbon plus dealkylation of the carboxamide on the other, which left `C1=c2ccccc2=[N+]=c2ccccc2=C1` and `NC=O`. Status: approved. Test: `test_carbamazepine_does_not_emit_two_double_nitrogen`.
- A dearomatizing pair edit must change the system that was kekulized. Other aromatic systems may stay. If every aromatic atom of that system is still aromatic and no non-aromatic double or triple bond touches it, the product re-aromatized and is dropped. Alprazolam site `{1, 4}` is that cation: `Cc1nnc2cnc(-c3ccccc3)c3cc(Cl)ccc3[n+]1-2`. A quinone, quinone-imine, or quinodimethane keeps a localized double bond, so it stays. Status: approved. Test: `test_alprazolam_re_aromatized_cation_is_not_a_quinone`.
- Methide is always on the rule data (Dehydrogenation `methide_end`, QuinoneFormation alkyl branch). Opt-in `pathways=("methide",)` is dropped (DROPPED.md). Both ends may be methide: a para-quinodimethane is two alkyl single-to-double ends. The one-side skip was a mistake. There is no separate "no methides" mode. `merge_effects` sets `methide` when either end has it. Search already skips alkyl `partner=="C"` ends unless `_alkyl_bond_raises`. Tests: `test_dh_methide_pathways.py`.
- **Adds-H vs removes-H filters (Hydrogenation ≠ Dehydrogenation).** `filter_sites` already refused removes-H effects unless some site atom has `h_delta < 0` (or intersects `loses_aromaticity`). Symmetric: effects whose `adds` contains H (and that do not also add O / cleave) are refused unless some site atom has `h_delta > 0`. `filter_rules` mirrors this on span.adds. **Hydrogenation adds H** (reduction / saturation); **Dehydrogenation removes H**. Do not conflate. `dearomatizes` alone is not enough toward a quinone target: both H and DH can dearomatize, but only DH/QF match oxidative `h_delta`. Status: approved. Tests: `test_hydrogenation_pattern_info.py`.
- **Hydrogenation `path_end` PatternInfo.** Declares `adds="H"` and `dearomatizes=True` (capability); `merge_effects` resolves dearomatizes against `system_aromatic`. Alkene SMARTS keep span `dearomatizes=False` so aliphatic C=C→CC is not refused by the pattern-level “all dearomatizes” check when `loses_aromaticity` is empty. The same check **skips** patterns whose span also adds H (path_end / alkene): those are gated by the adds-H filter instead, so aliphatic carbonyl path reduction (`CC=O`→`CCO`) is not refused. Status: approved.

## Schema proposals (not yet in PatternInfo)

- **when → name.** Optional: a resolved `When` refines the pattern label so a globbed SMARTS OR (e.g. `[#6h2,#6h3]`) can still report `h2` vs `h3` without a second overlapping pattern. Candidates: `When` grows optional `name=`; or `PatternInfo.name` is a base and the chosen branch supplies a suffix. `_unique_csmi_key` / trace `pattern` would need a defined resolved token (static `PatternInfo.name` vs effect-time name). Prefer data on `When` / `Effect`, not a code branch in metabolize. Status: not decided — do not invent the field until a caller needs branch-specific names after partition. Test hook ready: `emitable_pattern_names` / `optional_when_name` in `tests/forest/pattern_info_inventory.py` already fold an optional `When["name"]` into per-rule uniqueness (`test_emitable_names_unique_within_each_reaction_rule`).

## Dedup check vs yield (approved)

``unique_csmi`` controls **yield** only (suppress duplicate *emissions* whose
fragment-CSMI frozenset already appeared under the same rule/pattern). **Check**
(emission frozenset of product CSMIs + site ranks → ``SiteDeduplicationWarning``
on unique-edit miss) runs whenever site unique-edit is on — including with yield
off, so a RuleSet can force children ``unique_csmi=False`` for cross-rule
bubble-up without losing within-rule miss detection. Site unique-edit has no
public off switch today; if one lands, check-off with it. ``product.xf.csmi``
stays cached after first read. Emission identity for check/yield is
``frozenset(p.xf.csmi for p in products)`` computed internally — not a
``csmi`` field on ``info`` (callers use ``xf.csmi``). Status: approved.

## Schema in PatternInfo (approved)

- **``swap_group: str``** on `PatternInfo` (optional When override). Pair unique-edit resolves as When → PatternInfo → ``name``. Same non-empty resolved group on both ends → unordered; unequal → ordered by ``name``. Omit the field when it would equal ``name``; set it only when grouping should differ. Status: approved. Tests: `test_pair_signatures.py`.

## Site kind taxonomy (approved)

`RuleSiteKind`: ``"atom"`` | ``"bond"`` | ``"directed_bond"`` | ``"atom_pair"``.

- **``atom``** — singleton Site; unique-edit uses directed `MapRankKey` `((mapno, rank), …)`.
- **``bond``** — undirected bond frozenset (two adjacent atoms). Unique-edit first field is sorted site ranks (`bond_rank_key`). Epoxidation, AzoSplitting. ResonancePair **SMARTS** (Hydrogenation alkene/alkyne, Dehydrogenation one-bond) also use `bond_rank_key` via `site_signature` when `site_kind="atom_pair"` — pair-path emissions still use `pair_site_signature`.
- **``directed_bond``** — unique-edit uses ordered map ranks (`MapRankKey` / ``site_map`` order) because map 1 is chemically distinct (Dealkylation / NDealkylation / Benzodioxole / Nitroaromatic). **Public API:** ``info["site"]`` is always a frozenset; orientation is on ``info["discovered_site"]`` as the ordered map-order tuple (same ``Site`` union — no extra field). With ``canonical_emitted_sites``, ``site`` is the frozenset of the lex representative and ``discovered_site`` remains the directed discovery tuple. Coercion is yield-only — unique-edit `seen` never keys on the frozenset.
- **``atom_pair``** — ResonancePair ends only (DH / QF / Hydrogenation / TautomerRule stub). Never on plain SMARTS. Pair unique-edit stays `pair_site_signature` + `swap_group`. One-bond SMARTS on these rules share undirected bond ranks with ``bond`` (see above).

Status: approved. Tests: `test_epoxidation_unique_edit.py`, `test_site_kind.py`, `test_parity.py` (anisole dealk).

## Antipattern: global `_kekule_forms` on SmirksReactionRule

Status: **not approved**. Materializing every Kekulé form inside plain `SmirksReactionRule.metabolites` is combinatorial expansion — that is why `ResonanceRule` exists. Match once on the aromatic parent; react on a cached Kekulé parent selected by SMARTS-implied bond order. Do not reintroduce an all-forms loop or an easy opt-in that restores the tax. `_kekule_forms` may remain for tests / helpers that need the list explicitly.

## Narrowing SMARTS vs ``swap_group`` (ResonancePair)

Status: **not approved** as a replacement for ``swap_group`` / name-default groups. Hydroxylation-style H-count partitions fix *nested same-atom* SMARTS that double-emit under ``unique_csmi``. Pair same-role couples are different: two path ends that both match the same edit (hydroquinone ``phenol_end``×2, benzene ``add_carbonyl_o``×2, diene ``path_end``×2, QF ``dealkylate``×2, …) are real chemistry and stay unordered via resolved group (= ``name``). Narrowing SMARTS so those couples never co-apply would drop pathways. Keep resolved groups + nauty ordered/unordered. Residual: map ranks still distinguish dealkylate embeddings.

## TautomerRule (stub only)

Status: not decided for chemistry / unique-edit; stub lands first.

Archived ``Tautomerization`` extended alternating paths by one H-bearing
neighbor and flipped bonds (net heavy-atom formula and H count unchanged;
path swap, not ``RunReactants``). Live ``TautomerRule`` subclasses
``ResonancePairRule``, is patternless, raises ``NotImplementedError`` from
``metabolites``, and is not in PhaseOne. Orthogonal deferred work: tautomer
*SMARTS matching* via RDKit ``TautomerQuery`` (preferred direction only;
Status: not decided). When the rule is implemented, put ends/effects on
``PatternInfo`` rather than a silent search branch. Docstring on the class
and TODO.md carry the same split. Test: `test_tautomer_rule_stub.py`.
