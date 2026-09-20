# Find-path heuristics and correctness

Consult this file before changing `find_path` behavior, and again when a search miss or a slow pair shows up. Add an entry when a new idea appears. Do not delete entries. Mark one done when the test that names it passes without `xfail`.

Follow `.cursor/rules/data-not-branches.mdc`. Prefer a generic decision that reads `PatternInfo` over a special case in the search.
`leave_count` and `breaks_ring` on `Effect` are that data: the same `cleaves` bit, with the leaving piece named or open, and ring opening on the bond.

## Exceptions

An exception to that rule is allowed only with a good reason. Record it here. Each entry has `Status: approved`, `not approved`, or `not decided`. A branch in the code without an `approved` entry is not an exception. None are approved.

A rule may later carry data that it is less likely, and a filter may read that. The schema is not decided. Do not invent the field, and do not treat the missing field as an exception.

## Migration (not executed)

When the POC replaces the live library: archive old `xenosite.forest` for a while, then move POC modules into the `xenosite.forest` tree. Aim for a seamless import/call-surface swap (`AtomTracker` facade names, `xf` where forest had tags). Status: done (archived to xenosite._archive_forest; POC promoted to xenosite.forest).

## Correctness, not yet done

- Do not drop a child only because `atom_diff` cost is not strictly lower. The `closer` check in `find_path` refuses a sideways step. Some real routes do not move the score down on every hop. Tried equal-cost sideways when the site touched cleavage / dearomatization / oxygen; PhCH2OH and multi-oxidation quinones blew up (`mol_edits` 21→68). Status: not decided — needs a tighter predicate than “touches the diff”.
- Do not walk a step plan that is the same steps in another order. Once a `Deps` has been yielded, a later walk that is only a reordering of those steps is not a new path.

## Speed, not yet done

- Put the isotope-bearing atom first in the reaction SMARTS. `SubstructMatch` starts at query atom 0 and does not reorder, so a common atom first is walked before the isotope is tested.
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
- Two-site unique-edit identity is a pair-orbit signature from `graph_isomorphism` / `records` (beside the rank key, not in it): same-kind → ``AtomPairOrbitSignature`` / ``BondPairOrbitSignature``; DH ResonancePair → ``BondAtomPairOrbitSignature`` with a **joint second-order** ``pair_group`` (orbit of ``((b1,a1),(b2,a2))`` — first-order end signatures alone are not enough). One-site rules still collapse by atom rank. **Nauty six-family API** (`all_site_pair_orbits_nauty`) is the source of truth: `atom_atom` / `bond_bond` each have unordered and ordered families; `atom_bond` ordered ≡ unordered (types already distinguish — do not invent a split). Same-kind unordered = combinations; ordered = permutations. **Ordered vs unordered for pair unique-edit** reads resolved ``swap_group`` (When → PatternInfo → ``name``; absent/None means use ``name``): unordered iff both ends share the same non-empty group; else ordered with ends sorted by ``PatternInfo.name``. Map-rank embeddings are part of the signature (dealkylate methyl vs benzyl). Status: approved. **QF vs DH vs H:** distinct endpoint ``name``s are distinct groups by default — no cross-role unordered; same-name pairs (`add_carbonyl_o`×2, `phenol_end`×2, Hydrogenation `path_end`×2) unordered. DH `phenol_end` / `amine_end` / `methide_end` stay ordered across roles despite shared `edit=single_to_double`. Explicit ``swap_group`` only when grouping differs from ``name``. Dehydrogenation uses `unique_orbit = "bond_atom"`. **Surprising but intentional:** a directed bond–atom signature is a unique-edit key (not a Site) and captures *both* topological site identity *and* edit direction (bond vs atom role) in one isomorphism class — direction is not a separate Site field. Atom–atom / bond–bond have no type-level direction; they use ``swap_group`` ordered/unordered only. Narrowing SMARTS cannot replace ``swap_group`` for these couples (HEURISTICS: not approved). **Opt-in ``canonical_emitted_sites``** (default off; explicit kwargs on metabolize / find_path / bfs / dfs only — no env): ``filter_sites`` sees **discovery**; after accept, chemistry + emitted ``site`` remapped to lex orbit rep (singleton atom and pair/bond–atom kinds); ``discovered_site`` escape hatch when they differ; product forest last-layer restamp only; lex-rep cache read from **parent** after product ``clear_structure``. Filters must not assume ``site`` == discovery. Status: approved. Classical group orbits on pairs are not a novelty claim — see ``PAIR_ORBITS.md`` References. Legacy `site_pair_orbits_nauty` maps to unordered same-kind + atom_bond. Tests: `test_pair_signatures.py`, `test_nauty_six_family_orbits.py`, `test_dh_bond_atom_unique_edit.py`, `test_graph_isomorphism.py`, `test_canonical_emitted_sites.py`.
- Two dedup layers (do not conflate): (1) **Site unique-edit** — topological ranks + pair orbit signature; keeps site topology so equivalent embeddings are one edit before `RunReactants`. (2) **Product csmi dedup** (`unique_csmi`) — drops site topology; `ReactionRule.metabolize` and `RuleSet.metabolize` both key ``(rule name, PatternInfo.name | SMARTS, product csmi)`` via `_unique_csmi_key` (pair emissions without a pattern use ``None`` for the middle field). Different child rules that share a product structure both emit. Status: approved.
- Pattern SMARTS that implement the same edit must not nest applicability ranges (H-count nesting that double-emits under distinct names). Partition like Dealkylation `#6H0` vs `#6h`, Hydroxylation `#6h1` vs `#6h2,#6h3`, Glutathionation epoxide/aziridine `#6H0` for the substituted-carbon row. Residual multi-name same-csmi cases that are chemically distinct edits are listed under DIVERGENCES.md. Status: approved for partitions done; residuals documented.
- `find_path` and `ReactionRule.metabolize` / `RuleSet.metabolize` raise `ValueError` on `None` input (invalid parse must not soft-pass as any-path).
- Lazy heap / stale-priority handling on the find_path frontier. Status: approved. The queue had no walk priority key (plain deque + appendleft). Minimal key matches old behavior: `(0 if target_hit else 1, seq)` — hits before non-hits, then FIFO. On pop, cheaply rescore; if the fresh key differs and is worse than peek, `heappush` back instead of expanding (compare-to-peek; no O(n) resort). Not a forest ranker. Test: `test_stale_priority_is_pushed_back`.
- `leave_count` on a cleaving effect is read by `filter_sites`: when it is an int, the smaller fragment across the cleaved bond must have that many heavy atoms. Test: `test_leave_count_one_refuses_a_larger_leaving_fragment`.
- Methide is always on the rule data (Dehydrogenation `methide_end`, QuinoneFormation alkyl branch). Opt-in `pathways=("methide",)` is dropped (DROPPED.md). Pair resolve skips both-methide ends; `filter_sites` that reads `options["methide"]` refuses methide when the caller wants none. Search already skips alkyl `partner=="C"` ends unless `_alkyl_bond_raises`. Tests: `test_dh_methide_pathways.py`.

## Schema proposals (not yet in PatternInfo)

- **when → name.** Optional: a resolved `When` refines the pattern label so a globbed SMARTS OR (e.g. `[#6h2,#6h3]`) can still report `h2` vs `h3` without a second overlapping pattern. Candidates: `When` grows optional `name=`; or `PatternInfo.name` is a base and the chosen branch supplies a suffix. `_unique_csmi_key` / trace `pattern` would need a defined resolved token (static `PatternInfo.name` vs effect-time name). Prefer data on `When` / `Effect`, not a code branch in metabolize. Status: not decided — do not invent the field until a caller needs branch-specific names after partition. Test hook ready: `emitable_pattern_names` / `optional_when_name` in `tests/forest/pattern_info_inventory.py` already fold an optional `When["name"]` into per-rule uniqueness (`test_emitable_names_unique_within_each_reaction_rule`).

## Schema in PatternInfo (approved)

- **``swap_group: str``** on `PatternInfo` (optional When override). Pair unique-edit resolves as When → PatternInfo → ``name``. Same non-empty resolved group on both ends → unordered; unequal → ordered by ``name``. Omit the field when it would equal ``name``; set it only when grouping should differ. Status: approved. Tests: `test_pair_signatures.py`.

## Narrowing SMARTS vs ``swap_group`` (ResonancePair)

Status: **not approved** as a replacement for ``swap_group`` / name-default groups. Hydroxylation-style H-count partitions fix *nested same-atom* SMARTS that double-emit under ``unique_csmi``. Pair same-role couples are different: two path ends that both match the same edit (hydroquinone ``phenol_end``×2, benzene ``add_carbonyl_o``×2, diene ``path_end``×2, QF ``dealkylate``×2, …) are real chemistry and stay unordered via resolved group (= ``name``). Narrowing SMARTS so those couples never co-apply would drop pathways. Keep resolved groups + nauty ordered/unordered. Residual: map ranks still distinguish dealkylate embeddings.

## TautomerRule (stub only)

Status: not decided for chemistry / unique-edit; stub lands first.

Forest ``Tautomerization`` extended alternating paths by one H-bearing
neighbor and flipped bonds (net heavy-atom formula and H count unchanged;
path swap, not ``RunReactants``). POC ``TautomerRule`` subclasses
``ResonancePairRule``, is patternless, raises ``NotImplementedError`` from
``metabolites``, and is not in PhaseOne. Orthogonal deferred work: tautomer
*SMARTS matching* via RDKit ``TautomerQuery`` (preferred direction only;
Status: not decided). When the rule is implemented, put ends/effects on
``PatternInfo`` rather than a silent search branch. Docstring on the class
and TODO.md carry the same split. Test: `test_tautomer_rule_stub.py`.
