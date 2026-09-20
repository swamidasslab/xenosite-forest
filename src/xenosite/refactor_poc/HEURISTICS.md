# Find-path heuristics and correctness

Consult this file before changing `find_path` behavior, and again when a search miss or a slow pair shows up. Add an entry when a new idea appears. Do not delete entries. Mark one done when the test that names it passes without `xfail`.

Follow `.cursor/rules/data-not-branches.mdc`. Prefer a generic decision that reads `PatternInfo` over a special case in the search.
`leave_count` and `breaks_ring` on `Effect` are that data: the same `cleaves` bit, with the leaving piece named or open, and ring opening on the bond.

## Exceptions

An exception to that rule is allowed only with a good reason. Record it here. Each entry has `Status: approved`, `not approved`, or `not decided`. A branch in the code without an `approved` entry is not an exception. None are approved.

A rule may later carry data that it is less likely, and a filter may read that. The schema is not decided. Do not invent the field, and do not treat the missing field as an exception.

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
- Two-site unique-edit identity is ``AtomPairOrbitSignature(groups=(ga, gb), pair_group=PairGroupId)`` from `graph_isomorphism` / `records` (beside the rank key, not in it). One-site rules still collapse by atom rank. Unified recipes `site_pair_orbits_smiles` / `site_pair_orbits_nauty` compute `atom_atom`, `bond_bond`, and `bond_atom` together; each orbit membership is a CIP-sorted tuple of index pairs. Forest: `structure["site_pair_orbits_nauty"]` up front when `pynauty` imports; `site_pair_orbits_smiles` on demand only when both ends have topeqiv size `> 1` (else `TRIVIAL_PAIR_GROUP = PairGroupId(-1)`). ``pair_group`` is sequential ``PairGroupId`` ``0..n-1`` from sorting CIP membership keys (atom site ``("atom", cip)``; bond site ``("bond", cip_lo, cip_hi)``; atom before bond; index membership tie-breaks equal CIP). CIP from ``cip_ids`` / ``CanonicalRankAtoms(..., breakTies=False)``. QuinoneFormation / two-atom SMARTS use `atom_pair_orbit_key`. Bond–atom orbits for Dehydrogenation: **Status: decided (approved)** — recipe and `bond_atom_orbit_key` exist; Dehydrogenation unique-edit not wired yet (see TODO). Recipe 2 unused. Status: approved for atom pairs and for bond–atom / Dehydrogenation. Tests: `test_graph_isomorphism.py` (invariants + pynauty oracle on partitions and canonical ids).
- Two dedup layers (do not conflate): (1) **Site unique-edit** — topological ranks + pair orbit signature; keeps site topology so equivalent embeddings are one edit before `RunReactants`. (2) **Product csmi dedup** (`unique_csmi`) — drops site topology; `ReactionRule.metabolize` and `RuleSet.metabolize` both key ``(rule name, PatternInfo.name | SMARTS, product csmi)`` via `_unique_csmi_key` (pair emissions without a pattern use ``None`` for the middle field). Different child rules that share a product structure both emit. Status: approved.
- `find_path` and `ReactionRule.metabolize` / `RuleSet.metabolize` raise `ValueError` on `None` input (invalid parse must not soft-pass as any-path).
- Lazy heap / stale-priority handling on the find_path frontier. Status: approved. The queue had no walk priority key (plain deque + appendleft). Minimal key matches old behavior: `(0 if target_hit else 1, seq)` — hits before non-hits, then FIFO. On pop, cheaply rescore; if the fresh key differs and is worse than peek, `heappush` back instead of expanding (compare-to-peek; no O(n) resort). Not a forest ranker. Test: `test_stale_priority_is_pushed_back`.
- `leave_count` on a cleaving effect is read by `filter_sites`: when it is an int, the smaller fragment across the cleaved bond must have that many heavy atoms. Test: `test_leave_count_one_refuses_a_larger_leaving_fragment`.
