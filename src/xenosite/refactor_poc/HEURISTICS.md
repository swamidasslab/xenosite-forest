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
- A full priority heap instead of the queue. Order transforms by how likely they are to help, and let that priority change after a path has already been found. Partial work is under Already in (target-hit prepend; rule order from span+diff). Sibling cost-sort before enqueue looked like a heap but preferred locally cheap dead ends on PhCH2OH→quinone; do not revive without a better key.
- Pattern frequency. Count how often each pattern matches, and how often that match shows up in a real metabolite. Deprioritize or filter patterns that match almost everything and are rarely the one that was used. That number belongs on `PatternInfo`, not in a table inside `find_path`.
- Cleavage-first may not need its own branch once a heap prefers cleaving spans when the target is smaller. Today's `order_key` still reads `cleaves` on the span; do not add another cleavage predicate until a heap has been tried.

## Already in the search

- Cleavage patterns run first when the target is smaller or a mapped bond is missing. Speed, not a new rule class. `order_key` also prefers dearomatizing spans when `loses_aromaticity`, and oxygen-adding spans when `needs_oxygen`.
- Oxygen is required only at atoms the diff says need it. Speed. Filters also see the current mol: a non-cleaving oxygen-only pattern is refused when the live oxygen count already meets the target.
- Every full-size target placement stays open to filters and cost. `_mappings` / `_merge_views` keep each placement (not one union core; not only the smaller remainder). Tests: `test_multi_mcs_union_would_prune_both_benzyls`, `test_smaller_secondary_embedding_alone_misses_naphthaldehyde`.
- Ring-membership disagreement on a mapped atom adds the bridge into the ring to `cleavage_bonds`. MCS may keep an exocyclic atom on a ring image (PhCH2OH CH2 onto a quinone carbon); the cut is still that bridge. Reads the mapping, not a rule-name branch. Unblocks benzylic dealkylation on PhCH2OH→quinone.
- Target-hit early stop: a child that is already the target is prepended, and further site edits on that node stop once `max_paths` hits are queued (TBA).
- `leave_count` on a cleaving effect is read by `filter_sites`: when it is an int, the smaller fragment across the cleaved bond must have that many heavy atoms. Test: `test_leave_count_one_refuses_a_larger_leaving_fragment`.
- `find_path` and `ReactionRule.metabolize` / `RuleSet.metabolize` raise `ValueError` on `None` input (invalid parse must not soft-pass as any-path).
