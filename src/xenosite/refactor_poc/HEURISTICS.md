# Find-path heuristics and correctness

Consult this file before changing `find_path` behavior, and again when a search miss or a slow pair shows up. Add an entry when a new idea appears. Do not delete entries. Mark one done when the test that names it passes without `xfail`.

Follow `.cursor/rules/data-not-branches.mdc`. Prefer a generic decision that reads `PatternInfo` over a special case in the search.
`leave_count` and `breaks_ring` on `Effect` are that data: the same `cleaves` bit, with the leaving piece named or open, and ring opening on the bond.

## Exceptions

An exception to that rule is allowed only with a good reason. Record it here. Each entry has `Status: approved`, `not approved`, or `not decided`. A branch in the code without an `approved` entry is not an exception. None are approved.

A rule may later carry data that it is less likely, and a filter may read that. The schema is not decided. Do not invent the field, and do not treat the missing field as an exception.

## Correctness, not yet done

- Expand every full-size match of the target onto the start. `_best_mapping` in `find_path.py` keeps one embedding. Two placements (bis-benzyl to benzaldehyde, naphthyl-benzyl to naphthaldehyde) must both be searched. The union of those matches is not one conserved core, and the smaller leftover match is not the only one kept. Tests that go green when this is fixed: `test_multi_mcs_union_would_prune_both_benzyls` and `test_smaller_secondary_embedding_alone_misses_naphthaldehyde` in `tests/refactor_poc/test_ported.py`.
- Do not drop a child only because `atom_diff` cost is not strictly lower. The `closer` check in `find_path` refuses a sideways step. Some real routes do not move the score down on every hop.
- Filters must see the molecule under consideration, not only the diff. Today `filter_rules` and `filter_sites` are built from `atom_diff` (`find_path.py`, the call beside `parent_cost`). The current mol has information the diff does not: what the last transform already changed, which atoms are new, what the formula is now. The smallest fix is to pass that mol into the filter call. Do not add a one-off `if rule is Dealkylation` instead.
- Do not walk a step plan that is the same steps in another order. Once a `Deps` has been yielded, a later walk that is only a reordering of those steps is not a new path.

## Speed, not yet done

- Deferred metabolites. `metabolites` should be able to yield a light object instead of a sanitized mol. The object carries the anticipated formula (the same counts as `atom_trace["formula"]` and `delta_formula`) and says whether it is sure materialization can succeed. `metabolize` either returns that object or materializes it, so today's callers still get mols. A later hop materializes only when a filter or the plan needs the real mol. The same object is how a `StepPlan` can check that a later site was `added_by` an earlier step without replaying the edit.
- A priority heap instead of the queue. Order transforms by how likely they are to help, and let that priority change after a path has already been found.
- Pattern frequency. Count how often each pattern matches, and how often that match shows up in a real metabolite. Deprioritize or filter patterns that match almost everything and are rarely the one that was used. That number belongs on `PatternInfo`, not in a table inside `find_path`.
- Cleavage-first may not need its own branch. If the pattern says it cleaves, and the filter can see the mol, a heap that prefers those patterns when the target is smaller replaces the `order_key` special case. Do not add another cleavage predicate until that has been tried.

## Already in the search

- Cleavage patterns run first when the target is smaller or a mapped bond is missing. Speed, not a new rule class.
- Oxygen is required only at atoms the diff says need it. Speed. The filter is still diff-only; the correctness entry above is the gap.
