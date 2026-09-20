This file lists what live ``xenosite.forest`` left out relative to the pre-swap archive. A listing is not approval.
Each entry has `Status: approved`, `not approved`, or `not decided`.
`approved` means the cut is accepted for the long term.
`not approved` means the capability comes back. The old shape may still change.
`not decided` means it is absent for now and nobody has accepted that.
Do not delete entries. Change the status in the same commit that decides.

## PathContext and its MCS bookkeeping

Status: not decided.

Lived in `src/xenosite/forest/path_context.py`. Multi-embedding rounds,
chemical disagreement, frontier versus deep interior, safe-to-drop
thresholds, and oxygen slack. The atom diff in `find_path.py` is the part
search needs: one pairing, then the local change at each atom.

Old tests: `tests/test_cleavage_site_prune.py`, `tests/test_guided_path_fuzz.py`
(the context object, not the metabolites).

## Formula hint algebra

Status: not decided.

`FormulaHint`, `FormulaAny`, `FormulaMatch`, and `expand_formula_effects`
in `path_context.py`. Local atom differences cover compatibility. The spec
language does not.

Old tests: `tests/test_formula_smarts_hints.py`.

## and_cleave_plan SMARTS replay

Status: not decided.

`and_cleave_plan` and the apply-replay correction in `guided_path.py`.
Precedence comes from `added_by` and the resolved `Effect`.

Old tests: `tests/test_and_cleave_plan_fuzz.py`,
`tests/test_guided_path_cleavage.py`.

## And / Or plan trees and plan JSON

Status: approved.

`And` / `Or` in `step_plan.py`, plan JSON canonicalization, and exact
linearization matching against the old emission. Dependency step plans
(`Deps`) replace them. `Step` and `Deps` are imported, not reimplemented.
Callers do not get the old JSON shape.

Old tests: `tests/test_step_plan.py`, `tests/test_phase1_steps_fuzz.py`
(the JSON and `And`/`Or` assertions, not the chemistry).

## AtomTracker as a base class, and phase1_steps on each rule

Status: not decided for `AtomTracker` as a base class.
Status: not approved for dropping the composite look-ahead. It does not
have to be a `phase1_steps` method on every rule, and the name need not
stay "phase I". See the phase-one entry.

`AtomTracker` in `base.py`, and `phase1_steps` on each rule in
`forest/rules.py`. Live forest rules do not subclass `AtomTracker`.

Old tests: `tests/test_phase1_steps_fuzz.py`, `tests/test_atom_trace_char.py`.

## Phase-one reconstruction module

Status: not approved.

The archived module (`phaseone.py`, `PhaseOneRS`, `PhaseOneQF`) is not what live
forest calls. The capability stays. Quinone rules are composite look-aheads.
A ruleset is still a rule. Epoxidation sits inside StableOxygenation among
peers; NDealkylation is its own ruleset (not UnstableOxygenation —
that group is Dealkylation + OxidativeDehalogenation). Archive
``phase1_steps`` for both leaves are degenerate singletons naming the leaf
(``Epoxidation`` / ``NDealkylation``), not the group. Live
``canonical_plan`` matches that identity. Status: not approved for a
group-name look-ahead on those leaves (misremembered; fixed).

Old tests: `tests/test_phaseone.py`.

## Search controls that only bound the old expansion

Status: not decided.

`max_expansions` as a public mode and the old depth counters. A billed
counter stays, as a measurement (`mol_edits + nodes`). Live `bfs` and `dfs`
enumerate one ruleset up to a depth cap. They take `filter_rules` and
`filter_sites`, and they do not use the atom diff or the closer drop.
Adding them is not approval to drop enumeration, and it is not a port of
the old budget machinery.

Old tests: `tests/test_classic_path_counters.py`, `tests/test_bfs_fuzz.py`,
`tests/test_path_outcome_hard.py` (the budget flag, not the metabolite).

## Conjugate star-label modes

Status: not approved as a general drop. They belong only to the conjugation rules.

Bare `*` versus a full glucuronide or glutathione structure. Do not build
this as a mode on every rule.

Old tests: `tests/test_conjugates.py`, `tests/test_star_expand.py`.

## Dehydrogenation methide opt-in

Status: approved to drop the opt-in. Methide itself is not dropped.

The rule's data always offers methide. `PatternInfo` resolves at most one
methide site. A caller who wants none uses `filter_sites`. There is no
`pathways=("methide",)` flag.

Old tests: `tests/test_dh_methide_pathways.py`.

## Glutathionation thiol opt-out

Status: not decided.

`include_thiol=False` removed the substrate-thiol SMARTS. The pattern stays
on the rule. A caller who wants no thiol sites uses `filter_sites` and reads
`symbol`.

Old tests: `tests/test_conjugates.py`.

## Pair info `path` and merge_effects `partners` / `aromatic`

Status: not decided.

`ResonancePairRule` used to stash the alternating-path atom list on
`ProductsOfReaction.info["path"]` after `filter_sites`. `merge_effects`
also wrote `partners` and `aromatic` onto the merged options. Nothing in
live forest read those keys. Closed `SiteInfo` / `Effect` TypedDicts reject
them; the writes are gone. Per-end partners stay on `info["ends"]`.


## Network module / find_network_paths (xenosite.forest.net)

Status: approved to drop from the live forest plan.

Archived `net.py` builds a NetworkX reaction graph from one-step
`find_path` metabolites (Rxn records, SMILES-order atom maps). That
network API / `find_network_paths` scaffolding is not in live forest.
MS1 / MS2 path finding (if pursued) stays on `find_path`-style
APIs, not a revived `xenosite.net` package.

Old code: `src/xenosite/_archive_forest/net.py`.
