"""Enumerate metabolite structures with Metabolic Forest reaction rules.

This tree is the promoted refactor (formerly ``xenosite.refactor_poc``).
The pre-swap implementation lives in ``xenosite._archive_forest`` for a
while (parity / H2H / StepPlan apply). Prefer ``mol.xf`` over ``AtomTracker``.
"""

try:
    from ._version import __version__
except ImportError:  # pragma: no cover - missing only in incomplete checkouts
    __version__ = "0.0.0"

# StepPlan apply layer still lives on the archived tree until re-homed.
from xenosite._archive_forest.step_plan import (  # noqa: E402
    And,
    AtomRef,
    Deps,
    Linearization,
    Or,
    Step,
    StepPlan,
    pathway_from_json,
    pathway_to_json,
)

from . import rules
from .atom_tracker import AtomTracker
from .find_path import (
    CleavageSide,
    Maybe,
    PathCounters,
    PathOutcome,
    bfs,
    dfs,
    find_path,
)
from .phaseone import PhaseOneQF, PhaseOneRS, metabolize, reaction_labels
from .rulesets import PhaseOne, RuleSet

# Compatibility alias for callers that used the old counter name.
PathSearchCounters = PathCounters

__all__ = [
    "__version__",
    "AtomTracker",
    "bfs",
    "dfs",
    "find_path",
    "PathOutcome",
    "PathCounters",
    "PathSearchCounters",
    "Maybe",
    "CleavageSide",
    "PhaseOne",
    "PhaseOneRS",
    "PhaseOneQF",
    "RuleSet",
    "rules",
    "metabolize",
    "reaction_labels",
    "And",
    "Or",
    "Deps",
    "AtomRef",
    "Linearization",
    "Step",
    "StepPlan",
    "pathway_from_json",
    "pathway_to_json",
]
