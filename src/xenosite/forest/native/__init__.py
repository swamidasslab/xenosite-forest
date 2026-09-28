"""RDKit Python reference / validation engine.

Not the recommended product API. Prefer :mod:`xenosite.forest` (Rust).
Call signatures may diverge from the Rust door over time.
"""

from __future__ import annotations

import warnings

try:
    import rdkit  # noqa: F401
except ImportError as e:  # pragma: no cover - optional extra
    raise ImportError(
        "xenosite.forest.native requires RDKit; install with "
        "pip install 'xenosite-forest[rdkit]' "
        "(or uv add 'xenosite-forest[rdkit]')."
    ) from e

warnings.warn(
    "xenosite.forest.native is the RDKit reference/validation engine; "
    "the recommended API is xenosite.forest (Rust).",
    UserWarning,
    stacklevel=2,
)

from xenosite.forest.legacy.step_plan import (  # noqa: E402
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
