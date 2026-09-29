"""RDKit Python reference / validation engine — **feature-frozen**.

Not the product API. Prefer :mod:`xenosite.forest` (Rust chematic door).
**New features** (APIs, rules, search, sampling, conjugation helpers, …) go in
``crates/xenosite-forest`` and the ``xenosite.forest`` wrapper — **not here**.
See ``docs/forest/NATIVE.md``.

Call signatures may diverge from the Rust door over time. Allowed native edits
are parity bugfixes, CI, and docs that restate this freeze.
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
    "xenosite.forest.native is a frozen RDKit reference/validation engine; "
    "new features go in xenosite.forest (Rust) / crates/xenosite-forest, "
    "not here. See docs/forest/NATIVE.md.",
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

from . import rules  # noqa: E402
from .atom_tracker import AtomTracker  # noqa: E402
from .find_path import (  # noqa: E402
    CleavageSide,
    Maybe,
    PathCounters,
    PathOutcome,
    bfs,
    dfs,
    find_path,
)
from .phaseone import PhaseOneQF, PhaseOneRS, metabolize, reaction_labels  # noqa: E402
from .rulesets import PhaseOne, RuleSet  # noqa: E402

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
