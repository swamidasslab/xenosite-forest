"""Frozen 0.6.x Metabolic Forest archive (not the public API).

Import modules explicitly, e.g. ``xenosite.forest.legacy.rules``.
Public callers should use ``xenosite.forest`` (Rust). See README.md and
``docs/forest/LEGACY.md`` in this tree.
"""

from __future__ import annotations

import warnings

try:
    import rdkit  # noqa: F401
except ImportError as e:  # pragma: no cover - optional extra
    raise ImportError(
        "xenosite.forest.legacy requires RDKit; install with "
        "pip install 'xenosite-forest[rdkit]' "
        "(or uv add 'xenosite-forest[rdkit]')."
    ) from e

warnings.warn(
    "xenosite.forest.legacy is a frozen 0.6.x archive and is not the "
    "supported API; use xenosite.forest (Rust).",
    DeprecationWarning,
    stacklevel=2,
)

from .bfs import bfs, dfs  # noqa: E402
from .guided_path import (  # noqa: E402
    PathOutcome,
    PathSearchCounters,
    find_path,
    find_path_guided,
)
from .path_context import (  # noqa: E402
    ADD_O,
    CLEAVE,
    CLEAVE_OR_ADD_O,
    NEUTRAL,
    REMOVE_O,
    FormulaHint,
)
from .phaseone import PhaseOneQF, PhaseOneRS  # noqa: E402
from . import rules  # noqa: E402
from .rulesets import RULESETS, RuleSet, load_ruleset  # noqa: E402
from .step_plan import And, AtomRef, Deps, Linearization, Or, Step, StepPlan  # noqa: E402
from .step_plan import pathway_from_json, pathway_to_json  # noqa: E402
from .trace import AtomTrace  # noqa: E402

__all__ = [
    "bfs",
    "dfs",
    "find_path",
    "find_path_guided",
    "PathOutcome",
    "PathSearchCounters",
    "FormulaHint",
    "ADD_O",
    "CLEAVE",
    "CLEAVE_OR_ADD_O",
    "NEUTRAL",
    "REMOVE_O",
    "PhaseOneRS",
    "PhaseOneQF",
    "RULESETS",
    "RuleSet",
    "load_ruleset",
    "rules",
    "And",
    "Or",
    "Deps",
    "AtomRef",
    "Linearization",
    "Step",
    "StepPlan",
    "pathway_from_json",
    "pathway_to_json",
    "AtomTrace",
]
