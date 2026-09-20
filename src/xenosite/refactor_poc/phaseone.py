"""Phase I as rules, and the names already on each addition.

Call :data:`PhaseOne` or one of the grouped sets. Each product records the
rules that ran under ``atom_trace["additions"]``. :func:`reaction_labels`
reads that chain. The label is each rule's initialization name. A ruleset
with no name stays on the chain and contributes no label.

Epoxidation and N-dealkylation keep identity canonical plans (forest
``phase1_equivalent`` singletons). StableOxygenation / UnstableOxygenation
here are group RuleSets that contain those leaves among peers.
Tautomerization is not in these sets.
"""

from __future__ import annotations

from collections.abc import Generator, Sequence
from typing import Any

from xenosite.refactor_poc.rdkitutil import Mol, TracingMol
from xenosite.refactor_poc.records import Addition, ProductInfo, TraceAddition
from xenosite.refactor_poc.rules import ReactionRule

from . import rules
from .rulesets import PhaseOne, RuleSet

StableOxygenation_PhaseOne = RuleSet(
    (
        rules.Hydroxylation,
        rules.Epoxidation,
        rules.SulfurOxidation,
        rules.NitrogenOxidation,
    ),
    name="SO",
    longname="StableOxygenation",
)
Dehydrogenation_PhaseOne = RuleSet(
    (rules.Dehydrogenation,),
    name="DH",
    longname="Dehydrogenation",
)
Hydrolysis_PhaseOne = RuleSet(
    (rules.Dephosphorylation, rules.EpoxideOpening, rules.Hydrolysis),
    name="HD",
    longname="Hydrolysis",
)
Reduction_PhaseOne = RuleSet(
    (
        rules.Dehydration,
        rules.Hydrogenation,
        rules.NitrogenReduction,
        rules.OxygenReduction,
        rules.ReductiveDehalogenation,
        rules.SulfurReduction,
    ),
    name="RD",
    longname="Reduction",
)
UnstableOxygenation_PhaseOne = RuleSet(
    (rules.Dealkylation, rules.OxidativeDehalogenation),
    name="UO",
    longname="UnstableOxygenation",
)

# The grouped phase-I reactions, without quinone. :data:`PhaseOne` is the
# ruleset already named in rulesets, and it includes QuinoneFormation.
PhaseOneRS = RuleSet(
    (
        Dehydrogenation_PhaseOne,
        Hydrolysis_PhaseOne,
        Reduction_PhaseOne,
        StableOxygenation_PhaseOne,
        UnstableOxygenation_PhaseOne,
    ),
)

# Quinone formation is one rule in front of phase I. It is not a second
# copy of the reaction already inside :data:`PhaseOne`.
PhaseOneQF = RuleSet(
    (rules.QuinoneFormation, PhaseOneRS),
    name="PhaseOneQF",
    longname="Phase I + Quinone Formation",
)


def metabolize(
    mol: Mol, **kwargs: Any
) -> Generator[tuple[TracingMol, ProductInfo], None, None]:
    """Run :data:`PhaseOne`. Yields ``(product, info)``."""

    yield from PhaseOne.metabolize(mol, **kwargs)


def reaction_labels(addition: Addition | TraceAddition) -> tuple[str, ...]:
    """Initialization name of each named rule on the addition's chain.

    ``addition`` is the dict at ``atom_trace["additions"][id]`` or an
    :class:`~xenosite.refactor_poc.records.Addition`. Both store the chain
    on ``rules``. A ruleset with no name stays on the chain and is skipped
    here. The pattern that fired is on ``addition["pattern"]``, not here.
    """

    chain: Sequence[ReactionRule]
    if isinstance(addition, Addition):
        chain = addition.rules
    else:
        chain = addition["rules"]
    labels: list[str] = []
    for rule in chain:
        name = getattr(rule, "name", None)
        if name:
            labels.append(name)
    return tuple(labels)
