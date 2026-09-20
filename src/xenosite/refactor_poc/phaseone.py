"""Phase I as rules, and the names already on each addition.

Call :data:`PhaseOne` or one of the grouped sets. Each product records the
rules that ran under ``atom_trace["additions"]``. :func:`reaction_labels`
reads that chain. The label is each rule's initialization name. A ruleset
with no name stays on the chain and contributes no label.

Epoxidation and N-dealkylation phase-I look-aheads stay deferred.
Tautomerization is not in these sets.
"""

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


def metabolize(mol, **kwargs):
    """Run :data:`PhaseOne`. Yields ``(product, info)``."""

    yield from PhaseOne.metabolize(mol, **kwargs)


def reaction_labels(addition) -> tuple[str, ...]:
    """Initialization name of each named rule on the addition's chain.

    ``addition`` is the dict at ``atom_trace["additions"][id]`` or an
    :class:`~xenosite.refactor_poc.records.Addition`. Both store the chain
    on ``rules``. A ruleset with no name stays on the chain and is skipped
    here. The pattern that fired is on ``addition["pattern"]``, not here.
    """

    chain = getattr(addition, "rules", None)
    if chain is None:
        chain = addition["rules"]
    labels = []
    for rule in chain:
        name = getattr(rule, "name", None)
        if name:
            labels.append(name)
    return tuple(labels)
