"""Proof-of-concept find_path. Earlier cases stay; later increments append."""

import pytest
from rdkit import Chem

from xenosite.forest.step_plan import Deps
from xenosite.refactor_poc.find_path import PathCounters, atom_diff, canon_smiles, find_path
from xenosite.refactor_poc.records import AtomRef
from xenosite.refactor_poc.rules import Hydroxylation, ReactionRule, RuleSet


def _billed(counters):
    return (
        "billed=%s mol_edits=%s nodes=%s skipped=%s considered=%s"
        % (
            counters.billed,
            counters.mol_edits,
            counters.nodes,
            counters.sites_skipped,
            counters.sites_considered,
        )
    )


def test_ruleset_runs_children_and_filters_see_them():
    ruleset = RuleSet((Hydroxylation,), name="Poc")
    assert isinstance(ruleset, ReactionRule)
    assert [type(rule) for rule in ruleset] == [Hydroxylation]
    seen = []

    def filter_rules(rule, info):
        seen.append((type(rule).__name__, "span" in info))
        return False

    products = list(
        ruleset.metabolize(Chem.MolFromSmiles("CC"), filter_rules=filter_rules)
    )
    assert products == []
    assert seen
    assert all(name == "Hydroxylation" and has_span for name, has_span in seen)


def test_counters_start_at_zero_and_count_one_hydroxylation_edit():
    counters = PathCounters()
    assert counters.rule_expansions == 0
    assert counters.sites_considered == 0
    assert counters.sites_skipped == 0
    assert counters.mol_edits == 0
    assert counters.sanitize_dropped == 0
    assert counters.nodes == 0
    assert counters.billed == 0

    products = list(
        Hydroxylation().metabolize(Chem.MolFromSmiles("CC"), counters=counters)
    )
    assert products
    assert counters.mol_edits == 1


def test_atom_diff_ethane_needs_one_oxygen_and_cleaves_nothing():
    diff = atom_diff("CC", "CCO")
    assert len(diff.needs_oxygen) == 1
    assert not diff.cleaved
    mol = Chem.MolFromSmiles("CC")
    carbon = next(iter(diff.needs_oxygen))
    assert mol.GetAtomWithIdx(carbon).GetAtomicNum() == 6


def test_atom_diff_butylbenzene_ring_carbons_are_not_marked():
    reactant = "c1ccc(CCCC)cc1"
    diff = atom_diff(reactant, "OCCCCc1ccccc1")
    mol = Chem.MolFromSmiles(reactant)
    ring = {atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing()}
    assert diff.needs_oxygen
    assert diff.needs_oxygen.isdisjoint(ring)
    assert not diff.cleaved


def test_find_path_ethane_to_ethanol_is_one_hydroxylation():
    counters = PathCounters()
    hits = list(find_path("CC", "CCO", counters=counters))
    message = _billed(counters)
    assert hits, message
    outcome = hits[0]
    assert isinstance(outcome.plan, Deps), message
    steps = outcome.plan.children
    assert len(steps) == 1, message
    assert steps[0].rule == "Hydroxylation", message
    assert outcome.smiles == canon_smiles("CCO"), message
    assert not outcome.maybe, message


def test_butylbenzene_chain_alcohol_skips_ring_sites():
    counters = PathCounters()
    hits = list(
        find_path("c1ccc(CCCC)cc1", "OCCCCc1ccccc1", counters=counters)
    )
    message = _billed(counters)
    assert hits, message
    assert hits[0].smiles == canon_smiles("OCCCCc1ccccc1"), message
    assert hits[0].plan.children[0].rule == "Hydroxylation", message
    assert counters.sites_skipped > 0, message
    assert counters.mol_edits == 1, message
    assert counters.mol_edits < counters.sites_considered, message


def test_anisole_to_phenol_cleaves_and_keeps_the_methyl_side():
    counters = PathCounters()
    hits = list(find_path("COc1ccccc1", "Oc1ccccc1", counters=counters))
    message = _billed(counters)
    assert hits, message
    outcome = hits[0]
    steps = outcome.plan.children
    assert len(steps) == 1, message
    assert steps[0].rule == "Dealkylation", message
    assert outcome.smiles == canon_smiles("Oc1ccccc1"), message
    assert outcome.maybe, message
    side = outcome.maybe.entries[0].side
    assert outcome.allows(side=side), message
    # The methyl fragment is not another search node.
    assert counters.nodes == 2, message


def test_benzene_to_quinone_is_two_hydroxylations_then_dehydrogenation():
    quinone = "O=C1C=CC(=O)C=C1"
    on = PathCounters()
    hits = list(find_path("c1ccccc1", quinone, counters=on))
    message = _billed(on)
    assert hits, message
    steps = hits[0].plan.children
    names = [step.rule for step in steps]
    assert names.count("Hydroxylation") == 2, message
    assert names.count("Dehydrogenation") == 1, message
    hydroxyl = [i for i, name in enumerate(names) if name == "Hydroxylation"]
    dehydrogenation = names.index("Dehydrogenation")
    edges = set(hits[0].plan.precedes)
    assert (hydroxyl[0], hydroxyl[1]) not in edges, message
    assert (hydroxyl[1], hydroxyl[0]) not in edges, message
    assert (hydroxyl[0], dehydrogenation) in edges, message
    assert (hydroxyl[1], dehydrogenation) in edges, message
    assert hits[0].smiles == canon_smiles(quinone), message

    off = PathCounters()
    off_hits = list(find_path("c1ccccc1", quinone, counters=off, use_filters=False))
    assert off_hits, _billed(off)
    assert on.billed < off.billed, "on %s; off %s" % (_billed(on), _billed(off))


def test_phenol_to_quinone_plan_has_no_quinone_formation_step():
    quinone = "O=C1C=CC(=O)C=C1"
    counters = PathCounters()
    hits = list(find_path("Oc1ccccc1", quinone, counters=counters))
    message = _billed(counters)
    assert hits, message
    names = [step.rule for step in hits[0].plan.children]
    assert "QuinoneFormation" not in names, message
    assert names.count("Hydroxylation") == 1, message
    assert names.count("Dehydrogenation") == 1, message
    hydroxyl = names.index("Hydroxylation")
    dehydrogenation = names.index("Dehydrogenation")
    assert (hydroxyl, dehydrogenation) in set(hits[0].plan.precedes), message
    assert hits[0].smiles == canon_smiles(quinone), message


def test_quinone_oxygen_ref_resolves_to_the_atom_hydroxylation_adds():
    quinone = "O=C1C=CC(=O)C=C1"
    hits = list(find_path("Oc1ccccc1", quinone))
    assert hits
    steps = hits[0].plan.children
    names = [step.rule for step in steps]
    hydroxyl = names.index("Hydroxylation")
    dehydrogenation = names.index("Dehydrogenation")
    assert (hydroxyl, dehydrogenation) in set(hits[0].plan.precedes)

    origin = next(ref.origin for ref in steps[hydroxyl].site if ref.origin is not None)
    ref = AtomRef(origin, "O")
    mol = Chem.MolFromSmiles("Oc1ccccc1")
    product = next(
        candidate
        for candidate, info in Hydroxylation().metabolize(mol)
        if origin in info["site"]
    )
    idx = ref.resolve(product)
    atom = product.GetAtomWithIdx(idx)
    assert atom.GetAtomicNum() == 8
    added = product._forest["atom_trace"]["records"][atom.GetProp("forestLabel")]["added_by"]
    addition = product._forest["atom_trace"]["additions"][added]
    assert addition["name"] == "Hydroxylation"
    assert origin in addition["site"]


def test_butylbenzene_to_quinone_and_chain_alcohol_bills_those_sites():
    reactant = "c1ccc(CCCC)cc1"
    target = "OCCCCc1cc(=O)ccc1=O"
    counters = PathCounters()
    hits = list(find_path(reactant, target, counters=counters))
    message = _billed(counters)
    assert hits, message
    names = [step.rule for step in hits[0].plan.children]
    assert "QuinoneFormation" not in names, message
    assert "Hydroxylation" in names, message
    assert "Dehydrogenation" in names, message
    assert hits[0].smiles == canon_smiles(target), message
    # Chain oxygenation plus one ring dearomatization, not every ring carbon.
    assert counters.mol_edits <= 4, message
    assert counters.sites_skipped > counters.mol_edits, message


def _ruleset(*names):
    import xenosite.refactor_poc.rules as poc

    missing = [name for name in names if not hasattr(poc, name)]
    if missing:
        pytest.skip("not ported yet: %s" % ", ".join(missing))
    return RuleSet(tuple(getattr(poc, name) for name in names), name="Poc")


def _old_site_applies(reactant, target, *, ruleset=None, depth=None, ceiling=40):
    """Old guided search's ``site_applies``. Not the budget this search must beat."""

    from xenosite.forest import PathSearchCounters
    from xenosite.forest import find_path as old_find_path

    counters = PathSearchCounters()
    kwargs = {"max_paths": 1, "max_expansions": ceiling, "counters": counters}
    if ruleset is not None:
        kwargs["ruleset"] = ruleset
    if depth is not None:
        kwargs["depth"] = depth
    list(old_find_path(reactant, target, **kwargs))
    return counters.site_applies


def test_tba_mol_edits_within_old_expansion_budget():
    """Terbinafine to the enyne aldehyde. Ceiling is the old ``max_expansions``."""

    reactant = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
    target = "CC(C)(C)C#CC=CC=O"
    ceiling = 40
    old = _old_site_applies(reactant, target, depth=2, ceiling=ceiling)
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            target,
            ruleset=_ruleset("Dealkylation"),
            counters=counters,
            max_nodes=ceiling,
        )
    )
    message = "%s; old site_applies=%s" % (_billed(counters), old)
    assert hits, message
    assert hits[0].smiles == canon_smiles(target), message
    assert counters.mol_edits <= ceiling, message


def test_acetate_to_catechol_mol_edits_within_budget():
    reactant = "CC(=O)Oc1ccc(OC)cc1"
    target = "Oc1ccc(O)cc1"
    ceiling = 80
    old = _old_site_applies(
        reactant, target, ruleset="PhaseOneRS", depth=3, ceiling=ceiling
    )
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            target,
            ruleset=_ruleset("Hydrolysis", "Dealkylation"),
            counters=counters,
            max_nodes=ceiling,
        )
    )
    message = "%s; old site_applies=%s" % (_billed(counters), old)
    assert hits, message
    assert hits[0].smiles == canon_smiles(target), message
    assert counters.mol_edits <= ceiling, message


def test_impossible_targets_stay_under_the_old_ceiling():
    cases = [
        ("c1ccccc1O", "C", 40),
        ("c1cccc2ccccc12", "C", 40),
        ("c1ccccc1", "FC(F)(F)F", 40),
        ("c1ccccc1", "[Fe]", 50),
    ]
    for reactant, target, ceiling in cases:
        counters = PathCounters()
        hits = list(
            find_path(reactant, target, counters=counters, max_nodes=ceiling)
        )
        message = "%s → %s %s" % (reactant, target, _billed(counters))
        assert not hits, message
        assert counters.mol_edits <= ceiling, message


def test_ndealkylation_hard_pairs_when_the_class_exists():
    """Macrocycle and dialdehyde on the N-dealkylation ruleset."""

    ruleset = _ruleset("NDealkylation")
    cases = [
        ("C1CCCCCCNC2CCCC(CC2)NCCCC1", "NC1CCCC(=O)CC1", 120),
        ("CN(C)Cc1ccc(CN(C)CC)cc1", "O=Cc1ccc(C=O)cc1", 80),
        (
            "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
            "O=Cc1ccc(C=O)cc1",
            35,
        ),
    ]
    for reactant, target, ceiling in cases:
        counters = PathCounters()
        hits = list(
            find_path(
                reactant,
                target,
                ruleset=ruleset,
                counters=counters,
                max_nodes=ceiling,
            )
        )
        message = "%s → %s %s" % (reactant, target, _billed(counters))
        assert hits, message
        assert counters.mol_edits <= ceiling, message
