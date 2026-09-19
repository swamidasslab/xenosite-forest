"""Proof-of-concept find_path. Earlier cases stay; later increments append."""

from rdkit import Chem

from xenosite.forest.step_plan import Deps
from xenosite.refactor_poc.find_path import PathCounters, atom_diff, canon_smiles, find_path
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
