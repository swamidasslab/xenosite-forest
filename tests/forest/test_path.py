"""Proof-of-concept find_path. Earlier cases stay; later increments append."""

import pytest
from rdkit import Chem

from xenosite._archive_forest.step_plan import Deps
from xenosite.forest.find_path import PathCounters, atom_diff, canon_smiles, find_path
from xenosite.forest.records import AtomRef
from xenosite.forest.rules import Hydroxylation, ReactionRule
from xenosite.forest.rulesets import RuleSet


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


def test_find_path_and_metabolize_reject_none():
    """Invalid parse → None must not soft-pass; typing is Mol in, not Mol | None."""

    with pytest.raises(ValueError, match="required"):
        list(find_path(None, "CCO"))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="required"):
        list(find_path("CC", None))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="required"):
        list(Hydroxylation().metabolize(None))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="required"):
        list(RuleSet((Hydroxylation,), name="H").metabolize(None))  # type: ignore[arg-type]


def test_ruleset_runs_children_and_filters_see_them():
    ruleset = RuleSet((Hydroxylation,), name="Forest")
    assert isinstance(ruleset, ReactionRule)
    assert [type(rule) for rule in ruleset] == [Hydroxylation]
    seen = []

    def filter_rules(mol, rule, info):
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
    assert isinstance(outcome.plan, Deps), message  # pyright: ignore[reportArgumentType]
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
    record = product._forest["atom_trace"]["records"][atom.GetProp("forestLabel")]
    added = record.get("added_by")
    assert added is not None
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
    import xenosite.forest.rules as forest_rules

    missing = [name for name in names if not hasattr(forest_rules, name)]
    if missing:
        pytest.skip("not ported yet: %s" % ", ".join(missing))
    return RuleSet(tuple(getattr(forest_rules, name) for name in names), name="Forest")


def _old_site_applies(reactant, target, *, ruleset=None, depth=None, ceiling=40):
    """Old guided search's ``site_applies``. Not the budget this search must beat."""

    from xenosite._archive_forest import PathSearchCounters
    from xenosite._archive_forest.guided_path import find_path as old_find_path

    # forest/ is outside pyright; these calls are runtime-only against that API.
    counters = PathSearchCounters()  # pyright: ignore[reportCallIssue]
    kwargs = {"max_paths": 1, "max_expansions": ceiling, "counters": counters}
    if ruleset is not None:
        kwargs["ruleset"] = ruleset
    if depth is not None:
        kwargs["depth"] = depth
    list(old_find_path(reactant, target, **kwargs))  # pyright: ignore[reportCallIssue, reportArgumentType]
    return counters.site_applies  # pyright: ignore[reportAttributeAccessIssue]


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


def test_phch2oh_to_quinone_via_benzylic_dealkylation():
    """MCS may map CH2 onto a quinone carbon; the ring bridge is still a cut."""

    from xenosite.forest.rulesets import PhaseOne

    reactant = "OCc1ccccc1"
    target = "O=C1C=CC(=O)C=C1"
    diff = atom_diff(reactant, target)
    assert diff.site_is_cleavage({1, 2}), diff.cleavage_bonds

    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            target,
            ruleset=PhaseOne,
            counters=counters,
            max_nodes=80,
        )
    )
    assert hits, _billed(counters)
    assert counters.mol_edits <= 40, _billed(counters)
    names = [step.rule for step in hits[0].plan.children]
    assert "Dealkylation" in names


def test_hard_multi_oxidation_quinones_that_forest_struggled_on():
    """Large / multi-oxidation quinones: forest finds a path under a node ceiling.

    Forest PhaseOneQF exhausts budget on 4-methoxyphenol → hydroxyquinone
    (see tests/test_find_path_phase1_plan_fuzz.py xfail). Orthocarbonate and
    naphthalene→1,4-NQ are the other multi-oxidation anchors.
    """

    from xenosite.forest.rulesets import PhaseOne

    cases = [
        ("COc1ccc(O)cc1", "O=C1C=C(O)C(=O)C(O)=C1", 200),
        ("COc1ccc(O)cc1", "O=C1C=CC(OC(O)O)=CC1=O", 80),
        ("c1ccc2ccccc2c1", "O=C1C=CC(=O)c2ccccc12", 40),
    ]
    for reactant, target, ceiling in cases:
        counters = PathCounters()
        hits = list(
            find_path(
                reactant,
                target,
                ruleset=PhaseOne,
                counters=counters,
                max_nodes=ceiling,
            )
        )
        message = "%s → %s %s" % (reactant, target, _billed(counters))
        assert hits, message
        assert counters.nodes <= ceiling, message


def test_leave_count_one_refuses_a_larger_leaving_fragment():
    """``leave_count`` on the effect is the named leaving size; filters read it."""

    from xenosite.forest.find_path import _leaving_heavy_counts, _site_could_help
    from xenosite.forest.rules import NDealkylation

    mol = Chem.MolFromSmiles("CCN(C)C")
    assert mol is not None
    # Ethyl vs methyl on the same nitrogen: leave_count 1 keeps only methyl.
    ethyl = _leaving_heavy_counts(mol, {1, 2})
    methyl = _leaving_heavy_counts(mol, {2, 3})
    assert ethyl is not None and min(ethyl) == 2
    assert methyl is not None and min(methyl) == 1

    diff = atom_diff(mol, "CNC")
    info = {
        "site": (1, 2),
        "rule": NDealkylation(),
        "options": {"cleaves": True, "leave_count": 1},
        "rxn_num": 0,
        "pattern": NDealkylation().smarts[0][1],
    }
    assert not _site_could_help((1, 2), info, diff, mol)
    info_ok = dict(info, site=(2, 3), options={"cleaves": True, "leave_count": 1})
    assert diff.site_is_cleavage({2, 3})
    assert _site_could_help((2, 3), info_ok, diff, mol)


def test_stale_priority_is_pushed_back():
    """Lazy heap: a rescored key worse than peek is deferred, not expanded."""

    from xenosite.forest.find_path import _stale_vs_peek, _walk_priority

    best = _walk_priority(target_hit=True, seq=1)
    heap = [(best, 1, object())]
    worse = _walk_priority(target_hit=False, seq=0)
    assert _stale_vs_peek(worse, heap)
    assert not _stale_vs_peek(_walk_priority(target_hit=True, seq=0), heap)
    assert not _stale_vs_peek(worse, [])


def test_filter_sites_receives_live_mol_without_closure():
    """``FilterSites`` takes the traced mol first; callers need not close over it."""

    seen = []

    def filter_sites(mol, site, info):
        seen.append(mol)
        return True

    list(Hydroxylation().metabolize(Chem.MolFromSmiles("CC"), filter_sites=filter_sites))
    assert seen
    assert all(getattr(item, "_forest", None) is not None for item in seen)
