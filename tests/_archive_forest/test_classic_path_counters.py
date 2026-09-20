"""Classic RuleSet.find_path uses the same PathSearchCounters schema as guided."""

from rdkit import Chem

from xenosite._archive_forest import PathSearchCounters, find_path, load_ruleset


def _mols(r, t):
    return Chem.MolFromSmiles(r), Chem.MolFromSmiles(t)


def test_classic_counters_share_guided_schema_and_budget():
    classic = PathSearchCounters()
    r, t = _mols("CN(C)Cc1ccccc1", "O=Cc1ccccc1")
    hits = list(
        load_ruleset("ND").find_path(
            r, t, depth=2, max_expansions=20, counters=classic
        )
    )
    assert hits
    assert classic.rule_expansions >= 1
    assert classic.site_applies == classic.rule_expansions  # 1 bill / metabolize
    assert classic.sites_considered >= 1
    assert classic.mol_edits >= 1
    assert classic.billed() == classic.site_applies
    assert not classic.budget_exhausted

    guided = PathSearchCounters()
    list(
        find_path(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            ruleset="ND",
            depth=2,
            max_paths=1,
            max_expansions=20,
            counters=guided,
        )
    )
    assert set(classic.as_dict()) == set(guided.as_dict())


def test_classic_dict_counters_mirror_full_fields():
    d = {}
    r, t = _mols("CN(C)Cc1ccccc1", "O=Cc1ccccc1")
    hits = list(
        load_ruleset("ND").find_path(r, t, depth=2, max_expansions=10, counters=d)
    )
    assert hits
    assert "billed" in d
    assert "site_applies" in d
    assert "mol_edits" in d
    assert "sanitize_dropped" in d
    assert d["rule_expansions"] == d["site_applies"]
    assert d["budget_exhausted"] is False


def test_sanitize_dropped_counted_under_path_scope():
    """Constructed-but-unsanitizable metabolites increment sanitize_dropped."""
    from xenosite._archive_forest.utils import clean, path_counter_scope

    c = PathSearchCounters()
    bad = Chem.MolFromSmiles("C#C(C)C", sanitize=False)
    assert bad is not None
    with path_counter_scope(c):
        assert clean(bad) == []
    assert c.sanitize_dropped >= 1


def test_classic_tight_budget_exhausts():
    c = PathSearchCounters()
    r, t = _mols("c1ccccc1", "O=C1C=CC(=O)C=C1")
    list(
        load_ruleset("PhaseOneRS").find_path(
            r, t, depth=3, max_expansions=1, counters=c
        )
    )
    assert c.budget_exhausted
    assert c.billed() == 1
