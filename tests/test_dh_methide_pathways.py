"""Optional query_smarts pathway options (methide DH) + formula recommend."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest import PathSearchCounters, RuleSet, StepPlan, find_path
from xenosite.forest.guided_path import _canon
from xenosite.forest.rules import Dehydrogenation, Hydroxylation, QuinoneFormation
from xenosite.forest.utils import canon_smi, unmapped_smiles


def _products(rule, smi, **kw):
    mol = Chem.MolFromSmiles(smi)
    out = set()
    for _site, prods in rule.metabolize(
        mol,
        tag_atoms=False,
        only_emit_topologically_distinct_sites=True,
        **kw,
    ):
        for p in prods or []:
            if p:
                out.add(canon_smi(unmapped_smiles(p)))
    return out


def test_query_smarts_options_roundtrip_in_match_queries():
    """Third-field options appear on match_queries hits."""
    dh = Dehydrogenation()
    mol = Chem.MolFromSmiles("Oc1ccccc1C")
    hits = dh.match_queries(dh.standardize(mol))
    methide_opts = [
        opt
        for _atom, lst in hits.items()
        for _maps, _mods, opt in lst
        if opt.get("pathway") == "methide"
    ]
    assert methide_opts
    assert methide_opts[0].get("one_sided") is True
    assert "recommend" in methide_opts[0]


def test_dh_methide_off_by_default():
    assert "C=C1C=CC=CC1=O" not in _products(Dehydrogenation(), "Oc1ccccc1C")


def test_dh_methide_opt_in_instance_and_kwarg():
    want = "C=C1C=CC=CC1=O"
    assert want in _products(Dehydrogenation(pathways={"methide"}), "Oc1ccccc1C")
    assert want in _products(
        Dehydrogenation(), "Oc1ccccc1C", pathways=("methide",)
    )


def test_dh_methide_one_sided_never_both_ends():
    """Xylene has two alkyl ends; bis-methide pairs must not emit."""
    # No hetero end → only methide+methide pairs possible → empty when pathway on.
    assert _products(Dehydrogenation(pathways={"methide"}), "Cc1ccccc1C") == set()


def test_dh_pathways_toward_uses_formula_and_recommend():
    dh = Dehydrogenation()
    r = Chem.MolFromSmiles("Oc1ccccc1C")
    t = Chem.MolFromSmiles("C=C1C=CC=CC1=O")
    assert "methide" in dh.pathways_toward(r, t)
    # Formula mismatch: more O in target than reactant without prep — still
    # equal heavy formulas here; use a target that differs in O count.
    t_rich = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    assert "methide" not in dh.pathways_toward(
        Chem.MolFromSmiles("c1ccccc1"), t_rich
    )


def test_qf_phase1_stamps_methide_and_applies():
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    qf = QuinoneFormation()
    found = False
    for _site, prods in qf.metabolize(
        mol, tag_atoms=False, attach_phase1_steps=True
    ):
        for p in prods or []:
            if not p:
                continue
            ps = canon_smi(unmapped_smiles(p))
            if ps != "C=C1C=CC=CC1=O":
                continue
            plan = StepPlan.from_mol(p)
            dh = [s for s in plan.steps if s.rule == "Dehydrogenation"]
            assert dh and "methide" in dh[0].pathways
            assert any(
                ps in [canon_smi(unmapped_smiles(x)) for x in (lin.apply(mol) or [])]
                for lin in plan.linearizations()
            )
            found = True
    assert found


def test_guided_find_toluene_omethide_via_recommend():
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "Cc1ccccc1",
            "C=C1C=CC=CC1=O",
            ruleset=RuleSet([Hydroxylation(), Dehydrogenation()], name="t"),
            depth=3,
            maybe_prefixes=False,
            max_paths=1,
            max_expansions=40,
            counters=counters,
        )
    )
    assert hits
    assert _canon(hits[0][0][-1]) == _canon("C=C1C=CC=CC1=O")
    assert counters.rule_expansions <= 40
