"""Hypothesis fuzz: ``and_cleave_plan`` on multi-step quinone phase1 walks.

Multi-step ``QuinoneFormation.phase1_steps`` branches (e.g. ``And(OH, OH) → DH``,
``And(Dealk, OH) → DH``) are expanded to leaf Steps, replayed on a metabolite
chain, and re-emitted as ``Deps``. Emission must allow every phase1
linearization (not be stricter); apply-replay may widen vs phase1 And/Seq.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from hypothesis.database import DirectoryBasedExampleDatabase
from rdkit import Chem

from xenosite._archive_forest.guided_path import and_cleave_plan, _canon
from xenosite._archive_forest.rules import QuinoneFormation
from xenosite._archive_forest.step_plan import (
    And,
    AtomRef,
    Deps,
    Linearization,
    Step,
    StepPlan,
)
from xenosite._archive_forest.utils import unmapped_smiles

_HYPOTHESIS_DIR = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
_HYPOTHESIS_DB = DirectoryBasedExampleDatabase(str(_HYPOTHESIS_DIR))


def _deps_from_plan(plan: StepPlan) -> Deps:
    """Test helper: flat Deps view of a Seq-of-And / Deps plan."""
    if isinstance(plan, Deps):
        return plan
    steps = tuple(plan.steps)
    return Deps(steps, plan.precedes) if steps else Deps(())

_CORPUS = (
    "c1ccccc1",
    "Oc1ccccc1",
    "Oc1ccc(O)cc1",
    "Nc1ccc(O)cc1",
    "COc1ccc(O)cc1",
    "Cc1ccc(O)cc1",
    "CC(=O)Nc1ccc(O)cc1",
    "Clc1ccc(O)cc1",
    "c1ccc2ccccc2c1",
)

_MAX_HEAVY = 22


def _mol(smi: str):
    return Chem.MolFromSmiles(smi)


def _target_of(mol, steps) -> str | None:
    products = Linearization(tuple(steps)).apply(mol)
    if not products:
        return None
    return _canon(unmapped_smiles(products[0]))


def _reaches(mol, steps, target: str) -> bool:
    products = Linearization(tuple(steps)).apply(mol)
    if not products:
        return False
    return target in {_canon(unmapped_smiles(m)) for m in products}


def _mol_chain(mol, steps) -> list | None:
    """Apply ``steps`` in order; return ``[R, …, T]`` or None on failure."""
    if not Linearization(tuple(steps)).apply(mol):
        return None
    chain = [Chem.Mol(mol)]
    cur = Chem.Mol(mol)
    for step in steps:
        products = step.apply(cur)
        if not products:
            return None
        cur = max(products, key=lambda m: m.GetNumHeavyAtoms())
        chain.append(cur)
    if len(chain) != len(steps) + 1:
        return None
    return chain


def _unique_steps(steps) -> bool:
    return len(steps) == len(set(steps))


def _lin_set(plan: StepPlan) -> set:
    return set(plan.iter_linearizations())


def _assert_lins_reach(mol, plan, target: str):
    lin_orders = _lin_set(plan)
    assert lin_orders, plan
    for order in lin_orders:
        assert _reaches(mol, order, target), (
            "linearization missed target",
            plan,
            [str(s) for s in order],
            target,
        )
    return lin_orders


def _assert_plan_matches_chemistry(mol, steps, target: str, mols=None, *, exact=True):
    """All plan lins reach ``target``.

    When ``exact`` is True, every other permutation must miss (strict). For
    quinone phase1 spines, prefer comparing lin sets to the phase1 branch —
    some non-lins may still reach T even when phase1 requires Dealk before DH.
    """
    if mols is None:
        mols = _mol_chain(mol, steps)
    assert mols is not None
    plan = and_cleave_plan(steps, mols=mols)
    lin_orders = _assert_lins_reach(mol, plan, target)
    if exact:
        for order in permutations(steps):
            if order in lin_orders:
                continue
            assert not _reaches(mol, order, target), (
                "non-linearization reached target",
                plan,
                [str(s) for s in order],
                target,
            )
    return plan


def _assert_emitted_respects_phase1_deps(mol, branch: StepPlan):
    """Emitted Deps matches apply-replay gold (closure identity)."""
    from xenosite._archive_forest.guided_path import _correct_plan_by_apply_replay

    lin = next(branch.linearizations())
    steps = list(lin.steps)
    mols = _mol_chain(mol, steps)
    assert mols is not None
    target = _canon(unmapped_smiles(mols[-1]))
    emitted = _assert_plan_matches_chemistry(
        mol, steps, target, mols=mols, exact=False
    )
    assert isinstance(emitted, Deps), emitted

    # Correct identification: same lin set as always-before Deps from apply.
    if 2 <= len(steps) <= 4:
        gold = _correct_plan_by_apply_replay(list(emitted.steps), mol, target)
        assert gold is not None
        assert emitted.same_linearizations(gold), (emitted, gold)

    # Must not be stricter than the stamped phase1 And/Seq branch.
    assert set(branch.iter_linearizations()) <= set(
        emitted.iter_linearizations()
    ), (branch, emitted)

    dh = steps[-1]
    for prep in steps[:-1]:
        if any(
            getattr(r, "added_by", None)
            and r.added_by[0] == prep.rule
            and set(r.added_by[1])
            & {int(x.origin) for x in prep.site if x.origin is not None}
            for r in dh.site
        ):
            for order in emitted.iter_linearizations():
                assert order.index(prep) < order.index(dh), (emitted, order)
    return emitted


def json_dumps_plan(plan: StepPlan) -> str:
    import json

    return json.dumps(plan.to_json(), sort_keys=True)


def _is_multistep_quinone_branch(branch: StepPlan) -> bool:
    """True for And(preps)→DH or prep→DH (len ≥ 2 with a final DH)."""
    steps = branch.steps
    if len(steps) < 2:
        return False
    return steps[-1].rule == "Dehydrogenation"


def _phase1_quinone_branches(mol, *, min_steps=2, max_steps=4):
    """Multi-step phase1 branches from ``QuinoneFormation.phase1_steps``.

    Uses the public phase1 API so QF → leaf Steps is what emission re-globs.
    """
    qf = QuinoneFormation()
    out = []
    try:
        stream = qf.metabolites(mol, attach_phase1_steps=True, tag_atoms=False)
    except Exception:
        return out
    seen = set()
    for site, _products in stream:
        atoms = site[1] if isinstance(site, tuple) else site
        public = qf.phase1_steps(mol, atoms)
        if not public:
            continue
        for branch in public.branches():
            n = len(branch)
            if not (min_steps <= n <= max_steps):
                continue
            key = json_dumps_plan(branch)
            if key in seen:
                continue
            ok_lin = None
            for lin in branch.linearizations():
                steps = list(lin.steps)
                if not _unique_steps(steps):
                    continue
                if _target_of(mol, steps) is None:
                    continue
                ok_lin = steps
                break
            if ok_lin is None:
                continue
            seen.add(key)
            out.append(
                {
                    "branch": branch,
                    "steps": ok_lin,
                    "site": frozenset(atoms),
                }
            )
    return out


@st.composite
def multistep_quinone_case(draw):
    """Draw a multi-step QF phase1 branch (2–4 leaf steps)."""
    smi = draw(st.sampled_from(_CORPUS))
    mol = _mol(smi)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)

    branches = [
        b
        for b in _phase1_quinone_branches(mol, min_steps=2, max_steps=4)
        if _is_multistep_quinone_branch(b["branch"])
    ]
    # Prefer length ≥ 3 (And(preps)→DH) when available.
    long = [b for b in branches if len(b["branch"]) >= 3]
    pool = long or branches
    assume(pool)
    item = draw(st.sampled_from(pool[:32]))
    branch = item["branch"]
    lin = draw(st.sampled_from(list(branch.linearizations())))
    steps = list(lin.steps)
    assume(_unique_steps(steps))
    mols = _mol_chain(mol, steps)
    assume(mols is not None)
    target = _canon(unmapped_smiles(mols[-1]))
    return smi, steps, target, mols, branch


@given(case=multistep_quinone_case())
@settings(
    max_examples=40,
    deadline=30_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_fuzz_multistep_quinone_phase1_and_cleave(case):
    """Multi-step QF phase1 → leaf steps → and_cleave_plan Deps allows phase1 lins."""
    smi, steps, target, mols, branch = case
    mol = _mol(smi)
    assert len(steps) >= 2
    assert steps[-1].rule == "Dehydrogenation"
    emitted = _assert_plan_matches_chemistry(
        mol, steps, target, mols=mols, exact=False
    )
    # Emitted Deps must match apply-replay gold (closure), not be stricter than phase1.
    assert _lin_set(branch) <= _lin_set(emitted), (branch, emitted)
    if isinstance(emitted, Deps) and 2 <= len(steps) <= 4:
        from xenosite._archive_forest.guided_path import _correct_plan_by_apply_replay

        gold = _correct_plan_by_apply_replay(list(emitted.steps), mol, target)
        assert gold is not None and emitted.same_linearizations(gold), (
            emitted,
            gold,
        )


def test_benzene_qf_phase1_and_oh_oh_dh_matches_emission():
    """Benzene para QF: emitted Deps ≡ phase1 And(OH,OH)→DH (closure)."""
    mol = _mol("c1ccccc1")
    branch = next(
        b
        for b in QuinoneFormation().phase1_steps(mol, frozenset({0, 3})).branches()
        if len(b) == 3
    )
    assert isinstance(branch.children[0], And)
    emitted = _assert_emitted_respects_phase1_deps(mol, branch)
    assert emitted.same_linearizations(_deps_from_plan(branch))


def test_anisole_qf_phase1_dealk_globbed_with_oh():
    """Anisole phase1 globs Dealkylation & Hydroxylation as And peers before DH."""
    mol = _mol("COc1ccc(O)cc1")
    branches = [
        b["branch"]
        for b in _phase1_quinone_branches(mol, min_steps=3, max_steps=3)
        if any(s.rule == "Dealkylation" for s in b["branch"].steps)
        and any(s.rule == "Hydroxylation" for s in b["branch"].steps)
    ]
    assert branches
    branch = branches[0]
    assert isinstance(branch.children[0], And)
    prep_rules = {s.rule for s in branch.children[0].steps}
    assert prep_rules == {"Dealkylation", "Hydroxylation"}
    emitted = _assert_emitted_respects_phase1_deps(mol, branch)
    assert isinstance(emitted, Deps)
    assert {s.rule for s in emitted.steps} >= prep_rules

    # Opposite walk order still allows every phase1 lin.
    other = list(branch.linearizations())
    assert len(other) >= 2
    steps = list(other[1].steps)
    mols = _mol_chain(mol, steps)
    assert mols is not None
    again = and_cleave_plan(steps, mols=mols)
    assert isinstance(again, Deps)
    assert again.same_linearizations(emitted)
    assert set(branch.iter_linearizations()) <= set(again.iter_linearizations())


def test_apap_qf_phase1_multistep_emission():
    """APAP multi-step QF phase1 branch emits exact lin set."""
    mol = _mol("CC(=O)Nc1ccc(O)cc1")
    branches = [
        b["branch"]
        for b in _phase1_quinone_branches(mol, min_steps=3, max_steps=3)
    ]
    assert branches
    _assert_emitted_respects_phase1_deps(mol, branches[0])


def test_phase1_steps_not_opaque_quinoneformation_composite():
    """Leaf steps are Phase I rules — never a QuinoneFormation Step."""
    mol = _mol("c1ccccc1")
    for item in _phase1_quinone_branches(mol, min_steps=2, max_steps=4):
        for step in item["branch"].steps:
            assert step.rule != "QuinoneFormation", item["branch"]
            assert step.rule in {
                "Hydroxylation",
                "Dehydrogenation",
                "Dealkylation",
                "OxidativeDehalogenation",
            }


def test_and_cleave_plan_dual_oh_and_then_dh_exact():
    """And(OH, OH) → DH linearizations hit BQ; other orders miss."""
    mol = _mol("c1ccccc1")
    steps = [
        Step("Hydroxylation", frozenset({0})),
        Step("Hydroxylation", frozenset({3})),
        Step(
            "Dehydrogenation",
            frozenset(
                {
                    AtomRef(added_by=("Hydroxylation", frozenset({0}))),
                    AtomRef(added_by=("Hydroxylation", frozenset({3}))),
                }
            ),
        ),
    ]
    target = _target_of(mol, steps)
    assert target == _canon("O=C1C=CC(=O)C=C1")
    _assert_plan_matches_chemistry(mol, steps, target, exact=True)
