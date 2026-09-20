"""Shared assertions for Phase1-equivalent StepPlan vs Forest emissions."""

from __future__ import annotations

from xenosite._archive_forest.rules import QuinoneFormation
from xenosite._archive_forest.utils import canon_smi, unmapped_smiles


def _canon(smi: str) -> str:
    return canon_smi(smi)


def quinone_plan_reaches_product(mol, plan, target_smi: str) -> str:
    """How ``plan`` reaches the quinone product SMILES.

    Prefer full Forest ``Linearization.apply`` when it fires (rare on aromatics:
    Forest ``Dehydrogenation`` often cannot form quinones). Otherwise apply prep
    steps only (``drop_last=1``) and require ``QuinoneFormation`` on a mid to
    emit ``target_smi`` — the main check that the plan's prep matches the
    quinone rule's product.

    Returns ``"full-forest"`` or ``"prep+quinone"``.
    Raises ``AssertionError`` if the product is not reached.
    """
    qf = QuinoneFormation()
    target = _canon(target_smi)

    for lin in plan.linearizations():
        out = lin.apply(mol)
        if out and target in {_canon(unmapped_smiles(m)) for m in out}:
            return "full-forest"

    reached = False
    for lin in plan.linearizations():
        mids = lin.apply(mol, drop_last=1)
        assert mids, "prep apply empty for %s" % (plan,)
        for mid in mids:
            for _site, products in qf.metabolize(mid, tag_atoms=False):
                if any(_canon(unmapped_smiles(p)) == target for p in products):
                    reached = True
                    break
            if reached:
                break
        if reached:
            break

    assert reached, (
        "quinone product %r not reached via plan %s (Forest DH apply empty; "
        "prep+QuinoneFormation also missed)" % (target_smi, plan)
    )
    return "prep+quinone"
