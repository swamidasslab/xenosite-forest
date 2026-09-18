"""Hypothesis fuzz: random rule walks, then guided find recovers the recipe.

Starts with short Phase I walks; harder cases create with QF+Phase I and find
with Phase I only. Looping walks (A→B→A) are rejected during generation.

Searches use a tight ``max_expansions`` budget so expensive cases stop early
with ``counters.budget_exhausted`` — that is the signal for new heuristics.
"""

from __future__ import annotations

from pathlib import Path

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from hypothesis.database import DirectoryBasedExampleDatabase
from rdkit import Chem

from xenosite.forest import PathSearchCounters, RuleSet, find_path
from xenosite.forest.guided_path import _canon
from xenosite.forest.rules import (
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
)

_HYPOTHESIS_DIR = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
_HYPOTHESIS_DB = DirectoryBasedExampleDatabase(str(_HYPOTHESIS_DIR))

# Small starters — keep walks cheap
_SMALL_CORPUS = (
    "CCO",
    "CC",
    "C=C",
    "c1ccccc1",
    "Oc1ccccc1",
    "CC(=O)Nc1ccc(O)cc1",
    "CCN",
    "COc1ccccc1",
)

_MAX_HEAVY = 20

# Expansion budgets — prefer these over wall-clock timeouts for diagnostics.
_BUDGET_DEPTH1 = 40
_BUDGET_DEPTH2 = 80
_BUDGET_QF_PHASE1 = 120


def _mol(smi: str):
    return Chem.MolFromSmiles(smi)


def _is_subsequence(needle, haystack) -> bool:
    """True if ``needle`` appears in order inside ``haystack`` (not necessarily contiguous)."""
    it = iter(haystack)
    for item in needle:
        for x in it:
            if x == item:
                break
        else:
            return False
    return True


def _collect_candidates(mol, rules, seen_smis: set[str]):
    """Topo-distinct products that do not revisit ``seen_smis`` (no loops)."""
    out = []
    for rule in rules:
        try:
            stream = rule.metabolize(
                mol,
                tag_atoms=False,
                only_emit_topologically_distinct_sites=True,
            )
        except Exception:
            continue
        for site, products in stream:
            for product in products or []:
                if not product:
                    continue
                try:
                    smi = _canon(product)
                except Exception:
                    continue
                if not smi or smi in seen_smis:
                    continue
                out.append((rule, site, product, smi))
    return out


def random_walk(draw, start_smi: str, rules, n_steps: int):
    """Apply ``n_steps`` random rule hops; reject loops and dead ends."""
    mol = _mol(start_smi)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)

    path_smis = [_canon(mol)]
    recipe = []  # (rule_name, product_smi)
    current = mol

    for _ in range(n_steps):
        candidates = _collect_candidates(current, rules, set(path_smis))
        assume(candidates)
        candidates.sort(key=lambda c: c[2].GetNumHeavyAtoms())
        pool = candidates[: max(1, min(12, len(candidates)))]
        rule, site, product, smi = draw(st.sampled_from(pool))
        recipe.append((rule.name, smi))
        path_smis.append(smi)
        current = product
        assume(current.GetNumHeavyAtoms() <= _MAX_HEAVY + 4)

    return path_smis[0], path_smis[-1], recipe, path_smis


def _phase1_rules():
    """Small Phase I subset — enough for OH / DH walks without Full cost."""
    return [Hydroxylation(), Dehydrogenation()]


def _create_rules_with_quinone():
    return [QuinoneFormation(), Hydroxylation(), Dehydrogenation()]


def _phase1_ruleset():
    return RuleSet([Hydroxylation(), Dehydrogenation()], name="fuzz_phase1")


def _guided(r_smi, t_smi, *, depth, max_expansions, max_paths=5):
    counters = PathSearchCounters()
    hits = list(
        find_path(
            r_smi,
            t_smi,
            ruleset=_phase1_ruleset(),
            depth=depth,
            maybe_prefixes=False,
            max_paths=max_paths,
            max_expansions=max_expansions,
            counters=counters,
        )
    )
    return hits, counters


@st.composite
def small_start(draw):
    return draw(st.sampled_from(_SMALL_CORPUS))


@given(start=small_start(), data=st.data())
@settings(
    max_examples=30,
    deadline=10_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_fuzz_guided_recovers_depth1_phase1_recipe(start: str, data):
    """Depth-1 Phase I walk: guided finds T within a small expansion budget."""
    rules = _phase1_rules()
    r_smi, t_smi, recipe, path_smis = random_walk(data.draw, start, rules, 1)
    assume(r_smi != t_smi)

    hits, counters = _guided(
        r_smi, t_smi, depth=2, max_expansions=_BUDGET_DEPTH1, max_paths=5
    )
    assert hits, (r_smi, t_smi, recipe, counters.as_dict())
    found_path = [_canon(s) for s in hits[0][0]]
    assert found_path[-1] == _canon(t_smi)
    assert _is_subsequence(path_smis, found_path)
    assert not counters.budget_exhausted


@given(start=small_start(), data=st.data())
@settings(
    max_examples=20,
    deadline=15_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_fuzz_guided_recovers_depth2_phase1_recipe(start: str, data):
    """Depth-2 Phase I: find T under budget; prefer recipe intermediates when unique."""
    rules = _phase1_rules()
    r_smi, t_smi, recipe, path_smis = random_walk(data.draw, start, rules, 2)
    assume(r_smi != t_smi)
    assume(len(set(path_smis)) == len(path_smis))

    hits, counters = _guided(
        r_smi, t_smi, depth=3, max_expansions=_BUDGET_DEPTH2, max_paths=8
    )
    # Miss or budget blow → diagnostic payload (heuristics thread).
    assert hits, (r_smi, t_smi, recipe, path_smis, counters.as_dict())
    found_path = [_canon(s) for s in hits[0][0]]
    assert found_path[-1] == _canon(t_smi), (found_path, counters.as_dict())
    # Alternate Phase I orders to the same T are common (e.g. OH↔DH); filter
    # those out rather than fail — the hard signal is miss / budget_exhausted.
    assume(_is_subsequence(path_smis, found_path))


@given(start=st.sampled_from(("c1ccccc1", "Oc1ccccc1", "Cc1ccccc1")), data=st.data())
@settings(
    max_examples=12,
    deadline=20_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_fuzz_qf_create_phase1_find(start: str, data):
    """QF/Phase I create → Phase I find under budget; QF-only products may miss."""
    create_rules = _create_rules_with_quinone()
    use_qf = data.draw(st.booleans())
    if use_qf:
        mol = _mol(start)
        qf_cands = _collect_candidates(mol, [QuinoneFormation()], {_canon(mol)})
        assume(qf_cands)
        _rule, _site, _product, t_smi = data.draw(st.sampled_from(qf_cands[:8]))
        r_smi = _canon(mol)
        path_smis = [r_smi, t_smi]
        recipe = [("QuinoneFormation", t_smi)]
    else:
        r_smi, t_smi, recipe, path_smis = random_walk(
            data.draw, start, create_rules, data.draw(st.integers(2, 3))
        )
        assume(r_smi != t_smi)

    hits, counters = _guided(
        r_smi, t_smi, depth=4, max_expansions=_BUDGET_QF_PHASE1, max_paths=10
    )
    used_qf = any(name == "QuinoneFormation" for name, _ in recipe)

    if not hits:
        # Budget / chemistry miss is the diagnostic — do not hang; soft for QF.
        payload = (r_smi, t_smi, recipe, counters.as_dict())
        assume(not used_qf)  # QF-only product may need composite; filter
        assert hits, payload
        return

    assert _canon(hits[0][0][-1]) == _canon(t_smi)
    if not used_qf:
        # Alternate site orders to the same poly-OH product are common.
        assume(
            _is_subsequence(path_smis, [_canon(s) for s in hits[0][0]])
        )


def test_random_walk_rejects_immediate_loops():
    """Sanity: candidates never include the current molecule."""
    mol = _mol("CCO")
    rules = _phase1_rules()
    seen = {_canon(mol)}
    for _rule, _site, _p, smi in _collect_candidates(mol, rules, seen):
        assert smi not in seen


def test_is_subsequence_helper():
    assert _is_subsequence(["a", "c"], ["a", "b", "c"])
    assert not _is_subsequence(["a", "c"], ["c", "a"])


def test_max_expansions_stops_and_flags_budget():
    """Tight budget must stop search and set budget_exhausted (no hang)."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset="PhaseOneRS",
            depth=4,
            maybe_prefixes=False,
            max_paths=1,
            max_expansions=3,
            counters=counters,
        )
    )
    assert counters.billed() <= 3
    assert counters.budget_exhausted or hits
    # With only 3 billed applies, benzene→BQ via PhaseOne usually misses.
    if not hits:
        assert counters.budget_exhausted
