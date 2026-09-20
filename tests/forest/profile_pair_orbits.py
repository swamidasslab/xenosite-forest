#!/usr/bin/env python3
"""cProfile harness for pair-orbit recipes (isotope vs pynauty, cache paths).

Writes durable text under artifacts/pair_orbit_profile.out (+ .pstats).

  uv run python tests/forest/profile_pair_orbits.py
"""

from __future__ import annotations

import cProfile
import importlib.util
import io
import pstats
import time
from pathlib import Path

from xenosite.forest import graph_isomorphism as gi
from xenosite.forest.rdkitutil import MolFromSmiles

ROOT = Path(__file__).resolve().parents[2]
ART_OUT = ROOT / "artifacts" / "pair_orbit_profile.out"
ART_PSTATS = ROOT / "artifacts" / "pair_orbit_profile.pstats"

# Multi-symmetric aromatics + larger drug-like / H2H / substrate-library mols.
CASES: list[tuple[str, str]] = [
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("phenol", "Oc1ccccc1"),
    ("hydroquinone", "Oc1ccc(O)cc1"),
    ("MeOPhOH", "COc1ccc(O)cc1"),
    ("anisole", "COc1ccccc1"),
    ("benzodioxole", "c1ccc2c(c1)OCO2"),
    ("biphenyl-ish", "c1ccccc1CCCCc2ccccc2"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("diphenhydramine-like", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
    ("chloramphenicol-like", "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl"),
]

TOP_N = 25
HAS_PYNAUTY = importlib.util.find_spec("pynauty") is not None


def _log(msg: str) -> None:
    print(msg, flush=True)


def _mol(smiles: str):
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return mol


def _clear_pair_caches(mol) -> None:
    if not getattr(mol, "_forest", None):
        return
    structure = mol._forest.setdefault("cache", {})
    structure.pop("site_pair_orbits_nauty", None)
    structure.pop("site_pair_orbits_smiles", None)
    structure.pop("bond_topeqiv", None)


def run_isotope_batch(mol) -> None:
    gi.site_pair_orbits_smiles(mol)


def run_nauty_batch(mol) -> None:
    gi.site_pair_orbits_nauty(mol)


def run_signature_walk(mol) -> None:
    """Per-pair signature lookups (materialize-on-demand / forest cache path)."""

    n = mol.GetNumAtoms()
    for i in range(n):
        for j in range(i + 1, n):
            gi.atom_pair_orbit_key(mol, frozenset({i, j}))
    nb = mol.GetNumBonds()
    for i in range(nb):
        for j in range(i + 1, nb):
            gi.bond_pair_orbit_key(mol, frozenset({i, j}))
        for a in range(n):
            gi.bond_atom_orbit_key(mol, i, a)


def run_named_single_pair(mol) -> None:
    """Named single-pair helpers (profiling entry points)."""

    n = min(mol.GetNumAtoms(), 8)
    for i in range(n):
        for j in range(i + 1, n):
            gi.atom_pair_orbit_isotope(mol, i, j)
            if HAS_PYNAUTY:
                gi.atom_pair_orbit_pynauty(mol, i, j)
    nb = min(mol.GetNumBonds(), 8)
    for i in range(nb):
        for j in range(i + 1, nb):
            gi.bond_pair_orbit_isotope(mol, i, j)
            if HAS_PYNAUTY:
                gi.bond_pair_orbit_pynauty(mol, i, j)
        for a in range(min(n, 6)):
            gi.bond_atom_orbit_isotope(mol, i, a)
            if HAS_PYNAUTY:
                gi.bond_atom_orbit_pynauty(mol, i, a)


def main() -> int:
    ART_OUT.parent.mkdir(parents=True, exist_ok=True)
    _log(f"pair_orbit profile  pynauty={HAS_PYNAUTY}  cases={len(CASES)}")
    _log(f"writing {ART_OUT}")

    case_lines: list[str] = []
    wall0 = time.perf_counter()
    profiler = cProfile.Profile()
    profiler.enable()

    for name, smiles in CASES:
        mol = _mol(smiles)
        _log(f"  case {name} atoms={mol.GetNumAtoms()} bonds={mol.GetNumBonds()}")
        t0 = time.perf_counter()
        _clear_pair_caches(mol)
        run_isotope_batch(mol)
        t_iso = time.perf_counter() - t0

        t0 = time.perf_counter()
        _clear_pair_caches(mol)
        if HAS_PYNAUTY:
            run_nauty_batch(mol)
        t_nauty = time.perf_counter() - t0

        t0 = time.perf_counter()
        _clear_pair_caches(mol)
        run_signature_walk(mol)
        t_sig = time.perf_counter() - t0

        t0 = time.perf_counter()
        _clear_pair_caches(mol)
        run_named_single_pair(mol)
        t_named = time.perf_counter() - t0

        line = (
            f"  {name}: isotope_batch={t_iso:.4f}s nauty_batch={t_nauty:.4f}s "
            f"signature_walk={t_sig:.4f}s named_singles={t_named:.4f}s"
        )
        case_lines.append(line)
        _log(f"   {line.strip()}")

    profiler.disable()
    wall = time.perf_counter() - wall0

    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream)
    stats.sort_stats("tottime")
    stats.print_stats(TOP_N)

    focus = io.StringIO()
    stats_focus = pstats.Stats(profiler, stream=focus)
    stats_focus.sort_stats("tottime")
    stats_focus.print_stats(
        "graph_isomorphism|site_pair_orbits|atom_pair_orbit|bond_pair_orbit|"
        "bond_atom_orbit|marked_site_pair"
    )

    body = (
        f"# pair_orbit profile\n"
        f"# pynauty={HAS_PYNAUTY}\n"
        f"# wall={wall:.4f}s cases={len(CASES)}\n\n"
        f"## Per-case wall\n"
        + "\n".join(case_lines)
        + f"\n\n## cProfile tottime top {TOP_N}\n"
        f"{stream.getvalue()}\n"
        f"## Focus: graph_isomorphism named functions\n"
        f"{focus.getvalue()}\n"
    )
    ART_OUT.write_text(body)
    profiler.dump_stats(ART_PSTATS)
    _log(f"done wall={wall:.4f}s -> {ART_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
