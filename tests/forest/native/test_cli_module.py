"""Smoke: native CLI module is importable and exposes ``main``."""

from __future__ import annotations

from xenosite.forest.native import bfs, dfs
from xenosite.forest.native.cli import main
from xenosite.forest.native.find_path import bfs as find_path_bfs


def test_cli_module_wires_bfs():
    assert callable(main)
    assert bfs is find_path_bfs
    assert callable(dfs)
