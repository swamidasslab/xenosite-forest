"""Session-wide mol edit guard so closed-mol mutations fail tests."""

import pytest

from xenosite.forest.edit_guard import edit_guard


@pytest.fixture(autouse=True)
def _mol_edit_guard():
    with edit_guard():
        yield
