"""Committed replay cases for multi-rule sequence parity.

``SEQUENCE_PARITY_REPLAY`` — must stay green (regression lock).

``SEQUENCE_PARITY_XFAIL`` — appended automatically when random fuzz finds a
new product-bag gap (auto-commit unless ``CI`` / ``XENOSITE_SEQ_AUTOCOMMIT=0``).
Promote a case to ``SEQUENCE_PARITY_REPLAY`` once leaf parity makes it pass.
"""

from __future__ import annotations

SEQUENCE_PARITY_REPLAY: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("CCO", ("Hydroxylation", "Dehydrogenation")),
    ("c1ccccc1", ("Hydroxylation", "Dehydrogenation")),
    ("COc1ccccc1", ("Dealkylation", "Hydroxylation")),
    ("CCN", ("NDealkylation", "Hydroxylation", "Dehydrogenation")),
    ("CC(=O)Nc1ccc(O)cc1", ("Dehydrogenation", "Dealkylation", "Dehydration")),
    ("CC(=O)Nc1ccc(O)cc1", ("Dehydration", "Epoxidation", "Dehydrogenation")),
)

# Product-bag gaps still open on leaf parity — keep for replay / XPASS watch.
SEQUENCE_PARITY_XFAIL: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("O=C=Nc1ccccc1", ("Hydrogenation", "Dehydration", "Dehydrogenation")),
)
