"""Committed replay cases for multi-rule sequence parity.

Appended automatically when :mod:`test_rule_sequence_parity` discovers a new
failure (write the tuple; auto-commit unless ``CI`` /
``XENOSITE_SEQ_AUTOCOMMIT=0``). Each entry is
``(start_smiles, (rule_name, ...))`` — apply the named leaves in order and
assert Rust↔Python product bags match after every hop.
"""

from __future__ import annotations

SEQUENCE_PARITY_REPLAY: tuple[tuple[str, tuple[str, ...]], ...] = (
    ('CCO', ('Hydroxylation', 'Dehydrogenation',)),
    ('c1ccccc1', ('Hydroxylation', 'Dehydrogenation',)),
    ('COc1ccccc1', ('Dealkylation', 'Hydroxylation',)),
    ('CCN', ('NDealkylation', 'Hydroxylation', 'Dehydrogenation',)),
    ('CC(=O)Nc1ccc(O)cc1', ('Dehydrogenation', 'Dealkylation', 'Dehydration',)),
    ('CC(=O)Nc1ccc(O)cc1', ('Dehydration', 'Epoxidation', 'Dehydrogenation',)),
    ('O=C=Nc1ccccc1', ('Hydrogenation', 'Dehydration', 'Dehydrogenation',)),
)

