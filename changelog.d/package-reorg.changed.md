Split the forest package into a Rust public stub (``xenosite.forest``),
``native`` (RDKit reference), and ``legacy`` (frozen 0.6.x). Default install no
longer requires RDKit; use the ``rdkit`` extra for native/legacy. Public stub
exports only ``find_path``, ``available``, ``PhaseOne``, ``Epoxidation``,
``QuinoneFormation``, ``EpoxideOpening``, and ``NDealkylation``. See
``docs/forest/MIGRATING_0.8.md``.
