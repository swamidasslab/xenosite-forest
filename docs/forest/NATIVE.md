# Native RDKit engine (feature freeze)

`xenosite.forest.native` is the **RDKit Python reference** engine. It is kept for
parity checks, benches, and historical validation. It is **not** where new
product features land.

## Where new work goes

| Layer | Path | Role |
|-------|------|------|
| Rust chemistry + search | `crates/xenosite-forest/` | Product implementation |
| Python public wrapper | `src/xenosite/forest/` (`__init__.py`, `_ext.py`) | Thin PyO3 door (`xenosite.forest`) |
| Native RDKit | `src/xenosite/forest/native/` | **Frozen for new features** — reference / validation only |
| Legacy 0.6.x | `src/xenosite/forest/legacy/` | Frozen archive — see [LEGACY.md](LEGACY.md) |

Importing `xenosite.forest.native` emits a `UserWarning` stating this policy.

## Allowed native edits

Do **not** add APIs, rules, search modes, or sampling helpers to `native/` to
match Rust. Allowed changes are narrow:

- Bugfixes that keep reference parity / CI green
- Docstring or warning text that reinforces this freeze
- Explicit, recorded exceptions (note them here with `Status: approved`)

Anything else belongs in the Rust crate and the `xenosite.forest` wrapper.
See [MIGRATING_0.8.md](MIGRATING_0.8.md) and [RUST.md](RUST.md).
