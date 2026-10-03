# Local test + lint entry points. Prefer these over ad-hoc cargo/pytest.
# Mirrors the language jobs in .github/workflows/test.yml.

.PHONY: help stubs check-stubs check check-rust check-python test test-python test-python-smoke \
	test-python-legacy test-rust test-rust-smoke test-rust-python test-smoke

PYTEST ?= uv run pytest
# Same as CI: PATH ``cargo`` after dtolnay/rust-toolchain / rust-toolchain.toml.
# Do not rewrite to an absolute rustup path — that diverges from CI/CD.
CARGO ?= cargo
CARGO_TEST ?= $(CARGO) test -p xenosite-forest
# PyO3 links libpython for ``--features python`` unit tests. Prefer the project
# uv interpreter locally; on CI (no uv in the rust job) fall back to python3
# from setup-python / PATH so PYO3_PYTHON is never empty.
PYO3_PYTHON ?= $(shell \
	(command -v uv >/dev/null 2>&1 && uv run python -c 'import sys; print(sys.executable)') \
	|| python3 -c 'import sys; print(sys.executable)')

# Product-door smoke paths (Python). Opt-in only — not wired into ``test`` /
# ``test-python``. Skips native-only and legacy archives; keeps Rust wrapper
# tests plus a thin native↔Rust parity gate (focused corpus + coverage pool).
PYTHON_SMOKE_PATHS = \
	tests/forest/rust \
	tests/forest/native/test_rule_parity_pairs.py \
	tests/forest/native/test_rule_parity_fuzz.py \
	tests/forest/native/test_rust_parity_coverage.py

# Frozen 0.6.x archive + native↔legacy product pairs. Opt-in only — not CI,
# not ``make test`` / ``make test-python``.
PYTHON_LEGACY_PATHS = \
	tests/forest/legacy \
	tests/forest/native/test_parity.py

# Full ``test-rust`` still runs ms1_apply_fuzz. Smoke omits that binary (heavy
# + anthraquinone Hydrogenation plan-replay xfail lives there).
RUST_SMOKE_INTEGRATION_TESTS = gold_csmi phase1_plan_fuzz sssom_coverage derisk

help:
	@echo "stubs              regenerate src/xenosite/forest/_rust.pyi from Rust"
	@echo "check-stubs        fail if _rust.pyi is stale (CI gate)"
	@echo "check              CI non-test gates (fmt/clippy/WASM/stubs/ruff/pyright)"
	@echo "check-rust         rustfmt + clippy (+ python feature) + WASM build"
	@echo "check-python       ruff + pyright (forest package + tests)"
	@echo "test               Rust then Python (no legacy)"
	@echo "test-rust          cargo test -p xenosite-forest"
	@echo "test-rust-smoke    opt-in: lib + integration except ms1_apply_fuzz"
	@echo "test-rust-python   cargo test -p xenosite-forest --lib --features stubs"
	@echo "test-python        pytest tests/forest minus legacy / native↔legacy"
	@echo "test-python-smoke  opt-in: rust wrapper + focused native↔Rust parity"
	@echo "test-python-legacy opt-in: tests/forest/legacy + native test_parity"
	@echo "test-smoke         opt-in: test-rust-smoke + test-rust-python + test-python-smoke"

# CI lint/compile gates that ``make test*`` does not run (see test.yml rust +
# lint jobs). Needs ``wasm32-unknown-unknown`` on the active toolchain.
# If Homebrew cargo/rustc is ahead of rustup on PATH, put
# ``$$(dirname $$(rustup which rustc))`` first or WASM will fail to find core.
check: check-rust check-stubs check-python

# ``_rust.pyi`` is generated (pyo3-stub-gen + ``typed_dict!``); never hand-edit.
# ``stubs`` is a dev-only feature: wheels build ``python,extension-module``.
# ``test-rust-python`` also fails on a stale stub (checked_in_stub_is_current).
STUB_GEN = PYO3_PYTHON="$(PYO3_PYTHON)" $(CARGO) run -q -p xenosite-forest \
	--features stubs --bin stub_gen --

stubs:
	$(STUB_GEN)

check-stubs:
	$(STUB_GEN) --check

check-rust:
	$(CARGO) fmt --all -- --check
	$(CARGO) clippy -p xenosite-forest --all-targets --no-deps -- -D warnings
	$(CARGO) clippy -p xenosite-forest --all-targets --features python --no-deps -- -D warnings
	$(CARGO) clippy -p xenosite-forest --all-targets --features stubs --no-deps -- -D warnings
	$(CARGO) build -p xenosite-forest --target wasm32-unknown-unknown --features wasm

check-python:
	uv run ruff check src/xenosite/forest tests/forest
	uv run pyright src/xenosite/forest

test: test-rust test-python

test-rust:
	$(CARGO_TEST)

test-rust-smoke:
	$(CARGO_TEST) --lib
	@for t in $(RUST_SMOKE_INTEGRATION_TESTS); do \
		$(CARGO_TEST) --test $$t || exit 1; \
	done

test-rust-python:
	# ``--lib``: PyO3 unit tests live in the crate; integration binaries need
	# not link libpython. PYO3_PYTHON forces a consistent interpreter/dylib.
	# ``stubs`` ⊃ ``python``: also runs the checked_in_stub_is_current gate.
	PYO3_PYTHON="$(PYO3_PYTHON)" $(CARGO_TEST) --lib --features stubs

# Default Python suite: product door + native RDKit reference. Excludes frozen
# legacy archive and native↔legacy ``test_parity`` (use ``test-python-legacy``).
test-python:
	$(PYTEST) tests/forest -n auto \
		--ignore=tests/forest/legacy \
		--ignore=tests/forest/native/test_parity.py \
		--cov=xenosite.forest --cov-report=term-missing

# Opt-in CI/main smoke. XENOSITE_PARITY_FULL=0 / --parity-focused shrinks
# rule_parity_fuzz to CoverIntent cases; default make test-python stays full
# native↔Rust cartesian (still without legacy).
test-python-smoke:
	XENOSITE_PARITY_FULL=0 $(PYTEST) $(PYTHON_SMOKE_PATHS) -n auto \
		--parity-focused \
		--cov=xenosite.forest --cov-report=term-missing

# Opt-in: frozen legacy suite + native↔legacy product pairs. Not on CI.
test-python-legacy:
	$(PYTEST) $(PYTHON_LEGACY_PATHS) -n auto \
		--cov=xenosite.forest --cov-report=term-missing

test-smoke: test-rust-smoke test-rust-python test-python-smoke
