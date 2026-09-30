# Local test entry points. Prefer these over ad-hoc cargo/pytest.
# Mirrors the language jobs in .github/workflows/test.yml.

.PHONY: help test test-python test-python-smoke test-rust test-rust-smoke \
	test-rust-python test-smoke

PYTEST ?= uv run pytest
CARGO_TEST ?= cargo test -p xenosite-forest
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

# Full ``test-rust`` still runs ms1_apply_fuzz. Smoke omits that binary (heavy
# + anthraquinone Hydrogenation plan-replay xfail lives there).
RUST_SMOKE_INTEGRATION_TESTS = gold_csmi phase1_plan_fuzz sssom_coverage derisk

help:
	@echo "test               Rust then full Python (default)"
	@echo "test-rust          cargo test -p xenosite-forest"
	@echo "test-rust-smoke    opt-in: lib + integration except ms1_apply_fuzz"
	@echo "test-rust-python   cargo test -p xenosite-forest --features python"
	@echo "test-python        full pytest tests/forest (-n auto, coverage)"
	@echo "test-python-smoke  opt-in: rust wrapper + focused native↔Rust parity"
	@echo "test-smoke         opt-in: test-rust-smoke + test-rust-python + test-python-smoke"

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
	PYO3_PYTHON="$(PYO3_PYTHON)" $(CARGO_TEST) --lib --features python

test-python:
	$(PYTEST) tests/forest -n auto --cov=xenosite.forest --cov-report=term-missing

# Opt-in CI/main smoke. XENOSITE_PARITY_FULL=0 / --parity-focused shrinks
# rule_parity_fuzz to CoverIntent cases; default make test-python stays full.
test-python-smoke:
	XENOSITE_PARITY_FULL=0 $(PYTEST) $(PYTHON_SMOKE_PATHS) -n auto \
		--cov=xenosite.forest --cov-report=term-missing

test-smoke: test-rust-smoke test-rust-python test-python-smoke
