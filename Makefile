# Local test entry points. Prefer these over ad-hoc cargo/pytest.
# Mirrors the language jobs in .github/workflows/test.yml.

.PHONY: help test test-python test-rust test-rust-python

PYTEST ?= uv run pytest
CARGO_TEST ?= cargo test -p xenosite-forest

help:
	@echo "test              Rust then Python (default)"
	@echo "test-rust         cargo test -p xenosite-forest"
	@echo "test-rust-python  cargo test -p xenosite-forest --features python"
	@echo "test-python       pytest tests/forest (-n auto, coverage)"

test: test-rust test-python

test-rust:
	$(CARGO_TEST)

test-rust-python:
	$(CARGO_TEST) --features python

test-python:
	$(PYTEST) tests/forest -n auto --cov=xenosite.forest --cov-report=term-missing
