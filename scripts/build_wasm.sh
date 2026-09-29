#!/usr/bin/env bash
# Build wasm-bindgen package into js/src/wasm (browser + Node via fs init).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

export RUSTUP_TOOLCHAIN="${RUSTUP_TOOLCHAIN:-1.88.0}"
TC="$HOME/.rustup/toolchains/${RUSTUP_TOOLCHAIN}-aarch64-apple-darwin/bin"

# Toolchain cargo/rustc first; keep Homebrew rustup + ~/.cargo (wasm-pack) available.
export PATH="${TC}:/opt/homebrew/bin:${HOME}/.cargo/bin:${PATH}"

mkdir -p js/src/wasm
echo "==> wasm-pack → js/src/wasm (toolchain=$RUSTUP_TOOLCHAIN)"
echo "    rustc=$(command -v rustc) $(rustc --version)"
echo "    cargo=$(command -v cargo) $(cargo --version)"
echo "    rustup=$(command -v rustup)"

if ! rustup target list --installed --toolchain "$RUSTUP_TOOLCHAIN" | grep -q '^wasm32-unknown-unknown$'; then
  rustup target add wasm32-unknown-unknown --toolchain "$RUSTUP_TOOLCHAIN"
fi

wasm-pack build crates/xenosite-forest \
  --target web \
  --out-dir "$ROOT/js/src/wasm" \
  --out-name xenosite_forest \
  --features wasm

rm -f js/src/wasm/.gitignore js/src/wasm/package.json js/src/wasm/README.md
echo "Done: js/src/wasm"
