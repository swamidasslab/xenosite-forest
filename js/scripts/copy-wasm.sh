#!/usr/bin/env bash
# Copy wasm-pack artifacts into dist/ for the published package.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/dist/wasm"
cp -R "$ROOT/src/wasm/." "$ROOT/dist/wasm/"
# Drop wasm-pack's stub package.json if present.
rm -f "$ROOT/dist/wasm/package.json" "$ROOT/dist/wasm/README.md" "$ROOT/dist/wasm/.gitignore"
echo "Copied wasm → dist/wasm"
