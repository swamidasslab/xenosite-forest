#!/usr/bin/env bash
# Sparse-checkout chematic @v1.0.21 and apply our tiny local patch.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
sub="$root/vendor/chematic"
patch="$root/patches/chematic-v1.0.21-atom-tag-visit-order.patch"
pin="v1.0.21"

if [[ ! -e "$sub/.git" && ! -f "$sub/.git" ]]; then
  git -C "$root" submodule update --init --depth 1 vendor/chematic
fi

git -C "$sub" fetch --tags --depth 1 origin tag "$pin" 2>/dev/null || true
git -C "$sub" checkout -q "$pin"

git -C "$sub" sparse-checkout init --cone
git -C "$sub" sparse-checkout set \
  crates/chematic \
  crates/chematic-core \
  crates/chematic-smiles \
  crates/chematic-smarts \
  crates/chematic-rxn \
  crates/chematic-perception

# Drop any previous apply, then re-apply so the script is idempotent.
git -C "$sub" checkout -q -- .
git -C "$sub" clean -fdq
git -C "$sub" apply "$patch"

echo "chematic $pin sparse + patch applied -> $sub"
