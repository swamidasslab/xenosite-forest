#!/usr/bin/env python3
"""Set package versions from a release tag (vX.Y.Z) or bare X.Y.Z.

Rewrites:
  - crates/xenosite-forest/Cargo.toml  (maturin / PyPI source of truth)
  - Cargo.lock                         (workspace package entry)
  - js/package.json                    (kept in sync for the next tag)

Usage:
  python scripts/set_release_version.py           # reads GITHUB_REF_NAME
  python scripts/set_release_version.py v0.9.0
  python scripts/set_release_version.py 0.9.0
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CARGO_TOML = ROOT / "crates" / "xenosite-forest" / "Cargo.toml"
CARGO_LOCK = ROOT / "Cargo.lock"
JS_PACKAGE = ROOT / "js" / "package.json"

# Match release tags we publish; allow optional pre-release / build suffix.
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+(?:[.-][0-9A-Za-z.-]+)?$")


def normalize(raw: str) -> str:
    ver = raw.strip()
    if ver.startswith("v") or ver.startswith("V"):
        ver = ver[1:]
    if not VERSION_RE.fullmatch(ver):
        raise SystemExit(f"bad version {raw!r}; expected vX.Y.Z or X.Y.Z")
    return ver


def set_cargo_toml(ver: str) -> None:
    text = CARGO_TOML.read_text(encoding="utf-8")
    new, n = re.subn(
        r'(?m)^version = "[^"]*"',
        f'version = "{ver}"',
        text,
        count=1,
    )
    if n != 1:
        raise SystemExit(f"failed to rewrite version in {CARGO_TOML}")
    CARGO_TOML.write_text(new, encoding="utf-8")


def set_cargo_lock(ver: str) -> None:
    text = CARGO_LOCK.read_text(encoding="utf-8")
    new, n = re.subn(
        r'(name = "xenosite-forest"\n)version = "[^"]*"',
        rf'\1version = "{ver}"',
        text,
        count=1,
    )
    if n != 1:
        raise SystemExit(f"failed to rewrite xenosite-forest version in {CARGO_LOCK}")
    CARGO_LOCK.write_text(new, encoding="utf-8")


def set_js_package(ver: str) -> None:
    data = json.loads(JS_PACKAGE.read_text(encoding="utf-8"))
    data["version"] = ver
    JS_PACKAGE.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str]) -> None:
    if len(argv) > 1:
        raw = argv[1]
    else:
        raw = os.environ.get("GITHUB_REF_NAME", "")
        if not raw:
            raise SystemExit(
                "usage: set_release_version.py <vX.Y.Z>  "
                "(or set GITHUB_REF_NAME)"
            )
    ver = normalize(raw)
    set_cargo_toml(ver)
    set_cargo_lock(ver)
    set_js_package(ver)
    print(f"set release version {ver}")


if __name__ == "__main__":
    main(sys.argv)
