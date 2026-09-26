"""Shared HTTP + candidate record helpers (stdlib only)."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Iterable


USER_AGENT = "xenosite-xrm-harvest/0.1 (+https://github.com/swamidasslab/xenosite-forest)"


def http_get(url: str, timeout: float = 45.0, retries: int = 3) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:  # noqa: BLE001 — offline harvester; surface later
            last = e
            time.sleep(1.5 * (attempt + 1))
    assert last is not None
    raise last


def http_get_json(url: str, **kwargs: Any) -> Any:
    return json.loads(http_get(url, **kwargs).decode("utf-8", errors="replace"))


def http_get_text(url: str, **kwargs: Any) -> str:
    return http_get(url, **kwargs).decode("utf-8", errors="replace")


def candidate(
    *,
    candidate_id: str,
    pref_label: str,
    synonyms: Iterable[str] | None = None,
    definition: str | None = None,
    sources: list[dict[str, str]] | None = None,
    examples: list[dict[str, Any]] | None = None,
    suggested_spines: list[str] | None = None,
    notes: str | None = None,
) -> dict[str, Any]:
    syn = sorted({s.strip() for s in (synonyms or []) if s and s.strip() and s.strip() != pref_label})
    return {
        "candidate_id": candidate_id,
        "pref_label": pref_label,
        "synonyms": syn,
        "definition": definition,
        "sources": sources or [],
        "examples": examples or [],
        "suggested_spines": suggested_spines or [],
        "notes": notes,
    }


def merge_candidates(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Merge by lowercased pref_label; union synonyms/examples/sources."""
    by: dict[str, dict[str, Any]] = {}
    for row in rows:
        key = row["pref_label"].strip().lower()
        if key not in by:
            by[key] = dict(row)
            by[key]["synonyms"] = list(row.get("synonyms") or [])
            by[key]["sources"] = list(row.get("sources") or [])
            by[key]["examples"] = list(row.get("examples") or [])
            by[key]["suggested_spines"] = list(row.get("suggested_spines") or [])
            continue
        cur = by[key]
        cur["synonyms"] = sorted(
            set(cur.get("synonyms") or []) | set(row.get("synonyms") or [])
        )
        # de-dupe sources/examples by json
        for field in ("sources", "examples"):
            seen = {json.dumps(x, sort_keys=True) for x in cur.get(field) or []}
            for x in row.get(field) or []:
                k = json.dumps(x, sort_keys=True)
                if k not in seen:
                    cur[field].append(x)
                    seen.add(k)
        spines = set(cur.get("suggested_spines") or []) | set(row.get("suggested_spines") or [])
        cur["suggested_spines"] = sorted(spines)
        if not cur.get("definition") and row.get("definition"):
            cur["definition"] = row["definition"]
        if row.get("notes"):
            cur["notes"] = ((cur.get("notes") or "") + " | " + row["notes"]).strip(" |")
    return [by[k] for k in sorted(by)]
