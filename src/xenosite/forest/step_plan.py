"""Partial-order plans over reaction steps (domain-agnostic)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterator, Sequence


@dataclass(frozen=True)
class Step:
    """One named reaction application at a site of 0-based atom indices."""

    rule: str
    site: frozenset[int]

    def __post_init__(self):
        if not isinstance(self.site, frozenset):
            object.__setattr__(self, "site", frozenset(self.site))


class StepPlan:
    """Partial order over :class:`Step` nodes.

    Edges in ``precedes`` are pairs of indices ``(i, j)`` meaning
    ``steps[i]`` must occur before ``steps[j]``. Within a layer from
    :meth:`layers`, order is undefined.
    """

    __slots__ = ("_steps", "_precedes")

    def __init__(
        self,
        steps: Sequence[Step],
        precedes: Sequence[tuple[int, int]] = (),
    ):
        self._steps = tuple(steps)
        self._precedes = tuple((int(a), int(b)) for a, b in precedes)
        n = len(self._steps)
        for a, b in self._precedes:
            if not (0 <= a < n and 0 <= b < n):
                raise ValueError("precedes index out of range")
            if a == b:
                raise ValueError("precedes cannot be reflexive")

    @classmethod
    def singleton(cls, rule: str, site) -> StepPlan:
        """One step, no ordering constraints."""
        return cls((Step(rule, frozenset(site)),), ())

    @classmethod
    def layers(cls, layers: Sequence[Sequence[Step]]) -> StepPlan:
        """Stack layers: every node in layer k precedes every node in layer k+1.

        Order within a layer is undefined (no edges).
        """
        steps: list[Step] = []
        ranges: list[tuple[int, int]] = []
        for layer in layers:
            layer = list(layer)
            if not layer:
                continue
            start = len(steps)
            steps.extend(layer)
            ranges.append((start, len(steps)))
        precedes: list[tuple[int, int]] = []
        for (a0, a1), (b0, b1) in zip(ranges, ranges[1:]):
            for i in range(a0, a1):
                for j in range(b0, b1):
                    precedes.append((i, j))
        return cls(steps, precedes)

    def __len__(self) -> int:
        return len(self._steps)

    def __eq__(self, other) -> bool:
        if not isinstance(other, StepPlan):
            return NotImplemented
        return self._steps == other._steps and self._precedes == other._precedes

    def __hash__(self) -> int:
        return hash((self._steps, self._precedes))

    def __repr__(self) -> str:
        return "StepPlan(steps=%r, precedes=%r)" % (self._steps, self._precedes)

    @property
    def steps(self) -> tuple[Step, ...]:
        return self._steps

    @property
    def precedes(self) -> tuple[tuple[int, int], ...]:
        return self._precedes

    def to_json(self) -> dict:
        return {
            "steps": [
                {"rule": s.rule, "site": sorted(s.site)} for s in self._steps
            ],
            "precedes": [list(p) for p in self._precedes],
        }

    @classmethod
    def from_json(cls, data: dict) -> StepPlan:
        steps = [
            Step(item["rule"], frozenset(item["site"])) for item in data["steps"]
        ]
        precedes = [tuple(p) for p in data.get("precedes", ())]
        return cls(steps, precedes)

    @classmethod
    def from_mol(cls, mol) -> StepPlan:
        """Read the ``phase1_steps`` mol property."""
        if not mol.HasProp("phase1_steps"):
            raise ValueError("mol has no phase1_steps property")
        return cls.from_json(json.loads(mol.GetProp("phase1_steps")))

    def attach_to_mol(self, mol) -> None:
        """Stamp ``phase1_steps`` JSON onto ``mol``."""
        mol.SetProp("phase1_steps", json.dumps(self.to_json(), separators=(",", ":")))

    def iter_linearizations(self) -> Iterator[tuple[Step, ...]]:
        """Lazily yield every total order consistent with ``precedes``."""
        n = len(self._steps)
        if n == 0:
            yield ()
            return

        successors = [[] for _ in range(n)]
        indegree = [0] * n
        for a, b in self._precedes:
            successors[a].append(b)
            indegree[b] += 1

        def search(ready: list[int], remaining_indegree: list[int], path: list[int]):
            if len(path) == n:
                yield tuple(self._steps[i] for i in path)
                return
            # Stable choice order for determinism among equal-priority nodes.
            for idx, node in enumerate(list(ready)):
                path.append(node)
                next_ready = ready[:idx] + ready[idx + 1 :]
                next_indegree = list(remaining_indegree)
                for succ in successors[node]:
                    next_indegree[succ] -= 1
                    if next_indegree[succ] == 0:
                        next_ready.append(succ)
                yield from search(sorted(next_ready), next_indegree, path)
                path.pop()

        initial = sorted(i for i, d in enumerate(indegree) if d == 0)
        if not initial and n:
            raise ValueError("StepPlan precedes graph has a cycle")
        yield from search(initial, indegree, [])
