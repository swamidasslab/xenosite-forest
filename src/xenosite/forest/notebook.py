"""Notebook / REPL rich displays for forest pyclasses.

Attaches ``_repr_html_`` on Rust types after the extension loads. Text
``__str__`` / ``__repr__`` live in Rust (truncated). HTML stays here so RDKit
drawings and light markup do not pull into the extension.

Prefer **ForestMol tag tracing** (survivors vs born, stamp_end, atom notes)
over search hop lists. PathOutcome shows the tagged hit mol + StepPlan; it
does not dump per-hop intermediates unless tracing is unavailable.
"""

from __future__ import annotations

import html
from typing import Any

_MAX_LIN_NOTE = 4
_SVG_W = 280
_SVG_H = 180


def _esc(s: Any) -> str:
    return html.escape(str(s), quote=True)


def _try_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        return Chem, rdMolDraw2D
    except ImportError:
        return None, None


def _mol_svg(
    smiles: str,
    *,
    highlight: list[int] | None = None,
    tags: list[int | None] | None = None,
    born_cutoff: int | None = None,
) -> str | None:
    Chem, rdMolDraw2D = _try_rdkit()
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if tags is not None:
        for i, tag in enumerate(tags):
            if i >= mol.GetNumAtoms():
                break
            if tag is None:
                continue
            # Survivor tags: plain id; born atoms: id* (minted after root stamp).
            note = str(tag)
            if born_cutoff is not None and tag >= born_cutoff:
                note = f"{tag}*"
            mol.GetAtomWithIdx(i).SetProp("atomNote", note)
    drawer = rdMolDraw2D.MolDraw2DSVG(_SVG_W, _SVG_H)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    highlight = [i for i in (highlight or []) if 0 <= i < mol.GetNumAtoms()]
    # Born atoms get a soft highlight when no explicit SOM list is given.
    if not highlight and tags is not None and born_cutoff is not None:
        highlight = [
            i
            for i, tag in enumerate(tags)
            if tag is not None and tag >= born_cutoff and i < mol.GetNumAtoms()
        ]
    if highlight:
        drawer.DrawMolecule(mol, highlightAtoms=highlight)
    else:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _pre(text: str) -> str:
    return f"<pre style='margin:0.25em 0;white-space:pre-wrap'>{_esc(text)}</pre>"


def _trace_bits(mol: Any) -> tuple[str, list[int | None], int]:
    tags = list(mol.atom_tags())
    stamp_end, survivors, born, untagged = mol.trace_counts()
    bits = (
        f"trace stamp_end={stamp_end} survivors={survivors} "
        f"born={born} untagged={untagged} "
        f"(born marked with * on drawing)"
    )
    return bits, tags, int(stamp_end)


def _forest_mol_html(self: Any) -> str:
    csmi = self.csmi
    trace, tags, stamp_end = _trace_bits(self)
    svg = _mol_svg(csmi, tags=tags, born_cutoff=stamp_end)
    body = svg if svg else _pre(str(self))
    return (
        "<div class='xenosite-forest-mol'>"
        f"<div><b>ForestMol</b> <code>{_esc(csmi)}</code></div>"
        f"<div>{_esc(trace)}</div>"
        f"{body}"
        f"{_pre(str(self))}"
        "</div>"
    )


def _ruleset_html(self: Any) -> str:
    return f"<div class='xenosite-ruleset'>{_pre(str(self))}</div>"


def _network_html(self: Any) -> str:
    root = self.root_csmi
    targets = list(self.target_csmis())
    reached = any(self.reaches(t) for t in targets) if targets else None
    bits = [
        f"<b>MetabolicNetwork</b> nodes={self.n_nodes()} edges={self.n_edges()}",
        f"root: <code>{_esc(root)}</code>" if root else "root: (none)",
    ]
    if targets:
        shown = ", ".join(f"<code>{_esc(t)}</code>" for t in targets[:4])
        if len(targets) > 4:
            shown += f" … (+{len(targets) - 4})"
        bits.append(f"targets: {shown}")
        bits.append(f"target_reached: <b>{'yes' if reached else 'no'}</b>")
    else:
        bits.append("targets: (none marked)")
    root_svg = _mol_svg(root) if root else None
    tgt_svg = _mol_svg(targets[0]) if targets else None
    imgs = ""
    if root_svg or tgt_svg:
        imgs = "<div style='display:flex;gap:1em;flex-wrap:wrap'>"
        if root_svg:
            imgs += f"<div><div>start</div>{root_svg}</div>"
        if tgt_svg:
            imgs += f"<div><div>target</div>{tgt_svg}</div>"
        imgs += "</div>"
    return (
        "<div class='xenosite-network'>"
        + "<br/>".join(bits)
        + imgs
        + _pre(str(self))
        + "</div>"
    )


def _step_plan_html(self: Any) -> str:
    n_lin = self.n_linearizations()
    d = self.to_dict()
    steps = d.get("steps") or []
    lines = [f"StepPlan · {len(self)} steps · ~{n_lin} linearization(s)"]
    for i, step in enumerate(steps[:16]):
        rule = step.get("rule", "?")
        site = step.get("site", [])
        site_s = ",".join(str(x) for x in site[:6])
        if len(site) > 6:
            site_s += ",…"
        lines.append(f"  {i}: {rule} @ [{site_s}]")
    if len(steps) > 16:
        lines.append(f"  … ({len(steps) - 16} more steps)")
    precedes = d.get("precedes") or []
    if precedes:
        edges = ", ".join(f"{a}→{b}" for a, b in precedes[:12])
        if len(precedes) > 12:
            edges += " …"
        lines.append(f"  precedes: {edges}")
    note = ""
    if n_lin > _MAX_LIN_NOTE:
        note = (
            f"<div><i>~{n_lin} linearizations (not expanded; "
            "use StepPlan text / to_dict for edges).</i></div>"
        )
    return f"<div class='xenosite-stepplan'>{note}{_pre(chr(10).join(lines))}</div>"


def _path_outcome_html(self: Any) -> str:
    """Hit mol via ForestMol tracing + elementary StepPlan — not hop dumps."""

    mol = self.mol
    plan = self.plan
    trace, tags, stamp_end = _trace_bits(mol)
    svg = _mol_svg(self.smiles, tags=tags, born_cutoff=stamp_end)
    head = (
        f"<div><b>PathOutcome</b> <code>{_esc(self.smiles)}</code>"
        f" · plan {len(plan)} elementary step(s)"
        f" · ~{plan.n_linearizations()} lin</div>"
        f"<div>{_esc(trace)}</div>"
    )
    plan_block = _step_plan_html(plan)
    return (
        "<div class='xenosite-path'>"
        + head
        + (svg or "")
        + plan_block
        + _pre(str(self))
        + "</div>"
    )


def _emission_html(self: Any) -> str:
    products = list(self.products())
    site = int(self.site)
    blocks = [
        f"<div><b>Emission</b> {_esc(self.pattern_name)} site={site}</div>",
        _pre(str(self)),
    ]
    # Prefer tagged ForestMol products when available.
    for i, prod in enumerate(products[:4]):
        if hasattr(prod, "atom_tags") and hasattr(prod, "trace_counts"):
            trace, tags, stamp_end = _trace_bits(prod)
            # Highlight SOM site when it is still in range; else born atoms.
            highlight = [site] if site < len(tags) else None
            svg = _mol_svg(
                prod.csmi,
                highlight=highlight,
                tags=tags,
                born_cutoff=stamp_end,
            )
            blocks.append(
                f"<div>product {i}: <code>{_esc(prod.csmi)}</code> · {_esc(trace)}</div>"
            )
            if svg:
                blocks.append(svg)
        else:
            smi = self.product_csmis()[i]
            svg = _mol_svg(smi, highlight=[site])
            blocks.append(f"<div>product {i}: <code>{_esc(smi)}</code></div>")
            if svg:
                blocks.append(svg)
    if len(products) > 4:
        blocks.append(f"<div>… ({len(products) - 4} more products)</div>")
    return "<div class='xenosite-emission'>" + "".join(blocks) + "</div>"


def _random_path_html(self: Any) -> str:
    svg = _mol_svg(self.smiles)
    return (
        "<div class='xenosite-random-path'>"
        f"<div><b>RandomPathOutcome</b> <code>{_esc(self.smiles)}</code></div>"
        + (svg or "")
        + _pre(str(self))
        + "</div>"
    )


def install(rust: Any) -> None:
    """Attach ``_repr_html_`` callables onto loaded ``_rust`` pyclasses."""

    mapping = {
        "ForestMol": _forest_mol_html,
        "RuleSet": _ruleset_html,
        "BoundPattern": lambda self: _pre(str(self)),
        "MetabolicNetwork": _network_html,
        "StepPlan": _step_plan_html,
        "PathOutcome": _path_outcome_html,
        "PartialOutcome": lambda self: _pre(str(self)),
        "Emission": _emission_html,
        "RandomPathOutcome": _random_path_html,
        "PathCounters": lambda self: _pre(str(self.to_dict())),
        "GraphNode": lambda self: _pre(
            f"GraphNode[{self.index}] {self.csmi} sealed={self.sealed} expanded={self.expanded}"
        ),
        "GraphEdge": lambda self: _pre(
            f"GraphEdge {self.rule} site={self.site} parent={self.parent_index}→{self.child_index}"
        ),
    }
    for name, fn in mapping.items():
        cls = getattr(rust, name, None)
        if cls is not None:
            cls._repr_html_ = fn
