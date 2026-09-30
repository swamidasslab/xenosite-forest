"""Notebook / REPL rich displays for forest pyclasses.

Attaches ``_repr_html_`` on Rust types after the extension loads. Text
``__str__`` / ``__repr__`` live in Rust (truncated). HTML stays here so RDKit
drawings and light markup do not pull into the extension.

Pathway displays draw each hop's **reactant** with the SOM site highlighted
(the mol the rule was applied to), then the product. Final tagged ForestMol
trace is secondary. Do not mark sites on the terminal product.
"""

from __future__ import annotations

import html
from typing import Any, Sequence

_MAX_LIN_NOTE = 4
_MAX_HOPS = 8
_SVG_W = 220
_SVG_H = 150


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


def _sites(site: Any) -> list[int]:
    if site is None:
        return []
    if isinstance(site, int):
        return [site]
    return [int(x) for x in site]


def _trace_bits(mol: Any) -> tuple[str, list[int | None], int]:
    tags = list(mol.atom_tags())
    stamp_end, survivors, born, untagged = mol.trace_counts()
    bits = (
        f"trace stamp_end={stamp_end} survivors={survivors} "
        f"born={born} untagged={untagged} "
        f"(born marked with * on drawing)"
    )
    return bits, tags, int(stamp_end)


def _traced_mol_block(mol: Any, *, highlight: list[int] | None = None) -> str:
    csmi = mol.csmi
    trace, tags, stamp_end = _trace_bits(mol)
    svg = _mol_svg(csmi, highlight=highlight, tags=tags, born_cutoff=stamp_end)
    parts = [
        f"<div><code>{_esc(csmi)}</code></div>",
        f"<div>{_esc(trace)}</div>",
    ]
    if svg:
        parts.append(svg)
    return "".join(parts)


def _panel(label: str, smiles: str, *, highlight: list[int] | None = None) -> str:
    svg = _mol_svg(smiles, highlight=highlight)
    body = svg if svg else _pre(smiles)
    return (
        "<div style='min-width:11em'>"
        f"<div style='font-size:0.85em'>{label}</div>"
        f"{body}"
        f"<div style='font-size:0.75em'><code>{_esc(smiles)}</code></div>"
        "</div>"
    )


def _arrow(caption: str) -> str:
    return (
        "<div style='align-self:center;padding:0 0.4em;text-align:center;"
        "font-size:0.85em;max-width:9em'>"
        f"{_esc(caption)}<div style='font-size:1.4em'>→</div></div>"
    )


def _pathway_row(hops: Sequence[Any], *, final_smiles: str | None = None) -> str:
    """Horizontal reactant(SOM) → product trail from hop dicts / objects."""

    if not hops:
        if final_smiles:
            return _panel("<b>product</b>", final_smiles)
        return "<div><i>(no hops)</i></div>"

    parts: list[str] = [
        "<div style='display:flex;gap:0.25em;flex-wrap:wrap;align-items:flex-start'>"
    ]
    shown = list(hops[:_MAX_HOPS])
    for i, hop in enumerate(shown):
        if isinstance(hop, dict):
            reactant = hop.get("reactant") or ""
            product = hop.get("product") or ""
            rule = hop.get("rule") or hop.get("pattern_name") or "?"
            site = _sites(hop.get("site"))
        else:
            reactant = getattr(hop, "reactant", "") or ""
            product = getattr(hop, "product", "") or ""
            rule = getattr(hop, "rule", None) or getattr(hop, "pattern_name", "?")
            site = _sites(getattr(hop, "site", None))
        site_s = ",".join(str(s) for s in site[:6])
        label_r = f"<b>hop {i}</b> reactant"
        if i == 0:
            label_r = "<b>start</b> (SOM)"
        parts.append(_panel(label_r, str(reactant), highlight=site))
        parts.append(_arrow(f"{rule} @ [{site_s}]"))
        if i == len(shown) - 1:
            end = final_smiles or product
            parts.append(_panel("<b>end</b>", str(end)))
        else:
            # Intermediate product is next hop's reactant; omit duplicate draw.
            pass
    if len(hops) > _MAX_HOPS:
        parts.append(
            f"<div style='align-self:center'>… (+{len(hops) - _MAX_HOPS} hops)</div>"
        )
    parts.append("</div>")
    return "".join(parts)


def _forest_mol_html(self: Any) -> str:
    body = _traced_mol_block(self)
    return (
        "<div class='xenosite-forest-mol'>"
        f"<div><b>ForestMol</b></div>"
        f"{body}"
        f"{_pre(str(self))}"
        "</div>"
    )


def _ruleset_html(self: Any) -> str:
    return f"<div class='xenosite-ruleset'>{_pre(str(self))}</div>"


def _bound_pattern_html(self: Any) -> str:
    return (
        "<div class='xenosite-bound-pattern'>"
        f"<div><b>BoundPattern</b> <code>{_esc(self.curie)}</code>"
        f" · {len(self)} pattern(s)</div>"
        f"{_pre(str(self))}"
        "</div>"
    )


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
    if 0 < n_lin <= _MAX_LIN_NOTE:
        lines.append("  linearizations:")
        for li, lin in enumerate(self.linearizations()):
            order = []
            for step in lin:
                site = step.get("site") or []
                site_s = ",".join(str(x) for x in site[:4])
                if len(site) > 4:
                    site_s += ",…"
                order.append(f"{step.get('rule', '?')}@[{site_s}]")
            lines.append(f"    [{li}] {' → '.join(order)}")
    elif n_lin > _MAX_LIN_NOTE:
        note = (
            f"<div><i>~{n_lin} linearizations (not expanded; "
            "use StepPlan.linearizations() / to_dict).</i></div>"
        )
    return f"<div class='xenosite-stepplan'>{note}{_pre(chr(10).join(lines))}</div>"


def _path_outcome_html(self: Any) -> str:
    """Multistep walk: SOM on each reactant; end mol; StepPlan summary."""

    hops = list(self.hops())
    plan = self.plan
    head = (
        f"<div><b>PathOutcome</b> {len(hops)} hop(s)"
        f" · plan {len(plan)} elementary step(s)"
        f" · ~{plan.n_linearizations()} lin</div>"
    )
    return (
        "<div class='xenosite-path'>"
        + head
        + _pathway_row(hops, final_smiles=self.smiles)
        + "<details><summary>tagged end mol + StepPlan</summary>"
        + _traced_mol_block(self.mol)
        + _step_plan_html(plan)
        + _pre(str(self))
        + "</details>"
        + "</div>"
    )


def _partial_outcome_html(self: Any) -> str:
    hops = list(self.hops())
    plan = self.plan
    head = (
        f"<div><b>PartialOutcome</b> residual_cost={self.residual_cost}"
        f" · {len(hops)} hop(s)"
        f" · plan {len(plan)} step(s)</div>"
    )
    return (
        "<div class='xenosite-partial'>"
        + head
        + _pathway_row(hops, final_smiles=self.smiles)
        + "<details><summary>tagged end mol + StepPlan</summary>"
        + _traced_mol_block(self.mol)
        + _step_plan_html(plan)
        + _pre(str(self))
        + "</details>"
        + "</div>"
    )


def _emission_html(self: Any) -> str:
    """SOM on the reactant; products without site marks."""

    site = int(self.site)
    reactant = self.reactant
    products = list(self.products())
    blocks = [
        f"<div><b>Emission</b> {_esc(self.pattern_name)} site={site}</div>",
        "<div style='display:flex;gap:0.5em;flex-wrap:wrap;align-items:flex-start'>",
        _panel("<b>reactant</b> (SOM)", reactant.csmi, highlight=[site]),
        _arrow(f"{self.pattern_name} @ [{site}]"),
    ]
    for i, prod in enumerate(products[:4]):
        smi = prod.csmi if hasattr(prod, "csmi") else self.product_csmis()[i]
        blocks.append(_panel(f"product {i}", smi))
    if len(products) > 4:
        blocks.append(f"<div>… (+{len(products) - 4})</div>")
    blocks.append("</div>")
    blocks.append(_pre(str(self)))
    return "<div class='xenosite-emission'>" + "".join(blocks) + "</div>"


def _random_path_html(self: Any) -> str:
    path = list(self.path)
    steps = list(self.steps()) if callable(getattr(self, "steps", None)) else []
    # Prefer structured steps when present; else CSMI trail only.
    hops: list[dict[str, Any]] = []
    if steps and len(path) >= 2:
        for i, step in enumerate(steps):
            if isinstance(step, dict):
                rule = step.get("rule") or step.get("pattern") or "?"
                site = step.get("site")
                chosen = int(step.get("chosen", 0))
                products = step.get("products") or []
                product = products[chosen] if products else path[min(i + 1, len(path) - 1)]
            else:
                rule = getattr(step, "rule", "?")
                site = getattr(step, "site", [])
                chosen = int(getattr(step, "chosen", 0))
                products = list(getattr(step, "products", []))
                product = products[chosen] if products else path[min(i + 1, len(path) - 1)]
            reactant = path[i] if i < len(path) else ""
            hops.append(
                {
                    "reactant": reactant,
                    "product": product,
                    "rule": rule,
                    "site": site,
                }
            )
    head = f"<div><b>RandomPathOutcome</b> {max(0, len(path) - 1)} hop(s)</div>"
    body = (
        _pathway_row(hops, final_smiles=self.smiles)
        if hops
        else _panel("<b>end</b>", self.smiles)
    )
    return (
        "<div class='xenosite-random-path'>"
        + head
        + body
        + _pre(str(self))
        + "</div>"
    )


def _path_counters_html(self: Any) -> str:
    return f"<div class='xenosite-counters'>{_pre(str(self))}</div>"


def _graph_node_html(self: Any) -> str:
    mol = self.mol
    head = (
        f"<div><b>GraphNode</b>[{self.index}]"
        f" sealed={self.sealed} expanded={self.expanded}"
        f" inbound={self.n_inbound()}</div>"
    )
    return (
        "<div class='xenosite-graph-node'>"
        + head
        + _traced_mol_block(mol)
        + _pre(str(self))
        + "</div>"
    )


def _graph_edge_html(self: Any) -> str:
    """SOM on parent (reactant); kept child without site marks."""

    site = int(self.site)
    parent = self.parent_mol
    kept = self.kept_mol
    head = (
        f"<div><b>GraphEdge</b> {_esc(self.rule)} @ site={site}"
        f" · parent={self.parent_index}→{self.child_index}</div>"
    )
    hop = {
        "reactant": parent.csmi,
        "product": kept.csmi,
        "rule": self.rule,
        "site": site,
    }
    return (
        "<div class='xenosite-graph-edge'>"
        + head
        + _pathway_row([hop], final_smiles=kept.csmi)
        + _pre(str(self))
        + "</div>"
    )


def install(rust: Any) -> None:
    """Attach ``_repr_html_`` callables onto loaded ``_rust`` pyclasses."""

    mapping = {
        "ForestMol": _forest_mol_html,
        "RuleSet": _ruleset_html,
        "BoundPattern": _bound_pattern_html,
        "MetabolicNetwork": _network_html,
        "StepPlan": _step_plan_html,
        "PathOutcome": _path_outcome_html,
        "PartialOutcome": _partial_outcome_html,
        "Emission": _emission_html,
        "RandomPathOutcome": _random_path_html,
        "PathCounters": _path_counters_html,
        "GraphNode": _graph_node_html,
        "GraphEdge": _graph_edge_html,
    }
    for name, fn in mapping.items():
        cls = getattr(rust, name, None)
        if cls is not None:
            cls._repr_html_ = fn
