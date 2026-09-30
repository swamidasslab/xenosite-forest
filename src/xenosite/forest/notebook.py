"""Notebook / REPL rich displays for forest pyclasses.

Attaches ``_repr_html_`` on Rust types after the extension loads. Text
``__str__`` / ``__repr__`` live in Rust (truncated). HTML stays here so RDKit
drawings and light markup do not pull into the extension.

PathOutcome / PartialOutcome lead with **StepPlan** (elementary steps + Maybe
bags). Mol drawings always show atom labels (ForestMol tags when available,
else atom indices). SOM sites are one highlight color; Maybe formation /
span sites are another.
"""

from __future__ import annotations

import html
from collections.abc import Sequence
from typing import Any

_MAX_LIN_NOTE = 4
_MAX_HOPS = 8
_SVG_W = 240
_SVG_H = 160

# RDKit highlightAtomColors are 0–1 RGB tuples.
_COLOR_SOM = (0.95, 0.35, 0.2)  # coral — applied site
_COLOR_MAYBE = (0.25, 0.45, 0.95)  # blue — Maybe bag / span


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
    som: list[int] | None = None,
    maybe: list[int] | None = None,
    tags: list[int | None] | None = None,
    born_cutoff: int | None = None,
    label_indices: bool = True,
) -> str | None:
    Chem, rdMolDraw2D = _try_rdkit()
    if Chem is None or rdMolDraw2D is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n = mol.GetNumAtoms()
    if tags is not None:
        for i, tag in enumerate(tags):
            if i >= n:
                break
            if tag is None:
                if label_indices:
                    mol.GetAtomWithIdx(i).SetProp("atomNote", str(i))
                continue
            note = str(tag)
            if born_cutoff is not None and tag >= born_cutoff:
                note = f"{tag}*"
            mol.GetAtomWithIdx(i).SetProp("atomNote", note)
    elif label_indices:
        for i in range(n):
            mol.GetAtomWithIdx(i).SetProp("atomNote", str(i))

    som_l = [i for i in (som or []) if 0 <= i < n]
    maybe_l = [i for i in (maybe or []) if 0 <= i < n and i not in som_l]
    highlight = som_l + maybe_l
    colors = {i: _COLOR_SOM for i in som_l}
    colors.update({i: _COLOR_MAYBE for i in maybe_l})

    drawer = rdMolDraw2D.MolDraw2DSVG(_SVG_W, _SVG_H)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    if highlight:
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight,
            highlightAtomColors=colors,
        )
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
    out: list[int] = []
    for x in site:
        try:
            out.append(int(x))
        except (TypeError, ValueError):
            # PlanAtom labels like "6" or "WillAdd(...)" — keep digits-only.
            s = str(x)
            if s.isdigit():
                out.append(int(s))
    return out


def _maybe_atoms(maybe_bags: Sequence[Any]) -> list[int]:
    atoms: list[int] = []
    for bag in maybe_bags or []:
        if isinstance(bag, dict):
            for key in ("site", "span_sites"):
                val = bag.get(key)
                if val is None:
                    continue
                if val and isinstance(val[0], (list, tuple, set)):
                    for span in val:
                        atoms.extend(_sites(span))
                else:
                    atoms.extend(_sites(val))
            for span in bag.get("opens") or []:
                atoms.extend(_sites(span))
        else:
            atoms.extend(_sites(getattr(bag, "site", None)))
    # stable unique
    seen: set[int] = set()
    out: list[int] = []
    for a in atoms:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def _rule_path_label(hop: Any) -> str:
    if isinstance(hop, dict):
        path = hop.get("rule_path")
        rule = hop.get("rule") or hop.get("pattern_name") or "?"
    else:
        path = getattr(hop, "rule_path", None)
        rule = getattr(hop, "rule", None) or getattr(hop, "pattern_name", "?")
    if path:
        parts = [p for p in path if p]
        if parts:
            return "/".join(str(p) for p in parts)
    return str(rule)


def _trace_bits(mol: Any) -> tuple[str, list[int | None], int]:
    tags = list(mol.atom_tags())
    stamp_end, survivors, born, untagged = mol.trace_counts()
    bits = (
        f"trace stamp_end={stamp_end} survivors={survivors} "
        f"born={born} untagged={untagged} "
        f"(born marked with * · coral=SOM · blue=Maybe)"
    )
    return bits, tags, int(stamp_end)


def _traced_mol_block(
    mol: Any,
    *,
    som: list[int] | None = None,
    maybe: list[int] | None = None,
) -> str:
    csmi = mol.csmi
    trace, tags, stamp_end = _trace_bits(mol)
    svg = _mol_svg(
        csmi,
        som=som,
        maybe=maybe,
        tags=tags,
        born_cutoff=stamp_end,
    )
    parts = [
        f"<div><code>{_esc(csmi)}</code></div>",
        f"<div>{_esc(trace)}</div>",
    ]
    if svg:
        parts.append(svg)
    return "".join(parts)


def _panel(
    label: str,
    smiles: str,
    *,
    som: list[int] | None = None,
    maybe: list[int] | None = None,
    tags: list[int | None] | None = None,
    born_cutoff: int | None = None,
) -> str:
    svg = _mol_svg(
        smiles,
        som=som,
        maybe=maybe,
        tags=tags,
        born_cutoff=born_cutoff,
        label_indices=True,
    )
    body = svg if svg else _pre(smiles)
    return (
        "<div style='min-width:12em'>"
        f"<div style='font-size:0.85em'>{label}</div>"
        f"{body}"
        f"<div style='font-size:0.75em'><code>{_esc(smiles)}</code></div>"
        "</div>"
    )


def _arrow(caption: str) -> str:
    return (
        "<div style='align-self:center;padding:0 0.4em;text-align:center;"
        "font-size:0.85em;max-width:10em'>"
        f"{_esc(caption)}<div style='font-size:1.4em'>→</div></div>"
    )


def _legend() -> str:
    return (
        "<div style='font-size:0.8em;margin:0.25em 0'>"
        "<span style='color:#c24'>■</span> SOM / step site &nbsp; "
        "<span style='color:#247'>■</span> Maybe bag &nbsp; "
        "atom notes = tags (or indices)"
        "</div>"
    )


def _pathway_row(
    hops: Sequence[Any],
    *,
    final_smiles: str | None = None,
    final_tags: list[int | None] | None = None,
    final_born: int | None = None,
    maybe_atoms: list[int] | None = None,
) -> str:
    """Horizontal reactant(SOM) → end trail; atom labels on every mol."""

    maybe_atoms = maybe_atoms or []
    if not hops:
        if final_smiles:
            return _panel(
                "<b>product</b>",
                final_smiles,
                maybe=maybe_atoms,
                tags=final_tags,
                born_cutoff=final_born,
            )
        return "<div><i>(no hops)</i></div>"

    parts: list[str] = [
        "<div style='display:flex;gap:0.25em;flex-wrap:wrap;align-items:flex-start'>"
    ]
    shown = list(hops[:_MAX_HOPS])
    for i, hop in enumerate(shown):
        if isinstance(hop, dict):
            reactant = hop.get("reactant") or ""
            product = hop.get("product") or ""
            site = _sites(hop.get("site"))
            sides = hop.get("sides") or []
        else:
            reactant = getattr(hop, "reactant", "") or ""
            product = getattr(hop, "product", "") or ""
            site = _sites(getattr(hop, "site", None))
            sides = getattr(hop, "sides", None) or []
        # Cleavage hops: Maybe sites belong on the reactant this rule hit.
        hop_maybe = maybe_atoms if sides else ([] if i < len(shown) - 1 else maybe_atoms)
        if i == 0 and maybe_atoms and not hop_maybe:
            hop_maybe = maybe_atoms
        path_label = _rule_path_label(hop)
        site_s = ",".join(str(s) for s in site[:6])
        label_r = "<b>start</b>" if i == 0 else f"<b>hop {i}</b>"
        parts.append(
            _panel(
                f"{label_r} (SOM)",
                str(reactant),
                som=site,
                maybe=hop_maybe if sides or i == 0 else None,
            )
        )
        parts.append(_arrow(f"{path_label} @ [{site_s}]"))
        if i == len(shown) - 1:
            end = final_smiles or product
            parts.append(
                _panel(
                    "<b>end</b>",
                    str(end),
                    maybe=maybe_atoms,
                    tags=final_tags,
                    born_cutoff=final_born,
                )
            )
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
    maybe = d.get("maybe") or list(self.maybe())
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
    if maybe:
        lines.append(f"  maybe: {len(maybe)} bag(s)")
        for i, bag in enumerate(maybe[:8]):
            site = bag.get("site") if isinstance(bag, dict) else getattr(bag, "site", [])
            side = bag.get("side") if isinstance(bag, dict) else getattr(bag, "side", "")
            site_s = ",".join(str(x) for x in _sites(site)[:8])
            lines.append(f"    [{i}] site={{{site_s}}} side={side}")
        if len(maybe) > 8:
            lines.append(f"    … (+{len(maybe) - 8} bags)")
    return f"<div class='xenosite-stepplan'>{note}{_pre(chr(10).join(lines))}</div>"


def _path_outcome_html(self: Any) -> str:
    """StepPlan (+ Maybe) first; tagged pathway with SOM vs Maybe colors."""

    plan = self.plan
    hops = list(self.hops())
    maybe_bags = list(self.maybe())
    maybe_atoms = _maybe_atoms(maybe_bags)
    # Plan-step sites as SOM for the end-mol overlay when no hop site applies.
    plan_som: list[int] = []
    for step in plan.to_dict().get("steps") or []:
        plan_som.extend(_sites(step.get("site")))
    tags = list(self.mol.atom_tags())
    stamp_end = int(self.mol.trace_counts()[0])
    head = (
        f"<div><b>PathOutcome</b> plan {len(plan)} step(s)"
        f" · ~{plan.n_linearizations()} lin"
        f" · maybe {len(maybe_bags)} · {len(hops)} walk hop(s)</div>"
    )
    return (
        "<div class='xenosite-path'>"
        + head
        + _legend()
        + _step_plan_html(plan)
        + _pathway_row(
            hops,
            final_smiles=self.smiles,
            final_tags=tags,
            final_born=stamp_end,
            maybe_atoms=maybe_atoms,
        )
        + "<details><summary>tagged end mol</summary>"
        + _traced_mol_block(self.mol, som=plan_som[:12], maybe=maybe_atoms)
        + _pre(str(self))
        + "</details>"
        + "</div>"
    )


def _partial_outcome_html(self: Any) -> str:
    plan = self.plan
    hops = list(self.hops())
    maybe_bags = list(self.maybe())
    maybe_atoms = _maybe_atoms(maybe_bags)
    tags = list(self.mol.atom_tags())
    stamp_end = int(self.mol.trace_counts()[0])
    head = (
        f"<div><b>PartialOutcome</b> residual_cost={self.residual_cost}"
        f" · plan {len(plan)} · maybe {len(maybe_bags)}</div>"
    )
    return (
        "<div class='xenosite-partial'>"
        + head
        + _legend()
        + _step_plan_html(plan)
        + _pathway_row(
            hops,
            final_smiles=self.smiles,
            final_tags=tags,
            final_born=stamp_end,
            maybe_atoms=maybe_atoms,
        )
        + "<details><summary>tagged end mol</summary>"
        + _traced_mol_block(self.mol, maybe=maybe_atoms)
        + _pre(str(self))
        + "</details>"
        + "</div>"
    )


def _emission_html(self: Any) -> str:
    """SOM on the reactant; products with atom labels, no SOM marks."""

    site = int(self.site)
    reactant = self.reactant
    products = list(self.products())
    path = list(self.rule_path())
    # Outer→leaf for display caption.
    path_s = "/".join(p for p in reversed(path) if p) or self.pattern_name
    r_tags = list(reactant.atom_tags()) if hasattr(reactant, "atom_tags") else None
    r_born = int(reactant.trace_counts()[0]) if hasattr(reactant, "trace_counts") else None
    blocks = [
        f"<div><b>Emission</b> {_esc(self.pattern_name)} site={site}"
        f" · <code>{_esc(path_s)}</code></div>",
        _legend(),
        "<div style='display:flex;gap:0.5em;flex-wrap:wrap;align-items:flex-start'>",
        _panel(
            "<b>reactant</b> (SOM)",
            reactant.csmi,
            som=[site],
            tags=r_tags,
            born_cutoff=r_born,
        ),
        _arrow(f"{path_s} @ [{site}]"),
    ]
    for i, prod in enumerate(products[:4]):
        smi = prod.csmi if hasattr(prod, "csmi") else self.product_csmis()[i]
        p_tags = list(prod.atom_tags()) if hasattr(prod, "atom_tags") else None
        p_born = int(prod.trace_counts()[0]) if hasattr(prod, "trace_counts") else None
        blocks.append(
            _panel(f"product {i}", smi, tags=p_tags, born_cutoff=p_born)
        )
    if len(products) > 4:
        blocks.append(f"<div>… (+{len(products) - 4})</div>")
    blocks.append("</div>")
    blocks.append(_pre(str(self)))
    return "<div class='xenosite-emission'>" + "".join(blocks) + "</div>"


def _random_path_html(self: Any) -> str:
    path = list(self.path)
    steps = list(self.steps()) if callable(getattr(self, "steps", None)) else []
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
        + _legend()
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
    """SOM on parent (reactant); kept child with atom labels."""

    site = int(self.site)
    parent = self.parent_mol
    kept = self.kept_mol
    path = self.rule  # leaf; GraphEdge has no full rule_path today
    head = (
        f"<div><b>GraphEdge</b> {_esc(path)} @ site={site}"
        f" · parent={self.parent_index}→{self.child_index}</div>"
    )
    p_tags = list(parent.atom_tags())
    p_born = int(parent.trace_counts()[0])
    k_tags = list(kept.atom_tags())
    k_born = int(kept.trace_counts()[0])
    return (
        "<div class='xenosite-graph-edge'>"
        + head
        + _legend()
        + "<div style='display:flex;gap:0.5em;flex-wrap:wrap;align-items:flex-start'>"
        + _panel(
            "<b>parent</b> (SOM)",
            parent.csmi,
            som=[site],
            tags=p_tags,
            born_cutoff=p_born,
        )
        + _arrow(f"{path} @ [{site}]")
        + _panel("<b>kept</b>", kept.csmi, tags=k_tags, born_cutoff=k_born)
        + "</div>"
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
