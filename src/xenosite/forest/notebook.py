"""Notebook / REPL rich displays for forest pyclasses.

Attaches ``_repr_html_`` on Rust types after the extension loads. Text
``__str__`` / ``__repr__`` live in Rust (truncated). HTML stays here so RDKit
drawings and light markup do not pull into the extension.

PathOutcome / PartialOutcome lead with **StepPlan** (elementary steps + Maybe
bags). Mol drawings omit atom numbers (RDKit cannot show ForestMol trace
labels correctly). SOM / Maybe are color highlights only.
"""

from __future__ import annotations

import html
from collections.abc import Sequence
from typing import Any

_MAX_LIN_NOTE = 4
_MAX_HOPS = 8
_SVG_W = 240
_SVG_H = 160

# RDKit highlightAtomColors are 0–1 RGB tuples (optional alpha).
_COLOR_SOM = (0.95, 0.35, 0.2, 1.0)  # coral — applied site
_COLOR_MAYBE = (0.25, 0.45, 0.95, 0.30)  # translucent blue — full Maybe bag


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
    site_kind: str | None = None,
    tags: list[int | None] | None = None,
    born_cutoff: int | None = None,
) -> str | None:
    """Draw ``smiles`` with SOM / Maybe highlights. No atom numbers.

    ``tags`` / ``born_cutoff`` are accepted for call-site compatibility but
    unused — RDKit atom notes cannot carry ForestMol trace labels reliably.
    """
    del tags, born_cutoff  # numbering off
    Chem, rdMolDraw2D = _try_rdkit()
    if Chem is None or rdMolDraw2D is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n = mol.GetNumAtoms()

    som_l = [i for i in (som or []) if 0 <= i < n]
    maybe_l = [i for i in (maybe or []) if 0 <= i < n and i not in som_l]

    highlight_atoms = list(som_l)
    atom_colors = {i: _COLOR_SOM for i in som_l}
    # Maybe: all bag atoms, translucent so they do not obscure the drawing.
    for i in maybe_l:
        highlight_atoms.append(i)
        atom_colors[i] = _COLOR_MAYBE

    highlight_bonds: list[int] = []
    bond_colors: dict[int, Any] = {}
    kind = (site_kind or "atom").lower()
    if kind in ("bond", "directed_bond") and len(som_l) >= 2:
        for a, b in zip(som_l[::2], som_l[1::2], strict=False):
            bond = mol.GetBondBetweenAtoms(int(a), int(b))
            if bond is None and len(som_l) == 2:
                bond = mol.GetBondBetweenAtoms(int(som_l[0]), int(som_l[1]))
            if bond is not None:
                bi = bond.GetIdx()
                highlight_bonds.append(bi)
                bond_colors[bi] = _COLOR_SOM
        # Bond sites: prefer bond highlight; keep endpoint atoms lightly if needed.
        if highlight_bonds:
            highlight_atoms = [i for i in highlight_atoms if i not in som_l] + list(som_l)
    elif kind == "atom_pair":
        # Pair ends: atom highlights only (not the intervening path as a bond).
        pass

    drawer = rdMolDraw2D.MolDraw2DSVG(_SVG_W, _SVG_H)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    opts.addStereoAnnotation = False
    if highlight_atoms or highlight_bonds:
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atoms or None,
            highlightAtomColors=atom_colors or None,
            highlightBonds=highlight_bonds or None,
            highlightBondColors=bond_colors or None,
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
    """All atoms in Maybe bags (formation site + opens / span_sites)."""
    atoms: list[int] = []
    for bag in maybe_bags or []:
        if isinstance(bag, dict):
            # Prefer span_sites (site + opens); else site + opens separately.
            spans = bag.get("span_sites")
            if spans:
                for span in spans:
                    atoms.extend(_sites(span))
            else:
                atoms.extend(_sites(bag.get("site")))
                for span in bag.get("opens") or []:
                    atoms.extend(_sites(span))
        else:
            atoms.extend(_sites(getattr(bag, "site", None)))
    seen: set[int] = set()
    out: list[int] = []
    for a in atoms:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def _path_parts(path: Any) -> list[str]:
    """Outer→leaf path segments, omitting empty and the ``Default`` catalog root."""
    if not path:
        return []
    parts: list[str] = []
    for p in path:
        if not p:
            continue
        s = str(p)
        if s == "Default":
            continue
        parts.append(s)
    return parts


def _rule_path_label(hop: Any) -> str:
    if isinstance(hop, dict):
        path = hop.get("rule_path")
        rule = hop.get("rule") or hop.get("pattern_name") or "?"
    else:
        path = getattr(hop, "rule_path", None)
        rule = getattr(hop, "rule", None) or getattr(hop, "pattern_name", "?")
    parts = _path_parts(path)
    if parts:
        return "/".join(parts)
    return str(rule)


def _trace_bits(mol: Any) -> tuple[str, list[int | None], int]:
    tags = list(mol.atom_tags())
    stamp_end, survivors, born, untagged = mol.trace_counts()
    bits = (
        f"trace stamp_end={stamp_end} survivors={survivors} "
        f"born={born} untagged={untagged} "
        f"(coral=SOM · translucent blue=Maybe)"
    )
    return bits, tags, int(stamp_end)


def _traced_mol_block(
    mol: Any,
    *,
    som: list[int] | None = None,
    maybe: list[int] | None = None,
    site_kind: str | None = None,
) -> str:
    csmi = mol.csmi
    trace, _tags, _stamp_end = _trace_bits(mol)
    svg = _mol_svg(
        csmi,
        som=som,
        maybe=maybe,
        site_kind=site_kind,
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
    site_kind: str | None = None,
    tags: list[int | None] | None = None,
    born_cutoff: int | None = None,
) -> str:
    del tags, born_cutoff
    svg = _mol_svg(
        smiles,
        som=som,
        maybe=maybe,
        site_kind=site_kind,
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
        "<span style='color:#247;opacity:0.45'>■</span> Maybe bag (translucent) "
        "· no atom numbers"
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
    """Horizontal reactant(SOM) → end trail. No atom numbers on drawings."""
    del final_tags, final_born

    maybe_atoms = maybe_atoms or []
    if not hops:
        if final_smiles:
            return _panel(
                "<b>product</b>",
                final_smiles,
                maybe=maybe_atoms,
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
            site_kind = hop.get("site_kind")
        else:
            reactant = getattr(hop, "reactant", "") or ""
            product = getattr(hop, "product", "") or ""
            site = _sites(getattr(hop, "site", None))
            sides = getattr(hop, "sides", None) or []
            site_kind = getattr(hop, "site_kind", None)
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
                site_kind=site_kind,
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
    """SOM on the reactant by site_kind; products without SOM marks."""

    site_atoms = list(getattr(self, "site_atoms", None) or [int(self.site)])
    site_kind = getattr(self, "site_kind", None)
    reactant = self.reactant
    products = list(self.products())
    path_s = "/".join(_path_parts(reversed(list(self.rule_path())))) or self.pattern_name
    site_s = ",".join(str(s) for s in site_atoms[:6])
    blocks = [
        f"<div><b>Emission</b> {_esc(self.pattern_name)}"
        f" · { _esc(site_kind or '?') } [{_esc(site_s)}]"
        f" · <code>{_esc(path_s)}</code></div>",
        _legend(),
        "<div style='display:flex;gap:0.5em;flex-wrap:wrap;align-items:flex-start'>",
        _panel(
            "<b>reactant</b> (SOM)",
            reactant.csmi,
            som=site_atoms,
            site_kind=site_kind,
        ),
        _arrow(f"{path_s} @ [{site_s}]"),
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
    return (
        "<div class='xenosite-graph-edge'>"
        + head
        + _legend()
        + "<div style='display:flex;gap:0.5em;flex-wrap:wrap;align-items:flex-start'>"
        + _panel(
            "<b>parent</b> (SOM)",
            parent.csmi,
            som=[site],
        )
        + _arrow(f"{path} @ [{site}]")
        + _panel("<b>kept</b>", kept.csmi)
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
