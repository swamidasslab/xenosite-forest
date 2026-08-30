"""Pretty HTML/SVG drawings of metabolite networks for notebooks."""

from __future__ import annotations

import html
import itertools
import shutil
import subprocess
import tempfile
from collections import defaultdict, deque
from pathlib import Path
from typing import Iterable

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolfiles import MolFromSmiles

# Metabolic Rainbow (Wong / Okabe–Ito), plus extra classes from Full.
CLASS_STYLE = {
    "SO": ("#D55E00", "Stable oxygenation"),
    "UO": ("#E69F00", "Unstable oxygenation"),
    "DH": ("#009E73", "Dehydrogenation"),
    "HD": ("#56B4E9", "Hydrolysis"),
    "RD": ("#CC79A7", "Reduction"),
    "QF": ("#0072B2", "Quinone formation"),
    "CJ": ("#882255", "Conjugation"),
    "TT": ("#999999", "Tautomerization"),
    "TSO": ("#D55E00", "Thiophene S-oxidation"),
    "BA": ("#0072B2", "Bioactivation"),
}
_DEFAULT_STYLE = ("#666666", "Other")
_IDS = itertools.count()


class Drawing:
    """Notebook-friendly figure. Jupyter calls ``_repr_html_`` automatically.

    Use :meth:`save` to write SVG (always) or PDF/PS (via ``rsvg-convert`` or
    ImageMagick when installed).
    """

    def __init__(self, html: str, svg: str | None = None):
        self.html = html
        self.svg = svg

    def _repr_html_(self) -> str:
        return self.html

    def _repr_svg_(self) -> str | None:
        return self.svg

    def __str__(self) -> str:
        return self.html

    def __repr__(self) -> str:
        return f"<Drawing {len(self.html)} chars>"

    def save(self, path, *, format: str | None = None) -> Path:
        """Save the figure to ``path``.

        Formats:
          * ``.svg`` — written directly (no extra tools)
          * ``.pdf`` / ``.ps`` / ``.eps`` — need ``rsvg-convert`` or ImageMagick ``magick``
          * ``.png`` — same converters (raster)
        """
        path = Path(path)
        fmt = (format or path.suffix.lstrip(".")).lower()
        if fmt.startswith("."):
            fmt = fmt[1:]
        if fmt in {"svg", "pdf", "ps", "eps", "png"} and not self.svg:
            raise ValueError(
                "This drawing has no SVG payload (for example a grid-only "
                "fallback). Draw a smaller network or pass highlight=path."
            )
        if fmt == "svg":
            path.write_text(self.svg)
            return path
        if fmt in {"pdf", "ps", "eps", "png"}:
            return _convert_svg(self.svg, path, fmt)
        raise ValueError(
            f"Unsupported format {fmt!r}. Use svg, pdf, ps, eps, or png."
        )


def _uid(prefix: str = "mf") -> str:
    return f"{prefix}{next(_IDS)}"


def _escape(text: str) -> str:
    return html.escape(str(text), quote=True)


def _rgba(hex_color: str, alpha: float) -> str:
    """Bake alpha into an RGB color (SVG markers ignore parent opacity)."""
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return hex_color
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha:g})"


def _rule_class_map() -> dict[str, str]:
    from xenosite.forest.rulesets import RULESETS, Full

    mapping = {
        rule.name: getattr(rule, "ruleset_name", None) or "other"
        for rule in Full.rules
    }
    for extra in ("TSO", "BA"):
        rs = RULESETS.get(extra)
        if rs is None:
            continue
        for rule in rs.rules:
            mapping.setdefault(rule.name, extra)
    return mapping


_RULE_CLASS = None


def rule_class(name: str) -> str:
    global _RULE_CLASS
    if _RULE_CLASS is None:
        _RULE_CLASS = _rule_class_map()
    key = name.split("_")[0]
    return _RULE_CLASS.get(name) or _RULE_CLASS.get(key) or "other"


def rule_color(name: str) -> str:
    return CLASS_STYLE.get(rule_class(name), _DEFAULT_STYLE)[0]


def edge_types(data: dict) -> list[str]:
    types: list[str] = []
    seen: set[str] = set()
    for rxn in data.get("rxn") or ():
        name = rxn.type.split("_")[0]
        if name not in seen:
            seen.add(name)
            types.append(name)
    return types or ["reaction"]


def _match_reaction(types: list[str], wanted: Iterable[str] | None) -> bool:
    if wanted is None:
        return True
    wanted_set = {w.split("_")[0] for w in wanted}
    return bool(wanted_set & set(types))


def generations(graph, roots: Iterable[str] | None = None) -> dict[str, int]:
    """BFS distance from source molecules (in-degree 0 by default)."""
    if roots is None:
        roots = [n for n in graph.nodes if graph.in_degree(n) == 0]
    roots = list(roots)
    if not roots and graph.number_of_nodes():
        roots = [next(iter(graph.nodes))]

    dist: dict[str, int] = {}
    queue = deque()
    for root in roots:
        if root not in graph:
            continue
        dist[root] = 0
        queue.append(root)

    while queue:
        node = queue.popleft()
        for child in graph.successors(node):
            nxt = dist[node] + 1
            if child not in dist or nxt < dist[child]:
                dist[child] = nxt
                queue.append(child)

    extra = (max(dist.values()) + 1) if dist else 0
    for node in graph.nodes:
        dist.setdefault(node, extra)
    return dist


def _canon_list(graph, items: Iterable[str] | None) -> list[str] | None:
    if items is None:
        return None
    if hasattr(graph, "canonize"):
        return [graph.canonize(x) for x in items]
    return list(items)


def filter_graph(
    graph,
    *,
    nodes: Iterable[str] | None = None,
    max_generation: int | None = None,
    reaction_types: Iterable[str] | None = None,
    path: Iterable[str] | None = None,
    roots: Iterable[str] | None = None,
):
    """Return a DiGraph copy filtered by generation, reaction type, path, or nodes."""
    import networkx as nx

    out = nx.DiGraph()
    roots_c = _canon_list(graph, roots)
    path_c = _canon_list(graph, path)
    nodes_c = _canon_list(graph, nodes)
    dist = generations(graph, roots=roots_c)

    keep = set(graph.nodes)
    if nodes_c is not None:
        keep &= set(nodes_c)
    if path_c is not None:
        keep &= set(path_c)
    if max_generation is not None:
        keep = {n for n in keep if dist.get(n, 0) <= max_generation}

    for n in keep:
        out.add_node(n, **dict(graph.nodes[n]))

    for u, v, data in graph.edges(data=True):
        if u not in out or v not in out:
            continue
        types = edge_types(data)
        if not _match_reaction(types, reaction_types):
            continue
        if path_c is not None:
            # keep only consecutive path edges when pruning to a path
            pairs = set(zip(path_c, path_c[1:]))
            if (u, v) not in pairs:
                continue
        out.add_edge(u, v, **dict(data))

    # Drop isolates that were not explicitly requested (except when path/nodes given).
    if nodes_c is None and path_c is None:
        isolates = [n for n in list(out.nodes) if out.degree(n) == 0 and dist.get(n, 0) > 0]
        out.remove_nodes_from(isolates)
    return out


def mol_svg(smiles: str, size: int = 140) -> str:
    mol = MolFromSmiles(smiles)
    if mol is None:
        empty = (
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{size}' height='{size}'>"
            f"<text x='50%' y='50%' text-anchor='middle'>{_escape(smiles)}</text></svg>"
        )
        return empty
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.14
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    text = drawer.GetDrawingText()
    start = text.find("<svg")
    end = text.rfind("</svg>")
    return text[start : end + 6]


def _wrap(inner: str, width: str = "100%") -> str:
    return (
        "<div style='font-family:Helvetica,Arial,sans-serif;font-size:13px;"
        f"color:#222;line-height:1.35;max-width:{width}'>"
        f"{inner}</div>"
    )


def _legend_html(classes: Iterable[str]) -> str:
    chips = []
    for cls in classes:
        color, label = CLASS_STYLE.get(cls, _DEFAULT_STYLE)
        chips.append(
            "<span style='display:inline-flex;align-items:center;gap:6px;"
            "margin:0 12px 6px 0'>"
            f"<span style='width:12px;height:12px;border-radius:50%;"
            f"background:{color};display:inline-block'></span>"
            f"<span>{_escape(label)}</span></span>"
        )
    if not chips:
        return ""
    return "<div style='margin:0 0 10px 0'>" + "".join(chips) + "</div>"


def _card(
    smiles: str,
    size: int,
    highlight: bool = False,
    show_smiles: bool = True,
    caption_max: int = 28,
    *,
    mol_background: str = "#ffffff",
    mol_border: str = "#cccccc",
    mol_border_width: float = 1.0,
    highlight_border: str = "#222222",
    highlight_border_width: float = 2.2,
) -> str:
    stroke = highlight_border if highlight else mol_border
    width = f"{highlight_border_width:g}px" if highlight else f"{mol_border_width:g}px"
    caption = ""
    if show_smiles:
        caption = smiles if len(smiles) <= caption_max else smiles[: caption_max - 1] + "…"
    cap_html = (
        f"<div style='font-size:11px;color:#555;margin-top:2px'>{_escape(caption)}</div>"
        if show_smiles
        else ""
    )
    return (
        "<div style='display:inline-block;text-align:center;"
        f"border:{width} solid {stroke};border-radius:10px;"
        f"background:{mol_background};"
        "padding:6px 6px 4px'>"
        f"{mol_svg(smiles, size)}"
        f"{cap_html}</div>"
    )


def _convert_svg(svg: str, path: Path, fmt: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "figure.svg"
        src.write_text(svg)
        rsvg = shutil.which("rsvg-convert")
        if rsvg:
            cmd = [rsvg, "-f", fmt, "-o", str(path), str(src)]
            subprocess.run(cmd, check=True, capture_output=True)
            return path
        magick = shutil.which("magick") or shutil.which("convert")
        if magick:
            cmd = [magick, str(src), str(path)]
            subprocess.run(cmd, check=True, capture_output=True)
            return path
    raise RuntimeError(
        f"Saving {fmt.upper()} requires librsvg (`rsvg-convert`) or ImageMagick "
        "(`magick`). Install one of those, or save as `.svg` instead."
    )


def draw_grid(
    graph,
    mols: Iterable[str] | None = None,
    mol_size: int = 120,
    columns: int = 5,
    *,
    max_generation: int | None = None,
    reaction_types: Iterable[str] | None = None,
    show_smiles: bool = True,
    show_generations: bool = True,
    caption_max: int = 28,
    roots: Iterable[str] | None = None,
    mol_background: str = "#ffffff",
    mol_border: str = "#cccccc",
    mol_border_width: float = 1.0,
) -> Drawing:
    """Grid of metabolite structures, ordered by generation."""
    view = filter_graph(
        graph,
        nodes=mols,
        max_generation=max_generation,
        reaction_types=reaction_types,
        roots=roots,
    )
    # If only filtering generations/types, still show all kept nodes even if isolated roots.
    if mols is None and reaction_types is None and max_generation is not None:
        # filter_graph already applied max_generation; re-add isolates at those gens
        dist_full = generations(graph, roots=_canon_list(graph, roots))
        for n, g in dist_full.items():
            if g <= max_generation and n not in view:
                view.add_node(n)

    dist = generations(view, roots=_canon_list(graph, roots))
    nodes = sorted(view.nodes, key=lambda n: (dist.get(n, 0), n))
    if mols is not None:
        nodes = _canon_list(graph, mols) or []

    cells = []
    for smiles in nodes:
        gen = dist.get(smiles)
        label = f"g{gen}" if (show_generations and gen is not None) else ""
        cells.append(
            "<div>"
            f"{_card(smiles, mol_size, show_smiles=show_smiles, caption_max=caption_max, mol_background=mol_background, mol_border=mol_border, mol_border_width=mol_border_width)}"
            + (
                f"<div style='font-size:10px;color:#888;margin-top:4px;text-align:center'>"
                f"{_escape(label)}</div>"
                if label
                else ""
            )
            + "</div>"
        )
    body = (
        f"<p style='margin:0 0 8px'><b>{len(nodes)}</b> metabolites</p>"
        "<div style='display:grid;"
        f"grid-template-columns:repeat({columns}, max-content);gap:12px'>"
        + "".join(cells)
        + "</div>"
    )
    return Drawing(_wrap(body))


def draw_path(
    graph,
    path: list[str],
    mol_size: int = 130,
    *,
    show_legend: bool = True,
    show_smiles: bool = True,
    caption_max: int = 28,
    mol_background: str = "#ffffff",
    mol_border: str = "#cccccc",
    mol_border_width: float = 1.0,
    highlight_border: str = "#222222",
    highlight_border_width: float = 2.2,
) -> Drawing:
    """Linear reaction scheme for one path through the network."""
    path = _canon_list(graph, path) or []
    if not path:
        return Drawing(_wrap("<i>Empty path.</i>"))

    card_kw = dict(
        show_smiles=show_smiles,
        caption_max=caption_max,
        mol_background=mol_background,
        mol_border=mol_border,
        mol_border_width=mol_border_width,
        highlight_border=highlight_border,
        highlight_border_width=highlight_border_width,
    )
    parts = [_card(path[0], mol_size, highlight=True, **card_kw)]
    used = []
    for a, b in zip(path, path[1:]):
        types = edge_types(graph.edges[a, b]) if graph.has_edge(a, b) else ["?"]
        used.extend(rule_class(t) for t in types)
        color = rule_color(types[0])
        label = ", ".join(types)
        parts.append(
            "<div style='display:flex;flex-direction:column;align-items:center;"
            "min-width:88px;padding:0 6px'>"
            f"<div style='font-size:11px;font-weight:600;color:{color};"
            f"white-space:nowrap'>{_escape(label)}</div>"
            f"<div style='font-size:22px;color:{color};line-height:1'>→</div>"
            "</div>"
        )
        parts.append(_card(b, mol_size, highlight=True, **card_kw))

    legend = _legend_html(dict.fromkeys(used)) if show_legend else ""
    body = (
        legend
        + "<div style='display:flex;align-items:center;flex-wrap:wrap;gap:4px'>"
        + "".join(parts)
        + "</div>"
    )
    return Drawing(_wrap(body))


def draw_network(
    graph,
    highlight: Iterable[str] | None = None,
    mol_size: int = 120,
    max_nodes: int = 48,
    *,
    max_generation: int | None = None,
    reaction_types: Iterable[str] | None = None,
    nodes: Iterable[str] | None = None,
    roots: Iterable[str] | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
    show_smiles: bool = True,
    show_generations: bool = True,
    show_title: bool = True,
    title: str | None = None,
    fade_others: bool = True,
    target_width: int = 980,
    background: str = "#ffffff",
    mol_background: str = "#ffffff",
    mol_border: str = "#cccccc",
    mol_border_width: float = 1.0,
    highlight_border: str = "#222222",
    highlight_border_width: float = 2.2,
    caption_max: int = 26,
) -> Drawing:
    """Layered SVG of the network: structures as nodes, Rainbow-colored edges.

    Defaults are publication-friendly (white figure and molecule panels).

    Parameters
    ----------
    highlight:
        Path (node SMILES list) to emphasize.
    mol_size:
        Pixel size of each structure drawing.
    max_nodes:
        If the graph is larger and nothing is highlighted, fall back to a grid.
    max_generation:
        Only include metabolites up to this BFS generation.
    reaction_types:
        Keep only edges whose reaction names intersect this set
        (for example ``["Dehydrogenation", "Epoxidation"]``).
    nodes:
        Explicit molecule allow-list (canonical SMILES or inputs to canonize).
    roots:
        Roots used to compute generations (default: in-degree 0).
    show_labels / show_legend / show_smiles / show_generations / show_title:
        Toggle figure chrome.
    title:
        Custom title string; ignored when ``show_title=False``.
    fade_others:
        When highlighting a path, dim off-path nodes and edges.
    target_width:
        Soft maximum width; molecule size shrinks to fit.
    background:
        Figure (SVG) background color. Default white for print/export.
    mol_background:
        Fill color of each molecule panel.
    mol_border / mol_border_width:
        Panel border color and stroke width for ordinary nodes.
    highlight_border / highlight_border_width:
        Panel border for nodes on a highlighted path.
    caption_max:
        Max SMILES characters under each node.
    """
    view = filter_graph(
        graph,
        nodes=nodes,
        max_generation=max_generation,
        reaction_types=reaction_types,
        roots=roots,
    )
    # When only max_generation is set, keep isolates in those generations.
    if nodes is None and reaction_types is None and max_generation is not None:
        dist_full = generations(graph, roots=_canon_list(graph, roots))
        for n, g in dist_full.items():
            if g <= max_generation and n not in view:
                view.add_node(n, **dict(graph.nodes[n]))

    n_nodes = view.number_of_nodes()
    if n_nodes == 0:
        return Drawing(_wrap("<i>Empty MetaboliteNetwork.</i>"))

    dist = generations(view, roots=_canon_list(graph, roots))
    layers: dict[int, list[str]] = defaultdict(list)
    for node, gen in dist.items():
        layers[gen].append(node)
    for gen in layers:
        layers[gen].sort()

    highlight_nodes = set()
    highlight_edges: set[tuple[str, str]] = set()
    if highlight is not None:
        path = _canon_list(graph, highlight) or []
        highlight_nodes = set(path)
        highlight_edges = set(zip(path, path[1:]))

    if n_nodes > max_nodes and not highlight_nodes:
        note = (
            f"<p><b>MetaboliteNetwork</b> has {n_nodes} molecules and "
            f"{view.number_of_edges()} reactions — too many to draw as a graph. "
            "Showing a structure grid instead. Pass <code>highlight=path</code>, "
            "<code>max_generation=...</code>, <code>reaction_types=...</code>, "
            "or prune the network first.</p>"
        )
        grid = draw_grid(
            view,
            mol_size=max(72, mol_size - 20),
            show_smiles=show_smiles,
            show_generations=show_generations,
            caption_max=caption_max,
            roots=roots,
        )
        return Drawing(_wrap(note) + grid.html)

    if highlight_nodes and n_nodes > max_nodes:
        keep = set(highlight_nodes)
        for u, v in view.edges:
            if u in highlight_nodes or v in highlight_nodes:
                keep.add(u)
                keep.add(v)
        layers = defaultdict(list)
        for node in keep:
            if node in dist:
                layers[dist[node]].append(node)
        for gen in layers:
            layers[gen].sort()

    pad = 8
    caption_h = 18 if show_smiles else 4
    hgap = 18
    vgap = 70 if show_labels else 48
    left = 52 if show_generations else 24
    top = 62 if (show_title or show_legend) else 24
    max_in_layer = max(len(v) for v in layers.values())
    est = max_in_layer * (mol_size + pad * 2 + hgap)
    if est > target_width:
        mol_size = max(70, int(target_width / max_in_layer) - pad * 2 - hgap)
    card_w = mol_size + pad * 2
    card_h = mol_size + caption_h + pad * 2
    width = left + max_in_layer * (card_w + hgap) - hgap + 24
    height = top + (max(layers) + 1) * (card_h + vgap) - vgap + 28

    pos: dict[str, tuple[float, float]] = {}
    for gen, group in layers.items():
        span = len(group) * (card_w + hgap) - hgap
        x0 = left + (width - left - 24 - span) / 2
        y = top + gen * (card_h + vgap)
        for i, node in enumerate(group):
            pos[node] = (x0 + i * (card_w + hgap), y)

    used_classes: list[str] = []
    marker_id = _uid("arr")
    colors_needed: dict[str, str] = {}
    faded_colors: set[str] = set()

    edge_svg = []
    for u, v, data in view.edges(data=True):
        if u not in pos or v not in pos:
            continue
        types = edge_types(data)
        cls = rule_class(types[0])
        color = rule_color(types[0])
        colors_needed[color] = color
        if cls not in used_classes:
            used_classes.append(cls)

        on_path = (u, v) in highlight_edges
        # Dim every edge that is not on the highlighted path (including
        # other reactions leaving a path molecule). Bake alpha into stroke
        # and marker fill — SVG markers often ignore parent opacity.
        faded = bool(highlight_nodes) and fade_others and not on_path
        alpha = 0.18 if faded else 1.0
        stroke = _rgba(color, alpha)
        stroke_w = 2.6 if on_path else 1.6
        color_key = color.lstrip("#")
        marker_key = f"{color_key}f" if faded else color_key
        if faded:
            faded_colors.add(color)

        x1, y1 = pos[u]
        x2, y2 = pos[v]
        sx = x1 + card_w / 2
        sy = y1 + card_h
        tx = x2 + card_w / 2
        ty = y2
        mid_y = (sy + ty) / 2
        path_d = (
            f"M {sx:.1f} {sy:.1f} C {sx:.1f} {mid_y:.1f}, "
            f"{tx:.1f} {mid_y:.1f}, {tx:.1f} {ty:.1f}"
        )
        edge_svg.append(
            f"<path d='{path_d}' fill='none' stroke='{stroke}' "
            f"stroke-width='{stroke_w}' "
            f"marker-end='url(#{marker_id}{marker_key})'/>"
        )
        if show_labels and not faded:
            label = ", ".join(types)
            lx, ly = tx, ty - 14
            edge_svg.append(
                f"<text x='{lx:.1f}' y='{ly:.1f}' text-anchor='middle' "
                f"font-size='11' font-family='Helvetica,Arial,sans-serif' "
                f"font-weight='600' fill='{color}' "
                f"stroke='#fff' stroke-width='3' paint-order='stroke'>"
                f"{_escape(label)}</text>"
            )

    markers = []
    for color in colors_needed:
        key = color.lstrip("#")
        markers.append(
            f"<marker id='{marker_id}{key}' viewBox='0 0 10 10' refX='8' refY='5' "
            f"markerWidth='7' markerHeight='7' orient='auto-start-reverse'>"
            f"<path d='M 0 0 L 10 5 L 0 10 z' fill='{color}'/></marker>"
        )
        if color in faded_colors:
            markers.append(
                f"<marker id='{marker_id}{key}f' viewBox='0 0 10 10' refX='8' refY='5' "
                f"markerWidth='7' markerHeight='7' orient='auto-start-reverse'>"
                f"<path d='M 0 0 L 10 5 L 0 10 z' fill='{_rgba(color, 0.18)}'/></marker>"
            )

    node_svg = []
    for node, (x, y) in pos.items():
        on = node in highlight_nodes if highlight_nodes else False
        faded = bool(highlight_nodes) and fade_others and not on
        opacity = 0.22 if faded else 1.0
        stroke = highlight_border if on else mol_border
        sw = highlight_border_width if on else mol_border_width
        caption = ""
        if show_smiles:
            caption = node if len(node) <= caption_max else node[: caption_max - 1] + "…"
        inner = mol_svg(node, mol_size)
        inner_body_start = inner.find(">")
        inner_content = inner[inner_body_start + 1 : inner.rfind("</svg>")]
        caption_svg = ""
        if show_smiles:
            caption_svg = (
                f"<text x='{x + card_w / 2:.1f}' y='{y + card_h - 6:.1f}' "
                f"text-anchor='middle' font-size='11' "
                f"font-family='Helvetica,Arial,sans-serif' fill='#555'>"
                f"{_escape(caption)}</text>"
            )
        node_svg.append(
            f"<g opacity='{opacity}'>"
            f"<rect x='{x:.1f}' y='{y:.1f}' width='{card_w}' height='{card_h}' "
            f"rx='10' fill='{mol_background}' stroke='{stroke}' stroke-width='{sw}'/>"
            f"<svg x='{x + pad:.1f}' y='{y + pad:.1f}' width='{mol_size}' "
            f"height='{mol_size}' viewBox='0 0 {mol_size} {mol_size}'>"
            f"{inner_content}</svg>"
            f"{caption_svg}</g>"
        )

    gen_labels = []
    if show_generations:
        for gen in sorted(layers):
            y = pos[layers[gen][0]][1] + card_h / 2
            gen_labels.append(
                f"<text x='16' y='{y:.1f}' font-size='11' fill='#888' "
                f"font-family='Helvetica,Arial,sans-serif' dominant-baseline='middle'>"
                f"g{gen}</text>"
            )

    order = list(CLASS_STYLE)
    used_classes.sort(key=lambda c: order.index(c) if c in CLASS_STYLE else 99)

    legend_bits = []
    if show_legend:
        lx = 16
        for cls in used_classes:
            color, label = CLASS_STYLE.get(cls, _DEFAULT_STYLE)
            legend_bits.append(
                f"<circle cx='{lx}' cy='40' r='5' fill='{color}'/>"
                f"<text x='{lx + 10}' y='44' font-size='12' "
                f"font-family='Helvetica,Arial,sans-serif' fill='#333'>"
                f"{_escape(label)}</text>"
            )
            lx += 7.2 * len(label) + 28

    if title is None:
        title = (
            f"{view.number_of_nodes()} metabolites · {view.number_of_edges()} reactions"
            if graph.number_of_nodes() == len(pos)
            else f"Showing {len(pos)} of {graph.number_of_nodes()} metabolites"
        )
    title_svg = ""
    if show_title and title:
        title_svg = (
            f"<text x='16' y='20' font-size='13' "
            f"font-family='Helvetica,Arial,sans-serif' fill='#444'>"
            f"{_escape(title)}</text>"
        )

    svg = (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{int(width)}' "
        f"height='{int(height)}' viewBox='0 0 {int(width)} {int(height)}'>"
        f"<rect width='100%' height='100%' fill='{background}' rx='8'/>"
        f"<defs>{''.join(markers)}</defs>"
        + title_svg
        + "".join(legend_bits)
        + "".join(gen_labels)
        + "".join(edge_svg)
        + "".join(node_svg)
        + "</svg>"
    )
    return Drawing(_wrap(svg), svg=svg)
