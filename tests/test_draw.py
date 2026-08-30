import pytest

nx = pytest.importorskip("networkx")

from xenosite.forest import rules, RuleSet
from xenosite.forest.draw import rule_class, rule_color, generations
from xenosite.forest.net import MetaboliteNetwork
from xenosite.forest.phaseone import PhaseOneRS


def test_rule_rainbow_colors():
    assert rule_class("Hydroxylation") == "SO"
    assert rule_class("Dealkylation") == "UO"
    assert rule_class("Dehydrogenation") == "DH"
    assert rule_color("Dehydrogenation") == "#009E73"
    assert rule_color("Hydroxylation") == "#D55E00"


def test_draw_network_and_path():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)

    drawing = net.draw()
    assert drawing.svg is not None
    assert "CCO" in drawing.svg
    assert "#009E73" in drawing.svg  # DH
    assert "#D55E00" in drawing.svg  # SO

    path = next(net.paths("CCO", "CC=O"))
    scheme = net.draw_path(path)
    assert "Dehydrogenation" in scheme.html
    assert "CC=O" in scheme.html

    grid = net.grid()
    assert "CCO" in grid.html
    assert "g0" in grid.html

    highlighted = net.draw(highlight=path)
    assert "stroke-width='2.6'" in highlighted.svg
    # Off-path reactions that leave a path molecule must still be dimmed.
    # Alpha is baked into stroke/marker colors (markers ignore parent opacity).
    assert "rgba(" in highlighted.svg and ",0.18)" in highlighted.svg

    assert "Empty" in MetaboliteNetwork().draw().html


def test_network_repr_html():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)
    html = net._repr_html_()
    assert "<svg" in html
    assert "CCO" in html


def test_draw_several_steps_deep():
    """Keep the ruleset tiny so multi-generation graphs stay drawable."""
    rs = RuleSet([rules.Dehydrogenation(), rules.Epoxidation()], name="tiny")
    net = MetaboliteNetwork("C=CC")
    net.expand(rs)
    net.expand(rs)

    dist = generations(net)
    assert max(dist.values()) >= 2
    assert net.number_of_nodes() < 48

    drawing = net.draw()
    assert drawing.svg is not None
    assert ">g0<" in drawing.svg
    assert ">g2<" in drawing.svg

    path = next(net.paths("C=CC", "C=C1CO1", extra_depth=2))
    assert len(path) >= 3
    assert path[0] == "C=CC"
    # Substrate epoxidation is an off-path side reaction.
    assert net.has_edge("C=CC", "CC1CO1")
    scheme = net.draw_path(path)
    assert scheme.html.count("→") == len(path) - 1
    highlighted = net.draw(highlight=path)
    assert highlighted.svg is not None
    assert "stroke-width='2.6'" in highlighted.svg
    assert ",0.18)" in highlighted.svg
    # Path epoxidation keeps a full-color marker; substrate epoxidation uses faded.
    assert f"fill='#D55E00'" in highlighted.svg or "fill='#D55E00'" in highlighted.svg
    assert "rgba(213,94,0,0.18)" in highlighted.svg


def test_draw_options():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)

    only_dh = net.draw(reaction_types=["Dehydrogenation"], title="DH only")
    assert only_dh.svg is not None
    assert "DH only" in only_dh.svg
    assert "Dehydrogenation" in only_dh.svg
    assert "Hydroxylation" not in only_dh.svg

    no_chrome = net.draw(
        show_labels=False,
        show_legend=False,
        show_smiles=False,
        show_generations=False,
        show_title=False,
        mol_size=80,
    )
    assert no_chrome.svg is not None
    assert "Stable oxygenation" not in no_chrome.svg
    assert ">g0<" not in no_chrome.svg

    g0 = net.draw(max_generation=0)
    assert "CCO" in g0.svg
    assert g0.svg.count("<rect ") >= 1

    # Publication defaults: white figure background.
    assert 'fill="#ffffff"' in net.draw().svg or "fill='#ffffff'" in net.draw().svg

    subset = net.draw(nodes=["CCO", "CC=O", "C=CO"], title="Three molecules")
    assert "CCO" in subset.svg and "CC=O" in subset.svg
    assert "OCCO" not in subset.svg

    styled = net.draw(
        nodes=["CCO", "CC=O"],
        background="#f7f7f7",
        mol_background="#fffaf0",
        mol_border="#888888",
        mol_border_width=2,
        highlight=["CCO", "CC=O"],
        highlight_border="#d62728",
        highlight_border_width=3,
    )
    assert "fill='#f7f7f7'" in styled.svg or 'fill="#f7f7f7"' in styled.svg
    assert "fill='#fffaf0'" in styled.svg or 'fill="#fffaf0"' in styled.svg
    assert "#d62728" in styled.svg


def test_prune_network():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)
    n0 = net.number_of_nodes()

    dh = net.prune(reaction_types=["Dehydrogenation"])
    assert dh.number_of_nodes() < n0
    assert all(
        "Dehydrogenation" in {r.type for r in data["rxn"]}
        for _, _, data in dh.edges(data=True)
    )

    path = next(net.paths("CCO", "CC=O"))
    slim = net.prune(path=path)
    assert set(slim.nodes) == set(path)
    assert slim.number_of_edges() == len(path) - 1

    g0 = net.prune(max_generation=0)
    assert set(g0.nodes) == {"CCO"}

    # copy is independent
    cloned = net.copy()
    cloned.prune(path=path, inplace=True)
    assert set(cloned.nodes) == set(path)
    assert net.number_of_nodes() == n0


def test_save_svg_and_pdf(tmp_path):
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)
    svg_path = tmp_path / "net.svg"
    net.save_draw(svg_path, reaction_types=["Dehydrogenation"])
    assert svg_path.exists()
    assert "Dehydrogenation" in svg_path.read_text()

    pdf_path = tmp_path / "net.pdf"
    try:
        net.save_draw(pdf_path, reaction_types=["Dehydrogenation"])
    except RuntimeError as err:
        pytest.skip(str(err))
    assert pdf_path.exists() and pdf_path.stat().st_size > 100

    ps_path = tmp_path / "net.ps"
    try:
        net.draw(reaction_types=["Dehydrogenation"]).save(ps_path)
    except RuntimeError as err:
        pytest.skip(str(err))
    assert ps_path.exists() and ps_path.stat().st_size > 100
