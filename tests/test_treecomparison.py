import pytest
from PIL import Image

from treecomparison.cli import main
from treecomparison.comparison import (
    compare_trees,
    parse_newick,
    tree_to_dot,
)


def test_parse_newick_preserves_topology_labels_and_branch_lengths():
    """Confirm the parser reads a labelled tree with branch lengths.

    :return: None.
    """
    root = parse_newick("((A:1,'Taxon B':2)AB:3,C:4);")

    assert len(root.children) == 2
    assert root.children[0].name == "AB"
    assert root.children[0].branch_length == 3.0
    assert [child.name for child in root.children[0].children] == ["A", "Taxon B"]
    assert root.children[1].name == "C"
    assert root.children[1].branch_length == 4.0


def test_tree_to_dot_uses_horizontal_root_to_leaf_layout():
    """Confirm generated Graphviz source lays trees out from left to right.

    :return: None.
    """
    dot = tree_to_dot(
        parse_newick("(A:1,B:1);"),
        "Source tree",
        "strict_clock",
    )

    assert "rankdir=LR" in dot
    assert "Source tree" in dot
    assert "strict_clock" in dot
    assert 'POINT-SIZE="11"' in dot
    assert 'label="A"' in dot
    assert 'label="B"' in dot
    assert 'fontcolor="white"' in dot


@pytest.mark.parametrize(
    "newick,error",
    [
        ("", "must not be empty"),
        ("(A:1);", "at least two children"),
        ("(A:bad,B:1);", "Invalid Newick branch length"),
        ("(A:1,B:1)", "Expected ';'"),
        ("(A:1,B:1); trailing", "Unexpected content"),
    ],
)
def test_parse_newick_rejects_invalid_input(newick, error):
    """Confirm malformed Newick input raises a descriptive error.

    :param newick: Invalid Newick document.
    :param error: Expected error-message fragment.
    :return: None.
    """
    with pytest.raises(ValueError, match=error):
        parse_newick(newick)


def test_compare_trees_writes_side_by_side_png(tmp_path):
    """Confirm comparison rendering creates a valid landscape-style PNG.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    source_path = tmp_path / "source.newick"
    reconstructed_path = tmp_path / "reconstructed.newick"
    output_path = tmp_path / "images" / "comparison.png"
    source_path.write_text("((A:1,B:1):1,C:2);\n", encoding="utf-8")
    reconstructed_path.write_text("(A:2,(B:1,C:1):1);\n", encoding="utf-8")

    assert compare_trees(source_path, reconstructed_path, output_path) == output_path
    with Image.open(output_path) as image:
        # Horizontal panel composition should be wider than either tree's typical height.
        assert image.format == "PNG"
        assert image.width > image.height
        assert image.mode == "RGB"
        # The comparison uses a six-pixel outer frame instead of Graphviz's wide transparency.
        assert image.getpixel((0, 0)) == (0, 0, 0)
        assert image.getpixel((5, 5)) == (0, 0, 0)
        assert image.getpixel((6, 6)) == (255, 255, 255)


def test_cli_accepts_required_short_arguments(tmp_path, capsys):
    """Confirm the CLI accepts source, reconstructed, and output short options.

    :param tmp_path: Temporary directory supplied by pytest.
    :param capsys: Pytest fixture used to capture command output.
    :return: None.
    """
    source_path = tmp_path / "source.newick"
    reconstructed_path = tmp_path / "reconstructed.newick"
    output_path = tmp_path / "comparison.png"
    source_path.write_text("(A:1,B:1);", encoding="utf-8")
    reconstructed_path.write_text("(A:1,B:1);", encoding="utf-8")

    assert main([
        "-s",
        str(source_path),
        "-r",
        str(reconstructed_path),
        "-o",
        str(output_path),
    ]) == 0
    assert output_path.exists()
    assert f"Wrote tree comparison: {output_path}" in capsys.readouterr().out


def test_compare_trees_requires_png_extension(tmp_path):
    """Confirm the output filename identifies the required PNG format.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    source_path = tmp_path / "source.newick"
    reconstructed_path = tmp_path / "reconstructed.newick"
    source_path.write_text("(A,B);", encoding="utf-8")
    reconstructed_path.write_text("(A,B);", encoding="utf-8")

    with pytest.raises(ValueError, match="must end with .png"):
        compare_trees(source_path, reconstructed_path, tmp_path / "comparison.jpg")
