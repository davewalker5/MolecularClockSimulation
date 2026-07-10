"""Tests for reusable and command-line tree calibration."""

import json

import pytest

from treecalibration.calibration import calibrate_tree, to_newick
from treecalibration.cli import main
from treecomparison.comparison import parse_newick
from molecular_clock_simulation.calibration_ui import clear_calibration_state
from molecular_clock_simulation.reconstruction import calibrated_tree_to_dot
from common import DARK_THEME


def test_calibrate_tree_scales_copy_from_mrca_tip_depth():
    """Confirm calibration scales every branch without mutating its input.

    :return: None.
    """
    tree = parse_newick("((A:0.04,B:0.04):0.06,(C:0.03,D:0.03):0.07);")
    calibrated, metadata = calibrate_tree(tree, ["A", "B"], 20.0)

    assert metadata["reconstructed_depth"] == pytest.approx(0.04)
    assert metadata["scale_factor"] == pytest.approx(500.0)
    assert calibrated.children[0].children[0].branch_length == pytest.approx(20.0)
    assert calibrated.children[0].branch_length == pytest.approx(30.0)
    assert tree.children[0].children[0].branch_length == pytest.approx(0.04)


@pytest.mark.parametrize(
    "taxa,age,error",
    [
        (["A"], 10, "at least two"),
        (["A", "A"], 10, "unique"),
        (["A", "missing"], 10, "not found"),
        (["A", "B"], 0, "greater than zero"),
    ],
)
def test_calibrate_tree_rejects_invalid_calibrations(taxa, age, error):
    """Confirm invalid calibration values produce informative errors.

    :param taxa: Calibration taxa under test.
    :param age: Calibration age under test.
    :param error: Expected error-message fragment.
    :return: None.
    """
    with pytest.raises(ValueError, match=error):
        calibrate_tree(parse_newick("(A:1,B:1);"), taxa, age)


def test_calibrate_tree_requires_all_branch_lengths():
    """Confirm incomplete branch-length trees are rejected.

    :return: None.
    """
    with pytest.raises(ValueError, match="Every non-root branch"):
        calibrate_tree(parse_newick("(A:1,B);"), ["A", "B"], 10)


def test_calibrate_tree_accepts_negative_neighbor_joining_branch_lengths():
    """Confirm finite signed NJ branches are retained during calibration.

    :return: None.
    """
    tree = parse_newick("((A:0.002,B:-0.00002):0.003,C:0.005);")

    calibrated, metadata = calibrate_tree(tree, ["A", "B"], 10)

    assert metadata["reconstructed_depth"] == pytest.approx(0.00099)
    assert calibrated.children[0].children[1].branch_length < 0


def test_cli_writes_calibrated_tree_and_report(tmp_path):
    """Confirm the CLI writes both specified output artifacts.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    tree_path = tmp_path / "tree.nwk"
    calibration_path = tmp_path / "calibration.json"
    output_path = tmp_path / "output" / "dated_tree"
    tree_path.write_text("((A:2,B:2):1,C:3);", encoding="utf-8")
    calibration_path.write_text(json.dumps({"node": {"taxa": ["A", "B"]}, "age_mya": 10}), encoding="utf-8")

    assert main(["-i", str(tree_path), "-c", str(calibration_path), "-o", str(output_path)]) == 0
    assert (tmp_path / "output" / "dated_tree.newick").read_text(encoding="utf-8") == "((A:10,B:10):5,C:15);\n"
    report = json.loads((tmp_path / "output" / "dated_tree_metadata.json").read_text(encoding="utf-8"))
    assert report["units"] == "million years"
    assert report["scale_factor"] == pytest.approx(5.0)


def test_to_newick_preserves_quoted_labels():
    """Confirm serialization protects labels containing delimiters.

    :return: None.
    """
    assert to_newick(parse_newick("('Taxon A':1,B:1);")) == "('Taxon A':1,B:1);\n"


def test_clear_calibration_state_removes_results_and_widget_values():
    """Confirm changing an upstream tree invalidates all calibration state.

    :return: None.
    """
    state = {
        "strict_calibration_newick": "(A:1,B:1);",
        "strict_calibration_taxon_a": "A",
        "strict_calibration_taxon_b": "B",
        "strict_calibration_age_mya": 10.0,
        "unrelated": "preserved",
    }

    clear_calibration_state(state, "strict")

    assert state == {"unrelated": "preserved"}


def test_calibrated_tree_dot_matches_explorer_tree_styling():
    """Confirm calibrated trees use the explorers' dark reconstruction format.

    :return: None.
    """
    dot = calibrated_tree_to_dot(
        parse_newick("((A:10,B:10):5,C:15);"),
        graph_name="calibrated_tree",
        colors=DARK_THEME,
    )

    assert "rankdir=LR" in dot
    assert DARK_THEME["page_bg"] in dot
    assert 'shape=box, style="rounded,filled"' in dot
    assert 'label="A"' in dot
    assert 'label="10"' in dot
