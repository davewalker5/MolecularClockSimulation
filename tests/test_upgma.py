import json

import pytest

from phylogeny.cli import main
from phylogeny.upgma import (
    Cluster,
    cluster_distance,
    load_distance_matrix,
    to_newick,
    upgma,
)


FOUR_TAXON_LABELS = ["A", "B", "C", "D"]
FOUR_TAXON_MATRIX = [
    [0, 2, 6, 8],
    [2, 0, 6, 8],
    [6, 6, 0, 4],
    [8, 8, 4, 0],
]


def test_two_taxon_tree_has_half_distance_branch_lengths():
    """Confirm a two-taxon distance is divided equally between both branches.

    :return: None.
    """
    root = upgma(["A", "B"], [[0, 2], [2, 0]])

    assert root.height == 1.0
    assert to_newick(root) == "(A:1.000000,B:1.000000);"


def test_four_taxon_tree_has_expected_clustering_and_branch_lengths():
    """Confirm the worked example produces the specified UPGMA tree.

    :return: None.
    """
    root = upgma(FOUR_TAXON_LABELS, FOUR_TAXON_MATRIX)

    assert root.members == ("A", "B", "C", "D")
    assert root.height == 3.5
    assert to_newick(root) == (
        "((A:1.000000,B:1.000000):2.500000,"
        "(C:2.000000,D:2.000000):1.500000);"
    )


def test_cluster_distance_averages_every_member_pair():
    """Confirm cluster distance is the arithmetic mean of cross-cluster pairs.

    :return: None.
    """
    first = Cluster(name="A+B", members=("A", "B"), height=1.0)
    second = Cluster(name="C+D", members=("C", "D"), height=2.0)

    # The four cross-cluster distances are 6, 8, 6, and 8.
    assert cluster_distance(
        first,
        second,
        FOUR_TAXON_LABELS,
        FOUR_TAXON_MATRIX,
    ) == 7.0


def test_reconstructed_tree_is_ultrametric():
    """Confirm every leaf lies at the same distance from the root.

    :return: None.
    """
    root = upgma(FOUR_TAXON_LABELS, FOUR_TAXON_MATRIX)

    def root_to_tip_lengths(cluster, distance=0.0):
        """Collect cumulative distances from the root to each descendant leaf.

        :param cluster: Current cluster in the traversal.
        :param distance: Accumulated distance from the root.
        :return: Mapping of leaf labels to root-to-tip distances.
        """
        # A leaf ends the current traversal path.
        if cluster.left is None and cluster.right is None:
            return {cluster.name: distance}
        return {
            **root_to_tip_lengths(
                cluster.left,
                distance + cluster.left_branch_length,
            ),
            **root_to_tip_lengths(
                cluster.right,
                distance + cluster.right_branch_length,
            ),
        }

    assert root_to_tip_lengths(root) == {
        "A": 3.5,
        "B": 3.5,
        "C": 3.5,
        "D": 3.5,
    }


def test_equal_distances_use_alphabetical_tie_breaking():
    """Confirm tied candidate pairs are selected in alphabetical order.

    :return: None.
    """
    root = upgma(
        ["C", "B", "A"],
        [
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0],
        ],
    )

    assert root.left.members == ("A", "B")
    assert to_newick(root) == "((A:1.000000,B:1.000000):0.000000,C:1.000000);"


@pytest.mark.parametrize(
    "labels,matrix,error",
    [
        (["A"], [[0]], "at least two taxa"),
        (["A", "B"], [[0, 1]], "square"),
        (["A", "B"], [[1, 2], [2, 0]], "Diagonal"),
        (["A", "B"], [[0, 1], [2, 0]], "symmetric"),
        (["A", "B"], [[0, -1], [-1, 0]], "non-negative"),
        (["A", "B"], [[0, float("inf")], [float("inf"), 0]], "finite"),
    ],
)
def test_invalid_matrices_raise_clear_errors(labels, matrix, error):
    """Confirm malformed or unusable matrices are rejected clearly.

    :param labels: Taxon labels supplied to reconstruction.
    :param matrix: Invalid matrix supplied to reconstruction.
    :param error: Expected validation-message fragment.
    :return: None.
    """
    with pytest.raises(ValueError, match=error):
        upgma(labels, matrix)


def test_load_distance_matrix_reads_calculator_json(tmp_path):
    """Confirm the loader accepts the existing calculator's JSON structure.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    input_path = tmp_path / "distance_matrix.json"
    input_path.write_text(
        json.dumps({
            "distance_metric": "hamming",
            "labels": ["A", "B"],
            "matrix": [[0, 2], [2, 0]],
        }),
        encoding="utf-8",
    )

    assert load_distance_matrix(input_path) == (["A", "B"], [[0, 2], [2, 0]])


def test_cli_writes_newick_output(tmp_path, capsys):
    """Confirm the CLI reconstructs and writes a tree from calculator JSON.

    :param tmp_path: Temporary directory supplied by pytest.
    :param capsys: Pytest fixture used to capture command output.
    :return: None.
    """
    input_path = tmp_path / "distance_matrix.json"
    output_path = tmp_path / "trees" / "upgma_tree.nwk"
    input_path.write_text(
        json.dumps({
            "labels": ["A", "B"],
            "matrix": [[0, 2], [2, 0]],
            "distance_metric": "hamming",
        }),
        encoding="utf-8",
    )

    assert main(["--input", str(input_path), "--output", str(output_path)]) == 0
    assert output_path.read_text(encoding="utf-8") == "(A:1.000000,B:1.000000);\n"
    assert f"Wrote Newick tree: {output_path}" in capsys.readouterr().out
