import json

import pytest

from molecular_clock_simulation.reconstruction import (
    reconstructed_tree_newick,
    reconstructed_tree_to_dot,
    reconstruct_tree,
)
from phylogeny.cli import main
from phylogeny.neighbor_joining import (
    final_branch_lengths,
    neighbor_branch_lengths,
    neighbor_joining,
    q_matrix,
    select_neighbor_pair,
    total_distances,
    updated_distance_map,
)
from phylogeny.upgma import Cluster, to_newick


PUBLISHED_LABELS = ["A", "B", "C", "D", "E"]
PUBLISHED_MATRIX = [
    [0, 5, 9, 9, 8],
    [5, 0, 10, 10, 9],
    [9, 10, 0, 8, 7],
    [9, 10, 8, 0, 3],
    [8, 9, 7, 3, 0],
]


def published_clusters() -> list[Cluster]:
    """Build active leaf clusters for the standard NJ worked example.

    :return: Deterministically ordered leaf clusters.
    """
    return [
        Cluster(name=label, members=(label,), height=0.0)
        for label in PUBLISHED_LABELS
    ]


def published_distances() -> dict[tuple[str, str], float]:
    """Build active distances for the standard NJ worked example.

    :return: Pairwise distance map keyed by sorted taxon pairs.
    """
    return {
        ("A", "B"): 5,
        ("A", "C"): 9,
        ("A", "D"): 9,
        ("A", "E"): 8,
        ("B", "C"): 10,
        ("B", "D"): 10,
        ("B", "E"): 9,
        ("C", "D"): 8,
        ("C", "E"): 7,
        ("D", "E"): 3,
    }


def test_q_matrix_matches_known_neighbor_joining_example():
    """Confirm Q-matrix values match the standard five-taxon NJ example.

    :return: None.
    """
    clusters = published_clusters()
    distances = published_distances()

    q_values = q_matrix(clusters, distances)

    assert total_distances(clusters, distances) == {
        "A": 31,
        "B": 34,
        "C": 34,
        "D": 30,
        "E": 27,
    }
    assert q_values[("A", "B")] == -50
    assert q_values[("D", "E")] == -48
    assert select_neighbor_pair(q_values, clusters)[0].name == "A"
    assert select_neighbor_pair(q_values, clusters)[1].name == "B"


def test_neighbor_joining_branch_lengths_match_known_example():
    """Confirm the first NJ join assigns unequal child branch lengths correctly.

    :return: None.
    """
    clusters = published_clusters()
    distances = published_distances()

    left_length, right_length = neighbor_branch_lengths(
        clusters[0],
        clusters[1],
        clusters,
        distances,
    )

    assert left_length == pytest.approx(2.0)
    assert right_length == pytest.approx(3.0)


def test_neighbor_joining_distance_update_does_not_weight_cluster_size():
    """Confirm NJ updates distances with the pairwise formula rather than UPGMA averaging.

    :return: None.
    """
    clusters = published_clusters()
    distances = published_distances()
    merged = Cluster(name="A+B", members=("A", "B"), height=0.0)

    updated = updated_distance_map(clusters[0], clusters[1], merged, clusters, distances)

    assert updated[("A+B", "C")] == pytest.approx(7.0)
    assert updated[("A+B", "D")] == pytest.approx(7.0)
    assert updated[("A+B", "E")] == pytest.approx(6.0)
    assert updated[("D", "E")] == pytest.approx(3.0)


def test_neighbor_joining_final_three_way_branch_lengths_match_known_example():
    """Confirm the final NJ node keeps three unrooted branches.

    :return: None.
    """
    ab = Cluster(name="A+B", members=("A", "B"), height=0.0)
    c = Cluster(name="C", members=("C",), height=0.0)
    de = Cluster(name="D+E", members=("D", "E"), height=0.0)
    distances = {
        ("A+B", "C"): 7.0,
        ("A+B", "D+E"): 5.0,
        ("C", "D+E"): 6.0,
    }

    assert final_branch_lengths(ab, c, de, distances) == pytest.approx((3.0, 4.0, 2.0))


def test_neighbor_joining_reconstructs_known_published_tree():
    """Confirm the standard five-taxon example produces deterministic NJ Newick.

    :return: None.
    """
    root = neighbor_joining(PUBLISHED_LABELS, PUBLISHED_MATRIX)

    assert root.members == ("A", "B", "C", "D", "E")
    assert root.third is not None
    assert to_newick(root) == (
        "(((A:2.000000,B:3.000000):3.000000,C:4.000000):2.000000,"
        "D:2.000000,E:1.000000);"
    )


def test_neighbor_joining_final_node_renders_in_explorer_dot():
    """Confirm three-way NJ roots render through the shared explorer DOT helper.

    :return: None.
    """
    root = neighbor_joining(PUBLISHED_LABELS, PUBLISHED_MATRIX)
    colors = {
        "page_bg": "#000000",
        "surface_elevated": "#111111",
        "border_strong": "#222222",
        "text": "#eeeeee",
        "text_subtle": "#999999",
        "text_muted_strong": "#cccccc",
    }

    dot = reconstructed_tree_to_dot(root, graph_name="nj_tree", colors=colors)

    assert dot.startswith("digraph nj_tree")
    assert '"A+B+C+D+E" -> "E" [label="1"]' in dot


def test_neighbor_joining_two_taxon_tree_splits_final_distance():
    """Confirm the final two-cluster join divides the remaining distance equally.

    :return: None.
    """
    root = neighbor_joining(["A", "B"], [[0, 2], [2, 0]])

    assert to_newick(root) == "(A:1.000000,B:1.000000);"


def test_neighbor_joining_uses_deterministic_tie_breaking():
    """Confirm tied Q values select the alphabetically earliest cluster pair.

    :return: None.
    """
    root = neighbor_joining(
        ["C", "B", "A"],
        [
            [0, 2, 2],
            [2, 0, 2],
            [2, 2, 0],
        ],
    )

    assert to_newick(root) == "(A:1.000000,B:1.000000,C:1.000000);"
    assert to_newick(root) == to_newick(neighbor_joining(["C", "B", "A"], [
        [0, 2, 2],
        [2, 0, 2],
        [2, 2, 0],
    ]))


@pytest.mark.parametrize(
    "labels,matrix,error",
    [
        (["A"], [[0]], "at least two taxa"),
        (["A", "B"], [[0, 1]], "square"),
        (["A", "B"], [[1, 2], [2, 0]], "Diagonal"),
        (["A", "B"], [[0, 1], [2, 0]], "symmetric"),
        (["A", "B"], [[0, -1], [-1, 0]], "non-negative"),
    ],
)
def test_neighbor_joining_reuses_distance_matrix_validation(labels, matrix, error):
    """Confirm malformed matrices fail before NJ reconstruction begins.

    :param labels: Taxon labels supplied to reconstruction.
    :param matrix: Invalid matrix supplied to reconstruction.
    :param error: Expected validation-message fragment.
    :return: None.
    """
    with pytest.raises(ValueError, match=error):
        neighbor_joining(labels, matrix)


def test_reconstruct_tree_dispatches_neighbor_joining_from_payload():
    """Confirm shared explorer reconstruction helpers can select NJ.

    :return: None.
    """
    root = reconstruct_tree(
        {
            "labels": ["A", "B"],
            "matrix": [[0, 2], [2, 0]],
            "distance_metric": "hamming",
        },
        method="nj",
    )

    assert reconstructed_tree_newick(root) == "(A:1.000000,B:1.000000);"


def test_cli_writes_neighbor_joining_output_and_default_name(tmp_path):
    """Confirm the CLI selects NJ and names default output by method and metric.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    input_path = tmp_path / "distance_matrix_hamming.json"
    expected_output = tmp_path / "nj_hamming.newick"
    input_path.write_text(
        json.dumps({
            "labels": ["A", "B"],
            "matrix": [[0, 2], [2, 0]],
            "distance_metric": "hamming",
        }),
        encoding="utf-8",
    )

    assert main(["--input", str(input_path), "--method", "nj"]) == 0
    assert expected_output.read_text(encoding="utf-8") == "(A:1.000000,B:1.000000);\n"
