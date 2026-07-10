"""Shared reconstruction helpers for the Streamlit explorers."""

from __future__ import annotations

from typing import Any

from common.graphviz import dot_escape as _dot_escape

from phylogeny.neighbor_joining import neighbor_joining
from phylogeny.upgma import Cluster, to_newick, upgma


def reconstruct_upgma_tree(matrix_payload: dict[str, Any]) -> Cluster:
    """Reconstruct a tree from a calculated distance matrix using UPGMA.

    :param matrix_payload: Matrix payload returned by calculate_distance_matrix.
    :return: Root cluster of the reconstructed tree.
    """
    # Retain the original named helper for existing callers and tests.
    return upgma(matrix_payload["labels"], matrix_payload["matrix"])


def reconstruct_tree(matrix_payload: dict[str, Any], *, method: str) -> Cluster:
    """Reconstruct a tree from a calculated distance matrix using the selected method.

    :param matrix_payload: Matrix payload returned by calculate_distance_matrix.
    :param method: Reconstruction method identifier, either upgma or nj.
    :return: Root cluster of the reconstructed tree.
    """
    # The explorer stores user-facing labels, while the CLI uses lowercase identifiers.
    method_key = method.strip().lower().replace("neighbor joining", "nj")
    reconstructors = {
        "upgma": upgma,
        "nj": neighbor_joining,
    }
    if method_key not in reconstructors:
        raise ValueError(f"Unknown reconstruction method: {method}")
    return reconstructors[method_key](matrix_payload["labels"], matrix_payload["matrix"])


def reconstructed_tree_to_dot(
    root: Cluster,
    *,
    graph_name: str,
    colors: dict[str, str],
) -> str:
    """Render a reconstructed cluster tree as Graphviz DOT.

    :param root: Root cluster returned by a reconstruction algorithm.
    :param graph_name: DOT graph identifier.
    :param colors: Theme color mapping used by the explorers.
    :return: DOT source that can be displayed by Streamlit.
    """
    lines = [
        f"digraph {graph_name} {{",
        f"  graph [rankdir=LR, bgcolor=\"{colors['page_bg']}\", margin=0.08];",
        "  node [shape=box, style=\"rounded,filled\", "
        f"fillcolor=\"{colors['surface_elevated']}\", "
        f"color=\"{colors['border_strong']}\", "
        f"fontcolor=\"{colors['text']}\", "
        "fontname=\"Helvetica\", fontsize=11];",
        f"  edge [color=\"{colors['text_subtle']}\", "
        f"fontcolor=\"{colors['text_muted_strong']}\", "
        "fontname=\"Helvetica\", fontsize=10];",
    ]
    for node in _walk_clusters(root):
        lines.append(
            f'  "{_dot_escape(_cluster_id(node))}" '
            f'[label="{_dot_escape(_cluster_label(node))}"];'
        )
        for child, branch_length in _cluster_children(node):
            lines.append(
                f'  "{_dot_escape(_cluster_id(node))}" -> "{_dot_escape(_cluster_id(child))}" '
                f'[label="{branch_length:.4g}"];'
            )
    lines.append("}")
    return "\n".join(lines)


def reconstructed_tree_newick(root: Cluster) -> str:
    """Return Newick text for a reconstructed cluster tree.

    :param root: Root cluster returned by a reconstruction algorithm.
    :return: Newick tree ending with a semicolon.
    """
    return to_newick(root)


def _walk_clusters(root: Cluster) -> list[Cluster]:
    """Return clusters in preorder for stable graph rendering.

    :param root: Root cluster to traverse.
    :return: List of clusters in parent-before-child order.
    """
    # Start with the parent so DOT node declarations are stable and easy to inspect.
    nodes = [root]
    if root.left is not None:
        nodes.extend(_walk_clusters(root.left))
    if root.right is not None:
        nodes.extend(_walk_clusters(root.right))
    if root.third is not None:
        nodes.extend(_walk_clusters(root.third))
    return nodes


def _cluster_children(root: Cluster) -> list[tuple[Cluster, float]]:
    """Return existing child clusters with incoming branch lengths.

    :param root: Cluster whose children should be inspected.
    :return: Existing child clusters paired with branch lengths from the parent.
    """
    children: list[tuple[Cluster, float]] = []
    # Only include edges that have both a child node and a calculated branch length.
    if root.left is not None and root.left_branch_length is not None:
        children.append((root.left, root.left_branch_length))
    if root.right is not None and root.right_branch_length is not None:
        children.append((root.right, root.right_branch_length))
    if root.third is not None and root.third_branch_length is not None:
        children.append((root.third, root.third_branch_length))
    return children


def _cluster_id(cluster: Cluster) -> str:
    """Return a stable DOT node identifier for one cluster.

    :param cluster: Cluster to identify in DOT output.
    :return: Stable identifier derived from descendant labels.
    """
    # Member labels are unique and sorted by the reconstruction algorithms.
    return "+".join(cluster.members)


def _cluster_label(cluster: Cluster) -> str:
    """Return a concise display label for one reconstructed cluster.

    :param cluster: Cluster to label in the explorer graph.
    :return: Human-readable label for a leaf or internal cluster.
    """
    if cluster.left is None and cluster.right is None and cluster.third is None:
        return cluster.name
    # UPGMA heights are meaningful; NJ internal nodes keep height at zero as a placeholder.
    if cluster.height == 0.0:
        return "cluster"
    return f"cluster\nheight {cluster.height:.4g}"
