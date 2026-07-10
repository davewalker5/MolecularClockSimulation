"""Shared download filename helpers for the interactive explorers."""

from __future__ import annotations

from typing import Any

from common.constants import (
    DOWNLOAD_DISTANCE_MATRIX_CSV,
    DOWNLOAD_DISTANCE_MATRIX_JSON,
    DOWNLOAD_RECONSTRUCTED_TREE_NEWICK,
    DOWNLOAD_RECONSTRUCTED_TREE_PNG,
    DOWNLOAD_SIMULATION_METADATA,
    DOWNLOAD_TERMINAL_SEQUENCES,
    DOWNLOAD_TRUE_TREE_NEWICK,
    DOWNLOAD_TRUE_TREE_PNG,
)


def default_download_stem(
    selection: str,
    distance_matrix: dict[str, Any] | None = None,
) -> str:
    """Return the default file stem for a download selection.

    :param selection: User-selected download option.
    :param distance_matrix: Current calculated distance matrix, when available.
    :return: Default filename stem without a path or extension.
    """
    # Most exports have a fixed descriptive stem shared by both explorer variants.
    fixed_stems = {
        DOWNLOAD_TERMINAL_SEQUENCES: "terminal_sequences",
        DOWNLOAD_TRUE_TREE_NEWICK: "true_tree",
        DOWNLOAD_TRUE_TREE_PNG: "true_tree",
        DOWNLOAD_SIMULATION_METADATA: "simulation_metadata",
        DOWNLOAD_RECONSTRUCTED_TREE_NEWICK: "reconstructed_tree",
        DOWNLOAD_RECONSTRUCTED_TREE_PNG: "reconstructed_tree",
    }
    if selection in fixed_stems:
        return fixed_stems[selection]

    if selection in {DOWNLOAD_DISTANCE_MATRIX_JSON, DOWNLOAD_DISTANCE_MATRIX_CSV}:
        # The payload records the method that actually produced the current matrix.
        method = distance_matrix.get("distance_metric") if distance_matrix else None
        return f"distance_matrix_{method}" if method else "distance_matrix"

    raise ValueError(f"Unknown download selection: {selection}")


def download_unavailable_warning(
    selection: str,
    *,
    distance_matrix: dict[str, Any] | None = None,
    reconstructed_newick: str | None = None,
    reconstructed_dot: str | None = None,
) -> str | None:
    """Return a warning when the selected download data is not available.

    :param selection: User-selected download option.
    :param distance_matrix: Current calculated distance matrix, when available.
    :param reconstructed_newick: Current reconstructed Newick text, when available.
    :param reconstructed_dot: Current reconstructed tree DOT source, when available.
    :return: Warning text when unavailable, otherwise None.
    """
    # Matrix downloads are unavailable until a distance calculation has completed.
    if (
        selection in {DOWNLOAD_DISTANCE_MATRIX_JSON, DOWNLOAD_DISTANCE_MATRIX_CSV}
        and distance_matrix is None
    ):
        return "You must calculate a distance matrix before downloading it."

    # Check the representation required by the selected reconstruction format.
    if selection == DOWNLOAD_RECONSTRUCTED_TREE_NEWICK and reconstructed_newick is None:
        return "You must reconstruct a tree before downloading it."
    if selection == DOWNLOAD_RECONSTRUCTED_TREE_PNG and reconstructed_dot is None:
        return "You must reconstruct a tree before downloading it."
    return None
