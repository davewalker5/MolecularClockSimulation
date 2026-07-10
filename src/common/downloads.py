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
