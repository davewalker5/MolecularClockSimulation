"""Shared download filename helpers for the interactive explorers."""

from __future__ import annotations

from collections.abc import MutableMapping
from pathlib import PurePath
from typing import Any

from common.constants import (
    DOWNLOAD_CALIBRATED_TREE_NEWICK,
    DOWNLOAD_CALIBRATED_TREE_PNG,
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
    reconstruction_algorithm: str | None = None,
) -> str:
    """Return the default file stem for a download selection.

    :param selection: User-selected download option.
    :param distance_matrix: Current calculated distance matrix, when available.
    :param reconstruction_algorithm: Algorithm used for the current reconstructed tree.
    :return: Default filename stem without a path or extension.
    """
    # Most exports have a fixed descriptive stem shared by both explorer variants.
    fixed_stems = {
        DOWNLOAD_TERMINAL_SEQUENCES: "terminal_sequences",
        DOWNLOAD_TRUE_TREE_NEWICK: "true_tree",
        DOWNLOAD_TRUE_TREE_PNG: "true_tree",
        DOWNLOAD_SIMULATION_METADATA: "simulation_metadata",
    }
    if selection in fixed_stems:
        return fixed_stems[selection]

    if selection in {DOWNLOAD_DISTANCE_MATRIX_JSON, DOWNLOAD_DISTANCE_MATRIX_CSV}:
        # The payload records the method that actually produced the current matrix.
        method = distance_matrix.get("distance_metric") if distance_matrix else None
        return f"distance_matrix_{method}" if method else "distance_matrix"

    if selection in {DOWNLOAD_CALIBRATED_TREE_NEWICK, DOWNLOAD_CALIBRATED_TREE_PNG}:
        # Calibration inherits the distance method used to reconstruct its input tree.
        method = distance_matrix.get("distance_metric") if distance_matrix else None
        suffix = "_".join(part for part in (method, reconstruction_algorithm) if part)
        return f"calibrated_tree_{suffix}" if suffix else "calibrated_tree"

    if selection in {DOWNLOAD_RECONSTRUCTED_TREE_NEWICK, DOWNLOAD_RECONSTRUCTED_TREE_PNG}:
        # A reconstructed tree is identified by both stages that produced it.
        method = distance_matrix.get("distance_metric") if distance_matrix else None
        suffix = "_".join(part for part in (method, reconstruction_algorithm) if part)
        return f"reconstructed_tree_{suffix}" if suffix else "reconstructed_tree"

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


def validate_download_stem(value: str) -> tuple[str | None, str | None]:
    """Validate a user-entered download file stem.

    :param value: Raw user-entered file stem.
    :return: Cleaned stem and error message, one of which is None.
    """
    stem = value.strip()
    if not stem:
        return None, "Enter a file stem before downloads are available."
    # Streamlit download names must be plain filenames rather than paths.
    if "/" in stem or "\\" in stem:
        return None, "Enter a file stem only, without a folder or path."
    if stem in {".", ".."}:
        return None, "Enter a file stem, not a relative path marker."
    if PurePath(stem).suffix:
        return None, "Enter the name without a file extension."
    return stem, None


def download_filename(stem: str, extension: str) -> str:
    """Build a download filename from a validated stem and extension.

    :param stem: Validated file stem without path or extension.
    :param extension: File extension without a leading dot.
    :return: Complete filename for a download.
    """
    # Keep extension assembly consistent for every explorer download.
    return f"{stem}.{extension}"


def synchronize_download_state(
    state: MutableMapping[str, Any],
    *,
    selection_key: str,
    stem_key: str,
    unavailable_key: str,
    distance_matrix_key: str,
    reconstructed_newick_key: str,
    reconstructed_dot_key: str,
    reset_stem: bool = False,
) -> str | None:
    """Synchronize download availability and the editable default stem.

    :param state: Mutable session-state mapping used by the explorer.
    :param selection_key: State key containing the selected download type.
    :param stem_key: State key containing the editable filename stem.
    :param unavailable_key: State key recording prior unavailability.
    :param distance_matrix_key: State key containing the calculated matrix.
    :param reconstructed_newick_key: State key containing reconstructed Newick text.
    :param reconstructed_dot_key: State key containing reconstructed DOT source.
    :param reset_stem: Whether a selection change must overwrite the current stem.
    :return: Warning text when the selected download is unavailable, otherwise None.
    """
    selection = state[selection_key]
    distance_matrix = state.get(distance_matrix_key)
    warning = download_unavailable_warning(
        selection,
        distance_matrix=distance_matrix,
        reconstructed_newick=state.get(reconstructed_newick_key),
        reconstructed_dot=state.get(reconstructed_dot_key),
    )
    was_unavailable = state.get(unavailable_key, False)

    # Clear unavailable downloads and restore defaults on selection or availability changes.
    if warning:
        state[stem_key] = ""
    elif reset_stem or was_unavailable:
        state[stem_key] = default_download_stem(selection, distance_matrix)
    state[unavailable_key] = warning is not None
    return warning
