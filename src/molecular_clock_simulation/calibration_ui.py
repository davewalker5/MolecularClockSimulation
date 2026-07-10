"""Shared Streamlit controls for calibrating reconstructed trees."""

from __future__ import annotations

from typing import Any

from treecalibration.calibration import calibrate_tree, to_newick
from treecomparison.comparison import parse_newick, tree_to_dot


def calibration_state_keys(state_prefix: str) -> tuple[str, ...]:
    """Return session-state keys owned by one explorer's calibration panel.

    :param state_prefix: Explorer-specific prefix such as ``strict`` or ``relaxed``.
    :return: Calibration result and error session-state keys.
    """
    # Keeping the key list central makes stale calibration results easy to clear.
    return (
        f"{state_prefix}_calibration_newick",
        f"{state_prefix}_calibration_dot",
        f"{state_prefix}_calibration_metadata",
        f"{state_prefix}_calibration_error",
    )


def clear_calibration_state(session_state: Any, state_prefix: str) -> None:
    """Remove cached calibration results after an upstream tree changes.

    :param session_state: Streamlit-compatible mutable session state.
    :param state_prefix: Explorer-specific session-state prefix.
    :return: None.
    """
    # A calibrated tree is meaningful only for the exact reconstruction it used.
    for key in calibration_state_keys(state_prefix):
        session_state.pop(key, None)
    # Reset widget selections too because a new simulation may have different taxa.
    for suffix in ("taxon_a", "taxon_b", "age_mya"):
        session_state.pop(f"{state_prefix}_calibration_{suffix}", None)


def render_calibration_tab(
    reconstructed_newick: str | None,
    state_prefix: str,
) -> None:
    """Render calibration controls and the current calibrated tree.

    :param reconstructed_newick: Reconstructed Newick text, when available.
    :param state_prefix: Explorer-specific prefix for Streamlit widget and state keys.
    :return: None.
    """
    import streamlit as st

    st.subheader("Tree Calibration")
    with st.sidebar:
        st.header("Tree Calibration")
    if reconstructed_newick is None:
        # Calibration depends on branch lengths estimated during reconstruction.
        st.info("Complete tree reconstruction before calibrating the tree.")
        return

    try:
        tree = parse_newick(reconstructed_newick)
        taxa = sorted(
            node.name for node in tree.walk() if node.is_leaf and node.name is not None
        )
    except ValueError as error:
        st.warning(f"The reconstructed tree cannot be calibrated: {error}")
        return
    if len(taxa) < 2:
        st.warning("The reconstructed tree must contain at least two terminal taxa.")
        return

    with st.sidebar:
        with st.form(f"{state_prefix}_calibration_form"):
            # Keep workflow inputs together in the sidebar like the other stages.
            taxon_a = st.selectbox(
                "First taxon",
                options=taxa,
                key=f"{state_prefix}_calibration_taxon_a",
            )
            taxon_b = st.selectbox(
                "Second taxon",
                options=taxa,
                index=1,
                key=f"{state_prefix}_calibration_taxon_b",
            )
            age_mya = st.number_input(
                "MRCA age (million years)",
                min_value=0.000001,
                value=10.0,
                step=1.0,
                format="%.6f",
                key=f"{state_prefix}_calibration_age_mya",
            )
            submitted = st.form_submit_button(
                "Calibrate Tree",
                type="primary",
                width="stretch",
            )

    if submitted:
        try:
            # Reparse on submission so each result starts from the reconstructed tree.
            calibrated, metadata = calibrate_tree(
                parse_newick(reconstructed_newick),
                [taxon_a, taxon_b],
                float(age_mya),
            )
        except ValueError as error:
            st.session_state[f"{state_prefix}_calibration_error"] = str(error)
            for key in calibration_state_keys(state_prefix)[:3]:
                st.session_state.pop(key, None)
        else:
            st.session_state.pop(f"{state_prefix}_calibration_error", None)
            st.session_state[f"{state_prefix}_calibration_newick"] = to_newick(calibrated)
            st.session_state[f"{state_prefix}_calibration_dot"] = tree_to_dot(
                calibrated,
                "Calibrated Phylogeny",
                f"branch lengths in million years; scale factor {metadata['scale_factor']:.6g}",
            )
            st.session_state[f"{state_prefix}_calibration_metadata"] = metadata

    error = st.session_state.get(f"{state_prefix}_calibration_error")
    calibrated_dot = st.session_state.get(f"{state_prefix}_calibration_dot")
    if error:
        st.warning(error)
    elif calibrated_dot:
        st.graphviz_chart(calibrated_dot, width="stretch")
        with st.expander("Calibrated Newick and metadata"):
            st.code(
                st.session_state[f"{state_prefix}_calibration_newick"],
                language="text",
            )
            st.json(st.session_state[f"{state_prefix}_calibration_metadata"])
    else:
        st.info("Choose two taxa and their MRCA age, then submit the calibration.")
