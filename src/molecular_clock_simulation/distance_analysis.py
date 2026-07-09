"""Shared Streamlit distance-analysis UI for molecular clock explorers."""

from __future__ import annotations

import math
from typing import Any

from distancematrix.calculator import (
    calculate_distance_matrix,
    calculate_nucleotide_frequencies,
    count_hky85_substitutions,
    count_k80_substitutions,
    hamming_distance,
    proportional_distance,
)

DISTANCE_MODEL_OPTIONS = {
    "Hamming Distance": "hamming",
    "Proportional Distance (p-distance)": "proportional",
    "JC69": "jc69",
    "K80": "k80",
    "F81": "f81",
    "HKY85": "hky85",
}

MODEL_LABELS = {
    model_key: model_label
    for model_label, model_key in DISTANCE_MODEL_OPTIONS.items()
}

MODEL_EXPLANATIONS = {
    "hamming": (
        "Counts the aligned sites where two sequences differ. It makes no correction "
        "for sequence length or hidden repeated mutations, so it is best understood as "
        "the raw observed difference count."
    ),
    "proportional": (
        "Extends Hamming distance by dividing the observed difference count by sequence "
        "length. This makes comparisons easier across equal-length alignments, but it "
        "still treats observed differences as the full evolutionary distance."
    ),
    "jc69": (
        "Assumes all nucleotide substitutions are equally likely and all bases have "
        "equal long-term frequencies. It extends p-distance by correcting for multiple "
        "mutations occurring at the same site."
    ),
    "k80": (
        "Extends JC69 by recognising that transitions generally occur more frequently "
        "than transversions. It corrects the two observed substitution classes separately."
    ),
    "f81": (
        "Extends JC69 by allowing unequal nucleotide frequencies. This matters when the "
        "sequence composition is biased, because the chance of observing a difference "
        "depends on the bases that are common in the alignment."
    ),
    "hky85": (
        "Combines the K80 transition/transversion distinction with the F81 unequal-base "
        "frequency correction. It is useful when both substitution bias and nucleotide "
        "composition bias affect the observed differences."
    ),
}

NO_CORRECTION_MODELS = {"hamming", "proportional"}


def render_distance_analysis_tab(
    sequences: dict[str, str],
    *,
    state_key_prefix: str,
) -> None:
    """Render an interactive distance-analysis tab for generated terminal sequences.

    :param sequences: Mapping of taxon labels to aligned terminal sequences.
    :param state_key_prefix: Prefix used to isolate Streamlit session-state keys.
    :return: None.
    """
    # Import UI dependencies lazily so tests can import helper functions without Streamlit.
    import streamlit as st

    labels = list(sequences)
    model_state_key = f"{state_key_prefix}_distance_model"
    matrix_state_key = f"{state_key_prefix}_distance_matrix"
    selected_model = selected_distance_model(state_key_prefix)

    active_model = st.session_state.get(model_state_key)
    matrix_payload = st.session_state.get(matrix_state_key)
    if matrix_payload is None or active_model is None:
        st.info(
            "Choose a distance model in the sidebar, then calculate distances "
            "to analyse this simulation."
        )
        return

    if active_model != selected_model:
        st.caption(
            f"Showing the last calculated {MODEL_LABELS[active_model]} matrix. "
            "Press Calculate Distances in the sidebar to update the analysis."
        )

    st.subheader("Pairwise Sequence Comparison")
    first_column, second_column = st.columns(2)
    with first_column:
        taxon_a = st.selectbox(
            "Taxon A",
            options=labels,
            index=0,
            key=f"{state_key_prefix}_taxon_a",
        )
    with second_column:
        taxon_b = st.selectbox(
            "Taxon B",
            options=labels,
            index=1 if len(labels) > 1 else 0,
            key=f"{state_key_prefix}_taxon_b",
        )

    summary = pairwise_summary(
        sequences,
        taxon_a,
        taxon_b,
        active_model,
        matrix_payload,
    )
    render_pairwise_summary(summary)
    render_distance_correction_card(summary)

    st.subheader("Distance Matrix")
    heatmap = build_distance_heatmap(matrix_payload, title=MODEL_LABELS[active_model])
    st.plotly_chart(heatmap, width="stretch")

    with st.expander("Numeric matrix"):
        st.dataframe(
            matrix_table_rows(matrix_payload),
            width="stretch",
            hide_index=True,
        )


def render_distance_analysis_controls(
    sequences: dict[str, str],
    *,
    state_key_prefix: str,
) -> None:
    """Render sidebar controls for distance-analysis calculations.

    :param sequences: Mapping of taxon labels to aligned terminal sequences.
    :param state_key_prefix: Prefix used to isolate Streamlit session-state keys.
    :return: None.
    """
    import streamlit as st

    model_state_key = f"{state_key_prefix}_distance_model"
    matrix_state_key = f"{state_key_prefix}_distance_matrix"
    selected_label = st.selectbox(
        "Distance model",
        options=list(DISTANCE_MODEL_OPTIONS),
        key=f"{state_key_prefix}_distance_model_selection",
    )
    selected_model = DISTANCE_MODEL_OPTIONS[selected_label]

    if st.button("Calculate Distances", key=f"{state_key_prefix}_calculate_distances"):
        # Store the model together with the matrix so model changes do not recalculate implicitly.
        st.session_state[model_state_key] = selected_model
        st.session_state[matrix_state_key] = calculate_distance_matrix(
            sequences,
            distance_type=selected_model,
        )

    st.subheader("Model Explanation")
    st.markdown(model_explanation(selected_model))


def selected_distance_model(state_key_prefix: str) -> str:
    """Return the currently selected distance model for a Streamlit session.

    :param state_key_prefix: Prefix used to isolate Streamlit session-state keys.
    :return: Selected distance model key.
    """
    import streamlit as st

    selected_label = st.session_state.get(
        f"{state_key_prefix}_distance_model_selection",
        next(iter(DISTANCE_MODEL_OPTIONS)),
    )
    return DISTANCE_MODEL_OPTIONS[selected_label]


def pairwise_summary(
    sequences: dict[str, str],
    taxon_a: str,
    taxon_b: str,
    distance_type: str,
    matrix_payload: dict[str, Any],
) -> dict[str, Any]:
    """Summarize observed and model-estimated differences for two selected taxa.

    :param sequences: Mapping of taxon labels to aligned terminal sequences.
    :param taxon_a: Label for the first selected taxon.
    :param taxon_b: Label for the second selected taxon.
    :param distance_type: Calculated distance model key.
    :param matrix_payload: Matrix payload returned by calculate_distance_matrix.
    :return: Dictionary of display-ready pairwise comparison values.
    """
    first = sequences[taxon_a]
    second = sequences[taxon_b]
    labels = matrix_payload["labels"]
    row_index = labels.index(taxon_a)
    column_index = labels.index(taxon_b)

    # Pull the model estimate from the full matrix so pairwise and matrix values match exactly.
    estimated_distance = matrix_payload["matrix"][row_index][column_index]
    differing_sites = hamming_distance(first, second)
    p_distance = proportional_distance(first, second)
    summary: dict[str, Any] = {
        "model": distance_type,
        "sequence_length": len(first),
        "differing_sites": differing_sites,
        "hamming_distance": differing_sites,
        "proportional_distance": p_distance,
        "estimated_distance": estimated_distance,
    }

    if distance_type in {"k80", "hky85"}:
        # K80 and HKY85 expose transition/transversion counts as their key observation.
        counter = count_hky85_substitutions if distance_type == "hky85" else count_k80_substitutions
        transitions, transversions = counter(first, second)
        summary.update({
            "transitions": transitions,
            "transversions": transversions,
            "transition_transversion_ratio": transition_transversion_ratio(
                transitions,
                transversions,
            ),
        })

    if distance_type in {"f81", "hky85"}:
        # F81 and HKY85 estimate base composition from the selected pair.
        summary["nucleotide_frequencies"] = calculate_nucleotide_frequencies(first, second)

    return summary


def transition_transversion_ratio(transitions: int, transversions: int) -> float:
    """Calculate the observed transition/transversion ratio.

    :param transitions: Number of observed transition differences.
    :param transversions: Number of observed transversion differences.
    :return: Transition count divided by transversion count, or infinity when no transversions exist.
    """
    # A transition-only comparison has an undefined finite ratio, represented as infinity.
    if transversions == 0:
        return float("inf") if transitions else 0.0
    return transitions / transversions


def render_pairwise_summary(summary: dict[str, Any]) -> None:
    """Render pairwise comparison metrics and relevant intermediate statistics.

    :param summary: Pairwise summary returned by pairwise_summary.
    :return: None.
    """
    import streamlit as st

    metric_columns = st.columns(5)
    metric_columns[0].metric("Sequence length", summary["sequence_length"])
    metric_columns[1].metric("Differing sites", summary["differing_sites"])
    metric_columns[2].metric("Hamming distance", summary["hamming_distance"])
    metric_columns[3].metric("p-distance", format_distance(summary["proportional_distance"]))
    metric_columns[4].metric("Model distance", format_distance(summary["estimated_distance"]))

    if summary["model"] in {"k80", "hky85"}:
        substitution_columns = st.columns(3)
        substitution_columns[0].metric("Transitions", summary["transitions"])
        substitution_columns[1].metric("Transversions", summary["transversions"])
        substitution_columns[2].metric(
            "Ti/Tv ratio",
            format_distance(summary["transition_transversion_ratio"]),
        )

    if "nucleotide_frequencies" in summary:
        frequency_text = " | ".join(
            f"{base}: {frequency:.3f}"
            for base, frequency in summary["nucleotide_frequencies"].items()
        )
        st.caption(f"Observed nucleotide frequencies: {frequency_text}")


def correction_card_values(summary: dict[str, Any]) -> dict[str, int | float | str]:
    """Calculate display values for the pairwise distance-correction card.

    :param summary: Pairwise summary returned by pairwise_summary.
    :return: Dictionary containing observed, corrected, and correction display values.
    """
    model = summary["model"]
    sequence_length = summary["sequence_length"]
    observed_distance = summary["proportional_distance"]

    # Hamming is a count, so p-distance is the comparable per-site value for this card.
    corrected_distance = observed_distance if model == "hamming" else summary["estimated_distance"]
    correction_amount = (
        float("inf")
        if is_infinite_distance(corrected_distance)
        else corrected_distance - observed_distance
    )
    correction_percent = (
        float("inf")
        if is_infinite_distance(correction_amount)
        else correction_amount * 100
    )
    hidden_substitutions = (
        float("inf")
        if is_infinite_distance(correction_amount)
        else correction_amount * sequence_length
    )

    if model in NO_CORRECTION_MODELS:
        interpretation = "No model correction is applied; observed differences are shown directly."
    elif is_infinite_distance(corrected_distance):
        interpretation = (
            "The selected model treats this comparison as saturated, so a finite "
            "corrected distance cannot be estimated."
        )
    elif correction_amount == 0:
        interpretation = "The correction is zero because no hidden substitutions are inferred."
    else:
        interpretation = (
            "The corrected estimate is larger than the observed p-distance because "
            "the model accounts for substitutions hidden by repeated changes at the same site."
        )

    return {
        "observed_distance": observed_distance,
        "corrected_distance": corrected_distance,
        "correction_amount": correction_amount,
        "correction_percent": correction_percent,
        "hidden_substitutions": hidden_substitutions,
        "interpretation": interpretation,
    }


def render_distance_correction_card(summary: dict[str, Any]) -> None:
    """Render a compact card explaining the model correction for one pairwise comparison.

    :param summary: Pairwise summary returned by pairwise_summary.
    :return: None.
    """
    import streamlit as st

    values = correction_card_values(summary)

    with st.container(border=True):
        st.markdown("#### Distance Correction")
        correction_columns = st.columns(4)
        correction_columns[0].metric(
            "Observed p-distance",
            format_distance(values["observed_distance"]),
        )
        correction_columns[1].metric(
            "Corrected distance",
            format_distance(values["corrected_distance"]),
        )
        correction_columns[2].metric(
            "Distance correction",
            format_signed_percentage(values["correction_percent"]),
        )
        correction_columns[3].metric(
            "Estimated hidden substitutions",
            format_distance(values["hidden_substitutions"]),
        )
        st.caption(values["interpretation"])


def build_distance_heatmap(matrix_payload: dict[str, Any], *, title: str) -> Any:
    """Build a Plotly heatmap for a calculated distance matrix.

    :param matrix_payload: Matrix payload returned by calculate_distance_matrix.
    :param title: Human-readable model title used in the chart heading.
    :return: Plotly Figure containing the distance heatmap.
    """
    import plotly.graph_objects as go

    labels = matrix_payload["labels"]
    matrix = matrix_payload["matrix"]
    text_matrix = [
        [format_distance(value) for value in row]
        for row in matrix
    ]
    z_matrix = [
        [finite_or_none(value) for value in row]
        for row in matrix
    ]

    # Plotly scales the colour range from the finite z values in the current matrix.
    figure = go.Figure(
        data=go.Heatmap(
            z=z_matrix,
            x=labels,
            y=labels,
            text=text_matrix,
            hovertemplate="Taxon A: %{y}<br>Taxon B: %{x}<br>Distance: %{text}<extra></extra>",
            colorscale="YlOrBr",
            colorbar={"title": "Distance"},
        )
    )
    figure.update_layout(
        title=f"{title} distance matrix",
        xaxis_title="Taxon B",
        yaxis_title="Taxon A",
        margin={"l": 40, "r": 20, "t": 50, "b": 40},
    )
    return figure


def matrix_table_rows(matrix_payload: dict[str, Any]) -> list[dict[str, str]]:
    """Convert a distance matrix payload into rows for Streamlit dataframe display.

    :param matrix_payload: Matrix payload returned by calculate_distance_matrix.
    :return: List of row dictionaries with formatted distance values.
    """
    labels = matrix_payload["labels"]
    rows = []
    for label, values in zip(labels, matrix_payload["matrix"], strict=True):
        # Include a Taxon column so the row label remains visible in Streamlit's dataframe.
        row = {"Taxon": label}
        row.update({
            column_label: format_distance(value)
            for column_label, value in zip(labels, values, strict=True)
        })
        rows.append(row)
    return rows


def model_explanation(distance_type: str) -> str:
    """Return an educational explanation for one distance model.

    :param distance_type: Distance model key.
    :return: Markdown text describing the selected model.
    """
    label = MODEL_LABELS[distance_type]
    return f"**{label}**\n\n{MODEL_EXPLANATIONS[distance_type]}"


def format_distance(value: int | float) -> str:
    """Format a numeric distance for compact display.

    :param value: Distance value to display.
    :return: Human-readable distance string.
    """
    # Saturated model estimates are represented internally as infinity.
    if isinstance(value, float) and math.isinf(value):
        return "infinity"
    if isinstance(value, int):
        return str(value)
    return f"{value:.6g}"


def format_signed_percentage(value: int | float) -> str:
    """Format a correction percentage with an explicit sign.

    :param value: Percentage value to display.
    :return: Human-readable signed percentage string.
    """
    # Infinite corrections are possible when the selected model reports saturation.
    if is_infinite_distance(value):
        return "+infinity"
    if value == 0:
        return "0%"
    return f"{value:+.6g}%"


def is_infinite_distance(value: int | float) -> bool:
    """Return whether a distance value is positive or negative infinity.

    :param value: Numeric distance value to inspect.
    :return: True when the value is infinite, otherwise False.
    """
    # Keep the infinity check reusable across display helpers and card calculations.
    return isinstance(value, float) and math.isinf(value)


def finite_or_none(value: int | float) -> int | float | None:
    """Return finite numeric values and replace infinities with None.

    :param value: Distance value destined for Plotly colour scaling.
    :return: Original finite value, or None when the value is infinite.
    """
    # Plotly cannot colour an infinite value, but hover text still shows the exact estimate.
    if is_infinite_distance(value):
        return None
    return value
