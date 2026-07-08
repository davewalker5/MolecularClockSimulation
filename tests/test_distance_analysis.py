import pytest

from distancematrix.calculator import calculate_distance_matrix, hky85_distance, kimura_distance
from molecular_clock_simulation.distance_analysis import (
    finite_or_none,
    format_distance,
    matrix_table_rows,
    model_explanation,
    pairwise_summary,
    transition_transversion_ratio,
)


def test_pairwise_summary_uses_matrix_distance_for_selected_taxa():
    """Confirm pairwise model distance is read from the calculated matrix.

    :return: None.
    """
    sequences = {
        "Species_A": "AAAA",
        "Species_B": "GAAA",
        "Species_C": "ACAA",
    }
    matrix = calculate_distance_matrix(sequences, distance_type="k80")

    summary = pairwise_summary(sequences, "Species_B", "Species_C", "k80", matrix)

    assert summary["sequence_length"] == 4
    assert summary["differing_sites"] == 2
    assert summary["hamming_distance"] == 2
    assert summary["proportional_distance"] == 0.5
    assert summary["estimated_distance"] == pytest.approx(kimura_distance("GAAA", "ACAA"))
    assert summary["transitions"] == 1
    assert summary["transversions"] == 1
    assert summary["transition_transversion_ratio"] == 1.0


def test_pairwise_summary_includes_hky85_base_frequencies():
    """Confirm HKY85 summaries expose both substitution classes and base frequencies.

    :return: None.
    """
    sequences = {
        "Species_A": "AAAACCCC",
        "Species_B": "GGAACCCC",
    }
    matrix = calculate_distance_matrix(sequences, distance_type="hky85")

    summary = pairwise_summary(sequences, "Species_A", "Species_B", "hky85", matrix)

    assert summary["estimated_distance"] == pytest.approx(hky85_distance("AAAACCCC", "GGAACCCC"))
    assert summary["transitions"] == 2
    assert summary["transversions"] == 0
    assert summary["transition_transversion_ratio"] == float("inf")
    assert summary["nucleotide_frequencies"] == {
        "A": 0.375,
        "C": 0.5,
        "G": 0.125,
        "T": 0.0,
    }


def test_matrix_table_rows_formats_taxon_label_and_distances():
    """Confirm matrix payloads become Streamlit-friendly row dictionaries.

    :return: None.
    """
    rows = matrix_table_rows({
        "labels": ["Species_A", "Species_B"],
        "matrix": [[0.0, 0.125], [0.125, float("inf")]],
        "distance_metric": "jc69",
    })

    assert rows == [
        {"Taxon": "Species_A", "Species_A": "0", "Species_B": "0.125"},
        {"Taxon": "Species_B", "Species_A": "0.125", "Species_B": "infinity"},
    ]


def test_format_distance_and_heatmap_values_handle_infinity():
    """Confirm saturated distances display clearly and do not enter colour scaling.

    :return: None.
    """
    assert format_distance(float("inf")) == "infinity"
    assert finite_or_none(float("inf")) is None
    assert finite_or_none(0.25) == 0.25


def test_transition_transversion_ratio_handles_zero_transversions():
    """Confirm transition-only and identical comparisons have useful display ratios.

    :return: None.
    """
    assert transition_transversion_ratio(2, 0) == float("inf")
    assert transition_transversion_ratio(0, 0) == 0.0
    assert transition_transversion_ratio(3, 2) == 1.5


def test_model_explanation_names_selected_model():
    """Confirm model explanations are keyed by reusable distance metric names.

    :return: None.
    """
    explanation = model_explanation("jc69")

    assert explanation.startswith("**JC69**")
    assert "multiple mutations" in explanation
