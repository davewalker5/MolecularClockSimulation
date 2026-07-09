import csv
import json
import math

import pytest

from distancematrix import (
    calculate_nucleotide_frequencies,
    calculate_distance_matrix,
    f81_distance,
    hky85_distance,
    jc69_distance,
    kimura_distance,
    read_fasta,
    write_csv,
    write_json,
)
from distancematrix.calculator import count_hky85_substitutions, count_k80_substitutions, hamming_distance, write_outputs
from distancematrix.calculator import proportional_distance
from distancematrix.cli import main


def test_read_fasta_parses_multiline_sequences_and_preserves_labels(tmp_path):
    """Confirm aligned FASTA input is parsed into reusable sequence data.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">Species_A\nACGT\nACGT\n>Species_B\nACGT\nTCGT\n",
        encoding="utf-8",
    )

    assert read_fasta(fasta) == {
        "Species_A": "ACGTACGT",
        "Species_B": "ACGTTCGT",
    }


def test_calculate_distance_matrix_returns_hamming_counts():
    """Confirm pairwise Hamming counts are returned in a square matrix.

    :return: None.
    """
    matrix = calculate_distance_matrix({
        "Species_A": "ACGTACGT",
        "Species_B": "ACGTTCGT",
        "Species_C": "ACGGACGA",
    })

    assert matrix == {
        "labels": ["Species_A", "Species_B", "Species_C"],
        "matrix": [
            [0, 1, 2],
            [1, 0, 3],
            [2, 3, 0],
        ],
        "distance_metric": "hamming",
    }


def test_calculate_distance_matrix_returns_proportional_distances():
    """Confirm proportional distances divide Hamming counts by sequence length.

    :return: None.
    """
    matrix = calculate_distance_matrix(
        {
            "Species_A": "ACGTACGT",
            "Species_B": "ACGTTCGT",
            "Species_C": "ACGGACGA",
        },
        distance_type="proportional",
    )

    assert matrix == {
        "labels": ["Species_A", "Species_B", "Species_C"],
        "matrix": [
            [0.0, 0.125, 0.25],
            [0.125, 0.0, 0.375],
            [0.25, 0.375, 0.0],
        ],
        "distance_metric": "proportional",
    }


def test_jc69_distance_returns_zero_for_identical_sequences():
    """Confirm identical DNA sequences have no corrected evolutionary distance.

    :return: None.
    """
    assert jc69_distance("ACGTACGT", "ACGTACGT") == 0.0


def test_jc69_distance_matches_known_example():
    """Confirm JC69 correction matches the expected formula for a known p-distance.

    :return: None.
    """
    expected = -0.75 * math.log(1 - (4 * 0.25 / 3))

    assert jc69_distance("ACGTACGT", "ACGGACGA") == pytest.approx(expected)


def test_jc69_distance_is_finite_below_saturation_limit():
    """Confirm JC69 can calculate a finite distance just below saturation.

    :return: None.
    """
    expected = -0.75 * math.log(1 - (4 * 0.7 / 3))

    assert jc69_distance("AAAAAAAAAA", "CCCCCCCAAA") == pytest.approx(expected)


@pytest.mark.parametrize(
    "first,second",
    [
        ("AAAAAAAA", "CCCCCCAA"),
        ("AAAAAAAA", "CCCCCCCC"),
    ],
)
def test_jc69_distance_returns_infinity_at_and_above_saturation(first, second):
    """Confirm saturated JC69 comparisons return infinity instead of raising math errors.

    :param first: First aligned DNA sequence.
    :param second: Second aligned DNA sequence.
    :return: None.
    """
    assert jc69_distance(first, second) == float("inf")


def test_jc69_distance_rejects_invalid_dna_symbols():
    """Confirm JC69 validates DNA symbols before applying the substitution model.

    :return: None.
    """
    with pytest.raises(ValueError, match="invalid DNA bases"):
        jc69_distance("ACGTNN", "ACGTAA")


def test_kimura_distance_returns_zero_for_identical_sequences():
    """Confirm identical DNA sequences have no K80 corrected evolutionary distance.

    :return: None.
    """
    assert kimura_distance("AAAA", "AAAA") == 0.0


def test_k80_substitution_counter_classifies_transition_only_differences():
    """Confirm K80 counts transition substitutions independently of transversions.

    :return: None.
    """
    assert count_k80_substitutions("AAGC", "GGGC") == (2, 0)


def test_kimura_distance_matches_transition_only_example():
    """Confirm K80 correction matches the expected formula for transitions only.

    :return: None.
    """
    assert kimura_distance("AAGC", "GGGC") == float("inf")

    expected = -0.5 * math.log(1 - (2 * 0.25) - 0) - 0.25 * math.log(1 - 0)
    assert kimura_distance("AACC", "GACC") == pytest.approx(expected)


def test_k80_substitution_counter_classifies_transversion_only_differences():
    """Confirm K80 counts transversion substitutions independently of transitions.

    :return: None.
    """
    assert count_k80_substitutions("AAAA", "CCCC") == (0, 4)


def test_kimura_distance_matches_transversion_only_example():
    """Confirm K80 correction matches the expected formula for transversions only.

    :return: None.
    """
    expected = -0.5 * math.log(1 - 0 - 0.25) - 0.25 * math.log(1 - (2 * 0.25))

    assert kimura_distance("AAAA", "ACAA") == pytest.approx(expected)


def test_kimura_distance_matches_mixed_substitution_example():
    """Confirm K80 correction handles transition and transversion substitutions together.

    :return: None.
    """
    expected = -0.5 * math.log(1 - (2 * 0.25) - 0.25) - 0.25 * math.log(1 - (2 * 0.25))

    assert count_k80_substitutions("AAAA", "GCAA") == (1, 1)
    assert kimura_distance("AAAA", "GCAA") == pytest.approx(expected)


def test_kimura_distance_rejects_unequal_lengths():
    """Confirm K80 requires aligned sequences with equal lengths.

    :return: None.
    """
    with pytest.raises(ValueError, match="same length"):
        kimura_distance("ACGT", "ACG")


def test_kimura_distance_rejects_invalid_dna_symbols():
    """Confirm K80 validates DNA symbols before applying the substitution model.

    :return: None.
    """
    with pytest.raises(ValueError, match="invalid DNA bases"):
        kimura_distance("ACGTNN", "ACGTAA")


@pytest.mark.parametrize(
    "first,second",
    [
        ("AAAA", "GGGG"),
        ("AAAA", "CCCC"),
    ],
)
def test_kimura_distance_returns_infinity_when_logarithm_terms_are_undefined(first, second):
    """Confirm saturated K80 comparisons return infinity instead of raising math errors.

    :param first: First aligned DNA sequence.
    :param second: Second aligned DNA sequence.
    :return: None.
    """
    assert kimura_distance(first, second) == float("inf")


def test_calculate_nucleotide_frequencies_pools_pairwise_sequences():
    """Confirm F81 empirical frequencies are estimated from both sequences.

    :return: None.
    """
    assert calculate_nucleotide_frequencies("AAAA", "ACGT") == {
        "A": 0.625,
        "C": 0.125,
        "G": 0.125,
        "T": 0.125,
    }


def test_f81_distance_returns_zero_for_identical_sequences():
    """Confirm identical DNA sequences have no F81 corrected evolutionary distance.

    :return: None.
    """
    assert f81_distance("AAAA", "AAAA") == 0.0


def test_f81_distance_matches_unequal_frequency_example():
    """Confirm F81 correction uses empirical nucleotide frequencies.

    :return: None.
    """
    frequencies = calculate_nucleotide_frequencies("AAAACCCC", "AAAAGGCC")
    frequency_factor = 1 - sum(frequency ** 2 for frequency in frequencies.values())
    observed_distance = 0.25
    expected = -frequency_factor * math.log(1 - (observed_distance / frequency_factor))

    assert f81_distance("AAAACCCC", "AAAAGGCC") == pytest.approx(expected)


def test_f81_distance_matches_jc69_when_frequencies_are_equal():
    """Confirm F81 collapses to JC69 when empirical base frequencies are equal.

    :return: None.
    """
    first = "ACGTACGT"
    second = "TCGTACGA"

    assert calculate_nucleotide_frequencies(first, second) == {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    }
    assert f81_distance(first, second) == pytest.approx(jc69_distance(first, second))


def test_f81_distance_rejects_invalid_dna_symbols():
    """Confirm F81 validates DNA symbols before applying the substitution model.

    :return: None.
    """
    with pytest.raises(ValueError, match="invalid DNA bases"):
        f81_distance("ACGTNN", "ACGTAA")


def test_f81_distance_returns_infinity_when_logarithm_term_is_undefined():
    """Confirm saturated F81 comparisons return infinity instead of raising math errors.

    :return: None.
    """
    assert f81_distance("AAAA", "CCCC") == float("inf")


def test_hky85_distance_returns_zero_for_identical_sequences():
    """Confirm identical DNA sequences have no HKY85 corrected evolutionary distance.

    :return: None.
    """
    assert hky85_distance("ACGTACGT", "ACGTACGT") == 0.0


def test_hky85_substitution_counter_classifies_transition_and_transversion_differences():
    """Confirm HKY85 observes transitions and transversions separately.

    :return: None.
    """
    assert count_hky85_substitutions("AACC", "GATG") == (2, 1)


def test_hky85_distance_matches_transition_only_example():
    """Confirm HKY85 correction combines transition bias with observed base frequencies.

    :return: None.
    """
    frequencies = calculate_nucleotide_frequencies("AAAACCCC", "GGAACCCC")
    frequency_factor = 2 * (frequencies["A"] + frequencies["G"]) * (frequencies["C"] + frequencies["T"])
    transition_proportion = 0.25
    transversion_proportion = 0.0
    expected = (
        -frequency_factor
        * math.log(1 - (transition_proportion / frequency_factor) - transversion_proportion)
        - 0.5 * (1 - frequency_factor) * math.log(1 - (2 * transversion_proportion))
    )

    assert hky85_distance("AAAACCCC", "GGAACCCC") == pytest.approx(expected)
    assert hky85_distance("AAAACCCC", "GGAACCCC") != pytest.approx(
        f81_distance("AAAACCCC", "GGAACCCC")
    )


def test_hky85_distance_matches_k80_when_purine_pyrimidine_frequencies_are_balanced():
    """Confirm HKY85 reduces to K80 when the empirical frequency factor is balanced.

    :return: None.
    """
    first = "AAAACCCC"
    second = "GAAATCCC"

    assert calculate_nucleotide_frequencies(first, second) == {
        "A": 0.4375,
        "C": 0.4375,
        "G": 0.0625,
        "T": 0.0625,
    }
    assert hky85_distance(first, second) == pytest.approx(kimura_distance(first, second))


@pytest.mark.parametrize(
    "first,second",
    [
        ("AAAA", "GGGG"),
        ("AAAA", "CCCC"),
    ],
)
def test_hky85_distance_returns_infinity_when_logarithm_terms_are_undefined(first, second):
    """Confirm saturated HKY85 comparisons return infinity instead of raising math errors.

    :param first: First aligned DNA sequence.
    :param second: Second aligned DNA sequence.
    :return: None.
    """
    assert hky85_distance(first, second) == float("inf")


def test_hky85_distance_rejects_invalid_dna_symbols():
    """Confirm HKY85 validates DNA symbols before applying the substitution model.

    :return: None.
    """
    with pytest.raises(ValueError, match="invalid DNA bases"):
        hky85_distance("ACGTNN", "ACGTAA")


def test_calculate_distance_matrix_returns_jc69_distances():
    """Confirm JC69 distances are returned in a symmetric square matrix.

    :return: None.
    """
    matrix = calculate_distance_matrix(
        {
            "Species_A": "ACGTACGT",
            "Species_B": "ACGTTCGT",
            "Species_C": "ACGGACGA",
        },
        distance_type="jc69",
    )

    assert matrix["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert matrix["distance_metric"] == "jc69"
    assert matrix["matrix"][0][0] == 0.0
    assert matrix["matrix"][1][1] == 0.0
    assert matrix["matrix"][2][2] == 0.0
    assert matrix["matrix"][0][1] == matrix["matrix"][1][0]
    assert matrix["matrix"][0][2] == matrix["matrix"][2][0]
    assert matrix["matrix"][1][2] == matrix["matrix"][2][1]
    assert matrix["matrix"][0][1] == pytest.approx(-0.75 * math.log(1 - (4 * 0.125 / 3)))
    assert matrix["matrix"][0][2] == pytest.approx(-0.75 * math.log(1 - (4 * 0.25 / 3)))
    assert matrix["matrix"][1][2] == pytest.approx(-0.75 * math.log(1 - (4 * 0.375 / 3)))


def test_calculate_distance_matrix_returns_k80_distances():
    """Confirm K80 distances are returned in a symmetric square matrix.

    :return: None.
    """
    matrix = calculate_distance_matrix(
        {
            "Species_A": "AAAA",
            "Species_B": "GAAA",
            "Species_C": "ACAA",
        },
        distance_type="k80",
    )

    expected_transition = kimura_distance("AAAA", "GAAA")
    expected_transversion = kimura_distance("AAAA", "ACAA")
    expected_mixed = kimura_distance("GAAA", "ACAA")

    assert matrix["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert matrix["distance_metric"] == "k80"
    assert matrix["matrix"][0][0] == 0.0
    assert matrix["matrix"][1][1] == 0.0
    assert matrix["matrix"][2][2] == 0.0
    assert matrix["matrix"][0][1] == pytest.approx(expected_transition)
    assert matrix["matrix"][1][0] == pytest.approx(expected_transition)
    assert matrix["matrix"][0][2] == pytest.approx(expected_transversion)
    assert matrix["matrix"][2][0] == pytest.approx(expected_transversion)
    assert matrix["matrix"][1][2] == pytest.approx(expected_mixed)
    assert matrix["matrix"][2][1] == pytest.approx(expected_mixed)


def test_calculate_distance_matrix_returns_f81_distances():
    """Confirm F81 distances are returned in a symmetric square matrix.

    :return: None.
    """
    matrix = calculate_distance_matrix(
        {
            "Species_A": "AAAACCCC",
            "Species_B": "AAAAGGCC",
            "Species_C": "AAAATTCC",
        },
        distance_type="f81",
    )

    expected_ab = f81_distance("AAAACCCC", "AAAAGGCC")
    expected_ac = f81_distance("AAAACCCC", "AAAATTCC")
    expected_bc = f81_distance("AAAAGGCC", "AAAATTCC")

    assert matrix["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert matrix["distance_metric"] == "f81"
    assert matrix["matrix"][0][0] == 0.0
    assert matrix["matrix"][1][1] == 0.0
    assert matrix["matrix"][2][2] == 0.0
    assert matrix["matrix"][0][1] == pytest.approx(expected_ab)
    assert matrix["matrix"][1][0] == pytest.approx(expected_ab)
    assert matrix["matrix"][0][2] == pytest.approx(expected_ac)
    assert matrix["matrix"][2][0] == pytest.approx(expected_ac)
    assert matrix["matrix"][1][2] == pytest.approx(expected_bc)
    assert matrix["matrix"][2][1] == pytest.approx(expected_bc)


def test_calculate_distance_matrix_returns_hky85_distances():
    """Confirm HKY85 distances are returned in a symmetric square matrix.

    :return: None.
    """
    matrix = calculate_distance_matrix(
        {
            "Species_A": "AAAACCCC",
            "Species_B": "GGAACCCC",
            "Species_C": "AAAATCCC",
        },
        distance_type="hky85",
    )

    expected_ab = hky85_distance("AAAACCCC", "GGAACCCC")
    expected_ac = hky85_distance("AAAACCCC", "AAAATCCC")
    expected_bc = hky85_distance("GGAACCCC", "AAAATCCC")

    assert matrix["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert matrix["distance_metric"] == "hky85"
    assert matrix["matrix"][0][0] == 0.0
    assert matrix["matrix"][1][1] == 0.0
    assert matrix["matrix"][2][2] == 0.0
    assert matrix["matrix"][0][1] == pytest.approx(expected_ab)
    assert matrix["matrix"][1][0] == pytest.approx(expected_ab)
    assert matrix["matrix"][0][2] == pytest.approx(expected_ac)
    assert matrix["matrix"][2][0] == pytest.approx(expected_ac)
    assert matrix["matrix"][1][2] == pytest.approx(expected_bc)
    assert matrix["matrix"][2][1] == pytest.approx(expected_bc)


def test_hamming_distance_rejects_unequal_lengths():
    """Confirm direct distance calculation requires aligned sequences.

    :return: None.
    """
    with pytest.raises(ValueError, match="same length"):
        hamming_distance("ACGT", "ACG")


def test_proportional_distance_rejects_unequal_lengths():
    """Confirm proportional distance also requires aligned sequences.

    :return: None.
    """
    with pytest.raises(ValueError, match="same length"):
        proportional_distance("ACGT", "ACG")


def test_calculate_distance_matrix_rejects_unknown_distance_type():
    """Confirm callers must request one of the supported distance metrics.

    :return: None.
    """
    with pytest.raises(ValueError, match="distance_type"):
        calculate_distance_matrix({"Species_A": "ACGT", "Species_B": "ACGA"}, "jukes-cantor")


def test_read_fasta_rejects_missing_file(tmp_path):
    """Confirm missing input files fail with a clear validation error.

    :return: None.
    """
    with pytest.raises(ValueError, match="does not exist"):
        read_fasta(tmp_path / "missing.fasta")


def test_read_fasta_rejects_duplicate_labels(tmp_path):
    """Confirm FASTA labels must be unique.

    :return: None.
    """
    fasta = tmp_path / "duplicate.fasta"
    fasta.write_text(">Species_A\nACGT\n>Species_A\nACGA\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unique"):
        read_fasta(fasta)


def test_read_fasta_rejects_unequal_lengths(tmp_path):
    """Confirm FASTA records must already be aligned.

    :return: None.
    """
    fasta = tmp_path / "unaligned.fasta"
    fasta.write_text(">Species_A\nACGT\n>Species_B\nACGTA\n", encoding="utf-8")

    with pytest.raises(ValueError, match="same length"):
        read_fasta(fasta)


def test_read_fasta_rejects_single_sequence(tmp_path):
    """Confirm at least two sequences are required for pairwise distances.

    :return: None.
    """
    fasta = tmp_path / "single.fasta"
    fasta.write_text(">Species_A\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="At least two"):
        read_fasta(fasta)


def test_write_json_and_csv_create_expected_outputs(tmp_path):
    """Confirm matrix writers create Python-friendly and spreadsheet-friendly files.

    :return: None.
    """
    matrix = calculate_distance_matrix({
        "Species_A": "ACGT",
        "Species_B": "ACGA",
    })
    json_path = tmp_path / "distance_matrix.json"
    csv_path = tmp_path / "distance_matrix.csv"

    write_json(matrix, json_path)
    write_csv(matrix, csv_path)

    assert json.loads(json_path.read_text(encoding="utf-8")) == matrix
    with csv_path.open(encoding="utf-8", newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["", "Species_A", "Species_B"],
            ["Species_A", "0", "1"],
            ["Species_B", "1", "0"],
        ]


def test_write_outputs_creates_both_matrix_files(tmp_path):
    """Confirm the output helper includes the distance type in both filenames.

    :return: None.
    """
    matrix = calculate_distance_matrix({
        "Species_A": "ACGT",
        "Species_B": "ACGA",
    })
    written = write_outputs(matrix, tmp_path / "output")

    assert written["json"].name == "distance_matrix_hamming.json"
    assert written["csv"].name == "distance_matrix_hamming.csv"
    assert written["json"].exists()
    assert written["csv"].exists()


def test_cli_writes_distance_matrix_files(tmp_path):
    """Confirm the command-line wrapper runs the full calculator workflow.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nACGTACGT\n>Species_B\nACGTTCGT\n>Species_C\nACGGACGA\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output)]) == 0

    assert json.loads((output / "distance_matrix_hamming.json").read_text(encoding="utf-8"))["matrix"] == [
        [0, 1, 2],
        [1, 0, 3],
        [2, 3, 0],
    ]
    assert (output / "distance_matrix_hamming.csv").exists()


def test_cli_writes_proportional_distance_matrix_files(tmp_path):
    """Confirm the CLI can write proportional matrix outputs.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nACGTACGT\n>Species_B\nACGTTCGT\n>Species_C\nACGGACGA\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output), "-dt", "proportional"]) == 0

    payload = json.loads((output / "distance_matrix_proportional.json").read_text(encoding="utf-8"))
    assert payload["distance_metric"] == "proportional"
    assert payload["matrix"] == [
        [0.0, 0.125, 0.25],
        [0.125, 0.0, 0.375],
        [0.25, 0.375, 0.0],
    ]
    with (output / "distance_matrix_proportional.csv").open(encoding="utf-8", newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["", "Species_A", "Species_B", "Species_C"],
            ["Species_A", "0.0", "0.125", "0.25"],
            ["Species_B", "0.125", "0.0", "0.375"],
            ["Species_C", "0.25", "0.375", "0.0"],
        ]


def test_cli_writes_jc69_distance_matrix_files(tmp_path):
    """Confirm the CLI can write JC69 matrix outputs without changing file formats.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nACGTACGT\n>Species_B\nACGTTCGT\n>Species_C\nACGGACGA\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output), "-dt", "jc69"]) == 0

    payload = json.loads((output / "distance_matrix_jc69.json").read_text(encoding="utf-8"))
    assert payload["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert payload["distance_metric"] == "jc69"
    assert payload["matrix"][0][0] == 0.0
    assert payload["matrix"][0][1] == pytest.approx(jc69_distance("ACGTACGT", "ACGTTCGT"))
    assert payload["matrix"][1][0] == pytest.approx(jc69_distance("ACGTTCGT", "ACGTACGT"))
    with (output / "distance_matrix_jc69.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["", "Species_A", "Species_B", "Species_C"]
    assert rows[1][0] == "Species_A"
    assert float(rows[1][1]) == 0.0
    assert float(rows[1][2]) == pytest.approx(jc69_distance("ACGTACGT", "ACGTTCGT"))


def test_cli_writes_k80_distance_matrix_files(tmp_path):
    """Confirm the CLI can write K80 matrix outputs without changing file formats.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nAAAA\n>Species_B\nGAAA\n>Species_C\nACAA\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output), "-dt", "k80"]) == 0

    payload = json.loads((output / "distance_matrix_k80.json").read_text(encoding="utf-8"))
    assert payload["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert payload["distance_metric"] == "k80"
    assert payload["matrix"][0][0] == 0.0
    assert payload["matrix"][0][1] == pytest.approx(kimura_distance("AAAA", "GAAA"))
    assert payload["matrix"][1][0] == pytest.approx(kimura_distance("GAAA", "AAAA"))
    with (output / "distance_matrix_k80.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["", "Species_A", "Species_B", "Species_C"]
    assert rows[1][0] == "Species_A"
    assert float(rows[1][1]) == 0.0
    assert float(rows[1][2]) == pytest.approx(kimura_distance("AAAA", "GAAA"))


def test_cli_writes_f81_distance_matrix_files(tmp_path):
    """Confirm the CLI can write F81 matrix outputs without changing file formats.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nAAAACCCC\n>Species_B\nAAAAGGCC\n>Species_C\nAAAATTCC\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output), "-dt", "f81"]) == 0

    payload = json.loads((output / "distance_matrix_f81.json").read_text(encoding="utf-8"))
    assert payload["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert payload["distance_metric"] == "f81"
    assert payload["matrix"][0][0] == 0.0
    assert payload["matrix"][0][1] == pytest.approx(f81_distance("AAAACCCC", "AAAAGGCC"))
    assert payload["matrix"][1][0] == pytest.approx(f81_distance("AAAAGGCC", "AAAACCCC"))
    with (output / "distance_matrix_f81.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["", "Species_A", "Species_B", "Species_C"]
    assert rows[1][0] == "Species_A"
    assert float(rows[1][1]) == 0.0
    assert float(rows[1][2]) == pytest.approx(f81_distance("AAAACCCC", "AAAAGGCC"))


def test_cli_writes_hky85_distance_matrix_files(tmp_path):
    """Confirm the CLI can write HKY85 matrix outputs without changing file formats.

    :return: None.
    """
    fasta = tmp_path / "sequences.fasta"
    output = tmp_path / "output"
    fasta.write_text(
        ">Species_A\nAAAACCCC\n>Species_B\nGGAACCCC\n>Species_C\nAAAATCCC\n",
        encoding="utf-8",
    )

    assert main(["-i", str(fasta), "-o", str(output), "-dt", "hky85"]) == 0

    payload = json.loads((output / "distance_matrix_hky85.json").read_text(encoding="utf-8"))
    assert payload["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert payload["distance_metric"] == "hky85"
    assert payload["matrix"][0][0] == 0.0
    assert payload["matrix"][0][1] == pytest.approx(hky85_distance("AAAACCCC", "GGAACCCC"))
    assert payload["matrix"][1][0] == pytest.approx(hky85_distance("GGAACCCC", "AAAACCCC"))
    with (output / "distance_matrix_hky85.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["", "Species_A", "Species_B", "Species_C"]
    assert rows[1][0] == "Species_A"
    assert float(rows[1][1]) == 0.0
    assert float(rows[1][2]) == pytest.approx(hky85_distance("AAAACCCC", "GGAACCCC"))


def test_cli_reports_validation_errors(tmp_path):
    """Confirm CLI validation failures exit with a clear message.

    :return: None.
    """
    with pytest.raises(SystemExit, match="Error: Input FASTA file does not exist"):
        main(["-i", str(tmp_path / "missing.fasta"), "-o", str(tmp_path)])
