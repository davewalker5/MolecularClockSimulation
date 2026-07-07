import csv
import json
import math

import pytest

from distancematrix import (
    calculate_distance_matrix,
    jc69_distance,
    read_fasta,
    write_csv,
    write_json,
)
from distancematrix.calculator import hamming_distance, write_outputs
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
    """Confirm the output helper writes the brief's two expected filenames.

    :return: None.
    """
    matrix = calculate_distance_matrix({
        "Species_A": "ACGT",
        "Species_B": "ACGA",
    })
    written = write_outputs(matrix, tmp_path / "output")

    assert written["json"].name == "distance_matrix.json"
    assert written["csv"].name == "distance_matrix.csv"
    assert written["json"].exists()
    assert written["csv"].exists()


def test_cli_writes_distance_matrix_files(tmp_path, capsys):
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

    captured = capsys.readouterr()
    assert "Wrote JSON" in captured.out
    assert json.loads((output / "distance_matrix.json").read_text(encoding="utf-8"))["matrix"] == [
        [0, 1, 2],
        [1, 0, 3],
        [2, 3, 0],
    ]
    assert (output / "distance_matrix.csv").exists()


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

    payload = json.loads((output / "distance_matrix.json").read_text(encoding="utf-8"))
    assert payload["distance_metric"] == "proportional"
    assert payload["matrix"] == [
        [0.0, 0.125, 0.25],
        [0.125, 0.0, 0.375],
        [0.25, 0.375, 0.0],
    ]
    with (output / "distance_matrix.csv").open(encoding="utf-8", newline="") as handle:
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

    payload = json.loads((output / "distance_matrix.json").read_text(encoding="utf-8"))
    assert payload["labels"] == ["Species_A", "Species_B", "Species_C"]
    assert payload["distance_metric"] == "jc69"
    assert payload["matrix"][0][0] == 0.0
    assert payload["matrix"][0][1] == pytest.approx(jc69_distance("ACGTACGT", "ACGTTCGT"))
    assert payload["matrix"][1][0] == pytest.approx(jc69_distance("ACGTTCGT", "ACGTACGT"))
    with (output / "distance_matrix.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    assert rows[0] == ["", "Species_A", "Species_B", "Species_C"]
    assert rows[1][0] == "Species_A"
    assert float(rows[1][1]) == 0.0
    assert float(rows[1][2]) == pytest.approx(jc69_distance("ACGTACGT", "ACGTTCGT"))


def test_cli_reports_validation_errors(tmp_path):
    """Confirm CLI validation failures exit with a clear message.

    :return: None.
    """
    with pytest.raises(SystemExit, match="Error: Input FASTA file does not exist"):
        main(["-i", str(tmp_path / "missing.fasta"), "-o", str(tmp_path)])
