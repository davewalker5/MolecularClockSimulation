"""Reusable distance matrix calculation for aligned FASTA sequences."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read aligned molecular sequences from a FASTA file.

    :param path: Path to the FASTA file containing two or more labelled sequences.
    :return: Mapping of sequence labels to sequence strings.
    """
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise ValueError(f"Input FASTA file does not exist: {fasta_path}")
    if not fasta_path.is_file():
        raise ValueError(f"Input FASTA path is not a file: {fasta_path}")

    # Build the mapping manually because the initial format requirements are small.
    sequences: dict[str, str] = {}
    current_label: str | None = None
    current_lines: list[str] = []

    with fasta_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            # Blank lines are ignored so exported FASTA files with spacing still parse.
            if not line:
                continue

            if line.startswith(">"):
                _store_sequence(sequences, current_label, current_lines)
                current_label = line[1:].strip()
                current_lines = []

                # Empty labels cannot be used as stable row and column identifiers.
                if not current_label:
                    raise ValueError(f"FASTA label on line {line_number} is empty")
                if current_label in sequences:
                    raise ValueError(f"FASTA labels must be unique: {current_label}")
                continue

            if current_label is None:
                raise ValueError("FASTA sequence data appeared before the first label")

            # Uppercase sequence data so lower-case FASTA input compares consistently.
            current_lines.append(line.upper())

    _store_sequence(sequences, current_label, current_lines)
    _validate_sequences(sequences)
    return sequences


def calculate_distance_matrix(sequences: dict[str, str]) -> dict[str, Any]:
    """Calculate a pairwise Hamming distance matrix for aligned sequences.

    :param sequences: Mapping of unique sequence labels to aligned sequence strings.
    :return: Dictionary containing labels, a square distance matrix, and metric name.
    """
    _validate_sequences(sequences)

    # Preserve caller insertion order so labels and rows remain predictable.
    labels = list(sequences)
    matrix: list[list[int]] = []

    for row_label in labels:
        row: list[int] = []
        for column_label in labels:
            # Hamming distance is the count of positions with different characters.
            distance = hamming_distance(sequences[row_label], sequences[column_label])
            row.append(distance)
        matrix.append(row)

    return {
        "labels": labels,
        "matrix": matrix,
        "distance_metric": "hamming",
    }


def hamming_distance(first: str, second: str) -> int:
    """Count differing positions between two aligned sequences.

    :param first: First aligned sequence.
    :param second: Second aligned sequence.
    :return: Number of positions where the two sequences differ.
    """
    if len(first) != len(second):
        raise ValueError("Sequences must have the same length to calculate Hamming distance")

    # zip is safe after the length check and keeps the comparison position-by-position.
    return sum(left != right for left, right in zip(first, second, strict=True))


def write_json(matrix: dict[str, Any], path: str | Path) -> None:
    """Write a distance matrix payload to a JSON file.

    :param matrix: Distance matrix payload returned by calculate_distance_matrix.
    :param path: Destination JSON file path.
    :return: None.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Indented JSON keeps the output easy to inspect and simple for Python to reload.
    output_path.write_text(
        json.dumps(matrix, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_csv(matrix: dict[str, Any], path: str | Path) -> None:
    """Write a distance matrix payload to a CSV file.

    :param matrix: Distance matrix payload returned by calculate_distance_matrix.
    :param path: Destination CSV file path.
    :return: None.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    labels = matrix["labels"]
    rows = matrix["matrix"]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)

        # The first empty cell leaves room for row labels in spreadsheet views.
        writer.writerow(["", *labels])
        for label, row in zip(labels, rows, strict=True):
            writer.writerow([label, *row])


def write_outputs(matrix: dict[str, Any], output_dir: str | Path) -> dict[str, Path]:
    """Write JSON and CSV distance matrix files to an output directory.

    :param matrix: Distance matrix payload returned by calculate_distance_matrix.
    :param output_dir: Directory where distance_matrix.json and distance_matrix.csv are written.
    :return: Mapping of output format names to written file paths.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    json_path = output_path / "distance_matrix.json"
    csv_path = output_path / "distance_matrix.csv"

    # Delegate file-format details to the public writer helpers.
    write_json(matrix, json_path)
    write_csv(matrix, csv_path)

    return {"json": json_path, "csv": csv_path}


def _store_sequence(
    sequences: dict[str, str],
    label: str | None,
    sequence_lines: list[str],
) -> None:
    """Store a parsed FASTA record when a complete record is available.

    :param sequences: Mapping being populated with parsed FASTA records.
    :param label: Label for the record, or None before the first record is found.
    :param sequence_lines: Sequence lines collected for the current record.
    :return: None.
    """
    if label is None:
        return
    if not sequence_lines:
        raise ValueError(f"FASTA record has no sequence data: {label}")

    # Multi-line FASTA records are concatenated into one aligned sequence string.
    sequences[label] = "".join(sequence_lines)


def _validate_sequences(sequences: dict[str, str]) -> None:
    """Validate sequence count, labels, and alignment length.

    :param sequences: Mapping of sequence labels to sequence strings.
    :return: None.
    """
    if len(sequences) < 2:
        raise ValueError("At least two FASTA sequences are required")
    if any(not label for label in sequences):
        raise ValueError("All FASTA sequences must have labels")
    if len(set(sequences)) != len(sequences):
        raise ValueError("All FASTA labels must be unique")
    if any(not sequence for sequence in sequences.values()):
        raise ValueError("All FASTA sequences must contain sequence data")

    # Aligned input requires one shared sequence length across every record.
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("All FASTA sequences must be the same length")
