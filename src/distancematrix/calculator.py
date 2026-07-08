"""Reusable distance matrix calculation for aligned FASTA sequences."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

DNA_BASES = frozenset({"A", "C", "G", "T"})
TRANSITION_PAIRS = frozenset({
    frozenset({"A", "G"}),
    frozenset({"C", "T"}),
})
DISTANCE_TYPES = ("hamming", "proportional", "jc69", "k80")


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


def calculate_distance_matrix(
    sequences: dict[str, str],
    distance_type: str = "hamming",
) -> dict[str, Any]:
    """Calculate a pairwise distance matrix for aligned sequences.

    :param sequences: Mapping of unique sequence labels to aligned sequence strings.
    :param distance_type: Distance metric to calculate: hamming, proportional, jc69, or k80.
    :return: Dictionary containing labels, a square distance matrix, and metric name.
    """
    _validate_sequences(sequences)
    _validate_distance_type(distance_type)

    # Preserve caller insertion order so labels and rows remain predictable.
    labels = list(sequences)
    matrix: list[list[int | float]] = []

    for row_label in labels:
        row: list[int | float] = []
        for column_label in labels:
            # Dispatch through a helper so the matrix loop stays independent of the metric.
            distance = calculate_distance(
                sequences[row_label],
                sequences[column_label],
                distance_type,
            )
            row.append(distance)
        matrix.append(row)

    return {
        "labels": labels,
        "matrix": matrix,
        "distance_metric": distance_type,
    }


def calculate_distance(first: str, second: str, distance_type: str) -> int | float:
    """Calculate the selected distance metric between two aligned sequences.

    :param first: First aligned sequence.
    :param second: Second aligned sequence.
    :param distance_type: Distance metric to calculate: hamming, proportional, jc69, or k80.
    :return: Numeric distance for the two sequences under the selected metric.
    """
    _validate_distance_type(distance_type)

    # Keep metric selection data-driven so future substitution models can be registered here.
    distance_functions = {
        "hamming": hamming_distance,
        "proportional": proportional_distance,
        "jc69": jc69_distance,
        "k80": kimura_distance,
    }
    return distance_functions[distance_type](first, second)


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


def proportional_distance(first: str, second: str) -> float:
    """Calculate the proportion of differing positions between two aligned sequences.

    :param first: First aligned sequence.
    :param second: Second aligned sequence.
    :return: Differing positions divided by aligned sequence length.
    """
    if not first:
        raise ValueError("Sequences must not be empty to calculate proportional distance")

    # Reuse Hamming distance so both metrics share one definition of a difference.
    return hamming_distance(first, second) / len(first)


def jc69_distance(first: str, second: str) -> float:
    """Calculate the Jukes-Cantor corrected evolutionary distance.

    :param first: First aligned DNA sequence containing only A, C, G, and T.
    :param second: Second aligned DNA sequence containing only A, C, G, and T.
    :return: Estimated substitutions per site, or infinity when saturated.
    """
    _validate_dna_sequence(first, "first", "jc69")
    _validate_dna_sequence(second, "second", "jc69")

    # JC69 starts from the observed proportional distance between aligned DNA sequences.
    observed_distance = proportional_distance(first, second)
    if observed_distance == 0:
        return 0.0
    if observed_distance >= 0.75:
        return float("inf")

    # The correction accounts for unobserved repeated substitutions at the same site.
    return -0.75 * math.log(1 - (4 * observed_distance / 3))


def kimura_distance(first: str, second: str) -> float:
    """Calculate the Kimura two-parameter corrected evolutionary distance.

    :param first: First aligned DNA sequence containing only A, C, G, and T.
    :param second: Second aligned DNA sequence containing only A, C, G, and T.
    :return: Estimated substitutions per site, or infinity when saturated.
    """
    transitions, transversions = count_k80_substitutions(first, second)
    sequence_length = len(first)

    # Identical sequences have no observed substitutions and therefore zero distance.
    if transitions == 0 and transversions == 0:
        return 0.0

    # K80 separates transition and transversion proportions before correction.
    transition_proportion = transitions / sequence_length
    transversion_proportion = transversions / sequence_length

    # Both logarithm arguments must stay positive for the model to be defined.
    transition_term = 1 - (2 * transition_proportion) - transversion_proportion
    transversion_term = 1 - (2 * transversion_proportion)
    if transition_term <= 0 or transversion_term <= 0:
        return float("inf")

    # The two terms correct hidden substitutions while preserving the transition bias.
    return (
        -0.5 * math.log(transition_term)
        - 0.25 * math.log(transversion_term)
    )


def count_k80_substitutions(first: str, second: str) -> tuple[int, int]:
    """Count transition and transversion substitutions between aligned DNA sequences.

    :param first: First aligned DNA sequence containing only A, C, G, and T.
    :param second: Second aligned DNA sequence containing only A, C, G, and T.
    :return: Tuple containing transition count followed by transversion count.
    """
    _validate_dna_sequence(first, "first", "k80")
    _validate_dna_sequence(second, "second", "k80")
    if not first:
        raise ValueError("Sequences must not be empty to calculate K80 distance")
    if len(first) != len(second):
        raise ValueError("Sequences must have the same length to calculate K80 distance")

    transitions = 0
    transversions = 0

    for left, right in zip(first, second, strict=True):
        # Matching bases are unchanged sites and do not contribute to either count.
        if left == right:
            continue

        # A-G and C-T changes are transitions; every other DNA change is a transversion.
        substitution_pair = frozenset({left, right})
        if substitution_pair in TRANSITION_PAIRS:
            transitions += 1
        else:
            transversions += 1

    return transitions, transversions


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


def _validate_distance_type(distance_type: str) -> None:
    """Validate the requested distance metric.

    :param distance_type: Distance metric name supplied by the caller.
    :return: None.
    """
    # Restrict metrics to named implementations so unsupported requests fail clearly.
    if distance_type not in DISTANCE_TYPES:
        supported = "', '".join(DISTANCE_TYPES)
        raise ValueError(f"distance_type must be one of: '{supported}'")


def _validate_dna_sequence(sequence: str, sequence_name: str, model_name: str) -> None:
    """Validate that a sequence contains only model-compatible DNA bases.

    :param sequence: Sequence string to validate.
    :param sequence_name: Human-readable name used in validation errors.
    :param model_name: Distance model name used in validation errors.
    :return: None.
    """
    # Empty-sequence checks are handled by the distance models that require them.
    invalid_bases = sorted(set(sequence) - DNA_BASES)
    if invalid_bases:
        bases = ", ".join(invalid_bases)
        raise ValueError(f"{sequence_name} sequence contains invalid DNA bases for {model_name}: {bases}")
