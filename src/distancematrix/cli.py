"""Command-line interface for the distance matrix calculator."""

from __future__ import annotations

import argparse
from pathlib import Path

from distancematrix.calculator import (
    DISTANCE_TYPES,
    calculate_distance_matrix,
    read_fasta,
    write_outputs,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser.

    :return: Configured parser for the distance matrix calculator command.
    """
    # Keep the CLI surface intentionally small to match the project brief.
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Path to an aligned FASTA file")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Directory where distance_matrix.json and distance_matrix.csv will be written")
    parser.add_argument("-dt", "--distance-type", choices=DISTANCE_TYPES, default="hamming",
                        help="Distance metric to calculate: hamming, proportional, or jc69")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line distance matrix workflow.

    :param argv: Optional list of command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    # Parse paths before doing validation so argparse handles missing arguments.
    args = build_parser().parse_args(argv)

    try:
        # Keep the CLI thin by delegating parsing, calculation, and writing to core code.
        sequences = read_fasta(args.input)
        matrix = calculate_distance_matrix(sequences, args.distance_type)
        written = write_outputs(matrix, args.output)
    except ValueError as error:
        raise SystemExit(f"Error: {error}") from error
    except OSError as error:
        raise SystemExit(f"Error writing distance matrix outputs: {error}") from error

    print(f"Wrote JSON: {written['json']}")
    print(f"Wrote CSV: {written['csv']}")
    return 0
