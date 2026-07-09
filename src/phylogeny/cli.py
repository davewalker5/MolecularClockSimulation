"""Command-line interface for UPGMA tree reconstruction."""

from __future__ import annotations

import argparse
from pathlib import Path

from phylogeny.upgma import load_distance_matrix, upgma, write_newick


def build_parser() -> argparse.ArgumentParser:
    """Build the UPGMA command-line argument parser.

    :return: Configured parser for the UPGMA command.
    """
    # Keep the initial interface limited to one matrix input and one tree output.
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="Path to a distance_matrix_<type>.json file")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Path where the Newick tree will be written")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line UPGMA workflow.

    :param argv: Optional command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    args = build_parser().parse_args(argv)

    try:
        # Keep reconstruction reusable by delegating all substantive work to the core module.
        labels, matrix = load_distance_matrix(args.input)
        output_path = write_newick(upgma(labels, matrix), args.output)
    except (ValueError, OSError) as error:
        raise SystemExit(f"Error: {error}") from error

    print(f"Wrote Newick tree: {output_path}")
    return 0
