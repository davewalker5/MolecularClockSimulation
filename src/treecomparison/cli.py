"""Command-line interface for side-by-side tree comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

from treecomparison.comparison import compare_trees


def build_parser() -> argparse.ArgumentParser:
    """Build the tree comparison command-line parser.

    :return: Configured argument parser.
    """
    # All three paths are required because the utility performs one explicit comparison.
    parser = argparse.ArgumentParser(
        description="Render source and reconstructed Newick trees side by side.",
    )
    parser.add_argument("-s", "--source", required=True, type=Path, help="Source Newick tree")
    parser.add_argument("-r", "--reconstructed", required=True, type=Path, help="Reconstructed Newick tree")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output PNG file")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the tree comparison command-line workflow.

    :param argv: Optional command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    args = build_parser().parse_args(argv)
    try:
        # Keep parsing and rendering reusable by delegating to the core utility.
        output_path = compare_trees(args.source, args.reconstructed, args.output)
    except (OSError, RuntimeError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error

    print(f"Wrote tree comparison: {output_path}")
    return 0
