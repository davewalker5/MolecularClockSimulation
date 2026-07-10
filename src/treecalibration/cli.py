"""Command-line interface for single-point tree calibration."""

from __future__ import annotations

import argparse
from pathlib import Path

from treecalibration.calibration import calibrate_tree_files


def build_parser() -> argparse.ArgumentParser:
    """Build the tree-calibration command-line parser.

    :return: Configured argument parser.
    """
    # The utility intentionally accepts exactly one tree and one calibration file.
    parser = argparse.ArgumentParser(description="Calibrate Newick branch lengths to millions of years.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input Newick tree")
    parser.add_argument("-c", "--calibration", required=True, type=Path, help="Calibration JSON file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output path and shared filename stem",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the tree-calibration command-line workflow.

    :param argv: Optional command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    args = build_parser().parse_args(argv)
    try:
        # Keep file orchestration outside the CLI so it remains reusable.
        calibrate_tree_files(args.input, args.calibration, args.output)
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error
    return 0
