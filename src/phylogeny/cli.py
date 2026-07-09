"""Command-line interface for phylogenetic tree reconstruction."""

from __future__ import annotations

import argparse
from pathlib import Path

from phylogeny.neighbor_joining import neighbor_joining
from phylogeny.upgma import load_distance_matrix_with_metric, upgma, write_newick


def build_parser() -> argparse.ArgumentParser:
    """Build the phylogeny command-line argument parser.

    :return: Configured parser for the reconstruction command.
    """
    # Only the matrix is required because method and metric determine the default output name.
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="Path to a distance_matrix_<type>.json file")
    parser.add_argument("-o", "--output", type=Path, help="Newick output path")
    parser.add_argument(
        "--method",
        choices=("upgma", "nj"),
        required=True,
        help="Tree reconstruction method",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line reconstruction workflow.

    :param argv: Optional command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    args = build_parser().parse_args(argv)

    try:
        # Read the metric with the matrix so an omitted output has a deterministic filename.
        labels, matrix, distance_metric = load_distance_matrix_with_metric(args.input)
        method = args.method
        destination = (
            args.output
            if args.output is not None
            else args.input.parent / f"{method}_{distance_metric}.newick"
        )
        # Dispatch through a tiny method table so adding future reconstruction methods is local.
        reconstruct = {
            "upgma": upgma,
            "nj": neighbor_joining,
        }[method]
        _ = write_newick(reconstruct(labels, matrix), destination)
    except (ValueError, OSError) as error:
        raise SystemExit(f"Error: {error}") from error

    return 0
