"""Command-line interface for the molecular clock simulation generator."""

from __future__ import annotations

import argparse
from pathlib import Path

from strictclock.simulator import load_config, run_simulation, write_outputs

PROJECT_PATH = Path(__file__).parent.parent.parent
DATA_FOLDER = PROJECT_PATH / "data"
DEFAULT_OUTPUT_FOLDER = DATA_FOLDER / "output"


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser.

    :return: Configured parser for the simulation generator command.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, type=Path,
                        help="Path to a JSON simulation configuration file")
    parser.add_argument("-o", "--output-dir", type=Path, default=DEFAULT_OUTPUT_FOLDER,
                        help="Directory where the folder containing FASTA, Newick, and metadata files will be written")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line simulation workflow.

    :param argv: Optional list of command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    # Parse user input before loading the config so argparse can report usage errors.
    args = build_parser().parse_args(argv)

    # The config fully defines the run, including the random seed for reproducibility.
    config = load_config(args.config)
    result = run_simulation(config)

    # Output paths are printed so shell users can see exactly where artifacts landed.
    written = write_outputs(result, args.config, args.output_dir)

    print(f"Wrote FASTA: {written['fasta']}")
    print(f"Wrote Newick: {written['newick']}")
    print(f"Wrote metadata: {written['metadata']}")
    return 0
