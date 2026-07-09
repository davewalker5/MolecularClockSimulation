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

    try:
        # The loader checks that the file explicitly targets the strict clock simulator.
        config = load_config(args.config)
        result = run_simulation(config)
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error

    # Output paths are printed so shell users can see exactly where artifacts landed.
    _ = write_outputs(result, args.config, args.output_dir)

    return 0
