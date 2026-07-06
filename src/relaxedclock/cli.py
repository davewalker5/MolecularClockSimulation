"""Command-line interface for the relaxed molecular clock simulator."""

from __future__ import annotations

import argparse
from pathlib import Path

from relaxedclock.simulator import load_config, run_simulation, write_outputs

PROJECT_PATH = Path(__file__).parent.parent.parent
DATA_FOLDER = PROJECT_PATH / "data"
DEFAULT_OUTPUT_FOLDER = DATA_FOLDER / "output"


def build_parser() -> argparse.ArgumentParser:
    """Build the relaxed clock command-line argument parser.

    :return: Configured parser for relaxed clock simulation runs.
    """
    parser = argparse.ArgumentParser(
        description="Simulate sequence evolution under a relaxed molecular clock.",
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        type=Path,
        help="Path to a relaxed molecular clock JSON configuration file",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_FOLDER,
        help="Directory where relaxed clock output files will be written",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the relaxed clock command-line workflow.

    :param argv: Optional list of command-line arguments, used mainly by tests.
    :return: Process exit code, where zero indicates success.
    """
    # Parse arguments before reading files so argparse handles usage errors.
    args = build_parser().parse_args(argv)

    # The config controls reproducibility, topology, rate inheritance, and outputs.
    config = load_config(args.config)
    result = run_simulation(config)
    written = write_outputs(result, args.config, args.output_dir)

    # Print only files that were requested and written by the config.
    for label, path in written.items():
        print(f"Wrote {label}: {path}")
    return 0
