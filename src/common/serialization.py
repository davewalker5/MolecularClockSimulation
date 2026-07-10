"""Shared text serialization helpers for simulation outputs."""

from __future__ import annotations

import json
from typing import Any, Protocol


class SupportsMetadata(Protocol):
    """Describe results that expose JSON-serializable metadata."""

    def to_metadata(self) -> dict[str, Any]:
        """Return metadata describing a completed result.

        :return: JSON-serializable metadata dictionary.
        """
        ...


def format_fasta(sequences: dict[str, str], line_width: int = 80) -> str:
    """Format named sequences as FASTA text.

    :param sequences: Mapping of sequence name to sequence string.
    :param line_width: Maximum number of characters per FASTA sequence line.
    :return: FASTA-formatted text ending with a trailing newline.
    """
    records: list[str] = []
    for name, sequence in sequences.items():
        records.append(f">{name}")
        # Wrap long sequences to keep the FASTA readable and broadly compatible.
        records.extend(
            sequence[index : index + line_width]
            for index in range(0, len(sequence), line_width)
        )
    return "\n".join(records) + "\n"


def metadata_json(result: SupportsMetadata) -> str:
    """Serialize result metadata as stable formatted JSON.

    :param result: Completed result exposing a ``to_metadata`` method.
    :return: Pretty-printed metadata JSON ending with a newline.
    """
    # Sort keys so downloaded metadata is stable and easier to diff.
    return json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n"
