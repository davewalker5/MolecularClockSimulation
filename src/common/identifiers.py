"""Shared deterministic identifier generation."""


class IdentifierCounter:
    """Generate stable incrementing identifiers for simulation objects."""

    def __init__(self) -> None:
        """Create a counter starting before the first identifier.

        :return: None.
        """
        self.value = 0

    def next(self, prefix: str) -> str:
        """Return the next identifier with the requested prefix.

        :param prefix: Prefix describing the identifier type.
        :return: Stable identifier such as ``node_1`` or ``Taxon_2``.
        """
        # Increment before formatting so identifiers consistently begin at one.
        self.value += 1
        return f"{prefix}_{self.value}"
