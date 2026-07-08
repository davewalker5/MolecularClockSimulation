"""Shared transition/transversion-aware nucleotide mutation helpers."""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

DNA_BASES = ("A", "C", "G", "T")
DEFAULT_EQUILIBRIUM_FREQUENCIES = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25,
}
TRANSITION_PAIRS = frozenset({
    frozenset({"A", "G"}),
    frozenset({"C", "T"}),
})
PROJECT_PATH = Path(__file__).resolve().parents[2]
DEFAULT_BIOLOGY_CONFIG_PATH = PROJECT_PATH / "data" / "config" / "biology.json"


@dataclass(frozen=True)
class BiologySettings:
    transition_weight: float = 1.0
    transversion_weight: float = 1.0
    equilibrium_frequencies: dict[str, float] = field(
        default_factory=lambda: dict(DEFAULT_EQUILIBRIUM_FREQUENCIES)
    )

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "BiologySettings":
        """Create validated biology settings from decoded configuration data.

        :param data: Dictionary containing mutation weighting settings.
        :return: Validated biology settings.
        """
        # Default to equal weights so existing simulations keep their uniform behavior.
        settings = cls(
            transition_weight=float(data.get("transition_weight", 1.0)),
            transversion_weight=float(data.get("transversion_weight", 1.0)),
            equilibrium_frequencies=normalize_equilibrium_frequencies(
                data.get("equilibrium_frequencies", DEFAULT_EQUILIBRIUM_FREQUENCIES)
            ),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate mutation weighting settings.

        :return: None.
        """
        # Weights are relative probabilities, so zero or negative values are invalid.
        if self.transition_weight <= 0:
            raise ValueError("biology.transition_weight must be greater than zero")
        if self.transversion_weight <= 0:
            raise ValueError("biology.transversion_weight must be greater than zero")
        validate_equilibrium_frequencies(self.equilibrium_frequencies)

    def to_dict(self) -> dict[str, Any]:
        """Convert biology settings into JSON-serializable data.

        :return: Dictionary containing mutation weighting settings.
        """
        return {
            "transition_weight": self.transition_weight,
            "transversion_weight": self.transversion_weight,
            "equilibrium_frequencies": dict(self.equilibrium_frequencies),
        }


def load_biology_settings(path: str | Path = DEFAULT_BIOLOGY_CONFIG_PATH) -> BiologySettings:
    """Load mutation weighting settings from a biology JSON file.

    :param path: Path to the biology JSON configuration file.
    :return: Validated biology settings.
    """
    biology_path = Path(path)
    if not biology_path.exists():
        # Missing biology config preserves the historical uniform substitution model.
        return BiologySettings()

    with biology_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    # The biology file must be an object so named weights can be read safely.
    if not isinstance(data, dict):
        raise ValueError("Biology configuration JSON must contain an object")
    return BiologySettings.from_dict(data)


def mutate_base(
    base: str,
    transition_weight: float,
    transversion_weight: float,
    rng: random.Random | None = None,
    candidates: Iterable[str] | None = None,
    equilibrium_frequencies: Mapping[str, float] | None = None,
) -> str:
    """Select a weighted DNA substitution for one nucleotide base.

    :param base: Ancestral nucleotide to mutate.
    :param transition_weight: Relative weight assigned to transition substitutions.
    :param transversion_weight: Relative weight assigned to transversion substitutions.
    :param rng: Optional random number generator used for reproducible simulations.
    :param candidates: Optional replacement bases to consider after engine-specific filtering.
    :param equilibrium_frequencies: Optional target frequencies for A, C, G, and T.
    :return: Derived nucleotide selected from valid substitutions.
    """
    active_rng = rng or random
    settings = BiologySettings(
        transition_weight,
        transversion_weight,
        normalize_equilibrium_frequencies(
            equilibrium_frequencies or DEFAULT_EQUILIBRIUM_FREQUENCIES
        ),
    )
    settings.validate()
    normalized_base = base.upper()
    if normalized_base not in DNA_BASES:
        raise ValueError(f"Cannot apply DNA substitution model to base: {base}")

    # Candidate bases default to every DNA base except the ancestral base.
    if candidates is None:
        possible_replacements = [candidate for candidate in DNA_BASES if candidate != normalized_base]
    else:
        possible_replacements = [
            candidate.upper()
            for candidate in candidates
            if candidate.upper() in DNA_BASES and candidate.upper() != normalized_base
        ]
    if not possible_replacements:
        raise ValueError("At least one valid replacement nucleotide is required")

    weighted_replacements = [
        (
            candidate,
            substitution_weight(
                normalized_base,
                candidate,
                settings.transition_weight,
                settings.transversion_weight,
            ) * settings.equilibrium_frequencies[candidate],
        )
        for candidate in possible_replacements
    ]

    # Draw from the cumulative distribution so relative weights become probabilities.
    total_weight = sum(weight for _, weight in weighted_replacements)
    threshold = active_rng.random() * total_weight
    cumulative = 0.0
    for candidate, weight in weighted_replacements:
        cumulative += weight
        if threshold <= cumulative:
            return candidate

    # Floating-point roundoff can only reach this path at the very end of the interval.
    return weighted_replacements[-1][0]


def normalize_equilibrium_frequencies(frequencies: Mapping[str, Any]) -> dict[str, float]:
    """Normalize configured equilibrium frequencies to uppercase DNA base keys.

    :param frequencies: Mapping containing frequencies for A, C, G, and T.
    :return: Dictionary with uppercase DNA base keys and float frequencies.
    """
    if not isinstance(frequencies, Mapping):
        raise ValueError("biology.equilibrium_frequencies must be an object")

    # Uppercase keys allow JSON authors to use lower-case bases without changing behavior.
    normalized = {
        str(base).upper(): float(frequency)
        for base, frequency in frequencies.items()
    }
    validate_equilibrium_frequencies(normalized)
    return normalized


def validate_equilibrium_frequencies(frequencies: Mapping[str, float]) -> None:
    """Validate configured equilibrium nucleotide frequencies.

    :param frequencies: Mapping containing frequencies for A, C, G, and T.
    :return: None.
    """
    # F81-style DNA equilibrium settings require exactly one frequency per DNA base.
    if set(frequencies) != set(DNA_BASES):
        expected = ", ".join(DNA_BASES)
        raise ValueError(f"biology.equilibrium_frequencies must contain exactly: {expected}")

    for base, frequency in frequencies.items():
        # Positive values keep every base reachable during mutation generation.
        if frequency <= 0:
            raise ValueError(f"biology.equilibrium_frequencies.{base} must be greater than zero")

    # Use a small tolerance so ordinary JSON decimal values are accepted safely.
    total_frequency = sum(frequencies.values())
    if not math_isclose(total_frequency, 1.0):
        raise ValueError("biology.equilibrium_frequencies must sum to 1.0")


def math_isclose(left: float, right: float, tolerance: float = 1e-9) -> bool:
    """Report whether two floating-point values are effectively equal.

    :param left: First value to compare.
    :param right: Second value to compare.
    :param tolerance: Absolute comparison tolerance.
    :return: True when values differ by no more than the tolerance.
    """
    # Avoid importing math for one simple tolerance check in this small helper module.
    return abs(left - right) <= tolerance


def substitution_weight(
    ancestral_base: str,
    derived_base: str,
    transition_weight: float,
    transversion_weight: float,
) -> float:
    """Return the configured weight for one nucleotide substitution.

    :param ancestral_base: Original nucleotide.
    :param derived_base: Replacement nucleotide.
    :param transition_weight: Relative weight assigned to transitions.
    :param transversion_weight: Relative weight assigned to transversions.
    :return: Weight to use when sampling the replacement nucleotide.
    """
    # A substitution is a transition only for A-G or C-T changes.
    if is_transition(ancestral_base, derived_base):
        return transition_weight
    return transversion_weight


def is_transition(ancestral_base: str, derived_base: str) -> bool:
    """Report whether a nucleotide substitution is a transition.

    :param ancestral_base: Original nucleotide.
    :param derived_base: Replacement nucleotide.
    :return: True for A-G or C-T changes, otherwise false.
    """
    # frozenset makes the classification direction-independent.
    return frozenset({ancestral_base.upper(), derived_base.upper()}) in TRANSITION_PAIRS


def is_dna_alphabet(alphabet: Iterable[str]) -> bool:
    """Report whether an alphabet contains exactly the four DNA bases.

    :param alphabet: Iterable of sequence symbols.
    :return: True when the alphabet is A, C, G and T, otherwise false.
    """
    # Relaxed-clock configs support arbitrary alphabets, so DNA-specific bias is opt-in.
    return set(alphabet) == set(DNA_BASES)


__all__ = [
    "BiologySettings",
    "DNA_BASES",
    "DEFAULT_BIOLOGY_CONFIG_PATH",
    "DEFAULT_EQUILIBRIUM_FREQUENCIES",
    "is_dna_alphabet",
    "is_transition",
    "load_biology_settings",
    "mutate_base",
    "normalize_equilibrium_frequencies",
    "substitution_weight",
    "validate_equilibrium_frequencies",
]
