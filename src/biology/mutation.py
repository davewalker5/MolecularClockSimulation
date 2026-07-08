"""Shared nucleotide mutation helpers."""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

DNA_BASES = ("A", "C", "G", "T")
EXCHANGEABILITY_PAIRS = ("A_C", "A_G", "A_T", "C_G", "C_T", "G_T")
DEFAULT_EQUILIBRIUM_FREQUENCIES = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25,
}
PROJECT_PATH = Path(__file__).resolve().parents[2]
DEFAULT_BIOLOGY_CONFIG_PATH = PROJECT_PATH / "data" / "config" / "biology.json"


@dataclass(frozen=True)
class BiologySettings:
    transition_weight: float = 1.0
    transversion_weight: float = 1.0
    equilibrium_frequencies: dict[str, float] = field(
        default_factory=lambda: dict(DEFAULT_EQUILIBRIUM_FREQUENCIES)
    )
    exchangeability_rates: dict[str, float] | None = None

    def __post_init__(self) -> None:
        """Populate exchangeability rates from legacy weights when needed.

        :return: None.
        """
        if self.exchangeability_rates is None:
            # Frozen dataclasses need object.__setattr__ for post-init compatibility setup.
            object.__setattr__(
                self,
                "exchangeability_rates",
                derive_exchangeability_rates(
                    self.transition_weight,
                    self.transversion_weight,
                ),
            )
        else:
            # Store normalized pair keys so all later lookups share one representation.
            object.__setattr__(
                self,
                "exchangeability_rates",
                normalize_exchangeability_rates(self.exchangeability_rates),
            )

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "BiologySettings":
        """Create validated biology settings from decoded configuration data.

        :param data: Dictionary containing mutation weighting settings.
        :return: Validated biology settings.
        """
        transition_transversion = data.get("transition_transversion", {})
        if transition_transversion is None:
            transition_transversion = {}
        if not isinstance(transition_transversion, Mapping):
            raise ValueError("biology.transition_transversion must be an object")

        # Explicit exchangeabilities are the primary path; legacy weights only fill the gap.
        settings = cls(
            transition_weight=float(
                data.get(
                    "transition_weight",
                    transition_transversion.get("transition_weight", 1.0),
                )
            ),
            transversion_weight=float(
                data.get(
                    "transversion_weight",
                    transition_transversion.get("transversion_weight", 1.0),
                )
            ),
            equilibrium_frequencies=normalize_equilibrium_frequencies(
                data.get("equilibrium_frequencies", DEFAULT_EQUILIBRIUM_FREQUENCIES)
            ),
            exchangeability_rates=(
                normalize_exchangeability_rates(data["exchangeability_rates"])
                if "exchangeability_rates" in data
                else None
            ),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate mutation weighting settings.

        :return: None.
        """
        validate_equilibrium_frequencies(self.equilibrium_frequencies)
        validate_exchangeability_rates(self.exchangeability_rates or {})

    def to_dict(self) -> dict[str, Any]:
        """Convert biology settings into JSON-serializable data.

        :return: Dictionary containing mutation weighting settings.
        """
        return {
            "equilibrium_frequencies": dict(self.equilibrium_frequencies),
            "exchangeability_rates": dict(self.exchangeability_rates or {}),
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
    exchangeability_rates: Mapping[str, float] | None = None,
) -> str:
    """Select a weighted DNA substitution for one nucleotide base.

    :param base: Ancestral nucleotide to mutate.
    :param transition_weight: Legacy transition weight used only when exchangeability rates are omitted.
    :param transversion_weight: Legacy transversion weight used only when exchangeability rates are omitted.
    :param rng: Optional random number generator used for reproducible simulations.
    :param candidates: Optional replacement bases to consider after engine-specific filtering.
    :param equilibrium_frequencies: Optional target frequencies for A, C, G, and T.
    :param exchangeability_rates: Optional reversible nucleotide-pair exchangeability rates.
    :return: Derived nucleotide selected from valid substitutions.
    """
    active_rng = rng or random
    settings = BiologySettings(
        transition_weight,
        transversion_weight,
        normalize_equilibrium_frequencies(
            equilibrium_frequencies or DEFAULT_EQUILIBRIUM_FREQUENCIES
        ),
        (
            normalize_exchangeability_rates(exchangeability_rates)
            if exchangeability_rates is not None
            else None
        ),
    )
    settings.validate()
    return choose_mutation_target(
        base,
        settings.equilibrium_frequencies,
        settings.exchangeability_rates or {},
        active_rng,
        candidates,
    )


def choose_mutation_target(
    source: str,
    equilibrium_frequencies: Mapping[str, float],
    exchangeability_rates: Mapping[str, float],
    rng: random.Random,
    candidates: Iterable[str] | None = None,
) -> str:
    """Choose a mutation target using exchangeability and equilibrium frequencies.

    :param source: Current nucleotide before mutation.
    :param equilibrium_frequencies: Target frequencies for A, C, G, and T.
    :param exchangeability_rates: Reversible nucleotide-pair exchangeability rates.
    :param rng: Random number generator used for reproducible simulations.
    :param candidates: Optional replacement bases to consider after engine-specific filtering.
    :return: Derived nucleotide selected from valid substitutions.
    """
    normalized_frequencies = normalize_equilibrium_frequencies(equilibrium_frequencies)
    normalized_rates = normalize_exchangeability_rates(exchangeability_rates)
    normalized_source = normalize_dna_base(source, "source")

    # Candidate bases default to every DNA base except the source base.
    if candidates is None:
        possible_replacements = [
            candidate for candidate in DNA_BASES if candidate != normalized_source
        ]
    else:
        possible_replacements = normalize_candidate_bases(candidates, normalized_source)

    if not possible_replacements:
        raise ValueError("At least one valid replacement nucleotide is required")

    weighted_replacements = [
        (
            candidate,
            get_exchangeability_rate(
                normalized_source,
                candidate,
                normalized_rates,
            ) * normalized_frequencies[candidate],
        )
        for candidate in possible_replacements
    ]

    # Draw from the cumulative distribution so relative weights become probabilities.
    total_weight = sum(weight for _, weight in weighted_replacements)
    threshold = rng.random() * total_weight
    cumulative = 0.0
    for candidate, weight in weighted_replacements:
        cumulative += weight
        if threshold <= cumulative:
            return candidate

    # Floating-point roundoff can only reach this path at the very end of the interval.
    return weighted_replacements[-1][0]


def get_exchangeability_rate(
    source: str,
    target: str,
    exchangeability_rates: Mapping[str, float],
) -> float:
    """Return the reversible exchangeability rate for a nucleotide pair.

    :param source: Current nucleotide before mutation.
    :param target: Candidate nucleotide after mutation.
    :param exchangeability_rates: Reversible nucleotide-pair exchangeability rates.
    :return: Exchangeability rate shared by both substitution directions.
    """
    normalized_source = normalize_dna_base(source, "source")
    normalized_target = normalize_dna_base(target, "target")
    if normalized_source == normalized_target:
        raise ValueError("source and target nucleotides must not be identical")

    # Sorting the pair makes A_G and G_A resolve to the single configured A_G key.
    pair_key = "_".join(sorted((normalized_source, normalized_target)))
    normalized_rates = normalize_exchangeability_rates(exchangeability_rates)
    return normalized_rates[pair_key]


def derive_exchangeability_rates(
    transition_weight: float,
    transversion_weight: float,
) -> dict[str, float]:
    """Derive exchangeability rates from legacy transition/transversion weights.

    :param transition_weight: Relative weight assigned to transition substitutions.
    :param transversion_weight: Relative weight assigned to transversion substitutions.
    :return: Six reversible exchangeability rates keyed by nucleotide pair.
    """
    # Legacy weights are relative probabilities, so zero or negative values are invalid.
    if transition_weight <= 0:
        raise ValueError("biology.transition_weight must be greater than zero")
    if transversion_weight <= 0:
        raise ValueError("biology.transversion_weight must be greater than zero")

    # A-G and C-T were the old transition class; all other pairs were transversions.
    return {
        "A_C": float(transversion_weight),
        "A_G": float(transition_weight),
        "A_T": float(transversion_weight),
        "C_G": float(transversion_weight),
        "C_T": float(transition_weight),
        "G_T": float(transversion_weight),
    }


def normalize_exchangeability_rates(rates: Mapping[str, Any]) -> dict[str, float]:
    """Normalize configured exchangeability rates to canonical pair keys.

    :param rates: Mapping containing the six reversible DNA pair rates.
    :return: Dictionary with canonical uppercase pair keys and float rates.
    """
    if not isinstance(rates, Mapping):
        raise ValueError("biology.exchangeability_rates must be an object")

    normalized: dict[str, float] = {}
    for pair_key, rate in rates.items():
        # Pair keys are normalized through nucleotide validation so malformed keys fail clearly.
        normalized_key = normalize_exchangeability_pair_key(str(pair_key))
        if normalized_key in normalized:
            raise ValueError(f"Duplicate exchangeability pair key: {pair_key}")
        normalized[normalized_key] = float(rate)

    validate_exchangeability_rates(normalized)
    return normalized


def validate_exchangeability_rates(rates: Mapping[str, float]) -> None:
    """Validate configured reversible nucleotide exchangeability rates.

    :param rates: Mapping containing the six reversible DNA pair rates.
    :return: None.
    """
    # Exchangeability-aware mutation weighting requires one rate per reversible DNA pair.
    if set(rates) != set(EXCHANGEABILITY_PAIRS):
        expected = ", ".join(EXCHANGEABILITY_PAIRS)
        raise ValueError(f"biology.exchangeability_rates must contain exactly: {expected}")

    for pair_key, rate in rates.items():
        # Rates are relative probabilities, so zero or negative values are invalid.
        if rate <= 0:
            raise ValueError(f"biology.exchangeability_rates.{pair_key} must be greater than zero")


def normalize_exchangeability_pair_key(pair_key: str) -> str:
    """Normalize one reversible exchangeability pair key.

    :param pair_key: Pair key such as A_G or g_a.
    :return: Canonical uppercase pair key with alphabetically ordered nucleotides.
    """
    parts = pair_key.upper().split("_")
    if len(parts) != 2:
        raise ValueError(f"Invalid exchangeability pair key: {pair_key}")

    left = normalize_dna_base(parts[0], "exchangeability source")
    right = normalize_dna_base(parts[1], "exchangeability target")
    if left == right:
        raise ValueError("exchangeability pair nucleotides must not be identical")

    return "_".join(sorted((left, right)))


def normalize_candidate_bases(
    candidates: Iterable[str],
    normalized_source: str,
) -> list[str]:
    """Normalize mutation candidate bases while excluding the source base.

    :param candidates: Replacement bases supplied by a simulator.
    :param normalized_source: Uppercase source nucleotide.
    :return: Candidate replacement bases in caller-supplied order.
    """
    possible_replacements: list[str] = []
    for candidate in candidates:
        # Invalid candidates are ignored to preserve previous engine filtering behavior.
        normalized_candidate = str(candidate).upper()
        if normalized_candidate in DNA_BASES and normalized_candidate != normalized_source:
            possible_replacements.append(normalized_candidate)
    return possible_replacements


def normalize_dna_base(base: str, label: str) -> str:
    """Normalize and validate one DNA nucleotide symbol.

    :param base: Nucleotide symbol to normalize.
    :param label: Name used in validation errors.
    :return: Uppercase DNA nucleotide symbol.
    """
    normalized = str(base).upper()
    if normalized not in DNA_BASES:
        raise ValueError(f"Invalid {label} nucleotide: {base}")
    return normalized


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
    "EXCHANGEABILITY_PAIRS",
    "choose_mutation_target",
    "derive_exchangeability_rates",
    "get_exchangeability_rate",
    "is_dna_alphabet",
    "load_biology_settings",
    "mutate_base",
    "normalize_equilibrium_frequencies",
    "normalize_exchangeability_rates",
    "validate_equilibrium_frequencies",
    "validate_exchangeability_rates",
]
