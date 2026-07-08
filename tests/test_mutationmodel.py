import json
import random
from collections import Counter

import pytest

import biology
import relaxedclock.simulator as relaxed_simulator
import strictclock.simulator as strict_simulator
from biology import BiologySettings, load_biology_settings, mutate_base


def test_transition_weight_increases_transition_frequency():
    """Confirm transition-biased weights make A to G substitutions more frequent.

    :return: None.
    """
    rng = random.Random(123)
    draws = Counter(
        mutate_base("A", transition_weight=2.0, transversion_weight=1.0, rng=rng)
        for _ in range(12000)
    )
    total = sum(draws.values())

    assert draws["G"] / total == pytest.approx(0.5, abs=0.03)
    assert draws["C"] / total == pytest.approx(0.25, abs=0.03)
    assert draws["T"] / total == pytest.approx(0.25, abs=0.03)


def test_equal_weights_preserve_uniform_substitution_behavior():
    """Confirm equal weights reproduce the previous uniform substitution model.

    :return: None.
    """
    rng = random.Random(456)
    draws = Counter(
        mutate_base("A", transition_weight=1.0, transversion_weight=1.0, rng=rng)
        for _ in range(12000)
    )
    total = sum(draws.values())

    assert draws["G"] / total == pytest.approx(1 / 3, abs=0.03)
    assert draws["C"] / total == pytest.approx(1 / 3, abs=0.03)
    assert draws["T"] / total == pytest.approx(1 / 3, abs=0.03)


def test_transition_weight_sensitivity_changes_observed_transition_fraction():
    """Confirm larger transition weights increase observed transition substitutions.

    :return: None.
    """
    uniform_rng = random.Random(789)
    biased_rng = random.Random(789)
    uniform_draws = Counter(
        mutate_base("C", transition_weight=1.0, transversion_weight=1.0, rng=uniform_rng)
        for _ in range(10000)
    )
    biased_draws = Counter(
        mutate_base("C", transition_weight=5.0, transversion_weight=1.0, rng=biased_rng)
        for _ in range(10000)
    )

    uniform_transition_fraction = uniform_draws["T"] / sum(uniform_draws.values())
    biased_transition_fraction = biased_draws["T"] / sum(biased_draws.values())

    assert biased_transition_fraction > uniform_transition_fraction
    assert biased_transition_fraction == pytest.approx(5 / 7, abs=0.03)


def test_biology_settings_reject_non_positive_weights():
    """Confirm substitution weights must be positive.

    :return: None.
    """
    with pytest.raises(ValueError, match="transition_weight"):
        BiologySettings(0.0, 1.0).validate()
    with pytest.raises(ValueError, match="transversion_weight"):
        BiologySettings(1.0, -1.0).validate()


def test_load_biology_settings_reads_transition_transversion_weights(tmp_path):
    """Confirm biology.json weights are loaded into reusable settings.

    :param tmp_path: Temporary directory for an isolated biology config file.
    :return: None.
    """
    biology_path = tmp_path / "biology.json"
    biology_path.write_text(
        json.dumps({"transition_weight": 2.5, "transversion_weight": 0.5}),
        encoding="utf-8",
    )

    settings = load_biology_settings(biology_path)

    assert settings.transition_weight == 2.5
    assert settings.transversion_weight == 0.5


def test_both_simulation_engines_import_the_shared_mutation_routine():
    """Confirm strict and relaxed engines use the shared mutation implementation.

    :return: None.
    """
    assert strict_simulator.mutate_base is biology.mutate_base
    assert relaxed_simulator.mutate_base is biology.mutate_base


def test_strict_engine_mutate_sequence_uses_shared_mutation_routine(monkeypatch):
    """Confirm strict branch evolution delegates base replacement to the shared routine.

    :param monkeypatch: Pytest fixture used to replace the shared mutation routine.
    :return: None.
    """
    calls: list[tuple[str, float, float]] = []

    def fake_mutate_base(
        base: str,
        transition_weight: float,
        transversion_weight: float,
        rng: random.Random | None = None,
        candidates: list[str] | None = None,
    ) -> str:
        """Record strict-engine mutation parameters and return a valid substitution.

        :param base: Ancestral nucleotide to mutate.
        :param transition_weight: Relative transition weight.
        :param transversion_weight: Relative transversion weight.
        :param rng: Optional random generator supplied by the simulator.
        :param candidates: Optional candidate replacements supplied by the simulator.
        :return: Deterministic replacement nucleotide.
        """
        # Store the values passed by the engine so the assertion can verify delegation.
        calls.append((base, transition_weight, transversion_weight))
        return "G"

    monkeypatch.setattr(strict_simulator, "mutate_base", fake_mutate_base)

    sequence, events = strict_simulator.mutate_sequence(
        "AAAA",
        branch_length=1.0,
        mutation_rate=1_000_000.0,
        biology=BiologySettings(3.0, 1.0),
        rng=random.Random(1),
    )

    assert sequence == "GGGG"
    assert len(events) == 4
    assert calls == [("A", 3.0, 1.0)] * 4


def test_relaxed_engine_mutate_sequence_uses_shared_mutation_routine(monkeypatch):
    """Confirm relaxed branch evolution delegates DNA replacement to the shared routine.

    :param monkeypatch: Pytest fixture used to replace the shared mutation routine.
    :return: None.
    """
    calls: list[tuple[str, float, float, tuple[str, ...]]] = []

    def fake_mutate_base(
        base: str,
        transition_weight: float,
        transversion_weight: float,
        rng: random.Random | None = None,
        candidates: list[str] | None = None,
    ) -> str:
        """Record relaxed-engine mutation parameters and return a valid substitution.

        :param base: Ancestral nucleotide to mutate.
        :param transition_weight: Relative transition weight.
        :param transversion_weight: Relative transversion weight.
        :param rng: Optional random generator supplied by the simulator.
        :param candidates: Candidate replacements allowed by relaxed-clock settings.
        :return: Deterministic replacement nucleotide.
        """
        # Candidate recording confirms relaxed-clock filtering is preserved.
        calls.append((base, transition_weight, transversion_weight, tuple(candidates or ())))
        return "G"

    monkeypatch.setattr(relaxed_simulator, "mutate_base", fake_mutate_base)

    sequence, events = relaxed_simulator.mutate_sequence(
        "AAAA",
        root_sequence="CCCC",
        branch_duration=1.0,
        lineage_rate=1_000_000.0,
        alphabet=("A", "C", "G", "T"),
        allow_back_mutation=False,
        biology=BiologySettings(4.0, 1.0),
        rng=random.Random(2),
    )

    assert sequence == "GGGG"
    assert len(events) == 4
    assert calls == [("A", 4.0, 1.0, ("G", "T"))] * 4
