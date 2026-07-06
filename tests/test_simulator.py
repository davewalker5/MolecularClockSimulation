import json

import pytest

from molecular_clock_simulation import SimulationConfig, run_simulation, write_outputs


def make_config(**overrides):
    data = {
        "sequence_length": 50,
        "number_of_taxa": 8,
        "random_seed": 12345,
        "tree": {"topology": "balanced", "root_age": 1.0},
        "clock": {"model": "strict", "mutation_rate": 0.5},
        "substitution": {"model": "simple"},
    }
    data.update(overrides)
    return SimulationConfig.from_dict(data)


def test_simulation_is_reproducible_with_same_seed():
    first = run_simulation(make_config())
    second = run_simulation(make_config())

    assert first.ancestral_sequence == second.ancestral_sequence
    assert first.newick == second.newick
    assert first.terminal_sequences == second.terminal_sequences
    assert first.to_metadata() == second.to_metadata()


def test_terminal_taxa_are_ultrametric():
    result = run_simulation(make_config())

    distances = {
        node.name: distance_from_root(node, result.root)
        for node in result.root.walk()
        if node.is_leaf
    }

    assert len(distances) == 8
    assert set(round(distance, 10) for distance in distances.values()) == {1.0}


def test_random_topology_keeps_all_taxa_and_records_mutations():
    result = run_simulation(
        make_config(
            number_of_taxa=6,
            tree={"topology": "random", "root_age": 2.0},
            clock={"model": "strict", "mutation_rate": 2.0},
        )
    )

    assert sorted(result.terminal_sequences) == [
        "taxon_1",
        "taxon_2",
        "taxon_3",
        "taxon_4",
        "taxon_5",
        "taxon_6",
    ]
    assert all(len(sequence) == 50 for sequence in result.terminal_sequences.values())
    assert any(
        node.mutations_from_parent for node in result.root.walk() if node.parent_id is not None
    )


def test_write_outputs_creates_fasta_newick_and_metadata(tmp_path):
    result = run_simulation(make_config(sequence_length=12, number_of_taxa=4))
    written = write_outputs(result, "example_config.json", tmp_path)

    fasta = written["fasta"].read_text(encoding="utf-8")
    metadata = json.loads(written["metadata"].read_text(encoding="utf-8"))

    assert written["newick"].read_text(encoding="utf-8").strip().endswith(";")
    assert fasta.count(">") == 4
    assert metadata["ancestral_sequence"] == result.ancestral_sequence
    assert metadata["tree_newick"] == result.newick
    assert len(metadata["nodes"]) == 7


@pytest.mark.parametrize(
    "bad_config,error",
    [
        ({"sequence_length": 0}, "sequence_length"),
        ({"number_of_taxa": 1}, "number_of_taxa"),
        ({"tree": {"topology": "star", "root_age": 1.0}}, "tree.topology"),
        ({"clock": {"model": "relaxed", "mutation_rate": 0.1}}, "clock.model"),
    ],
)
def test_config_validation(bad_config, error):
    data = {
        "sequence_length": 10,
        "number_of_taxa": 2,
        "random_seed": 1,
        "tree": {"topology": "balanced", "root_age": 1.0},
        "clock": {"model": "strict", "mutation_rate": 0.1},
        "substitution": {"model": "simple"},
    }
    data.update(bad_config)

    with pytest.raises(ValueError, match=error):
        SimulationConfig.from_dict(data)


def distance_from_root(node, root):
    by_id = {candidate.id: candidate for candidate in root.walk()}
    distance = 0.0
    current = node
    while current.parent_id is not None:
        distance += current.branch_length
        current = by_id[current.parent_id]
    return distance
