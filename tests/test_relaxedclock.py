import json

import pytest

from relaxedclock.cli import main
from relaxedclock import RelaxedClockConfig, run_simulation, write_outputs


def make_config(**overrides):
    """Build a valid relaxed clock config with optional top-level overrides.

    :param overrides: Top-level config sections to replace.
    :return: Validated relaxed clock configuration.
    """
    data = {
        "simulation": {
            "name": "test-relaxed-clock",
            "random_seed": 12345,
        },
        "sequence": {
            "length": 80,
            "alphabet": ["A", "C", "G", "T"],
            "root_sequence": None,
        },
        "tree": {
            "max_depth": 3,
            "branching_mode": "binary",
            "branch_duration": 1.0,
            "duration_jitter": 0.0,
        },
        "clock": {
            "model": "autocorrelated_relaxed",
            "root_rate": 0.2,
            "rate_distribution": "lognormal",
            "rate_sigma": 0.35,
            "minimum_rate": 0.01,
            "maximum_rate": 0.8,
        },
        "mutation": {
            "model": "simple_substitution",
            "allow_back_mutation": True,
        },
        "outputs": {
            "write_fasta": True,
            "write_newick": True,
            "write_metadata": True,
            "newick_branch_lengths": "genetic_change",
        },
    }
    data.update(overrides)
    return RelaxedClockConfig.from_dict(data)


def test_relaxed_simulation_is_reproducible_with_same_seed():
    """Confirm seeded relaxed clock runs reproduce exactly.

    :return: None.
    """
    first = run_simulation(make_config())
    second = run_simulation(make_config())

    assert first.root_sequence == second.root_sequence
    assert first.newick == second.newick
    assert first.terminal_sequences == second.terminal_sequences
    assert first.to_metadata() == second.to_metadata()


def test_zero_duration_jitter_preserves_time_ultrametric_tree_with_variable_rates():
    """Confirm time depths remain equal while genetic depths vary.

    :return: None.
    """
    result = run_simulation(make_config())

    time_distances = {
        node.name: distance_from_root(node, result.root, "branch_duration")
        for node in result.root.walk()
        if node.is_leaf
    }
    genetic_distances = {
        node.name: distance_from_root(node, result.root, "genetic_change")
        for node in result.root.walk()
        if node.is_leaf
    }
    lineage_rates = {
        round(node.lineage_rate, 8)
        for node in result.root.walk()
        if node.parent_id is not None
    }

    assert len(time_distances) == 8
    assert set(round(distance, 10) for distance in time_distances.values()) == {3.0}
    assert len({round(distance, 8) for distance in genetic_distances.values()}) > 1
    assert len(lineage_rates) > 1


def test_metadata_includes_branch_rates_substitutions_and_summary():
    """Confirm metadata carries the relaxed-clock fields from the brief.

    :return: None.
    """
    result = run_simulation(make_config(sequence={"length": 20, "alphabet": ["A", "C", "G", "T"], "root_sequence": "A" * 20}))
    metadata = result.to_metadata()
    branch_node = next(node for node in metadata["nodes"] if node["parent"] is not None)

    assert metadata["configuration"]["clock"]["model"] == "autocorrelated_relaxed"
    assert metadata["root_sequence"] == "A" * 20
    assert metadata["generated_taxa"] == sorted(metadata["terminal_sequences"])
    assert metadata["tree_newick"] == result.newick
    assert "branch_duration" in branch_node
    assert "lineage_rate" in branch_node
    assert "genetic_change" in branch_node
    assert "expected_substitutions" in branch_node
    assert "observed_substitutions" in branch_node
    assert metadata["summary"]["terminal_taxa"] == 8
    assert metadata["summary"]["newick_branch_lengths"] == "genetic_change"


def test_write_outputs_honors_output_flags(tmp_path):
    """Confirm output writer creates only requested relaxed clock files.

    :return: None.
    """
    config = make_config(
        outputs={
            "write_fasta": True,
            "write_newick": False,
            "write_metadata": True,
            "newick_branch_lengths": "genetic_change",
        }
    )
    result = run_simulation(config)
    written = write_outputs(result, "relaxed_example.json", tmp_path)

    assert set(written) == {"FASTA", "metadata"}
    assert written["FASTA"] == tmp_path / "relaxed_example" / "terminal_sequences.fasta"
    assert written["metadata"] == tmp_path / "relaxed_example" / "simulation_metadata.json"
    assert not (tmp_path / "relaxed_example" / "true_tree.newick").exists()
    assert json.loads(written["metadata"].read_text(encoding="utf-8"))["summary"]


def test_newick_can_use_time_branch_lengths():
    """Confirm Newick export can represent elapsed time instead of genetic change.

    :return: None.
    """
    result = run_simulation(
        make_config(
            outputs={
                "write_fasta": True,
                "write_newick": True,
                "write_metadata": True,
                "newick_branch_lengths": "time",
            }
        )
    )

    assert ":1" in result.newick
    assert result.to_metadata()["summary"]["newick_branch_lengths"] == "time"


def test_module_style_cli_accepts_config_option(tmp_path):
    """Confirm the relaxed CLI mirrors strictclock's --config style.

    :param tmp_path: Temporary directory for config and generated output files.
    :return: None.
    """
    config_path = tmp_path / "relaxed_cli.json"
    output_dir = tmp_path / "output"
    config_path.write_text(
        json.dumps(
            {
                "simulation": {"name": "cli-test", "random_seed": 7},
                "sequence": {
                    "length": 12,
                    "alphabet": ["A", "C", "G", "T"],
                    "root_sequence": None,
                },
                "tree": {
                    "max_depth": 1,
                    "branching_mode": "binary",
                    "branch_duration": 1.0,
                    "duration_jitter": 0.0,
                },
                "clock": {
                    "model": "autocorrelated_relaxed",
                    "root_rate": 0.1,
                    "rate_distribution": "lognormal",
                    "rate_sigma": 0.1,
                    "minimum_rate": 0.01,
                    "maximum_rate": 0.5,
                },
                "mutation": {
                    "model": "simple_substitution",
                    "allow_back_mutation": True,
                },
                "outputs": {
                    "write_fasta": True,
                    "write_newick": True,
                    "write_metadata": True,
                    "newick_branch_lengths": "genetic_change",
                },
            }
        ),
        encoding="utf-8",
    )

    exit_code = main(["--config", str(config_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "relaxed_cli" / "terminal_sequences.fasta").exists()
    assert (output_dir / "relaxed_cli" / "true_tree.newick").exists()
    assert (output_dir / "relaxed_cli" / "simulation_metadata.json").exists()


@pytest.mark.parametrize(
    "bad_config,error",
    [
        ({"tree": {"max_depth": 0, "branching_mode": "binary", "branch_duration": 1.0, "duration_jitter": 0.0}}, "tree.max_depth"),
        ({"clock": {"model": "strict", "root_rate": 0.1, "rate_distribution": "lognormal", "rate_sigma": 0.1, "minimum_rate": 0.01, "maximum_rate": 1.0}}, "clock.model"),
        ({"outputs": {"write_fasta": True, "write_newick": True, "write_metadata": True, "newick_branch_lengths": "height"}}, "outputs.newick_branch_lengths"),
        ({"sequence": {"length": 4, "alphabet": ["A"], "root_sequence": None}}, "sequence.alphabet"),
    ],
)
def test_relaxed_config_validation(bad_config, error):
    """Confirm invalid relaxed clock configuration fails clearly.

    :param bad_config: Top-level section override that should fail validation.
    :param error: Expected error message fragment.
    :return: None.
    """
    with pytest.raises(ValueError, match=error):
        make_config(**bad_config)


def distance_from_root(node, root, attribute):
    """Measure a terminal node distance by summing one branch attribute.

    :param node: Descendant node whose distance should be measured.
    :param root: Root node of the tree containing the descendant.
    :param attribute: Branch attribute name to sum while walking to root.
    :return: Total distance from root to the descendant node.
    """
    # Build an id lookup once so parent links can be followed directly.
    by_id = {candidate.id: candidate for candidate in root.walk()}
    distance = 0.0
    current = node
    while current.parent_id is not None:
        distance += getattr(current, attribute)
        current = by_id[current.parent_id]
    return distance
