import json

import pytest

from molecular_clock_simulation import SimulationConfig, run_simulation, write_outputs
from strictclock.cli import main as strict_main
from strictclock.explorer import (
    DARK_THEME,
    count_mutations,
    download_payload,
    download_filename,
    fasta_text,
    dark_theme_css,
    metadata_json,
    summarize_result,
    tree_png_bytes,
    tree_to_dot,
    validate_download_stem,
)


def make_config(**overrides):
    data = {
        "clock_model": "strict",
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


def test_explorer_helpers_render_current_simulation_result():
    result = run_simulation(make_config(sequence_length=20, number_of_taxa=4))

    summary = summarize_result(result)
    dot = tree_to_dot(result.root)
    fasta = fasta_text(result)
    metadata = json.loads(metadata_json(result))

    assert summary["number_of_taxa"] == 4
    assert summary["sequence_length"] == 20
    assert summary["number_of_mutations"] == count_mutations(result.root)
    assert dot.startswith("digraph strict_clock_tree")
    assert "rankdir=LR" in dot
    assert DARK_THEME["page_bg"] in dot
    assert DARK_THEME["surface_elevated"] in dot
    assert "taxon_1" in dot
    assert fasta.count(">") == 4
    assert metadata["tree_newick"] == result.newick


def test_explorer_can_render_tree_png():
    result = run_simulation(make_config(sequence_length=20, number_of_taxa=4))
    png = tree_png_bytes(tree_to_dot(result.root))

    assert png.startswith(b"\x89PNG\r\n\x1a\n")


def test_dark_theme_css_uses_dark_palette():
    css = dark_theme_css()

    assert DARK_THEME["page_bg"] in css
    assert DARK_THEME["text"] in css
    assert DARK_THEME["link"] in css


def test_download_stem_validation_builds_expected_filenames():
    stem, error = validate_download_stem(" example_run ")

    assert stem == "example_run"
    assert error is None
    assert download_filename(stem, "fasta") == "example_run.fasta"
    assert download_filename(stem, "newick") == "example_run.newick"
    assert download_filename(stem, "json") == "example_run.json"
    assert download_filename(stem, "png") == "example_run.png"


@pytest.mark.parametrize(
    "selection,extension,mime",
    [
        ("FASTA", "fasta", "text/plain"),
        ("Newick", "newick", "text/plain"),
        ("Metadata JSON", "json", "application/json"),
        ("Tree PNG", "png", "image/png"),
    ],
)
def test_download_payload_matches_selected_format(selection, extension, mime):
    result = run_simulation(make_config(sequence_length=20, number_of_taxa=4))
    fasta = fasta_text(result)
    metadata = metadata_json(result)
    tree_dot = tree_to_dot(result.root)

    data, actual_extension, actual_mime = download_payload(
        selection,
        result,
        fasta=fasta,
        metadata=metadata,
        tree_dot=tree_dot,
    )

    assert actual_extension == extension
    assert actual_mime == mime
    if selection == "Tree PNG":
        assert data.startswith(b"\x89PNG\r\n\x1a\n")
    else:
        assert isinstance(data, str)
        assert data


@pytest.mark.parametrize("value", ["", "   ", "folder/run", "folder\\run", "run.json", ".", ".."])
def test_download_stem_validation_rejects_missing_paths_and_extensions(value):
    stem, error = validate_download_stem(value)

    assert stem is None
    assert error


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
        "clock_model": "strict",
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


@pytest.mark.parametrize("clock_model", [None, "relaxed"])
def test_strict_config_rejects_wrong_intended_clock_model(clock_model):
    """Confirm strict configs must explicitly identify the strict simulator.

    :param clock_model: Missing or mismatched top-level simulator identifier.
    :return: None.
    """
    data = {
        "sequence_length": 10,
        "number_of_taxa": 2,
        "random_seed": 1,
        "tree": {"topology": "balanced", "root_age": 1.0},
        "clock": {"model": "strict", "mutation_rate": 0.1},
        "substitution": {"model": "simple"},
    }
    if clock_model is not None:
        # Omit the key for the None case to test that it is required.
        data["clock_model"] = clock_model

    with pytest.raises(ValueError, match="clock_model must be 'strict'"):
        SimulationConfig.from_dict(data)


def test_strict_cli_exits_for_relaxed_clock_config(tmp_path):
    """Confirm the strict CLI exits clearly when given a relaxed clock file.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    config_path = tmp_path / "relaxed.json"
    config_path.write_text('{"clock_model": "relaxed"}', encoding="utf-8")

    with pytest.raises(SystemExit, match="clock_model must be 'strict'"):
        strict_main(["--config", str(config_path), "--output-dir", str(tmp_path)])


def distance_from_root(node, root):
    by_id = {candidate.id: candidate for candidate in root.walk()}
    distance = 0.0
    current = node
    while current.parent_id is not None:
        distance += current.branch_length
        current = by_id[current.parent_id]
    return distance
