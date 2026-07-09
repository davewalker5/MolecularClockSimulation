import json

import pytest

from relaxedclock.cli import main
from relaxedclock import RelaxedClockConfig, run_simulation, write_outputs
from relaxedclock.explorer import (
    DARK_THEME,
    build_config,
    count_mutations,
    download_filename,
    download_payload,
    fasta_text,
    dark_theme_css,
    metadata_json,
    summarize_result,
    tree_png_bytes,
    tree_to_dot,
    tree_to_svg,
    validate_download_stem,
)


def make_config(**overrides):
    """Build a valid relaxed clock config with optional top-level overrides.

    :param overrides: Top-level config sections to replace.
    :return: Validated relaxed clock configuration.
    """
    data = {
        "clock_model": "relaxed",
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


def test_relaxed_explorer_helpers_render_current_simulation_result():
    """Confirm relaxed explorer helpers expose tree, sequence, and rate summaries.

    :return: None.
    """
    result = run_simulation(
        make_config(
            sequence={
                "length": 20,
                "alphabet": ["A", "C", "G", "T"],
                "root_sequence": None,
            }
        )
    )

    summary = summarize_result(result)
    dot = tree_to_dot(result.root)
    time_dot = tree_to_dot(result.root, "time")
    svg = tree_to_svg(result.root)
    time_svg = tree_to_svg(result.root, "time")
    fasta = fasta_text(result)
    metadata = json.loads(metadata_json(result))

    assert summary["terminal_taxa"] == 8
    assert summary["sequence_length"] == 20
    assert summary["total_observed_substitutions"] == count_mutations(result.root)
    assert summary["minimum_lineage_rate"] <= summary["mean_lineage_rate"]
    assert summary["mean_lineage_rate"] <= summary["maximum_lineage_rate"]
    assert dot.startswith("digraph relaxed_clock_tree")
    assert "rankdir=LR" in dot
    assert DARK_THEME["page_bg"] in dot
    assert DARK_THEME["surface_elevated"] in dot
    assert "Taxon_1" in dot
    assert "Newick genetic" in dot
    assert "Newick time" in time_dot
    assert time_dot != dot
    assert svg.startswith("<svg")
    assert "genetic" in svg
    assert "time" in time_svg
    assert time_svg != svg
    assert fasta.count(">") == 8
    assert metadata["tree_newick"] == result.newick


def test_relaxed_explorer_build_config_uses_public_validator():
    """Confirm explorer control values build a valid relaxed clock config.

    :return: None.
    """
    config = build_config(
        sequence_length=100,
        max_depth=2,
        random_seed=9,
        branch_duration=1.0,
        duration_jitter=0.1,
        root_rate=0.05,
        rate_sigma=0.25,
        minimum_rate=0.001,
        maximum_rate=0.2,
        newick_branch_lengths="time",
        allow_back_mutation=False,
    )

    assert config.sequence.length == 100
    assert config.tree.max_depth == 2
    assert config.outputs.newick_branch_lengths == "time"
    assert config.mutation.allow_back_mutation is False


def test_relaxed_explorer_can_render_tree_png():
    """Confirm the relaxed explorer can export Graphviz trees as PNG.

    :return: None.
    """
    result = run_simulation(
        make_config(
            sequence={
                "length": 20,
                "alphabet": ["A", "C", "G", "T"],
                "root_sequence": None,
            }
        )
    )
    dot_png = tree_png_bytes(tree_to_dot(result.root))
    scaled_png = tree_png_bytes(result.root, "genetic_change")
    time_png = tree_png_bytes(result.root, "time")

    assert dot_png.startswith(b"\x89PNG\r\n\x1a\n")
    assert scaled_png.startswith(b"\x89PNG\r\n\x1a\n")
    assert time_png.startswith(b"\x89PNG\r\n\x1a\n")
    assert time_png != scaled_png


def test_relaxed_dark_theme_css_uses_dark_palette():
    """Confirm the relaxed explorer uses the same dark visual palette.

    :return: None.
    """
    css = dark_theme_css()

    assert DARK_THEME["page_bg"] in css
    assert DARK_THEME["text"] in css
    assert DARK_THEME["link"] in css


def test_relaxed_download_stem_validation_builds_expected_filenames():
    """Confirm relaxed explorer download filenames are validated consistently.

    :return: None.
    """
    stem, error = validate_download_stem(" relaxed_run ")

    assert stem == "relaxed_run"
    assert error is None
    assert download_filename(stem, "fasta") == "relaxed_run.fasta"
    assert download_filename(stem, "newick") == "relaxed_run.newick"
    assert download_filename(stem, "json") == "relaxed_run.json"
    assert download_filename(stem, "png") == "relaxed_run.png"


@pytest.mark.parametrize(
    "selection,extension,mime",
    [
        ("FASTA", "fasta", "text/plain"),
        ("Newick", "newick", "text/plain"),
        ("Metadata JSON", "json", "application/json"),
        ("Tree PNG", "png", "image/png"),
    ],
)
def test_relaxed_download_payload_matches_selected_format(selection, extension, mime):
    """Confirm relaxed explorer download options return the expected payload types.

    :return: None.
    """
    result = run_simulation(
        make_config(
            sequence={
                "length": 20,
                "alphabet": ["A", "C", "G", "T"],
                "root_sequence": None,
            }
        )
    )
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


@pytest.mark.parametrize(
    "value",
    ["", "   ", "folder/run", "folder\\run", "run.json", ".", ".."],
)
def test_relaxed_download_stem_validation_rejects_missing_paths_and_extensions(value):
    """Confirm invalid relaxed explorer download file stems are rejected.

    :return: None.
    """
    stem, error = validate_download_stem(value)

    assert stem is None
    assert error


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
                "clock_model": "relaxed",
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


@pytest.mark.parametrize("clock_model", [None, "strict"])
def test_relaxed_config_rejects_wrong_intended_clock_model(clock_model):
    """Confirm relaxed configs must explicitly identify the relaxed simulator.

    :param clock_model: Missing or mismatched top-level simulator identifier.
    :return: None.
    """
    data = {
        "simulation": {"name": "test", "random_seed": 1},
        "sequence": {"length": 4, "alphabet": ["A", "C"], "root_sequence": None},
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
        "mutation": {"model": "simple_substitution", "allow_back_mutation": True},
        "outputs": {
            "write_fasta": True,
            "write_newick": True,
            "write_metadata": True,
            "newick_branch_lengths": "genetic_change",
        },
    }
    if clock_model is not None:
        # Omit the key for the None case to test that it is required.
        data["clock_model"] = clock_model

    with pytest.raises(ValueError, match="clock_model must be 'relaxed'"):
        RelaxedClockConfig.from_dict(data)


def test_relaxed_cli_exits_for_strict_clock_config(tmp_path):
    """Confirm the relaxed CLI exits clearly when given a strict clock file.

    :param tmp_path: Temporary directory supplied by pytest.
    :return: None.
    """
    config_path = tmp_path / "strict.json"
    config_path.write_text('{"clock_model": "strict"}', encoding="utf-8")

    with pytest.raises(SystemExit, match="clock_model must be 'relaxed'"):
        main(["--config", str(config_path), "--output-dir", str(tmp_path)])


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
