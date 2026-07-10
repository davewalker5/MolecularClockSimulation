"""Streamlit explorer for strict molecular clock simulations."""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from typing import Any

from molecular_clock_simulation.distance_analysis import (
    distance_matrix_csv_text,
    render_distance_analysis_controls,
    render_distance_analysis_tab,
)
from molecular_clock_simulation.reconstruction import (
    reconstructed_tree_newick,
    reconstructed_tree_to_dot,
    reconstruct_tree,
)
from molecular_clock_simulation.calibration_ui import (
    clear_calibration_state,
    render_calibration_tab,
)
from strictclock.simulator import (
    SimulationConfig,
    SimulationResult,
    TreeNode,
    run_simulation,
)

from common import (
    DOWNLOAD_TERMINAL_SEQUENCES,
    DOWNLOAD_TRUE_TREE_NEWICK,
    DOWNLOAD_TRUE_TREE_PNG,
    DOWNLOAD_SIMULATION_METADATA,
    DOWNLOAD_DISTANCE_MATRIX_JSON,
    DOWNLOAD_DISTANCE_MATRIX_CSV,
    DOWNLOAD_RECONSTRUCTED_TREE_NEWICK,
    DOWNLOAD_RECONSTRUCTED_TREE_PNG,
    DOWNLOAD_OPTIONS,
    DARK_THEME,
    dark_theme_css,
    default_download_stem,
    download_filename,
    dot_escape,
    format_fasta,
    metadata_json,
    render_dot_png,
    validate_download_stem,
)


@dataclass(frozen=True)
class ExplorerDefaults:
    """Default values used to initialise the interactive explorer.

    :return: Immutable group of default UI control values.
    """

    sequence_length: int = 1000
    number_of_taxa: int = 8
    random_seed: int = 734635
    tree_topology: str = "random"
    root_age: float = 1.0
    mutation_rate: float = 0.01


def build_config(
    *,
    sequence_length: int,
    number_of_taxa: int,
    random_seed: int | None,
    tree_topology: str,
    root_age: float,
    mutation_rate: float,
) -> SimulationConfig:
    """Create a simulator config from explorer controls.

    :param sequence_length: Number of nucleotide bases to generate.
    :param number_of_taxa: Number of terminal taxa in the generated tree.
    :param random_seed: Optional random seed used to reproduce the simulation.
    :param tree_topology: Branching model, either balanced or random.
    :param root_age: Age assigned to the root of the tree.
    :param mutation_rate: Strict-clock substitution rate per site per unit time.
    :return: Validated simulator configuration.
    """
    # Keep the Streamlit controls mapped through the public config validator.
    return SimulationConfig.from_dict(
        {
            "clock_model": "strict",
            "sequence_length": sequence_length,
            "number_of_taxa": number_of_taxa,
            "random_seed": random_seed,
            "tree": {"topology": tree_topology, "root_age": root_age},
            "clock": {"model": "strict", "mutation_rate": mutation_rate},
            "substitution": {"model": "simple"},
        }
    )


def summarize_result(result: SimulationResult) -> dict[str, Any]:
    """Build compact summary values for display in the explorer.

    :param result: Completed simulation result to summarize.
    :return: Dictionary of scalar values suitable for Streamlit metrics.
    """
    return {
        "number_of_taxa": result.config.number_of_taxa,
        "sequence_length": result.config.sequence_length,
        "mutation_rate": result.config.mutation_rate,
        "random_seed": result.config.random_seed,
        "tree_topology": result.config.tree_topology,
        "root_age": result.config.root_age,
        "number_of_mutations": count_mutations(result.root),
    }


def count_mutations(root: TreeNode) -> int:
    """Count mutation events recorded across all non-root branches.

    :param root: Root node of the simulated tree.
    :return: Total number of mutation events across the tree.
    """
    # Mutation events are stored on child nodes as events from their parent branch.
    return sum(len(node.mutations_from_parent) for node in root.walk())


def tree_to_dot(root: TreeNode) -> str:
    """Render a simulated tree as Graphviz DOT for Streamlit.

    :param root: Root node of the simulated tree.
    :return: DOT source that can be displayed or rendered by Graphviz.
    """
    # Use left-to-right rank direction so terminal taxa sit on the right.
    colors = DARK_THEME
    lines = [
        "digraph strict_clock_tree {",
        f"  graph [rankdir=LR, bgcolor=\"{colors['page_bg']}\", margin=0.08];",
        "  node [shape=box, style=\"rounded,filled\", "
        f"fillcolor=\"{colors['surface_elevated']}\", "
        f"color=\"{colors['border_strong']}\", "
        f"fontcolor=\"{colors['text']}\", "
        "fontname=\"Helvetica\", fontsize=11];",
        f"  edge [color=\"{colors['text_subtle']}\", "
        f"fontcolor=\"{colors['text_muted_strong']}\", "
        "fontname=\"Helvetica\", fontsize=10];",
    ]
    for node in root.walk():
        label = node_label(node)
        lines.append(f'  "{dot_escape(node.id)}" [label="{dot_escape(label)}"];')
        for child in node.children:
            # Branch length labels expose the ultrametric timing used by the simulator.
            edge_label = f"{child.branch_length:.4g}"
            lines.append(
                f'  "{dot_escape(node.id)}" -> "{dot_escape(child.id)}" '
                f'[label="{dot_escape(edge_label)}"];'
            )
    lines.append("}")
    return "\n".join(lines)


def node_label(node: TreeNode) -> str:
    """Return a concise display label for one tree node.

    :param node: Tree node to label.
    :return: Display label containing node identity and age.
    """
    name = node.name if node.is_leaf else f"internal {node.id}"
    return f"{name}\\nage {node.age:.4g}"


def fasta_text(result: SimulationResult) -> str:
    """Return terminal sequences as FASTA text.

    :param result: Completed simulation result.
    :return: FASTA-formatted terminal sequences.
    """
    return format_fasta(result.terminal_sequences)


tree_png_bytes = render_dot_png


def download_payload(
    selection: str,
    result: SimulationResult,
    *,
    fasta: str,
    metadata: str,
    tree_dot: str,
    distance_matrix: dict[str, Any] | None = None,
    reconstructed_newick: str | None = None,
    reconstructed_dot: str | None = None,
) -> tuple[str | bytes, str, str]:
    """Return data, extension, and MIME type for one download option.

    :param selection: User-selected download option.
    :param result: Completed simulation result.
    :param fasta: Pre-rendered FASTA text for the result.
    :param metadata: Pre-rendered metadata JSON for the result.
    :param tree_dot: DOT source for the generated tree.
    :param distance_matrix: Calculated distance matrix payload, when available.
    :param reconstructed_newick: Reconstructed Newick text, when available.
    :param reconstructed_dot: Reconstructed tree DOT source, when available.
    :return: Tuple of download data, filename extension, and MIME type.
    """
    # Centralize option handling so the UI and tests share one source of truth.
    if selection == DOWNLOAD_TERMINAL_SEQUENCES:
        return fasta, "fasta", "text/plain"
    if selection == DOWNLOAD_TRUE_TREE_NEWICK:
        return result.newick + "\n", "newick", "text/plain"
    if selection == DOWNLOAD_SIMULATION_METADATA:
        return metadata, "json", "application/json"
    if selection == DOWNLOAD_TRUE_TREE_PNG:
        return tree_png_bytes(tree_dot), "png", "image/png"
    if selection == DOWNLOAD_DISTANCE_MATRIX_JSON:
        if distance_matrix is None:
            raise ValueError("You must calculate a distance matrix before downloading it.")
        return json.dumps(distance_matrix, indent=2) + "\n", "json", "application/json"
    if selection == DOWNLOAD_DISTANCE_MATRIX_CSV:
        if distance_matrix is None:
            raise ValueError("You must calculate a distance matrix before downloading it.")
        return distance_matrix_csv_text(distance_matrix), "csv", "text/csv"
    if selection == DOWNLOAD_RECONSTRUCTED_TREE_NEWICK:
        if reconstructed_newick is None:
            raise ValueError("You must reconstruct a tree before downloading it.")
        # Add a trailing newline so the downloaded Newick is a complete text file.
        return reconstructed_newick.rstrip("\n") + "\n", "newick", "text/plain"
    if selection == DOWNLOAD_RECONSTRUCTED_TREE_PNG:
        if reconstructed_dot is None:
            raise ValueError("You must reconstruct a tree before downloading it.")
        # Render the same DOT source displayed in the reconstruction results tab.
        return tree_png_bytes(reconstructed_dot), "png", "image/png"
    raise ValueError(f"Unknown download selection: {selection}")


def main() -> int:
    """Launch the Streamlit explorer via the installed console script.

    :return: Process exit code from the Streamlit command.
    """
    # Delegate to Streamlit's runner so installed users can launch a console script.
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        __file__,
        "--server.headless=true",
        *sys.argv[1:],
    ]
    return subprocess.call(command)


def render_app() -> None:
    """Render the strict explorer as one explicit, conditionally rendered workflow."""
    import streamlit as st

    defaults = ExplorerDefaults()
    st.set_page_config(
        page_title="Strict Molecular Clock Explorer",
        page_icon="",
        layout="wide",
    )
    st.markdown(dark_theme_css(), unsafe_allow_html=True)
    st.title("Strict Molecular Clock Explorer")

    stage = st.segmented_control(
        "Workflow",
        options=("Simulation", "Distance Matrix", "Reconstruction", "Calibration", "Downloads"),
        default="Simulation",
        key="strict_workflow_stage",
        label_visibility="collapsed",
        width="stretch",
    )

    if "result" not in st.session_state:
        initial_config = build_config(
            sequence_length=defaults.sequence_length,
            number_of_taxa=defaults.number_of_taxa,
            random_seed=defaults.random_seed,
            tree_topology=defaults.tree_topology,
            root_age=defaults.root_age,
            mutation_rate=defaults.mutation_rate,
        )
        st.session_state.result = run_simulation(initial_config)

    if stage == "Simulation":
        with st.sidebar:
            st.header("Simulation")
            number_of_taxa = st.slider(
                "Number of taxa", 2, 64, defaults.number_of_taxa, step=1
            )
            sequence_length = st.slider(
                "Sequence length", 10, 5000, defaults.sequence_length, step=10
            )
            mutation_rate = st.number_input(
                "Mutation rate",
                0.00001,
                1.0,
                defaults.mutation_rate,
                step=0.001,
                format="%.5f",
            )
            random_seed = st.number_input(
                "Random seed", 0, 2_147_483_647, defaults.random_seed, step=1
            )
            tree_topology = st.selectbox(
                "Branching model",
                options=("balanced", "random"),
                index=0 if defaults.tree_topology == "balanced" else 1,
            )
            root_age = st.number_input(
                "Root age", 0.01, 100.0, defaults.root_age, step=0.1, format="%.2f"
            )
            generate = st.button("Generate", type="primary", width="stretch")

        if generate:
            config = build_config(
                sequence_length=sequence_length,
                number_of_taxa=number_of_taxa,
                random_seed=int(random_seed),
                tree_topology=tree_topology,
                root_age=root_age,
                mutation_rate=mutation_rate,
            )
            st.session_state.result = run_simulation(config)
            for key in (
                "strict_distance_matrix",
                "strict_distance_model",
                "strict_reconstruction_dot",
                "strict_reconstruction_newick",
                "strict_reconstruction_error",
            ):
                st.session_state.pop(key, None)
            clear_calibration_state(st.session_state, "strict")

    result: SimulationResult = st.session_state.result
    summary = summarize_result(result)
    fasta = fasta_text(result)
    metadata = metadata_json(result)
    tree_dot = tree_to_dot(result.root)

    if stage == "Simulation":
        st.subheader("Simulated Phylogeny")
        st.graphviz_chart(tree_dot, width="stretch")
        columns = st.columns(4)
        columns[0].metric("Taxa", summary["number_of_taxa"])
        columns[1].metric("Sequence length", summary["sequence_length"])
        columns[2].metric("Mutations", summary["number_of_mutations"])
        columns[3].metric("Seed", summary["random_seed"])
        st.subheader("Terminal Sequences")
        st.code(fasta, language="text")
        with st.expander("True tree Newick"):
            st.code(result.newick, language="text")

    elif stage == "Distance Matrix":
        with st.sidebar:
            st.header("Distance Matrix")
            render_distance_analysis_controls(
                result.terminal_sequences,
                state_key_prefix="strict",
            )
        st.subheader("Distance Analysis")
        render_distance_analysis_tab(
            result.terminal_sequences,
            state_key_prefix="strict",
        )

    elif stage == "Reconstruction":
        with st.sidebar:
            st.header("Reconstruction")
            algorithm = st.selectbox(
                "Algorithm",
                options=("UPGMA", "Neighbor Joining"),
                key="strict_workflow_reconstruction_algorithm",
            )
            reconstruct = st.button("Reconstruct Tree", type="primary", width="stretch")

        if reconstruct:
            clear_calibration_state(st.session_state, "strict")
            matrix_payload = st.session_state.get("strict_distance_matrix")
            if matrix_payload is None:
                st.warning("Calculate a distance matrix before reconstructing a tree.")
            else:
                try:
                    reconstructed = reconstruct_tree(matrix_payload, method=algorithm)
                except ValueError as error:
                    st.session_state["strict_reconstruction_error"] = str(error)
                    st.session_state.pop("strict_reconstruction_dot", None)
                    st.session_state.pop("strict_reconstruction_newick", None)
                else:
                    st.session_state.pop("strict_reconstruction_error", None)
                    st.session_state["strict_reconstruction_dot"] = reconstructed_tree_to_dot(
                        reconstructed,
                        graph_name="strict_reconstructed_tree",
                        colors=DARK_THEME,
                    )
                    st.session_state["strict_reconstruction_newick"] = (
                        reconstructed_tree_newick(reconstructed)
                    )

        st.subheader("Reconstructed Phylogeny")
        if st.session_state.get("strict_reconstruction_error"):
            st.warning(st.session_state["strict_reconstruction_error"])
        elif st.session_state.get("strict_reconstruction_dot"):
            st.graphviz_chart(st.session_state["strict_reconstruction_dot"], width="stretch")
            with st.expander("Reconstructed Newick"):
                st.code(st.session_state["strict_reconstruction_newick"], language="text")
        else:
            st.info("Calculate a distance matrix, then reconstruct a tree.")

    elif stage == "Calibration":
        render_calibration_tab(
            st.session_state.get("strict_reconstruction_newick"),
            "strict",
        )

    else:
        st.subheader("Downloads")
        selection = st.selectbox(
            "Download",
            options=DOWNLOAD_OPTIONS,
            key="strict_workflow_download_selection",
        )
        default_stem = default_download_stem(
            selection, st.session_state.get("strict_distance_matrix")
        )
        stem_input = st.text_input("File stem", value=default_stem)
        stem, stem_error = validate_download_stem(stem_input)
        try:
            data, extension, mime = download_payload(
                selection,
                result,
                fasta=fasta,
                metadata=metadata,
                tree_dot=tree_dot,
                distance_matrix=st.session_state.get("strict_distance_matrix"),
                reconstructed_newick=st.session_state.get("strict_reconstruction_newick"),
                reconstructed_dot=st.session_state.get("strict_reconstruction_dot"),
            )
        except (ValueError, RuntimeError, subprocess.CalledProcessError) as error:
            st.warning(str(error))
            st.download_button(
                "Download", data="", file_name="unavailable", disabled=True, width="stretch"
            )
        else:
            if stem_error:
                st.warning(stem_error)
            st.download_button(
                "Download",
                data=data,
                file_name=download_filename(stem, extension),
                mime=mime,
                disabled=stem_error is not None,
                width="stretch",
            )


if __name__ == "__main__":
    render_app()
