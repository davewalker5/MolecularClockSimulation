"""Streamlit explorer for strict molecular clock simulations."""

from __future__ import annotations

import html
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import PurePath
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
from strictclock.simulator import (
    SimulationConfig,
    SimulationResult,
    TreeNode,
    format_fasta,
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


def dot_escape(value: str) -> str:
    """Escape a value for use in a quoted DOT string.

    :param value: Raw string value to embed in DOT source.
    :return: Escaped value safe for a quoted DOT string.
    """
    return value.replace("\\", "\\\\").replace('"', '\\"')


def fasta_text(result: SimulationResult) -> str:
    """Return terminal sequences as FASTA text.

    :param result: Completed simulation result.
    :return: FASTA-formatted terminal sequences.
    """
    return format_fasta(result.terminal_sequences)


def metadata_json(result: SimulationResult) -> str:
    """Return simulation metadata as formatted JSON text.

    :param result: Completed simulation result.
    :return: Pretty-printed metadata JSON ending with a newline.
    """
    # Sort keys so downloaded metadata is stable and easier to diff.
    return json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n"


def tree_png_bytes(dot_source: str) -> bytes:
    """Render Graphviz DOT source to PNG bytes.

    :param dot_source: DOT source describing the tree to render.
    :return: PNG image bytes generated by Graphviz.
    """
    # Streamlit can display DOT directly, but downloads need a concrete image file.
    if shutil.which("dot") is None:
        raise RuntimeError("Graphviz 'dot' command is required for PNG export")

    completed = subprocess.run(
        ["dot", "-Tpng"],
        input=dot_source.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return completed.stdout


def validate_download_stem(value: str) -> tuple[str | None, str | None]:
    """Validate a user-entered download file stem.

    :param value: Raw user-entered file stem.
    :return: Tuple containing the cleaned stem and error message, one of which is None.
    """
    stem = value.strip()
    if not stem:
        return None, "Enter a file stem before downloads are available."
    # Reject path-like values because Streamlit download names should be filenames only.
    if "/" in stem or "\\" in stem:
        return None, "Enter a file stem only, without a folder or path."
    if stem in {".", ".."}:
        return None, "Enter a file stem, not a relative path marker."
    if PurePath(stem).suffix:
        return None, "Enter the name without a file extension."
    return stem, None


def download_filename(stem: str, extension: str) -> str:
    """Build a download filename from a validated stem and extension.

    :param stem: Validated file stem without path or extension.
    :param extension: File extension without a leading dot.
    :return: Complete filename for a download.
    """
    return f"{stem}.{extension}"


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
    ]
    return subprocess.call(command)


def render_app() -> None:
    """Render the Streamlit application.

    :return: None.
    """
    # Import Streamlit lazily so helper functions remain testable without app startup.
    import streamlit as st

    defaults = ExplorerDefaults()

    # Configure the page before emitting any visible Streamlit elements.
    st.set_page_config(
        page_title="Strict Molecular Clock Explorer",
        page_icon="",
        layout="wide",
    )
    st.markdown(dark_theme_css(), unsafe_allow_html=True)
    st.title("Strict Molecular Clock Explorer")

    sidebar_tab_key = "strict_sidebar_tab"
    main_tab_key = "strict_main_tab"
    reconstruction_warning_key = "strict_reconstruction_warning"
    reconstruction_error_key = "strict_reconstruction_error"
    reconstruction_dot_key = "strict_reconstruction_dot"
    reconstruction_newick_key = "strict_reconstruction_newick"
    download_selection_key = "strict_download_selection"
    download_stem_key = "strict_download_stem"
    main_tab_labels = [
        "FASTA Sequences",
        "Newick Output",
        "Distance Analysis",
        "Reconstruction",
        "Downloads",
    ]

    def sync_main_tab_from_sidebar() -> None:
        """Select the matching main output tab when the sidebar tab changes.

        :return: None.
        """
        # Read the sidebar's persisted selection before mapping it to an output tab.
        sidebar_tab = st.session_state.get(sidebar_tab_key)
        if sidebar_tab == "Simulation":
            st.session_state[main_tab_key] = "FASTA Sequences"
        elif sidebar_tab == "Distance":
            st.session_state[main_tab_key] = "Distance Analysis"
        elif sidebar_tab == "Reconstruction":
            st.session_state[main_tab_key] = "Reconstruction"

    def reset_download_stem() -> None:
        """Reset the editable file stem after the download selection changes.

        :return: None.
        """
        # Selection changes intentionally replace any file stem edited by the user.
        st.session_state[download_stem_key] = default_download_stem(
            st.session_state[download_selection_key],
            st.session_state.get("strict_distance_matrix"),
        )
    def sync_sidebar_tab_from_main() -> None:
        """Select the matching sidebar tab when the main output tab changes.

        :return: None.
        """
        # Map output tabs back to the sidebar section that controls their content.
        main_tab = st.session_state.get(main_tab_key)
        if main_tab in {"FASTA Sequences", "Newick Output"}:
            st.session_state[sidebar_tab_key] = "Simulation"
        elif main_tab == "Distance Analysis":
            st.session_state[sidebar_tab_key] = "Distance"
        elif main_tab == "Reconstruction":
            st.session_state[sidebar_tab_key] = "Reconstruction"

    with st.sidebar:
        simulation_sidebar_tab, distance_sidebar_tab, reconstruction_sidebar_tab = st.tabs(
            ["Simulation", "Distance", "Reconstruction"],
            key=sidebar_tab_key,
            on_change=sync_main_tab_from_sidebar,
        )

    with simulation_sidebar_tab:
        st.header("Simulation")
        number_of_taxa = st.slider(
            "Number of taxa",
            min_value=2,
            max_value=64,
            value=defaults.number_of_taxa,
            step=1,
        )
        sequence_length = st.slider(
            "Sequence length",
            min_value=10,
            max_value=5000,
            value=defaults.sequence_length,
            step=10,
        )
        mutation_rate = st.number_input(
            "Mutation rate",
            min_value=0.0,
            max_value=10.0,
            value=defaults.mutation_rate,
            step=0.01,
            format="%.5f",
        )
        random_seed = st.number_input(
            "Random seed",
            min_value=0,
            max_value=2_147_483_647,
            value=defaults.random_seed,
            step=1,
        )
        tree_topology = st.selectbox(
            "Branching model",
            options=("balanced", "random"),
            index=0 if defaults.tree_topology == "balanced" else 1,
        )
        root_age = st.number_input(
            "Root age",
            min_value=0.01,
            max_value=100.0,
            value=defaults.root_age,
            step=0.1,
            format="%.2f",
        )
        generate = st.button("Generate", type="primary", width="stretch")

    # Store the latest result in session state so tab changes do not rerun simulation.
    if "result" not in st.session_state or generate:
        config = build_config(
            sequence_length=sequence_length,
            number_of_taxa=number_of_taxa,
            random_seed=int(random_seed),
            tree_topology=tree_topology,
            root_age=root_age,
            mutation_rate=mutation_rate,
        )
        st.session_state.result = run_simulation(config)

    result: SimulationResult = st.session_state.result
    summary = summarize_result(result)
    fasta = fasta_text(result)
    metadata = metadata_json(result)
    tree_dot = tree_to_dot(result.root)

    with distance_sidebar_tab:
        st.header("Distance")
        render_distance_analysis_controls(
            result.terminal_sequences,
            state_key_prefix="strict",
        )

    with reconstruction_sidebar_tab:
        st.header("Reconstruction")
        st.selectbox(
            "Algorithm",
            options=("UPGMA", "Neighbor Joining"),
            key="strict_reconstruction_algorithm",
        )
        if st.button("Reconstruct Tree", key="strict_reconstruct_tree", width="stretch"):
            st.session_state[main_tab_key] = "Reconstruction"
            matrix_payload = st.session_state.get("strict_distance_matrix")
            if matrix_payload is None:
                st.session_state[reconstruction_warning_key] = True
                st.session_state.pop(reconstruction_error_key, None)
                st.session_state.pop(reconstruction_dot_key, None)
                st.session_state.pop(reconstruction_newick_key, None)
            else:
                try:
                    reconstructed_tree = reconstruct_tree(
                        matrix_payload,
                        method=st.session_state["strict_reconstruction_algorithm"],
                    )
                except ValueError as error:
                    st.session_state[reconstruction_warning_key] = False
                    st.session_state[reconstruction_error_key] = str(error)
                    st.session_state.pop(reconstruction_dot_key, None)
                    st.session_state.pop(reconstruction_newick_key, None)
                else:
                    st.session_state[reconstruction_warning_key] = False
                    st.session_state.pop(reconstruction_error_key, None)
                    st.session_state[reconstruction_dot_key] = reconstructed_tree_to_dot(
                        reconstructed_tree,
                        graph_name="strict_reconstructed_tree",
                        colors=DARK_THEME,
                    )
                    st.session_state[reconstruction_newick_key] = reconstructed_tree_newick(
                        reconstructed_tree
                    )

    # The tree is the primary visual output for this explorer release.
    st.subheader("Phylogenetic Tree")
    st.graphviz_chart(tree_dot, width="stretch")

    st.subheader("Simulation Summary")
    metric_columns = st.columns(4)
    metric_columns[0].metric("Taxa", summary["number_of_taxa"])
    metric_columns[1].metric("Sequence length", summary["sequence_length"])
    metric_columns[2].metric("Mutations", summary["number_of_mutations"])
    metric_columns[3].metric("Seed", summary["random_seed"])
    st.caption(
        f"Branching model: {html.escape(summary['tree_topology'])} | "
        f"root age: {summary['root_age']} | "
        f"mutation rate: {summary['mutation_rate']}"
    )

    sequence_tab, newick_tab, distance_tab, reconstruction_tab, download_tab = st.tabs(
        main_tab_labels,
        key=main_tab_key,
        on_change=sync_sidebar_tab_from_main,
    )
    with sequence_tab:
        st.code(fasta, language="text")
    with newick_tab:
        st.code(result.newick, language="text")
    with distance_tab:
        render_distance_analysis_tab(
            result.terminal_sequences,
            state_key_prefix="strict",
        )
    with reconstruction_tab:
        if (
            st.session_state.get(reconstruction_warning_key)
            and st.session_state.get("strict_distance_matrix") is None
        ):
            st.warning("Calculate a distance matrix before reconstructing a tree.")
        elif st.session_state.get(reconstruction_error_key):
            st.warning(st.session_state[reconstruction_error_key])
        elif st.session_state.get(reconstruction_dot_key):
            st.graphviz_chart(st.session_state[reconstruction_dot_key], width="stretch")
            with st.expander("Reconstructed Newick"):
                st.code(st.session_state[reconstruction_newick_key], language="text")
    with download_tab:
        # Initialize the editable stem before Streamlit creates its stateful widget.
        st.session_state.setdefault(
            download_stem_key,
            default_download_stem(
                st.session_state.get(download_selection_key, DOWNLOAD_OPTIONS[0]),
                st.session_state.get("strict_distance_matrix"),
            ),
        )
        download_selection = st.selectbox(
            "Download",
            options=DOWNLOAD_OPTIONS,
            key=download_selection_key,
            on_change=reset_download_stem,
        )
        download_stem_input = st.text_input(
            "File stem",
            key=download_stem_key,
            help="Enter the name to use for downloads, without a path or extension.",
        )
        download_stem, download_stem_error = validate_download_stem(download_stem_input)
        if download_stem_error:
            # Keep the button visible but disabled so the required action is clear.
            st.warning(download_stem_error)
            st.download_button(
                "Download",
                data="",
                file_name="enter_file_stem",
                mime="application/octet-stream",
                width="stretch",
                disabled=True,
            )
        else:
            try:
                # PNG generation can fail when Graphviz is missing or cannot render.
                download_data, extension, mime = download_payload(
                    download_selection,
                    result,
                    fasta=fasta,
                    metadata=metadata,
                    tree_dot=tree_dot,
                    distance_matrix=st.session_state.get("strict_distance_matrix"),
                    reconstructed_newick=st.session_state.get(reconstruction_newick_key),
                    reconstructed_dot=st.session_state.get(reconstruction_dot_key),
                )
            except ValueError as error:
                if download_selection in {
                    DOWNLOAD_DISTANCE_MATRIX_JSON,
                    DOWNLOAD_DISTANCE_MATRIX_CSV,
                    DOWNLOAD_RECONSTRUCTED_TREE_NEWICK,
                    DOWNLOAD_RECONSTRUCTED_TREE_PNG,
                }:
                    if st.button("Download", width="stretch"):
                        st.warning(str(error))
                else:
                    st.warning(str(error))
            except (RuntimeError, subprocess.CalledProcessError) as error:
                st.warning(f"{download_selection} export is unavailable: {error}")
            else:
                st.download_button(
                    "Download",
                    data=download_data,
                    file_name=download_filename(download_stem, extension),
                    mime=mime,
                    width="stretch",
                )


if __name__ == "__main__":
    render_app()
