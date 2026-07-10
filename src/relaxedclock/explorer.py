"""Streamlit explorer for relaxed molecular clock simulations."""

from __future__ import annotations

import html
import io
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
from relaxedclock.simulator import (
    RelaxedClockConfig,
    RelaxedClockResult,
    RelaxedTreeNode,
    run_simulation,
    summary_statistics,
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
    synchronize_download_state,
    validate_download_stem,
)

BRANCH_LENGTH_LABELS = {
    "time": "time",
    "rate": "rate",
    "genetic_change": "genetic",
}


@dataclass(frozen=True)
class TreeLayout:
    """Computed coordinates for a scaled phylogram.

    :return: Immutable coordinate data for tree rendering.
    """

    x_by_id: dict[str, float]
    y_by_id: dict[str, float]
    width: int
    height: int
    left_margin: int
    right_margin: int


@dataclass(frozen=True)
class ExplorerDefaults:
    """Default values used to initialise the interactive explorer.

    :return: Immutable group of default UI control values.
    """

    sequence_length: int = 500
    max_depth: int = 3
    random_seed: int = 42
    branch_duration: float = 1.0
    duration_jitter: float = 0.0
    root_rate: float = 0.01
    rate_sigma: float = 0.25
    minimum_rate: float = 0.001
    maximum_rate: float = 0.05
    newick_branch_lengths: str = "genetic_change"
    allow_back_mutation: bool = True


def build_config(
    *,
    sequence_length: int,
    max_depth: int,
    random_seed: int | None,
    branch_duration: float,
    duration_jitter: float,
    root_rate: float,
    rate_sigma: float,
    minimum_rate: float,
    maximum_rate: float,
    newick_branch_lengths: str,
    allow_back_mutation: bool,
) -> RelaxedClockConfig:
    """Create a simulator config from explorer controls.

    :param sequence_length: Number of nucleotide bases to generate.
    :param max_depth: Number of binary branching levels from root to terminal taxa.
    :param random_seed: Optional random seed used to reproduce the simulation.
    :param branch_duration: Base elapsed time assigned to each branch.
    :param duration_jitter: Maximum absolute perturbation applied to branch durations.
    :param root_rate: Substitution rate assigned to the root lineage.
    :param rate_sigma: Lognormal rate inheritance sigma.
    :param minimum_rate: Lower bound for lineage-specific rates.
    :param maximum_rate: Upper bound for lineage-specific rates.
    :param newick_branch_lengths: Branch value used when exporting Newick.
    :param allow_back_mutation: Whether substitutions may return to the root state.
    :return: Validated relaxed clock simulator configuration.
    """
    # Keep the Streamlit controls mapped through the public config validator.
    return RelaxedClockConfig.from_dict(
        {
            "clock_model": "relaxed",
            "simulation": {
                "name": "relaxed-clock-explorer",
                "random_seed": random_seed,
            },
            "sequence": {
                "length": sequence_length,
                "alphabet": ["A", "C", "G", "T"],
                "root_sequence": None,
            },
            "tree": {
                "max_depth": max_depth,
                "branching_mode": "binary",
                "branch_duration": branch_duration,
                "duration_jitter": duration_jitter,
            },
            "clock": {
                "model": "autocorrelated_relaxed",
                "root_rate": root_rate,
                "rate_distribution": "lognormal",
                "rate_sigma": rate_sigma,
                "minimum_rate": minimum_rate,
                "maximum_rate": maximum_rate,
            },
            "mutation": {
                "model": "simple_substitution",
                "allow_back_mutation": allow_back_mutation,
            },
            "outputs": {
                "write_fasta": True,
                "write_newick": True,
                "write_metadata": True,
                "newick_branch_lengths": newick_branch_lengths,
            },
        }
    )


def summarize_result(result: RelaxedClockResult) -> dict[str, Any]:
    """Build compact summary values for display in the explorer.

    :param result: Completed simulation result to summarize.
    :return: Dictionary of scalar values suitable for Streamlit metrics.
    """
    summary = summary_statistics(result)
    return {
        **summary,
        "max_depth": result.config.tree.max_depth,
        "random_seed": result.config.simulation.random_seed,
        "branch_duration": result.config.tree.branch_duration,
        "duration_jitter": result.config.tree.duration_jitter,
        "rate_sigma": result.config.clock.rate_sigma,
        "allow_back_mutation": result.config.mutation.allow_back_mutation,
    }


def count_mutations(root: RelaxedTreeNode) -> int:
    """Count mutation events recorded across all non-root branches.

    :param root: Root node of the simulated tree.
    :return: Total number of mutation events across the tree.
    """
    # Mutation events are stored on child nodes as events from their parent branch.
    return sum(len(node.mutations_from_parent) for node in root.walk())


def tree_to_dot(root: RelaxedTreeNode, branch_lengths: str = "genetic_change") -> str:
    """Render a simulated tree as Graphviz DOT for Streamlit.

    :param root: Root node of the simulated tree.
    :param branch_lengths: Branch value emphasized by the Newick output setting.
    :return: DOT source that can be displayed or rendered by Graphviz.
    """
    # Use left-to-right rank direction so terminal taxa sit on the right.
    colors = DARK_THEME
    lines = [
        "digraph relaxed_clock_tree {",
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
            edge_label = branch_label(child, branch_lengths)
            lines.append(
                f'  "{dot_escape(node.id)}" -> "{dot_escape(child.id)}" '
                f'[label="{dot_escape(edge_label)}"];'
            )
    lines.append("}")
    return "\n".join(lines)


def tree_to_svg(root: RelaxedTreeNode, branch_lengths: str = "genetic_change") -> str:
    """Render a scaled phylogram as SVG.

    :param root: Root node of the simulated tree.
    :param branch_lengths: Branch value used to scale horizontal distances.
    :return: SVG markup for display in Streamlit.
    """
    colors = DARK_THEME
    layout = build_tree_layout(root, branch_lengths)
    lines = [
        f'<svg class="relaxed-clock-tree" viewBox="0 0 {layout.width} {layout.height}" '
        'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="Relaxed molecular clock tree scaled by {html.escape(branch_lengths)}">',
        f'<rect width="100%" height="100%" fill="{colors["page_bg"]}" />',
        '<style>'
        '.relaxed-clock-tree .branch{fill:none;stroke:#94a3b8;stroke-width:2.2;stroke-linecap:round;}'
        '.relaxed-clock-tree .node{fill:#e5e7eb;stroke:#0b1220;stroke-width:1.5;}'
        '.relaxed-clock-tree .leaf{fill:#0b1220;stroke:#e5e7eb;stroke-width:1.8;}'
        '.relaxed-clock-tree .leaf-label{fill:#e5e7eb;font:14px Helvetica,Arial,sans-serif;}'
        '.relaxed-clock-tree .node-label{fill:#94a3b8;font:10px Helvetica,Arial,sans-serif;}'
        '.relaxed-clock-tree .branch-label{fill:#cbd5e1;font:10px Helvetica,Arial,sans-serif;}'
        '</style>',
    ]

    for node in root.walk():
        parent_x = layout.x_by_id[node.id]
        parent_y = layout.y_by_id[node.id]
        for child in node.children:
            child_x = layout.x_by_id[child.id]
            child_y = layout.y_by_id[child.id]
            label_x = parent_x + max((child_x - parent_x) / 2, 8)
            label_y = child_y - 7
            lines.append(
                f'<path class="branch" d="M {parent_x:.2f} {parent_y:.2f} '
                f'V {child_y:.2f} H {child_x:.2f}" />'
            )
            lines.append(
                f'<text class="branch-label" x="{label_x:.2f}" y="{label_y:.2f}">'
                f'{html.escape(compact_branch_label(child, branch_lengths))}</text>'
            )

    for node in root.walk():
        x = layout.x_by_id[node.id]
        y = layout.y_by_id[node.id]
        if node.is_leaf:
            label = node.name or node.id
            lines.append(f'<circle class="leaf" cx="{x:.2f}" cy="{y:.2f}" r="5" />')
            lines.append(
                f'<text class="leaf-label" x="{x + 12:.2f}" y="{y + 5:.2f}">'
                f'{html.escape(label)}</text>'
            )
        else:
            lines.append(f'<circle class="node" cx="{x:.2f}" cy="{y:.2f}" r="4.5" />')
            lines.append(
                f'<text class="node-label" x="{x + 8:.2f}" y="{y - 8:.2f}">'
                f'{html.escape(node.id)}</text>'
            )

    lines.append("</svg>")
    return "\n".join(lines)


def build_tree_layout(
    root: RelaxedTreeNode,
    branch_lengths: str,
    *,
    width: int = 1500,
    row_height: int = 72,
    top_margin: int = 36,
    bottom_margin: int = 36,
    left_margin: int = 28,
    right_margin: int = 220,
) -> TreeLayout:
    """Calculate scaled x/y positions for a rectangular phylogram.

    :param root: Root node of the simulated tree.
    :param branch_lengths: Branch value used to scale horizontal distances.
    :param width: SVG or bitmap width in pixels.
    :param row_height: Vertical spacing between terminal taxa.
    :param top_margin: Empty space above the first terminal taxon.
    :param bottom_margin: Empty space below the last terminal taxon.
    :param left_margin: Empty space before the root.
    :param right_margin: Space reserved for terminal labels.
    :return: Computed tree layout.
    """
    leaves = [node for node in root.walk() if node.is_leaf]
    height = max(top_margin + bottom_margin + row_height * max(len(leaves) - 1, 1), 180)
    distances: dict[str, float] = {root.id: 0.0}

    def assign_distances(node: RelaxedTreeNode) -> None:
        for child in node.children:
            distances[child.id] = distances[node.id] + branch_value(child, branch_lengths)
            assign_distances(child)

    assign_distances(root)
    max_distance = max(distances.values()) or 1.0
    drawable_width = width - left_margin - right_margin
    x_by_id = {
        node_id: left_margin + (distance / max_distance) * drawable_width
        for node_id, distance in distances.items()
    }
    y_by_id: dict[str, float] = {}
    for index, leaf in enumerate(leaves):
        y_by_id[leaf.id] = top_margin + index * row_height

    def assign_internal_y(node: RelaxedTreeNode) -> float:
        if node.is_leaf:
            return y_by_id[node.id]
        child_y_values = [assign_internal_y(child) for child in node.children]
        y_by_id[node.id] = sum(child_y_values) / len(child_y_values)
        return y_by_id[node.id]

    assign_internal_y(root)
    return TreeLayout(
        x_by_id=x_by_id,
        y_by_id=y_by_id,
        width=width,
        height=height,
        left_margin=left_margin,
        right_margin=right_margin,
    )


def branch_value(node: RelaxedTreeNode, branch_lengths: str) -> float:
    """Return the selected branch value for visual scaling.

    :param node: Tree node whose incoming branch should be measured.
    :param branch_lengths: Branch value type, either time, rate, or genetic_change.
    :return: Numeric branch length for the selected mode.
    """
    if branch_lengths == "time":
        return node.branch_duration
    if branch_lengths == "rate":
        return node.lineage_rate
    if branch_lengths == "genetic_change":
        return node.genetic_change
    raise ValueError(f"Unknown branch length mode: {branch_lengths}")


def compact_branch_label(node: RelaxedTreeNode, branch_lengths: str) -> str:
    """Return a compact branch label for the scaled tree.

    :param node: Tree node whose incoming branch should be labelled.
    :param branch_lengths: Branch value emphasized by the Newick output setting.
    :return: Single-line label for display on the branch.
    """
    label = BRANCH_LENGTH_LABELS[branch_lengths]
    return f"{label}: {branch_value(node, branch_lengths):.4g} | obs: {node.observed_substitutions}"


def branch_label(node: RelaxedTreeNode, branch_lengths: str) -> str:
    """Return a graph edge label for one incoming relaxed-clock branch.

    :param node: Tree node whose incoming branch should be labelled.
    :param branch_lengths: Branch value emphasized by the Newick output setting.
    :return: Multiline DOT label for the branch.
    """
    labels = {
        "time": ("Newick time", node.branch_duration),
        "rate": ("Newick rate", node.lineage_rate),
        "genetic_change": ("Newick genetic", node.genetic_change),
    }
    title, value = labels[branch_lengths]
    supporting = [
        ("time", node.branch_duration),
        ("rate", node.lineage_rate),
        ("genetic", node.genetic_change),
    ]
    return "\\n".join(
        [f"{title} {value:.4g}"]
        + [
            f"{label} {supporting_value:.4g}"
            for label, supporting_value in supporting
            if label not in title
        ]
        + [f"obs {node.observed_substitutions}"]
    )


def node_label(node: RelaxedTreeNode) -> str:
    """Return a concise display label for one tree node.

    :param node: Tree node to label.
    :return: Display label containing node identity and depth.
    """
    name = node.name if node.is_leaf else f"internal {node.id}"
    return f"{name}\\ndepth {node.depth}"


def fasta_text(result: RelaxedClockResult) -> str:
    """Return terminal sequences as FASTA text.

    :param result: Completed simulation result.
    :return: FASTA-formatted terminal sequences.
    """
    return format_fasta(result.terminal_sequences)


def tree_png_bytes(
    tree: str | RelaxedTreeNode,
    branch_lengths: str = "genetic_change",
) -> bytes:
    """Render tree graphics to PNG bytes.

    :param tree: DOT source or relaxed tree root to render.
    :param branch_lengths: Branch value used to scale tree graphics.
    :return: PNG image bytes.
    """
    if isinstance(tree, RelaxedTreeNode):
        return scaled_tree_png_bytes(tree, branch_lengths)

    # DOT-based reconstruction graphics use the shared Graphviz renderer.
    return render_dot_png(tree)


def scaled_tree_png_bytes(
    root: RelaxedTreeNode,
    branch_lengths: str = "genetic_change",
) -> bytes:
    """Render the scaled phylogram to PNG bytes.

    :param root: Root node of the simulated tree.
    :param branch_lengths: Branch value used to scale horizontal distances.
    :return: PNG image bytes generated with Pillow.
    """
    from PIL import Image, ImageDraw, ImageFont

    colors = DARK_THEME
    layout = build_tree_layout(root, branch_lengths)
    image = Image.new("RGB", (layout.width, layout.height), colors["page_bg"])
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    small_font = ImageFont.load_default()

    for node in root.walk():
        parent_x = layout.x_by_id[node.id]
        parent_y = layout.y_by_id[node.id]
        for child in node.children:
            child_x = layout.x_by_id[child.id]
            child_y = layout.y_by_id[child.id]
            label_x = parent_x + max((child_x - parent_x) / 2, 8)
            label_y = child_y - 14
            draw.line(
                [(parent_x, parent_y), (parent_x, child_y), (child_x, child_y)],
                fill=colors["text_subtle"],
                width=2,
            )
            draw.text(
                (label_x, label_y),
                compact_branch_label(child, branch_lengths),
                fill=colors["text_muted_strong"],
                font=small_font,
            )

    for node in root.walk():
        x = layout.x_by_id[node.id]
        y = layout.y_by_id[node.id]
        radius = 5 if node.is_leaf else 4
        fill = colors["page_bg"] if node.is_leaf else colors["text"]
        outline = colors["text"] if node.is_leaf else colors["page_bg"]
        draw.ellipse(
            [(x - radius, y - radius), (x + radius, y + radius)],
            fill=fill,
            outline=outline,
            width=2,
        )
        if node.is_leaf:
            draw.text(
                (x + 12, y - 7),
                node.name or node.id,
                fill=colors["text"],
                font=font,
            )
        else:
            draw.text(
                (x + 8, y - 16),
                node.id,
                fill=colors["text_subtle"],
                font=small_font,
            )

    output = io.BytesIO()
    image.save(output, format="PNG")
    return output.getvalue()


def download_payload(
    selection: str,
    result: RelaxedClockResult,
    *,
    fasta: str,
    metadata: str,
    tree_dot: str | RelaxedTreeNode,
    branch_lengths: str = "genetic_change",
    distance_matrix: dict[str, Any] | None = None,
    reconstructed_newick: str | None = None,
    reconstructed_dot: str | None = None,
) -> tuple[str | bytes, str, str]:
    """Return data, extension, and MIME type for one download option.

    :param selection: User-selected download option.
    :param result: Completed simulation result.
    :param fasta: Pre-rendered FASTA text for the result.
    :param metadata: Pre-rendered metadata JSON for the result.
    :param tree_dot: DOT source or relaxed tree root for the generated tree.
    :param branch_lengths: Branch value used to scale tree image downloads.
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
        return tree_png_bytes(tree_dot, branch_lengths), "png", "image/png"
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
        page_title="Relaxed Molecular Clock Explorer",
        page_icon="",
        layout="wide",
    )
    st.markdown(dark_theme_css(), unsafe_allow_html=True)
    st.title("Relaxed Molecular Clock Explorer")

    sidebar_tab_key = "relaxed_sidebar_tab"
    main_tab_key = "relaxed_main_tab"
    reconstruction_warning_key = "relaxed_reconstruction_warning"
    reconstruction_error_key = "relaxed_reconstruction_error"
    reconstruction_dot_key = "relaxed_reconstruction_dot"
    reconstruction_newick_key = "relaxed_reconstruction_newick"
    download_selection_key = "relaxed_download_selection"
    download_stem_key = "relaxed_download_stem"
    download_unavailable_key = "relaxed_download_unavailable"
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
        # Selection changes intentionally overwrite any user-edited file stem.
        synchronize_download_state(
            st.session_state,
            selection_key=download_selection_key,
            stem_key=download_stem_key,
            unavailable_key=download_unavailable_key,
            distance_matrix_key="relaxed_distance_matrix",
            reconstructed_newick_key=reconstruction_newick_key,
            reconstructed_dot_key=reconstruction_dot_key,
            reset_stem=True,
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
        st.subheader("Tree")
        max_depth = st.slider(
            "Tree depth",
            min_value=1,
            max_value=6,
            value=defaults.max_depth,
            step=1,
            help="A full binary tree produces 2^depth terminal taxa.",
        )
        branch_duration = st.number_input(
            "Branch duration",
            min_value=0.01,
            max_value=100.0,
            value=defaults.branch_duration,
            step=0.1,
            format="%.2f",
        )
        duration_jitter = st.number_input(
            "Duration jitter",
            min_value=0.0,
            max_value=max(0.0, branch_duration - 0.01),
            value=min(defaults.duration_jitter, max(0.0, branch_duration - 0.01)),
            step=0.01,
            format="%.2f",
        )
        st.divider()

        st.subheader("Sequence")
        sequence_length = st.slider(
            "Sequence length",
            min_value=10,
            max_value=5000,
            value=defaults.sequence_length,
            step=10,
        )
        st.divider()

        st.subheader("Clock")
        root_rate = st.number_input(
            "Root rate",
            min_value=0.00001,
            max_value=10.0,
            value=defaults.root_rate,
            step=0.001,
            format="%.5f",
        )
        rate_sigma = st.number_input(
            "Rate sigma",
            min_value=0.0,
            max_value=5.0,
            value=defaults.rate_sigma,
            step=0.05,
            format="%.2f",
        )
        minimum_rate = st.number_input(
            "Minimum rate",
            min_value=0.00001,
            max_value=root_rate,
            value=min(defaults.minimum_rate, root_rate),
            step=0.001,
            format="%.5f",
        )
        maximum_rate = st.number_input(
            "Maximum rate",
            min_value=root_rate,
            max_value=10.0,
            value=max(defaults.maximum_rate, root_rate),
            step=0.001,
            format="%.5f",
        )
        st.divider()

        st.subheader("Outputs")
        newick_branch_lengths = st.selectbox(
            "Newick branch lengths",
            options=("genetic_change", "time", "rate"),
            index=("genetic_change", "time", "rate").index(defaults.newick_branch_lengths),
        )
        st.divider()

        st.subheader("Mutation")
        allow_back_mutation = st.checkbox(
            "Allow back mutation",
            value=defaults.allow_back_mutation,
        )
        st.divider()

        st.subheader("Reproducibility")
        random_seed = st.number_input(
            "Random seed",
            min_value=0,
            max_value=2_147_483_647,
            value=defaults.random_seed,
            step=1,
        )
        generate = st.button("Generate", type="primary", width="stretch")

    # Store the latest result in session state so tab changes do not rerun simulation.
    if "relaxed_result" not in st.session_state or generate:
        config = build_config(
            sequence_length=sequence_length,
            max_depth=max_depth,
            random_seed=int(random_seed),
            branch_duration=branch_duration,
            duration_jitter=duration_jitter,
            root_rate=root_rate,
            rate_sigma=rate_sigma,
            minimum_rate=minimum_rate,
            maximum_rate=maximum_rate,
            newick_branch_lengths=newick_branch_lengths,
            allow_back_mutation=allow_back_mutation,
        )
        st.session_state.relaxed_result = run_simulation(config)

    result: RelaxedClockResult = st.session_state.relaxed_result
    summary = summarize_result(result)
    fasta = fasta_text(result)
    metadata = metadata_json(result)
    branch_lengths = result.config.outputs.newick_branch_lengths
    tree_dot = tree_to_dot(result.root, branch_lengths)
    tree_svg = tree_to_svg(result.root, branch_lengths)

    with distance_sidebar_tab:
        st.header("Distance")
        render_distance_analysis_controls(
            result.terminal_sequences,
            state_key_prefix="relaxed",
        )

    with reconstruction_sidebar_tab:
        st.header("Reconstruction")
        st.selectbox(
            "Algorithm",
            options=("UPGMA", "Neighbor Joining"),
            key="relaxed_reconstruction_algorithm",
        )
        if st.button("Reconstruct Tree", key="relaxed_reconstruct_tree", width="stretch"):
            st.session_state[main_tab_key] = "Reconstruction"
            matrix_payload = st.session_state.get("relaxed_distance_matrix")
            if matrix_payload is None:
                st.session_state[reconstruction_warning_key] = True
                st.session_state.pop(reconstruction_error_key, None)
                st.session_state.pop(reconstruction_dot_key, None)
                st.session_state.pop(reconstruction_newick_key, None)
            else:
                try:
                    reconstructed_tree = reconstruct_tree(
                        matrix_payload,
                        method=st.session_state["relaxed_reconstruction_algorithm"],
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
                        graph_name="relaxed_reconstructed_tree",
                        colors=DARK_THEME,
                    )
                    st.session_state[reconstruction_newick_key] = reconstructed_tree_newick(
                        reconstructed_tree
                    )

    # The tree is the primary visual output for this explorer release.
    st.subheader("Phylogenetic Tree")
    st.markdown(tree_svg, unsafe_allow_html=True)

    st.subheader("Simulation Summary")
    metric_columns = st.columns(4)
    metric_columns[0].metric("Taxa", summary["terminal_taxa"])
    metric_columns[1].metric("Sequence length", summary["sequence_length"])
    metric_columns[2].metric("Mutations", summary["total_observed_substitutions"])
    metric_columns[3].metric("Seed", summary["random_seed"])
    rate_columns = st.columns(4)
    rate_columns[0].metric("Mean rate", f"{summary['mean_lineage_rate']:.5g}")
    rate_columns[1].metric("Min rate", f"{summary['minimum_lineage_rate']:.5g}")
    rate_columns[2].metric("Max rate", f"{summary['maximum_lineage_rate']:.5g}")
    rate_columns[3].metric("Expected substitutions", f"{summary['total_expected_substitutions']:.4g}")
    st.caption(
        f"Tree depth: {summary['max_depth']} | "
        f"branch duration: {summary['branch_duration']} | "
        f"duration jitter: {summary['duration_jitter']} | "
        f"root rate: {summary['root_rate']} | "
        f"rate sigma: {summary['rate_sigma']} | "
        f"Newick branch lengths: {html.escape(summary['newick_branch_lengths'])}"
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
            state_key_prefix="relaxed",
        )
    with reconstruction_tab:
        if (
            st.session_state.get(reconstruction_warning_key)
            and st.session_state.get("relaxed_distance_matrix") is None
        ):
            st.warning("You must calculate a distance matrix before reconstructing a tree")
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
                st.session_state.get("relaxed_distance_matrix"),
            ),
        )
        download_selection = st.selectbox(
            "Download",
            options=DOWNLOAD_OPTIONS,
            key=download_selection_key,
            on_change=reset_download_stem,
        )
        unavailable_warning = synchronize_download_state(
            st.session_state,
            selection_key=download_selection_key,
            stem_key=download_stem_key,
            unavailable_key=download_unavailable_key,
            distance_matrix_key="relaxed_distance_matrix",
            reconstructed_newick_key=reconstruction_newick_key,
            reconstructed_dot_key=reconstruction_dot_key,
        )
        download_stem_input = st.text_input(
            "File stem",
            key=download_stem_key,
            help="Enter the name to use for downloads, without a path or extension.",
            disabled=unavailable_warning is not None,
        )
        download_stem, download_stem_error = validate_download_stem(download_stem_input)
        if unavailable_warning:
            st.warning(unavailable_warning)
            st.download_button(
                "Download",
                data="",
                file_name="download_unavailable",
                mime="application/octet-stream",
                width="stretch",
                disabled=True,
            )
        elif download_stem_error:
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
                    tree_dot=result.root,
                    branch_lengths=branch_lengths,
                    distance_matrix=st.session_state.get("relaxed_distance_matrix"),
                    reconstructed_newick=st.session_state.get(reconstruction_newick_key),
                    reconstructed_dot=st.session_state.get(reconstruction_dot_key),
                )
            except ValueError as error:
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
