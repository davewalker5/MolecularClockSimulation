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

from strictclock.simulator import (
    SimulationConfig,
    SimulationResult,
    TreeNode,
    format_fasta,
    run_simulation,
)

DOWNLOAD_OPTIONS = (
    "FASTA",
    "Newick",
    "Metadata JSON",
    "Tree PNG",
)


@dataclass(frozen=True)
class ExplorerDefaults:
    """Default values used to initialise the interactive explorer."""

    sequence_length: int = 1000
    number_of_taxa: int = 16
    random_seed: int = 734635
    tree_topology: str = "balanced"
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
    """Create a simulator config from explorer controls."""
    return SimulationConfig.from_dict(
        {
            "sequence_length": sequence_length,
            "number_of_taxa": number_of_taxa,
            "random_seed": random_seed,
            "tree": {"topology": tree_topology, "root_age": root_age},
            "clock": {"model": "strict", "mutation_rate": mutation_rate},
            "substitution": {"model": "simple"},
        }
    )


def summarize_result(result: SimulationResult) -> dict[str, Any]:
    """Build compact summary values for display in the explorer."""
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
    """Count mutation events recorded across all non-root branches."""
    return sum(len(node.mutations_from_parent) for node in root.walk())


def tree_to_dot(root: TreeNode) -> str:
    """Render a simulated tree as Graphviz DOT for Streamlit."""
    lines = [
        "digraph strict_clock_tree {",
        "  graph [rankdir=TB, bgcolor=transparent, margin=0.08];",
        "  node [shape=box, style=\"rounded,filled\", fillcolor=\"#f8fafc\", color=\"#64748b\", fontname=\"Helvetica\", fontsize=11];",
        "  edge [color=\"#64748b\", fontname=\"Helvetica\", fontsize=10];",
    ]
    for node in root.walk():
        label = node_label(node)
        lines.append(f'  "{dot_escape(node.id)}" [label="{dot_escape(label)}"];')
        for child in node.children:
            edge_label = f"{child.branch_length:.4g}"
            lines.append(
                f'  "{dot_escape(node.id)}" -> "{dot_escape(child.id)}" '
                f'[label="{dot_escape(edge_label)}"];'
            )
    lines.append("}")
    return "\n".join(lines)


def node_label(node: TreeNode) -> str:
    """Return a concise display label for one tree node."""
    name = node.name if node.is_leaf else f"internal {node.id}"
    return f"{name}\\nage {node.age:.4g}"


def dot_escape(value: str) -> str:
    """Escape a value for use in a quoted DOT string."""
    return value.replace("\\", "\\\\").replace('"', '\\"')


def fasta_text(result: SimulationResult) -> str:
    """Return terminal sequences as FASTA text."""
    return format_fasta(result.terminal_sequences)


def metadata_json(result: SimulationResult) -> str:
    """Return simulation metadata as formatted JSON text."""
    return json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n"


def tree_png_bytes(dot_source: str) -> bytes:
    """Render Graphviz DOT source to PNG bytes."""
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
    """Validate a user-entered download file stem."""
    stem = value.strip()
    if not stem:
        return None, "Enter a file stem before downloads are available."
    if "/" in stem or "\\" in stem:
        return None, "Enter a file stem only, without a folder or path."
    if stem in {".", ".."}:
        return None, "Enter a file stem, not a relative path marker."
    if PurePath(stem).suffix:
        return None, "Enter the name without a file extension."
    return stem, None


def download_filename(stem: str, extension: str) -> str:
    """Build a download filename from a validated stem and extension."""
    return f"{stem}.{extension}"


def download_payload(
    selection: str,
    result: SimulationResult,
    *,
    fasta: str,
    metadata: str,
    tree_dot: str,
) -> tuple[str | bytes, str, str]:
    """Return data, extension, and MIME type for one download option."""
    if selection == "FASTA":
        return fasta, "fasta", "text/plain"
    if selection == "Newick":
        return result.newick + "\n", "newick", "text/plain"
    if selection == "Metadata JSON":
        return metadata, "json", "application/json"
    if selection == "Tree PNG":
        return tree_png_bytes(tree_dot), "png", "image/png"
    raise ValueError(f"Unknown download selection: {selection}")


def main() -> int:
    """Launch the Streamlit explorer via the installed console script."""
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
    """Render the Streamlit application."""
    import streamlit as st

    defaults = ExplorerDefaults()

    st.set_page_config(
        page_title="Strict Molecular Clock Explorer",
        page_icon="",
        layout="wide",
    )
    st.title("Strict Molecular Clock Explorer")

    with st.sidebar:
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

    sequence_tab, newick_tab, download_tab = st.tabs(
        ["FASTA Sequences", "Newick Output", "Downloads"]
    )
    with sequence_tab:
        st.code(fasta, language="text")
    with newick_tab:
        st.code(result.newick, language="text")
    with download_tab:
        download_selection = st.selectbox(
            "Download",
            options=DOWNLOAD_OPTIONS,
        )
        download_stem_input = st.text_input(
            "File stem",
            value="",
            placeholder="example_run",
            help="Enter the name to use for downloads, without a path or extension.",
        )
        download_stem, download_stem_error = validate_download_stem(download_stem_input)
        if download_stem_error:
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
                download_data, extension, mime = download_payload(
                    download_selection,
                    result,
                    fasta=fasta,
                    metadata=metadata,
                    tree_dot=tree_dot,
                )
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
