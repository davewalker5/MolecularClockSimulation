"""Reusable single-point calibration for reconstructed phylogenetic trees."""

from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
from typing import Any

from treecomparison.comparison import NewickNode, read_newick


def _leaf_names(node: NewickNode) -> set[str]:
    """Collect the labels of all leaves below a node.

    :param node: Root of the subtree to inspect.
    :return: Set of descendant leaf labels.
    """
    # A leaf contributes its own label; internal nodes combine their child sets.
    if node.is_leaf:
        return {node.name} if node.name is not None else set()
    names: set[str] = set()
    for child in node.children:
        names.update(_leaf_names(child))
    return names


def _find_mrca(node: NewickNode, taxa: set[str]) -> NewickNode:
    """Find the smallest subtree containing every requested taxon.

    :param node: Root of the complete tree or current search subtree.
    :param taxa: Taxon labels whose most recent common ancestor is required.
    :return: Node representing the taxa's unique MRCA.
    """
    # Descend while a single child still contains every calibration taxon.
    for child in node.children:
        if taxa.issubset(_leaf_names(child)):
            return _find_mrca(child, taxa)
    return node


def _distance_to_taxon(node: NewickNode, taxon: str, distance: float = 0.0) -> float | None:
    """Measure the branch-length path from a node to one descendant taxon.

    :param node: Current node in the recursive search.
    :param taxon: Descendant leaf label to locate.
    :param distance: Accumulated distance from the starting node.
    :return: Path length when found, otherwise None.
    """
    if node.is_leaf:
        return distance if node.name == taxon else None
    for child in node.children:
        # Every non-root branch was validated before this traversal.
        result = _distance_to_taxon(child, taxon, distance + child.branch_length)  # type: ignore[operator]
        if result is not None:
            return result
    return None


def _validate_tree(root: NewickNode) -> dict[str, NewickNode]:
    """Validate branch lengths and build a unique leaf-name lookup.

    :param root: Root node of the reconstructed tree.
    :return: Mapping from taxon label to leaf node.
    """
    leaves: dict[str, NewickNode] = {}
    for node in root.walk():
        if node is not root:
            if node.branch_length is None:
                raise ValueError("Every non-root branch must have a branch length")
            if not math.isfinite(node.branch_length):
                raise ValueError("Branch lengths must be finite")
            # Neighbor Joining can produce small negative branches when the
            # observed distances are not perfectly additive. Retain these
            # signed estimates so calibration does not alter the topology or
            # silently impose a separate branch-correction procedure.
        if node.is_leaf:
            if node.name is None:
                raise ValueError("Every leaf must have a taxon label")
            if node.name in leaves:
                raise ValueError(f"Taxon label is not unique in the tree: {node.name}")
            leaves[node.name] = node
    return leaves


def calibrate_tree(
    tree: NewickNode,
    calibration_taxa: list[str],
    calibration_age_mya: float,
) -> tuple[NewickNode, dict[str, Any]]:
    """Scale all tree branches using the known age of one MRCA.

    :param tree: Root of a parsed reconstructed Newick tree.
    :param calibration_taxa: At least two unique taxa identifying the MRCA.
    :param calibration_age_mya: Positive calibration age in millions of years.
    :return: Independent calibrated tree and calibration metadata.
    """
    if len(calibration_taxa) < 2:
        raise ValueError("Calibration must specify at least two taxa")
    if len(set(calibration_taxa)) != len(calibration_taxa):
        raise ValueError("Calibration taxa must be unique")
    if isinstance(calibration_age_mya, bool) or not isinstance(calibration_age_mya, (int, float)):
        raise ValueError("Calibration age must be a number")
    if not math.isfinite(calibration_age_mya) or calibration_age_mya <= 0:
        raise ValueError("Calibration age must be finite and greater than zero")

    leaves = _validate_tree(tree)
    missing = [taxon for taxon in calibration_taxa if taxon not in leaves]
    if missing:
        raise ValueError(f"Calibration taxa not found in tree: {', '.join(missing)}")

    # The mean descendant path is exact for an ultrametric tree and remains a
    # simple, transparent approximation when reconstructed tip depths differ.
    mrca = _find_mrca(tree, set(calibration_taxa))
    distances = [_distance_to_taxon(mrca, taxon) for taxon in calibration_taxa]
    if any(distance is None for distance in distances):
        raise ValueError("No unique MRCA could be found for the calibration taxa")
    reconstructed_depth = sum(distance for distance in distances if distance is not None) / len(distances)
    if reconstructed_depth <= 0:
        raise ValueError("Calibration node reconstructed depth must be greater than zero")
    scale_factor = float(calibration_age_mya) / reconstructed_depth

    # Copy first so callers can retain and compare the uncalibrated tree.
    calibrated = deepcopy(tree)
    for node in calibrated.walk():
        if node.branch_length is not None:
            node.branch_length *= scale_factor
    metadata = {
        "calibration_taxa": list(calibration_taxa),
        "calibration_age_mya": float(calibration_age_mya),
        "reconstructed_depth": reconstructed_depth,
        "scale_factor": scale_factor,
        "units": "million years",
    }
    return calibrated, metadata


def _format_label(label: str) -> str:
    """Encode a node label safely for Newick output.

    :param label: Decoded node label.
    :return: Unquoted or apostrophe-quoted Newick label.
    """
    # Quote labels containing Newick delimiters or whitespace and double apostrophes.
    if any(character in "(),:;" or character.isspace() or character == "'" for character in label):
        return "'" + label.replace("'", "''") + "'"
    return label


def to_newick(node: NewickNode) -> str:
    """Serialize a parsed tree as a semicolon-terminated Newick document.

    :param node: Root node to serialize.
    :return: Newick text with compact numeric branch lengths.
    """
    def serialize(current: NewickNode) -> str:
        """Serialize one node and its descendants.

        :param current: Node currently being encoded.
        :return: Newick fragment for the subtree.
        """
        # Children precede the optional internal label in Newick grammar.
        prefix = "(" + ",".join(serialize(child) for child in current.children) + ")" if current.children else ""
        label = _format_label(current.name) if current.name is not None else ""
        length = f":{current.branch_length:.12g}" if current.branch_length is not None else ""
        return prefix + label + length

    return serialize(node) + ";\n"


def load_calibration(path: str | Path) -> tuple[list[str], float]:
    """Read and validate the structure of a single-point calibration file.

    :param path: Path to the calibration JSON document.
    :return: Calibration taxa and age in millions of years.
    """
    calibration_path = Path(path)
    try:
        data = json.loads(calibration_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid calibration JSON: {error}") from error
    if not isinstance(data, dict) or not isinstance(data.get("node"), dict):
        raise ValueError("Calibration JSON must contain a 'node' object")
    taxa = data["node"].get("taxa")
    if not isinstance(taxa, list) or any(not isinstance(taxon, str) or not taxon for taxon in taxa):
        raise ValueError("Calibration 'node.taxa' must be a list of non-empty strings")
    if "age_mya" not in data:
        raise ValueError("Calibration JSON must contain 'age_mya'")
    return taxa, data["age_mya"]


def calibrate_tree_files(
    input_path: str | Path,
    calibration_path: str | Path,
    output_stem: str | Path,
) -> tuple[Path, Path]:
    """Calibrate a Newick file and write its tree and JSON report.

    :param input_path: Reconstructed Newick tree path.
    :param calibration_path: Single-point calibration JSON path.
    :param output_stem: Path and shared filename stem for generated output files.
    :return: Paths to the calibrated Newick tree and metadata report.
    """
    # Parse both inputs before creating output, avoiding partial results on validation errors.
    tree = read_newick(input_path)
    taxa, age = load_calibration(calibration_path)
    calibrated, metadata = calibrate_tree(tree, taxa, age)
    output_path = Path(output_stem)
    if not output_path.name or output_path.name in {".", ".."}:
        raise ValueError("Output must include a filename stem")
    # Create only the parent directory; the final path component names both files.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree_path = output_path.with_name(f"{output_path.name}.newick")
    report_path = output_path.with_name(f"{output_path.name}_metadata.json")
    tree_path.write_text(to_newick(calibrated), encoding="utf-8")
    report_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return tree_path, report_path
