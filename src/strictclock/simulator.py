"""Core simulation engine for synthetic molecular clock datasets."""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from biology import BiologySettings, load_biology_settings, mutate_base

BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class SimulationConfig:
    sequence_length: int
    number_of_taxa: int
    random_seed: int | None
    tree_topology: str
    root_age: float
    clock_model: str
    mutation_rate: float
    substitution_model: str
    biology: BiologySettings
    raw: dict[str, Any]

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SimulationConfig":
        """Create a validated simulation config from decoded JSON data.

        :param data: Configuration dictionary loaded from a JSON file.
        :return: Validated configuration object ready for simulation.
        """
        # Pull nested sections up front so future config keys can be added cleanly.
        tree = data.get("tree", {})
        clock = data.get("clock", {})
        substitution = data.get("substitution", {})
        biology = data.get("biology", {})

        # Normalize scalar values into the types used by the simulator.
        config = cls(
            sequence_length=int(data.get("sequence_length", 0)),
            number_of_taxa=int(data.get("number_of_taxa", 0)),
            random_seed=data.get("random_seed"),
            tree_topology=str(tree.get("topology", "balanced")),
            root_age=float(tree.get("root_age", 1.0)),
            clock_model=str(clock.get("model", "strict")),
            mutation_rate=float(clock.get("mutation_rate", 0.0)),
            substitution_model=str(substitution.get("model", "simple")),
            biology=BiologySettings.from_dict(biology),
            raw=data,
        )
        config.validate()
        return config

    def validate(self) -> None:
        """Validate Version 1 simulation settings.

        :return: None.
        """
        # Keep validation explicit so unsupported future roadmap settings fail early.
        if self.sequence_length <= 0:
            raise ValueError("sequence_length must be greater than zero")
        if self.number_of_taxa < 2:
            raise ValueError("number_of_taxa must be at least 2")
        if self.random_seed is not None and not isinstance(self.random_seed, int):
            raise ValueError("random_seed must be an integer when provided")
        if self.tree_topology not in {"balanced", "random"}:
            raise ValueError("tree.topology must be 'balanced' or 'random'")
        if self.root_age <= 0:
            raise ValueError("tree.root_age must be greater than zero")
        if self.clock_model != "strict":
            raise ValueError("Version 1 supports only clock.model='strict'")
        if self.mutation_rate < 0:
            raise ValueError("clock.mutation_rate must be non-negative")
        if self.substitution_model != "simple":
            raise ValueError("Version 1 supports only substitution.model='simple'")


@dataclass
class MutationEvent:
    position: int
    ancestral_base: str
    derived_base: str

    def to_dict(self) -> dict[str, Any]:
        """Convert the mutation event into JSON-serializable data.

        :return: Dictionary representation of the mutation event.
        """
        return {
            "position": self.position,
            "ancestral_base": self.ancestral_base,
            "derived_base": self.derived_base,
        }


@dataclass
class TreeNode:
    id: str
    name: str | None = None
    children: list["TreeNode"] = field(default_factory=list)
    age: float = 0.0
    sequence: str = ""
    parent_id: str | None = None
    branch_length: float = 0.0
    mutations_from_parent: list[MutationEvent] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        """Report whether this node is a terminal taxon.

        :return: True when the node has no children, otherwise false.
        """
        return not self.children

    def walk(self) -> list["TreeNode"]:
        """Return this node and all descendants in preorder.

        :return: List containing this node followed by descendant nodes.
        """
        # Preorder traversal keeps parent records before child records in metadata.
        nodes = [self]
        for child in self.children:
            nodes.extend(child.walk())
        return nodes


@dataclass
class SimulationResult:
    config: SimulationConfig
    root: TreeNode
    ancestral_sequence: str
    newick: str

    @property
    def terminal_sequences(self) -> dict[str, str]:
        """Return final DNA sequences for terminal taxa only.

        :return: Mapping of taxon name to simulated DNA sequence.
        """
        # Sort by taxon label to make FASTA and metadata stable across traversals.
        leaves = sorted(
            (node for node in self.root.walk() if node.is_leaf),
            key=lambda n: n.name or "",
        )
        return {node.name or node.id: node.sequence for node in leaves}

    def to_metadata(self) -> dict[str, Any]:
        """Build the complete JSON metadata payload for a simulation run.

        :return: Dictionary containing configuration, tree, sequences, and mutations.
        """
        return {
            "configuration": self.config.raw,
            "random_seed": self.config.random_seed,
            "ancestral_sequence": self.ancestral_sequence,
            "tree_newick": self.newick,
            "mutation_rate": self.config.mutation_rate,
            "nodes": [node_to_dict(node) for node in self.root.walk()],
            "terminal_sequences": self.terminal_sequences,
        }


def load_config(path: str | Path) -> SimulationConfig:
    """Load and validate a JSON simulation configuration file.

    :param path: Path to a JSON configuration file.
    :return: Validated simulation configuration.
    """
    with Path(path).open(encoding="utf-8") as handle:
        data = json.load(handle)

    # The top-level config must be an object so named sections can be accessed safely.
    if not isinstance(data, dict):
        raise ValueError("Configuration JSON must contain an object")
    config = SimulationConfig.from_dict(data)
    biology = load_biology_settings()
    raw = {**config.raw, "biology": biology.to_dict()}

    # Attach shared biology settings without requiring legacy config files to change.
    return replace(config, biology=biology, raw=raw)


def run_simulation(config: SimulationConfig) -> SimulationResult:
    """Generate a complete synthetic molecular clock dataset.

    :param config: Validated simulation configuration.
    :return: Simulation result containing tree, sequences, and metadata helpers.
    """
    # Use a local RNG so identical seeds reproduce runs without touching global state.
    rng = random.Random(config.random_seed)
    ancestral_sequence = random_sequence(config.sequence_length, rng)
    taxa = [f"taxon_{index}" for index in range(1, config.number_of_taxa + 1)]

    # Build the full tree before evolving sequences down from the root.
    id_counter = _IdCounter()
    root = build_tree(taxa, config.tree_topology, config.root_age, rng, id_counter)
    root.sequence = ancestral_sequence
    evolve_children(root, config.mutation_rate, config.biology, rng)

    # Newick export is derived after evolution because branch lengths are already set.
    newick = to_newick(root) + ";"
    return SimulationResult(
        config=config,
        root=root,
        ancestral_sequence=ancestral_sequence,
        newick=newick,
    )


def write_outputs(
    result: SimulationResult,
    config: str,
    output_dir: str | Path
) -> dict[str, Path]:
    """Write FASTA, Newick, and metadata files for a simulation result.

    :param result: Completed simulation result to export.
    :param config: Configuration path whose file stem names the run output folder.
    :param output_dir: Parent directory for generated simulation artifacts.
    :return: Mapping of output type to written file path.
    """
    # Group files by config stem so multiple experiment configs can share a parent dir.
    output_path = Path(output_dir) / Path(config).stem
    output_path.mkdir(parents=True, exist_ok=True)

    fasta_path = output_path / "terminal_sequences.fasta"
    newick_path = output_path / "true_tree.newick"
    metadata_path = output_path / "simulation_metadata.json"

    # Write all outputs as plain text formats that downstream tools can consume.
    fasta_path.write_text(format_fasta(result.terminal_sequences), encoding="utf-8")
    newick_path.write_text(result.newick + "\n", encoding="utf-8")
    metadata_path.write_text(
        json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {"fasta": fasta_path, "newick": newick_path, "metadata": metadata_path}


def random_sequence(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence.

    :param length: Number of nucleotide bases to generate.
    :param rng: Random number generator used for reproducibility.
    :return: DNA sequence containing only A, C, G, and T.
    """
    return "".join(rng.choice(BASES) for _ in range(length))


def build_tree(
    taxa: list[str],
    topology: str,
    root_age: float,
    rng: random.Random,
    id_counter: "_IdCounter",
) -> TreeNode:
    """Build a rooted ultrametric binary tree.

    :param taxa: Terminal taxon labels to place in the tree.
    :param topology: Tree topology strategy, either balanced or random.
    :param root_age: Age assigned to the root node.
    :param rng: Random number generator used for random topologies.
    :param id_counter: Counter used to assign stable node identifiers.
    :return: Root node of the generated tree.
    """
    if topology == "random":
        # Shuffle a copy so callers do not see their original taxon list reordered.
        taxa = taxa[:]
        rng.shuffle(taxa)

    # Build topology first, then assign ages and branch lengths in separate passes.
    root = build_topology(taxa, topology, rng, id_counter)
    root.age = root_age
    assign_internal_ages(root, root_age)
    assign_branch_metadata(root)
    return root


def build_topology(
    taxa: list[str],
    topology: str,
    rng: random.Random,
    id_counter: "_IdCounter",
) -> TreeNode:
    """Recursively construct a binary tree topology.

    :param taxa: Taxa to place under the current subtree.
    :param topology: Tree topology strategy, either balanced or random.
    :param rng: Random number generator used for random split points.
    :param id_counter: Counter used to assign stable node identifiers.
    :return: Root node for the constructed subtree.
    """
    if len(taxa) == 1:
        # A single taxon is a terminal node and ends this recursive branch.
        return TreeNode(id=id_counter.next("leaf"), name=taxa[0])

    if topology == "balanced":
        # Balanced trees split as evenly as possible at every internal node.
        split_at = len(taxa) // 2
    else:
        # Random trees choose any non-empty left and right partition.
        split_at = rng.randint(1, len(taxa) - 1)

    left = build_topology(taxa[:split_at], topology, rng, id_counter)
    right = build_topology(taxa[split_at:], topology, rng, id_counter)
    return TreeNode(id=id_counter.next("node"), children=[left, right])


def assign_internal_ages(root: TreeNode, root_age: float) -> None:
    """Assign node ages so all terminal taxa are the same distance from the root.

    :param root: Root node of the tree.
    :param root_age: Age to assign to the root node.
    :return: None.
    """
    root_height = topology_height(root)

    def assign(node: TreeNode) -> None:
        """Assign an age to one node and recurse into its children.

        :param node: Node receiving an age.
        :return: None.
        """
        if node.is_leaf:
            # Present-day terminal taxa have age zero in this Version 1 model.
            node.age = 0.0
            return

        # Scale ages by subtree height so every root-to-leaf path totals root_age.
        node.age = root_age * topology_height(node) / root_height
        for child in node.children:
            assign(child)

    assign(root)


def topology_height(node: TreeNode) -> int:
    """Measure the number of edges on the longest path to a descendant leaf.

    :param node: Node whose subtree height should be measured.
    :return: Longest descendant path length in edges.
    """
    if node.is_leaf:
        return 0
    return 1 + max(topology_height(child) for child in node.children)


def assign_branch_metadata(node: TreeNode) -> None:
    """Populate parent identifiers and branch lengths throughout the tree.

    :param node: Node whose child branch metadata should be assigned.
    :return: None.
    """
    for child in node.children:
        # Branch length is elapsed time from parent age down to child age.
        child.parent_id = node.id
        child.branch_length = node.age - child.age
        assign_branch_metadata(child)


def evolve_children(
    node: TreeNode,
    mutation_rate: float,
    biology: BiologySettings,
    rng: random.Random,
) -> None:
    """Evolve child sequences from a parent node sequence.

    :param node: Parent node whose descendants should be evolved.
    :param mutation_rate: Strict-clock substitution rate per site per unit time.
    :param biology: Transition/transversion weighting used for substitutions.
    :param rng: Random number generator used for substitutions.
    :return: None.
    """
    for child in node.children:
        # Each child receives an independently mutated copy of its parent sequence.
        child.sequence, child.mutations_from_parent = mutate_sequence(
            node.sequence,
            child.branch_length,
            mutation_rate,
            biology,
            rng,
        )
        evolve_children(child, mutation_rate, biology, rng)


def mutate_sequence(
    sequence: str,
    branch_length: float,
    mutation_rate: float,
    biology: BiologySettings,
    rng: random.Random,
) -> tuple[str, list[MutationEvent]]:
    """Apply simple nucleotide substitutions along one branch.

    :param sequence: Parent sequence at the start of the branch.
    :param branch_length: Evolutionary time represented by the branch.
    :param mutation_rate: Strict-clock substitution rate per site per unit time.
    :param biology: Transition/transversion weighting used for substitutions.
    :param rng: Random number generator used for mutation decisions.
    :return: Mutated sequence and list of mutation events on the branch.
    """
    # Convert rate and time into the chance that a site changes on this branch.
    probability = 1.0 - math.exp(-mutation_rate * branch_length)
    sequence_chars = list(sequence)
    events: list[MutationEvent] = []

    for index, base in enumerate(sequence_chars):
        if rng.random() < probability:
            # The shared helper keeps transition/transversion weighting consistent.
            derived = mutate_base(
                base,
                biology.transition_weight,
                biology.transversion_weight,
                rng,
            )
            sequence_chars[index] = derived
            events.append(
                MutationEvent(
                    position=index,
                    ancestral_base=base,
                    derived_base=derived,
                )
            )

    return "".join(sequence_chars), events


def to_newick(node: TreeNode) -> str:
    """Render a tree node and descendants in Newick format.

    :param node: Node to render.
    :return: Newick fragment for the node subtree, without the final semicolon.
    """
    label = node.name if node.is_leaf else f"internal_{node.id}"
    branch = "" if node.parent_id is None else f":{node.branch_length:.10g}"
    if node.is_leaf:
        return f"{label}{branch}"

    # Internal Newick nodes wrap comma-separated child subtrees in parentheses.
    children = ",".join(to_newick(child) for child in node.children)
    return f"({children}){label}{branch}"


def format_fasta(sequences: dict[str, str], line_width: int = 80) -> str:
    """Format terminal sequences as FASTA text.

    :param sequences: Mapping of sequence name to DNA sequence.
    :param line_width: Maximum number of bases per FASTA sequence line.
    :return: FASTA-formatted text ending with a trailing newline.
    """
    records: list[str] = []
    for name, sequence in sequences.items():
        records.append(f">{name}")
        # Wrap long sequences to keep the FASTA readable and broadly compatible.
        records.extend(
            sequence[index : index + line_width]
            for index in range(0, len(sequence), line_width)
        )
    return "\n".join(records) + "\n"


def node_to_dict(node: TreeNode) -> dict[str, Any]:
    """Convert a tree node into JSON-serializable metadata.

    :param node: Tree node to serialize.
    :return: Dictionary describing node identity, sequence, branches, and mutations.
    """
    return {
        "id": node.id,
        "name": node.name,
        "type": "terminal" if node.is_leaf else "internal",
        "age": node.age,
        "parent_id": node.parent_id,
        "branch_length": node.branch_length,
        "children": [child.id for child in node.children],
        "sequence": node.sequence,
        "mutations_from_parent": [
            mutation.to_dict() for mutation in node.mutations_from_parent
        ],
    }


class _IdCounter:
    def __init__(self) -> None:
        """Create a counter for deterministic node identifiers.

        :return: None.
        """
        self.value = 0

    def next(self, prefix: str) -> str:
        """Return the next identifier with the requested prefix.

        :param prefix: Prefix describing the node type.
        :return: Stable identifier such as node_1 or leaf_2.
        """
        self.value += 1
        return f"{prefix}_{self.value}"
