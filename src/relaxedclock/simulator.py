"""Core simulation engine for relaxed molecular clock datasets."""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass, field, replace
from pathlib import Path
from statistics import mean
from typing import Any

from biology import BiologySettings, is_dna_alphabet, load_biology_settings, mutate_base

@dataclass(frozen=True)
class SimulationSettings:
    name: str
    random_seed: int | None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SimulationSettings":
        """Create simulation-level settings from decoded configuration data.

        :param data: Dictionary containing the simulation configuration section.
        :return: Validated simulation settings.
        """
        # Keep default names human-readable when callers omit the optional field.
        settings = cls(
            name=str(data.get("name", "relaxed-clock-simulation")),
            random_seed=data.get("random_seed"),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate simulation-level settings.

        :return: None.
        """
        # Empty run names make metadata harder to interpret later.
        if not self.name:
            raise ValueError("simulation.name must not be empty")
        if self.random_seed is not None and not isinstance(self.random_seed, int):
            raise ValueError("simulation.random_seed must be an integer when provided")


@dataclass(frozen=True)
class SequenceSettings:
    length: int
    alphabet: tuple[str, ...]
    root_sequence: str | None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SequenceSettings":
        """Create sequence settings from decoded configuration data.

        :param data: Dictionary containing the sequence configuration section.
        :return: Validated sequence settings.
        """
        # Convert the alphabet into an immutable tuple so simulations cannot alter it.
        settings = cls(
            length=int(data.get("length", 0)),
            alphabet=tuple(str(base) for base in data.get("alphabet", [])),
            root_sequence=data.get("root_sequence"),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate sequence generation settings.

        :return: None.
        """
        # Sequence simulation needs at least two symbols so substitutions can change.
        if self.length <= 0:
            raise ValueError("sequence.length must be greater than zero")
        if len(self.alphabet) < 2:
            raise ValueError("sequence.alphabet must contain at least two symbols")
        if len(set(self.alphabet)) != len(self.alphabet):
            raise ValueError("sequence.alphabet must not contain duplicate symbols")
        if any(len(base) != 1 for base in self.alphabet):
            raise ValueError("sequence.alphabet entries must be single characters")
        if self.root_sequence is not None:
            if len(self.root_sequence) != self.length:
                raise ValueError("sequence.root_sequence must match sequence.length")
            if any(base not in self.alphabet for base in self.root_sequence):
                raise ValueError("sequence.root_sequence contains bases outside sequence.alphabet")


@dataclass(frozen=True)
class TreeSettings:
    max_depth: int
    branching_mode: str
    branch_duration: float
    duration_jitter: float

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TreeSettings":
        """Create tree settings from decoded configuration data.

        :param data: Dictionary containing the tree configuration section.
        :return: Validated tree settings.
        """
        # The first version deliberately supports one simple full binary tree shape.
        settings = cls(
            max_depth=int(data.get("max_depth", 0)),
            branching_mode=str(data.get("branching_mode", "binary")),
            branch_duration=float(data.get("branch_duration", 0.0)),
            duration_jitter=float(data.get("duration_jitter", 0.0)),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate tree generation settings.

        :return: None.
        """
        # A max depth of one produces two terminal taxa and is the smallest tree.
        if self.max_depth < 1:
            raise ValueError("tree.max_depth must be at least 1")
        if self.branching_mode != "binary":
            raise ValueError("tree.branching_mode must be 'binary'")
        if self.branch_duration <= 0:
            raise ValueError("tree.branch_duration must be greater than zero")
        if self.duration_jitter < 0:
            raise ValueError("tree.duration_jitter must be non-negative")
        if self.duration_jitter >= self.branch_duration:
            raise ValueError("tree.duration_jitter must be less than tree.branch_duration")


@dataclass(frozen=True)
class ClockSettings:
    model: str
    root_rate: float
    rate_distribution: str
    rate_sigma: float
    minimum_rate: float
    maximum_rate: float

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ClockSettings":
        """Create relaxed clock settings from decoded configuration data.

        :param data: Dictionary containing the clock configuration section.
        :return: Validated relaxed clock settings.
        """
        # Store all rate limits explicitly so generated metadata is self-contained.
        settings = cls(
            model=str(data.get("model", "autocorrelated_relaxed")),
            root_rate=float(data.get("root_rate", 0.0)),
            rate_distribution=str(data.get("rate_distribution", "lognormal")),
            rate_sigma=float(data.get("rate_sigma", 0.0)),
            minimum_rate=float(data.get("minimum_rate", 0.0)),
            maximum_rate=float(data.get("maximum_rate", 0.0)),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate relaxed clock settings.

        :return: None.
        """
        # The initial engine intentionally implements one educational model.
        if self.model != "autocorrelated_relaxed":
            raise ValueError("clock.model must be 'autocorrelated_relaxed'")
        if self.rate_distribution != "lognormal":
            raise ValueError("clock.rate_distribution must be 'lognormal'")
        if self.root_rate <= 0:
            raise ValueError("clock.root_rate must be greater than zero")
        if self.rate_sigma < 0:
            raise ValueError("clock.rate_sigma must be non-negative")
        if self.minimum_rate <= 0:
            raise ValueError("clock.minimum_rate must be greater than zero")
        if self.maximum_rate < self.minimum_rate:
            raise ValueError("clock.maximum_rate must be greater than or equal to clock.minimum_rate")
        if not self.minimum_rate <= self.root_rate <= self.maximum_rate:
            raise ValueError("clock.root_rate must fall between clock.minimum_rate and clock.maximum_rate")


@dataclass(frozen=True)
class MutationSettings:
    model: str
    allow_back_mutation: bool

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MutationSettings":
        """Create mutation settings from decoded configuration data.

        :param data: Dictionary containing the mutation configuration section.
        :return: Validated mutation settings.
        """
        # Default to allowing reversions because that matches common substitution models.
        settings = cls(
            model=str(data.get("model", "simple_substitution")),
            allow_back_mutation=bool(data.get("allow_back_mutation", True)),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate mutation settings.

        :return: None.
        """
        # Keep the first relaxed engine focused on simple substitutions only.
        if self.model != "simple_substitution":
            raise ValueError("mutation.model must be 'simple_substitution'")


@dataclass(frozen=True)
class OutputSettings:
    write_fasta: bool
    write_newick: bool
    write_metadata: bool
    newick_branch_lengths: str

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "OutputSettings":
        """Create output settings from decoded configuration data.

        :param data: Dictionary containing the outputs configuration section.
        :return: Validated output settings.
        """
        # Each output can be toggled independently for lightweight experiment runs.
        settings = cls(
            write_fasta=bool(data.get("write_fasta", True)),
            write_newick=bool(data.get("write_newick", True)),
            write_metadata=bool(data.get("write_metadata", True)),
            newick_branch_lengths=str(data.get("newick_branch_lengths", "genetic_change")),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate output settings.

        :return: None.
        """
        # Newick branch lengths expose the main educational contrast in this engine.
        if self.newick_branch_lengths not in {"time", "rate", "genetic_change"}:
            raise ValueError(
                "outputs.newick_branch_lengths must be 'time', 'rate', or 'genetic_change'"
            )


@dataclass(frozen=True)
class RelaxedClockConfig:
    simulation: SimulationSettings
    sequence: SequenceSettings
    tree: TreeSettings
    clock: ClockSettings
    mutation: MutationSettings
    outputs: OutputSettings
    biology: BiologySettings
    raw: dict[str, Any]

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RelaxedClockConfig":
        """Create a validated relaxed clock config from decoded JSON data.

        :param data: Configuration dictionary loaded from a JSON file.
        :return: Validated relaxed clock configuration object.
        """
        # Reject missing or mismatched files before interpreting simulator-specific sections.
        intended_clock_model = data.get("clock_model")
        if intended_clock_model != "relaxed":
            raise ValueError(
                "clock_model must be 'relaxed' for the relaxed clock simulator "
                f"(received {intended_clock_model!r})"
            )

        # Validate each required section separately so error messages name the section.
        required_sections = ("simulation", "sequence", "tree", "clock", "mutation", "outputs")
        for section in required_sections:
            if section not in data:
                raise ValueError(f"Missing required configuration section: {section}")
            if not isinstance(data[section], dict):
                raise ValueError(f"Configuration section {section} must be an object")

        return cls(
            simulation=SimulationSettings.from_dict(data["simulation"]),
            sequence=SequenceSettings.from_dict(data["sequence"]),
            tree=TreeSettings.from_dict(data["tree"]),
            clock=ClockSettings.from_dict(data["clock"]),
            mutation=MutationSettings.from_dict(data["mutation"]),
            outputs=OutputSettings.from_dict(data["outputs"]),
            biology=BiologySettings.from_dict(data.get("biology", {})),
            raw=data,
        )


@dataclass
class RelaxedMutationEvent:
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
class RelaxedTreeNode:
    id: str
    name: str | None = None
    children: list["RelaxedTreeNode"] = field(default_factory=list)
    parent_id: str | None = None
    depth: int = 0
    sequence: str = ""
    branch_duration: float = 0.0
    lineage_rate: float = 0.0
    genetic_change: float = 0.0
    expected_substitutions: float = 0.0
    mutations_from_parent: list[RelaxedMutationEvent] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        """Report whether this node is a terminal taxon.

        :return: True when the node has no children, otherwise false.
        """
        return not self.children

    @property
    def observed_substitutions(self) -> int:
        """Count observed substitutions on the branch from the parent.

        :return: Number of mutation events recorded for this node's incoming branch.
        """
        return len(self.mutations_from_parent)

    def walk(self) -> list["RelaxedTreeNode"]:
        """Return this node and all descendants in preorder.

        :return: List containing this node followed by descendant nodes.
        """
        # Preorder traversal keeps parent records before child records in metadata.
        nodes = [self]
        for child in self.children:
            nodes.extend(child.walk())
        return nodes


@dataclass
class RelaxedClockResult:
    config: RelaxedClockConfig
    root: RelaxedTreeNode
    root_sequence: str
    newick: str

    @property
    def terminal_sequences(self) -> dict[str, str]:
        """Return final sequences for terminal taxa only.

        :return: Mapping of taxon name to simulated sequence.
        """
        # Sort labels so FASTA and metadata remain deterministic.
        leaves = sorted(
            (node for node in self.root.walk() if node.is_leaf),
            key=lambda node: node.name or "",
        )
        return {node.name or node.id: node.sequence for node in leaves}

    def to_metadata(self) -> dict[str, Any]:
        """Build the complete JSON metadata payload for a relaxed clock run.

        :return: Dictionary containing configuration, tree, sequences, and summaries.
        """
        # Metadata keeps both time and genetic branch values for later analysis.
        return {
            "configuration": self.config.raw,
            "simulation": {
                "name": self.config.simulation.name,
                "random_seed": self.config.simulation.random_seed,
            },
            "root_sequence": self.root_sequence,
            "generated_taxa": list(self.terminal_sequences),
            "tree_newick": self.newick,
            "nodes": [node_to_dict(node) for node in self.root.walk()],
            "terminal_sequences": self.terminal_sequences,
            "summary": summary_statistics(self),
        }


def load_config(path: str | Path) -> RelaxedClockConfig:
    """Load and validate a relaxed clock JSON configuration file.

    :param path: Path to a JSON configuration file.
    :return: Validated relaxed clock configuration.
    """
    with Path(path).open(encoding="utf-8") as handle:
        data = json.load(handle)

    # The top-level config must be an object so sections can be validated by name.
    if not isinstance(data, dict):
        raise ValueError("Configuration JSON must contain an object")
    config = RelaxedClockConfig.from_dict(data)
    biology = load_biology_settings()
    raw = {**config.raw, "biology": biology.to_dict()}

    # Attach shared biology settings without changing existing relaxed config files.
    return replace(config, biology=biology, raw=raw)


def run_simulation(config: RelaxedClockConfig) -> RelaxedClockResult:
    """Generate a complete relaxed molecular clock simulation result.

    :param config: Validated relaxed clock configuration.
    :return: Simulation result containing tree, terminal sequences, and metadata.
    """
    # Use a local RNG so seeded runs reproduce without changing global random state.
    rng = random.Random(config.simulation.random_seed)
    root_sequence = config.sequence.root_sequence or random_sequence(
        config.sequence.length,
        config.sequence.alphabet,
        rng,
    )

    # Build topology first, then assign branch durations, inherited rates, and sequences.
    id_counter = _IdCounter()
    taxon_counter = _IdCounter()
    root = build_time_tree(config.tree.max_depth, id_counter, taxon_counter)
    root.sequence = root_sequence
    root.lineage_rate = config.clock.root_rate
    assign_branch_values(root, config, rng)
    evolve_children(root, config, root_sequence, rng)

    # Newick rendering is derived last because branch values depend on evolution settings.
    newick = to_newick(root, config.outputs.newick_branch_lengths) + ";"
    return RelaxedClockResult(
        config=config,
        root=root,
        root_sequence=root_sequence,
        newick=newick,
    )


def write_outputs(
    result: RelaxedClockResult,
    config: str | Path,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write configured FASTA, Newick, and metadata files for a relaxed run.

    :param result: Completed relaxed clock simulation result to export.
    :param config: Configuration path whose file stem names the run output folder.
    :param output_dir: Parent directory for generated relaxed clock artifacts.
    :return: Mapping of output type to written file path.
    """
    # Group files by config stem so experiment configs can share a parent directory.
    output_path = Path(output_dir) / Path(config).stem
    output_path.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    if result.config.outputs.write_fasta:
        fasta_path = output_path / "terminal_sequences.fasta"
        fasta_path.write_text(format_fasta(result.terminal_sequences), encoding="utf-8")
        written["FASTA"] = fasta_path

    if result.config.outputs.write_newick:
        newick_path = output_path / "true_tree.newick"
        newick_path.write_text(result.newick + "\n", encoding="utf-8")
        written["Newick"] = newick_path

    if result.config.outputs.write_metadata:
        metadata_path = output_path / "simulation_metadata.json"
        metadata_path.write_text(
            json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        written["metadata"] = metadata_path

    return written


def random_sequence(length: int, alphabet: tuple[str, ...], rng: random.Random) -> str:
    """Generate a random root sequence from the configured alphabet.

    :param length: Number of symbols to generate.
    :param alphabet: Allowed sequence symbols.
    :param rng: Random number generator used for reproducibility.
    :return: Random sequence string.
    """
    # Draw each site independently to make the root sequence easy to reason about.
    return "".join(rng.choice(alphabet) for _ in range(length))


def build_time_tree(
    max_depth: int,
    id_counter: "_IdCounter",
    taxon_counter: "_IdCounter",
    depth: int = 0,
) -> RelaxedTreeNode:
    """Build a full binary time tree with terminal taxa at the same depth.

    :param max_depth: Number of branching levels from root to terminal taxa.
    :param id_counter: Counter used to assign stable node identifiers.
    :param taxon_counter: Counter used to assign terminal taxon labels.
    :param depth: Current depth during recursive construction.
    :return: Root node for the constructed subtree.
    """
    # Nodes at max depth are present-day terminal taxa.
    if depth == max_depth:
        taxon_id = taxon_counter.next("Taxon")
        return RelaxedTreeNode(
            id=id_counter.next("leaf"),
            name=taxon_id,
            depth=depth,
        )

    # Internal nodes always split into two daughters in the initial binary model.
    node = RelaxedTreeNode(id=id_counter.next("node"), depth=depth)
    node.children = [
        build_time_tree(max_depth, id_counter, taxon_counter, depth + 1),
        build_time_tree(max_depth, id_counter, taxon_counter, depth + 1),
    ]
    return node


def assign_branch_values(
    node: RelaxedTreeNode,
    config: RelaxedClockConfig,
    rng: random.Random,
) -> None:
    """Assign parent links, branch durations, lineage rates, and genetic lengths.

    :param node: Parent node whose child branch values should be assigned.
    :param config: Relaxed clock configuration controlling branch values.
    :param rng: Random number generator used for jitter and rate modifiers.
    :return: None.
    """
    for child in node.children:
        # Duration jitter can break time ultrametricity when intentionally configured.
        duration = jittered_duration(
            config.tree.branch_duration,
            config.tree.duration_jitter,
            rng,
        )
        rate = inherited_rate(node.lineage_rate, config.clock, rng)

        # Genetic change is the relaxed-clock branch length seen by downstream tools.
        child.parent_id = node.id
        child.branch_duration = duration
        child.lineage_rate = rate
        child.genetic_change = duration * rate
        child.expected_substitutions = config.sequence.length * child.genetic_change
        assign_branch_values(child, config, rng)


def jittered_duration(base_duration: float, jitter: float, rng: random.Random) -> float:
    """Draw a branch duration around the configured base duration.

    :param base_duration: Base elapsed time for a branch.
    :param jitter: Maximum absolute uniform perturbation around the base duration.
    :param rng: Random number generator used when jitter is non-zero.
    :return: Positive branch duration.
    """
    # Zero jitter preserves a time-ultrametric full binary tree exactly.
    if jitter == 0:
        return base_duration
    return base_duration + rng.uniform(-jitter, jitter)


def inherited_rate(
    parent_rate: float,
    clock: ClockSettings,
    rng: random.Random,
) -> float:
    """Draw an autocorrelated daughter lineage rate from a parent rate.

    :param parent_rate: Rate inherited from the parent lineage.
    :param clock: Relaxed clock settings including sigma and rate bounds.
    :param rng: Random number generator used for the rate modifier.
    :return: Daughter lineage rate after applying configured bounds.
    """
    # The lognormal mean is shifted so the modifier is centered near one.
    modifier = rng.lognormvariate(-0.5 * clock.rate_sigma**2, clock.rate_sigma)
    return clamp(parent_rate * modifier, clock.minimum_rate, clock.maximum_rate)


def clamp(value: float, lower: float, upper: float) -> float:
    """Constrain a numeric value to inclusive lower and upper bounds.

    :param value: Candidate numeric value.
    :param lower: Minimum allowed value.
    :param upper: Maximum allowed value.
    :return: Value constrained to the provided interval.
    """
    # Apply the lower bound first, then the upper bound.
    return min(max(value, lower), upper)


def evolve_children(
    node: RelaxedTreeNode,
    config: RelaxedClockConfig,
    root_sequence: str,
    rng: random.Random,
) -> None:
    """Evolve child sequences from a parent node sequence.

    :param node: Parent node whose descendants should be evolved.
    :param config: Relaxed clock configuration controlling substitutions.
    :param root_sequence: Original root sequence used when preventing back mutation.
    :param rng: Random number generator used for substitutions.
    :return: None.
    """
    for child in node.children:
        # Each child receives an independently mutated copy of its parent sequence.
        child.sequence, child.mutations_from_parent = mutate_sequence(
            node.sequence,
            root_sequence,
            child.branch_duration,
            child.lineage_rate,
            config.sequence.alphabet,
            config.mutation.allow_back_mutation,
            config.biology,
            rng,
        )
        evolve_children(child, config, root_sequence, rng)


def mutate_sequence(
    sequence: str,
    root_sequence: str,
    branch_duration: float,
    lineage_rate: float,
    alphabet: tuple[str, ...],
    allow_back_mutation: bool,
    biology: BiologySettings,
    rng: random.Random,
) -> tuple[str, list[RelaxedMutationEvent]]:
    """Apply simple substitutions along one relaxed-clock branch.

    :param sequence: Parent sequence at the start of the branch.
    :param root_sequence: Original root sequence used to avoid reversions when requested.
    :param branch_duration: Elapsed time represented by the branch.
    :param lineage_rate: Branch-specific substitution rate.
    :param alphabet: Allowed sequence symbols.
    :param allow_back_mutation: Whether a site may revert to its root state.
    :param biology: Biology settings used for DNA substitutions.
    :param rng: Random number generator used for mutation decisions.
    :return: Mutated sequence and list of mutation events on the branch.
    """
    # Convert branch duration and lineage rate into a per-site substitution probability.
    probability = 1.0 - math.exp(-lineage_rate * branch_duration)
    sequence_chars = list(sequence)
    events: list[RelaxedMutationEvent] = []

    for index, base in enumerate(sequence_chars):
        if rng.random() < probability:
            # Substitutions must change the current base; optional back mutation excludes root.
            excluded = {base}
            if not allow_back_mutation:
                excluded.add(root_sequence[index])
            alternatives = [candidate for candidate in alphabet if candidate not in excluded]
            if not alternatives:
                continue
            if is_dna_alphabet(alphabet):
                # DNA alphabets use the shared exchangeability and frequency-aware helper.
                derived = mutate_base(
                    base,
                    biology.transition_weight,
                    biology.transversion_weight,
                    rng,
                    candidates=alternatives,
                    equilibrium_frequencies=biology.equilibrium_frequencies,
                    exchangeability_rates=biology.exchangeability_rates,
                )
            else:
                # Non-DNA alphabets retain the historical uniform substitution behavior.
                derived = rng.choice(alternatives)
            sequence_chars[index] = derived
            events.append(
                RelaxedMutationEvent(
                    position=index,
                    ancestral_base=base,
                    derived_base=derived,
                )
            )

    return "".join(sequence_chars), events


def to_newick(node: RelaxedTreeNode, branch_lengths: str) -> str:
    """Render a tree node and descendants in Newick format.

    :param node: Node to render.
    :param branch_lengths: Branch value to use, either time, rate, or genetic_change.
    :return: Newick fragment for the node subtree, without the final semicolon.
    """
    label = node.name if node.is_leaf else f"internal_{node.id}"
    branch = "" if node.parent_id is None else f":{newick_branch_value(node, branch_lengths):.10g}"
    if node.is_leaf:
        return f"{label}{branch}"

    # Internal Newick nodes wrap comma-separated child subtrees in parentheses.
    children = ",".join(to_newick(child, branch_lengths) for child in node.children)
    return f"({children}){label}{branch}"


def newick_branch_value(node: RelaxedTreeNode, branch_lengths: str) -> float:
    """Return the configured branch length value for Newick export.

    :param node: Node whose incoming branch should be measured.
    :param branch_lengths: Branch value type, either time, rate, or genetic_change.
    :return: Numeric branch length for Newick output.
    """
    # Keep this mapping explicit so unsupported values fail during config validation.
    if branch_lengths == "time":
        return node.branch_duration
    if branch_lengths == "rate":
        return node.lineage_rate
    return node.genetic_change


def format_fasta(sequences: dict[str, str], line_width: int = 80) -> str:
    """Format terminal sequences as FASTA text.

    :param sequences: Mapping of sequence name to sequence string.
    :param line_width: Maximum number of characters per FASTA sequence line.
    :return: FASTA-formatted text ending with a trailing newline.
    """
    records: list[str] = []
    for name, sequence in sequences.items():
        records.append(f">{name}")
        # Wrap long sequences to keep the FASTA readable for downstream tools.
        records.extend(
            sequence[index : index + line_width]
            for index in range(0, len(sequence), line_width)
        )
    return "\n".join(records) + "\n"


def node_to_dict(node: RelaxedTreeNode) -> dict[str, Any]:
    """Convert a relaxed tree node into JSON-serializable metadata.

    :param node: Tree node to serialize.
    :return: Dictionary describing node identity, branch values, sequence, and mutations.
    """
    return {
        "id": node.id,
        "name": node.name,
        "parent": node.parent_id,
        "children": [child.id for child in node.children],
        "depth": node.depth,
        "is_terminal": node.is_leaf,
        "branch_duration": node.branch_duration,
        "lineage_rate": node.lineage_rate,
        "genetic_change": node.genetic_change,
        "expected_substitutions": node.expected_substitutions,
        "observed_substitutions": node.observed_substitutions,
        "sequence": node.sequence,
        "mutations_from_parent": [
            mutation.to_dict() for mutation in node.mutations_from_parent
        ],
    }


def summary_statistics(result: RelaxedClockResult) -> dict[str, Any]:
    """Calculate summary statistics for relaxed clock metadata.

    :param result: Completed relaxed clock simulation result.
    :return: Dictionary of summary statistics for inspection and comparison.
    """
    # Exclude the root from branch-specific rate and substitution summaries.
    branch_nodes = [node for node in result.root.walk() if node.parent_id is not None]
    rates = [node.lineage_rate for node in branch_nodes]
    observed = sum(node.observed_substitutions for node in branch_nodes)
    expected = sum(node.expected_substitutions for node in branch_nodes)

    return {
        "terminal_taxa": len(result.terminal_sequences),
        "sequence_length": result.config.sequence.length,
        "root_rate": result.config.clock.root_rate,
        "minimum_lineage_rate": min(rates) if rates else result.config.clock.root_rate,
        "maximum_lineage_rate": max(rates) if rates else result.config.clock.root_rate,
        "mean_lineage_rate": mean(rates) if rates else result.config.clock.root_rate,
        "total_expected_substitutions": expected,
        "total_observed_substitutions": observed,
        "newick_branch_lengths": result.config.outputs.newick_branch_lengths,
    }


class _IdCounter:
    def __init__(self) -> None:
        """Create a counter for deterministic identifiers.

        :return: None.
        """
        self.value = 0

    def next(self, prefix: str) -> str:
        """Return the next identifier with the requested prefix.

        :param prefix: Prefix describing the identifier type.
        :return: Stable identifier such as node_1 or Taxon_2.
        """
        self.value += 1
        return f"{prefix}_{self.value}"
