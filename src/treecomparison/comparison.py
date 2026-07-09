"""Parse Newick trees and render visual side-by-side comparisons."""

from __future__ import annotations

from dataclasses import dataclass, field
import html
from io import BytesIO
import math
from pathlib import Path
import shutil
import subprocess


@dataclass
class NewickNode:
    """Represent one parsed Newick tree node.

    :param name: Optional leaf or internal-node label.
    :param branch_length: Optional branch length from the parent to this node.
    :param children: Child nodes below this node.
    """

    name: str | None = None
    branch_length: float | None = None
    children: list["NewickNode"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        """Report whether this node has no descendants.

        :return: True for a leaf node, otherwise false.
        """
        # Newick leaves are exactly the nodes without a parenthesized child list.
        return not self.children

    def walk(self) -> list["NewickNode"]:
        """Return this node and all descendants in preorder.

        :return: Nodes ordered with each parent before its descendants.
        """
        # Preorder gives stable node identifiers when generating Graphviz source.
        nodes = [self]
        for child in self.children:
            nodes.extend(child.walk())
        return nodes


class _NewickParser:
    """Maintain parser state while recursively reading one Newick tree.

    :param text: Complete Newick document.
    """

    def __init__(self, text: str) -> None:
        """Initialize parsing at the first character.

        :param text: Complete Newick document.
        :return: None.
        """
        # Keep one shared cursor so recursive calls consume the same input stream.
        self.text = text
        self.position = 0

    def parse(self) -> NewickNode:
        """Parse the complete Newick document.

        :return: Root node of the parsed tree.
        """
        # A document contains exactly one subtree followed by a required semicolon.
        self._skip_whitespace()
        root = self._parse_subtree()
        self._skip_whitespace()
        self._expect(";")
        self._skip_whitespace()
        if self.position != len(self.text):
            raise ValueError(
                f"Unexpected content after Newick tree at character {self.position + 1}"
            )
        return root

    def _parse_subtree(self) -> NewickNode:
        """Parse one leaf or parenthesized internal subtree.

        :return: Parsed node at the root of this subtree.
        """
        self._skip_whitespace()
        if self._peek() == "(":
            # Internal nodes begin with two or more comma-separated child subtrees.
            self.position += 1
            children = [self._parse_subtree()]
            self._skip_whitespace()
            while self._peek() == ",":
                self.position += 1
                children.append(self._parse_subtree())
                self._skip_whitespace()
            self._expect(")")
            if len(children) < 2:
                raise ValueError("Newick internal nodes must contain at least two children")
            name = self._parse_optional_label()
            branch_length = self._parse_optional_branch_length()
            return NewickNode(name=name, branch_length=branch_length, children=children)

        # Leaves require a label so they can be identified in the comparison.
        name = self._parse_optional_label()
        if name is None:
            raise ValueError(f"Expected a leaf label at character {self.position + 1}")
        return NewickNode(
            name=name,
            branch_length=self._parse_optional_branch_length(),
        )

    def _parse_optional_label(self) -> str | None:
        """Parse an optional quoted or unquoted node label.

        :return: Decoded label, or None when no label is present.
        """
        self._skip_whitespace()
        if self._peek() == "'":
            # Newick escapes an apostrophe inside a quoted label by doubling it.
            self.position += 1
            characters: list[str] = []
            while self.position < len(self.text):
                character = self.text[self.position]
                self.position += 1
                if character != "'":
                    characters.append(character)
                    continue
                if self._peek() == "'":
                    self.position += 1
                    characters.append("'")
                    continue
                return "".join(characters)
            raise ValueError("Unterminated quoted label in Newick tree")

        start = self.position
        while self.position < len(self.text):
            character = self.text[self.position]
            if character in "(),:;" or character.isspace():
                break
            self.position += 1
        return self.text[start:self.position] or None

    def _parse_optional_branch_length(self) -> float | None:
        """Parse an optional numeric branch length.

        :return: Branch length, or None when the node has no length.
        """
        self._skip_whitespace()
        if self._peek() != ":":
            return None
        self.position += 1
        self._skip_whitespace()
        start = self.position
        while self.position < len(self.text):
            character = self.text[self.position]
            if character in ",();" or character.isspace():
                break
            self.position += 1
        raw_length = self.text[start:self.position]
        try:
            branch_length = float(raw_length)
        except ValueError as error:
            raise ValueError(f"Invalid Newick branch length: {raw_length!r}") from error
        if not math.isfinite(branch_length):
            raise ValueError("Newick branch lengths must be finite")
        return branch_length

    def _skip_whitespace(self) -> None:
        """Advance past whitespace between Newick tokens.

        :return: None.
        """
        # Whitespace outside labels has no structural meaning in Newick.
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def _peek(self) -> str | None:
        """Inspect the current character without consuming it.

        :return: Current character, or None at the end of input.
        """
        # Bounds checking here keeps the recursive parser methods concise.
        if self.position >= len(self.text):
            return None
        return self.text[self.position]

    def _expect(self, expected: str) -> None:
        """Consume one required structural character.

        :param expected: Character required at the current parser position.
        :return: None.
        """
        self._skip_whitespace()
        # Report a one-based location to make malformed input easier to fix.
        if self._peek() != expected:
            found = self._peek()
            raise ValueError(
                f"Expected {expected!r} at character {self.position + 1}, found {found!r}"
            )
        self.position += 1


def parse_newick(text: str) -> NewickNode:
    """Parse one Newick document into a reusable node tree.

    :param text: Newick text containing one semicolon-terminated tree.
    :return: Root node of the parsed tree.
    """
    # Reject empty input explicitly instead of reporting an obscure missing-label error.
    if not text.strip():
        raise ValueError("Newick input must not be empty")
    return _NewickParser(text).parse()


def read_newick(path: str | Path) -> NewickNode:
    """Read and parse a Newick tree file.

    :param path: Path to a Newick tree file.
    :return: Root node of the parsed tree.
    """
    input_path = Path(path)
    if not input_path.exists():
        raise ValueError(f"Newick file does not exist: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"Newick path is not a file: {input_path}")
    # UTF-8 covers ordinary taxon labels and preserves quoted Unicode labels.
    return parse_newick(input_path.read_text(encoding="utf-8"))


def tree_to_dot(root: NewickNode, title: str, subtitle: str | None = None) -> str:
    """Convert a parsed tree to left-to-right Graphviz DOT.

    :param root: Root node of the parsed Newick tree.
    :param title: Heading displayed above the tree.
    :param subtitle: Optional smaller heading displayed below the title.
    :return: Graphviz DOT source for the tree.
    """
    nodes = root.walk()
    identifiers = {id(node): f"node_{index}" for index, node in enumerate(nodes)}
    # HTML-like labels allow the filename subtitle to use a smaller font.
    graph_label = (
        f'<<FONT POINT-SIZE="18">{html.escape(title)}</FONT>'
        + (
            f'<BR/><FONT POINT-SIZE="11">{html.escape(subtitle)}</FONT>>'
            if subtitle is not None
            else ">"
        )
    )
    # rankdir=LR places the root on the left and descendant leaves on the right.
    lines = [
        "digraph newick_tree {",
        '  graph [rankdir=LR, bgcolor="white", margin=0.12, pad=0.15, '
        f"label={graph_label}, labelloc=t, fontname=\"Helvetica\"];",
        '  node [shape=circle, width=0.13, height=0.13, fixedsize=true, '
        'label="", color="#334155", fillcolor="#334155", style=filled];',
        '  edge [color="#64748b", fontcolor="#475569", fontsize=9, '
        'fontname="Helvetica"];',
    ]

    for node in nodes:
        node_id = identifiers[id(node)]
        if node.is_leaf:
            # Leaf labels sit to the right of a small endpoint marker.
            label = _dot_escape(node.name or "")
            lines.append(
                f'  "{node_id}" [shape=plaintext, fixedsize=false, '
                f'label="{label}", fontname="Helvetica", fontsize=11, fontcolor="white"];'
            )
        for child in node.children:
            child_id = identifiers[id(child)]
            length_label = (
                ""
                if child.branch_length is None
                else f' [label="{child.branch_length:g}"]'
            )
            lines.append(f'  "{node_id}" -> "{child_id}"{length_label};')
    lines.append("}")
    return "\n".join(lines)


def render_tree_png(
    root: NewickNode,
    title: str,
    subtitle: str | None = None,
) -> bytes:
    """Render one parsed tree as PNG bytes using Graphviz.

    :param root: Root node of the parsed Newick tree.
    :param title: Heading displayed above the tree.
    :param subtitle: Optional smaller heading displayed below the title.
    :return: PNG image bytes.
    """
    if shutil.which("dot") is None:
        raise RuntimeError("Graphviz 'dot' command is required for tree comparison")
    # Pass DOT through standard input so no intermediate files need cleanup.
    process = subprocess.run(
        ["dot", "-Tpng"],
        input=tree_to_dot(root, title, subtitle).encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        message = process.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"Graphviz could not render the tree: {message}")
    return process.stdout


def compare_trees(
    source_path: str | Path,
    reconstructed_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Render source and reconstructed Newick trees side by side.

    :param source_path: Newick source tree shown in the left panel.
    :param reconstructed_path: Newick reconstructed tree shown in the right panel.
    :param output_path: Destination PNG file path.
    :return: Destination path after writing the comparison image.
    """
    from PIL import Image

    destination = Path(output_path)
    if destination.suffix.lower() != ".png":
        raise ValueError("Tree comparison output path must end with .png")

    # Render independently so horizontal composition always preserves panel ordering.
    source_png = render_tree_png(
        read_newick(source_path),
        "Source tree",
        Path(source_path).stem,
    )
    reconstructed_png = render_tree_png(
        read_newick(reconstructed_path),
        "Reconstructed tree",
        Path(reconstructed_path).stem,
    )
    with Image.open(BytesIO(source_png)) as source_image:
        with Image.open(BytesIO(reconstructed_png)) as reconstructed_image:
            # Flatten Graphviz transparency onto white so only intentional borders are dark.
            source = Image.new("RGB", source_image.size, "white")
            source.paste(source_image, mask=source_image.getchannel("A"))
            reconstructed = Image.new("RGB", reconstructed_image.size, "white")
            reconstructed.paste(
                reconstructed_image,
                mask=reconstructed_image.getchannel("A"),
            )
            border_width = 6
            panel_height = max(source.height, reconstructed.height)
            canvas = Image.new(
                "RGB",
                (
                    source.width + reconstructed.width + (3 * border_width),
                    panel_height + (2 * border_width),
                ),
                "black",
            )
            # White equal-height panels leave one thin divider and a matching outer frame.
            source_x = border_width
            reconstructed_x = source_x + source.width + border_width
            canvas.paste(
                "white",
                (source_x, border_width, source_x + source.width, border_width + panel_height),
            )
            canvas.paste(
                "white",
                (
                    reconstructed_x,
                    border_width,
                    reconstructed_x + reconstructed.width,
                    border_width + panel_height,
                ),
            )
            canvas.paste(source, (source_x, border_width))
            canvas.paste(reconstructed, (reconstructed_x, border_width))
            destination.parent.mkdir(parents=True, exist_ok=True)
            canvas.save(destination, format="PNG")
    return destination


def _dot_escape(value: str) -> str:
    """Escape text embedded in a quoted Graphviz string.

    :param value: Raw label or title.
    :return: Graphviz-safe string content.
    """
    # Escape backslashes first so quote escaping is not itself altered.
    return value.replace("\\", "\\\\").replace('"', '\\"')
