"""Render `docs/architecture/website-flow.drawio` to its SVG preview.

The two files describe the same diagram, and until now both were maintained by hand. That is
the arrangement that let the promoter map and its legend drift apart until four of five
swatches were wrong, and there is no reason to expect a diagram to behave better. The `.drawio`
is the source; this produces the `.svg` from it, so the preview cannot say something the
diagram does not.

Deliberately not a general draw.io renderer. It understands the subset this one file uses:
rounded rectangles, text labels, and orthogonal edges between them. An element it does not
recognise makes it fail rather than silently omit something.

Run after editing the diagram::

    .venv/bin/python scripts/render_flow_diagram.py
"""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree

PAGE_GAP = 60
"""Blank space between the two stacked pages."""

ARROW = (
    '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
    'markerHeight="7" orient="auto-start-reverse">'
    '<path d="M0,0 L10,5 L0,10 z" fill="#748b87" /></marker>'
    '<marker id="arrow-remote" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
    'markerHeight="7" orient="auto-start-reverse">'
    '<path d="M0,0 L10,5 L0,10 z" fill="#c08a4e" /></marker></defs>'
)


@dataclass(frozen=True, slots=True)
class Box:
    """One rectangle, with the text inside it."""

    ident: str
    lines: tuple[str, ...]
    x: int
    y: int
    width: int
    height: int
    fill: str
    stroke: str
    ink: str
    size: int
    is_text: bool

    @property
    def centre_y(self) -> float:
        """Vertical centre, used to route edges."""
        return self.y + self.height / 2

    @property
    def centre_x(self) -> float:
        """Horizontal centre."""
        return self.x + self.width / 2


def style_value(style: str, key: str, default: str) -> str:
    """Return one `key=value;` from a draw.io style string."""
    found = re.search(rf"(?:^|;){re.escape(key)}=([^;]*)", style)
    return found.group(1) if found else default


def escape(text: str) -> str:
    """Escape text for XML content."""
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def read_pages(path: Path) -> list[tuple[str, list[Box], list[tuple[str, str, str, str]]]]:
    """Return each diagram page as its name, its boxes and its edges.

    Raises:
        ValueError: If a cell is neither a recognised vertex nor an edge between two of them.
    """
    tree = ElementTree.parse(path)
    pages = []
    for diagram in tree.getroot().iter("diagram"):
        boxes: list[Box] = []
        edges: list[tuple[str, str, str, str]] = []
        for cell in diagram.iter("mxCell"):
            style = cell.get("style") or ""
            if cell.get("vertex") == "1":
                geometry = cell.find("mxGeometry")
                if geometry is None:
                    raise ValueError(f"vertex {cell.get('id')} has no geometry")
                is_text = style.startswith("text;")
                boxes.append(Box(
                    ident=cell.get("id") or "",
                    lines=tuple((cell.get("value") or "").split("\n")),
                    x=int(float(geometry.get("x") or 0)),
                    y=int(float(geometry.get("y") or 0)),
                    width=int(float(geometry.get("width") or 0)),
                    height=int(float(geometry.get("height") or 0)),
                    fill=style_value(style, "fillColor", "none"),
                    stroke=style_value(style, "strokeColor", "none"),
                    ink=style_value(style, "fontColor", "#183932"),
                    size=int(style_value(style, "fontSize", "14")),
                    is_text=is_text,
                ))
            elif cell.get("edge") == "1":
                source, target = cell.get("source"), cell.get("target")
                if not source or not target:
                    raise ValueError(f"edge {cell.get('id')} is not between two cells")
                edges.append((source, target, cell.get("value") or "",
                              style_value(style, "strokeColor", "#748b87")))
        pages.append((diagram.get("name") or "", boxes, edges))
    return pages


def route(source: Box, target: Box) -> tuple[list[tuple[float, float]], tuple[float, float]]:
    """Return the polyline between two boxes, and where to put the label.

    Three cases, which is all this diagram's grid produces: side by side, stacked, or an elbow
    between the two. The elbow turns at the midpoint of the vertical gap, as a reader expects.
    """
    if abs(source.centre_y - target.centre_y) < 1:  # same row
        forward = target.x > source.x
        start = ((float(source.x + source.width), source.centre_y) if forward
                 else (float(source.x), source.centre_y))
        end = ((float(target.x), target.centre_y) if forward
               else (float(target.x + target.width), target.centre_y))
        return [start, end], ((start[0] + end[0]) / 2, start[1] - 10)
    if abs(source.centre_x - target.centre_x) < 1:  # same column
        down = target.y > source.y
        start = (source.centre_x, float(source.y + source.height if down else source.y))
        end = (target.centre_x, float(target.y if down else target.y + target.height))
        return [start, end], (start[0] + 8, (start[1] + end[1]) / 2)
    down = target.centre_y > source.centre_y
    start = (source.centre_x, float(source.y + source.height if down else source.y))
    end = (target.centre_x, float(target.y if down else target.y + target.height))
    turn = (start[1] + end[1]) / 2
    points = [start, (start[0], turn), (end[0], turn), end]
    return points, ((start[0] + end[0]) / 2, turn - 8)


def wrap(lines: tuple[str, ...], width: int, size: int) -> list[str]:
    """Break text to fit a declared width.

    SVG does not wrap text, so a caption written as one long string runs off the page and is
    simply lost. draw.io wraps it for you, which is how the two renderings diverged.
    """
    # 0.52 em per character is close enough for this sans-serif at these sizes; it only has to
    # keep the line inside the box.
    budget = max(10, int(width / (size * 0.52)))
    out: list[str] = []
    for line in lines:
        words, current = line.split(), ""
        for word in words:
            candidate = f"{current} {word}".strip()
            if len(candidate) > budget and current:
                out.append(current)
                current = word
            else:
                current = candidate
        out.append(current)
    return out


def render_page(boxes: list[Box], edges: list[tuple[str, str, str, str]], offset: int) -> str:
    """Return one page's SVG, shifted down by ``offset``."""
    by_id = {box.ident: box for box in boxes}
    pairs = {(source, target) for source, target, _label, _colour in edges}
    parts = [f'  <g transform="translate(0,{offset})" font-family="system-ui,sans-serif">']
    for source, target, label, colour in edges:
        if source not in by_id or target not in by_id:
            raise ValueError(f"edge {source}->{target} names a cell that is not on this page")
        points, (label_x, label_y) = route(by_id[source], by_id[target])
        drawn = " ".join(f"{x:.0f},{y:.0f}" for x, y in points)
        marker = "arrow-remote" if colour != "#748b87" else "arrow"
        dash = ' stroke-dasharray="6 4"' if colour != "#748b87" else ""
        parts.append(
            f'    <polyline points="{drawn}" fill="none" stroke="{colour}" stroke-width="1.7"'
            f'{dash} marker-end="url(#{marker})" />')
        if label:
            # Reciprocal edges share a midpoint, so their labels land on top of each other.
            # Nudge the second one clear.
            if (target, source) in pairs and source > target:
                label_y += 15
            width = len(label) * 5.6 + 8
            parts.append(
                f'    <rect x="{label_x - width / 2:.0f}" y="{label_y - 10:.0f}" '
                f'width="{width:.0f}" height="14" fill="#f6f9f7" />')
            parts.append(
                f'    <text x="{label_x:.0f}" y="{label_y:.0f}" fill="#617974" font-size="11" '
                f'text-anchor="middle">{escape(label)}</text>')
    for box in boxes:
        if not box.is_text:
            parts.append(
                f'    <rect x="{box.x}" y="{box.y}" width="{box.width}" height="{box.height}" '
                f'rx="10" fill="{box.fill}" stroke="{box.stroke}" stroke-width="1.5" />')
        weight = "600" if box.size >= 20 else "400"
        # A heading sits on its own baseline; a box's lines are centred in it, with the first
        # line carrying the box's name in bold.
        top = (box.y + box.size + 6 if box.is_text
               else box.centre_y - (len(box.lines) - 1) * 10 + 5)
        lines = wrap(box.lines, box.width, box.size) if box.is_text else list(box.lines)
        for index, line in enumerate(lines):
            if not line:
                continue
            bold = "600" if (not box.is_text and index == 0) else weight
            size = box.size if index == 0 else max(11, box.size - 2)
            ink = box.ink if index == 0 else "#617974" if not box.is_text else box.ink
            x = box.x if box.is_text else box.x + 18
            parts.append(
                f'    <text x="{x}" y="{top + index * 20:.0f}" fill="{ink}" '
                f'font-size="{size}" font-weight="{bold}">{escape(line)}</text>')
    parts.append("  </g>")
    return "\n".join(parts)


def render(source: Path) -> str:
    """Return the whole SVG for every page in the diagram."""
    pages = read_pages(source)
    rendered, offset, names = [], 0, []
    for name, boxes, edges in pages:
        names.append(name)
        rendered.append(render_page(boxes, edges, offset))
        offset += max(box.y + box.height for box in boxes) + PAGE_GAP
    title = " and ".join(names)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="{offset}" '
        f'viewBox="0 0 1280 {offset}" role="img" aria-labelledby="flow-title">\n'
        f"  <title id=\"flow-title\">saRNA Designer: {escape(title.lower())}</title>\n"
        f"  <desc>Generated from website-flow.drawio by scripts/render_flow_diagram.py. "
        f"Edit the .drawio and re-run that script; do not edit this file.</desc>\n"
        f"  {ARROW}\n"
        f'  <rect width="1280" height="{offset}" fill="#f6f9f7" />\n'
        + "\n".join(rendered)
        + "\n</svg>\n"
    )


def main() -> None:
    """Write the SVG beside the diagram, or check it is current."""
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent / "docs" / "architecture"
    parser.add_argument("--source", type=Path, default=root / "website-flow.drawio")
    parser.add_argument("--output", type=Path, default=root / "website-flow.svg")
    parser.add_argument("--check", action="store_true",
                        help="Exit non-zero if the SVG is out of date, writing nothing.")
    arguments = parser.parse_args()
    produced = render(arguments.source)
    if arguments.check:
        current = (arguments.output.read_text(encoding="utf-8")
                   if arguments.output.exists() else "")
        if current != produced:
            raise SystemExit(
                f"{arguments.output.name} is out of date; run {Path(__file__).name}")
        print(f"{arguments.output.name} is up to date")
        return
    arguments.output.write_text(produced, encoding="utf-8")
    print(f"wrote {arguments.output} ({len(produced)} bytes)")


if __name__ == "__main__":
    main()
