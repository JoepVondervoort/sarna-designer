"""Workspace tabs must own their controls and remain a valid accessible document."""

from html.parser import HTMLParser
from pathlib import Path


class WorkspaceLayout(HTMLParser):
    def __init__(self, source: str) -> None:
        """Collect IDs and panel ancestry while checking the static document structure."""
        super().__init__()
        self.stack: list[tuple[str, dict[str, str | None]]] = []
        self.nodes: dict[str, dict[str, str | None]] = {}
        self.panels: dict[str, str | None] = {}
        self.parents: dict[str, list[str | None]] = {}
        self.feed(source)
        assert not self.stack

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        identifier = values.get("id")
        if identifier:
            assert identifier not in self.nodes, f"Duplicate id: {identifier}"
            self.nodes[identifier] = values
            self.panels[identifier] = next((
                a.get("id") for _, a in reversed(self.stack) if a.get("role") == "tabpanel"
            ), None)
            self.parents[identifier] = [a.get("id") for _, a in self.stack]
        if tag not in {"input", "meta", "link", "br", "hr", "img"}:
            self.stack.append((tag, values))

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        pass

    def handle_endtag(self, tag: str) -> None:
        opened, _ = self.stack.pop()
        assert opened == tag, f"Mismatched HTML tags: {opened} / {tag}"


def test_four_accessible_workspaces_own_their_controls() -> None:
    layout = WorkspaceLayout(Path("src/sarna/app/static/index.html").read_text())
    tabs = {id: node for id, node in layout.nodes.items() if node.get("role") == "tab"}
    assert set(tabs) == {f"tab-{key}" for key in ("design", "candidates", "screening", "export")}
    for identifier, tab in tabs.items():
        panel = layout.nodes[str(tab["aria-controls"])]
        assert panel["role"] == "tabpanel"
        assert panel["aria-labelledby"] == identifier
        assert ("hidden" in panel) == (identifier != "tab-design")
    ownership = {
        "design": ["setup", "gene-overview", "empty-state", "evidence-guide"],
        "candidates": ["screened-shortlist", "candidate-layout", "design-comparison",
                       "shortlist-card", "promoter-landscape"],
        "screening": ["local-reference-panel", "offtarget-card", "offtarget-results",
                      "screening-comparisons"],
        "export": ["export-workspace", "export-selection", "reproducibility-panel"],
    }
    for key, ids in ownership.items():
        for identifier in ids:
            assert layout.panels[identifier] == f"workspace-{key}"
    for identifier in ["gc-min", "gc-max", "window-start", "window-end", "open-settings"]:
        assert "design-form" in layout.parents[identifier]
        assert "advanced-design-settings" in layout.parents[identifier]
    assert layout.panels["duplex-sheet"] is None  # Printing must not depend on the active tab.


def test_promoter_is_directly_visible_and_exports_have_format_descriptions() -> None:
    source = Path("src/sarna/app/static/index.html").read_text()
    layout = WorkspaceLayout(source)
    assert layout.parents["promoter-landscape"][-1] == "candidate-run"
    assert source.index('id="promoter-landscape"') < source.index('id="result-metrics"')
    for identifier in ["export-xlsx", "export-analysis", "export-order", "export-csv",
                       "export-json", "export-html", "export-plot"]:
        assert layout.panels[identifier] == "workspace-export"
