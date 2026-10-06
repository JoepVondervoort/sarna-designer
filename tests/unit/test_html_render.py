"""Tests for the self-contained HTML renderer.

A rendered page is a scientific record, so the checks are about the record surviving: that the
run's data is really embedded, that it cannot break out of its script block, and that a template
which has lost its placeholder fails loudly rather than producing an empty browser.
"""

import json
import re
from typing import Any

import pytest

from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import Gene
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.sequence.nucleic import NucleicSequence
from sarna.view.html import DATA_PLACEHOLDER, FALLBACK_TITLE, render
from sarna.view.promoter_view import PromoterView

BASES = ("GCTACGGTATCATGTCATAAGTTCAATTGACATTGGAATTACCATAGTTAACGTTTAGCATA" * 10)[:600]


def view(with_gene: bool = False) -> PromoterView:
    promoter = Promoter(
        sequence=NucleicSequence.parse(BASES),
        tss_start=TssPosition(-500),
        interval=GenomicInterval(Assembly.GRCH38, "2", 100_000, 100_600, Strand.MINUS),
    )
    gene = (
        Gene(
            ensembl_id="ENSG00000144285",
            symbol="SCN1A",
            interval=GenomicInterval(Assembly.GRCH38, "2", 100_000, 200_000),
            description="sodium voltage-gated channel alpha subunit 1",
        )
        if with_gene
        else None
    )
    return PromoterView(
        result=design_sarnas(promoter, WANG_2015.resolve()), gene=gene
    )


def embedded(html: str) -> dict[str, Any]:
    """Return the run's data as the page's own script block holds it."""
    block = re.search(
        r'<script id="view-data" type="application/json">(.*?)</script>', html, re.S
    )
    assert block is not None
    decoded: dict[str, Any] = json.loads(block.group(1).replace("<\\/", "</"))
    return decoded


class TestTheRecordIsSelfContained:
    def test_the_runs_data_is_embedded_in_the_page(self) -> None:
        # Self-contained on purpose: a record that needs a server or a sibling file stops being
        # reproducible as soon as either moves.
        data = embedded(render(view()))
        assert data["counts"]["generated"] > 0
        assert data["promoter"]["sequence"] == BASES

    def test_the_placeholder_is_gone(self) -> None:
        assert DATA_PLACEHOLDER not in render(view())

    def test_the_resolved_configuration_travels_with_the_page(self) -> None:
        data = embedded(render(view()))
        assert data["configuration"]["rules"]
        assert data["configuration"]["preset"]["id"] == "wang-2015"

    def test_ranks_survive_into_the_page(self) -> None:
        ranks = [c["rank"] for c in embedded(render(view()))["candidates"] if c["rank"]]
        assert ranks
        assert min(ranks) == 1

    def test_rejected_candidates_are_abbreviated_by_default_but_still_present(self) -> None:
        data = embedded(render(view()))
        assert len(data["candidates"]) == data["counts"]["generated"]
        assert any(c.get("abbreviated") for c in data["candidates"])

    def test_full_rules_can_be_requested(self) -> None:
        data = embedded(render(view(), full_rules_for_rejected=True))
        assert not any(c.get("abbreviated") for c in data["candidates"])


class TestTheDataCannotBreakOut:
    def test_no_closing_script_tag_appears_in_the_data_block(self) -> None:
        # A literal </script inside the block would close it early and the page would render
        # the rest of the run as text.
        html = render(view())
        block = re.search(
            r'<script id="view-data" type="application/json">(.*?)</script>', html, re.S
        )
        assert block is not None
        assert "</script" not in block.group(1)

    def test_javascript_line_terminators_are_escaped(self) -> None:
        # U+2028 and U+2029 terminate a line in JavaScript but not in JSON.
        html = render(view())
        assert "\u2028" not in html
        assert "\u2029" not in html

    def test_the_page_carries_no_document_skeleton(self) -> None:
        # The publish step wraps the file; a second doctype or <body> would nest documents.
        html = render(view())
        assert not re.search(r"<!doctype|<html|<head>|<body>", html, re.I)


class TestTitle:
    def test_a_page_for_a_gene_is_named_after_it(self) -> None:
        assert "<title>SCN1A Promoter Browser</title>" in render(view(with_gene=True))

    def test_a_page_without_a_gene_falls_back(self) -> None:
        assert f"<title>{FALLBACK_TITLE}</title>" in render(view())


class TestAMissingPlaceholderFailsLoudly:
    def test_a_template_without_the_data_placeholder_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Otherwise the page renders an empty browser and reports nothing wrong.
        monkeypatch.setattr(
            "sarna.view.html.TEMPLATE_FILENAME", "jaspar_core_selection.json"
        )
        monkeypatch.setattr("sarna.view.html.TEMPLATE_PACKAGE", "sarna.regulatory.data")
        with pytest.raises(ValueError, match="render empty"):
            render(view())
