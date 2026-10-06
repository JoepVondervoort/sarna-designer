"""Provenance is mandatory and self-validating.

ADR 0003: a rule's provenance is a citation *with a locator*, or `legacy-script`, or
`unsourced`. Silence is not permitted, and `unsourced` renders as a warning rather than
passing as established biology.
"""

import pytest
from pydantic import ValidationError

from sarna.rules.provenance import EvidenceClass, Provenance, ProvenanceKind


class TestCitation:
    def test_builds_with_source_and_locator(self):
        p = Provenance.citation("Wang et al. 2015", "Table 1")
        assert p.kind is ProvenanceKind.CITATION
        assert p.source == "Wang et al. 2015"
        assert p.locator == "Table 1"

    def test_renders_source_and_locator_together(self):
        assert str(Provenance.citation("Wang et al. 2015", "Table 1")) == (
            "Wang et al. 2015, Table 1"
        )

    def test_a_citation_without_a_locator_is_refused(self):
        # "Wang et al. 2015" alone is not checkable; the locator is the point.
        with pytest.raises(ValidationError, match="locator"):
            Provenance.citation("Wang et al. 2015", "")

    def test_a_citation_without_a_source_is_refused(self):
        with pytest.raises(ValidationError, match="source"):
            Provenance.citation("", "Table 1")

    def test_blank_locator_is_refused(self):
        with pytest.raises(ValidationError, match="locator"):
            Provenance.citation("Wang et al. 2015", "   ")

    def test_is_not_a_warning(self):
        assert Provenance.citation("Wang et al. 2015", "Table 1").is_warning is False


class TestLegacyScript:
    def test_builds_with_a_locator(self):
        p = Provenance.legacy_script("sarna_promoter_scanner.py:checkpos19")
        assert p.kind is ProvenanceKind.LEGACY_SCRIPT
        assert p.locator == "sarna_promoter_scanner.py:checkpos19"

    def test_requires_a_locator(self):
        with pytest.raises(ValidationError, match="locator"):
            Provenance.legacy_script("")

    def test_renders_as_inherited_from_the_legacy_script(self):
        rendered = str(Provenance.legacy_script("sarna_promoter_scanner.py:checkpos7"))
        assert "legacy script" in rendered
        assert "checkpos7" in rendered

    def test_is_not_a_warning_but_is_not_a_citation_either(self):
        p = Provenance.legacy_script("sarna_promoter_scanner.py:checkpos19")
        assert p.is_warning is False
        assert p.kind is not ProvenanceKind.CITATION


class TestUnsourced:
    def test_builds_with_a_note(self):
        p = Provenance.unsourced("Position 7 appears nowhere in Wang 2015.")
        assert p.kind is ProvenanceKind.UNSOURCED
        assert p.source is None
        assert p.locator is None

    def test_requires_a_note_explaining_what_was_looked_for(self):
        # "unsourced" must be a finding, not a shrug.
        with pytest.raises(ValidationError, match="note"):
            Provenance.unsourced("")

    def test_is_a_warning(self):
        assert Provenance.unsourced("no source found").is_warning is True

    def test_renders_as_a_warning_carrying_the_note(self):
        rendered = str(Provenance.unsourced("Position 7 appears nowhere in Wang 2015."))
        assert "unsourced" in rendered.lower()
        assert "Position 7" in rendered


class TestProvenanceIsData:
    def test_is_immutable(self):
        p = Provenance.citation("Wang et al. 2015", "Table 1")
        with pytest.raises(ValidationError):
            p.source = "something else"

    def test_round_trips_through_serialisation(self):
        # ADR 0003: exports carry the rule records, so they must serialise losslessly.
        p = Provenance.citation("Wang et al. 2015", "Table 1")
        assert Provenance.model_validate(p.model_dump()) == p

    def test_unsourced_round_trips_too(self):
        p = Provenance.unsourced("no source found")
        assert Provenance.model_validate(p.model_dump()) == p


class TestEvidenceClass:
    def test_the_three_classes_exist(self):
        assert {e.value for e in EvidenceClass} == {"established", "empirical", "exploratory"}

    def test_renders_as_its_label(self):
        assert str(EvidenceClass.EXPLORATORY) == "exploratory"


class TestProvenanceKind:
    def test_renders_as_its_label(self):
        assert str(ProvenanceKind.LEGACY_SCRIPT) == "legacy-script"
        assert str(ProvenanceKind.UNSOURCED) == "unsourced"
