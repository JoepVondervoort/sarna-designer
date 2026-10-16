"""Tests for what a rule record carries into an export.

ADR 0003 commits to exports that carry the rule *records*, not just the resolved values. They did
not. ``rationale`` and ``consequence`` are required non-blank on all 26 shipped rules (the largest
body of "why this choice" text in the project) and reached no export and no screen. Nor did a
parameter's unit, range or description, so a reader saw ``0.6`` with no way to know what it measured
or what values were legal. Nor did the options, each of which carries its own provenance and a
required statement of what choosing it costs.

Provenance was exported, but flattened to one display string, which means a consumer cannot tell a
citation from an unsourced rule without parsing English. ``pos7-u`` is the one shipped rule with no
source, and nothing downstream could mark it.

These tests pin the structure rather than the prose: a renderer that joins on ``rule_id`` must find
every field it needs to explain a verdict without this codebase to hand.
"""

import json
from typing import Any

import pytest

from sarna.design.presets import LEGACY_SCANNER, WANG_2015
from sarna.rules.catalog import CATALOGUE
from sarna.rules.provenance import Provenance, ProvenanceKind


def parameters_by_name(entry: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return a rule export entry's parameters keyed by name."""
    return {parameter["name"]: parameter for parameter in entry["parameters"]}


class TestProvenanceExport:
    def test_a_citation_keeps_its_source_and_locator_apart(self) -> None:
        exported = Provenance.citation("Wang et al. 2015", "Table 1").to_dict()
        assert exported["kind"] == "citation"
        assert exported["source"] == "Wang et al. 2015"
        assert exported["locator"] == "Table 1"
        assert exported["is_warning"] is False

    def test_the_display_string_is_kept_alongside_the_structure(self) -> None:
        # So a consumer that only wants one line does not have to reassemble it, and the line it
        # gets is the same one the rest of the project renders.
        provenance = Provenance.citation("Wang et al. 2015", "Table 1")
        assert provenance.to_dict()["display"] == str(provenance)

    def test_an_unsourced_rule_is_flagged_as_a_warning(self) -> None:
        exported = Provenance.unsourced("searched Table 1 and the text; not found").to_dict()
        assert exported["kind"] == "unsourced"
        assert exported["is_warning"] is True
        assert "not found" in exported["note"]

    def test_a_legacy_script_names_the_line_it_came_from(self) -> None:
        exported = Provenance.legacy_script("sarna_promoter_scanner.py:checkpos19").to_dict()
        assert exported["kind"] == "legacy-script"
        assert exported["locator"] == "sarna_promoter_scanner.py:checkpos19"
        assert exported["is_warning"] is False


class TestRuleExport:
    def test_the_reasoning_text_is_exported(self) -> None:
        exported = CATALOGUE["gc-content"].to_dict()
        assert exported["rationale"].strip()
        assert exported["consequence"].strip()

    def test_every_shipped_rule_exports_both(self) -> None:
        # Required on construction, so an empty one here would mean the export dropped it.
        for rule in CATALOGUE:
            exported = rule.to_dict()
            assert exported["rationale"].strip(), rule.id
            assert exported["consequence"].strip(), rule.id

    def test_a_parameter_carries_its_unit_range_and_description(self) -> None:
        parameter = parameters_by_name(CATALOGUE["gc-content"].to_dict())["maximum-gc-fraction"]
        assert parameter["unit"] == "fraction"
        assert parameter["minimum"] is not None
        assert parameter["maximum"] is not None
        assert parameter["description"].strip()

    def test_options_are_exported_with_their_own_provenance(self) -> None:
        options = {option["id"]: option for option in CATALOGUE["gc-content"].to_dict()["options"]}
        assert options["wang-2015"]["provenance"]["kind"] == "citation"
        assert options["legacy-scanner"]["provenance"]["kind"] == "legacy-script"
        assert options["reynolds-2004"]["provenance"]["source"] == "Reynolds et al. 2004"

    def test_an_option_carries_the_trade_off_it_costs(self) -> None:
        options = {option["id"]: option for option in CATALOGUE["gc-content"].to_dict()["options"]}
        assert options["legacy-scanner"]["note"].strip()

    def test_the_unsourced_rule_is_identifiable_in_its_export(self) -> None:
        exported = CATALOGUE["pos7-u"].to_dict()
        assert exported["provenance"]["kind"] == ProvenanceKind.UNSOURCED.value
        assert exported["provenance"]["is_warning"] is True
        assert exported["enabled_by_default"] is False

    def test_every_shipped_rule_exports_as_json(self) -> None:
        json.dumps([rule.to_dict() for rule in CATALOGUE])


class TestConfigurationExport:
    """The per-run record. A reader must be able to judge it without the codebase."""

    def exported(self) -> dict[str, Any]:
        return WANG_2015.resolve().to_dict()

    def rule(self, rule_id: str, exported: dict[str, Any] | None = None) -> dict[str, Any]:
        entries: dict[str, dict[str, Any]] = {
            entry["id"]: entry for entry in (exported or self.exported())["rules"]
        }
        return entries[rule_id]

    def test_a_parameter_shows_the_value_used_beside_its_declaration(self) -> None:
        # One merged list rather than two structures to join: the renderer needs both at once,
        # because "0.6" is meaningless without "of a possible 0.0 to 1.0, as a fraction".
        parameter = parameters_by_name(self.rule("gc-content"))["maximum-gc-fraction"]
        assert parameter["value"] == 0.60
        assert parameter["default"] == 0.60
        assert parameter["unit"] == "fraction"
        assert parameter["description"].strip()

    def test_a_preset_that_diverges_from_a_default_says_so(self) -> None:
        legacy = self.rule("gc-content", LEGACY_SCANNER.resolve().to_dict())
        parameter = parameters_by_name(legacy)["maximum-gc-fraction"]
        assert parameter["value"] == 0.65
        assert parameter["default"] == 0.60
        assert parameter["is_default"] is False

    def test_a_parameter_left_alone_is_marked_as_the_default(self) -> None:
        parameter = parameters_by_name(self.rule("gc-content"))["maximum-gc-fraction"]
        assert parameter["is_default"] is True

    def test_provenance_arrives_structured_not_as_prose(self) -> None:
        provenance = self.rule("gc-content")["provenance"]
        assert provenance["source"] == "Wang et al. 2015"
        assert provenance["locator"] == "Table 1"

    def test_the_reasoning_travels_with_the_run(self) -> None:
        entry = self.rule("gc-content")
        assert entry["rationale"].strip()
        assert entry["consequence"].strip()

    def test_the_alternatives_not_chosen_travel_too(self) -> None:
        # A reader judging a threshold needs to know what else was on offer and why it lost.
        option_ids = {option["id"] for option in self.rule("gc-content")["options"]}
        assert {"wang-2015", "legacy-scanner", "reynolds-2004"} <= option_ids

    def test_a_run_selecting_the_unsourced_rule_can_be_warned_about(self) -> None:
        legacy = LEGACY_SCANNER.resolve().to_dict()
        warned = [
            entry["id"] for entry in legacy["rules"] if entry["provenance"]["is_warning"]
        ]
        assert warned == ["pos7-u"]

    def test_the_default_preset_carries_no_warning(self) -> None:
        assert not [
            entry for entry in self.exported()["rules"] if entry["provenance"]["is_warning"]
        ]

    def test_the_whole_configuration_exports_as_json(self) -> None:
        json.dumps(self.exported())

    @pytest.mark.parametrize("field", ["id", "name", "kind", "evidence_class"])
    def test_the_fields_that_were_already_exported_are_still_there(self, field: str) -> None:
        assert self.rule("gc-content")[field]


def test_rule_evidence_distinguishes_source_scope_from_activation_validation() -> None:
    assert 'protocol guidance' in CATALOGUE['pos19-a'].to_dict()['evidence_note']
    sirna = CATALOGUE['au-rich-guide-5-prime'].to_dict()
    assert 'extrapolation' in sirna['evidence_note']
    assert sirna['provenance']['url'] == 'https://doi.org/10.1093/nar/gkh247'
    assert 'unsupported' in CATALOGUE['pos7-u'].to_dict()['evidence_note']
    assert CATALOGUE['pos7-u'].to_dict()['provenance']['url'] is None
    assert Provenance.citation('Unverified source', 'Table 1').to_dict()['url'] is None
    for rule in WANG_2015.resolve().to_dict()['rules']:
        assert 'not evidence strength' in rule['application_note']
