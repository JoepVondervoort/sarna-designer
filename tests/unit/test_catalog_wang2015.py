"""The Wang 2015 rule catalogue.

These tests pin the scientific content, not just the plumbing: the criteria that must exist,
the defaults that must match Table 1, and the honesty invariants that stop an unsourced rule
being presented as settled biology.
"""

import pytest

from sarna.rules.catalog import CATALOGUE
from sarna.rules.provenance import EvidenceClass, ProvenanceKind
from sarna.rules.rule import RuleKind

TABLE_1_MANDATORY = {
    "target-length-19",
    "gc-content",
    "no-homopolymer-run",
    "cpg-count",
    "not-in-cpg-island",
    "no-simple-repeats",
    "tss-window",
    "pos18-at",
    "pos19-a",
    "thermodynamic-asymmetry",
}

TABLE_1_PREFERENCES = {
    "start-1-gc",
    "start-2-gc",
    "flank-20-au",
    "flank-21-au",
    "flank-22-au",
    "flank-23-au",
    "tss-preferred-window",
}


class TestCoverageOfTable1:
    def test_every_mandatory_criterion_is_present(self):
        assert {r.id for r in CATALOGUE.mandatory} == TABLE_1_MANDATORY

    def test_the_seven_scored_preferences_are_present(self):
        enabled = {r.id for r in CATALOGUE.preferences if r.enabled_by_default}
        assert enabled == TABLE_1_PREFERENCES

    def test_there_is_exactly_one_tiebreak(self):
        assert [r.id for r in CATALOGUE.tiebreaks] == ["asymmetry-margin"]

    def test_the_default_preference_ceiling_is_seven(self):
        # The denominator in "meets 5 of 7 preferences" (ADR 0006).
        enabled = [r for r in CATALOGUE.preferences if r.enabled_by_default]
        assert sum(r.points for r in enabled) == 7


class TestEveryRuleIsExplained:
    @pytest.mark.parametrize("rule_id", sorted(TABLE_1_MANDATORY | TABLE_1_PREFERENCES))
    def test_carries_rationale_and_consequence(self, rule_id: str):
        rule = CATALOGUE[rule_id]
        assert rule.rationale.strip()
        assert rule.consequence.strip()

    def test_no_rule_is_silently_unsourced(self):
        # unsourced is permitted, but only when flagged as such and never as established.
        for rule in CATALOGUE:
            if rule.provenance.kind is ProvenanceKind.UNSOURCED:
                assert rule.needs_warning
                assert rule.evidence_class is not EvidenceClass.ESTABLISHED

    def test_cited_rules_name_a_locator(self):
        for rule in CATALOGUE:
            if rule.provenance.kind is ProvenanceKind.CITATION:
                assert rule.provenance.locator


class TestDefaultsMatchTable1:
    def test_gc_band_is_forty_to_sixty_percent(self):
        gc = CATALOGUE["gc-content"]
        assert gc.parameter("minimum-gc-fraction").default == 0.40
        # 0.65 is the legacy deviation; the default must be the paper's 0.60.
        assert gc.parameter("maximum-gc-fraction").default == 0.60

    def test_at_most_one_cpg_is_allowed(self):
        # The original script tolerated zero, which is stricter than the paper.
        assert CATALOGUE["cpg-count"].parameter("maximum-cpg-count").default == 1

    def test_cpg_islands_are_defined_by_the_takai_jones_thresholds(self):
        # Wang et al. say to avoid CpG islands but name no detection algorithm, so the
        # definition is this project's choice and is pinned here. Takai & Jones 2002 rather than
        # Gardiner-Garden & Frommer 1987, because this rule excludes candidates and the 1987
        # criteria call more than half their islands on Alu repeats.
        island = CATALOGUE["not-in-cpg-island"]
        assert island.parameter("island-minimum-length").default == 500
        assert island.parameter("island-minimum-gc-fraction").default == 0.55
        assert island.parameter("island-minimum-observed-expected-cpg").default == 0.65

    def test_both_cpg_island_definitions_are_offered_with_their_own_sources(self):
        island = CATALOGUE["not-in-cpg-island"]
        assert {option.id for option in island.options} == {
            "takai-jones-2002",
            "gardiner-garden-frommer-1987",
        }
        default = island.default_option
        assert default is not None
        assert default.id == "takai-jones-2002"
        # The looser definition must not be presented as unsourced: it is the original paper.
        looser = island.option("gardiner-garden-frommer-1987")
        assert looser.provenance.is_warning is False
        assert "1987" in str(looser.provenance)

    def test_homopolymer_runs_of_five_or_more_are_rejected(self):
        assert CATALOGUE["no-homopolymer-run"].parameter("maximum-run-length").default == 4

    def test_target_length_is_nineteen(self):
        assert CATALOGUE["target-length-19"].parameter("target-length").default == 19

    def test_the_mandatory_tss_window_is_minus_one_thousand_to_minus_one_hundred(self):
        window = CATALOGUE["tss-window"]
        assert window.parameter("window-start-bp").default == -1000
        assert window.parameter("window-end-bp").default == -100

    def test_the_preferred_tss_window_sits_inside_the_mandatory_one(self):
        preferred = CATALOGUE["tss-preferred-window"]
        mandatory = CATALOGUE["tss-window"]
        start = float(preferred.parameter("preferred-window-start-bp").default)
        end = float(preferred.parameter("preferred-window-end-bp").default)
        assert start == -500
        assert end == -200
        assert float(mandatory.parameter("window-start-bp").default) <= start
        assert end <= float(mandatory.parameter("window-end-bp").default)


class TestThePositionSevenRule:
    def test_is_present_but_off_by_default(self):
        # Deleting it would make historical legacy rankings unexplainable (ADR 0002).
        rule = CATALOGUE["pos7-u"]
        assert rule.enabled_by_default is False

    def test_is_unsourced_and_exploratory_and_warns(self):
        rule = CATALOGUE["pos7-u"]
        assert rule.provenance.kind is ProvenanceKind.UNSOURCED
        assert rule.evidence_class is EvidenceClass.EXPLORATORY
        assert rule.needs_warning is True

    def test_its_note_records_what_was_searched_for(self):
        assert "Wang" in (CATALOGUE["pos7-u"].provenance.note or "")

    def test_if_enabled_it_is_worth_one_point_like_any_other_preference(self):
        # Never a dominant term as it was in the original script.
        assert CATALOGUE["pos7-u"].kind is RuleKind.PREFERENCE
        assert CATALOGUE["pos7-u"].points == 1


class TestGrouping:
    def test_the_four_flank_preferences_share_a_group(self):
        groups = {CATALOGUE[f"flank-{p}-au"].group for p in (20, 21, 22, 23)}
        assert len(groups) == 1
        assert groups.pop()

    def test_the_two_start_preferences_share_a_group(self):
        assert CATALOGUE["start-1-gc"].group == CATALOGUE["start-2-gc"].group
        assert CATALOGUE["start-1-gc"].group is not None


class TestParametersAreUsable:
    def test_every_numeric_parameter_declares_a_unit(self):
        for rule in CATALOGUE:
            for parameter in rule.parameters:
                if not parameter.is_boolean:
                    assert parameter.unit, f"{rule.id}.{parameter.name} has no unit"

    def test_every_parameter_default_is_within_its_own_bounds(self):
        for rule in CATALOGUE:
            for parameter in rule.parameters:
                assert parameter.validated(parameter.default) == parameter.default

    def test_the_whole_catalogue_serialises(self):
        # ADR 0003: exports carry the rule records.
        assert all(r.model_dump(mode="json")["id"] == r.id for r in CATALOGUE)
