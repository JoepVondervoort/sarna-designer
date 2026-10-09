"""Tests for presets and the resolved design configuration (ADR 0010).

The point of these is less the mechanics than the scientific claims: that `wang-2015` matches the
paper, that `legacy-scanner` matches the script where it is meant to, and that the two documented
exceptions are present and visible rather than quietly absent.
"""

from typing import Any

import pytest

from sarna.design.config import DesignConfig
from sarna.design.presets import (
    EXPLORATORY,
    LEGACY_SCANNER,
    PRESETS,
    STRICT,
    WANG_2015,
    preset,
)
from sarna.rules.catalog import CATALOGUE
from sarna.rules.catalog.optional import (
    AU_RICH_GUIDE_5_PRIME,
    CENTRAL_STABILITY,
    NO_G_AT_13,
    OPTIONAL_RULES,
    POS3_A,
    POS10_G,
    POS10_U,
)
from sarna.rules.provenance import ProvenanceKind
from sarna.rules.rule import RuleKind
from sarna.validation.benchmark import benchmark
from sarna.validation.dataset import load_chapter4_dataset


class TestPresetRegistry:
    def test_the_four_presets_ship(self) -> None:
        assert set(PRESETS) == {"wang-2015", "strict", "exploratory", "legacy-scanner"}

    def test_exactly_one_preset_is_the_default(self) -> None:
        assert sum(1 for p in PRESETS.values() if p.is_default) == 1

    def test_the_default_is_the_paper(self) -> None:
        assert WANG_2015.is_default is True
        assert LEGACY_SCANNER.is_default is False

    def test_a_preset_is_retrievable_by_id(self) -> None:
        assert preset("wang-2015") is WANG_2015

    def test_an_unknown_preset_raises_and_lists_what_exists(self) -> None:
        with pytest.raises(KeyError, match="wang-2015"):
            preset("wang-2016")

    def test_every_preset_explains_itself(self) -> None:
        for named in PRESETS.values():
            assert named.name.strip()
            assert named.description.strip()

    def test_every_selected_rule_exists_in_the_catalogue(self) -> None:
        # A preset naming a rule that does not exist would silently select nothing.
        for named in PRESETS.values():
            assert named.resolve()  # raises if a rule id is unknown


class TestWang2015:
    CONFIG = WANG_2015.resolve()

    def test_the_gc_ceiling_is_the_papers(self) -> None:
        gc = self.CONFIG.rule("gc-content")
        assert gc.number("maximum-gc-fraction") == pytest.approx(0.60)

    def test_one_cpg_is_tolerated(self) -> None:
        assert self.CONFIG.rule("cpg-count").integer("maximum-cpg-count") == 1

    def test_the_rna_table_is_used(self) -> None:
        asymmetry = self.CONFIG.rule("thermodynamic-asymmetry")
        assert asymmetry.text("nearest-neighbour-table") == "rna-xia-1998"

    def test_the_unsourced_position_seven_rule_is_not_selected(self) -> None:
        assert "pos7-u" not in self.CONFIG

    def test_the_tss_rules_are_selected(self) -> None:
        assert "tss-window" in self.CONFIG
        assert "tss-preferred-window" in self.CONFIG

    def test_the_island_and_repeat_rules_are_selected(self) -> None:
        assert "not-in-cpg-island" in self.CONFIG
        assert "no-simple-repeats" in self.CONFIG

    def test_there_are_seven_available_preference_points(self) -> None:
        # ADR 0006: two start positions, four flank positions, one preferred window.
        assert self.CONFIG.max_preference_points == 7


class TestLegacyScanner:
    CONFIG = LEGACY_SCANNER.resolve()

    def test_the_gc_ceiling_is_the_scripts(self) -> None:
        gc = self.CONFIG.rule("gc-content")
        assert gc.number("maximum-gc-fraction") == pytest.approx(0.65)

    def test_no_cpg_is_tolerated(self) -> None:
        assert self.CONFIG.rule("cpg-count").integer("maximum-cpg-count") == 0

    def test_the_defective_legacy_table_is_used(self) -> None:
        # Defect B2 preserved deliberately: a corrected legacy preset reproduces nothing.
        assert (
            self.CONFIG.rule("thermodynamic-asymmetry").text("nearest-neighbour-table")
            == "legacy-dna-kj"
        )

    def test_the_tiebreak_orders_by_the_same_table_that_admitted(self) -> None:
        assert (
            self.CONFIG.rule("asymmetry-margin").text("nearest-neighbour-table")
            == "legacy-dna-kj"
        )

    def test_the_unsourced_position_seven_rule_is_selected(self) -> None:
        # One of the script's two dominant ranking terms, so historical rankings need it.
        assert "pos7-u" in self.CONFIG

    def test_the_tss_rules_are_not_selected(self) -> None:
        # The script scans a pasted string and has no concept of a TSS at all.
        assert "tss-window" not in self.CONFIG
        assert "tss-preferred-window" not in self.CONFIG

    def test_the_cpg_island_rule_is_not_selected(self) -> None:
        # Defect B8: the script's `inisland` call is commented out, so the rule never ran.
        assert "not-in-cpg-island" not in self.CONFIG

    def test_the_simple_repeat_rule_is_not_selected(self) -> None:
        # Absent from the script entirely.
        assert "no-simple-repeats" not in self.CONFIG


class TestTheStatedExceptionsToReproduction:
    """ADR 0010 records two ways legacy-scanner departs from the original. Both are asserted."""

    def test_position_one_gc_is_a_preference_not_a_gate(self) -> None:
        # The script's `endbase` makes it mandatory (finding B16). The decision is to follow
        # wang-2015 and keep it a preference, which means legacy-scanner admits candidates the
        # original rejected. Asserted so the divergence cannot vanish unnoticed.
        config = LEGACY_SCANNER.resolve()
        start = config.rule("start-1-gc")
        assert start.rule.kind is RuleKind.PREFERENCE
        assert not start.rule.is_gate

    def test_no_mandatory_position_one_rule_exists_anywhere(self) -> None:
        config = LEGACY_SCANNER.resolve()
        gates = {resolved.rule.id for resolved in config.mandatory}
        assert "pos1-gc-required" not in gates

    def test_the_preset_declares_its_known_divergences(self) -> None:
        # Stated in the preset itself, so it reaches an export and the UI rather than living
        # only in an ADR nobody reads at the bench.
        assert LEGACY_SCANNER.divergences
        joined = " ".join(LEGACY_SCANNER.divergences).lower()
        assert "position 1" in joined
        assert "rank" in joined

    def test_the_default_preset_declares_no_divergences(self) -> None:
        assert WANG_2015.divergences == ()


class TestOverrides:
    def test_a_user_override_layers_over_the_preset(self) -> None:
        config = WANG_2015.resolve({"gc-content": {"maximum-gc-fraction": 0.62}})
        assert config.rule("gc-content").number("maximum-gc-fraction") == pytest.approx(0.62)

    def test_an_override_does_not_disturb_other_parameters_of_the_same_rule(self) -> None:
        config = WANG_2015.resolve({"gc-content": {"maximum-gc-fraction": 0.62}})
        assert config.rule("gc-content").number("minimum-gc-fraction") == pytest.approx(0.40)

    def test_an_override_for_an_unselected_rule_is_an_error(self) -> None:
        # Silently ignoring it would let a user believe they had changed something.
        with pytest.raises(KeyError, match="not selected"):
            LEGACY_SCANNER.resolve({"tss-window": {"window-end-bp": -50}})

    def test_an_invalid_override_is_rejected_with_the_parameter_named(self) -> None:
        with pytest.raises(ValueError, match="maximum-gc-fraction"):
            WANG_2015.resolve({"gc-content": {"maximum-gc-fraction": 2.0}})

    def test_resolving_twice_does_not_accumulate_overrides(self) -> None:
        WANG_2015.resolve({"gc-content": {"maximum-gc-fraction": 0.62}})
        assert WANG_2015.resolve().rule("gc-content").number(
            "maximum-gc-fraction"
        ) == pytest.approx(0.60)


class TestConfigShape:
    CONFIG = WANG_2015.resolve()

    def test_rules_are_grouped_by_kind(self) -> None:
        assert all(r.rule.kind is RuleKind.MANDATORY for r in self.CONFIG.mandatory)
        assert all(r.rule.kind is RuleKind.PREFERENCE for r in self.CONFIG.preferences)
        assert all(r.rule.kind is RuleKind.TIEBREAK for r in self.CONFIG.tiebreaks)

    def test_the_groups_partition_the_selection(self) -> None:
        grouped = (
            len(self.CONFIG.mandatory)
            + len(self.CONFIG.preferences)
            + len(self.CONFIG.tiebreaks)
        )
        assert grouped == len(self.CONFIG)

    def test_the_target_length_is_available_to_the_generator(self) -> None:
        # The scan needs it before any rule runs, so it is surfaced rather than dug out.
        assert self.CONFIG.target_length == 19

    def test_the_island_criteria_come_from_the_island_rule(self) -> None:
        criteria = self.CONFIG.cpg_island_criteria
        assert criteria is not None
        assert criteria.minimum_length == 500
        assert criteria.minimum_gc_fraction == pytest.approx(0.55)

    def test_a_preset_without_the_island_rule_has_no_criteria(self) -> None:
        # None, so the caller knows not to run detection at all, and the rule reports
        # not-applicable rather than the engine guessing a definition.
        assert LEGACY_SCANNER.resolve().cpg_island_criteria is None

    def test_an_unknown_rule_lookup_raises(self) -> None:
        with pytest.raises(KeyError):
            self.CONFIG.rule("no-such-rule")

    def test_is_immutable(self) -> None:
        with pytest.raises((AttributeError, TypeError)):
            self.CONFIG.preset_id = "other"  # type: ignore[misc]


def resolved_values(entry: dict[str, Any]) -> dict[str, Any]:
    """Return a rule export entry's resolved parameter values, keyed by name.

    The export carries each parameter as a self-describing object (value beside its unit, legal
    range and description) so a value lookup is derived rather than the shape on the wire.
    """
    return {p["name"]: p["value"] for p in entry["parameters"]}


class TestExport:
    def test_the_configuration_serialises_to_the_rules_that_ran(self) -> None:
        # ADR 0010: a preset name is not a reproducible record, because its meaning can change
        # between versions. The resolved values are.
        exported = WANG_2015.resolve().to_dict()
        assert exported["preset"]["id"] == "wang-2015"
        by_id = {entry["id"]: entry for entry in exported["rules"]}
        assert resolved_values(by_id["gc-content"])["maximum-gc-fraction"] == 0.60

    def test_each_exported_rule_carries_its_provenance_and_evidence_class(self) -> None:
        exported = WANG_2015.resolve().to_dict()
        for entry in exported["rules"]:
            assert entry["provenance"]
            assert entry["evidence_class"]
            assert entry["kind"]

    def test_the_export_records_the_divergences(self) -> None:
        exported = LEGACY_SCANNER.resolve().to_dict()
        assert exported["preset"]["divergences"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        for named in PRESETS.values():
            json.dumps(named.resolve().to_dict())

    def test_two_presets_can_be_diffed_by_their_exports(self) -> None:
        wang = {
            entry["id"]: resolved_values(entry)
            for entry in WANG_2015.resolve().to_dict()["rules"]
        }
        legacy = {
            entry["id"]: resolved_values(entry)
            for entry in LEGACY_SCANNER.resolve().to_dict()["rules"]
        }
        differing = {
            rule_id
            for rule_id in set(wang) & set(legacy)
            if wang[rule_id] != legacy[rule_id]
        }
        assert differing == {
            "gc-content",
            "cpg-count",
            "thermodynamic-asymmetry",
            "asymmetry-margin",
        }


class TestDesignConfigDirectly:
    def test_a_config_can_be_built_without_a_preset(self) -> None:
        # Needed for tests and notebook use: a bare selection with no preset identity.
        config = DesignConfig.of(["gc-content", "pos19-a"])
        assert len(config) == 2
        assert config.preset_id is None

    def test_a_bare_config_still_exports(self) -> None:
        exported = DesignConfig.of(["gc-content"]).to_dict()
        assert exported["preset"]["id"] is None
        assert exported["rules"][0]["id"] == "gc-content"

    def test_a_duplicate_rule_id_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="more than once"):
            DesignConfig.of(["gc-content", "gc-content"])

    def test_an_unknown_rule_id_names_itself(self) -> None:
        with pytest.raises(KeyError, match="gc-contents"):
            DesignConfig.of(["gc-contents"])


class TestIteration:
    def test_a_config_iterates_its_resolved_rules(self) -> None:
        # Used by the design pipeline to walk the gates and preferences in order.
        config = DesignConfig.of(["gc-content", "pos19-a"])
        assert [resolved.rule.id for resolved in config] == ["gc-content", "pos19-a"]


class TestStrictPreset:
    """Strict promotes Wang's softer statements to requirements, and says that it has."""

    CONFIG = STRICT.resolve()

    def test_it_admits_on_the_region_the_paper_reports_activity_in(self) -> None:
        window = self.CONFIG.rule("tss-window")
        assert window.number("window-start-bp") == -500
        assert window.number("window-end-bp") == -200

    def test_it_uses_the_looser_island_definition_to_exclude_more(self) -> None:
        # Gardiner-Garden & Frommer detects more and shorter islands than Takai & Jones, so
        # excluding any overlap with one rejects more candidates. The more permissive detector
        # is the stricter gate.
        island = self.CONFIG.rule("not-in-cpg-island")
        assert island.number("island-minimum-length") == 200
        assert island.number("island-minimum-gc-fraction") == pytest.approx(0.50)

    def test_it_leaves_wangs_sequence_criteria_alone(self) -> None:
        # The measured activators support them exactly as published, so there is nothing to
        # tighten. Only location criteria, which this dataset cannot test, are narrowed.
        assert self.CONFIG.rule("cpg-count").integer("maximum-cpg-count") == 1
        assert self.CONFIG.rule("gc-content").number("maximum-gc-fraction") == pytest.approx(0.60)
        assert self.CONFIG.rule("thermodynamic-asymmetry").number(
            "minimum-asymmetry-kcal"
        ) == pytest.approx(0.0)

    def test_it_does_not_reject_more_known_activators_than_the_default(self) -> None:
        # The check that caught the first version of this preset, which excluded 5 of the 14
        # measured activators by requiring zero CpG and a 0.5 kcal/mol asymmetry margin.
        dataset = load_chapter4_dataset()
        strict = benchmark(dataset, self.CONFIG, permutations=1)
        default = benchmark(dataset, WANG_2015.resolve(), permutations=1)
        assert len(strict.rejected_activators) <= len(default.rejected_activators)

    def test_it_adds_no_preferences_so_its_score_stays_comparable_with_the_default(self) -> None:
        assert self.CONFIG.max_preference_points == WANG_2015.resolve().max_preference_points

    def test_it_declares_that_it_is_stricter_than_the_paper(self) -> None:
        # Being stricter than the source is a divergence from it, and must be stated as one.
        joined = " ".join(STRICT.divergences).lower()
        assert "stricter than the paper" in joined
        assert "narrative text" in joined

    def test_it_declares_the_tightenings_it_measured_and_rejected(self) -> None:
        # A preset that declined an obvious setting for a measured reason should say so, or the
        # next person will "fix" it back.
        joined = " ".join(STRICT.divergences).lower()
        assert "measured and rejected" in joined
        assert "4 of the 14" in joined


class TestExploratoryPreset:
    CONFIG = EXPLORATORY.resolve()

    def test_it_selects_the_best_supported_optional_rule(self) -> None:
        assert "au-rich-guide-5-prime" in self.CONFIG

    def test_it_selects_the_chapter_four_hypotheses(self) -> None:
        assert "pos10-g" in self.CONFIG
        assert "central-stability" in self.CONFIG

    def test_it_omits_the_rule_that_contradicts_one_it_selects(self) -> None:
        # pos10-u (Reynolds, U at position 10) and pos10-g (chapter 4, G at position 10) cannot
        # both be met. Shipping both in one preset would guarantee every candidate loses a point.
        assert "pos10-u" not in self.CONFIG
        assert "pos10-g" in self.CONFIG

    def test_it_adds_no_gates_so_the_candidate_set_is_unchanged(self) -> None:
        # Every added rule is a preference, so this preset reorders rather than filters.
        assert len(self.CONFIG.mandatory) == len(WANG_2015.resolve().mandatory)

    def test_its_denominator_differs_from_the_default(self) -> None:
        assert self.CONFIG.max_preference_points == 13
        assert WANG_2015.resolve().max_preference_points == 7

    def test_it_warns_that_its_scores_are_not_comparable(self) -> None:
        joined = " ".join(EXPLORATORY.divergences).lower()
        assert "not comparable" in joined

    def test_it_warns_that_the_hypotheses_do_not_survive_correction(self) -> None:
        joined = " ".join(EXPLORATORY.divergences).lower()
        assert "multiple comparisons" in joined

    def test_it_warns_about_overlap_with_the_asymmetry_gate(self) -> None:
        joined = " ".join(EXPLORATORY.divergences).lower()
        assert "overlap" in joined


class TestTheDefaultIsUnaffectedByTheNewRules:
    """Adding optional rules to the catalogue must not change what wang-2015 does (ADR 0011)."""

    def test_the_default_still_has_seven_preferences(self) -> None:
        assert WANG_2015.resolve().max_preference_points == 7

    def test_the_default_selects_no_optional_rule(self) -> None:
        config = WANG_2015.resolve()
        for rule_id in (
            "au-rich-guide-5-prime",
            "no-gc-stretch",
            "pos3-a",
            "pos10-u",
            "pos10-g",
            "no-g-at-13",
            "central-stability",
        ):
            assert rule_id not in config

    def test_every_optional_rule_is_off_by_default_in_the_catalogue(self) -> None:
        for rule in OPTIONAL_RULES:
            assert rule.enabled_by_default is False

    def test_every_optional_rule_is_a_preference_and_never_a_gate(self) -> None:
        # So selecting any of them can only reorder, never exclude.
        for rule in OPTIONAL_RULES:
            assert rule.kind is RuleKind.PREFERENCE
            assert not rule.is_gate


class TestOptionalRuleProvenance:
    """These rules exist to be judged, so each must publish enough to judge it by."""

    def test_every_optional_rule_cites_a_source(self) -> None:
        for rule in OPTIONAL_RULES:
            assert rule.provenance.kind is ProvenanceKind.CITATION
            assert rule.provenance.source
            assert rule.provenance.locator

    def test_the_knockdown_derived_rules_say_so(self) -> None:
        # The single most important caveat: these were derived for knockdown, not activation.
        for rule in (POS3_A, POS10_U, NO_G_AT_13):
            assert "knockdown" in rule.consequence.lower()

    def test_the_chapter_four_rules_state_that_they_fail_correction(self) -> None:
        for rule in (POS10_G, CENTRAL_STABILITY):
            assert "does not survive correction" in rule.consequence

    def test_the_contradiction_between_position_10_rules_is_stated_on_both(self) -> None:
        assert "pos10-g" in POS10_U.consequence
        assert "pos10-u" in POS10_G.consequence

    def test_the_best_supported_rule_names_its_three_agreeing_sources(self) -> None:
        rationale = AU_RICH_GUIDE_5_PRIME.rationale
        assert "Ui-Tei" in rationale
        assert "Reynolds" in rationale
        assert "Chapter 4" in rationale or "chapter 4" in rationale

    def test_the_overlap_with_the_asymmetry_rule_is_disclosed(self) -> None:
        # A preference that re-scores what a gate already measured must say so, or a reader
        # will read two independent lines of evidence where there is one.
        assert "asymmetry" in AU_RICH_GUIDE_5_PRIME.consequence
        assert "double-count" in CENTRAL_STABILITY.consequence


class TestTheReynoldsGcOption:
    def test_the_gc_rule_offers_the_reynolds_band(self) -> None:
        option = CATALOGUE["gc-content"].option("reynolds-2004")
        assert option.values == {"minimum-gc-fraction": 0.30, "maximum-gc-fraction": 0.52}

    def test_it_explains_that_it_is_a_different_band_not_a_tighter_one(self) -> None:
        # It is lower at both ends than Wang's, so "stricter" would be the wrong mental model.
        note = CATALOGUE["gc-content"].option("reynolds-2004").note.lower()
        assert "not a tightening" in note

    def test_the_wang_band_remains_the_default(self) -> None:
        default = CATALOGUE["gc-content"].default_option
        assert default is not None
        assert default.id == "wang-2015"
