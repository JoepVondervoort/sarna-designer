"""Custom scientific roles must agree in gates, ranks, provenance and archived pages."""

import json
import re

import pytest

from sarna.design.config import DesignConfig
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import STRICT, WANG_2015
from sarna.genome.promoter import Promoter
from sarna.rules.catalog import CATALOGUE
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleOutcome
from sarna.rules.rule import RuleKind
from sarna.sequence.nucleic import NucleicSequence
from sarna.view.html import render
from sarna.view.promoter_view import PromoterView


def three_windows() -> Promoter:
    # Three 19 nt windows end in A, T, A; their first bases are G, C, T.
    return Promoter(sequence=NucleicSequence.parse("GCTATGCTATGCTATGCT" + "ATA"))


@pytest.mark.parametrize(
    ("role", "ranked_offsets", "points", "denominator"),
    [
        (RuleKind.MANDATORY, [0, 2], [0, 0], 0),
        (RuleKind.PREFERENCE, [0, 2, 1], [1, 1, 0], 1),
        (None, [0, 1, 2], [0, 0, 0], 0),
    ],
)
def test_roles_change_gating_and_ranking_with_hand_computed_expectations(
    role: RuleKind | None, ranked_offsets: list[int], points: list[int], denominator: int,
) -> None:
    config = DesignConfig.of(
        ["target-length-19", "pos19-a"], role_overrides={"pos19-a": role},
    )
    run = design_sarnas(three_windows(), config)
    assert len(run.generated) == 3
    assert [c.offset for c in run.ranked] == ranked_offsets
    assert [c.preference_points for c in run.ranked] == points
    assert all(c.max_preference_points == denominator for c in run.ranked)
    if role is None:
        assert all(all(r.rule_id != "pos19-a" for r in c.results) for c in run.generated)
        assert config.to_dict()["disabled_rules"][0]["kind"] == "off"
    else:
        assert all(c.result("pos19-a").kind is role for c in run.generated)


def test_promoting_a_preference_to_a_gate_changes_exclusions() -> None:
    config = DesignConfig.of(
        ["target-length-19", "start-1-gc"],
        role_overrides={"start-1-gc": RuleKind.MANDATORY},
    )
    run = design_sarnas(three_windows(), config)
    assert [c.offset for c in run.ranked] == [0, 1]
    assert run.rejected[0].excluded_by[0].rule_id == "start-1-gc"
    assert config.max_preference_points == 0


def test_missing_inputs_remain_unevaluated_when_role_changes() -> None:
    config = DesignConfig.of(
        ["target-length-19", "tss-window"],
        role_overrides={"tss-window": RuleKind.PREFERENCE},
    )
    run = design_sarnas(three_windows(), config)
    assert len(run.ranked) == 3
    assert all(c.result("tss-window").outcome is RuleOutcome.NOT_APPLICABLE for c in run.ranked)
    assert all(c.evaluated_preference_points == 0 for c in run.ranked)


def test_custom_wang_retains_default_and_baseline_and_disabled_sources() -> None:
    before = WANG_2015.resolve().to_dict()
    custom = WANG_2015.resolve(role_overrides={
        "tss-window": RuleKind.PREFERENCE, "pos19-a": None, "pos3-a": RuleKind.PREFERENCE,
    })
    assert custom.is_custom
    assert "pos3-a" in custom
    assert "pos19-a" not in custom
    exported = custom.to_dict()
    assert exported["baseline_rules"] == before["rules"]
    assert exported["disabled_rules"][0]["provenance"] == CATALOGUE["pos19-a"].to_dict()[
        "provenance"
    ]
    assert {change["rule_id"] for change in custom.customizations} == {
        "tss-window", "pos19-a", "pos3-a",
    }
    assert WANG_2015.resolve().to_dict() == before
    assert WANG_2015.resolve().rule("tss-window").kind is RuleKind.MANDATORY
    assert not WANG_2015.resolve().is_custom


def test_preset_values_are_the_baseline_for_custom_parameters() -> None:
    config = STRICT.resolve({"tss-window": {"window-start-bp": -600}})
    assert list(config.customizations) == [{
        "rule_id": "tss-window", "setting": "window-start-bp", "before": -500, "after": -600,
    }]
    assert not STRICT.resolve().is_custom


def test_overlapping_window_preferences_are_explained_and_can_be_avoided() -> None:
    two = WANG_2015.resolve(role_overrides={"tss-window": RuleKind.PREFERENCE})
    one = WANG_2015.resolve(role_overrides={
        "tss-window": RuleKind.PREFERENCE, "tss-preferred-window": None,
    })
    assert two.max_preference_points == 8
    assert two.warnings
    assert two.to_dict()["warnings"] == list(two.warnings)
    assert one.max_preference_points == 7
    assert not one.warnings
    assert not WANG_2015.resolve().warnings


def test_explicit_option_is_distinct_from_manual_matching_values() -> None:
    manual = ResolvedRule(CATALOGUE["gc-content"], {"maximum-gc-fraction": 0.65}).to_dict()
    selected = ResolvedRule.from_option(CATALOGUE["gc-content"], "legacy-scanner").to_dict()
    assert manual["selected_option"] is None
    assert selected["selected_option"]["id"] == "legacy-scanner"
    assert selected["parameters"][1]["origin"] == "named-option"
    assert manual["parameters"][1]["origin"] == "override"


def test_named_option_replaces_conflicting_preset_value_before_manual_override() -> None:
    config = STRICT.resolve(option_ids={"tss-window": "wang-2015"})
    assert config.rule("tss-window").integer("window-start-bp") == -1000
    adjusted = STRICT.resolve(
        {"tss-window": {"window-start-bp": -900}},
        option_ids={"tss-window": "wang-2015"},
    )
    assert adjusted.rule("tss-window").integer("window-start-bp") == -900
    assert adjusted.rule("tss-window").option_id == "wang-2015"


@pytest.mark.parametrize("role", [None, RuleKind.PREFERENCE, RuleKind.TIEBREAK])
def test_generation_length_cannot_be_disabled_or_scored(role: RuleKind | None) -> None:
    with pytest.raises(ValueError):
        WANG_2015.resolve(role_overrides={"target-length-19": role})


@pytest.mark.parametrize("role", [RuleKind.MANDATORY, RuleKind.PREFERENCE])
def test_continuous_tiebreak_cannot_be_given_a_binary_role(role: RuleKind) -> None:
    with pytest.raises(ValueError, match="continuous tie-break"):
        WANG_2015.resolve(role_overrides={"asymmetry-margin": role})


def test_unknown_or_ineffective_changes_are_rejected() -> None:
    with pytest.raises(KeyError, match="unknown-rule"):
        WANG_2015.resolve(role_overrides={"unknown-rule": None})
    with pytest.raises(KeyError, match="not selected"):
        DesignConfig.of(["target-length-19"], {"gc-content": {"maximum-gc-fraction": 0.65}})
    with pytest.raises(KeyError, match="disabled"):
        WANG_2015.resolve(
            {"gc-content": {"maximum-gc-fraction": 0.65}}, role_overrides={"gc-content": None},
        )
    with pytest.raises(ValueError, match="RuleKind"):
        WANG_2015.resolve(role_overrides={"pos19-a": "preference"})  # type: ignore[dict-item]


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_nonfinite_settings_are_rejected_before_evaluation(value: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        WANG_2015.resolve({"gc-content": {"maximum-gc-fraction": value}})


def test_inverted_ranges_and_fractional_counts_are_rejected_at_resolution() -> None:
    with pytest.raises(ValueError, match="must not exceed"):
        WANG_2015.resolve({"gc-content": {"minimum-gc-fraction": 0.7}})
    with pytest.raises(ValueError, match="whole number"):
        WANG_2015.resolve({"cpg-count": {"maximum-cpg-count": 0.5}})


def test_frozen_html_contains_effective_and_baseline_roles_and_choices() -> None:
    config = WANG_2015.resolve(
        role_overrides={"pos19-a": None, "tss-window": RuleKind.PREFERENCE},
        option_ids={"gc-content": "legacy-scanner"},
    )
    run = design_sarnas(three_windows(), config)
    page = render(PromoterView(result=run))
    match = re.search(r'<script id="view-data" type="application/json">(.*?)</script>', page, re.S)
    assert match is not None
    exported = json.loads(match.group(1))["configuration"]
    assert exported["is_custom"]
    assert exported["disabled_rules"][0]["id"] == "pos19-a"
    window = next(r for r in exported["rules"] if r["id"] == "tss-window")
    assert window["kind"] == "preference"
    assert window["default_kind"] == "mandatory"
    gc = next(r for r in exported["rules"] if r["id"] == "gc-content")
    assert gc["selected_option"]["id"] == "legacy-scanner"
