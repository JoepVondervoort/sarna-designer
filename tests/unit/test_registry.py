"""The registry is an immutable catalogue of rules, keyed by id.

Presets are subsets of it, so the registry itself stays a plain collection with no notion of
defaults or configuration.
"""

import pytest

from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.registry import DuplicateRuleError, RuleNotFoundError, RuleRegistry
from sarna.rules.rule import Rule, RuleKind

WANG_TABLE_1 = Provenance.citation("Wang et al. 2015", "Table 1")


def rule(rule_id: str, kind: RuleKind = RuleKind.PREFERENCE, **overrides: object) -> Rule:
    defaults: dict[str, object] = {
        "id": rule_id,
        "name": rule_id.replace("-", " "),
        "kind": kind,
        "evidence_class": EvidenceClass.EMPIRICAL,
        "provenance": WANG_TABLE_1,
        "rationale": "because.",
        "consequence": "changes the candidate set.",
    }
    return Rule(**{**defaults, **overrides})  # type: ignore[arg-type]


@pytest.fixture
def registry() -> RuleRegistry:
    return RuleRegistry(
        [
            rule("gc-content", RuleKind.MANDATORY),
            rule("pos19-a", RuleKind.MANDATORY),
            rule("start-1-gc"),
            rule("start-2-gc"),
            rule("asymmetry-margin", RuleKind.TIEBREAK),
        ]
    )


class TestConstruction:
    def test_holds_every_rule_given(self, registry: RuleRegistry):
        assert len(registry) == 5

    def test_is_empty_when_given_nothing(self):
        assert len(RuleRegistry([])) == 0

    def test_duplicate_ids_are_refused(self):
        with pytest.raises(DuplicateRuleError, match="gc-content"):
            RuleRegistry([rule("gc-content"), rule("gc-content")])

    def test_preserves_the_order_rules_were_given_in(self, registry: RuleRegistry):
        # Rule order is the order the UI presents them in, so it must be deterministic.
        assert [r.id for r in registry] == [
            "gc-content",
            "pos19-a",
            "start-1-gc",
            "start-2-gc",
            "asymmetry-margin",
        ]


class TestLookup:
    def test_finds_a_rule_by_id(self, registry: RuleRegistry):
        assert registry["pos19-a"].id == "pos19-a"

    def test_unknown_id_raises_and_names_the_id(self, registry: RuleRegistry):
        with pytest.raises(RuleNotFoundError, match="no-such-rule"):
            registry["no-such-rule"]

    def test_membership_is_testable(self, registry: RuleRegistry):
        assert "gc-content" in registry
        assert "no-such-rule" not in registry


class TestPartitions:
    def test_separates_gates_from_preferences_from_tiebreaks(self, registry: RuleRegistry):
        assert [r.id for r in registry.mandatory] == ["gc-content", "pos19-a"]
        assert [r.id for r in registry.preferences] == ["start-1-gc", "start-2-gc"]
        assert [r.id for r in registry.tiebreaks] == ["asymmetry-margin"]

    def test_the_partitions_cover_every_rule_exactly_once(self, registry: RuleRegistry):
        partitioned = [*registry.mandatory, *registry.preferences, *registry.tiebreaks]
        assert sorted(r.id for r in partitioned) == sorted(r.id for r in registry)

    def test_the_maximum_preference_score_is_the_number_of_preferences(
        self, registry: RuleRegistry
    ):
        # This is the denominator in "meets 5 of 7 preferences".
        assert registry.max_preference_points == 2


class TestImmutability:
    def test_adding_a_rule_returns_a_new_registry(self, registry: RuleRegistry):
        extended = registry.with_rule(rule("flank-20-au"))
        assert len(extended) == 6
        assert len(registry) == 5, "the original registry must be untouched"

    def test_adding_a_duplicate_is_refused(self, registry: RuleRegistry):
        with pytest.raises(DuplicateRuleError, match="gc-content"):
            registry.with_rule(rule("gc-content"))

    def test_selecting_a_subset_returns_a_new_registry(self, registry: RuleRegistry):
        subset = registry.select(["gc-content", "start-1-gc"])
        assert [r.id for r in subset] == ["gc-content", "start-1-gc"]
        assert len(registry) == 5

    def test_selecting_an_unknown_id_raises(self, registry: RuleRegistry):
        with pytest.raises(RuleNotFoundError, match="no-such-rule"):
            registry.select(["gc-content", "no-such-rule"])


class TestWarnings:
    def test_surfaces_rules_whose_provenance_needs_a_warning(self):
        unsourced = rule(
            "pos7-u",
            evidence_class=EvidenceClass.EXPLORATORY,
            provenance=Provenance.unsourced("Position 7 appears nowhere in Wang 2015."),
        )
        registry = RuleRegistry([rule("gc-content", RuleKind.MANDATORY), unsourced])
        assert [r.id for r in registry.needing_warning] == ["pos7-u"]
