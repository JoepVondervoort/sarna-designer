"""A design rule is a record, not a predicate.

ADR 0003 fixes the minimum it must carry. ADR 0006 fixes what each kind contributes to a
ranking: mandatory rules gate, preference rules are worth exactly one point each, and
tie-breaks order without scoring.
"""

import pytest
from pydantic import ValidationError

from sarna.rules.provenance import EvidenceClass, Provenance
from sarna.rules.rule import Rule, RuleKind

WANG_TABLE_1 = Provenance.citation("Wang et al. 2015", "Table 1")


def make_rule(**overrides: object) -> Rule:
    defaults: dict[str, object] = {
        "id": "gc-content",
        "name": "GC content 40-60 %",
        "kind": RuleKind.MANDATORY,
        "evidence_class": EvidenceClass.EMPIRICAL,
        "provenance": WANG_TABLE_1,
        "rationale": "Duplex stability and RISC loading degrade outside this band.",
        "consequence": "Tightening this excludes GC-rich candidates near CpG islands.",
    }
    return Rule(**{**defaults, **overrides})  # type: ignore[arg-type]


class TestConstruction:
    def test_builds_a_complete_rule(self):
        rule = make_rule()
        assert rule.id == "gc-content"
        assert rule.kind is RuleKind.MANDATORY
        assert rule.provenance == WANG_TABLE_1

    def test_rules_are_enabled_by_default(self):
        assert make_rule().enabled_by_default is True

    def test_a_rule_can_ship_disabled(self):
        # The position-7 rule is carried for traceability but off by default (ADR 0002).
        assert make_rule(enabled_by_default=False).enabled_by_default is False

    def test_grouping_is_optional(self):
        assert make_rule().group is None
        assert make_rule(group="3' flank A/U richness").group == "3' flank A/U richness"


class TestRequiredMetadata:
    @pytest.mark.parametrize("field", ["rationale", "consequence"])
    def test_empty_explanatory_text_is_refused(self, field: str):
        with pytest.raises(ValidationError, match=field):
            make_rule(**{field: "   "})

    def test_provenance_cannot_be_omitted(self):
        # ADR 0003: no rule may be added without provenance.
        with pytest.raises(ValidationError, match="provenance"):
            Rule(  # type: ignore[call-arg]
                id="x",
                name="x",
                kind=RuleKind.MANDATORY,
                evidence_class=EvidenceClass.EMPIRICAL,
                rationale="r",
                consequence="c",
            )

    @pytest.mark.parametrize("bad_id", ["GC_content", "gc content", "gc--content", "-gc", ""])
    def test_ids_must_be_kebab_case(self, bad_id: str):
        with pytest.raises(ValidationError, match="kebab-case"):
            make_rule(id=bad_id)

    @pytest.mark.parametrize("good_id", ["gc-content", "pos19-a", "flank-20-au", "tss-window"])
    def test_valid_kebab_case_ids_are_accepted(self, good_id: str):
        assert make_rule(id=good_id).id == good_id


class TestScientificIntegrity:
    def test_an_unsourced_rule_cannot_claim_to_be_established(self):
        # This is the invariant that would have caught the position-7 rule.
        with pytest.raises(ValidationError, match="established"):
            make_rule(
                provenance=Provenance.unsourced("Position 7 appears nowhere in Wang 2015."),
                evidence_class=EvidenceClass.ESTABLISHED,
            )

    def test_an_unsourced_rule_may_be_exploratory(self):
        rule = make_rule(
            provenance=Provenance.unsourced("Position 7 appears nowhere in Wang 2015."),
            evidence_class=EvidenceClass.EXPLORATORY,
        )
        assert rule.needs_warning is True

    def test_a_cited_rule_needs_no_warning(self):
        assert make_rule().needs_warning is False


class TestPointsContributed:
    def test_a_preference_is_worth_exactly_one_point(self):
        # ADR 0006: all preferences are weighted equally. There is no other value.
        assert make_rule(kind=RuleKind.PREFERENCE).points == 1

    def test_a_mandatory_rule_scores_nothing_because_it_gates(self):
        assert make_rule(kind=RuleKind.MANDATORY).points == 0

    def test_a_tiebreak_scores_nothing(self):
        assert make_rule(kind=RuleKind.TIEBREAK).points == 0

    def test_gating_and_scoring_are_mutually_exclusive(self):
        assert make_rule(kind=RuleKind.MANDATORY).is_gate is True
        assert make_rule(kind=RuleKind.PREFERENCE).is_gate is False


class TestRuleIsData:
    def test_is_immutable(self):
        with pytest.raises(ValidationError):
            make_rule().id = "other"

    def test_round_trips_through_serialisation(self):
        # ADR 0003: exports carry the rule records, not just the resolved values.
        rule = make_rule()
        assert Rule.model_validate(rule.model_dump()) == rule

    def test_serialisation_includes_the_provenance_and_evidence_class(self):
        dumped = make_rule().model_dump(mode="json")
        assert dumped["evidence_class"] == "empirical"
        assert dumped["provenance"]["source"] == "Wang et al. 2015"
        assert dumped["provenance"]["locator"] == "Table 1"


class TestRuleKind:
    def test_renders_as_its_label(self):
        assert str(RuleKind.MANDATORY) == "mandatory"
        assert str(RuleKind.PREFERENCE) == "preference"
        assert str(RuleKind.TIEBREAK) == "tiebreak"


class TestRuleParameters:
    def test_a_rule_may_expose_no_parameters(self):
        assert make_rule().parameters == ()

    def test_parameters_are_reachable_by_name(self):
        from sarna.rules.parameters import RuleParameter

        bound = RuleParameter(
            name="minimum-gc-fraction",
            default=0.40,
            minimum=0.0,
            maximum=1.0,
            unit="fraction",
            description="Lowest GC fraction a target may have.",
        )
        rule = make_rule(parameters=(bound,))
        assert rule.parameter("minimum-gc-fraction").default == 0.40

    def test_an_unknown_parameter_name_raises(self):
        with pytest.raises(KeyError, match="no-such-parameter"):
            make_rule().parameter("no-such-parameter")
