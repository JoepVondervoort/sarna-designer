"""Tests for what a rule reports after it has been evaluated.

Every claim the UI and the exports make about a candidate passes through this type, so the
tests here are mostly about what it *refuses* to represent: a verdict with no stated reason, a
not-applicable result that silently counts as a pass, a tie-break folded into the score.
"""

import pytest
from pydantic import ValidationError

from sarna.rules.catalog.mandatory import GC_CONTENT, TSS_WINDOW
from sarna.rules.catalog.preferences import ASYMMETRY_MARGIN, PREFERENCE_RULES
from sarna.rules.results import RuleOutcome, RuleResult
from sarna.rules.rule import RuleKind

A_PREFERENCE = PREFERENCE_RULES[0]


class TestConstruction:
    def test_passed_records_the_rule_and_outcome(self) -> None:
        result = RuleResult.passed(GC_CONTENT, "GC content 47.4 %, inside the 40-60 % band")
        assert result.rule_id == GC_CONTENT.id
        assert result.kind is RuleKind.MANDATORY
        assert result.outcome is RuleOutcome.PASS

    def test_failed_records_the_outcome(self) -> None:
        result = RuleResult.failed(GC_CONTENT, "GC content 68.4 %, above the 60 % maximum")
        assert result.outcome is RuleOutcome.FAIL

    def test_not_applicable_records_the_outcome(self) -> None:
        result = RuleResult.not_applicable(
            TSS_WINDOW, "no TSS is known for a pasted promoter sequence"
        )
        assert result.outcome is RuleOutcome.NOT_APPLICABLE

    def test_measurements_are_carried_alongside_the_prose(self) -> None:
        # The prose is for a human; the measurements are what a CSV column or a plot reads,
        # so neither has to be recovered by parsing the other.
        result = RuleResult.passed(
            GC_CONTENT,
            "GC content 47.4 %, inside the 40-60 % band",
            measurements={"gc-fraction": 0.4736},
        )
        assert result.measurements == {"gc-fraction": 0.4736}

    def test_measurements_default_to_empty(self) -> None:
        assert RuleResult.passed(GC_CONTENT, "inside the band").measurements == {}

    def test_the_kind_is_taken_from_the_rule_not_the_caller(self) -> None:
        assert RuleResult.passed(A_PREFERENCE, "position 1 is G").kind is RuleKind.PREFERENCE


class TestDetailIsMandatory:
    @pytest.mark.parametrize("detail", ["", "   ", "\n"])
    def test_a_verdict_without_a_reason_is_rejected(self, detail: str) -> None:
        # A rule that says "fail" without saying why is the thing this project exists to stop
        # shipping (ADR 0003).
        with pytest.raises(ValidationError):
            RuleResult.failed(GC_CONTENT, detail)

    def test_a_pass_also_has_to_say_what_it_measured(self) -> None:
        with pytest.raises(ValidationError):
            RuleResult.passed(GC_CONTENT, "")

    def test_not_applicable_has_to_say_why_it_could_not_be_evaluated(self) -> None:
        with pytest.raises(ValidationError):
            RuleResult.not_applicable(TSS_WINDOW, "  ")


class TestExclusion:
    def test_a_failed_mandatory_rule_excludes_the_candidate(self) -> None:
        assert RuleResult.failed(GC_CONTENT, "GC content 68.4 %").excludes is True

    def test_a_passed_mandatory_rule_does_not_exclude(self) -> None:
        assert RuleResult.passed(GC_CONTENT, "GC content 47.4 %").excludes is False

    def test_an_unmet_preference_does_not_exclude(self) -> None:
        assert RuleResult.failed(A_PREFERENCE, "position 1 is A").excludes is False

    def test_a_mandatory_rule_that_could_not_be_evaluated_does_not_exclude(self) -> None:
        # ADR 0007: "we could not check this" must not become "this failed". A pasted promoter
        # with no coordinates would otherwise have every candidate excluded by tss-window.
        result = RuleResult.not_applicable(TSS_WINDOW, "no TSS is known")
        assert result.excludes is False
        assert result.is_unevaluated is True


class TestPoints:
    def test_a_met_preference_is_worth_one_point(self) -> None:
        assert RuleResult.passed(A_PREFERENCE, "position 1 is G").points == 1

    def test_an_unmet_preference_is_worth_nothing(self) -> None:
        assert RuleResult.failed(A_PREFERENCE, "position 1 is A").points == 0

    def test_a_passed_mandatory_rule_scores_nothing(self) -> None:
        # ADR 0006: an admission requirement is not a preference you can trade off, so passing
        # one must not inflate the score.
        assert RuleResult.passed(GC_CONTENT, "GC content 47.4 %").points == 0

    def test_an_unevaluated_preference_scores_nothing(self) -> None:
        result = RuleResult.not_applicable(A_PREFERENCE, "sequence too short")
        assert result.points == 0

    def test_a_tiebreak_scores_nothing(self) -> None:
        # ADR 0006: tie-breaks are reported as tie-breaks, never folded into the score.
        result = RuleResult.tiebreak(ASYMMETRY_MARGIN, 1.4, "asymmetry 1.4 kcal/mol")
        assert result.points == 0


class TestTiebreaks:
    def test_a_tiebreak_carries_its_ordering_value(self) -> None:
        result = RuleResult.tiebreak(ASYMMETRY_MARGIN, 1.4, "asymmetry 1.4 kcal/mol")
        assert result.ordering_value == pytest.approx(1.4)
        assert result.kind is RuleKind.TIEBREAK

    def test_a_tiebreak_without_an_ordering_value_is_rejected(self) -> None:
        with pytest.raises(ValidationError, match="ordering value"):
            RuleResult(
                rule_id=ASYMMETRY_MARGIN.id,
                kind=RuleKind.TIEBREAK,
                outcome=RuleOutcome.PASS,
                detail="asymmetry 1.4 kcal/mol",
            )

    def test_a_tiebreak_that_could_not_be_computed_needs_no_ordering_value(self) -> None:
        result = RuleResult.not_applicable(ASYMMETRY_MARGIN, "no thermodynamic table selected")
        assert result.ordering_value is None

    def test_a_non_tiebreak_rule_must_not_carry_an_ordering_value(self) -> None:
        # Allowing one would reintroduce the legacy weighted sum by the back door.
        with pytest.raises(ValidationError, match="ordering value"):
            RuleResult(
                rule_id=GC_CONTENT.id,
                kind=RuleKind.MANDATORY,
                outcome=RuleOutcome.PASS,
                detail="GC content 47.4 %",
                ordering_value=1.4,
            )


class TestValueSemantics:
    def test_is_immutable(self) -> None:
        result = RuleResult.passed(GC_CONTENT, "GC content 47.4 %")
        with pytest.raises(ValidationError):
            result.outcome = RuleOutcome.FAIL

    def test_outcome_renders_as_a_stable_label(self) -> None:
        assert str(RuleOutcome.NOT_APPLICABLE) == "not-applicable"
        assert str(RuleOutcome.PASS) == "pass"
        assert str(RuleOutcome.FAIL) == "fail"
