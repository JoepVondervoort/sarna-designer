"""Tests for the optional siRNA-literature and chapter 4 evaluators.

The positional rules are thin, so the substance here is the window arithmetic: the A/T-rich rule
counts from the target's **3' end** and the central-stability rule counts from its 5' end, and
confusing the two would silently score the wrong end of the molecule. Each is therefore checked
against a sequence built so that the right answer and the wrong answer differ.
"""

import pytest

from sarna.design.config import DesignConfig
from sarna.design.context import CandidateContext
from sarna.design.evaluators import evaluate
from sarna.genome.promoter import Promoter
from sarna.rules.results import RuleOutcome
from sarna.sequence.nucleic import NucleicSequence

#                       1234567890123456789
# positions:            1         10      19
AT_RICH_3_PRIME = "GCGCGCGCGCGCATATAAA"
"""G/C through position 12, then A/T from 13 to 19. Seven weak bases at the 3' end."""

GC_RICH_3_PRIME = "ATATATATATATGCGCGCG"
"""The mirror image: weak bases at the 5' end, strong at the 3'. Fails the same rule."""

STRONG_CENTRE = "ATATATATAGCGCATATAT"
"""G/C at positions 10-13 and nowhere else, so the central region is the only strong part."""


def context_for(target: str, offset: int = 0) -> CandidateContext:
    """Build a context whose target is exactly ``target``, padded so flanks exist."""
    promoter = Promoter(sequence=NucleicSequence.parse("A" * offset + target + "TTTT"))
    return CandidateContext.at(promoter, offset=offset, length=len(target))


def result_for(rule_id: str, target: str, **overrides: float | int | str):
    config = DesignConfig.of([rule_id], {rule_id: overrides} if overrides else None)
    return evaluate(context_for(target), config.rule(rule_id))


class TestAuRichGuide5Prime:
    def test_a_weak_3_prime_end_passes(self) -> None:
        result = result_for("au-rich-guide-5-prime", AT_RICH_3_PRIME)
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["au-count"] == 7
        assert result.measurements["window-bases"] == "ATATAAA"

    def test_the_window_is_counted_from_the_3_prime_end(self) -> None:
        # The decisive test: the same bases arranged the other way round must fail. If the
        # window were read from the 5' end, GC_RICH_3_PRIME would pass instead.
        result = result_for("au-rich-guide-5-prime", GC_RICH_3_PRIME)
        assert result.outcome is RuleOutcome.FAIL
        assert result.measurements["au-count"] == 0

    def test_the_window_spans_the_last_seven_positions_by_default(self) -> None:
        result = result_for("au-rich-guide-5-prime", AT_RICH_3_PRIME)
        assert result.measurements["window-start"] == 13
        assert result.measurements["window-end"] == 19

    def test_the_verdict_explains_the_strand_relationship(self) -> None:
        # A reader must be told that a weak target 3' end is a weak guide 5' end, since that
        # is the whole mechanism and the positions run in opposite directions.
        detail = result_for("au-rich-guide-5-prime", AT_RICH_3_PRIME).detail
        assert "3' end on the target" in detail
        assert "5' end on the guide" in detail

    def test_a_shorter_window_can_be_requested(self) -> None:
        result = result_for("au-rich-guide-5-prime", AT_RICH_3_PRIME, **{"window-length-nt": 3})
        assert result.measurements["window-start"] == 17
        assert result.measurements["au-count"] == 3

    def test_the_reynolds_threshold_is_more_permissive_than_ui_teis(self) -> None:
        # Three of the last five, rather than five of the last seven.
        borderline = "GCGCGCGCGCGCGCGATAA"
        assert result_for("au-rich-guide-5-prime", borderline).outcome is RuleOutcome.FAIL
        relaxed = result_for(
            "au-rich-guide-5-prime",
            borderline,
            **{"window-length-nt": 5, "minimum-au-count": 3},
        )
        assert relaxed.outcome is RuleOutcome.PASS

    def test_a_window_longer_than_the_target_is_unevaluable(self) -> None:
        # Not a failure: a window that does not fit was never measured.
        result = result_for("au-rich-guide-5-prime", "GCGCATATA", **{"window-length-nt": 10})
        assert result.outcome is RuleOutcome.NOT_APPLICABLE
        assert result.points == 0


class TestNoGcStretch:
    def test_a_target_within_the_limit_passes(self) -> None:
        result = result_for("no-gc-stretch", "ATATGCGCGCATATATATA")
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["longest-gc-stretch"] == 6

    def test_a_long_gc_tract_fails(self) -> None:
        result = result_for("no-gc-stretch", "AGCGCGCGCGCGCGATATA")
        assert result.outcome is RuleOutcome.FAIL
        assert result.measurements["longest-gc-stretch"] == 13

    def test_it_counts_g_and_c_together_not_repeats_of_one_base(self) -> None:
        # This is what distinguishes it from no-homopolymer-run: ten alternating G/C bases are
        # one G/C stretch and no homopolymer at all.
        result = result_for("no-gc-stretch", "AGCGCGCGCGCATATATA", **{"maximum-gc-stretch-nt": 5})
        assert result.outcome is RuleOutcome.FAIL
        assert result.measurements["longest-gc-stretch"] == 10

    def test_a_target_with_no_gc_at_all_passes(self) -> None:
        result = result_for("no-gc-stretch", "ATATATATATATATATATA")
        assert result.measurements["longest-gc-stretch"] == 0

    def test_the_limit_is_configurable(self) -> None:
        target = "ATATGCGCGCATATATATA"
        assert result_for("no-gc-stretch", target, **{"maximum-gc-stretch-nt": 5}).outcome is (
            RuleOutcome.FAIL
        )


class TestCentralStability:
    def test_a_strong_central_region_passes(self) -> None:
        # Positions 10-13 of this target are GCGC.
        result = result_for("central-stability", STRONG_CENTRE)
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["region-bases"] == "GCGC"
        assert result.measurements["strong-count"] == 4

    def test_a_weak_central_region_fails(self) -> None:
        result = result_for("central-stability", "GCGCGCGCGATATGCGCGC")
        assert result.outcome is RuleOutcome.FAIL
        assert result.measurements["region-bases"] == "ATAT"

    def test_the_region_is_counted_from_the_5_prime_end(self) -> None:
        # Mirror of the A/T-rich test: this rule's positions are literal, not from the 3' end.
        result = result_for("central-stability", STRONG_CENTRE)
        assert result.measurements["region-start"] == 10
        assert result.measurements["region-end"] == 13

    def test_the_region_is_configurable(self) -> None:
        result = result_for(
            "central-stability",
            "GCGCATATATATATATATA",
            **{"central-region-start": 1, "central-region-end": 4},
        )
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["region-bases"] == "GCGC"

    def test_three_of_four_strong_bases_is_enough_by_default(self) -> None:
        result = result_for("central-stability", "ATATATATAGCGAATATATA"[:19])
        assert result.measurements["strong-count"] == 3
        assert result.outcome is RuleOutcome.PASS

    def test_a_region_past_the_end_of_a_short_target_is_unevaluable(self) -> None:
        result = result_for("central-stability", "ATATATATAG")
        assert result.outcome is RuleOutcome.NOT_APPLICABLE

    def test_an_inverted_region_is_a_configuration_error(self) -> None:
        # Silently passing or failing every candidate would be worse than refusing: an empty
        # region is a setting mistake, not a property of any candidate.
        with pytest.raises(ValueError, match="must not exceed"):
            result_for(
                "central-stability",
                STRONG_CENTRE,
                **{"central-region-start": 13, "central-region-end": 10},
            )


class TestThePositionalOptionalRules:
    @pytest.mark.parametrize(
        ("rule_id", "passing", "failing"),
        [
            ("pos3-a", "GCATATATATATATATATA", "GCGTATATATATATATATA"),
            ("pos10-u", "GCGCGCGCGTATATATATA", "GCGCGCGCGGATATATATA"),
            ("pos10-g", "GCGCGCGCGGATATATATA", "GCGCGCGCGTATATATATA"),
            ("no-g-at-13", "GCGCGCGCGCGCATATATA", "GCGCGCGCGCGCGTATATA"),
        ],
    )
    def test_each_rule_passes_and_fails_the_expected_target(
        self, rule_id: str, passing: str, failing: str
    ) -> None:
        assert result_for(rule_id, passing).outcome is RuleOutcome.PASS
        assert result_for(rule_id, failing).outcome is RuleOutcome.FAIL

    def test_the_two_position_10_rules_are_mutually_exclusive(self) -> None:
        # Stated in both rules' consequence text; asserted here so it stays true.
        for target in ("GCGCGCGCGGATATATATA", "GCGCGCGCGTATATATATA"):
            outcomes = {
                result_for(rule_id, target).outcome for rule_id in ("pos10-u", "pos10-g")
            }
            assert outcomes == {RuleOutcome.PASS, RuleOutcome.FAIL}

    def test_no_g_at_13_is_satisfied_by_any_of_three_bases(self) -> None:
        for base in "ACT":
            target = "GCGCGCGCGCGC" + base + "ATATAT"
            assert result_for("no-g-at-13", target).outcome is RuleOutcome.PASS

    def test_a_position_beyond_a_short_target_is_unevaluable(self) -> None:
        assert result_for("no-g-at-13", "GCGCGC").outcome is RuleOutcome.NOT_APPLICABLE
