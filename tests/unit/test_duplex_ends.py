"""Tests for terminal end stability and guide-strand prediction (ADR 0009).

Two things are pinned here. First, the arithmetic must reproduce the legacy ``endstability``
exactly when the legacy table is selected, because that function's output is the scientific
ground truth for the ``legacy-scanner`` preset (ADR 0001). The expected values below are
computed by hand from the script's own constants rather than from this implementation.

Second, every prediction must explain itself in words naming the mechanism, since that text is
engine output and not UI decoration.
"""

import pytest

from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.ends import DuplexEnd, end_stability
from sarna.thermo.guide import GuideStrand, assess_asymmetry, compare_asymmetry
from sarna.thermo.tables import DNA_SANTALUCIA_1998, LEGACY_DNA_KJ, RNA_XIA_1998

# Legacy constants, copied from sarna_promoter_scanner.py lines 23-32 (kJ/mol).
AAE, ATE, TAE, CAE = -4.26, -3.67, -2.50, -6.12
GTE, CTE, GAE, CGE, GCE, GGE = -6.09, -5.40, -5.51, -9.07, -9.36, -7.66

# A 19-mer whose two ends differ clearly, so the sign of the difference is unambiguous.
# 5' end GCGC, 3' end TTAA.
TARGET = NucleicSequence.parse("GCGCACGGATCCAAGTTAA")


class TestEndStability:
    def test_sums_three_stacks_over_a_four_base_end(self) -> None:
        # end = 4 in the legacy script, sourced in its own comment to the helicase footprint
        # reported by Schwarz et al. 2003: three stacks span four bases.
        end = end_stability(TARGET, LEGACY_DNA_KJ, DuplexEnd.FIVE_PRIME, end_length=4)
        assert len(end.stacks) == 3
        assert [step for step, _ in end.stacks] == ["GC", "CG", "GC"]

    def test_the_five_prime_energy_matches_the_legacy_constants(self) -> None:
        end = end_stability(TARGET, LEGACY_DNA_KJ, DuplexEnd.FIVE_PRIME, end_length=4)
        assert end.energy == pytest.approx(GCE + CGE + GCE)

    def test_the_three_prime_energy_matches_the_legacy_constants(self) -> None:
        # Last four bases are TTAA, so the stacks are TT, TA, AA. The legacy table assigns TT
        # the AA value, which is correct: AA and TT are reverse-complement equivalent.
        end = end_stability(TARGET, LEGACY_DNA_KJ, DuplexEnd.THREE_PRIME, end_length=4)
        assert [step for step, _ in end.stacks] == ["TT", "TA", "AA"]
        assert end.energy == pytest.approx(AAE + TAE + AAE)

    def test_the_end_renders_as_a_prime_label(self) -> None:
        assert str(DuplexEnd.FIVE_PRIME) == "5'"
        assert str(DuplexEnd.THREE_PRIME) == "3'"

    def test_the_unit_comes_from_the_table(self) -> None:
        end = end_stability(TARGET, LEGACY_DNA_KJ, DuplexEnd.FIVE_PRIME, end_length=4)
        assert end.unit == "kJ/mol"
        rna = end_stability(TARGET, RNA_XIA_1998, DuplexEnd.FIVE_PRIME, end_length=4)
        assert rna.unit == "kcal/mol"

    def test_an_rna_table_is_applied_to_the_transcribed_target(self) -> None:
        # The target site is DNA; an RNA parameter set describes the duplex, so the sequence is
        # transcribed before lookup. T/U is a one-to-one relabelling, so the stacks correspond.
        end = end_stability(TARGET, RNA_XIA_1998, DuplexEnd.THREE_PRIME, end_length=4)
        assert [step for step, _ in end.stacks] == ["UU", "UA", "AA"]

    def test_a_dna_table_is_applied_to_the_dna_target(self) -> None:
        end = end_stability(TARGET, DNA_SANTALUCIA_1998, DuplexEnd.THREE_PRIME, end_length=4)
        assert [step for step, _ in end.stacks] == ["TT", "TA", "AA"]

    @pytest.mark.parametrize("end_length", [2, 3, 4, 19])
    def test_the_end_length_is_configurable(self, end_length: int) -> None:
        end = end_stability(TARGET, RNA_XIA_1998, DuplexEnd.FIVE_PRIME, end_length=end_length)
        assert len(end.stacks) == end_length - 1

    @pytest.mark.parametrize("end_length", [0, 1, 20])
    def test_rejects_an_end_length_outside_the_target(self, end_length: int) -> None:
        with pytest.raises(ValueError):
            end_stability(TARGET, RNA_XIA_1998, DuplexEnd.FIVE_PRIME, end_length=end_length)


class TestAsymmetry:
    def test_the_difference_is_three_prime_minus_five_prime(self) -> None:
        # The legacy sign convention, preserved: positive means the target's 5' end is the more
        # stable, which is the desirable condition.
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        expected = (AAE + TAE + AAE) - (GCE + CGE + GCE)
        assert assessment.difference == pytest.approx(expected)
        assert assessment.difference > 0

    def test_a_stable_five_prime_end_predicts_the_antisense_guide(self) -> None:
        # The antisense strand's 5' end pairs with the target's 3' end. Here the target's 3' end
        # is the weaker, so the antisense 5' end is the weaker and RISC should load it.
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert assessment.predicted_guide is GuideStrand.ANTISENSE

    def test_a_stable_three_prime_end_predicts_the_sense_guide(self) -> None:
        reversed_target = NucleicSequence.parse("TTAAACGGATCCAAGGCGC")
        assessment = assess_asymmetry(reversed_target, LEGACY_DNA_KJ, end_length=4)
        assert assessment.difference < 0
        assert assessment.predicted_guide is GuideStrand.SENSE

    def test_symmetrical_ends_are_ambiguous_rather_than_guessed(self) -> None:
        # Both ends are GCGC, so the difference is exactly zero. Picking a strand here would be
        # a coin flip presented as a prediction.
        symmetrical = NucleicSequence.parse("GCGCACGGATCCAAGGCGC")
        assessment = assess_asymmetry(symmetrical, LEGACY_DNA_KJ, end_length=4)
        assert assessment.difference == pytest.approx(0.0)
        assert assessment.predicted_guide is GuideStrand.AMBIGUOUS

    def test_the_legacy_comparison_tolerance_is_the_default(self) -> None:
        # The legacy script compares against `fudge = 1e-8`. Retained as the default so the
        # legacy preset's boundary behaviour is reproduced, but configurable rather than buried.
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert assessment.tolerance == pytest.approx(1e-8)

    def test_a_difference_inside_a_widened_tolerance_is_ambiguous(self) -> None:
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4, tolerance=100.0)
        assert assessment.predicted_guide is GuideStrand.AMBIGUOUS

    def test_rejects_a_negative_tolerance(self) -> None:
        with pytest.raises(ValueError, match="tolerance"):
            assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4, tolerance=-1.0)

    def test_the_assessment_names_its_table(self) -> None:
        assessment = assess_asymmetry(TARGET, RNA_XIA_1998, end_length=4)
        assert assessment.table_id == "rna-xia-1998"
        assert assessment.table_name == RNA_XIA_1998.name


class TestExplanation:
    """ADR 0009: a prediction must say why, in words, naming the mechanism."""

    def test_the_explanation_names_the_predicted_strand(self) -> None:
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert "antisense" in assessment.explanation

    def test_the_explanation_quotes_both_end_energies_and_the_unit(self) -> None:
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert "kJ/mol" in assessment.explanation
        assert f"{assessment.five_prime.energy:.2f}" in assessment.explanation
        assert f"{assessment.three_prime.energy:.2f}" in assessment.explanation

    def test_the_explanation_cites_the_mechanism(self) -> None:
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert "Schwarz" in assessment.explanation

    def test_a_sense_prediction_explains_that_it_is_the_undesirable_outcome(self) -> None:
        # The researcher needs to know not just which strand, but that this one is the wrong one:
        # the antisense strand is the one intended to engage the promoter.
        reversed_target = NucleicSequence.parse("TTAAACGGATCCAAGGCGC")
        assessment = assess_asymmetry(reversed_target, LEGACY_DNA_KJ, end_length=4)
        assert assessment.predicted_guide is GuideStrand.SENSE
        assert "sense strand is predicted to be the guide" in assessment.explanation
        assert "undesirable" in assessment.explanation

    def test_an_ambiguous_explanation_says_no_strand_is_favoured(self) -> None:
        symmetrical = NucleicSequence.parse("GCGCACGGATCCAAGGCGC")
        assessment = assess_asymmetry(symmetrical, LEGACY_DNA_KJ, end_length=4)
        assert "neither" in assessment.explanation.lower()

    def test_the_explanation_mentions_the_defective_table_by_name(self) -> None:
        # A reader seeing a legacy verdict needs to know it came from the defective parameters
        # without cross-referencing the table list.
        assessment = assess_asymmetry(TARGET, LEGACY_DNA_KJ, end_length=4)
        assert LEGACY_DNA_KJ.name in assessment.explanation


class TestComparison:
    def test_reports_one_assessment_per_table(self) -> None:
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert [a.table_id for a in comparison.assessments] == [
            "rna-xia-1998",
            "legacy-dna-kj",
        ]

    def test_agreeing_tables_report_agreement(self) -> None:
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert comparison.predictions == {GuideStrand.ANTISENSE}
        assert comparison.tables_agree is True

    def test_the_summary_names_the_agreed_strand(self) -> None:
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert "antisense" in comparison.summary
        assert "agree" in comparison.summary

    def test_requires_at_least_one_table(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            compare_asymmetry(TARGET, (), end_length=4)

    def test_rejects_duplicate_tables(self) -> None:
        with pytest.raises(ValueError, match="more than once"):
            compare_asymmetry(TARGET, (RNA_XIA_1998, RNA_XIA_1998), end_length=4)

    def test_a_single_table_trivially_agrees(self) -> None:
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998,), end_length=4)
        assert comparison.tables_agree is True

    def test_a_single_table_summary_names_that_table_rather_than_claiming_agreement(
        self,
    ) -> None:
        # "All 1 parameter sets agree" would be an absurd thing to show a user.
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998,), end_length=4)
        assert comparison.summary.startswith(RNA_XIA_1998.name)
        assert "agree" not in comparison.summary


class TestDisagreementIsSurfaced:
    """The case ADR 0009 exists for: defect B2 moving a candidate across the threshold."""

    def _disagreeing_target(self) -> NucleicSequence | None:
        """Find a 19-mer where the RNA and legacy tables predict different guide strands.

        Searched rather than asserted, because whether such a sequence exists at all is a
        factual claim about the two parameter sets, not something to assume.
        """
        from itertools import product

        for head in product("ACGT", repeat=4):
            for tail in product("ACGT", repeat=4):
                candidate = NucleicSequence.parse(
                    "".join(head) + "ACGGATCCAAG" + "".join(tail)
                )
                rna = assess_asymmetry(candidate, RNA_XIA_1998, end_length=4)
                legacy = assess_asymmetry(candidate, LEGACY_DNA_KJ, end_length=4)
                if rna.predicted_guide is not legacy.predicted_guide:
                    return candidate
        return None

    def test_such_a_candidate_exists(self) -> None:
        assert self._disagreeing_target() is not None, (
            "if no 19-mer makes the two tables disagree, ADR 0009's premise is wrong and the "
            "comparison machinery is unnecessary"
        )

    def test_disagreement_is_reported_rather_than_resolved(self) -> None:
        target = self._disagreeing_target()
        assert target is not None
        comparison = compare_asymmetry(target, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert comparison.tables_agree is False
        assert len(comparison.predictions) > 1

    def test_the_summary_says_the_tables_disagree_and_names_both(self) -> None:
        target = self._disagreeing_target()
        assert target is not None
        comparison = compare_asymmetry(target, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert "disagree" in comparison.summary
        assert RNA_XIA_1998.name in comparison.summary
        assert LEGACY_DNA_KJ.name in comparison.summary

    def test_no_single_prediction_is_offered_when_tables_disagree(self) -> None:
        # There is deliberately no `comparison.predicted_guide`: picking one would be the
        # silent resolution ADR 0009 forbids.
        comparison = compare_asymmetry(TARGET, (RNA_XIA_1998, LEGACY_DNA_KJ), end_length=4)
        assert not hasattr(comparison, "predicted_guide")


class TestThresholdUnitGuards:
    """The asymmetry evaluator's unit handling, tested at its boundaries."""

    def test_an_unknown_table_id_is_rejected_by_name(self) -> None:
        from sarna.design.evaluators.thermodynamic_rules import _selected_table
        from sarna.rules.catalog.mandatory import THERMODYNAMIC_ASYMMETRY
        from sarna.rules.resolved import ResolvedRule

        resolved = ResolvedRule(THERMODYNAMIC_ASYMMETRY)
        # Bypasses the parameter's own `choices`, which is the only way this can happen: it
        # means the parameter and the shipped tables have drifted apart.
        object.__setattr__(
            resolved, "values", {**resolved.values, "nearest-neighbour-table": "rna-turner-2004"}
        )
        with pytest.raises(KeyError, match="rna-turner-2004"):
            _selected_table(resolved)

    def test_an_unconvertible_unit_is_refused_rather_than_compared(self) -> None:
        from sarna.design.evaluators.thermodynamic_rules import _threshold_in

        assert _threshold_in("kcal/mol", 2.0) == pytest.approx(2.0)
        assert _threshold_in("kJ/mol", 2.0) == pytest.approx(8.368)
        with pytest.raises(ValueError, match="joules"):
            _threshold_in("joules", 2.0)


class TestTheEndsAreExportedAsNumbersNotOnlyAsProse:
    """The evidence for a guide call must survive serialisation (ADR 0009).

    Until this was added, ``AsymmetryAssessment`` exported only the signed difference and an
    English sentence. The two end energies and the stacks they were summed from existed in
    Python and reached no consumer that was not a human reading a paragraph, so an interface
    could state the verdict but could not draw the measurement behind it.
    """

    def test_an_end_exports_its_energy_stacks_and_unit(self) -> None:
        # Legacy constants above, summed by hand: GC + CG + GC for the 5' end of GCGC...
        end = end_stability(TARGET, LEGACY_DNA_KJ, DuplexEnd.FIVE_PRIME)
        exported = end.to_dict()
        assert exported["end"] == "5'"
        assert exported["energy"] == pytest.approx(GCE + CGE + GCE)
        assert exported["unit"] == "kJ/mol"
        assert exported["stacks"] == [
            ["GC", pytest.approx(GCE)], ["CG", pytest.approx(CGE)], ["GC", pytest.approx(GCE)]
        ]

    def test_the_stacks_sum_to_the_energy_so_the_two_cannot_drift(self) -> None:
        for end in (DuplexEnd.FIVE_PRIME, DuplexEnd.THREE_PRIME):
            exported = end_stability(TARGET, RNA_XIA_1998, end).to_dict()
            assert sum(energy for _step, energy in exported["stacks"]) == pytest.approx(
                exported["energy"]
            )

    def test_an_assessment_carries_both_ends_so_a_reader_can_check_the_subtraction(self) -> None:
        assessment = assess_asymmetry(TARGET, RNA_XIA_1998)
        exported = assessment.to_dict()
        # The sign convention is E(3') - E(5'); exporting the difference alone makes that
        # claim uncheckable, which is precisely what ADR 0009 refuses.
        assert exported["three_prime"]["energy"] - exported["five_prime"]["energy"] == (
            pytest.approx(exported["difference"])
        )
        assert exported["predicted_guide"] == str(assessment.predicted_guide)
        assert exported["table"] == assessment.table_id
