"""Tests for the empirical benchmark harness.

The AUC is the number that would be quoted, so it is checked against hand-computed values and
against planted perfect and inverted rankings. Tie handling gets its own tests: preference points
are small integers, so ties are the common case rather than an edge case, and an AUC that ignored
them would read as signal.
"""

from collections.abc import Sequence
from typing import ClassVar

import pytest

from sarna.design.config import DesignConfig
from sarna.design.presets import EXPLORATORY, LEGACY_SCANNER, STRICT, WANG_2015
from sarna.validation.benchmark import (
    UNINFORMATIVE_AUC,
    BenchmarkResult,
    _auc,
    benchmark,
    score_sequence,
)
from sarna.validation.dataset import ValidationSequence, load_chapter4_dataset

DATASET = load_chapter4_dataset()


class TestAuc:
    def test_a_perfect_ranking_scores_one(self) -> None:
        assert _auc([3, 4, 5], [0, 1, 2]) == pytest.approx(1.0)

    def test_a_perfectly_inverted_ranking_scores_zero(self) -> None:
        assert _auc([0, 1, 2], [3, 4, 5]) == pytest.approx(0.0)

    def test_identical_values_score_one_half(self) -> None:
        # Every comparison is a tie, so the ranking is exactly uninformative.
        assert _auc([2, 2], [2, 2]) == pytest.approx(UNINFORMATIVE_AUC)

    def test_ties_count_as_half_a_win(self) -> None:
        # One win, one tie, out of two comparisons: 1.5 / 2.
        assert _auc([2], [1, 2]) == pytest.approx(0.75)

    def test_a_hand_computed_mixed_case(self) -> None:
        # group 3 beats 1 and 2 (2 wins); group 1 ties 1 (0.5) and loses to 2 (0).
        # 2.5 wins out of 4 comparisons.
        assert _auc([3, 1], [1, 2]) == pytest.approx(0.625)


class TestScoringOneSequence:
    def test_a_sequence_is_scored_as_a_candidate(self) -> None:
        scored = score_sequence(DATASET.by_name("P21-322"), WANG_2015.resolve())
        assert scored.candidate.id == "P21-322"
        assert scored.candidate.target.bases == "CCAACTCATTCTCCAAGTA"
        assert scored.is_activating

    def test_rules_needing_promoter_context_are_unevaluable_not_failed(self) -> None:
        # The dataset records no coordinates. Reporting these as failures would silently
        # penalise every sequence in it for a property of the dataset.
        scored = score_sequence(DATASET.by_name("P21-322"), WANG_2015.resolve())
        for rule_id in ("tss-window", "tss-preferred-window", "flank-20-au"):
            assert scored.candidate.result(rule_id).is_unevaluated

    def test_sequence_intrinsic_rules_are_evaluated(self) -> None:
        scored = score_sequence(DATASET.by_name("P21-322"), WANG_2015.resolve())
        for rule_id in ("gc-content", "cpg-count", "pos19-a", "no-homopolymer-run"):
            assert not scored.candidate.result(rule_id).is_unevaluated

    def test_an_unevaluable_rule_contributes_no_points(self) -> None:
        scored = score_sequence(DATASET.by_name("P21-322"), WANG_2015.resolve())
        assert scored.candidate.result("flank-20-au").points == 0


class TestBenchmarkMechanics:
    def test_both_classes_are_required(self) -> None:
        from sarna.validation.dataset import ValidationDataset

        only_positive = ValidationDataset(
            id="x",
            name="x",
            sequences=DATASET.activating,
            provenance=DATASET.provenance,
            caveats=(),
        )
        with pytest.raises(ValueError, match="both classes"):
            benchmark(only_positive, WANG_2015.resolve(), permutations=10)

    def test_the_seed_makes_the_p_value_reproducible(self) -> None:
        first = benchmark(DATASET, WANG_2015.resolve(), permutations=200, seed=5)
        second = benchmark(DATASET, WANG_2015.resolve(), permutations=200, seed=5)
        assert first.auc_p_value == second.auc_p_value
        assert first.seed == 5

    def test_every_sequence_is_scored(self) -> None:
        result = benchmark(DATASET, WANG_2015.resolve(), permutations=10)
        assert len(result.scored) == 125
        assert len(result.activating) == 14
        assert len(result.not_activating) == 111

    def test_the_result_records_which_preset_produced_it(self) -> None:
        result = benchmark(DATASET, WANG_2015.resolve(), permutations=10)
        assert result.preset_id == "wang-2015"

    def test_a_bare_config_has_no_preset_id(self) -> None:
        result = benchmark(DATASET, DesignConfig.of(["gc-content"]), permutations=10)
        assert result.preset_id is None

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps(benchmark(DATASET, WANG_2015.resolve(), permutations=10).to_dict())

    def test_the_export_carries_the_activators_the_gate_lost(self) -> None:
        exported = benchmark(DATASET, LEGACY_SCANNER.resolve(), permutations=10).to_dict()
        assert exported["counts"]["activators_rejected_by_a_gate"] == 5


class TestTheMeasuredResults:
    """The findings in docs/scientific/empirical-validation.md, pinned.

    These are the project's only empirical claims about its own rules. A change that altered any
    of them without a documented reason would be a regression in the science, not in the code.
    """

    WANG: ClassVar = benchmark(DATASET, WANG_2015.resolve(), permutations=400, seed=3)
    LEGACY: ClassVar = benchmark(DATASET, LEGACY_SCANNER.resolve(), permutations=400, seed=3)

    def test_wangs_gates_exclude_only_one_known_activator(self) -> None:
        assert len(self.WANG.rejected_activators) == 1

    def test_the_one_activator_wang_excludes_is_cebpa_on_position_19(self) -> None:
        # MTL-CEBPA, the first saRNA to reach human clinical trials, fails Table 1.
        rejected = self.WANG.rejected_activators[0]
        assert rejected.sequence.name == "CEBPA"
        assert [r.rule_id for r in rejected.candidate.excluded_by] == ["pos19-a"]

    def test_the_legacy_gates_exclude_five_known_activators(self) -> None:
        names = {s.sequence.name for s in self.LEGACY.rejected_activators}
        assert names == {"P53", "SCN1A-10", "CEBPA", "DEPDC5-4", "KLF4"}

    def test_the_legacy_cpg_gate_is_what_loses_most_of_them(self) -> None:
        by_cpg = {
            s.sequence.name
            for s in self.LEGACY.rejected_activators
            if any(r.rule_id == "cpg-count" for r in s.candidate.excluded_by)
        }
        assert by_cpg == {"P53", "SCN1A-10", "CEBPA", "KLF4"}

    def test_activators_carry_more_cpg_than_non_activators(self) -> None:
        # The direction that removes any empirical case for tightening the gate to zero.
        def mean_cpg(sequences: Sequence[ValidationSequence]) -> float:
            counts = [s.target_dna.cpg_count for s in sequences]
            return sum(counts) / len(counts)

        assert mean_cpg(DATASET.activating) > mean_cpg(DATASET.not_activating)

    def test_the_unsourced_position_7_rule_is_anti_correlated_with_activity(self) -> None:
        # Finding B18. The legacy script weighted this at 1000.
        pos7 = next(r for r in self.LEGACY.rule_discrimination if r.rule_id == "pos7-u")
        assert pos7.difference < -15

    def test_no_preset_ranking_beats_chance(self) -> None:
        for preset in (WANG_2015, STRICT, EXPLORATORY, LEGACY_SCANNER):
            result = benchmark(DATASET, preset.resolve(), permutations=400, seed=3)
            assert not result.discriminates, f"{preset.id} unexpectedly discriminates"

    def test_the_summary_says_the_ranking_is_uninformative(self) -> None:
        assert "uninformative" in self.WANG.summary

    def test_the_chapter_four_rules_separate_best_because_they_were_fitted_here(self) -> None:
        # Circular by construction: pos10-g was derived from this dataset. Asserted so the
        # circularity is visible in the test suite rather than only in the write-up.
        result = benchmark(DATASET, EXPLORATORY.resolve(), permutations=10)
        best = result.rule_discrimination[0]
        assert best.rule_id == "pos10-g"
        assert best.difference == pytest.approx(25.0, abs=1.0)

    def test_strict_lost_no_more_activators_than_the_default(self) -> None:
        strict = benchmark(DATASET, STRICT.resolve(), permutations=10)
        assert len(strict.rejected_activators) == len(self.WANG.rejected_activators)


class TestRuleDiscriminationReporting:
    RESULT: ClassVar = benchmark(DATASET, WANG_2015.resolve(), permutations=10)

    def test_every_selected_rule_is_reported(self) -> None:
        reported = {r.rule_id for r in self.RESULT.rule_discrimination}
        assert reported == {resolved.rule.id for resolved in WANG_2015.resolve()}

    def test_rules_are_ordered_by_absolute_separation(self) -> None:
        magnitudes = [abs(r.difference) for r in self.RESULT.rule_discrimination]
        assert magnitudes == sorted(magnitudes, reverse=True)

    def test_a_rule_without_promoter_context_is_marked_unevaluable(self) -> None:
        flank = next(r for r in self.RESULT.rule_discrimination if r.rule_id == "flank-20-au")
        assert flank.was_never_evaluated
        assert flank.unevaluated_count == 125
        assert "not evaluable without promoter context" in str(flank)

    def test_an_evaluated_rule_renders_both_pass_rates(self) -> None:
        gc = next(r for r in self.RESULT.rule_discrimination if r.rule_id == "gc-content")
        assert "of activators pass vs" in str(gc)


class TestSummaryRendering:
    """Both branches of the verdict wording, since only one occurs on the real dataset."""

    @staticmethod
    def _result(auc: float, p_value: float) -> BenchmarkResult:
        return BenchmarkResult(
            preset_id="x",
            scored=(),
            auc=auc,
            auc_p_value=p_value,
            permutations=100,
            seed=1,
            rule_discrimination=(),
            max_preference_points=7,
        )

    def test_a_discriminating_ranking_says_it_carries_signal(self) -> None:
        result = self._result(0.82, 0.004)
        assert result.discriminates
        assert "carries signal" in result.summary
        assert "82%" in result.summary

    def test_an_uninformative_ranking_names_the_coin_flip(self) -> None:
        result = self._result(0.46, 0.71)
        assert not result.discriminates
        assert "uninformative" in result.summary
        assert "coin flip" in result.summary

    def test_an_auc_above_one_half_is_described_as_above(self) -> None:
        assert "above the 50%" in self._result(0.53, 0.40).summary

    def test_an_auc_below_one_half_is_described_as_below(self) -> None:
        assert "below the 50%" in self._result(0.41, 0.90).summary
