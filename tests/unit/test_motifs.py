"""Tests for positional motif analysis and its significance test.

The statistics here decide whether a reported enrichment becomes a design rule, so they are
verified against hand-computed values and against a planted signal, not only against the real
dataset. A permutation test that always returned 0.9 would look correct on data with no signal.
"""

from typing import ClassVar

import pytest

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence
from sarna.validation.dataset import Activity, load_chapter4_dataset
from sarna.validation.motifs import (
    PermutationResult,
    frequency_differences,
    motifs_at,
    permutation_test,
    positional_frequencies,
)


def rna(bases: str) -> NucleicSequence:
    return NucleicSequence.parse(bases, Nucleic.RNA)


class TestMotifsAt:
    def test_a_single_base_is_read_1_based(self) -> None:
        assert motifs_at(rna("ACGU"), 1, 1) == "A"
        assert motifs_at(rna("ACGU"), 4, 1) == "U"

    def test_a_dinucleotide_starts_at_the_given_position(self) -> None:
        assert motifs_at(rna("ACGU"), 2, 2) == "CG"

    def test_a_motif_running_past_the_end_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="does not fit"):
            motifs_at(rna("ACGU"), 4, 2)

    def test_position_zero_is_an_error(self) -> None:
        # Positions inside a sequence are 1-based throughout this project.
        with pytest.raises(ValueError, match="does not fit"):
            motifs_at(rna("ACGU"), 0, 1)


class TestPositionalFrequencies:
    def test_frequencies_are_percentages_of_the_group(self) -> None:
        # Two sequences, one with A at position 1, one with C: 50% each.
        frequencies = positional_frequencies([rna("AC"), rna("CC")])
        assert frequencies[(1, "A")] == pytest.approx(50.0)
        assert frequencies[(1, "C")] == pytest.approx(50.0)
        assert frequencies[(2, "C")] == pytest.approx(100.0)

    def test_an_unobserved_motif_is_absent_rather_than_zero(self) -> None:
        frequencies = positional_frequencies([rna("AA")])
        assert (1, "C") not in frequencies

    def test_dinucleotide_positions_stop_one_short_of_the_end(self) -> None:
        frequencies = positional_frequencies([rna("ACGU")], motif_length=2)
        assert set(frequencies) == {(1, "AC"), (2, "CG"), (3, "GU")}

    def test_ragged_sequences_are_refused(self) -> None:
        # Comparing position 19 of a 19-mer against position 19 of a 21-mer compares different
        # parts of different molecules. The caller must decide how to align them.
        with pytest.raises(ValueError, match="equal-length"):
            positional_frequencies([rna("AC"), rna("ACG")])

    def test_an_empty_group_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no sequences"):
            positional_frequencies([])


class TestFrequencyDifferences:
    def test_differences_are_group_minus_reference_in_points(self) -> None:
        differences = frequency_differences([rna("AA")], [rna("CC")])
        by_cell = {(d.position, d.motif): d.difference for d in differences}
        assert by_cell[(1, "A")] == pytest.approx(100.0)
        assert by_cell[(1, "C")] == pytest.approx(-100.0)

    def test_a_motif_absent_from_one_group_counts_as_zero_there(self) -> None:
        # So a depletion is reported as readily as an enrichment.
        differences = frequency_differences([rna("AA")], [rna("CC")])
        assert any(d.motif == "C" and d.difference < 0 for d in differences)

    def test_results_are_ordered_by_absolute_difference(self) -> None:
        differences = frequency_differences(
            [rna("AAA"), rna("AAC")], [rna("CCC"), rna("ACC")]
        )
        magnitudes = [abs(d.difference) for d in differences]
        assert magnitudes == sorted(magnitudes, reverse=True)

    def test_a_difference_renders_readably(self) -> None:
        difference = frequency_differences([rna("A")], [rna("C")])[0]
        assert str(difference) == "pos 1 A: +100.0 pp"


class TestPermutationTest:
    def test_a_planted_perfect_signal_is_detected(self) -> None:
        # Every member of one group has G at position 1, none of the other does. If the test
        # cannot find this, it cannot find anything.
        group = [rna("GAAAA") for _ in range(8)]
        reference = [rna("CAAAA") for _ in range(8)]
        result = permutation_test(group, reference, permutations=500, seed=1)
        assert result.observed_maximum == pytest.approx(100.0)
        assert result.is_significant
        assert result.exceeds_chance

    def test_identical_groups_show_no_signal(self) -> None:
        sequences = [rna("ACGUA"), rna("CGUAC"), rna("GUACG"), rna("UACGU")]
        result = permutation_test(sequences, list(sequences), permutations=500, seed=1)
        assert not result.is_significant

    def test_the_p_value_is_plus_one_corrected(self) -> None:
        # With no shuffle reaching the observed maximum the honest answer is 1/(n+1), never 0.
        group = [rna("GAAAA") for _ in range(8)]
        reference = [rna("CAAAA") for _ in range(8)]
        result = permutation_test(group, reference, permutations=99, seed=1)
        assert result.family_wise_p >= 1 / 100
        assert result.family_wise_p > 0

    def test_the_seed_is_reported_so_the_number_can_be_reproduced(self) -> None:
        first = permutation_test(
            [rna("GAAAA")] * 4, [rna("CAAAA")] * 4, permutations=200, seed=7
        )
        second = permutation_test(
            [rna("GAAAA")] * 4, [rna("CAAAA")] * 4, permutations=200, seed=7
        )
        assert first.seed == 7
        assert first.family_wise_p == second.family_wise_p

    def test_a_different_seed_is_allowed_to_differ(self) -> None:
        sequences = [rna("ACGUA"), rna("CGUAC"), rna("GUACG"), rna("UACGU")]
        one = permutation_test(sequences[:2], sequences[2:], permutations=200, seed=1)
        two = permutation_test(sequences[:2], sequences[2:], permutations=200, seed=2)
        assert one.permutations == two.permutations == 200

    def test_the_cell_count_records_how_large_the_grid_was(self) -> None:
        result = permutation_test(
            [rna("ACGUA")] * 3, [rna("UGCAU")] * 3, permutations=100, seed=1
        )
        assert result.cells > 1

    def test_an_empty_group_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no sequences"):
            permutation_test([], [rna("ACGUA")], permutations=10)


class TestTheChapterFourResult:
    """The finding that governs how the thesis's motif analysis may be used (ADR 0011).

    Pinned as a test because it is the reason the positional rules are off by default. If a
    future change made these enrichments look significant, that would be a bug in the analysis
    rather than a discovery, and it should break here.
    """

    DATASET: ClassVar = load_chapter4_dataset()
    ACTIVATING: ClassVar = [s.target for s in DATASET.activating]
    INACTIVE: ClassVar = [s.target for s in DATASET.not_activating]

    def test_the_chapters_headline_enrichment_reproduces(self) -> None:
        differences = {
            (d.position, d.motif): d.difference
            for d in frequency_differences(self.ACTIVATING, self.INACTIVE)
        }
        assert differences[(10, "G")] == pytest.approx(24.8, abs=0.1)
        assert differences[(15, "U")] == pytest.approx(23.0, abs=0.1)
        assert differences[(9, "C")] == pytest.approx(-21.7, abs=0.1)

    def test_the_chapters_dinucleotide_claims_reproduce(self) -> None:
        differences = {
            (d.position, d.motif): d.difference
            for d in frequency_differences(self.ACTIVATING, self.INACTIVE, 2)
        }
        assert differences[(15, "UG")] == pytest.approx(24.1, abs=0.1)
        assert differences[(12, "GG")] == pytest.approx(20.5, abs=0.1)
        assert differences[(17, "AA")] == pytest.approx(18.6, abs=0.1)

    def test_the_chapters_ac_at_position_9_is_really_ag(self) -> None:
        # The chapter lists "AC (pos 9)" at +19.6 pp. The motif at that value is AG; AC is
        # +0.8 pp. Recorded because the 1-mer result (G enriched, C depleted at 9) predicts AG.
        differences = {
            (d.position, d.motif): d.difference
            for d in frequency_differences(self.ACTIVATING, self.INACTIVE, 2)
        }
        assert differences[(9, "AG")] == pytest.approx(19.6, abs=0.1)
        assert differences[(9, "AC")] == pytest.approx(0.8, abs=0.1)

    def test_the_positional_signal_does_not_survive_the_grid_it_was_found_in(self) -> None:
        # Fewer permutations than the published run so the suite stays fast; the verdict is
        # stable well below 20 000 because the observed maximum is below the null median.
        result = permutation_test(self.ACTIVATING, self.INACTIVE, 1, permutations=400, seed=11)
        assert not result.is_significant
        assert not result.exceeds_chance
        assert "smaller than chance" in result.summary

    def test_the_grid_searched_is_76_single_nucleotide_cells(self) -> None:
        # 19 positions x 4 bases, all of them observed in 125 sequences.
        differences = frequency_differences(self.ACTIVATING, self.INACTIVE)
        assert len(differences) == 76


class TestResultRendering:
    """The summary wording is what a reader sees, so each of its three branches is pinned."""

    @staticmethod
    def _result(observed: float, median: float, p: float) -> PermutationResult:
        return PermutationResult(
            observed_maximum=observed,
            null_median=median,
            null_95th_percentile=median + 10,
            family_wise_p=p,
            permutations=1000,
            seed=1,
            cells=76,
        )

    def test_a_significant_result_says_larger_than_chance(self) -> None:
        summary = self._result(90.0, 30.0, 0.001).summary
        assert "larger than chance" in summary
        assert "family-wise p = 0.001" in summary

    def test_a_non_significant_result_above_the_null_median_says_within_chance(self) -> None:
        # The middle case: more positional structure than a typical shuffle, but not enough to
        # survive the grid. It must not read as either a finding or an absence of one.
        result = self._result(35.0, 30.0, 0.30)
        assert result.exceeds_chance
        assert not result.is_significant
        assert "within what chance" in result.summary

    def test_a_result_below_the_null_median_says_smaller_than_chance(self) -> None:
        result = self._result(24.8, 31.1, 0.94)
        assert not result.exceeds_chance
        assert "smaller than chance" in result.summary
        assert "no positional signal" in result.summary


class TestActivityLabel:
    def test_an_activity_renders_as_its_label(self) -> None:
        assert str(Activity.ACTIVATING) == "activating"
        assert str(Activity.NOT_ACTIVATING) == "not-activating"
