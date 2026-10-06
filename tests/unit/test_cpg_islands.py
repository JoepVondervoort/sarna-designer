"""Tests for CpG-island detection.

The expected island boundaries below are derived by hand from the criteria, not read off the
implementation. The worked derivation is in the test that uses it, because an island boundary
that nobody can re-derive is an island boundary nobody can check.
"""

from itertools import pairwise

import pytest

from sarna.sequence.cpg import (
    GARDINER_GARDEN_FROMMER_1987,
    TAKAI_JONES_2002,
    CpGIslandCriteria,
    find_cpg_islands,
    overlaps_any_island,
)
from sarna.sequence.nucleic import NucleicSequence

# "CG" repeated is the strongest possible island signal: 100 % GC, maximal CpG density.
ALL_ISLAND = NucleicSequence.parse("CG" * 300)
# "AT" repeated has no G, no C and no CpG at all.
NO_ISLAND = NucleicSequence.parse("AT" * 300)


class TestCriteria:
    def test_the_takai_jones_thresholds(self) -> None:
        assert TAKAI_JONES_2002.minimum_length == 500
        assert TAKAI_JONES_2002.minimum_gc_fraction == pytest.approx(0.55)
        assert TAKAI_JONES_2002.minimum_observed_expected == pytest.approx(0.65)

    def test_the_gardiner_garden_frommer_thresholds(self) -> None:
        assert GARDINER_GARDEN_FROMMER_1987.minimum_length == 200
        assert GARDINER_GARDEN_FROMMER_1987.minimum_gc_fraction == pytest.approx(0.50)
        assert GARDINER_GARDEN_FROMMER_1987.minimum_observed_expected == pytest.approx(0.60)

    def test_rejects_a_non_positive_length(self) -> None:
        with pytest.raises(ValueError, match="length"):
            CpGIslandCriteria(
                minimum_length=0,
                minimum_gc_fraction=0.55,
                minimum_observed_expected=0.65,
            )

    @pytest.mark.parametrize("fraction", [-0.1, 1.1])
    def test_rejects_a_gc_fraction_outside_zero_to_one(self, fraction: float) -> None:
        with pytest.raises(ValueError, match="fraction"):
            CpGIslandCriteria(
                minimum_length=500,
                minimum_gc_fraction=fraction,
                minimum_observed_expected=0.65,
            )

    def test_rejects_a_negative_cpg_count(self) -> None:
        with pytest.raises(ValueError, match="count"):
            CpGIslandCriteria(
                minimum_length=500,
                minimum_gc_fraction=0.55,
                minimum_observed_expected=0.65,
                minimum_cpg_count=-1,
            )

    def test_rejects_a_negative_observed_expected_ratio(self) -> None:
        with pytest.raises(ValueError, match="ratio"):
            CpGIslandCriteria(
                minimum_length=500,
                minimum_gc_fraction=0.55,
                minimum_observed_expected=-1.0,
            )


class TestDetection:
    def test_a_sequence_of_pure_cpg_is_entirely_one_island(self) -> None:
        islands = find_cpg_islands(ALL_ISLAND, GARDINER_GARDEN_FROMMER_1987)
        assert len(islands) == 1
        assert (islands[0].start, islands[0].end) == (0, 600)

    def test_rejects_an_rna_sequence(self) -> None:
        # CpG content is a property of the DNA target site, so asking of RNA is a mistake
        # rather than a question with an answer.
        with pytest.raises(ValueError, match="DNA"):
            find_cpg_islands(ALL_ISLAND.as_rna(), TAKAI_JONES_2002)

    def test_a_sequence_with_no_g_or_c_has_no_islands(self) -> None:
        assert find_cpg_islands(NO_ISLAND, GARDINER_GARDEN_FROMMER_1987) == ()

    def test_a_sequence_shorter_than_the_window_has_no_islands(self) -> None:
        # Not an error and not a warning printed to stdout, as the legacy detector does: a short
        # promoter simply cannot support a positive call, and the caller decides what to say.
        short = NucleicSequence.parse("CG" * 50)
        assert find_cpg_islands(short, TAKAI_JONES_2002) == ()

    def test_a_sequence_exactly_the_window_length_can_be_an_island(self) -> None:
        # The legacy detector tests `len(sequence) > 200` and so discards a sequence of exactly
        # the window length. The criteria say "at least".
        exact = NucleicSequence.parse("CG" * 100)
        islands = find_cpg_islands(exact, GARDINER_GARDEN_FROMMER_1987)
        assert len(islands) == 1
        assert (islands[0].start, islands[0].end) == (0, 200)

    def test_the_stricter_criteria_find_no_island_in_a_short_one(self) -> None:
        # 300 bp of perfect CpG is an island under the 1987 criteria and not under the 2002
        # ones, which is the whole substance of the choice between them.
        three_hundred = NucleicSequence.parse("CG" * 150)
        assert find_cpg_islands(three_hundred, GARDINER_GARDEN_FROMMER_1987) != ()
        assert find_cpg_islands(three_hundred, TAKAI_JONES_2002) == ()


class TestBoundaries:
    """A worked case with boundaries derived by hand from the Takai & Jones criteria.

    The sequence is 1000 bp of ``AT`` repeat followed by 1000 bp of ``CG`` repeat, so the
    GC-rich block occupies ``[1000, 2000)``.

    A 500 bp window starting at ``s`` overlaps the block by ``s + 500 - 1000`` bases while
    ``s < 1000``. The GC criterion needs at least ``0.55 * 500 = 275`` GC bases, so
    ``s - 500 >= 275``, giving ``s >= 775``. At ``s = 775`` the window holds 275 block bases,
    about 138 CpGs against an expected ``138 * 137 / 500 = 37.8``, a ratio of 3.65 — far above
    0.65 — so the GC content is the binding constraint and the first qualifying window starts at
    775. The last window that fits starts at 1500. Merging gives ``[775, 2000)``.
    """

    SEQUENCE = NucleicSequence.parse("AT" * 500 + "CG" * 500)

    def test_the_island_spans_the_hand_derived_boundaries(self) -> None:
        islands = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)
        assert len(islands) == 1
        assert (islands[0].start, islands[0].end) == (775, 2000)

    def test_the_island_reports_its_length(self) -> None:
        island = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)[0]
        assert island.length == 1225

    def test_the_island_statistics_cover_the_merged_span_not_the_first_window(self) -> None:
        # Reporting the first qualifying window's numbers would understate a long island's GC
        # content, since the window that first qualifies is the least GC-rich part of it.
        island = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)[0]
        assert island.gc_fraction == pytest.approx(1000 / 1225, abs=1e-6)
        assert island.observed_cpg == 500

    def test_the_island_reports_its_observed_expected_ratio(self) -> None:
        island = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)[0]
        expected = (500 * 500) / 1225
        assert island.expected_cpg == pytest.approx(expected)
        assert island.observed_expected_ratio == pytest.approx(500 / expected)


class TestTwoIslands:
    def test_two_separated_blocks_are_reported_separately(self) -> None:
        sequence = NucleicSequence.parse("CG" * 400 + "AT" * 600 + "CG" * 400)
        islands = find_cpg_islands(sequence, TAKAI_JONES_2002)
        assert len(islands) == 2
        assert islands[0].start == 0
        assert islands[1].end == len(sequence)

    def test_the_islands_are_returned_in_sequence_order(self) -> None:
        sequence = NucleicSequence.parse("CG" * 400 + "AT" * 600 + "CG" * 400)
        islands = find_cpg_islands(sequence, TAKAI_JONES_2002)
        assert islands[0].end <= islands[1].start

    def test_the_islands_do_not_overlap(self) -> None:
        sequence = NucleicSequence.parse("CG" * 400 + "AT" * 600 + "CG" * 400)
        islands = find_cpg_islands(sequence, TAKAI_JONES_2002)
        assert all(first.end <= second.start for first, second in pairwise(islands))


class TestCpGCountGuard:
    def test_a_window_with_too_few_cpgs_is_not_an_island(self) -> None:
        # A GC-rich window with its Cs and Gs segregated has a tiny expected CpG count, so the
        # observed-to-expected ratio can pass on almost no CpGs at all. The count guard is what
        # stops that arithmetic artefact being reported as an island.
        segregated = NucleicSequence.parse("C" * 300 + "G" * 300)
        criteria = CpGIslandCriteria(
            minimum_length=200,
            minimum_gc_fraction=0.50,
            minimum_observed_expected=0.60,
            minimum_cpg_count=7,
        )
        assert find_cpg_islands(segregated, criteria) == ()

    def test_disabling_the_guard_admits_the_artefact(self) -> None:
        # Shown explicitly so the guard's effect is visible rather than assumed.
        segregated = NucleicSequence.parse("C" * 300 + "G" * 300)
        criteria = CpGIslandCriteria(
            minimum_length=200,
            minimum_gc_fraction=0.50,
            minimum_observed_expected=0.60,
            minimum_cpg_count=0,
        )
        assert find_cpg_islands(segregated, criteria) != ()


class TestOverlapQuery:
    SEQUENCE = NucleicSequence.parse("AT" * 500 + "CG" * 500)

    def test_a_target_wholly_inside_an_island_overlaps(self) -> None:
        islands = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)
        assert overlaps_any_island(islands, start=1000, end=1019) is True

    def test_a_target_wholly_outside_any_island_does_not_overlap(self) -> None:
        islands = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)
        assert overlaps_any_island(islands, start=0, end=19) is False

    def test_a_target_partly_overlapping_an_island_counts_as_overlapping(self) -> None:
        # The island starts at 775, so a target at 770-789 has 14 of its 19 bases outside it.
        # It still counts: this rule excludes candidates, so the cautious reading is the right
        # one. A methylatable CpG anywhere in the target site can block engagement.
        islands = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)
        assert overlaps_any_island(islands, start=770, end=789) is True

    def test_a_target_ending_exactly_where_an_island_starts_does_not_overlap(self) -> None:
        # Half-open spans: end 775 is the base before the island's first base.
        islands = find_cpg_islands(self.SEQUENCE, TAKAI_JONES_2002)
        assert overlaps_any_island(islands, start=756, end=775) is False

    def test_no_islands_means_no_overlap(self) -> None:
        assert overlaps_any_island((), start=0, end=19) is False

    def test_rejects_an_empty_or_inverted_span(self) -> None:
        with pytest.raises(ValueError):
            overlaps_any_island((), start=19, end=0)


class TestIncrementalCountingMatchesBruteForce:
    """The sliding counts are maintained incrementally; this checks them against rescanning.

    The incremental CpG update is the subtlest arithmetic in the module — a CpG leaves when the
    departing base was the C of a pair, and arrives when the entering base completes one — and an
    error in it would shift island boundaries silently. Comparing against a naive recount over
    randomised sequences is the cheapest way to keep it honest.
    """

    CRITERIA = CpGIslandCriteria(
        minimum_length=50,
        minimum_gc_fraction=0.50,
        minimum_observed_expected=0.60,
        minimum_cpg_count=3,
    )

    @staticmethod
    def _brute_force(bases: str, criteria: CpGIslandCriteria) -> list[tuple[int, int]]:
        """Find islands by rescanning each window from scratch."""
        window = criteria.minimum_length
        spans: list[list[int]] = []
        for start in range(len(bases) - window + 1):
            piece = bases[start : start + window]
            c_count, g_count = piece.count("C"), piece.count("G")
            expected = (c_count * g_count) / window
            qualifies = (
                piece.count("CG") >= criteria.minimum_cpg_count
                and (c_count + g_count) / window >= criteria.minimum_gc_fraction
                and expected > 0
                and piece.count("CG") / expected >= criteria.minimum_observed_expected
            )
            if not qualifies:
                continue
            if spans and start <= spans[-1][1]:
                spans[-1][1] = start + window
            else:
                spans.append([start, start + window])
        return [(start, end) for start, end in spans]

    def test_agrees_with_brute_force_over_randomised_sequences(self) -> None:
        import random

        # Seeded, so a failure is reproducible. The alphabet is weighted towards G and C so that
        # islands actually occur rather than the test passing on empty results.
        generator = random.Random(20261006)
        for _ in range(300):
            bases = "".join(
                generator.choices("ACGTCGCG", k=generator.randint(50, 300))
            )
            sequence = NucleicSequence.parse(bases)
            found = [
                (island.start, island.end)
                for island in find_cpg_islands(sequence, self.CRITERIA)
            ]
            assert found == self._brute_force(bases, self.CRITERIA), bases

    def test_the_brute_force_reference_actually_finds_islands(self) -> None:
        # Guards the guard: if the reference found nothing the comparison above would be vacuous.
        assert self._brute_force("CG" * 100, self.CRITERIA) == [(0, 200)]


class TestPerformance:
    def test_detection_is_linear_enough_for_a_long_promoter(self) -> None:
        # Defect B7: the legacy detector recomputes the expected CpG count with str.count on
        # every step, making it O(n * window). A 10 kb promoter is an ordinary request, so this
        # asserts the counts are maintained incrementally rather than rescanned.
        import time

        long_promoter = NucleicSequence.parse(("ACGTTGCA" * 1250)[:10_000])
        started = time.perf_counter()
        find_cpg_islands(long_promoter, TAKAI_JONES_2002)
        assert time.perf_counter() - started < 1.0


class TestUniformRandomSequenceIsCalledAsAnIsland:
    """Documents why a synthetic promoter fixture must be CpG-depleted.

    A CpG island is defined by CpG being *more frequent than composition predicts*. Uniform
    random sequence has observed/expected around 1.0 by construction, so any correct detector
    calls it an island. Real DNA outside islands is depleted to around 0.2 because methylated
    CpG deaminates to TpG, and that depletion is what makes an island stand out.

    This is kept as a test because the first end-to-end design run used a uniform random
    promoter and returned zero candidates, which looked like a bug in the design rules.
    See docs/scientific/synthetic-promoter-fixtures.md.
    """

    @staticmethod
    def _random(bases: str, weights: list[int], seed: int, length: int = 1200) -> str:
        import random

        return "".join(random.Random(seed).choices(bases, weights=weights, k=length))

    def test_uniform_random_sequence_is_called_as_an_island(self) -> None:
        uniform = NucleicSequence.parse(self._random("ACGT", [22, 28, 28, 22], seed=7))
        assert find_cpg_islands(uniform, TAKAI_JONES_2002) != ()

    def test_cpg_depleted_sequence_of_the_same_kind_is_not(self) -> None:
        import random

        generator = random.Random(7)
        out: list[str] = []
        for _ in range(1200):
            base = generator.choices("ACGT", weights=[26, 24, 24, 26])[0]
            # Suppress most CpG, as deamination has done to real genomic DNA.
            if out and out[-1] == "C" and base == "G" and generator.random() > 0.2:
                base = generator.choice("AT")
            out.append(base)
        depleted = NucleicSequence.parse("".join(out))
        # The observed/expected ratio is the quantity the criteria actually test, and the one
        # that distinguishes real genomic baseline (around 0.2) from random sequence (1.0).
        bases = depleted.bases
        expected = (bases.count("C") * bases.count("G")) / len(bases)
        assert depleted.cpg_count / expected < TAKAI_JONES_2002.minimum_observed_expected
        assert find_cpg_islands(depleted, TAKAI_JONES_2002) == ()

    def test_uniform_random_sequence_has_a_ratio_near_one(self) -> None:
        # The reason the test above works: expected is derived from the same composition that
        # generated the sequence, so the ratio is 1.0 by construction.
        uniform = NucleicSequence.parse(self._random("ACGT", [22, 28, 28, 22], seed=7))
        bases = uniform.bases
        expected = (bases.count("C") * bases.count("G")) / len(bases)
        assert uniform.cpg_count / expected == pytest.approx(1.0, abs=0.1)
