"""Tests for the per-candidate evaluation context (ADR 0007).

The context is the single input every rule evaluator reads. Two behaviours matter most: the
1-based position accessors must agree with the way Wang 2015 Table 1 numbers nucleotides, and
every genomic field must be absent rather than guessed when the promoter is unanchored, so that
coordinate-dependent rules can report not-applicable.
"""

import pytest

from sarna.design.context import CandidateContext
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.sequence.nucleic import NucleicSequence

# 30 bp, so a 19 nt target at offset 0 has a full 4 nt flank and room to slide.
PROMOTER_BASES = "GCTTACGGATCCAAGTTGACAGGCCTATAA"
INTERVAL = GenomicInterval(Assembly.GRCH38, "1", 1000, 1030, Strand.PLUS)
MINUS_INTERVAL = GenomicInterval(Assembly.GRCH38, "1", 1000, 1030, Strand.MINUS)


def anchored(interval: GenomicInterval = INTERVAL) -> Promoter:
    return Promoter(
        sequence=NucleicSequence.parse(PROMOTER_BASES),
        tss_start=TssPosition(-1000),
        interval=interval,
    )


def unanchored() -> Promoter:
    return Promoter(sequence=NucleicSequence.parse(PROMOTER_BASES))


class TestTarget:
    def test_the_target_is_the_window_at_the_offset(self) -> None:
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.target.bases == PROMOTER_BASES[:19]

    def test_a_later_offset_slides_the_window(self) -> None:
        context = CandidateContext.at(unanchored(), offset=3, length=19)
        assert context.target.bases == PROMOTER_BASES[3:22]

    def test_the_length_is_the_target_length(self) -> None:
        assert CandidateContext.at(unanchored(), offset=0, length=19).length == 19

    def test_the_target_is_dna(self) -> None:
        # The promoter target site is DNA, so CpG content remains answerable.
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.target.cpg_count >= 0

    def test_rejects_a_window_outside_the_promoter(self) -> None:
        with pytest.raises(ValueError):
            CandidateContext.at(unanchored(), offset=20, length=19)


class TestPositions:
    def test_base_at_uses_table_one_numbering(self) -> None:
        # Position 1 is the target's first base and position 19 its last, exactly as Table 1
        # numbers them.
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.base_at(1) == PROMOTER_BASES[0]
        assert context.base_at(19) == PROMOTER_BASES[18]

    def test_base_at_is_relative_to_the_target_not_the_promoter(self) -> None:
        context = CandidateContext.at(unanchored(), offset=3, length=19)
        assert context.base_at(1) == PROMOTER_BASES[3]

    @pytest.mark.parametrize("position", [0, 20, -1])
    def test_base_at_rejects_a_position_outside_the_target(self, position: int) -> None:
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        with pytest.raises(IndexError):
            context.base_at(position)


class TestFlank:
    def test_the_flank_follows_the_target(self) -> None:
        # Positions 20-23 in Table 1's numbering: the four bases immediately 3' of a 19 nt
        # target, which lie in the promoter rather than in the duplex.
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.flank is not None
        assert context.flank.bases == PROMOTER_BASES[19:23]

    @pytest.mark.parametrize(
        ("position", "promoter_index"),
        [(20, 19), (21, 20), (22, 21), (23, 22)],
    )
    def test_flank_base_at_uses_table_one_numbering(
        self, position: int, promoter_index: int
    ) -> None:
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.flank_base_at(position) == PROMOTER_BASES[promoter_index]

    def test_a_truncated_flank_reports_the_bases_it_has(self) -> None:
        # Offset 9 leaves only 2 bases after a 19 nt target in a 30 bp promoter.
        context = CandidateContext.at(unanchored(), offset=9, length=19)
        assert context.flank is not None
        assert len(context.flank) == 2
        assert context.flank_base_at(20) == PROMOTER_BASES[28]
        assert context.flank_base_at(22) is None

    def test_a_target_at_the_promoter_end_has_no_flank(self) -> None:
        # Not an error, and not silently treated as a failed preference: the flank preferences
        # must report not-applicable, which is why this is None rather than an empty sequence.
        context = CandidateContext.at(unanchored(), offset=11, length=19)
        assert context.flank is None
        assert context.flank_base_at(20) is None

    def test_flank_base_at_rejects_a_position_inside_the_target(self) -> None:
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        with pytest.raises(ValueError, match="19"):
            context.flank_base_at(19)


class TestTssSpan:
    def test_an_anchored_promoter_gives_the_span(self) -> None:
        context = CandidateContext.at(anchored(), offset=0, length=19)
        assert context.tss_span == (TssPosition(-1000), TssPosition(-982))

    def test_an_unanchored_promoter_has_no_span(self) -> None:
        assert CandidateContext.at(unanchored(), offset=0, length=19).tss_span is None

    def test_has_tss_context_reports_whether_the_window_rules_can_run(self) -> None:
        assert CandidateContext.at(anchored(), offset=0, length=19).has_tss_context is True
        assert CandidateContext.at(unanchored(), offset=0, length=19).has_tss_context is False


class TestTssWindowContainment:
    """ADR 0008: every base of the target must lie inside the window."""

    def test_a_target_wholly_inside_the_window_passes(self) -> None:
        promoter = Promoter(
            sequence=NucleicSequence.parse(PROMOTER_BASES), tss_start=TssPosition(-500)
        )
        context = CandidateContext.at(promoter, offset=0, length=19)
        assert context.tss_window_contains(-1000, -100) is True

    def test_a_target_extending_past_the_downstream_bound_fails(self) -> None:
        # 5' end at -105, so the target runs to -87 and overlaps the transcript. A 5'-anchored
        # test would wrongly pass this.
        promoter = Promoter(
            sequence=NucleicSequence.parse(PROMOTER_BASES), tss_start=TssPosition(-105)
        )
        context = CandidateContext.at(promoter, offset=0, length=19)
        assert context.tss_window_contains(-1000, -100) is False

    def test_a_target_ending_exactly_on_the_bound_passes(self) -> None:
        promoter = Promoter(
            sequence=NucleicSequence.parse(PROMOTER_BASES), tss_start=TssPosition(-118)
        )
        context = CandidateContext.at(promoter, offset=0, length=19)
        assert context.tss_span == (TssPosition(-118), TssPosition(-100))
        assert context.tss_window_contains(-1000, -100) is True

    def test_a_target_starting_before_the_upstream_bound_fails(self) -> None:
        promoter = Promoter(
            sequence=NucleicSequence.parse(PROMOTER_BASES), tss_start=TssPosition(-1001)
        )
        context = CandidateContext.at(promoter, offset=0, length=19)
        assert context.tss_window_contains(-1000, -100) is False

    def test_an_unanchored_promoter_cannot_answer(self) -> None:
        # None, not False: the rule must say it could not be evaluated, not that it failed.
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        assert context.tss_window_contains(-1000, -100) is None

    def test_rejects_a_window_whose_bounds_are_the_wrong_way_round(self) -> None:
        context = CandidateContext.at(anchored(), offset=0, length=19)
        with pytest.raises(ValueError, match="start"):
            context.tss_window_contains(-100, -1000)


class TestGenomicInterval:
    def test_an_anchored_promoter_gives_the_interval(self) -> None:
        context = CandidateContext.at(anchored(), offset=0, length=19)
        assert context.interval == GenomicInterval(
            Assembly.GRCH38, "1", 1000, 1019, Strand.PLUS
        )

    def test_the_interval_follows_the_strand(self) -> None:
        context = CandidateContext.at(anchored(MINUS_INTERVAL), offset=0, length=19)
        assert context.interval is not None
        assert context.interval.end == 1030

    def test_an_unanchored_promoter_has_no_interval(self) -> None:
        assert CandidateContext.at(unanchored(), offset=0, length=19).interval is None

    def test_has_coordinates_reports_whether_annotation_can_run(self) -> None:
        assert CandidateContext.at(anchored(), offset=0, length=19).has_coordinates is True
        assert CandidateContext.at(unanchored(), offset=0, length=19).has_coordinates is False


class TestValueSemantics:
    def test_is_immutable(self) -> None:
        context = CandidateContext.at(unanchored(), offset=0, length=19)
        with pytest.raises((AttributeError, TypeError)):
            context.offset = 5  # type: ignore[misc]

    def test_equal_contexts_compare_equal(self) -> None:
        assert CandidateContext.at(unanchored(), offset=0, length=19) == (
            CandidateContext.at(unanchored(), offset=0, length=19)
        )
