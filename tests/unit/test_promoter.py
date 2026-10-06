"""Tests for the promoter region and the window-to-coordinates mapping.

This is the defect-B1 surface. The legacy EpiMarker asked Ensembl for a minus-strand region,
received the reverse complement, then added the match index to the region start as though it
were a forward offset — so every minus-strand position and every TSS-relative distance it
reported was wrong. The mapping tested here is the one place that arithmetic happens, and the
minus-strand cases below are written out with explicit coordinates rather than derived, so that
a regression cannot hide behind the same formula it broke.
"""

import pytest

from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.sequence.nucleic import NucleicSequence

# A 20 nt promoter occupying 1000-1020 on contig 1, used for both strands so that the two
# mappings can be compared base for base.
PLUS_INTERVAL = GenomicInterval(Assembly.GRCH38, "1", 1000, 1020, Strand.PLUS)
MINUS_INTERVAL = GenomicInterval(Assembly.GRCH38, "1", 1000, 1020, Strand.MINUS)
TWENTY_MER = "ACGTACGTACGTACGTACGT"


def promoter(interval: GenomicInterval | None = PLUS_INTERVAL, **kwargs: object) -> Promoter:
    return Promoter(
        sequence=NucleicSequence.parse(TWENTY_MER),
        interval=interval,
        **kwargs,  # type: ignore[arg-type]
    )


def mapped_interval(region: Promoter, offset: int, length: int) -> GenomicInterval:
    """Return a window's genomic interval, asserting the promoter is anchored."""
    interval = region.genomic_interval(offset=offset, length=length)
    assert interval is not None, "the promoter under test should be anchored"
    return interval


class TestConstruction:
    def test_a_sequence_only_promoter_needs_no_coordinates(self) -> None:
        # Workflow B: a pasted promoter with no genomic context still designs candidates.
        region = Promoter(sequence=NucleicSequence.parse(TWENTY_MER))
        assert region.interval is None
        assert region.tss_start is None
        assert region.is_anchored is False

    def test_an_anchored_promoter_reports_itself_as_anchored(self) -> None:
        region = promoter(tss_start=TssPosition(-1000))
        assert region.is_anchored is True

    def test_rejects_an_rna_sequence(self) -> None:
        # A promoter is DNA. Accepting RNA here would make CpG content unanswerable later.
        with pytest.raises(ValueError, match="DNA"):
            Promoter(sequence=NucleicSequence.parse(TWENTY_MER).as_rna())

    def test_rejects_an_interval_whose_length_differs_from_the_sequence(self) -> None:
        with pytest.raises(ValueError, match="20"):
            Promoter(
                sequence=NucleicSequence.parse(TWENTY_MER),
                interval=GenomicInterval(Assembly.GRCH38, "1", 1000, 1010, Strand.PLUS),
            )

    def test_length(self) -> None:
        assert promoter().length == 20


class TestReferenceFactory:
    def test_a_plus_strand_promoter_keeps_the_reference_sequence(self) -> None:
        region = Promoter.from_reference_sequence(
            NucleicSequence.parse(TWENTY_MER), PLUS_INTERVAL
        )
        assert region.sequence.bases == TWENTY_MER

    def test_a_minus_strand_promoter_reverse_complements_the_reference(self) -> None:
        # The engine always holds the promoter 5'->3' in the transcript's own orientation, so a
        # plus-strand reference sequence must be flipped for a minus-strand gene. Doing this at
        # the boundary is what stops the rest of the engine needing to know about strand.
        region = Promoter.from_reference_sequence(
            NucleicSequence.parse(TWENTY_MER), MINUS_INTERVAL
        )
        assert region.sequence == NucleicSequence.parse(TWENTY_MER).reverse_complement()


class TestTssSpan:
    def test_the_first_window_starts_at_the_promoter_start(self) -> None:
        region = promoter(tss_start=TssPosition(-1000))
        assert region.tss_span(offset=0, length=19) == (
            TssPosition(-1000),
            TssPosition(-982),
        )

    def test_a_later_window_is_shifted_by_its_offset(self) -> None:
        region = promoter(tss_start=TssPosition(-1000))
        assert region.tss_span(offset=1, length=19) == (
            TssPosition(-999),
            TssPosition(-981),
        )

    def test_a_span_crossing_the_tss_skips_zero(self) -> None:
        region = promoter(tss_start=TssPosition(-10))
        assert region.tss_span(offset=0, length=19) == (TssPosition(-10), TssPosition(9))

    def test_an_unanchored_promoter_has_no_tss_span(self) -> None:
        # Not an error: it is the signal that the TSS rules must report not-applicable.
        assert Promoter(sequence=NucleicSequence.parse(TWENTY_MER)).tss_span(0, 19) is None

    def test_the_span_is_the_same_on_both_strands(self) -> None:
        # TSS-relative positions are measured in the transcript's own frame, so they do not
        # depend on strand at all. Any difference here would mean the frame leaked.
        plus = promoter(PLUS_INTERVAL, tss_start=TssPosition(-1000))
        minus = promoter(MINUS_INTERVAL, tss_start=TssPosition(-1000))
        assert plus.tss_span(1, 19) == minus.tss_span(1, 19)


class TestGenomicMappingOnThePlusStrand:
    def test_the_first_window_starts_at_the_interval_start(self) -> None:
        window = mapped_interval(promoter(PLUS_INTERVAL), 0, 19)
        assert window == GenomicInterval(Assembly.GRCH38, "1", 1000, 1019, Strand.PLUS)

    def test_a_later_window_moves_forward(self) -> None:
        window = mapped_interval(promoter(PLUS_INTERVAL), 1, 19)
        assert window == GenomicInterval(Assembly.GRCH38, "1", 1001, 1020, Strand.PLUS)

    def test_a_single_base_window(self) -> None:
        window = mapped_interval(promoter(PLUS_INTERVAL), 5, 1)
        assert (window.start, window.end) == (1005, 1006)


class TestGenomicMappingOnTheMinusStrand:
    def test_the_first_window_ends_at_the_interval_end(self) -> None:
        # The promoter's first base is its 5' end in the transcript's orientation, which on the
        # minus strand is the *highest* genomic coordinate. A 19 nt window at offset 0 therefore
        # occupies 1001-1020, not 1000-1019. Getting this backwards is defect B1.
        window = mapped_interval(promoter(MINUS_INTERVAL), 0, 19)
        assert window == GenomicInterval(Assembly.GRCH38, "1", 1001, 1020, Strand.MINUS)

    def test_a_later_window_moves_backward_in_genomic_coordinates(self) -> None:
        window = mapped_interval(promoter(MINUS_INTERVAL), 1, 19)
        assert window == GenomicInterval(Assembly.GRCH38, "1", 1000, 1019, Strand.MINUS)

    def test_a_single_base_window(self) -> None:
        window = mapped_interval(promoter(MINUS_INTERVAL), 5, 1)
        assert (window.start, window.end) == (1014, 1015)

    def test_the_fifth_prime_end_of_the_window_is_its_highest_coordinate(self) -> None:
        window = mapped_interval(promoter(MINUS_INTERVAL), 0, 19)
        assert window.five_prime == 1019

    def test_the_mapped_strand_is_the_promoters_strand(self) -> None:
        window = mapped_interval(promoter(MINUS_INTERVAL), 0, 19)
        assert window.strand is Strand.MINUS


class TestGenomicMappingInvariants:
    @pytest.mark.parametrize("strand_interval", [PLUS_INTERVAL, MINUS_INTERVAL])
    @pytest.mark.parametrize("offset", [0, 1])
    def test_the_mapped_interval_always_has_the_window_length(
        self, strand_interval: GenomicInterval, offset: int
    ) -> None:
        assert mapped_interval(promoter(strand_interval), offset, 19).length == 19

    @pytest.mark.parametrize("strand_interval", [PLUS_INTERVAL, MINUS_INTERVAL])
    def test_the_mapped_interval_stays_inside_the_promoter(
        self, strand_interval: GenomicInterval
    ) -> None:
        region = promoter(strand_interval)
        for offset, _ in region.sequence.windows(19):
            inside = mapped_interval(region, offset, 19)
            assert inside.start >= strand_interval.start
            assert inside.end <= strand_interval.end

    def test_the_two_strands_map_offset_zero_to_opposite_ends(self) -> None:
        # The clearest statement of the asymmetry: the same offset is at opposite ends of the
        # same genomic span depending on strand.
        plus = mapped_interval(promoter(PLUS_INTERVAL), 0, 19)
        minus = mapped_interval(promoter(MINUS_INTERVAL), 0, 19)
        assert plus.start == 1000
        assert minus.end == 1020

    def test_an_unanchored_promoter_has_no_genomic_mapping(self) -> None:
        unanchored = Promoter(sequence=NucleicSequence.parse(TWENTY_MER))
        assert unanchored.genomic_interval(offset=0, length=19) is None

    @pytest.mark.parametrize(("offset", "length"), [(-1, 19), (0, 21), (2, 19), (0, 0)])
    def test_rejects_a_window_outside_the_promoter(self, offset: int, length: int) -> None:
        with pytest.raises(ValueError):
            promoter().genomic_interval(offset=offset, length=length)

    @pytest.mark.parametrize(("offset", "length"), [(-1, 19), (0, 21), (0, 0)])
    def test_the_tss_span_rejects_the_same_windows(self, offset: int, length: int) -> None:
        region = promoter(tss_start=TssPosition(-1000))
        with pytest.raises(ValueError):
            region.tss_span(offset=offset, length=length)


class TestSequenceAgreesWithCoordinates:
    """The mapping is correct only if the sequence found at the reported coordinates matches.

    These tests check that invariant directly against a stand-in reference, rather than
    re-deriving the same formula the implementation uses. This is the check the legacy pipeline
    never made, and making it would have caught defect B1 immediately.
    """

    REFERENCE = NucleicSequence.parse("TTGACCGGATCCTAGGCAAT")

    def _reference_slice(self, mapped: GenomicInterval) -> NucleicSequence:
        """Return the plus-strand reference sequence at ``mapped``, as a genome server would."""
        offset_into_reference = mapped.start - 1000
        return self.REFERENCE.subsequence(
            offset_into_reference + 1, offset_into_reference + mapped.length
        )

    @pytest.mark.parametrize("offset", [0, 1])
    def test_a_plus_strand_window_matches_the_reference_at_its_coordinates(
        self, offset: int
    ) -> None:
        region = Promoter.from_reference_sequence(self.REFERENCE, PLUS_INTERVAL)
        window = region.sequence.subsequence(offset + 1, offset + 19)
        assert window == self._reference_slice(mapped_interval(region, offset, 19))

    @pytest.mark.parametrize("offset", [0, 1])
    def test_a_minus_strand_window_matches_the_reverse_complement_of_the_reference(
        self, offset: int
    ) -> None:
        region = Promoter.from_reference_sequence(self.REFERENCE, MINUS_INTERVAL)
        window = region.sequence.subsequence(offset + 1, offset + 19)
        reference = self._reference_slice(mapped_interval(region, offset, 19))
        assert window == reference.reverse_complement()

    def test_the_minus_strand_window_is_not_the_plain_reference_slice(self) -> None:
        # Guards the guard: if the reference happened to be its own reverse complement, the test
        # above would pass under defect B1's arithmetic too.
        region = Promoter.from_reference_sequence(self.REFERENCE, MINUS_INTERVAL)
        reference = self._reference_slice(mapped_interval(region, 0, 19))
        assert reference != reference.reverse_complement()


class TestTssCoordinate:
    def test_locates_the_tss_on_the_plus_strand(self) -> None:
        # The promoter's first base is -20, so the TSS (+1) is 20 bp further on: 1000 + 20.
        region = promoter(PLUS_INTERVAL, tss_start=TssPosition(-20))
        assert region.tss_coordinate == 1020

    def test_locates_the_tss_on_the_minus_strand(self) -> None:
        # On the minus strand "further on" means a lower coordinate: 1019 - 20.
        region = promoter(MINUS_INTERVAL, tss_start=TssPosition(-20))
        assert region.tss_coordinate == 999

    def test_an_unanchored_promoter_cannot_locate_a_tss(self) -> None:
        assert Promoter(sequence=NucleicSequence.parse(TWENTY_MER)).tss_coordinate is None

    def test_a_promoter_with_no_interval_cannot_locate_a_tss(self) -> None:
        region = Promoter(
            sequence=NucleicSequence.parse(TWENTY_MER), tss_start=TssPosition(-20)
        )
        assert region.tss_coordinate is None


@pytest.mark.strand
class TestGenomicToPromoterMapping:
    """``offset_of`` and ``span_of`` are the inverse of ``genomic_interval``.

    Needed to place anything arriving in genomic coordinates onto a promoter held in transcript
    orientation. Verified as a round trip against the forward mapping rather than by restating
    the arithmetic, which is what would make a sign error invisible.
    """

    def plus(self) -> Promoter:
        return Promoter(
            sequence=NucleicSequence.parse("ACGT" * 25),
            tss_start=TssPosition(-100),
            interval=GenomicInterval(Assembly.GRCH38, "1", 1000, 1100, Strand.PLUS),
        )

    def minus(self) -> Promoter:
        return Promoter(
            sequence=NucleicSequence.parse("ACGT" * 25),
            tss_start=TssPosition(-100),
            interval=GenomicInterval(Assembly.GRCH38, "1", 1000, 1100, Strand.MINUS),
        )

    def test_on_the_plus_strand_offset_zero_is_the_lowest_coordinate(self) -> None:
        assert self.plus().offset_of(1000) == 0
        assert self.plus().offset_of(1099) == 99

    def test_on_the_minus_strand_offset_zero_is_the_highest_coordinate(self) -> None:
        # The whole point. Getting this backwards is defect B1 in its other direction.
        assert self.minus().offset_of(1099) == 0
        assert self.minus().offset_of(1000) == 99

    def test_the_mapping_round_trips_against_genomic_interval_on_both_strands(self) -> None:
        for promoter in (self.plus(), self.minus()):
            for offset in (0, 1, 19, 50, 80):
                window = promoter.genomic_interval(offset, 19)
                assert window is not None
                recovered = promoter.span_of(window)
                assert recovered == (offset, offset + 19)

    def test_a_coordinate_outside_the_promoter_is_none_not_clamped(self) -> None:
        # A feature starting before the promoter must be drawn truncated, not at zero.
        assert self.plus().offset_of(999) is None
        assert self.plus().offset_of(1100) is None

    def test_an_unanchored_promoter_cannot_map_anything(self) -> None:
        bare = Promoter(sequence=NucleicSequence.parse("ACGT" * 25))
        assert bare.offset_of(1000) is None
        interval = GenomicInterval(Assembly.GRCH38, "1", 1000, 1100)
        assert bare.span_of(interval) is None

    def test_a_span_overlapping_one_edge_is_clipped(self) -> None:
        overlapping = GenomicInterval(Assembly.GRCH38, "1", 990, 1010, Strand.PLUS)
        assert self.plus().span_of(overlapping) == (0, 10)

    def test_a_span_clipped_on_the_minus_strand_lands_at_the_other_end(self) -> None:
        # Low coordinates are the promoter's 3' end on the minus strand, so a span hanging off
        # the low side clips to the *end* of the promoter, not the start.
        overlapping = GenomicInterval(Assembly.GRCH38, "1", 990, 1010, Strand.MINUS)
        assert self.minus().span_of(overlapping) == (90, 100)

    def test_a_span_that_does_not_overlap_is_none(self) -> None:
        elsewhere = GenomicInterval(Assembly.GRCH38, "1", 2000, 2100)
        assert self.plus().span_of(elsewhere) is None

    def test_a_span_on_another_contig_is_an_error(self) -> None:
        other = GenomicInterval(Assembly.GRCH38, "2", 1000, 1100)
        with pytest.raises(ValueError, match="contig"):
            self.plus().span_of(other)

    def test_a_span_on_another_assembly_is_an_error(self) -> None:
        # Silently producing offsets across assemblies is defect B15.
        other = GenomicInterval(Assembly.GRCH37, "1", 1000, 1100)
        with pytest.raises(ValueError, match="liftover adapter"):
            self.plus().span_of(other)

    def test_a_tss_just_past_the_promoters_end_maps_to_none(self) -> None:
        # These fixtures start at -100 and are 100 bp long, so +1 is exactly one base past the
        # last offset. The forward and inverse mappings must agree that it is outside.
        for promoter in (self.plus(), self.minus()):
            coordinate = promoter.tss_coordinate
            assert coordinate is not None
            assert promoter.offset_of(coordinate) is None

    def test_a_tss_inside_the_promoter_maps_to_the_offset_tss_start_implies(self) -> None:
        # 120 bp starting at -100, so +1 sits at offset 100 and the window runs to +20.
        for strand in (Strand.PLUS, Strand.MINUS):
            promoter = Promoter(
                sequence=NucleicSequence.parse("ACGT" * 30),
                tss_start=TssPosition(-100),
                interval=GenomicInterval(Assembly.GRCH38, "1", 1000, 1120, strand),
            )
            coordinate = promoter.tss_coordinate
            assert coordinate is not None
            assert promoter.offset_of(coordinate) == 100
