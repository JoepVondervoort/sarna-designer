"""The candidate's own locus is not an off-target finding.

Every candidate is a stretch of a real human promoter, so a genomic search finds that promoter.
Counted, it makes "perfect matches elsewhere" read at least 1 for every candidate in the run,
which is the same as reading nothing at all. ADR 0023 removes it from the findings, keeps it
on the record, and states the count.

Identification is by coordinate and never by sequence or gene name. These tests work the
conversion out by hand, because this is the boundary where BLAST's 1-based inclusive subject
coordinates meet this project's 0-based half-open intervals, and an off-by-one here would
quietly remove a real off-target hit next door to the candidate.
"""

from sarna.offtarget.models import OffTargetHit
from sarna.offtarget.ontarget import (
    is_on_target,
    loci_from_candidates,
    split_on_target,
)

# chr2:100,000-100,019 as the run records it: 0-based, half-open, 19 bases.
LOCUS = {"saRNA-1": ("2", 100_000, 100_019)}


def hit(
    candidate: str = "saRNA-1",
    accession: str = "NC_000002.12",
    start: int = 100_001,
    end: int = 100_019,
) -> OffTargetHit:
    """Build a genomic hit. Subject coordinates are BLAST's: 1-based and inclusive."""
    return OffTargetHit(
        candidate_id=candidate,
        accession=accession,
        title="Homo sapiens chromosome 2",
        identities=19,
        alignment_length=19,
        gaps=0,
        query_start=1,
        query_end=19,
        subject_start=start,
        subject_end=end,
        subject_strand="plus",
        evalue=1e-5,
        bit_score=38.2,
    )


class TestIdentifyingTheCandidatesOwnLocus:
    def test_the_hit_at_the_candidates_own_coordinates_is_recognised(self) -> None:
        # 1-based 100,001 to 100,019 inclusive is 0-based half-open 100,000 to 100,019: the
        # candidate's own interval exactly.
        assert is_on_target(hit(), LOCUS)

    def test_a_hit_on_another_chromosome_is_not(self) -> None:
        assert not is_on_target(hit(accession="NC_000003.12"), LOCUS)

    def test_a_hit_elsewhere_on_the_same_chromosome_is_not(self) -> None:
        assert not is_on_target(hit(start=200_001, end=200_019), LOCUS)

    def test_the_base_immediately_after_the_locus_is_not_on_target(self) -> None:
        # 1-based 100,020 is 0-based 100,019, which is the first base *outside* a half-open
        # interval ending at 100,019. Treating it as inside would remove a genuine neighbouring
        # hit, and nothing downstream would ever show that it had happened.
        assert not is_on_target(hit(start=100_020, end=100_038), LOCUS)

    def test_the_last_base_of_the_locus_is_on_target(self) -> None:
        # The other side of the same boundary: 1-based 100,019 is 0-based 100,018, inside.
        assert is_on_target(hit(start=100_019, end=100_037), LOCUS)

    def test_a_minus_strand_alignment_reports_its_coordinates_descending(self) -> None:
        # BLAST writes subject-from > subject-to on the minus strand. Read literally that is an
        # empty interval and the own-locus hit would survive on half of all candidates.
        assert is_on_target(hit(start=100_019, end=100_001), LOCUS)

    def test_a_candidate_with_no_recorded_locus_keeps_every_hit(self) -> None:
        # A run from pasted sequence has no coordinates. Nothing can be identified, so nothing
        # is removed, and the screening says so rather than guessing.
        assert not is_on_target(hit(), {})

    def test_an_unplaced_contig_is_not_matched_to_a_chromosome(self) -> None:
        assert not is_on_target(hit(accession="NT_187361.1"), LOCUS)


class TestSplitting:
    def test_the_own_locus_hit_is_moved_rather_than_deleted(self) -> None:
        elsewhere = hit(start=500_001, end=500_019)
        off, own = split_on_target([hit(), elsewhere], LOCUS)
        assert off == (elsewhere,)
        assert own == (hit(),)

    def test_order_is_preserved_within_each_list(self) -> None:
        first = hit(start=500_001, end=500_019)
        second = hit(start=700_001, end=700_019)
        off, _ = split_on_target([first, hit(), second], LOCUS)
        assert off == (first, second)

    def test_one_candidates_locus_does_not_remove_another_candidates_hit(self) -> None:
        # Attribution runs through the query id. A hit of saRNA-2 that happens to land on
        # saRNA-1's locus is a real finding: saRNA-2 was not designed there.
        other = hit(candidate="saRNA-2")
        off, own = split_on_target([other], LOCUS)
        assert off == (other,)
        assert own == ()


class TestReadingLociFromARun:
    def test_coordinates_are_taken_from_the_candidate_record(self) -> None:
        loci = loci_from_candidates([
            {"id": "saRNA-1",
             "genomic": {"contig": "2", "start": 100_000, "end": 100_019, "strand": "+"}},
        ])
        assert loci == {"saRNA-1": ("2", 100_000, 100_019)}

    def test_a_candidate_without_coordinates_is_absent_rather_than_invented(self) -> None:
        loci = loci_from_candidates([
            {"id": "saRNA-1", "genomic": None},
            {"id": "saRNA-2"},
        ])
        assert loci == {}

    def test_a_malformed_record_is_skipped_rather_than_crashing_the_screening(self) -> None:
        # Snapshots are read back from storage and from frozen exports, so this is a boundary.
        loci = loci_from_candidates([
            {"id": "saRNA-1", "genomic": {"contig": "2", "start": "100000", "end": 100_019}},
            {"id": "saRNA-2", "genomic": {"contig": "7", "start": 10, "end": 29}},
        ])
        assert loci == {"saRNA-2": ("7", 10, 29)}
