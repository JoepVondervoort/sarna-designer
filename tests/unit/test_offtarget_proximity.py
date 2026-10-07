"""Placing genomic hits against annotated start sites, with no invented threshold.

The expected distances here are computed by hand, not taken from the implementation. The sign
convention is the one the rest of the engine uses: negative is upstream, in the transcript's own
direction, which means the arithmetic flips on the minus strand.
"""

from typing import Any

from sarna.offtarget.models import OffTargetHit
from sarna.offtarget.proximity import (
    contig_for,
    describe,
    place_hit,
    signed_distance,
    tss_of,
)


def hit(accession: str = "NC_000002.12", start: int = 166_128_000) -> OffTargetHit:
    """Build a genomic hit at a chosen coordinate."""
    return OffTargetHit(
        candidate_id="saRNA-1", accession=accession, title="Homo sapiens chromosome 2",
        identities=19, alignment_length=19, gaps=0, query_start=1, query_end=19,
        subject_start=start, subject_end=start + 19, subject_strand="plus",
        evalue=0.004, bit_score=38.2,
    )


def lookup_returning(*transcripts: dict[str, Any]):
    """Return a lookup that always answers with ``transcripts``."""
    return lambda _interval: list(transcripts)


class TestAccessionMapping:
    def test_refseq_chromosome_accessions_map_to_ensembl_contigs(self) -> None:
        assert contig_for("NC_000001.11") == "1"
        assert contig_for("NC_000002.12") == "2"
        assert contig_for("NC_000022.11") == "22"
        assert contig_for("NC_000023.11") == "X"
        assert contig_for("NC_000024.10") == "Y"

    def test_an_unplaced_scaffold_is_not_guessed_at(self) -> None:
        assert contig_for("NT_187361.1") is None
        assert contig_for("NW_009646201.1") is None

    def test_an_unmapped_hit_is_reported_without_a_placement(self) -> None:
        placed = place_hit(hit(accession="NT_187361.1"), lookup_returning(), window=None)
        assert placed.nearest_gene is None
        assert placed.distance_to_tss is None
        assert placed.inside_run_window is False
        assert "could not be established" in describe(placed)


class TestTheStartSiteOfATranscript:
    def test_a_plus_strand_transcript_starts_at_its_low_coordinate(self) -> None:
        assert tss_of({"start": 1000, "end": 5000, "strand": 1}) == 1000

    def test_a_minus_strand_transcript_starts_at_its_high_coordinate(self) -> None:
        # The invariant the rest of the engine keeps: on the minus strand the start site is the
        # highest coordinate, not the lowest.
        assert tss_of({"start": 1000, "end": 5000, "strand": -1}) == 5000


class TestSignedDistance:
    def test_upstream_is_negative_on_the_plus_strand(self) -> None:
        # Hand-computed: a hit 300 bp before a TSS at 10 000.
        assert signed_distance(position=9_700, tss=10_000, strand=1) == -300
        assert signed_distance(position=10_400, tss=10_000, strand=1) == 400

    def test_the_sign_flips_on_the_minus_strand(self) -> None:
        # On the minus strand a *lower* coordinate is downstream, so the same arithmetic
        # inverted. A hit at 9 700 against a minus-strand TSS at 10 000 is 300 bp downstream.
        assert signed_distance(position=9_700, tss=10_000, strand=-1) == 300
        assert signed_distance(position=10_400, tss=10_000, strand=-1) == -400


class TestPlacement:
    def test_the_nearest_transcript_wins_regardless_of_direction(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning(
                {"id": "ENST-far", "external_name": "FAR", "start": 90_000, "end": 95_000,
                 "strand": 1},
                {"id": "ENST-near", "external_name": "NEAR", "start": 100_600, "end": 110_000,
                 "strand": 1},
            ),
            window=None,
        )
        assert placed.nearest_transcript == "ENST-near"
        assert placed.nearest_gene == "NEAR"
        assert placed.distance_to_tss == -600

    def test_no_threshold_decides_whether_this_is_a_promoter_hit(self) -> None:
        # A hit 40 bp from a start site is reported as 40 bp from a start site. The tool does
        # not pronounce it "in a promoter": that definition has not been defended here.
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "G", "start": 100_040,
                              "end": 120_000, "strand": 1}),
            window=None,
        )
        assert placed.distance_to_tss == -40
        exported = placed.to_dict()
        assert "is_promoter" not in exported
        assert "in_promoter" not in exported

    def test_the_only_flag_refers_to_the_runs_own_window(self) -> None:
        # Wang's default window, -1000 to -100. A hit 600 bp upstream of another gene's start
        # site would have been a legal target for that gene under this very configuration.
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "OTHER", "start": 100_600,
                              "end": 120_000, "strand": 1}),
            window=(-1000, -100),
        )
        assert placed.distance_to_tss == -600
        assert placed.inside_run_window is True
        assert "same window this run designed against" in describe(placed)

    def test_a_hit_outside_the_runs_window_is_not_flagged(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "OTHER", "start": 105_000,
                              "end": 120_000, "strand": 1}),
            window=(-1000, -100),
        )
        assert placed.distance_to_tss == -5000
        assert placed.inside_run_window is False

    def test_a_run_with_no_window_flags_nothing(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "G", "start": 100_600,
                              "end": 120_000, "strand": 1}),
            window=None,
        )
        assert placed.inside_run_window is False

    def test_nothing_within_the_search_radius_leaves_the_distance_unknown(self) -> None:
        placed = place_hit(hit(start=100_000), lookup_returning(), window=(-1000, -100))
        assert placed.distance_to_tss is None
        assert placed.nearest_gene is None

    def test_a_transcript_without_coordinates_is_skipped_not_guessed(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "broken"},
                             {"id": "E1", "external_name": "G", "start": 100_200,
                              "end": 120_000, "strand": 1}),
            window=None,
        )
        assert placed.nearest_transcript == "E1"


class TestTheDescription:
    def test_it_states_the_distance_and_the_gene(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "SCN1A", "start": 100_600,
                              "end": 120_000, "strand": 1}),
            window=None,
        )
        text = describe(placed)
        assert "600 bp upstream of" in text
        assert "SCN1A" in text

    def test_it_never_claims_the_promoter_is_active_or_that_binding_occurs(self) -> None:
        placed = place_hit(
            hit(start=100_000),
            lookup_returning({"id": "E1", "external_name": "OTHER", "start": 100_600,
                              "end": 120_000, "strand": 1}),
            window=(-1000, -100),
        )
        assert "does not establish that the promoter is active" in describe(placed)
