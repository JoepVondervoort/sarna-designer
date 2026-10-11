"""Running a screening end to end, with the clock and the network both injected.

No test here sleeps or reaches the network. What they pin is the behaviour a researcher depends
on while waiting: that the whole run goes in one job, that each stage is reported, that waiting
is bounded, and that a candidate with nothing found is still present in the table rather than
missing from it.
"""

from typing import Any

import httpx
import pytest

from sarna.genome.coordinates import GenomicInterval
from sarna.offtarget.models import GENOMIC_SCOPE, TRANSCRIPT_SCOPE
from sarna.offtarget.ncbi import MINIMUM_POLL_SECONDS, BlastClient, BlastError
from sarna.offtarget.screening import MAX_CANDIDATES, run_screening, specificity_table
from tests.unit.test_offtarget_ncbi import (
    PUT_REPLY,
    READY_NO_HITS,
    READY_WITH_HITS,
    REPORT,
    WAITING,
)

CANDIDATES: list[dict[str, Any]] = [
    {"id": "saRNA-1549", "rank": 1, "target_site_dna": "GTGAAGCCATAGATCTGTA"},
    {"id": "saRNA-0002", "rank": 2, "target_site_dna": "GGGTAATTTTGACAGGTCA"},
    {"id": "saRNA-9999", "rank": None, "target_site_dna": "TTTTTTTTTTTTTTTTTTT"},
]


def scripted(replies: list[bytes]) -> tuple[BlastClient, list[httpx.Request]]:
    """Return a client answering with ``replies`` in order, holding the last one."""
    seen: list[httpx.Request] = []

    def handler(request: httpx.Request) -> httpx.Response:
        seen.append(request)
        return httpx.Response(200, content=replies[min(len(seen) - 1, len(replies) - 1)])

    return BlastClient(client=httpx.Client(transport=httpx.MockTransport(handler))), seen


def never_sleep(_seconds: float) -> None:
    """Stand in for time.sleep so a poll costs nothing."""


class TestTheHappyPath:
    def test_a_screening_returns_hits_with_its_own_provenance(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, REPORT])
        record = run_screening(
            "run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep
        )
        assert record.run_id == "run-1"
        assert record.request_id == "ABC123XYZ014"
        assert record.parameters.database == "refseq_rna"
        assert record.parameters.word_size == 7
        # The moment matters: NCBI's databases change, so a result without a date is unusable.
        assert record.searched_at is not None
        assert len(record.hits) == 2

    def test_only_candidates_that_passed_the_gates_are_screened(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep
        )
        assert record.candidates_searched == ("saRNA-1549", "saRNA-0002")

    def test_each_stage_is_reported_so_the_interface_need_not_spin(self) -> None:
        stages: list[str] = []
        client, _ = scripted([PUT_REPLY, WAITING, READY_WITH_HITS, REPORT])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      progress=stages.append, sleep=never_sleep)
        assert any("Submitting" in s for s in stages)
        assert any("accepted" in s for s in stages)
        assert any("still searching" in s for s in stages)
        assert any("alignments" in s for s in stages)

    def test_polling_never_goes_faster_than_ncbi_asks(self) -> None:
        intervals: list[float] = []
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append, poll_seconds=1.0)
        assert intervals and min(intervals) >= MINIMUM_POLL_SECONDS


class TestFailureIsNotAnEmptyResult:
    def test_waiting_for_ever_is_bounded_and_says_the_designs_are_unaffected(self) -> None:
        client, _ = scripted([PUT_REPLY, WAITING])
        clock = iter([0.0, 10.0, 99_999.0, 99_999.0])
        with pytest.raises(BlastError, match="had not finished"):
            run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                          sleep=never_sleep, now=lambda: next(clock), max_wait=60)

    def test_a_run_with_nothing_passing_is_refused_rather_than_screened_emptily(self) -> None:
        client, _ = scripted([PUT_REPLY])
        with pytest.raises(ValueError, match="nothing to screen"):
            run_screening("run-1", [{"id": "x", "rank": None, "target_site_dna": "ACGT"}],
                          TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep)

    def test_the_genomic_scope_refuses_to_run_without_a_way_to_place_its_hits(self) -> None:
        # A genomic hit with no nearest-TSS lookup is a coordinate with no meaning.
        client, _ = scripted([PUT_REPLY])
        with pytest.raises(ValueError, match="transcript lookup"):
            run_screening("run-1", CANDIDATES, GENOMIC_SCOPE, 19, client, sleep=never_sleep)


class TestTheGenomicScope:
    def test_hits_are_placed_against_annotated_start_sites(self) -> None:
        genomic = REPORT.replace(b"NM_006920", b"NC_000002.12")
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, genomic])

        def lookup(_interval: Any) -> list[dict[str, Any]]:  # noqa: ANN401
            return [{"id": "ENST1", "external_name": "OTHER", "start": 1055,
                     "end": 9000, "strand": 1}]

        record = run_screening(
            "run-1", CANDIDATES, GENOMIC_SCOPE, 19, client, transcripts=lookup,
            window=(-1000, -100), sleep=never_sleep,
        )
        assert record.parameters.database == "refseq_genomic"
        assert len(record.proximities) == len(record.hits)
        # Hit at 455, TSS at 1055 -> 600 bp upstream, inside Wang's -1000..-100 window.
        first = record.proximities[0]
        assert first.distance_to_tss == -600
        assert first.inside_run_window is True


class TestTheTable:
    def test_a_candidate_with_no_hits_is_present_with_zeroes(self) -> None:
        # Absent would read as "not screened", which is a different statement.
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, REPORT])
        record = run_screening(
            "run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep
        )
        rows = {r["candidate"]: r for r in specificity_table(record, 19)}
        assert set(rows) == {"saRNA-1549", "saRNA-0002"}
        assert rows["saRNA-0002"]["hit_count"] == 0
        assert rows["saRNA-0002"]["best_evalue"] is None
        assert rows["saRNA-1549"]["perfect_hits"] == 1

    def test_the_table_offers_columns_rather_than_a_score(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, REPORT])
        record = run_screening(
            "run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep
        )
        row = specificity_table(record, 19)[0]
        for column in ("hit_count", "best_evalue", "near_perfect_hits", "perfect_hits",
                       "longest_complementary_stretch", "in_run_window_hits"):
            assert column in row, column
        assert "score" not in row

    def test_a_truncated_screening_says_how_many_it_left_out(self) -> None:
        many = [
            {"id": f"saRNA-{i:04d}", "rank": i + 1, "target_site_dna": "ACGTACGTACGTACGTACG"}
            for i in range(MAX_CANDIDATES + 10)
        ]
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", many, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep,
            limit=MAX_CANDIDATES,
        )
        assert len(record.candidates_searched) == MAX_CANDIDATES
        # The unscreened must not read as screened clean, which is the whole point of saying so.
        assert record.notes and "50 of the 60" in record.notes[0]
        assert "not the same as having been searched clean" in record.notes[0]


class TestWaitingIsPacedByNcbisOwnEstimate:
    """NCBI replies to every submission with ``RTOE``: how long it expects the search to take.

    That number was parsed, exported and then ignored. The loop slept a flat fifteen seconds
    from the moment of submission, so for the seven-to-fifteen-minute searches ADR 0019
    measured it made twenty-eight to sixty requests, most of them before the search could
    possibly have finished. Every one is a request against a shared public service that the
    same ADR notes can block the whole institution's address.
    """

    def test_the_first_wait_is_ncbis_estimate_rather_than_the_floor(self) -> None:
        # The stubbed submission answers RTOE = 27.
        intervals: list[float] = []
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append)
        assert intervals[0] == pytest.approx(27.0)

    def test_an_estimate_below_the_floor_does_not_beat_the_rate_limit(self) -> None:
        impatient = PUT_REPLY.replace(b"RTOE = 27", b"RTOE = 1")
        intervals: list[float] = []
        client, _ = scripted([impatient, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append)
        assert min(intervals) >= MINIMUM_POLL_SECONDS

    def test_a_wild_estimate_cannot_consume_the_whole_budget(self) -> None:
        # If NCBI says an hour, waiting an hour before looking once would turn a search that
        # finishes early into a search that reports nothing until the deadline.
        wild = PUT_REPLY.replace(b"RTOE = 27", b"RTOE = 36000")
        intervals: list[float] = []
        client, _ = scripted([wild, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append, max_wait=600)
        assert intervals[0] <= 600 / 2

    def test_later_polls_return_to_the_floor_so_completion_is_noticed_promptly(self) -> None:
        # Honouring the estimate must not become a backoff: a job that overruns its estimate
        # should still be collected within one poll of finishing.
        intervals: list[float] = []
        client, _ = scripted([PUT_REPLY, WAITING, WAITING, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append)
        assert intervals[0] == pytest.approx(27.0)
        assert intervals[1:] == [MINIMUM_POLL_SECONDS] * (len(intervals) - 1)

    def test_no_estimate_falls_back_to_the_floor(self) -> None:
        silent = PUT_REPLY.replace(b"RTOE = 27", b"")
        intervals: list[float] = []
        client, _ = scripted([silent, READY_NO_HITS])
        run_screening("run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client,
                      sleep=intervals.append)
        assert intervals[0] == pytest.approx(MINIMUM_POLL_SECONDS)

    def test_the_requests_saved_are_counted_against_the_old_flat_interval(self) -> None:
        # The measurement this change exists for, on NCBI's own estimate for this submission.
        estimate, floor = 27.0, MINIMUM_POLL_SECONDS
        polls_before = estimate / floor
        assert polls_before > 1, "a flat interval polls before the search can have finished"


MANY: list[dict[str, Any]] = [
    {"id": f"saRNA-{i:04d}", "rank": i + 1, "target_site_dna": "ACGTACGTACGTACGTACG"}
    for i in range(40)
]


class TestScreeningInBatches:
    """ADR 0023: a fixed top 50 becomes the researcher's choice, resolved to one id list.

    Every button the interface offers, a count, the next fifteen, the ticked selection, comes
    through the same ``only`` list, so there is one selection path to get right rather than
    one per button.
    """

    def test_a_count_screens_the_best_that_many(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", MANY, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep, limit=15
        )
        assert len(record.candidates_searched) == 15
        assert record.candidates_searched[0] == "saRNA-0000"

    def test_an_explicit_selection_screens_exactly_those(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", MANY, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep,
            only=["saRNA-0030", "saRNA-0007"],
        )
        # Rank order regardless of the order asked in, because the FASTA and the table should
        # not depend on which checkbox was ticked first.
        assert record.candidates_searched == ("saRNA-0007", "saRNA-0030")

    def test_the_next_batch_is_just_another_selection(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        already = {c["id"] for c in MANY[:15]}
        rest = [c["id"] for c in MANY if c["id"] not in already][:15]
        record = run_screening(
            "run-1", MANY, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep, only=rest
        )
        assert record.candidates_searched[0] == "saRNA-0015"
        assert len(record.candidates_searched) == 15

    def test_one_job_stays_capped_however_much_is_asked_for(self) -> None:
        # The per-job ceiling is a courtesy to NCBI and is not negotiable from the interface.
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", MANY * 3, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep, limit=10_000
        )
        assert len(record.candidates_searched) <= MAX_CANDIDATES

    def test_a_partial_batch_says_the_rest_were_not_searched_clean(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        record = run_screening(
            "run-1", MANY, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep, limit=15
        )
        assert any("not the same as having been searched clean" in n for n in record.notes)

    def test_a_selection_of_only_excluded_candidates_is_refused_with_a_reason(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_NO_HITS])
        with pytest.raises(ValueError, match="none of the chosen candidates"):
            run_screening(
                "run-1", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep,
                only=["saRNA-9999"],
            )


GENOMIC_CANDIDATES: list[dict[str, Any]] = [
    {"id": "saRNA-1549", "rank": 1, "target_site_dna": "GTGAAGCCATAGATCTGTA",
     "genomic": {"contig": "2", "start": 454, "end": 473, "strand": "+"}},
]

GENOMIC_REPORT = b"""<?xml version="1.0"?>
<BlastOutput><BlastOutput_iterations><Iteration>
  <Iteration_query-def>saRNA-1549</Iteration_query-def>
  <Iteration_hits>
    <Hit><Hit_accession>NC_000002.12</Hit_accession><Hit_def>Homo sapiens chromosome 2</Hit_def>
      <Hit_hsps><Hsp><Hsp_identity>19</Hsp_identity><Hsp_align-len>19</Hsp_align-len>
      <Hsp_gaps>0</Hsp_gaps><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>19</Hsp_query-to>
      <Hsp_hit-from>455</Hsp_hit-from><Hsp_hit-to>473</Hsp_hit-to><Hsp_hit-frame>1</Hsp_hit-frame>
      <Hsp_evalue>0.004</Hsp_evalue><Hsp_bit-score>38.2</Hsp_bit-score></Hsp></Hit_hsps></Hit>
    <Hit><Hit_accession>NC_000007.14</Hit_accession><Hit_def>Homo sapiens chromosome 7</Hit_def>
      <Hit_hsps><Hsp><Hsp_identity>17</Hsp_identity><Hsp_align-len>19</Hsp_align-len>
      <Hsp_gaps>0</Hsp_gaps><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>19</Hsp_query-to>
      <Hsp_hit-from>900</Hsp_hit-from><Hsp_hit-to>918</Hsp_hit-to><Hsp_hit-frame>1</Hsp_hit-frame>
      <Hsp_evalue>0.5</Hsp_evalue><Hsp_bit-score>30.1</Hsp_bit-score></Hsp></Hit_hsps></Hit>
  </Iteration_hits>
</Iteration></BlastOutput_iterations></BlastOutput>"""


class TestTheCandidatesOwnLocusIsNotAFinding:
    """Without this, every genomic candidate reports a perfect match and the column is noise."""

    def test_the_own_locus_hit_is_removed_from_the_findings(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        record = run_screening(
            "run-1", GENOMIC_CANDIDATES, GENOMIC_SCOPE, 19, client, sleep=never_sleep,
            transcripts=lambda _interval: [],
        )
        assert [h.accession for h in record.hits] == ["NC_000007.14"]
        assert [h.accession for h in record.on_target_hits] == ["NC_000002.12"]

    def test_the_perfect_match_column_is_no_longer_one_for_everything(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        record = run_screening(
            "run-1", GENOMIC_CANDIDATES, GENOMIC_SCOPE, 19, client, sleep=never_sleep,
            transcripts=lambda _interval: [],
        )
        row = specificity_table(record, 19)[0]
        assert row["perfect_hits"] == 0, "the only perfect match was the candidate's own locus"
        assert row["on_target_hits"] == 1, "and the row says it was taken out"

    def test_the_removal_is_stated_on_the_record(self) -> None:
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        record = run_screening(
            "run-1", GENOMIC_CANDIDATES, GENOMIC_SCOPE, 19, client, sleep=never_sleep,
            transcripts=lambda _interval: [],
        )
        assert any("own target loci" in n for n in record.notes)

    def test_a_transcript_search_removes_nothing(self) -> None:
        # A promoter is not usually transcribed, and the only handle in a transcript report
        # would be a gene name, which is a guess this project will not make silently.
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        record = run_screening(
            "run-1", GENOMIC_CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep
        )
        assert record.on_target_hits == ()
        assert len(record.hits) == 2

    def test_a_run_without_coordinates_says_it_could_not_tell(self) -> None:
        # Pasted sequence has no locus to compare against. Saying nothing would let a reader
        # assume the removal happened.
        pasted = [{"id": "saRNA-1549", "rank": 1, "target_site_dna": "GTGAAGCCATAGATCTGTA"}]
        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        record = run_screening(
            "run-1", pasted, GENOMIC_SCOPE, 19, client, sleep=never_sleep,
            transcripts=lambda _interval: [],
        )
        assert record.on_target_hits == ()
        assert any("no genomic coordinates" in n for n in record.notes)

    def test_placement_is_not_spent_on_a_hit_that_was_removed(self) -> None:
        placed: list[GenomicInterval] = []

        def lookup(interval: GenomicInterval) -> list[dict[str, Any]]:
            placed.append(interval)
            return []

        client, _ = scripted([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])
        run_screening(
            "run-1", GENOMIC_CANDIDATES, GENOMIC_SCOPE, 19, client, sleep=never_sleep,
            transcripts=lookup,
        )
        assert len(placed) == 1, "only the remaining off-target hit is placed"
