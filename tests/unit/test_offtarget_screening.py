"""Running a screening end to end, with the clock and the network both injected.

No test here sleeps or reaches the network. What they pin is the behaviour a researcher depends
on while waiting: that the whole run goes in one job, that each stage is reported, that waiting
is bounded, and that a candidate with nothing found is still present in the table rather than
missing from it.
"""

from typing import Any

import httpx
import pytest

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
        record = run_screening("run-1", many, TRANSCRIPT_SCOPE, 19, client, sleep=never_sleep)
        assert len(record.candidates_searched) == MAX_CANDIDATES
        assert record.notes and "out of" in record.notes[0]
