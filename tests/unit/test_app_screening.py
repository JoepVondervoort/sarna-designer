"""Off-target screening through the app: queued, polled, stored as its own dated record.

The design run is complete before any of this happens and is never modified by it. That is the
property these tests exist to hold: a screening is an *addition* to the evidence about a run,
with its own provenance, not an edit to the frozen record (ADR 0019).
"""

import json
import threading
import time
from pathlib import Path
from typing import Any

import httpx
import pytest

from sarna.app.service import DesignService
from sarna.app.web import AppRouter
from sarna.offtarget.models import GENOMIC_SCOPE, TRANSCRIPT_SCOPE
from sarna.offtarget.ncbi import BlastClient
from tests.unit.test_app import DNA, run_request, service_for
from tests.unit.test_offtarget_ncbi import PUT_REPLY, READY_WITH_HITS, REPORT


def blast_answering(replies: list[bytes]) -> BlastClient:
    """A BLAST client replying from a script, holding the last reply."""
    seen = {"n": 0}

    def handler(_request: httpx.Request) -> httpx.Response:
        seen["n"] += 1
        return httpx.Response(200, content=replies[min(seen["n"] - 1, len(replies) - 1)])

    return BlastClient(client=httpx.Client(transport=httpx.MockTransport(handler)))


def screen_now(service: DesignService, run_id: str, scope: str = TRANSCRIPT_SCOPE) -> str:
    """Run a screening synchronously, as the background worker would."""
    screening_id = service.store.start_screening(run_id, scope)
    # The real client sleeps between polls; this reply is ready immediately.
    service.screen(screening_id, run_id, scope,
                   blast_client=blast_answering([PUT_REPLY, READY_WITH_HITS, REPORT]),
                   sleep=lambda _seconds: None)
    return screening_id


def completed_run(tmp_path: Path) -> tuple[DesignService, str, dict[str, Any]]:
    """Return a service with one finished run."""
    service = service_for(tmp_path)
    saved = run_request(service)
    return service, saved["app"]["run_id"], saved


class TestTheScreeningIsASeparateRecord:
    def test_a_screening_does_not_touch_the_design_snapshot(self, tmp_path: Path) -> None:
        service, run_id, before = completed_run(tmp_path)
        screen_now(service, run_id)
        assert service.store.snapshot(run_id) == before

    def test_it_carries_its_own_date_database_and_parameters(self, tmp_path: Path) -> None:
        service, run_id, _ = completed_run(tmp_path)
        record = service.store.screening(screen_now(service, run_id))
        assert record["status"] == "complete"
        result = record["result"]
        assert result["parameters"]["database"] == "refseq_rna"
        assert result["parameters"]["word_size"] == 7
        assert result["searched_at"]
        assert result["request_id"] == "ABC123XYZ014"

    def test_screening_twice_keeps_both_records(self, tmp_path: Path) -> None:
        # NCBI's databases move. A later screening is a second observation, not a correction.
        service, run_id, _ = completed_run(tmp_path)
        first, second = screen_now(service, run_id), screen_now(service, run_id)
        assert first != second
        listed = service.store.screenings_for(run_id)
        assert {s["id"] for s in listed} == {first, second}
        assert all(s["status"] == "complete" for s in listed)

    def test_both_scopes_can_be_screened_for_one_run(self, tmp_path: Path) -> None:
        service, run_id, _ = completed_run(tmp_path)
        service.store.start_screening(run_id, TRANSCRIPT_SCOPE)
        service.store.start_screening(run_id, GENOMIC_SCOPE)
        assert {s["scope"] for s in service.store.screenings_for(run_id)} == {
            TRANSCRIPT_SCOPE, GENOMIC_SCOPE
        }


class TestTheResultTheInterfaceReads:
    def test_the_specificity_table_has_a_row_for_every_screened_candidate(
        self, tmp_path: Path
    ) -> None:
        service, run_id, saved = completed_run(tmp_path)
        result = service.store.screening(screen_now(service, run_id))["result"]
        passing = [c["id"] for c in saved["candidates"] if c["rank"] is not None]
        assert {r["candidate"] for r in result["specificity"]} == set(passing)

    def test_a_candidate_with_no_hits_is_reported_rather_than_omitted(
        self, tmp_path: Path
    ) -> None:
        service, run_id, _ = completed_run(tmp_path)
        result = service.store.screening(screen_now(service, run_id))["result"]
        empty = [r for r in result["specificity"] if r["hit_count"] == 0]
        assert empty, "the recorded report has hits for only one query"
        assert empty[0]["best_evalue"] is None
        assert "does not establish specificity" in empty[0]["summary"]

    def test_no_row_carries_a_combined_specificity_score(self, tmp_path: Path) -> None:
        service, run_id, _ = completed_run(tmp_path)
        result = service.store.screening(screen_now(service, run_id))["result"]
        for row in result["specificity"]:
            assert "score" not in row


class TestFailureKeepsTheDesigns:
    def test_an_unreachable_ncbi_fails_the_screening_and_nothing_else(
        self, tmp_path: Path
    ) -> None:
        service, run_id, before = completed_run(tmp_path)

        def refuse(_request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("no route to host")

        screening_id = service.store.start_screening(run_id, TRANSCRIPT_SCOPE)
        service.screen(screening_id, run_id, TRANSCRIPT_SCOPE,
                       blast_client=BlastClient(
                           client=httpx.Client(transport=httpx.MockTransport(refuse))),
                       sleep=lambda _seconds: None)
        record = service.store.screening(screening_id)
        assert record["status"] == "error"
        assert "could not reach NCBI" in record["error"]
        # The run is untouched and still readable.
        assert service.store.snapshot(run_id) == before
        assert service.store.status(run_id)["status"] == "complete"


class TestTheRoutes:
    def test_a_screening_is_queued_and_then_readable(self, tmp_path: Path) -> None:
        service, run_id, _ = completed_run(tmp_path)
        router = AppRouter(service)
        try:
            status, _, body = router.dispatch(
                "POST", f"/api/runs/{run_id}/screenings",
                json.dumps({"scope": "transcript"}).encode(),
            )
            assert status == 202
            screening_id = json.loads(body)["id"]
            status, _, body = router.dispatch("GET", f"/api/screenings/{screening_id}", b"")
            assert status == 200
            assert json.loads(body)["run_id"] == run_id
            listed = json.loads(
                router.dispatch("GET", f"/api/runs/{run_id}/screenings", b"")[2])
            assert [s["id"] for s in listed] == [screening_id]
        finally:
            router.close()

    def test_an_unknown_scope_is_refused_with_the_choices_named(self, tmp_path: Path) -> None:
        service, run_id, _ = completed_run(tmp_path)
        router = AppRouter(service)
        try:
            status, _, body = router.dispatch(
                "POST", f"/api/runs/{run_id}/screenings", b'{"scope":"everything"}')
            assert status == 422
            assert "transcript" in body.decode()
        finally:
            router.close()

    def test_screening_a_run_with_nothing_passing_is_refused(self, tmp_path: Path) -> None:
        service = service_for(tmp_path)
        # A GC floor no candidate in this promoter reaches leaves nothing to order, so
        # nothing to screen. Valid as a configuration, unsatisfiable on this sequence.
        saved = run_request(service, overrides={
            "gc-content": {"minimum-gc-fraction": 0.99, "maximum-gc-fraction": 1.0}})
        router = AppRouter(service)
        try:
            assert saved["counts"]["passing"] == 0
            status, _, body = router.dispatch(
                "POST", f"/api/runs/{saved['app']['run_id']}/screenings", b"{}")
            assert status == 422
            assert "nothing to screen" in body.decode()
        finally:
            router.close()

    def test_an_unknown_screening_is_a_clean_not_found(self, tmp_path: Path) -> None:
        service, _, _ = completed_run(tmp_path)
        router = AppRouter(service)
        try:
            assert router.dispatch("GET", "/api/screenings/nope", b"")[0] == 404
        finally:
            router.close()


def test_the_order_csv_still_says_screening_was_not_performed_by_itself(
    tmp_path: Path,
) -> None:
    # The order export's sanity report is generated from the run alone. A screening lives in its
    # own record, so until that is wired into the export the honest answer stays "not performed".
    service, run_id, _ = completed_run(tmp_path)
    screen_now(service, run_id)
    router = AppRouter(service)
    try:
        body = router.dispatch("GET", f"/api/runs/{run_id}/export.order.csv", b"")[2]
        assert "off-target-screening [not-performed]" in body.decode()
    finally:
        router.close()


@pytest.mark.parametrize("scope", [TRANSCRIPT_SCOPE, GENOMIC_SCOPE])
def test_every_supported_scope_has_a_database(scope: str) -> None:
    from sarna.offtarget.ncbi import parameters_for

    assert parameters_for(scope).database


class TestAScreeningMustNotStarveTheRestOfTheApp:
    """A screening waits on NCBI for minutes; a design takes under a second.

    ADR 0019's whole premise is that the researcher reads their candidates while the search
    runs. That only holds if the search is not occupying the capacity the reading needs. Both
    scopes are supported and the user asked for both, so two screenings in flight is the
    expected state, not an edge case.
    """

    @staticmethod
    def _block_screening(
        monkeypatch: pytest.MonkeyPatch, service: DesignService, release: threading.Event
    ) -> list[str]:
        """Replace screening with one that blocks until released, recording what started."""
        started: list[str] = []

        def blocked(*args: object, **kwargs: object) -> None:
            started.append(str(kwargs.get("scope", args[2] if len(args) > 2 else "?")))
            release.wait(timeout=60)

        monkeypatch.setattr(service, "screen", blocked)
        return started

    def test_two_screenings_in_flight_do_not_delay_a_new_design(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        service, run_id, _ = completed_run(tmp_path)
        release = threading.Event()
        started = self._block_screening(monkeypatch, service, release)
        router = AppRouter(service)
        try:
            for scope in (TRANSCRIPT_SCOPE, GENOMIC_SCOPE):
                status, _, _ = router.dispatch(
                    "POST", f"/api/runs/{run_id}/screenings",
                    json.dumps({"scope": scope}).encode())
                assert status == 202
            deadline = time.monotonic() + 10
            while len(started) < 2 and time.monotonic() < deadline:
                time.sleep(0.02)
            assert len(started) == 2, f"only {started} reached NCBI"

            # Both searches are now sitting in NCBI's queue. A design must still run.
            body = json.dumps({"source": "sequence", "sequence": DNA,
                               "regulatory": False, "motifs": False}).encode()
            status, _, payload = router.dispatch("POST", "/api/runs", body)
            assert status == 202
            new_run = json.loads(payload)["id"]
            deadline = time.monotonic() + 20
            while time.monotonic() < deadline:
                state = service.store.status(new_run)["status"]
                if state != "running":
                    break
                if service.store.status(new_run)["stage"] != "Queued":
                    break
                time.sleep(0.05)
            stage = service.store.status(new_run)["stage"]
            assert stage != "Queued", (
                "a design submitted while two screenings waited on NCBI never left the queue; "
                "the screening work is occupying the pool the interactive work needs"
            )
        finally:
            release.set()
            router.close()
