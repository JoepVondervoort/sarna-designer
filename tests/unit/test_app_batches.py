"""Verify separate scientific contexts, partial failure, archival groups and atomic saving."""

import json
import time
from pathlib import Path
from typing import Any

import pytest

from sarna.app.models import BatchRequest, DesignRequest
from sarna.app.service import DesignService
from sarna.app.store import RunStore
from sarna.app.web import AppRouter
from sarna.genome.ensembl import EnsemblClient
from tests.unit.test_app import service_for
from tests.unit.test_ensembl import (
    AUTO_SEQUENCE,
    GENE_PAYLOAD,
    GENE_ROUTE,
    TRANSCRIPTS_PAYLOAD,
    TRANSCRIPTS_ROUTE,
    client_for,
    requests_made,
)


def child_request(**changes: Any) -> dict[str, Any]:  # noqa: ANN401
    """Use a short recorded window; relax four gates so local ranks can be inspected."""
    return {
        "symbol": "SCN1A", "transcript_id": "ENST00000674923",
        "upstream": 25, "downstream": 5, "regulatory": False, "motifs": False,
        "roles": {"tss-window": "off", "gc-content": "off", "pos19-a": "off",
                  "thermodynamic-asymmetry": "off"},
        **changes,
    }


def finished(router: AppRouter, batch_id: str) -> dict[str, Any]:
    """Wait only for the actual small jobs, with bounded time and no external sources."""
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        code, _, body = router.dispatch("GET", f"/api/batches/{batch_id}")
        assert code == 200
        data: dict[str, Any] = json.loads(body)
        if data["status"] != "running":
            return data
        time.sleep(0.01)
    raise AssertionError("Recorded-source batch did not finish")


def submit(router: AppRouter, requests: list[dict[str, Any]], label: str = "") -> dict[str, Any]:
    """Submit the same full requests used by the browser, then inspect persistent progress."""
    code, _, body = router.dispatch("POST", "/api/batches", json.dumps({
        "label": label, "requests": requests,
    }).encode())
    assert code == 202, body.decode()
    return finished(router, json.loads(body)["id"])


@pytest.mark.strand
def test_batch_keeps_two_selected_tss_contexts_and_independent_ranks(tmp_path: Path):
    service = service_for(tmp_path)
    router = AppRouter(service)
    try:
        batch = submit(router, [child_request(), child_request(
            transcript_id="ENST00000409050", upstream=26,
        )], label="Two SCN1A references")
        assert batch["status"] == "complete"
        assert len({run["id"] for run in batch["runs"]}) == 2
        first, second = [service.store.snapshot(run["id"]) for run in batch["runs"]]
        # Minus strand: upstream is toward higher genomic coordinates, once at the boundary.
        assert (first["promoter"]["start"], first["promoter"]["end"]) == (
            166128015, 166128045,
        )
        assert (second["promoter"]["start"], second["promoter"]["end"]) == (
            166149209, 166149240,
        )
        assert first["counts"]["generated"] == 12
        assert second["counts"]["generated"] == 13
        for result, run in zip((first, second), batch["runs"], strict=True):
            ranks = sorted(c["rank"] for c in result["candidates"] if c["rank"] is not None)
            assert ranks and ranks == list(range(1, len(ranks) + 1))
            assert run["summary"]["counts"] == result["counts"]
            assert run["summary"]["anchor"] == result["app"]["anchor"]
            assert run["summary"]["best_rating"] == next(
                c["rating"] for c in result["candidates"] if c["rank"] == 1
            )
        assert "rank" not in batch  # There is no pooled cross-promoter rank.
        assert "ENST00000409050" in batch["runs"][1]["label"]
        assert first["app"]["request"]["upstream"] == 25
        assert second["app"]["request"]["upstream"] == 26
        assert json.loads(router.dispatch("GET", "/api/batches")[2])[0]["id"] == batch["id"]
    finally:
        router.close()


@pytest.mark.parametrize("invalid", [
    {"requests": []},
    {"requests": [child_request(), child_request(symbol="CHD2")]},
    {"requests": [child_request(), child_request(transcript_id=None)]},
    {"requests": [child_request(), child_request(overrides={
        "gc-content": {"minimum-gc-fraction": 0.8},
    }, roles={})]},
    {"requests": [child_request(), child_request(preset_id="unknown-preset")]},
    {"requests": [child_request()] * 33},
])
def test_invalid_batch_creates_neither_parent_nor_child_jobs(
    tmp_path: Path, invalid: dict[str, Any],
):
    router = AppRouter(service_for(tmp_path))
    try:
        assert router.dispatch("POST", "/api/batches", json.dumps(invalid).encode())[0] == 422
        assert router.service.store.list_batches() == []
        assert router.service.store.list_runs() == []
    finally:
        router.close()


def test_partial_batch_retry_reuses_success_and_never_rewrites_previous_results(tmp_path: Path):
    cache = tmp_path / "cache"
    routes = {GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
              "/sequence/region/": {"seq": AUTO_SEQUENCE}}
    clients: list[EnsemblClient] = []

    def factory(offline: bool) -> EnsemblClient:
        client = client_for(routes, cache_dir=cache, offline=offline)
        clients.append(client)
        return client

    service = DesignService(tmp_path / "runs", factory)
    # Warm only the first reference; the second window is deliberately absent.
    request = DesignRequest.model_validate(child_request())
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    router = AppRouter(service)
    try:
        original = submit(router, [child_request(offline=True), child_request(
            transcript_id="ENST00000409050", offline=True,
        )])
        assert original["status"] == "partial"
        success, failure = original["runs"]
        assert success["status"] == "complete"
        assert failure["error_code"] == "offline-cache-miss"
        assert all(requests_made(client) == [] for client in clients[1:])
        saved = service.store.snapshot(success["id"])
        before = len(clients)
        code, _, body = router.dispatch("POST", f"/api/batches/{original['id']}/retry",
                                        b'{"offline":false}')
        assert code == 202, body.decode()
        retry = finished(router, json.loads(body)["id"])
        assert retry["status"] == "complete"
        assert retry["retry_of"] == original["id"]
        assert retry["runs"][0]["id"] == success["id"]
        assert retry["runs"][1]["id"] != failure["id"]
        assert len(clients) == before + 1  # Only the failed child ran again.
        assert len(requests_made(clients[-1])) == 1  # Gene/transcripts reused; new window fetched.
        assert service.store.snapshot(success["id"]) == saved
        assert service.store.batch(original["id"])["status"] == "partial"
        assert service.store.status(failure["id"])["status"] == "error"
    finally:
        router.close()


def test_batch_archive_contains_all_verdicts_and_failure_metadata_without_source_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    service = service_for(tmp_path)
    router = AppRouter(service)
    try:
        batch = submit(router, [child_request(), child_request(transcript_id="NOT_A_TRANSCRIPT")],
                       label="</script><script>injected()</script>")
        assert batch["status"] == "partial"

        def forbidden(*args: Any, **kwargs: Any) -> None:  # noqa: ANN401
            raise AssertionError("Saved groups must not query sources or rerun science")

        monkeypatch.setattr(service, "client_factory", forbidden)
        monkeypatch.setattr(service, "run", forbidden)
        code, _, raw = router.dispatch("GET", f"/api/batches/{batch['id']}/export.json")
        assert code == 200
        archived = json.loads(raw)
        success, failure = archived["runs"]
        assert all(candidate["rules"] for candidate in success["result"]["candidates"])
        assert failure["status"] == "error" and failure["error"]
        assert "result" not in failure
        code, _, raw = router.dispatch("GET", f"/api/batches/{batch['id']}/export.html")
        assert code == 200
        assert b"window.SARNA_FROZEN=" in raw and b'"batch":' in raw
        assert b"</script><script>injected()" not in raw
        assert b'<script src="/app.js"' not in raw
        assert b'<link rel="stylesheet"' not in raw
        assert b'class="brand" href="/"' not in raw
        assert router.dispatch("GET", f"/api/batches/{batch['id']}")[0] == 200
        assert router.dispatch("GET", "/api/batches/missing")[0] == 404
    finally:
        router.close()


def test_store_group_creation_rolls_back_if_reuse_is_not_a_completed_matching_run(tmp_path: Path):
    store = RunStore(tmp_path / "runs.sqlite3")
    request = BatchRequest.model_validate({"requests": [child_request(), child_request()]})
    with pytest.raises(ValueError, match="unchanged completed"):
        store.create_batch(request.model_dump(mode="json"), reuse={1: "missing"})
    assert store.list_batches() == []
    assert store.list_runs() == []


def test_interrupted_batch_reopens_with_failure_metadata(tmp_path: Path):
    store = RunStore(tmp_path / "runs.sqlite3")
    request = BatchRequest.model_validate({"requests": [child_request()]})
    batch_id, _, _ = store.create_batch(request.model_dump(mode="json"))
    reopened = RunStore(store.path)
    batch = reopened.batch(batch_id)
    assert batch["status"] == "error"
    assert batch["runs"][0]["stage"] == "Interrupted"
    assert "Application stopped" in batch["runs"][0]["error"]
