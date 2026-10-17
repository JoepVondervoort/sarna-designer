"""Screening actions through real HTTP: include the media/origin checks router mocks bypass."""

from pathlib import Path
from threading import Event, Thread

import httpx

from sarna.app.web import AppRouter, make_server
from sarna.offtarget.models import TRANSCRIPT_SCOPE
from tests.unit.test_app_screening import completed_run


def test_resume_crosses_the_http_boundary_without_resubmitting(tmp_path: Path) -> None:
    service, run_id, _ = completed_run(tmp_path)
    screening_id, _ = service.store.reserve_screening(run_id, TRANSCRIPT_SCOPE, ["saRNA-001"])
    service.store.pause_screening(screening_id, "EXISTING-RID", "Still waiting", ["saRNA-001"])
    collected = Event()
    calls = []

    def collect(identifier: str) -> None:
        calls.append(identifier)
        service.store.finish_screening(identifier, {"request_id": "EXISTING-RID",
                                                    "candidates_searched": ["saRNA-001"]})
        collected.set()

    service.resume = collect  # type: ignore[method-assign]
    router = AppRouter(service)
    server = make_server(router, 0)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        with httpx.Client(base_url=base, trust_env=False, timeout=5) as client:
            # This is the browser's previous request. It must explain the rejection in JSON.
            rejected = client.post(f"/api/screenings/{screening_id}/resume")
            assert rejected.status_code == 415
            assert rejected.json() == {"detail": "Submit a JSON request", "code": "http-415"}
            assert service.store.screening(screening_id)["status"] == "waiting"
            assert not calls
            # Fixed browser request: existing job is collected, not created again.
            accepted = client.post(f"/api/screenings/{screening_id}/resume", json={})
            assert accepted.status_code == 202
            assert accepted.headers["Location"] == f"/api/screenings/{screening_id}"
            assert accepted.headers["Retry-After"] == "5"
            assert collected.wait(2)
            result = client.get(accepted.headers["Location"])
            assert result.json()["status"] == "complete"
            assert result.json()["request_id"] == "EXISTING-RID"
            assert calls == [screening_id]
            assert len(service.store.screenings_for(run_id)) == 1
            denied = client.post(f"/api/screenings/{screening_id}/resume", json={},
                                 headers={"Origin": "https://unrelated.example"})
            assert denied.status_code == 403
            assert denied.headers["Content-Type"].startswith("application/json")
            assert calls == [screening_id]
    finally:
        server.shutdown()
        thread.join(timeout=5)
        server.server_close()
        router.close()
