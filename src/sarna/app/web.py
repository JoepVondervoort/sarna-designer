"""Loopback-only HTTP adapter; the Python engine owns every scientific result."""

import csv
import io
import json
import logging
import re
from base64 import b64encode
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from typing import Any, cast
from urllib.parse import parse_qs, unquote, urlsplit

from pydantic import ValidationError

from sarna import __version__, settings
from sarna.app.models import BatchRequest, DesignRequest
from sarna.app.service import DesignService
from sarna.design.order import (
    CSV_COLUMNS,
    DEFAULT_MODEL,
    csv_row,
    order_from_record,
)
from sarna.design.presets import PRESETS
from sarna.genome.ensembl import (
    EnsemblCacheMissError,
    EnsemblUnavailableError,
    GeneNotFoundError,
)
from sarna.offtarget.models import GENOMIC_SCOPE, TRANSCRIPT_SCOPE
from sarna.regulatory.evidence import CELL_TYPES
from sarna.rules.catalog import CATALOGUE
from sarna.rules.resolved import ResolvedRule

STATIC = resources.files("sarna.app").joinpath("static")
MAX_BODY = 200_000


def frozen_page(data: dict[str, Any], snapshot: dict[str, Any] | None) -> bytes:
    """Embed a single run or whole group in the same UI with all local assets and verdicts."""
    page = STATIC.joinpath("index.html").read_text(encoding="utf-8")
    css = STATIC.joinpath("app.css").read_text(encoding="utf-8")
    js = STATIC.joinpath("app.js").read_text(encoding="utf-8")
    data["catalogue"] = {
        "engine_version": snapshot["app"]["engine_version"] if snapshot else __version__,
        "rules": (snapshot["app"].get("available_rules") or list({
            rule["id"]: rule for rule in (
                *snapshot["configuration"]["baseline_rules"],
                *snapshot["configuration"]["rules"],
            )
        }.values())) if snapshot else [],
        "presets": [], "cells": [],
    }
    encoded = json.dumps(data).replace("</", "<\\/")
    page = page.replace('<link rel="stylesheet" href="/app.css">', f"<style>{css}</style>")
    page = page.replace('<script src="/app.js" defer></script>',
                        f"<script>window.SARNA_FROZEN={encoded};</script><script>{js}</script>")
    icon = b64encode(STATIC.joinpath("favicon.svg").read_bytes()).decode()
    page = page.replace('href="/favicon.svg"', f'href="data:image/svg+xml;base64,{icon}"')
    page = page.replace('<a class="brand" href="/"', '<a class="brand" href="#"')
    return page.encode()


def _selected_for_export(
    snapshot: dict[str, Any], query: dict[str, list[str]]
) -> list[dict[str, Any]]:
    """Return the candidates a CSV export should write, **in rank order**.

    A snapshot's ``candidates`` are in *position* order carrying their ranks, which is what the
    promoter map draws along and what the browser table sorts for itself. A CSV is neither: it is
    a flat file read top-down and pasted into an order form, so leaving it in position order puts
    rank 8 on the first line and makes "take the top three" select the wrong three molecules.

    Explicitly selected ids are honoured as given, including rejected candidates, which have no
    rank and sort last.

    Raises:
        ValueError: If an id does not belong to this run. Silently dropping it would hand back a
            short file that looks complete.
    """
    ids = set(query.get("ids", [""])[0].split(",")) - {""}
    selected = [c for c in snapshot["candidates"] if c["id"] in ids] if ids else [
        c for c in snapshot["candidates"] if c["rank"] is not None
    ]
    if ids - {c["id"] for c in selected}:
        raise ValueError("A selected candidate does not belong to this run.")
    return sorted(
        selected, key=lambda c: (c["rank"] is None, c["rank"] or 0, c["id"])
    )


def catalogue() -> dict[str, Any]:
    """Publish editable rule declarations and frozen preset baselines, without recalculation."""
    return {
        "engine_version": __version__,
        "rules": [ResolvedRule(rule).to_dict() for rule in CATALOGUE],
        "presets": [{"id": p.id, "name": p.name, "description": p.description,
                     "configuration": p.resolve().to_dict()} for p in PRESETS.values()],
        "cells": [cell.to_dict() for cell in CELL_TYPES.values()],
        "data_access": {"offline_forced": settings.is_offline()},
    }


class AppRouter:
    """Route local requests independently of sockets, for direct offline integration tests."""

    def __init__(self, service: DesignService) -> None:
        """Hold a bounded worker pool and the local run service."""
        self.service = service
        self.workers = ThreadPoolExecutor(max_workers=2, thread_name_prefix="sarna")
        # Screenings get their own pool because they are a different kind of work. A design
        # finishes in under a second; a screening waits on NCBI's queue for seven to fifteen
        # minutes and may wait thirty. Sharing one two-worker pool meant that starting both
        # scopes -- which is the supported and expected use -- occupied every worker, and the
        # next design sat in "Queued" until a search came back. ADR 0019's premise is that the
        # researcher reads their candidates while the search runs, so the search must not hold
        # the capacity that reading needs. Two workers also caps this app at two concurrent
        # NCBI jobs, which is the other half of the same ADR's rate-limit concern.
        self.screeners = ThreadPoolExecutor(max_workers=2, thread_name_prefix="sarna-screen")

    def close(self) -> None:
        """Finish workers when the application or a test shuts down."""
        self.workers.shutdown(wait=True, cancel_futures=True)
        # Queued screenings are dropped and a running one is *not* waited for. A screening sits
        # in NCBI's queue for minutes and may sit for thirty, so waiting would make closing the
        # app -- or a test closing its router -- block for exactly that long. The screening's
        # record stays "running", which is the honest state: the job was abandoned, not failed,
        # and NCBI still holds the request id.
        self.screeners.shutdown(wait=False, cancel_futures=True)

    def dispatch(
        self, method: str, target: str, body: bytes = b"",
    ) -> tuple[int, str, bytes]:
        """Validate a request and return status, media type and response bytes."""
        try:
            return self._dispatch(method, target, body)
        except ValidationError as error:
            return self._json(422, {"detail": [{"loc": list(e["loc"]), "msg": e["msg"]}
                                             for e in error.errors()]})
        except (ValueError, TypeError) as error:
            return self._json(422, {"detail": str(error)})
        except KeyError as error:
            return self._json(404, {"detail": str(error)})
        except EnsemblCacheMissError as error:
            return self._json(409, {"detail": str(error), "code": "offline-cache-miss"})
        except GeneNotFoundError as error:
            return self._json(404, {"detail": str(error), "code": "gene-not-found"})
        except EnsemblUnavailableError as error:
            return self._json(502, {"detail": str(error), "code": "ensembl-unavailable"})
        except Exception as error:
            logging.getLogger(__name__).exception("Local request failed")
            return self._json(502, {"detail": str(error) or type(error).__name__})

    @staticmethod
    def _json(status: int, data: object) -> tuple[int, str, bytes]:
        return status, "application/json; charset=utf-8", json.dumps(data).encode()

    def _dispatch(self, method: str, target: str, body: bytes) -> tuple[int, str, bytes]:
        parsed = urlsplit(target)
        path = unquote(parsed.path)
        query = parse_qs(parsed.query)
        if method == "GET" and path in ("/", "/app.js", "/app.css", "/favicon.svg", "/favicon.ico"):
            name = {"/": "index.html", "/favicon.ico": "favicon.svg"}.get(path, path[1:])
            media = {"html": "text/html", "js": "text/javascript", "css": "text/css",
                     "svg": "image/svg+xml"}[name.rsplit(".", 1)[1]]
            return 200, media + "; charset=utf-8", STATIC.joinpath(name).read_bytes()
        if method == "GET" and path == "/api/catalogue":
            return self._json(200, catalogue())
        if method == "GET" and path.startswith("/api/genes/"):
            symbol = path.removeprefix("/api/genes/")
            if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", symbol):
                raise ValueError("Use a human gene symbol or Ensembl gene ID, without spaces.")
            return self._json(200, self.service.lookup(symbol, query.get("offline") == ["true"]))
        if path == "/api/batches":
            if method == "GET":
                return self._json(200, self.service.store.list_batches())
            if method == "POST":
                if len(body) > MAX_BODY:
                    raise ValueError("Request is too large.")
                return self._queue_batch(BatchRequest.model_validate_json(body))
        if path == "/api/runs":
            if method == "GET":
                return self._json(200, self.service.store.list_runs())
            if method == "POST":
                if len(body) > MAX_BODY:
                    raise ValueError("Request is too large.")
                request = DesignRequest.model_validate_json(body)
                try:
                    request.resolve()  # No job or lookup for an invalid custom configuration.
                except KeyError as error:
                    raise ValueError(str(error)) from error
                run_id = self.service.store.create(request.model_dump(mode="json"))
                self.workers.submit(self.service.run, run_id, request)
                return self._json(202, {"id": run_id})
        parts = path.strip("/").split("/")
        if method == "GET" and len(parts) == 3 and parts[:2] == ["api", "screenings"]:
            return self._json(200, self.service.store.screening(parts[2]))
        if (method == "POST" and len(parts) == 4
                and parts[:2] == ["api", "runs"] and parts[3] == "screenings"):
            return self._start_screening(parts[2], body)
        if len(parts) >= 3 and parts[:2] == ["api", "batches"]:
            batch = self.service.store.batch(parts[2])
            if method == "GET" and len(parts) == 3:
                return self._json(200, batch)
            if method == "POST" and len(parts) == 4 and parts[3] == "retry":
                if len(body) > MAX_BODY:
                    raise ValueError("Request is too large.")
                retry = json.loads(body)
                if not isinstance(retry, dict) or set(retry) != {"offline"} or (
                    retry["offline"] is not False
                ):
                    raise ValueError(
                        "Retry requires {offline: false} to allow reference retrieval."
                    )
                if settings.is_offline():
                    return self._json(409, {
                        "detail": "SARNA_OFFLINE forces cache-only mode. Unset it and restart "
                                  "the app before retrying online.",
                        "code": "offline-forced",
                    })
                if batch["status"] == "running":
                    raise ValueError("Wait for this group to finish before retrying its failures.")
                failed = [run for run in batch["runs"] if run["status"] == "error"]
                if not failed:
                    raise ValueError("This group has no failed references to retry.")
                requests = [run["request"] for run in batch["runs"]]
                reuse = {}
                for ordinal, run in enumerate(batch["runs"]):
                    if run["status"] == "complete":
                        reuse[ordinal] = run["id"]
                    else:
                        requests[ordinal] = {**requests[ordinal], "offline": False}
                retry_request = BatchRequest.model_validate({
                    "label": ("Retry · " + batch["label"])[:150], "requests": requests,
                })
                return self._queue_batch(retry_request, reuse=reuse, retry_of=batch["id"])
            if method == "GET" and len(parts) == 4 and parts[3] in (
                "export.json", "export.html",
            ):
                if batch["status"] == "running":
                    return self._json(409, {"detail": "Wait for this group to finish exporting."})
                for run in batch["runs"]:
                    if run["status"] == "complete":
                        run["result"] = self.service.store.snapshot(run["id"])
                if parts[3] == "export.json":
                    return self._json(200, batch)
                first = next((run["result"] for run in batch["runs"]
                              if run["status"] == "complete"), None)
                return 200, "text/html; charset=utf-8", frozen_page({"batch": batch}, first)
        if method == "GET" and len(parts) >= 3 and parts[:2] == ["api", "runs"]:
            run_id = parts[2]
            if len(parts) == 3:
                status = self.service.store.status(run_id)
                if status["status"] == "complete":
                    snapshot = self.service.store.snapshot(run_id)
                    # Large verdict arrays load only when a candidate is selected.
                    snapshot["candidates"] = [{k: v for k, v in c.items() if k != "rules"}
                                              for c in snapshot["candidates"]]
                    status["result"] = snapshot
                return self._json(200, status)
            snapshot = self.service.store.snapshot(run_id)
            if len(parts) == 5 and parts[3] == "candidates":
                candidate = next((c for c in snapshot["candidates"] if c["id"] == parts[4]), None)
                if candidate is None:
                    raise KeyError("Candidate not in this run.")
                return self._json(200, candidate)
            if len(parts) == 4 and parts[3] == "export.json":
                return self._json(200, snapshot)
            if len(parts) == 4 and parts[3] == "export.html":
                return 200, "text/html; charset=utf-8", frozen_page({"result": snapshot}, snapshot)
            if len(parts) == 4 and parts[3] == "export.csv":
                selected = _selected_for_export(snapshot, query)
                output = io.StringIO()
                writer = csv.writer(output)
                writer.writerow(["run_id", "candidate", "rank", "target_dna_5to3",
                                 "sense_rna_core_5to3",
                                 "antisense_rna_core_5to3", "predicted_guide", "guide_model",
                                 "preferences_met", "preferences_evaluated", "core_length_nt"])
                for candidate in selected:
                    rating = candidate["rating"]
                    writer.writerow([run_id, candidate["id"], candidate["rank"],
                                     candidate["target_site_dna"], candidate["sense_strand_rna"],
                                     candidate["antisense_strand_rna"], rating["guide"],
                                     rating["guide_model"], rating["met"], rating["evaluated"],
                                     len(candidate["target_site_dna"])])
                return 200, "text/csv; charset=utf-8", output.getvalue().encode()
            if len(parts) == 4 and parts[3] == "export.order.csv":
                return self._order_csv(run_id, snapshot, query)
            if len(parts) == 4 and parts[3] == "screenings":
                return self._json(200, self.service.store.screenings_for(run_id))
        return self._json(404, {"detail": "No such local app route."})


    def _start_screening(self, run_id: str, body: bytes) -> tuple[int, str, bytes]:
        """Queue an off-target screening for a finished run.

        Opt-in and asynchronous by design (ADR 0019): the run is already complete and readable,
        and starting a search sends this run's candidate sequences to NCBI. Returns at once with
        the screening's id; the interface polls it.
        """
        if len(body) > MAX_BODY:
            raise ValueError("Request is too large.")
        requested = json.loads(body) if body.strip() else {}
        scope = requested.get("scope", TRANSCRIPT_SCOPE)
        if scope not in (TRANSCRIPT_SCOPE, GENOMIC_SCOPE):
            raise ValueError(
                f"Unknown screening scope {scope!r}. Choose {TRANSCRIPT_SCOPE!r} for transcript "
                f"off-targets or {GENOMIC_SCOPE!r} for genomic ones."
            )
        snapshot = self.service.store.snapshot(run_id)
        if not any(c["rank"] is not None for c in snapshot["candidates"]):
            raise ValueError(
                "No candidate in this run passed the gates, so there is nothing to screen."
            )
        screening_id = self.service.store.start_screening(run_id, scope)
        self.screeners.submit(self.service.screen, screening_id, run_id, scope)
        return self._json(202, {"id": screening_id, "scope": scope})

    @staticmethod
    def _order_csv(
        run_id: str, snapshot: dict[str, Any], query: dict[str, list[str]]
    ) -> tuple[int, str, bytes]:
        """Write the order-ready duplexes: 19 nt RNA cores plus dTdT, with their sanity report.

        Distinct from ``export.csv``, which is deliberately a core-sequence export and not order
        chemistry (ADR 0015). Rejected candidates are never included by default: a duplex that
        failed a gate is not something to order, and its rule verdicts remain the record of why.
        """
        selected = _selected_for_export(snapshot, query)
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(CSV_COLUMNS)
        for candidate in selected:
            # The model this run actually displayed, taken from the frozen record rather than
            # today's default, so a reopened run reports the strand prediction it was read with.
            model = (candidate.get("rating") or {}).get("guide_model") or DEFAULT_MODEL
            writer.writerow(csv_row(order_from_record(candidate, model), run_id))
        return 200, "text/csv; charset=utf-8", output.getvalue().encode()

    def _queue_batch(
        self, request: BatchRequest, reuse: dict[int, str] | None = None,
        retry_of: str | None = None,
    ) -> tuple[int, str, bytes]:
        """Validate every child's configuration before atomically saving and queuing the group."""
        for child in request.requests:
            try:
                child.resolve()
            except KeyError as error:
                raise ValueError(str(error)) from error
        batch_id, run_ids, queued_ids = self.service.store.create_batch(
            request.model_dump(mode="json"), reuse=reuse, retry_of=retry_of,
        )
        queued = set(queued_ids)
        for run_id, child in zip(run_ids, request.requests, strict=True):
            if run_id in queued:
                self.workers.submit(self.service.run, run_id, child)
        return self._json(202, {"id": batch_id, "run_ids": run_ids, "retry_of": retry_of})


def make_server(router: AppRouter, port: int) -> ThreadingHTTPServer:
    """Create a same-origin server bound only to this machine's loopback interface."""
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            """Serve a resource or read a saved run."""
            self.respond("GET")

        def do_POST(self) -> None:
            """Accept a bounded JSON draft from the local workbench."""
            self.respond("POST")

        def respond(self, method: str) -> None:
            """Enforce local host/origin and body limits before routing."""
            host = self.headers.get("Host", "")
            origin = self.headers.get("Origin")
            local_port = cast("ThreadingHTTPServer", self.server).server_port
            if host not in (f"127.0.0.1:{local_port}", f"localhost:{local_port}") or (
                origin is not None and origin != f"http://{host}"
            ):
                self.send_error(403, "Use the local application's own origin.")
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                self.send_error(400, "Invalid Content-Length")
                return
            if not 0 <= length <= MAX_BODY:
                self.send_error(413, "Request too large")
                return
            if method == "POST" and self.headers.get_content_type() != "application/json":
                self.send_error(415, "Submit a JSON request")
                return
            status, media, payload = router.dispatch(method, self.path, self.rfile.read(length))
            self.send_response(status)
            self.send_header("Content-Type", media)
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(payload)

    return ThreadingHTTPServer(("127.0.0.1", port), Handler)
