"""Loopback-only HTTP adapter; the Python engine owns every scientific result."""

import csv
import io
import json
import logging
import re
from base64 import b64encode
from collections.abc import Mapping, Sequence
from collections.abc import Set as AbstractSet
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from typing import Any, cast
from urllib.parse import parse_qs, unquote, urlsplit

from pydantic import ValidationError

from sarna import __version__, settings
from sarna.app.models import BatchRequest, DesignRequest, ShortlistRequest
from sarna.app.service import DesignService
from sarna.app.workbook import Sheet, write_workbook
from sarna.design.duplicates import identical_duplex_groups
from sarna.design.order import (
    CSV_COLUMNS,
    DEFAULT_MODEL,
    ScreeningFinding,
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
from sarna.offtarget.ncbi import candidate_sequences
from sarna.offtarget.screening import DEFAULT_CANDIDATES, MAX_CANDIDATES
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
    offtarget = STATIC.joinpath("offtarget.js").read_text(encoding="utf-8")
    workflow = STATIC.joinpath("workflow.js").read_text(encoding="utf-8")
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
    # Order matters: the off-target panel is a second classic script sharing app.js's globals,
    # and a frozen page that dropped it would lose the whole screening card without saying so.
    page = page.replace('<script src="/offtarget.js" defer></script>',
                        f"<script>{offtarget}</script>")
    page = page.replace('<script src="/workflow.js" defer></script>',
                        f"<script>{workflow}</script>")
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
        if method == "GET" and path in (
            "/", "/app.js", "/offtarget.js", "/workflow.js", "/app.css",
            "/favicon.svg", "/favicon.ico",
        ):
            name = {"/": "index.html", "/favicon.ico": "favicon.svg"}.get(path, path[1:])
            media = {"html": "text/html", "js": "text/javascript", "css": "text/css",
                     "svg": "image/svg+xml"}[name.rsplit(".", 1)[1]]
            return 200, media + "; charset=utf-8", STATIC.joinpath(name).read_bytes()
        if method == "GET" and path == "/api/catalogue":
            return self._json(200, catalogue())
        if method == "GET" and path.startswith("/api/genes/"):
            rest = path.removeprefix("/api/genes/")
            # Expression is a second, optional request against a second service, so it gets
            # its own route rather than making every gene lookup wait on GTEx.
            symbol, _, leaf = rest.partition("/")
            if leaf not in ("", "expression"):
                raise KeyError(f"No such gene route: {leaf!r}")
            if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", symbol):
                raise ValueError("Use a human gene symbol or Ensembl gene ID, without spaces.")
            offline = query.get("offline") == ["true"]
            if leaf == "expression":
                return self._json(200, self.service.expression(symbol, offline))
            return self._json(200, self.service.lookup(symbol, offline))
        if path == "/api/batches":
            if method == "GET":
                return self._json(200, self.service.store.list_batches())
            if method == "POST":
                if len(body) > MAX_BODY:
                    raise ValueError("Request is too large.")
                return self._queue_batch(BatchRequest.model_validate_json(body))
        if path == "/api/runs":
            if method == "GET":
                return self._json(200, self.service.store.list_runs(
                    search=query.get("q", [""])[0], status=query.get("status", [""])[0],
                    offset=int(query.get("offset", ["0"])[0]),
                    limit=int(query.get("limit", ["100"])[0]),
                ))
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
        if (method == "POST" and len(parts) == 4
                and parts[:2] == ["api", "runs"] and parts[3] == "shortlist"):
            if len(body) > MAX_BODY:
                raise ValueError("Request is too large.")
            draft = ShortlistRequest.model_validate_json(body)
            snapshot = self.service.store.snapshot(parts[2])
            known = {candidate["id"] for candidate in snapshot["candidates"]}
            if any(entry.candidate not in known for entry in draft.entries):
                raise ValueError("A shortlisted candidate does not belong to this run.")
            saved = self.service.store.save_shortlist(
                parts[2], draft.revision, [entry.model_dump() for entry in draft.entries],
            )
            return (self._json(200, saved) if saved is not None else self._json(409, {
                "detail": "This shortlist changed in another window. Reload saved notes "
                          "before replacing them.", "code": "shortlist-conflict",
            }))
        if method == "GET" and len(parts) == 3 and parts[:2] == ["api", "screenings"]:
            return self._json(200, self.service.store.screening(parts[2]))
        if (method == "POST" and len(parts) == 4
                and parts[:2] == ["api", "screenings"] and parts[3] == "annotations"):
            if len(body) > MAX_BODY:
                raise ValueError("Request is too large.")
            requested = json.loads(body)
            accessions = requested.get("accessions")
            if (not isinstance(accessions, list) or not 1 <= len(accessions) <= 25
                    or not all(isinstance(a, str) for a in accessions)):
                raise ValueError("Request between 1 and 25 hit accessions.")
            return self._json(200, self.service.screening_annotations(parts[2], accessions))
        if (method == "POST" and len(parts) == 4
                and parts[:2] == ["api", "screenings"] and parts[3] == "resume"):
            return self._resume_screening(parts[2])
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
                        run["result"] = self._workspace_snapshot(run["id"])
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
                    snapshot = self._workspace_snapshot(run_id)
                    # Large verdict arrays load only when a candidate is selected.
                    snapshot["candidates"] = [{k: v for k, v in c.items() if k != "rules"}
                                              for c in snapshot["candidates"]]
                    status["result"] = snapshot
                return self._json(200, status)
            snapshot = self.service.store.snapshot(run_id)
            if len(parts) == 4 and parts[3] == "shortlist" and method == "GET":
                return self._json(200, self.service.store.shortlist(run_id))
            if len(parts) == 5 and parts[3] == "candidates":
                candidate = next((c for c in snapshot["candidates"] if c["id"] == parts[4]), None)
                if candidate is None:
                    raise KeyError("Candidate not in this run.")
                return self._json(200, candidate)
            if len(parts) == 4 and parts[3] == "export.json":
                return self._json(200, self._workspace_snapshot(run_id))
            if len(parts) == 4 and parts[3] == "export.html":
                # The screenings travel with the archive. They are separate dated records and
                # stay separate, but an offline copy of a screened run that showed no screening
                # would be a record quietly missing its own evidence.
                return 200, "text/html; charset=utf-8", frozen_page(
                    {"result": self._workspace_snapshot(run_id),
                     "screenings": self._complete_screenings(run_id)},
                    snapshot,
                )
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
            if len(parts) == 4 and parts[3] == "export.xlsx":
                return self._export_xlsx(run_id, snapshot, query)
            if len(parts) == 4 and parts[3] == "screenings":
                return self._json(200, self.service.store.screenings_for(run_id))
        return self._json(404, {"detail": "No such local app route."})

    def _workspace_snapshot(self, run_id: str) -> dict[str, Any]:
        """Attach current researcher notes to exports without changing the scientific snapshot."""
        return {**self.service.store.snapshot(run_id),
                "workspace": {"shortlist": self.service.store.shortlist(run_id)}}


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
        passing = {c["id"] for c in snapshot["candidates"] if c["rank"] is not None}
        if not passing:
            raise ValueError(
                "No candidate in this run passed the gates, so there is nothing to screen."
            )
        only = self._requested_candidates(requested, passing)
        limit = self._requested_limit(requested)
        chosen = list(candidate_sequences(snapshot["candidates"], limit=limit, only=only))
        screening_id, reused = self.service.store.reserve_screening(run_id, scope, chosen)
        if not reused:
            self.screeners.submit(
                self.service.screen, screening_id, run_id, scope, only, limit
            )
        return self._json(202, {
            "id": screening_id,
            "reused": reused,
            "scope": scope,
            "candidates": len(only) if only is not None else min(limit, len(passing)),
        })

    def _resume_screening(self, screening_id: str) -> tuple[int, str, bytes]:
        """Collect a search that outran its wait, rather than running it again.

        NCBI holds a request id for about a day. Resubmitting would queue a second identical
        job against a shared service and throw away the first (ADR 0024).
        """
        held = self.service.store.screening(screening_id)
        if held["status"] == "complete":
            return self._json(409, {"detail": "This screening has already finished."})
        if held["status"] == "running":
            return self._json(409, {"detail": "This screening is already being collected."})
        if held["status"] != "waiting" or not held.get("request_id"):
            raise ValueError(
                "This screening has no NCBI request id to collect, so it has to be run again "
                "rather than resumed."
            )
        # Back to 'running' before the worker starts, or the interface offers to collect the
        # same search again while it is already being collected.
        if not self.service.store.resume_screening(screening_id):
            return self._json(409, {"detail": "This screening is already being collected."})
        self.screeners.submit(self.service.resume, screening_id)
        return self._json(202, {"id": screening_id, "request_id": held["request_id"]})

    @staticmethod
    def _requested_candidates(
        requested: Mapping[str, Any], passing: AbstractSet[str]
    ) -> list[str] | None:
        """Return the explicit candidate selection, validated, or ``None`` for "the best N".

        Every batch the interface offers resolves to a list of identifiers (ADR 0023), so this
        is the one place the selection is checked. An unknown or excluded identifier is named
        rather than dropped: silently screening fewer candidates than were asked for would be
        indistinguishable from screening them and finding nothing.
        """
        chosen = requested.get("candidates")
        if chosen is None:
            return None
        if not isinstance(chosen, list) or not all(isinstance(c, str) for c in chosen):
            raise ValueError("'candidates' must be a list of candidate identifiers.")
        if not chosen:
            raise ValueError("No candidates were chosen, so there is nothing to screen.")
        if len(chosen) > MAX_CANDIDATES:
            raise ValueError(
                f"{len(chosen)} candidates were chosen; one search covers at most "
                f"{MAX_CANDIDATES}. Screening runs in batches, so start the next one when this "
                f"finishes."
            )
        unknown = [c for c in chosen if c not in passing]
        if unknown:
            raise ValueError(
                f"These are not candidates of this run that passed its gates: "
                f"{', '.join(sorted(unknown)[:5])}."
            )
        return list(dict.fromkeys(chosen))

    @staticmethod
    def _requested_limit(requested: Mapping[str, Any]) -> int:
        """Return how many to screen when no explicit selection was made."""
        count = requested.get("count", DEFAULT_CANDIDATES)
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError("'count' must be a whole number of candidates, at least 1.")
        return min(count, MAX_CANDIDATES)

    def _complete_screenings(self, run_id: str) -> list[dict[str, Any]]:
        """Return every finished screening of one run, oldest first, with its full result.

        Oldest first so a reader merging them sees later searches win, which is the order the
        interface applies. A running or failed screening is left out: a search that did not
        finish is not a search that found nothing.
        """
        listed = self.service.store.screenings_for(run_id)
        return [
            self.service.store.screening(entry["id"])
            for entry in reversed(listed)
            if entry["status"] == "complete"
        ]

    def _screenings_by_candidate(
        self, run_id: str
    ) -> dict[str, list[ScreeningFinding]]:
        """Return the newest complete finding per candidate and scope.

        ADR 0019 left this join open, so the order export reported "not-performed" even after
        a real screening, honest, and useless to somebody deciding what to synthesise. A
        screening stays its own dated record; this reads it rather than writing into the run.

        Newest per candidate and scope, because batches cover different candidates and the
        export should carry the most recent observation. An incomplete or failed screening
        contributes nothing: a search that did not finish is not a search that found nothing.
        """
        findings: dict[str, list[ScreeningFinding]] = {}
        seen: set[tuple[str, str]] = set()
        for listed in self.service.store.screenings_for(run_id):
            if listed["status"] != "complete":
                continue
            result = (self.service.store.screening(listed["id"]) or {}).get("result") or {}
            parameters = result.get("parameters") or {}
            for row in result.get("specificity") or []:
                key = (row["candidate"], listed["scope"])
                if key in seen:
                    continue
                seen.add(key)
                findings.setdefault(row["candidate"], []).append(ScreeningFinding(
                    scope=listed["scope"],
                    database=str(parameters.get("database", "unknown")),
                    searched_at=str(result.get("searched_at", "")),
                    row=row,
                ))
        return findings

    CANDIDATE_COLUMNS = (
        "rank", "candidate", "tss_start", "tss_end", "contig", "genomic_start", "genomic_end",
        "strand", "target_dna_5to3", "sense_core_rna_5to3", "antisense_core_rna_5to3",
        "sense_order_5to3", "antisense_order_5to3", "gc_percent", "core_length_nt",
        "predicted_guide", "guide_model", "models_agree", "preferences_met",
        "preferences_evaluated", "preferences_missed", "unchecked_gates", "excluded_by",
        "off_target_perfect", "off_target_wang80", "off_target_near_perfect",
        "off_target_alignments", "off_target_list_capped", "off_target_scopes",
    )
    """Everything the interface knows about a candidate, in one row.

    Wider than ``export.csv`` on purpose: that one is deliberately a short core-sequence export
    (ADR 0015), and this is the "give me all of them with their values" file.
    """

    @staticmethod
    def _candidate_row(
        candidate: dict[str, Any], findings: Sequence[ScreeningFinding],
        rule_names: Mapping[str, str],
    ) -> tuple[Any, ...]:
        """Return one candidate as a spreadsheet row.

        ``rule_names`` is passed rather than read from anywhere shared: this server is
        threaded, and two exports of different runs can be in flight at once.
        """
        def named(ids: Sequence[str]) -> str:
            """Rule ids are stable keys; a reader wants the names."""
            return "; ".join(rule_names.get(identifier, identifier) for identifier in ids)
        rating = candidate.get("rating") or {}
        genomic = candidate.get("genomic") or {}
        span = candidate.get("tss_span") or (None, None)
        order = order_from_record(
            candidate, rating.get("guide_model") or DEFAULT_MODEL, screenings=tuple(findings)
        ) if candidate.get("rank") is not None else None
        total = {key: sum(int(f.row.get(key) or 0) for f in findings)
                 for key in ("perfect_hits", "wang_homology_hits", "near_perfect_hits",
                             "hit_count")}
        # A capped list makes every count above it a floor, which has to travel with the row
        # rather than stay in the workbench (ADR 0023).
        capped = any(f.row.get("alignments_capped") for f in findings)
        return (
            candidate.get("rank"), candidate["id"], span[0], span[1],
            genomic.get("contig"), genomic.get("start"), genomic.get("end"),
            genomic.get("strand"),
            candidate["target_site_dna"], candidate["sense_strand_rna"],
            candidate["antisense_strand_rna"],
            order.sense.annotated if order else "", order.antisense.annotated if order else "",
            candidate.get("gc_percent"), len(candidate["target_site_dna"]),
            rating.get("guide"), rating.get("guide_model"), rating.get("models_agree"),
            rating.get("met"), rating.get("evaluated"),
            named(rating.get("missed") or []),
            named(rating.get("unchecked_gates") or []),
            named(candidate.get("excluded_by") or []),
            total["perfect_hits"] if findings else "",
            total["wang_homology_hits"] if findings else "",
            total["near_perfect_hits"] if findings else "",
            total["hit_count"] if findings else "",
            ("yes" if capped else "no") if findings else "",
            "; ".join(f.scope for f in findings) or "not screened",
        )

    def _export_xlsx(
        self, run_id: str, snapshot: dict[str, Any], query: dict[str, list[str]]
    ) -> tuple[int, str, bytes]:
        """Write the whole run as a spreadsheet: candidates, the run, and its rules.

        An ``.xlsx`` rather than a third CSV because Excel resolves a CSV's separator against
        the machine's locale: on a Dutch or German Windows a comma-separated file opens with
        every column in column A. See :mod:`sarna.app.workbook`.
        """
        selected = _selected_for_export(snapshot, query)
        screenings = self._screenings_by_candidate(run_id)
        rules = snapshot.get("configuration", {}).get("rules", [])
        rule_names = {rule["id"]: rule["name"] for rule in rules}
        candidates = Sheet(
            "saRNA candidates", self.CANDIDATE_COLUMNS,
            tuple(self._candidate_row(c, screenings.get(c["id"], ()), rule_names)
                  for c in selected),
        )
        performance = snapshot.get("ranking_performance") or {}
        reproducibility = snapshot.get("reproducibility") or {}
        software = reproducibility.get("software") or {}
        app = snapshot.get("app", {})
        gene = snapshot.get("gene") or {}
        anchor = app.get("anchor") or {}
        no_reference = "none (designed against a supplied sequence)"
        run = Sheet("Run", ("field", "value"), (
            ("run id", run_id),
            ("created", app.get("created", "")),
            ("gene", gene.get("symbol", "pasted sequence")),
            ("Ensembl gene id", gene.get("ensembl_id", "")),
            # A pasted sequence genuinely has no annotated reference. Saying so beats an
            # empty cell, which reads as a field somebody forgot to fill in.
            ("promoter reference", anchor.get("label", "") or no_reference),
            ("reference kind", anchor.get("reference_kind", "") or no_reference),
            ("reference provenance", anchor.get("provenance", "") or no_reference),
            ("reference caveat", anchor.get("caveat") or "none recorded"),
            ("promoter span (bp)", snapshot.get("promoter", {}).get("length")),
            ("cell line", app.get("cell_line", "")),
            ("configuration", snapshot.get("configuration", {}).get("preset", {}).get("name")),
            ("customised", snapshot.get("configuration", {}).get("is_custom")),
            ("targets evaluated", snapshot.get("counts", {}).get("generated")),
            ("passing every evaluable gate", snapshot.get("counts", {}).get("passing")),
            ("engine version", app.get("engine_version", "")),
            ("Git revision", software.get("git_revision") or "not recorded"),
            ("uncommitted changes", software.get("git_dirty")
             if software.get("git_dirty") is not None else "not recorded"),
            ("package SHA-256", software.get("package_sha256") or "not recorded"),
            ("promoter sequence SHA-256", reproducibility.get("sequence_sha256") or "not recorded"),
            ("resolved configuration SHA-256", reproducibility.get("configuration_sha256")
             or "not recorded"),
            ("replay scope", reproducibility.get("replay_scope") or "not recorded"),
            ("ranking performance", performance.get("headline", "not measured for this run")),
            ("what a rating is", "criteria met, counted equally. Not a probability that the "
                                 "duplex activates the gene."),
            ("off-target screening", "; ".join(sorted({
                f.scope for rows in screenings.values() for f in rows})) or "none performed"),
        ))
        rule_sheet = Sheet("Rules", ("rule", "name", "role", "evidence", "source",
                                     "parameters", "source URL", "evidence scope",
                                     "application note"), tuple(
            (rule["id"], rule["name"], rule["kind"], rule.get("evidence_class", ""),
             rule.get("provenance", {}).get("display", ""),
             "; ".join(f"{p['name']}={p['value']}" for p in rule.get("parameters", [])),
             rule.get("provenance", {}).get("url") or "not recorded",
             rule.get("evidence_note", "Evidence scope not recorded in this saved run"),
             rule.get("application_note", "Application note not recorded in this saved run"))
            for rule in rules
        ))
        source_sheet = Sheet("Reference sources", (
            "source", "URL", "response SHA-256", "cache hit", "observed at",
            "original retrieval", "annotation release", "note",
        ), tuple(
            (source["source"], source["url"], source["response_sha256"], source["cache_hit"],
             source["observed_at"], source["retrieved_at"] or "unknown",
             source["annotation_release"] or "unknown", source["note"])
            for source in reproducibility.get("sources", [])
        ))
        shortlist = self.service.store.shortlist(run_id)
        by_id = {candidate["id"]: candidate for candidate in snapshot["candidates"]}
        shortlist_sheet = Sheet("Shortlist", (
            "candidate", "researcher reason", "checklist rank", "excluded by gates",
            "saved revision", "updated at",
        ), tuple(
            (entry["candidate"], entry["reason"], by_id[entry["candidate"]]["rank"],
             "; ".join(by_id[entry["candidate"]]["excluded_by"]),
             shortlist["revision"], shortlist["updated"])
            for entry in shortlist["entries"]
        ))
        data = write_workbook([candidates, run, rule_sheet, source_sheet, shortlist_sheet])
        return 200, ("application/vnd.openxmlformats-officedocument."
                     "spreadsheetml.sheet"), data

    def _order_csv(
        self, run_id: str, snapshot: dict[str, Any], query: dict[str, list[str]]
    ) -> tuple[int, str, bytes]:
        """Write the order-ready duplexes: 19 nt RNA cores plus dTdT, with their sanity report.

        Distinct from ``export.csv``, which is deliberately a core-sequence export and not order
        chemistry (ADR 0015). Rejected candidates are never included by default: a duplex that
        failed a gate is not something to order, and its rule verdicts remain the record of why.
        """
        selected = _selected_for_export(snapshot, query)
        screenings = self._screenings_by_candidate(run_id)
        # Twins are looked for across the whole run, not the selection, so ticking one of eight
        # identical rows still says so. Which of them are in *this* file is the separate
        # question, and the only one that means a product is about to be bought twice.
        twins = identical_duplex_groups(
            (c["id"], c["target_site_dna"]) for c in snapshot["candidates"]
        )
        exported = {c["id"] for c in selected}
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(CSV_COLUMNS)
        for candidate in selected:
            # The model this run actually displayed, taken from the frozen record rather than
            # today's default, so a reopened run reports the strand prediction it was read with.
            model = (candidate.get("rating") or {}).get("guide_model") or DEFAULT_MODEL
            in_run = twins.get(candidate["id"], ())
            writer.writerow(csv_row(
                order_from_record(
                    candidate, model,
                    screenings=tuple(screenings.get(candidate["id"], ())),
                    identical_in_run=in_run,
                    identical_in_export=tuple(i for i in in_run if i in exported),
                ),
                run_id))
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
