"""SQLite run snapshots; reopening never consults today's presets or remote sources."""

import json
import sqlite3
import zlib
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from uuid import uuid4


def run_summary(snapshot: dict[str, Any]) -> dict[str, Any]:
    """Keep authoritative per-reference metrics without pooling or re-ranking candidates."""
    first = next((c for c in snapshot["candidates"] if c["rank"] == 1), None)
    return {
        "anchor": snapshot["app"]["anchor"], "counts": snapshot["counts"],
        "best_rating": first["rating"] if first else None,
        "configuration": {"preset_name": snapshot["configuration"]["preset"]["name"],
                          "is_custom": snapshot["configuration"]["is_custom"]},
        "promoter": {key: snapshot["promoter"][key] for key in (
            "contig", "start", "end", "strand", "length"
        )},
        "cell_line": snapshot["app"]["cell_line"],
    }


def batch_status(runs: list[dict[str, Any]]) -> str:
    """Separate unfinished, wholly complete, partly failed and wholly failed groups."""
    states = {run["status"] for run in runs}
    if "running" in states:
        return "running"
    if states == {"complete"}:
        return "complete"
    return "partial" if "complete" in states else "error"


class RunStore:
    """Keep request, status and compressed complete records in one local database."""

    def __init__(self, path: Path) -> None:
        """Create the versioned local schema if it is absent."""
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        with self.connect() as connection:
            connection.execute("PRAGMA journal_mode=WAL")
            connection.execute(
                "CREATE TABLE IF NOT EXISTS runs ("
                "id TEXT PRIMARY KEY, created TEXT NOT NULL, label TEXT NOT NULL, "
                "status TEXT NOT NULL, stage TEXT NOT NULL, request TEXT NOT NULL, "
                "snapshot BLOB, error TEXT, error_code TEXT, summary TEXT)"
            )
            columns = {row["name"] for row in connection.execute("PRAGMA table_info(runs)")}
            if "error_code" not in columns:
                connection.execute("ALTER TABLE runs ADD COLUMN error_code TEXT")
            if "summary" not in columns:
                connection.execute("ALTER TABLE runs ADD COLUMN summary TEXT")
            connection.execute(
                "CREATE TABLE IF NOT EXISTS batches (id TEXT PRIMARY KEY, "
                "created TEXT NOT NULL, label TEXT NOT NULL, symbol TEXT NOT NULL, "
                "request TEXT NOT NULL, retry_of TEXT)"
            )
            connection.execute(
                "CREATE TABLE IF NOT EXISTS batch_runs (batch_id TEXT NOT NULL, "
                "run_id TEXT NOT NULL, ordinal INTEGER NOT NULL, "
                "PRIMARY KEY(batch_id,ordinal), "
                "FOREIGN KEY(batch_id) REFERENCES batches(id), "
                "FOREIGN KEY(run_id) REFERENCES runs(id))"
            )
            # Off-target screenings are separate, dated records rather than columns on a run:
            # NCBI's databases change week to week, so a result is only meaningful with its own
            # date, database and parameters, and re-screening must add a record rather than
            # overwrite one (ADR 0019). The design snapshot stays frozen.
            connection.execute(
                "CREATE TABLE IF NOT EXISTS screenings ("
                "id TEXT PRIMARY KEY, run_id TEXT NOT NULL, scope TEXT NOT NULL, "
                "created TEXT NOT NULL, status TEXT NOT NULL, stage TEXT NOT NULL, "
                "request_id TEXT, result BLOB, error TEXT, "
                "FOREIGN KEY(run_id) REFERENCES runs(id))"
            )
            screening_columns = {
                row["name"] for row in connection.execute("PRAGMA table_info(screenings)")
            }
            if "progress" not in screening_columns:
                # Added for the progress bar. A stage sentence says what is happening; this
                # says how far along it is, which a sentence cannot. Nullable, so screenings
                # saved before it exists stay readable.
                connection.execute("ALTER TABLE screenings ADD COLUMN progress TEXT")
            connection.execute(
                "CREATE INDEX IF NOT EXISTS screenings_by_run ON screenings(run_id, created)"
            )
            connection.execute(
                "CREATE TABLE IF NOT EXISTS shortlists (run_id TEXT PRIMARY KEY, "
                "revision INTEGER NOT NULL, updated TEXT, entries TEXT NOT NULL, "
                "FOREIGN KEY(run_id) REFERENCES runs(id))"
            )
            connection.execute("PRAGMA user_version=6")
            connection.execute(
                "UPDATE screenings SET status='waiting', stage='Report retrieval needs retry' "
                "WHERE status='error' AND request_id IS NOT NULL AND result IS NOT NULL "
                "AND error LIKE 'could not parse NCBI''s XML reply:%'"
            )
            connection.execute(
                "UPDATE screenings SET status='waiting', stage='Ready to check NCBI again', "
                "error='The app restarted. The submitted search can still be collected.' "
                "WHERE status='running' AND request_id IS NOT NULL AND result IS NOT NULL"
            )
            connection.execute(
                # Only 'running' is interrupted by a restart. A 'waiting' screening is not
                # running here at all: it is queued at NCBI under a request id that outlives
                # this process, so it stays collectable (ADR 0024).
                "UPDATE screenings SET status='error', stage='Interrupted', "
                "error='The application stopped while this screening was running. Start it "
                "again.' WHERE status='running'"
            )
            connection.execute(
                "UPDATE runs SET status='error', error='Application stopped before this run "
                "completed. Submit it again.', stage='Interrupted' WHERE status='running'"
            )

    @contextmanager
    def connect(self) -> Iterator[sqlite3.Connection]:
        """Open one connection per operation so background jobs never share connections."""
        connection = sqlite3.connect(self.path, timeout=30)
        connection.execute("PRAGMA foreign_keys=ON")
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def create(self, request: dict[str, Any]) -> str:
        """Store the submitted draft and return a unique run id."""
        run_id = uuid4().hex
        with self.connect() as connection:
            connection.execute(
                "INSERT INTO runs(id,created,label,status,stage,request) VALUES(?,?,?,?,?,?)",
                (run_id, datetime.now(UTC).isoformat(), self.request_label(request),
                 "running", "Queued", json.dumps(request)),
            )
        return run_id

    @staticmethod
    def request_label(request: dict[str, Any]) -> str:
        """Identify reference-specific runs in history instead of repeating only a gene name."""
        label = request.get("symbol") or "Pasted DNA"
        custom = request.get("custom_reference")
        reference = (custom["name"] if custom else
                     request.get("anchor_id") or request.get("transcript_id"))
        return f"{label} · {reference}" if reference else str(label)

    def create_batch(
        self, request: dict[str, Any], reuse: dict[int, str] | None = None,
        retry_of: str | None = None,
    ) -> tuple[str, list[str], list[str]]:
        """Create a group and child jobs atomically, optionally reusing completed records."""
        batch_id, created = uuid4().hex, datetime.now(UTC).isoformat()
        requests = request["requests"]
        symbol = requests[0]["symbol"]
        label = request.get("label") or f"{symbol} · {len(requests)} promoter references"
        run_ids, queued_ids = [], []
        with self.connect() as connection:
            connection.execute("INSERT INTO batches VALUES(?,?,?,?,?,?)", (
                batch_id, created, label, symbol, json.dumps(request), retry_of,
            ))
            for ordinal, child in enumerate(requests):
                run_id = (reuse or {}).get(ordinal)
                if run_id is not None:
                    existing = connection.execute(
                        "SELECT status,request FROM runs WHERE id=?", (run_id,),
                    ).fetchone()
                    if existing is None or existing["status"] != "complete" or (
                        json.loads(existing["request"]) != child
                    ):
                        raise ValueError("Only an unchanged completed run can be reused.")
                else:
                    run_id = uuid4().hex
                    connection.execute(
                        "INSERT INTO runs(id,created,label,status,stage,request) "
                        "VALUES(?,?,?,?,?,?)",
                        (run_id, created, self.request_label(child), "running", "Queued",
                         json.dumps(child)),
                    )
                    queued_ids.append(run_id)
                connection.execute("INSERT INTO batch_runs VALUES(?,?,?)", (
                    batch_id, run_id, ordinal,
                ))
                run_ids.append(run_id)
        return batch_id, run_ids, queued_ids

    def batch(self, batch_id: str) -> dict[str, Any]:
        """Read group progress and saved summaries without opening full candidate arrays."""
        with self.connect() as connection:
            row = connection.execute("SELECT * FROM batches WHERE id=?", (batch_id,)).fetchone()
            if row is None:
                raise KeyError(batch_id)
            runs = [self.decode_status(child) for child in connection.execute(
                "SELECT r.id,r.created,r.label,r.status,r.stage,r.error,r.error_code,"
                "r.request,r.summary FROM runs r JOIN batch_runs b ON b.run_id=r.id "
                "WHERE b.batch_id=? ORDER BY b.ordinal", (batch_id,),
            )]
        result = dict(row)
        result["request"] = json.loads(result["request"])
        result["runs"] = runs
        result["status"] = batch_status(runs)
        return result

    def list_batches(self) -> list[dict[str, Any]]:
        """List recent saved groups, including partial failures and retry lineage."""
        with self.connect() as connection:
            rows = list(connection.execute(
                "SELECT id FROM batches ORDER BY created DESC,rowid DESC LIMIT 100"
            ))
        return [{key: batch[key] for key in (
            "id", "created", "label", "symbol", "status", "retry_of"
        )} for row in rows for batch in (self.batch(row["id"]),)]

    def progress(self, run_id: str, stage: str) -> None:
        """Publish a stage while calculations continue in a worker thread."""
        with self.connect() as connection:
            connection.execute("UPDATE runs SET stage=? WHERE id=?", (stage, run_id))

    def complete(self, run_id: str, snapshot: dict[str, Any]) -> None:
        """Freeze a complete scientific record, compressed to avoid repeated prose on disk."""
        blob = zlib.compress(json.dumps(snapshot, separators=(",", ":")).encode(), level=6)
        with self.connect() as connection:
            connection.execute(
                "UPDATE runs SET status='complete', stage='Complete', snapshot=?,summary=? "
                "WHERE id=?",
                (blob, json.dumps(run_summary(snapshot)), run_id),
            )

    def fail(self, run_id: str, reason: str, error_code: str | None = None) -> None:
        """Retain a failed request with its explanation instead of a misleading empty run."""
        with self.connect() as connection:
            connection.execute(
                "UPDATE runs SET status='error',stage='Failed',error=?,error_code=? WHERE id=?",
                (reason, error_code, run_id),
            )

    def status(self, run_id: str) -> dict[str, Any]:
        """Return progress without loading the snapshot."""
        with self.connect() as connection:
            row = connection.execute(
                "SELECT id,created,label,status,stage,error,error_code,request,summary "
                "FROM runs WHERE id=?",
                (run_id,),
            ).fetchone()
        if row is None:
            raise KeyError(run_id)
        return self.decode_status(row)

    @staticmethod
    def decode_status(row: sqlite3.Row) -> dict[str, Any]:
        """Decode stored request and summary metadata while leaving snapshots compressed."""
        result = dict(row)
        result["request"] = json.loads(result["request"])
        result["summary"] = json.loads(result["summary"]) if result.get("summary") else None
        return result

    def snapshot(self, run_id: str) -> dict[str, Any]:
        """Read the saved record without rerunning any scientific code."""
        with self.connect() as connection:
            row = connection.execute("SELECT snapshot FROM runs WHERE id=?", (run_id,)).fetchone()
        if row is None or row["snapshot"] is None:
            raise KeyError(run_id)
        result: dict[str, Any] = json.loads(zlib.decompress(row["snapshot"]))
        return result


    def reserve_screening(
        self, run_id: str, scope: str, candidates: list[str], backend: str = "ncbi",
        reference_release: str | None = None,
    ) -> tuple[str, bool]:
        """Atomically reuse an identical active search, otherwise reserve a new queue slot."""
        self.status(run_id)
        wanted = sorted(candidates)
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            rows = connection.execute(
                "SELECT id,result FROM screenings WHERE run_id=? AND scope=? "
                "AND status IN ('running','waiting')", (run_id, scope),
            ).fetchall()
            for row in rows:
                pending = json.loads(zlib.decompress(row["result"])) if row["result"] else {}
                if (sorted(pending.get("candidates_searched", [])) == wanted
                        and pending.get("backend", "ncbi") == backend
                        and pending.get("reference_release") == reference_release):
                    return row["id"], True
            screening_id = uuid4().hex
            pending = zlib.compress(json.dumps(
                {"candidates_searched": candidates, "backend": backend,
                 "reference_release": reference_release}
            ).encode())
            connection.execute(
                "INSERT INTO screenings(id,run_id,scope,created,status,stage,result) "
                "VALUES(?,?,?,?,?,?,?)",
                (screening_id, run_id, scope, datetime.now(UTC).isoformat(),
                 "running", "Queued locally", pending),
            )
        return screening_id, False

    def start_screening(self, run_id: str, scope: str) -> str:
        """Record that a screening has begun and return its id.

        Raises:
            KeyError: If the run does not exist. A screening with no run would be unattributable.
        """
        self.status(run_id)
        screening_id = uuid4().hex
        with self.connect() as connection:
            connection.execute(
                "INSERT INTO screenings(id,run_id,scope,created,status,stage) "
                "VALUES(?,?,?,?,?,?)",
                (screening_id, run_id, scope, datetime.now(UTC).isoformat(),
                 "running", "Queued locally"),
            )
        return screening_id

    def screening_progress(
        self,
        screening_id: str,
        stage: str,
        request_id: str | None = None,
        progress: dict[str, Any] | None = None,
    ) -> None:
        """Record how far a screening has got, so the interface can say rather than spin.

        ``progress`` carries the countable parts: which phase, how long NCBI said it would
        take, how long it has actually taken, and how many hits have been placed. A bar needs
        numbers; ``stage`` is a sentence and cannot supply them.
        """
        fields = ["stage=?"]
        values: list[Any] = [stage]
        request_id = request_id or (progress or {}).get("request_id")
        if request_id is not None:
            fields.append("request_id=?")
            values.append(request_id)
        if progress is not None:
            fields.append("progress=?")
            values.append(json.dumps(progress, separators=(",", ":")))
        if progress and progress.get("candidates_searched"):
            fields.append("result=?")
            values.append(zlib.compress(json.dumps({
                "candidates_searched": progress["candidates_searched"],
                "parameters": progress.get("parameters"),
                "request_id": request_id,
                "submitted_at": datetime.now(UTC).isoformat(),
            }).encode()))
        values.append(screening_id)
        with self.connect() as connection:
            connection.execute(
                f"UPDATE screenings SET {', '.join(fields)} WHERE id=?", tuple(values))

    def finish_screening(self, screening_id: str, result: dict[str, Any]) -> None:
        """Store a completed screening."""
        with self.connect() as connection:
            connection.execute(
                "UPDATE screenings SET status='complete', stage='Complete', error=NULL, "
                "request_id=?, result=? WHERE id=?",
                (result.get("request_id"),
                 zlib.compress(json.dumps(result, separators=(",", ":")).encode(), level=6),
                 screening_id),
            )

    def pause_screening(
        self,
        screening_id: str,
        request_id: str,
        reason: str,
        searched: list[str],
        waited_seconds: float = 0.0,
        estimate_seconds: int = 0,
        stage: str = "Still running at NCBI",
    ) -> None:
        """Record that a search is still running at NCBI and can be collected later.

        A fourth status beside running, complete and error, because "we stopped waiting" is
        none of those: the job exists, nothing failed, and the request id is the handle to
        finish it (ADR 0024). What is stored is only enough to resume, never a partial result,
        which would be a finding nobody made.
        """
        held = self.screening(screening_id).get("result") or {}
        waiting = {
            "parameters": held.get("parameters"),
            "submitted_at": held.get("submitted_at"),
            "candidates_searched": searched,
            "waited_seconds": waited_seconds,
            "estimate_seconds": estimate_seconds,
            "request_id": request_id,
        }
        with self.connect() as connection:
            connection.execute(
                "UPDATE screenings SET status='waiting', stage=?, "
                "request_id=?, error=?, result=? WHERE id=?",
                (stage, request_id, reason,
                 zlib.compress(json.dumps(waiting, separators=(",", ":")).encode(), level=6),
                 screening_id),
            )

    def resume_screening(self, screening_id: str) -> bool:
        """Mark a waiting screening as being collected again.

        Without this the record stays ``waiting`` while a worker is already collecting it, so
        the interface offers to collect it a second time and a researcher can queue the same
        work twice over.
        """
        with self.connect() as connection:
            updated = connection.execute(
                "UPDATE screenings SET status='running', stage='Checking NCBI again', "
                "error=NULL WHERE id=? AND status='waiting'",
                (screening_id,),
            )
            return updated.rowcount == 1

    def fail_screening(self, screening_id: str, reason: str) -> None:
        """Record that a screening failed, keeping the reason for the researcher."""
        with self.connect() as connection:
            connection.execute(
                "UPDATE screenings SET status='error', stage='Failed', error=? WHERE id=?",
                (reason, screening_id),
            )

    def screening(self, screening_id: str) -> dict[str, Any]:
        """Return one screening, with its result if it has finished.

        Raises:
            KeyError: If there is no such screening.
        """
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM screenings WHERE id=?", (screening_id,)).fetchone()
        if row is None:
            raise KeyError(f"No screening {screening_id!r} in this local store.")
        result = json.loads(zlib.decompress(row["result"])) if row["result"] else None
        local = result is not None and result.get("backend") == "local"
        record = {
            "id": row["id"], "run_id": row["run_id"], "scope": row["scope"],
            "created": row["created"], "status": row["status"], "stage": row["stage"],
            "request_id": row["request_id"], "error": row["error"],
            "progress": json.loads(row["progress"]) if row["progress"] else None,
            "poll_after_seconds": (1 if local else 5) if row["status"] == "running" else 60,
        }
        if row["result"] is not None:
            record["result"] = json.loads(zlib.decompress(row["result"]))
        return record

    def save_screening_annotations(
        self, screening_id: str, annotations: dict[str, Any],
    ) -> None:
        """Merge optional gene mappings without replacing the original search provenance."""
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute("SELECT result FROM screenings WHERE id=? AND "
                                     "status='complete'", (screening_id,)).fetchone()
            if row is None:
                raise ValueError("Annotations require a completed screening.")
            result = json.loads(zlib.decompress(row["result"]))
            result["gene_annotations"] = {**result.get("gene_annotations", {}), **annotations}
            connection.execute("UPDATE screenings SET result=? WHERE id=?", (
                zlib.compress(json.dumps(result).encode()), screening_id,
            ))

    def screenings_for(self, run_id: str) -> list[dict[str, Any]]:
        """Return every screening of one run, newest first, without their full results.

        Several screenings of one run are expected: a transcript search and a genomic one, and
        a repeat months later when the databases have moved on. None replaces another.
        """
        with self.connect() as connection:
            rows = connection.execute(
                "SELECT id,run_id,scope,created,status,stage,request_id,error,"
                "       (result IS NOT NULL) AS has_result "
                "FROM screenings WHERE run_id=? ORDER BY created DESC,rowid DESC", (run_id,)
            ).fetchall()
        return [dict(row) for row in rows]

    def list_runs(
        self, search: str = "", status: str = "", offset: int = 0, limit: int = 100,
    ) -> list[dict[str, Any]]:
        """Search saved labels, references and identifiers with bounded pagination."""
        if status not in ("", "complete", "error", "running"):
            raise ValueError("Unknown run status.")
        if len(search) > 300 or offset < 0 or not 1 <= limit <= 100:
            raise ValueError("Use a search up to 300 characters and a valid page range.")
        with self.connect() as connection:
            rows = connection.execute(
                "SELECT id,created,label,status,stage,error,summary FROM runs "
                "WHERE (?='' OR status=?) AND instr(lower(label || ' ' || id || ' ' || "
                "COALESCE(summary,'')),lower(?))>0 "
                "ORDER BY created DESC,rowid DESC LIMIT ? OFFSET ?",
                (status, status, search, limit, offset),
            ).fetchall()
        return [{**dict(row), "summary": json.loads(row["summary"]) if row["summary"] else None}
                for row in rows]

    def shortlist(self, run_id: str) -> dict[str, Any]:
        """Return researcher notes separately from the frozen scientific snapshot."""
        self.status(run_id)
        with self.connect() as connection:
            row = connection.execute(
                "SELECT revision,updated,entries FROM shortlists WHERE run_id=?", (run_id,),
            ).fetchone()
        return ({"revision": row["revision"], "updated": row["updated"],
                 "entries": json.loads(row["entries"])} if row else
                {"revision": 0, "updated": None, "entries": []})

    def save_shortlist(
        self, run_id: str, revision: int, entries: list[dict[str, str]],
    ) -> dict[str, Any] | None:
        """Use an atomic revision check to prevent one browser overwriting another's notes."""
        updated = datetime.now(UTC).isoformat()
        with self.connect() as connection:
            connection.execute(
                "INSERT OR IGNORE INTO shortlists VALUES(?,0,NULL,'[]')", (run_id,),
            )
            changed = connection.execute(
                "UPDATE shortlists SET revision=revision+1,updated=?,entries=? "
                "WHERE run_id=? AND revision=?",
                (updated, json.dumps(entries), run_id, revision),
            ).rowcount
        return ({"revision": revision + 1, "updated": updated, "entries": entries}
                if changed else None)
