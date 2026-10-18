"""Analysis exports preserve complete counts, retained subsets, provenance and unknowns."""

import csv
import io
import json
from pathlib import Path
from typing import Any

from sarna.app.screening_exports import screening_sheets, screening_summary
from sarna.app.web import AppRouter
from tests.unit.test_app_screening import completed_run
from tests.unit.test_workbook import sheets_of


def record(identifier: str, candidate: str, total: int, *, reference: str = "pack-A",
           created: str = "2026-10-18T10:00:00Z") -> dict[str, Any]:
    """Return a complete search with one detail but independently enumerated counts."""
    return {"id": identifier, "scope": "genomic", "status": "complete", "created": created,
            "result": {"backend": "local", "reference": {"release": reference},
                       "parameters": {"database": reference}, "searched_at": created,
                       "software": {"version": "test"}, "notes": ["Full counts, retained subset"],
                       "candidates_searched": [candidate], "proximities": [],
                       "specificity": [{"candidate": candidate, "hit_count": total,
                                        "identity_counts": {"19": 0, "18": 0, "17": 0,
                                                            "16": total},
                                        "details_retained": 1 if total else 0,
                                        "enumeration_complete": True, "alignments_capped": False,
                                        "perfect_hits": 0, "near_perfect_hits": 0,
                                        "wang_homology_hits": total, "on_target_hits": 1}],
                       "hits": [{"candidate": candidate, "accession": "NC_test",
                                 "title": "reference & retained detail", "identities": 16,
                                 "alignment_length": 19, "gaps": 0, "subject_start": 21,
                                 "subject_end": 39, "subject_strand": "plus", "query_start": 1,
                                 "query_end": 19, "mismatch_positions": [2, 5, 8]}]
                       if total else []}}


def test_reference_groups_latest_records_unknown_and_zero_are_distinct() -> None:
    candidates = [{"id": "c1", "rank": 1, "rating": {"met": 6, "evaluated": 7}},
                  {"id": "c2", "rank": 2}, {"id": "c3", "rank": 3}]
    old = record("old", "c1", 900, created="2026-10-17T10:00:00Z")
    newest = record("new", "c1", 100)
    alternate = record("alternate", "c1", 30, reference="pack-B")
    zero = record("zero", "c2", 0)
    pending = {"id": "pending", "scope": "genomic", "status": "running",
               "result": {"backend": "local", "reference_release": "pack-A",
                          "candidates_searched": ["c3"]}}
    records = [newest, old, alternate, zero, pending]
    rows = screening_summary("run", candidates, records)
    genomic = [row for row in rows if row["scope"] == "genomic"]
    assert [(r["reference"], r["reported_alignments"]) for r in genomic
            if r["candidate"] == "c1"] == [("pack-A", 100), ("pack-B", 30)]
    row = next(r for r in genomic if r["candidate"] == "c2")
    assert row["reported_alignments"] == 0 and row["matches_16_of_19"] == 0
    row = next(r for r in genomic if r["candidate"] == "c3")
    assert row["status"] == "running" and row["reported_alignments"] is None
    assert row["matches_19_of_19"] is None
    assert all(r["status"] == "not_screened" and r["reported_alignments"] is None
               for r in rows if r["scope"] == "transcript")
    sheets = {sheet.name: sheet for sheet in screening_sheets("run", candidates, records)}
    assert len(sheets["Retained alignments"].rows) == 3  # Historical searches remain traceable.
    assert len(sheets["Screening searches"].rows) == 5
    assert newest["result"]["specificity"][0]["hit_count"] == 100


def test_csv_excel_json_and_html_include_the_recorded_analysis(tmp_path: Path) -> None:
    service, run_id, snapshot = completed_run(tmp_path)
    candidate = next(c for c in snapshot["candidates"] if c["rank"] == 1)["id"]
    search_id, _ = service.store.reserve_screening(
        run_id, "genomic", [candidate], backend="local", reference_release="pack-A")
    evidence = record(search_id, candidate, 100)["result"]
    service.store.finish_screening(search_id, evidence)
    router = AppRouter(service)
    try:
        code, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.json")
        assert code == 200
        exported = json.loads(data)
        assert exported["screenings"][0]["id"] == search_id
        assert exported["screenings"][0]["result"]["hits"] == evidence["hits"]
        assert exported["screening_export_notes"]
        assert exported["screening_order"]["rows"]
        assert any(r["candidate"] == candidate and r["reported_alignments"] == 100
                   for r in exported["screening_summary"])
        _, _, response = router.dispatch("GET", f"/api/runs/{run_id}")
        assert "screenings" not in json.loads(response)["result"]
        code, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.screening.csv"
                                       f"?ids={candidate}")
        assert code == 200
        rows = list(csv.DictReader(io.StringIO(data.decode())))
        assert len(rows) == 2 and {r["candidate"] for r in rows} == {candidate}
        genomic = next(r for r in rows if r["scope"] == "genomic")
        assert genomic["reported_alignments"] == "100"
        assert genomic["matches_16_of_19"] == "100"
        assert genomic["details_retained"] == "1"
        transcript = next(r for r in rows if r["scope"] == "transcript")
        assert transcript["status"] == "not_screened"
        assert transcript["reported_alignments"] == ""
        code, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.xlsx?ids={candidate}")
        assert code == 200
        sheets = sheets_of(data)
        assert len(sheets["Retained alignments"]) == 2
        assert len(sheets["Off-target summary"]) == 3
        assert len(sheets["Screening searches"]) == 2
        assert sheets["Screening readme"]
        code, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.html")
        assert code == 200
        assert b'"screening_records"' in data and b'"screening_summary"' in data
        assert search_id.encode() in data
        assert service.store.snapshot(run_id) == snapshot
    finally:
        router.close()
