"""Exercise real local jobs and exports using recorded sources, without network or sockets."""

import csv
import io
import json
import sqlite3
import zlib
from copy import deepcopy
from pathlib import Path
from typing import Any

import httpx
import pytest

from sarna.app.models import DesignRequest
from sarna.app.service import DesignService
from sarna.app.store import RunStore
from sarna.app.web import AppRouter
from sarna.genome.ensembl import EnsemblCacheMissError, EnsemblClient
from sarna.offtarget.ncbi import BlastClient
from sarna.offtarget.screening import DEFAULT_CANDIDATES, MAX_CANDIDATES
from tests.unit.test_design_pipeline import PROMOTER_BASES
from tests.unit.test_ensembl import (
    AUTO_SEQUENCE,
    GENE_PAYLOAD,
    GENE_ROUTE,
    TRANSCRIPTS_PAYLOAD,
    TRANSCRIPTS_ROUTE,
    client_for,
    requests_made,
)
from tests.unit.test_offtarget_ncbi import READY_NO_HITS
from tests.unit.test_workbook import sheets_of

# Known 19-nt target from the cached SCN1A run; deliberately CpG-depleted flanks.
DNA = "GCTCTACTCTCAGTTCCAAACTATGACCTACTATGATACCATAC"


def offline_blast() -> BlastClient:
    """A BLAST client that answers from memory and can never reach NCBI.

    Injected by default, not only where a test is about screening. The route hands ``screen``
    to a worker pool without a client, so any test that merely *queues* a screening would
    otherwise build a real one and submit a live job: the suite used to do exactly that, and
    the seam only helps if nothing has to remember to use it.
    """
    def handler(_request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=READY_NO_HITS)

    return BlastClient(client=httpx.Client(transport=httpx.MockTransport(handler)))


def service_for(tmp_path: Path, routes: dict[str, Any] | None = None) -> DesignService:
    client = client_for(routes or {
        GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
        "/sequence/region/": {"seq": AUTO_SEQUENCE}, "/overlap/region/": [],
    })
    service = DesignService(tmp_path, client_factory=lambda _: client)
    service.blast_factory = offline_blast
    service.poll_sleep = lambda _seconds: None
    return service


def run_request(service: DesignService, **kwargs: Any) -> dict[str, Any]:  # noqa: ANN401
    values = {"source": "sequence", "sequence": DNA, "regulatory": False, "motifs": False}
    values.update(kwargs)
    request = DesignRequest.model_validate(values)
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    status = service.store.status(run_id)
    assert status["status"] == "complete", status["error"]
    return service.store.snapshot(run_id)


def test_gene_lookup_has_transcripts_and_uncertain_literature_references(tmp_path: Path):
    lookup = service_for(tmp_path).lookup("SCN1A")
    assert lookup["default_transcript"] == "ENST00000674923"
    assert len(lookup["transcripts"]) == 3
    distal = next(a for a in lookup["anchors"] if "distal" in a["id"])
    assert distal["reference_kind"] == "inferred"
    assert distal["caveat"]
    assert "SCN1A" in distal["label"]


@pytest.mark.strand
def test_gene_design_uses_requested_minus_strand_promoter_and_correct_coordinates(tmp_path: Path):
    result = run_request(service_for(tmp_path), source="gene", symbol="SCN1A",
                         sequence="", upstream=25, downstream=5)
    # MANE TSS is 166128019 (0-based). Minus-strand 25 upstream + 5 downstream.
    assert result["promoter"]["start"] == 166128015
    assert result["promoter"]["end"] == 166128045
    assert result["promoter"]["strand"] == "-"
    assert result["promoter"]["sequence"] == "TTGTTTGTTTGTTTGTTTGTTTGTTTGTTT"
    first = next(c for c in result["candidates"] if c["offset"] == 0)
    assert first["tss_span"] == [-25, -7]
    assert first["genomic"]["start"] == 166128026
    assert first["genomic"]["end"] == 166128045
    assert result["app"]["anchor"]["gene"] == "SCN1A"


def test_second_gene_is_retrieved_instead_of_reusing_a_demo(tmp_path: Path):
    gene = {**GENE_PAYLOAD, "display_name": "CDKN1A", "id": "ENSG_TEST", "strand": 1}
    transcripts = deepcopy(TRANSCRIPTS_PAYLOAD)
    transcripts["id"] = "ENSG_TEST"
    for t in transcripts["Transcript"]:
        t["strand"] = 1
    service = service_for(tmp_path, {
        "/lookup/symbol/homo_sapiens/CDKN1A": gene,
        "/lookup/id/ENSG_TEST": transcripts,
        "/sequence/region/": {"seq": AUTO_SEQUENCE},
    })
    assert service.lookup("CDKN1A")["symbol"] == "CDKN1A"
    result = run_request(service, source="gene", symbol="CDKN1A", sequence="",
                         upstream=25, downstream=5)
    assert result["gene"]["symbol"] == "CDKN1A"
    assert result["promoter"]["sequence"] == "AAACAAACAAACAAACAAACAAACAAACAA"
    assert result["promoter"]["strand"] == "+"


@pytest.mark.parametrize("symbol", ["CHD2", "KCNQ2"])
def test_new_gene_fetches_then_reopens_offline_regardless_of_case(tmp_path: Path, symbol: str):
    """Recorded-source integration: cache misses never fetch, online mode populates the cache."""
    gene = {**GENE_PAYLOAD, "display_name": symbol}
    clients: list[EnsemblClient] = []

    def factory(offline: bool) -> EnsemblClient:
        client = client_for({
            f"/lookup/symbol/homo_sapiens/{symbol}": gene,
            TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
            "/sequence/region/": {"seq": AUTO_SEQUENCE},
        }, cache_dir=tmp_path / "cache", offline=offline)
        clients.append(client)
        return client

    service = DesignService(tmp_path / "runs", factory)
    router = AppRouter(service)
    try:
        code, _, body = router.dispatch("GET", f"/api/genes/{symbol.lower()}?offline=true")
        assert code == 409
        assert json.loads(body)["code"] == "offline-cache-miss"
        assert "SARNA_OFFLINE" not in json.loads(body)["detail"]
        assert requests_made(clients[-1]) == []
        code, _, body = router.dispatch("GET", f"/api/genes/{symbol.lower()}?offline=false")
        assert code == 200
        assert json.loads(body)["symbol"] == symbol
        assert len(requests_made(clients[-1])) == 2
        assert f"/homo_sapiens/{symbol}?" in requests_made(clients[-1])[0]
        online = run_request(service, source="gene", sequence="", symbol=symbol.lower(),
                             upstream=25, downstream=5, offline=False)
        assert online["app"]["request"]["symbol"] == symbol
        assert online["gene"]["symbol"] == symbol
        assert len(requests_made(clients[-1])) == 1  # Only the missing promoter window.
        assert service.lookup(symbol.lower(), offline=True)["symbol"] == symbol
        assert requests_made(clients[-1]) == []
        cached = run_request(service, source="gene", sequence="", symbol=symbol.lower(),
                             upstream=25, downstream=5, offline=True)
        assert cached["candidates"] == online["candidates"]
        assert cached["app"]["data_access"]["offline"] is True
        assert requests_made(clients[-1]) == []
        # Cached gene/transcripts do not imply that every adjustable window is cached.
        request = DesignRequest(symbol=symbol.lower(), upstream=26, downstream=5,
                                offline=True, regulatory=False, motifs=False)
        run_id = service.store.create(request.model_dump(mode="json"))
        service.run(run_id, request)
        status = service.store.status(run_id)
        assert status["status"] == "error"
        assert status["error_code"] == "offline-cache-miss"
        assert "sequence/region" in status["error"]
        assert requests_made(clients[-1]) == []
    finally:
        router.close()


def test_environment_forced_offline_is_visible_and_cannot_be_bypassed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setenv("SARNA_OFFLINE", "1")
    monkeypatch.setenv("SARNA_CACHE_DIR", str(tmp_path / "empty-cache"))
    service = DesignService(tmp_path / "runs")
    router = AppRouter(service)
    try:
        _, _, body = router.dispatch("GET", "/api/catalogue")
        assert json.loads(body)["data_access"]["offline_forced"] is True
        code, _, body = router.dispatch("GET", "/api/genes/chd2?offline=false")
        assert code == 409
        assert "SARNA_OFFLINE forces this mode" in json.loads(body)["detail"]
        assert "restart" in json.loads(body)["detail"]
        with pytest.raises(EnsemblCacheMissError):
            service.lookup("chd2", offline=False)
    finally:
        router.close()


def test_unknown_gene_is_reported_as_missing_instead_of_an_outage(tmp_path: Path):
    router = AppRouter(service_for(tmp_path))
    try:
        code, _, body = router.dispatch("GET", "/api/genes/UNKNOWN")
        assert code == 404
        assert json.loads(body)["code"] == "gene-not-found"
    finally:
        router.close()


def test_gene_without_mane_requires_explicit_transcript(tmp_path: Path):
    transcripts = deepcopy(TRANSCRIPTS_PAYLOAD)
    transcripts["Transcript"][0].pop("MANE")
    service = service_for(tmp_path, {GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: transcripts,
                                   "/sequence/region/": {"seq": AUTO_SEQUENCE}})
    assert service.lookup("SCN1A")["default_transcript"] is None
    request = DesignRequest(symbol="SCN1A", upstream=25, downstream=5,
                            regulatory=False, motifs=False)
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    assert "Choose a transcript" in service.store.status(run_id)["error"]
    explicit = run_request(service, source="gene", symbol="SCN1A", sequence="",
                           transcript_id="ENST00000409050", upstream=25, downstream=5)
    assert explicit["transcript"]["ensembl_id"] == "ENST00000409050"


def test_pasted_fasta_preserves_missing_context_and_tristate_checks(tmp_path: Path):
    result = run_request(service_for(tmp_path), sequence=">my promoter\n" + DNA)
    assert result["promoter"]["sequence"] == DNA
    assert result["promoter"]["contig"] is None
    c = next(c for c in result["candidates"] if c["offset"] == 0)
    assert "tss-window" in c["rating"]["unchecked_gates"]
    assert c["rating"]["unavailable"] == 1  # TSS preference, not a failed preference.
    assert "tss-window" not in c["excluded_by"]
    assert c["rating"]["guide"] == "antisense"
    assert c["antisense_strand_rna"] == "UUGGAACUGAGAGUAGAGC"
    assert any("No gene-expression" in n for n in result["app"]["notices"])


def test_relative_pasted_reference_draws_tss_with_no_position_zero(tmp_path: Path):
    result = run_request(service_for(tmp_path), tss_start=-2)
    assert result["promoter"]["tss_offset"] == 2
    first = next(c for c in result["candidates"] if c["offset"] == 0)
    assert first["tss_span"] == [-2, 17]
    assert next(r for r in first["rules"] if r["id"] == "tss-window")["outcome"] == "fail"


def test_changed_rule_role_and_numeric_threshold_change_actual_results(tmp_path: Path):
    service = service_for(tmp_path)
    original = run_request(service, tss_start=-1500)
    changed = run_request(service, tss_start=-1500,
                          roles={"tss-window": "preference", "tss-preferred-window": "off"},
                          overrides={"gc-content": {"maximum-gc-fraction": 0.65}})
    assert original["counts"]["passing"] == 0
    assert changed["counts"]["passing"] > 0
    assert changed["configuration"]["is_custom"]
    assert changed["configuration"]["baseline_rules"] == original["configuration"]["rules"]
    assert service.store.snapshot(original["app"]["run_id"])["counts"]["passing"] == 0


def test_named_parameter_alternative_retains_option_provenance(tmp_path: Path):
    result = run_request(service_for(tmp_path), options={"gc-content": "reynolds-2004"})
    rule = next(r for r in result["configuration"]["rules"] if r["id"] == "gc-content")
    assert rule["selected_option"]["id"] == "reynolds-2004"
    assert rule["parameters"][0]["value"] == 0.3
    assert rule["parameters"][0]["origin"] == "named-option"


def test_custom_end_length_controls_displayed_guide_and_all_model_comparisons(tmp_path: Path):
    # Hand-calculated RNA dinucleotide energies for two-base ends of this target:
    # GC = -3.42; AA = -0.93 kcal/mol. E3-E5 = +2.49, hence antisense.
    result = run_request(service_for(tmp_path),
                         overrides={"thermodynamic-asymmetry": {"end-length-nt": 2}})
    first = next(c for c in result["candidates"] if c["offset"] == 0)
    assessment = first["guide_prediction"]["assessments"][0]
    assert first["guide_prediction"]["end_length_nt"] == 2
    assert assessment["difference"] == pytest.approx(2.49)
    assert first["rating"]["guide"] == "antisense"
    gate = next(r for r in first["rules"] if r["id"] == "thermodynamic-asymmetry")
    assert gate["measurements"]["asymmetry"] == pytest.approx(2.49)


def test_source_failure_retains_failed_job_instead_of_empty_success(tmp_path: Path):
    service = service_for(tmp_path, {"/lookup/symbol/": 503})
    request = DesignRequest(symbol="SCN1A", motifs=False, regulatory=False)
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    status = service.store.status(run_id)
    assert status["status"] == "error"
    assert status["error"]
    with pytest.raises(KeyError):
        service.store.snapshot(run_id)


@pytest.mark.parametrize("change", [
    {"overrides": {"gc-content": {"minimum-gc-fraction": 0.8}}},
    {"overrides": {"target-length-19": {"target-length": 21}}},
    {"overrides": {"tss-window": {"window-end-bp": 0}}},
    {"roles": {"asymmetry-margin": "mandatory"}},
    {"roles": {"target-length-19": "off"}},
])
def test_invalid_scientific_configuration_never_queues_a_job(
    tmp_path: Path, change: dict[str, Any],
):
    router = AppRouter(service_for(tmp_path))
    try:
        status, _, _ = router.dispatch("POST", "/api/runs", json.dumps({
            "source": "sequence", "sequence": DNA, **change,
        }).encode())
        assert status == 422
        assert router.service.store.list_runs() == []
    finally:
        router.close()


def test_saved_run_and_candidate_details_reopen_without_engine_or_source(tmp_path: Path):
    service = service_for(tmp_path)
    saved = run_request(service)
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        status, _, data = router.dispatch("GET", f"/api/runs/{run_id}")
        assert status == 200
        result = json.loads(data)["result"]
        assert result["counts"] == saved["counts"]
        assert "rules" not in result["candidates"][0]
        originals = {c["id"]: c for c in saved["candidates"]}
        for candidate in result["candidates"]:
            if candidate["rank"] is not None and candidate["rank"] <= 15:
                assert candidate["preference_verdicts"] == [
                    rule for rule in originals[candidate["id"]]["rules"]
                    if rule["kind"] == "preference"
                ]
            else:
                assert "preference_verdicts" not in candidate
        assert "preference_verdicts" not in service.store.snapshot(run_id)["candidates"][0]
        cid = saved["candidates"][0]["id"]
        _, _, detail = router.dispatch("GET", f"/api/runs/{run_id}/candidates/{cid}")
        assert json.loads(detail)["rules"] == saved["candidates"][0]["rules"]
        with pytest.raises(KeyError):
            service.store.snapshot("unknown-run")
    finally:
        router.close()


def test_csv_exports_selected_core_sequences_even_for_rejected_candidates(tmp_path: Path):
    service = service_for(tmp_path)
    saved = run_request(service, tss_start=-1500)
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        status, media, body = router.dispatch("GET", f"/api/runs/{run_id}/export.csv?ids=saRNA-001")
        assert (status, media) == (200, "text/csv; charset=utf-8")
        rows = list(csv.DictReader(io.StringIO(body.decode())))
        assert len(rows) == 1
        assert rows[0]["rank"] == ""
        assert rows[0]["sense_rna_core_5to3"] == "GCUCUACUCUCAGUUCCAA"
        assert rows[0]["antisense_rna_core_5to3"] == "UUGGAACUGAGAGUAGAGC"
        assert rows[0]["core_length_nt"] == "19"
        assert router.dispatch("GET", f"/api/runs/{run_id}/export.csv?ids=unknown")[0] == 422
    finally:
        router.close()


def test_both_csv_exports_are_written_in_rank_order(tmp_path: Path):
    # A snapshot's candidates are in *position* order carrying their ranks, which is right for
    # the promoter map and for a table the browser sorts. A CSV is neither: it is a flat file a
    # researcher reads top-down and pastes into an order form. Emitting rank 8 on the first row
    # means "take the top three" orders the wrong three molecules.
    service = service_for(tmp_path)
    # This promoter ranks its survivors in the order 3, 4, 5, 2, 1 by position, so a CSV that
    # follows the snapshot rather than the rank is visibly wrong. The default fixture has two
    # passing candidates that are already in rank order and cannot show the difference.
    saved = run_request(service, sequence=PROMOTER_BASES, tss_start=-200)
    assert [c["rank"] for c in saved["candidates"] if c["rank"] is not None] != sorted(
        c["rank"] for c in saved["candidates"] if c["rank"] is not None
    ), "fixture no longer scrambles rank against position, so this test proves nothing"
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        for route in ("export.csv", "export.order.csv"):
            body = router.dispatch("GET", f"/api/runs/{run_id}/{route}")[2]
            ranks = [int(r["rank"]) for r in csv.DictReader(io.StringIO(body.decode()))]
            assert ranks, f"{route} exported no rows"
            assert ranks == sorted(ranks), f"{route} is not in rank order: {ranks}"
            assert ranks[0] == 1
    finally:
        router.close()


def test_the_order_csv_writes_the_dtdt_duplex_and_its_caveats(tmp_path: Path):
    # The order export is distinct from the core-sequence CSV: it carries the 21-mers that are
    # actually synthesised, and the caveats travel in the same row as the sequences.
    service = service_for(tmp_path)
    saved = run_request(service, tss_start=-1500)
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        status, media, body = router.dispatch(
            "GET", f"/api/runs/{run_id}/export.order.csv?ids=saRNA-001"
        )
        assert (status, media) == (200, "text/csv; charset=utf-8")
        rows = list(csv.DictReader(io.StringIO(body.decode())))
        assert len(rows) == 1
        assert rows[0]["sense_order_5to3"] == "GCUCUACUCUCAGUUCCAA[dT][dT]"
        assert rows[0]["antisense_order_5to3"] == "UUGGAACUGAGAGUAGAGC[dT][dT]"
        assert rows[0]["ordered_length_nt"] == "21"
        assert rows[0]["core_length_nt"] == "19"
        assert "off-target-screening [not-performed]" in rows[0]["sanity_report"]
        assert router.dispatch("GET", f"/api/runs/{run_id}/export.order.csv?ids=unknown")[0] == 422
    finally:
        router.close()


def test_every_run_carries_the_rankings_measured_performance(tmp_path: Path):
    # A reader takes a rank for a prediction unless told otherwise. The run record therefore
    # carries the engine's own measurement, so an archived run shows the same caveat the live
    # page does rather than the caveat living only in the UI.
    saved = run_request(service_for(tmp_path))
    measured = saved["ranking_performance"]
    assert measured["preset"] == "wang-2015"
    assert measured["beats_chance"] is False
    assert measured["dataset"]["assayed"] == 125
    assert measured["applies"] is True
    assert "chance" in measured["headline"].lower()
    assert "range-restricted" in measured["caveat"].lower()


def test_a_custom_configuration_does_not_inherit_the_presets_measurement(tmp_path: Path):
    # The AUC describes the unmodified preset. Once a researcher changes a rule they are not
    # running what was measured, and claiming the number for their configuration would attach
    # evidence to something nobody assessed.
    saved = run_request(service_for(tmp_path), roles={"tss-window": "preference"})
    measured = saved["ranking_performance"]
    assert saved["configuration"]["is_custom"] is True
    assert measured["applies"] is False
    assert "custom configuration" in measured["custom_note"].lower()


def test_the_record_carries_the_same_order_rows_the_route_serves(tmp_path: Path):
    # The frozen archive has no backend, so the duplexes it hands over are the ones stored here
    # by Python. If these two ever diverge, an archived run and a live run would name different
    # molecules for the same candidate. JavaScript joins these strings; it never builds them.
    service = service_for(tmp_path)
    saved = run_request(service, sequence=PROMOTER_BASES, tss_start=-200)
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        body = router.dispatch("GET", f"/api/runs/{run_id}/export.order.csv")[2]
        served = list(csv.reader(io.StringIO(body.decode())))
        assert served[0] == saved["order"]["columns"]
        assert served[1:] == saved["order"]["rows"]
        assert all(row[4].endswith("[dT][dT]") for row in saved["order"]["rows"])
    finally:
        router.close()


def test_archive_has_inline_ui_complete_verdicts_and_escaped_data(tmp_path: Path):
    service = service_for(tmp_path)
    saved = run_request(service, cell_line="</script><script>alert('injected')</script>")
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        status, _, page = router.dispatch("GET", f"/api/runs/{run_id}/export.html")
        assert status == 200
        assert b"window.SARNA_FROZEN=" in page
        assert b'<script src="/app.js"' not in page
        assert b'<link rel="stylesheet"' not in page
        assert b'href="/favicon.svg"' not in page
        assert b'class="brand" href="/"' not in page
        assert b"</script><script>alert('injected')" not in page
        _, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.json")
        archived = json.loads(data)
        assert archived.pop("workspace")["shortlist"]["entries"] == []
        assert archived == saved
    finally:
        router.close()


def test_thesis_option_labels_and_evaluates_its_actual_upstream_extent(tmp_path: Path):
    result = run_request(service_for(tmp_path), tss_start=-45,
                         options={"tss-window": "vondervoort-2026"})
    first = next(c for c in result["candidates"] if c["offset"] == 0)
    gate = next(r for r in first["rules"] if r["id"] == "tss-window")
    assert gate["outcome"] == "pass"
    rule = next(r for r in result["configuration"]["rules"] if r["id"] == "tss-window")
    assert rule["parameters"][1]["value"] == -1
    assert "upstream part" in rule["selected_option"]["label"]
    assert "does not reproduce the full thesis" in rule["selected_option"]["note"]


def test_store_recovers_interrupted_jobs_without_changing_complete_records(tmp_path: Path):
    service = service_for(tmp_path)
    saved = run_request(service)
    interrupted = service.store.create({"symbol": "TP53"})
    reopened = RunStore(service.store.path)
    assert reopened.status(interrupted)["status"] == "error"
    assert reopened.snapshot(saved["app"]["run_id"]) == saved


def test_store_migrates_existing_runs_and_preserves_their_snapshot(tmp_path: Path):
    path = tmp_path / "runs.sqlite3"
    snapshot = {"saved": "unchanged"}
    with sqlite3.connect(path) as connection:
        connection.execute(
            "CREATE TABLE runs (id TEXT PRIMARY KEY, created TEXT NOT NULL, "
            "label TEXT NOT NULL, status TEXT NOT NULL, stage TEXT NOT NULL, "
            "request TEXT NOT NULL, snapshot BLOB, error TEXT)"
        )
        connection.execute("PRAGMA user_version=1")
        connection.execute("INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)", (
            "old", "2026-10-07", "CHD2", "complete", "Complete", '{"symbol":"CHD2"}',
            zlib.compress(json.dumps(snapshot).encode()), None,
        ))
    store = RunStore(path)
    assert store.snapshot("old") == snapshot
    assert store.status("old")["request"] == {"symbol": "CHD2"}
    assert store.status("old")["error_code"] is None
    with store.connect() as connection:
        # Researcher shortlists are separate from frozen runs in schema 6.
        assert connection.execute("PRAGMA user_version").fetchone()[0] == 6
        screening_columns = {
            row[1] for row in connection.execute("PRAGMA table_info(screenings)")
        }
        assert "progress" in screening_columns
        tables = {
            row[0]
            for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }
        assert {"runs", "batches", "batch_runs", "screenings"} <= tables


def test_regulatory_unavailability_does_not_change_sequence_rank(tmp_path: Path):
    service = service_for(tmp_path, {GENE_ROUTE: GENE_PAYLOAD,
                                   TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
                                   "/sequence/region/": {"seq": AUTO_SEQUENCE},
                                   "/overlap/region/": 503})
    base = {"source": "gene", "symbol": "SCN1A", "sequence": "",
            "upstream": 25, "downstream": 5}
    missing = run_request(service, **base, regulatory=True)
    skipped = run_request(service, **base, regulatory=False)
    assert missing["candidates"] == skipped["candidates"]
    assert any("Regulatory annotations unavailable" in n for n in missing["app"]["notices"])


def test_local_routes_serve_real_setup_and_reject_unknown_files(tmp_path: Path):
    router = AppRouter(service_for(tmp_path))
    try:
        assert b'"wang-2015"' in router.dispatch("GET", "/api/catalogue")[2]
        assert b'id="gene-symbol"' in router.dispatch("GET", "/")[2]
        assert b"Design saRNAs" in router.dispatch("GET", "/")[2]
        assert router.dispatch("GET", "/../CLAUDE.md")[0] == 404
        assert router.dispatch("GET", "/favicon.ico")[0] == 200
    finally:
        router.close()


def test_the_app_publishes_the_end_energies_behind_each_guide_call(tmp_path: Path):
    # Hand-calculated as in the end-length test above: for a two-base end of this target the
    # 5' stack is GC = -3.42 and the 3' stack AA = -0.93 kcal/mol (RNA, Xia 1998).
    # An interface cannot draw the asymmetry from the signed difference alone.
    result = run_request(service_for(tmp_path),
                         overrides={"thermodynamic-asymmetry": {"end-length-nt": 2}})
    first = next(c for c in result["candidates"] if c["offset"] == 0)
    assessment = next(a for a in first["guide_prediction"]["assessments"]
                      if a["table"] == "rna-xia-1998")
    assert assessment["five_prime"]["energy"] == pytest.approx(-3.42)
    assert assessment["five_prime"]["stacks"] == [["GC", pytest.approx(-3.42)]]
    assert assessment["three_prime"]["energy"] == pytest.approx(-0.93)
    assert assessment["three_prime"]["stacks"] == [["AA", pytest.approx(-0.93)]]
    assert assessment["unit"] == "kcal/mol"


class TestTheSpreadsheetExport:
    """"Download all the saRNAs with their values", as a file Excel opens correctly.

    A CSV is not a file format, it is a convention Excel resolves against the machine's locale:
    on a Dutch or German Windows the list separator is ``;`` and a comma-separated export
    arrives with every column in column A. The first thing a colleague would see.
    """

    @staticmethod
    def book(service: DesignService, run_id: str, query: str = "") -> dict[str, list[list[Any]]]:
        router = AppRouter(service)
        try:
            status, content_type, body = router.dispatch(
                "GET", f"/api/runs/{run_id}/export.xlsx{query}", b"")
            assert status == 200, body[:200]
            assert "spreadsheetml.sheet" in content_type
            return sheets_of(body)
        finally:
            router.close()

    def test_the_workbook_has_a_sheet_for_the_candidates_the_run_and_the_rules(
        self, tmp_path: Path
    ) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        sheets = self.book(service, saved["app"]["run_id"])
        assert list(sheets) == ["saRNA candidates", "Run", "Rules", "Reference sources",
                                "Shortlist"]

    def test_every_candidate_appears_in_rank_order_with_its_values(
        self, tmp_path: Path
    ) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        rows = self.book(service, saved["app"]["run_id"])["saRNA candidates"]
        header = rows[0]
        assert rows[1:], "the export has no candidates"
        ranks = [row[header.index("rank")] for row in rows[1:]]
        assert ranks == sorted(r for r in ranks if r is not None) or len(ranks) == 1
        first = dict(zip(header, rows[1], strict=True))
        # The sequences a researcher came for, and the ordered 21-mers beside the cores.
        assert set(str(first["target_dna_5to3"])) <= set("ACGT")
        assert str(first["sense_order_5to3"]).endswith("[dT][dT]")
        assert str(first["antisense_order_5to3"]).endswith("[dT][dT]")
        assert first["predicted_guide"] in {"antisense", "sense", "ambiguous"}

    def test_sequences_stay_text_so_excel_cannot_reinterpret_them(
        self, tmp_path: Path
    ) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        rows = self.book(service, saved["app"]["run_id"])["saRNA candidates"]
        column = rows[0].index("target_dna_5to3")
        assert all(isinstance(row[column], str) for row in rows[1:])

    def test_numbers_stay_numbers_so_a_column_can_be_sorted(self, tmp_path: Path) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        rows = self.book(service, saved["app"]["run_id"])["saRNA candidates"]
        column = rows[0].index("gc_percent")
        assert all(isinstance(row[column], float) for row in rows[1:])

    def test_the_run_sheet_carries_the_caveat_the_interface_shows(
        self, tmp_path: Path
    ) -> None:
        # The spreadsheet leaves the application. The measured ranking performance has to
        # leave with it, or the numbers travel without the one statement that qualifies them.
        service = service_for(tmp_path)
        saved = run_request(service)
        run = {row[0]: row[1] for row in self.book(service, saved["app"]["run_id"])["Run"][1:]}
        assert "Not a probability" in str(run["what a rating is"])
        assert run["off-target screening"] == "none performed"
        assert run["targets evaluated"]
        # A promoter reference that does not say where it came from is a coordinate with no
        # standing, and the spreadsheet is what leaves the application.
        assert run["reference provenance"], "the reference reached Excel without its source"
        assert run["reference kind"]

    def test_the_rules_sheet_names_every_rule_with_its_source(self, tmp_path: Path) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        rows = self.book(service, saved["app"]["run_id"])["Rules"]
        assert rows[0] == ["rule", "name", "role", "evidence", "source", "parameters",
                           "source URL", "evidence scope", "application note"]
        by_id = {row[0]: row for row in rows[1:]}
        assert by_id["pos19-a"][6] == "https://doi.org/10.14440/jbm.2015.39"
        assert "protocol guidance" in by_id["pos19-a"][7]
        assert "not evidence strength" in by_id["pos19-a"][8]
        sources = [row[4] for row in rows[1:]]
        assert all(sources), "a rule reached the spreadsheet without its provenance"
        assert any("Wang" in str(s) for s in sources)

    def test_a_selection_exports_only_what_was_selected(self, tmp_path: Path) -> None:
        service = service_for(tmp_path)
        saved = run_request(service)
        wanted = saved["candidates"][0]["id"]
        rows = self.book(service, saved["app"]["run_id"], f"?ids={wanted}")["saRNA candidates"]
        assert len(rows) == 2
        assert rows[1][rows[0].index("candidate")] == wanted


def test_identical_duplexes_are_named_in_the_order_csv(tmp_path: Path):
    # A promoter with the same window twice produces two candidates and one product. The route
    # must say so: "Order CSV, all 8" on a repeat-rich promoter otherwise buys one molecule
    # eight times, and nothing in the file would mention it.
    service = service_for(tmp_path)
    # A 62 bp tandem repeat, anchored so several copies fall inside Wang's -1000..-100 window
    # and therefore pass together. With a shorter anchor the window gate happens to keep only
    # one copy of each duplex, which is the quiet case this test is not about.
    repeated = ("GCTACGGTATCATGTCATAAGTTCAATTGACATTGGAATTACCATAGTTAACGTTTAGCATA" * 10)[:600]
    saved = run_request(service, sequence=repeated, tss_start=-600)
    router = AppRouter(service)
    try:
        run_id = saved["app"]["run_id"]
        body = router.dispatch("GET", f"/api/runs/{run_id}/export.order.csv")[2]
        rows = list(csv.DictReader(io.StringIO(body.decode())))
        exported = {r["candidate"] for r in rows}
        twinned = [
            r for r in rows
            if any(t in exported for t in r["identical_duplex_in_run"].split())
        ]
        assert twinned, "the repeated promoter put no two identical duplexes in one export"
        first = twinned[0]
        # Every id it names is a real candidate of this run, and never the row itself. The
        # twins span the whole run, so an excluded position is named too: the same molecule
        # sitting somewhere the gates rejected is worth knowing before ordering it.
        named = first["identical_duplex_in_run"].split()
        assert first["candidate"] not in named
        every_candidate = {c["id"] for c in saved["candidates"]}
        assert all(name in every_candidate for name in named)
        # Both twins are in this file, so it is a concern rather than a note.
        assert "identical-duplex [concern]" in first["sanity_report"]

        # Exporting one of the pair alone is not buying anything twice.
        alone = router.dispatch(
            "GET", f"/api/runs/{run_id}/export.order.csv?ids={first['candidate']}")[2]
        single = next(iter(csv.DictReader(io.StringIO(alone.decode()))))
        assert single["identical_duplex_in_run"] == first["identical_duplex_in_run"]
        assert "identical-duplex [ok]" in single["sanity_report"]
    finally:
        router.close()


class TestChoosingHowManyCandidatesToScreen:
    """ADR 0023: the fixed top 50 becomes a choice, validated here at the boundary.

    None of these reaches NCBI. The job is queued on a pool whose work is never run, which is
    enough to see what the route accepted and what it refused.
    """

    @staticmethod
    def _router(tmp_path: Path) -> tuple[AppRouter, str, list[str]]:
        service = service_for(tmp_path)
        saved = run_request(service)
        passing = [c["id"] for c in saved["candidates"] if c["rank"] is not None]
        return AppRouter(service), saved["app"]["run_id"], passing

    def _start(self, router: AppRouter, run_id: str, body: dict[str, Any]) -> tuple[int, Any]:
        status, _, data = router.dispatch(
            "POST", f"/api/runs/{run_id}/screenings", json.dumps(body).encode()
        )
        return status, json.loads(data)

    def test_the_default_is_fifteen_rather_than_everything(self, tmp_path: Path):
        router, run_id, passing = self._router(tmp_path)
        try:
            status, answer = self._start(router, run_id, {"scope": "transcript"})
            assert status == 202
            assert answer["candidates"] == min(DEFAULT_CANDIDATES, len(passing))
        finally:
            router.close()

    def test_a_count_is_honoured_and_capped_at_one_jobs_worth(self, tmp_path: Path):
        router, run_id, passing = self._router(tmp_path)
        try:
            assert self._start(router, run_id, {"count": 1})[1]["candidates"] == 1
            # A courtesy to NCBI that the interface cannot talk its way out of.
            asked = self._start(router, run_id, {"count": 10_000})[1]["candidates"]
            assert asked <= min(MAX_CANDIDATES, len(passing))
        finally:
            router.close()

    def test_an_explicit_selection_is_screened_as_given(self, tmp_path: Path):
        router, run_id, passing = self._router(tmp_path)
        try:
            chosen = passing[:2]
            status, answer = self._start(router, run_id, {"candidates": chosen})
            assert status == 202
            assert answer["candidates"] == len(chosen)
        finally:
            router.close()

    def test_a_candidate_that_did_not_pass_the_gates_is_named_not_dropped(
        self, tmp_path: Path
    ):
        # Quietly screening fewer than were asked for is indistinguishable from screening them
        # and finding nothing, which is the confusion this whole change exists to remove.
        service = service_for(tmp_path)
        saved = run_request(service)
        excluded = next(c["id"] for c in saved["candidates"] if c["rank"] is None)
        router = AppRouter(service)
        try:
            status, answer = self._start(
                router, saved["app"]["run_id"], {"candidates": [excluded]}
            )
            assert status == 422
            assert excluded in answer["detail"]
        finally:
            router.close()

    def test_an_unknown_candidate_is_refused(self, tmp_path: Path):
        router, run_id, _ = self._router(tmp_path)
        try:
            status, answer = self._start(router, run_id, {"candidates": ["saRNA-nope"]})
            assert status == 422
            assert "saRNA-nope" in answer["detail"]
        finally:
            router.close()

    def test_more_than_one_job_may_hold_is_refused_with_the_batch_advice(
        self, tmp_path: Path
    ):
        router, run_id, _ = self._router(tmp_path)
        try:
            status, answer = self._start(
                router, run_id, {"candidates": [f"saRNA-{i:04d}" for i in range(60)]}
            )
            assert status == 422
            assert "batches" in answer["detail"]
        finally:
            router.close()

    @pytest.mark.parametrize("body", [
        {"candidates": []},
        {"candidates": "saRNA-001"},
        {"candidates": [1, 2]},
        {"count": 0},
        {"count": -1},
        {"count": "fifteen"},
        {"count": True},
    ])
    def test_a_malformed_selection_never_queues_a_job(
        self, tmp_path: Path, body: dict[str, Any]
    ):
        router, run_id, _ = self._router(tmp_path)
        try:
            assert self._start(router, run_id, body)[0] == 422
            assert router.service.store.screenings_for(run_id) == []
        finally:
            router.close()


class TestASlowSearchIsCollectedRatherThanLost:
    """NCBI holds a request id for about a day, so running out of patience is not failing.

    The 1800s cap used to raise, the screening was recorded as an error, and the job NCBI was
    still working on became unreachable. A researcher lost the whole genomic lane and the only
    way forward was to queue an identical second job (ADR 0024).
    """

    @staticmethod
    def _waiting(tmp_path: Path) -> tuple[AppRouter, str, str]:
        service = service_for(tmp_path)
        saved = run_request(service)
        run_id = saved["app"]["run_id"]
        screening_id = service.store.start_screening(run_id, "genomic")
        searched = [c["id"] for c in saved["candidates"] if c["rank"] is not None][:2]
        service.store.pause_screening(
            screening_id, "RID-ABC-1", "NCBI search RID-ABC-1 has run for 1800s without "
            "finishing.", searched, waited_seconds=1800.0,
        )
        return AppRouter(service), run_id, screening_id

    def test_a_paused_screening_keeps_its_request_id_and_its_candidates(
        self, tmp_path: Path
    ):
        router, _run_id, screening_id = self._waiting(tmp_path)
        try:
            held = router.service.store.screening(screening_id)
            assert held["status"] == "waiting"
            assert held["request_id"] == "RID-ABC-1"
            assert held["result"]["candidates_searched"]
            assert held["result"]["waited_seconds"] == 1800.0
        finally:
            router.close()

    def test_it_is_neither_complete_nor_failed(self, tmp_path: Path):
        # Both would be wrong, and in opposite directions: one invents a finding, the other
        # throws away a job that still exists.
        router, run_id, _screening_id = self._waiting(tmp_path)
        try:
            listed = json.loads(router.dispatch("GET", f"/api/runs/{run_id}/screenings")[1:][1])
            assert [s["status"] for s in listed] == ["waiting"]
            # Nothing incomplete may reach the order export.
            body = router.dispatch("GET", f"/api/runs/{run_id}/export.order.csv")[2]
            assert b"not-performed" in body
        finally:
            router.close()

    def test_resuming_queues_a_collection_rather_than_a_second_search(self, tmp_path: Path):
        router, _run_id, screening_id = self._waiting(tmp_path)
        try:
            status, _, body = router.dispatch(
                "POST", f"/api/screenings/{screening_id}/resume")
            assert status == 202
            assert json.loads(body)["request_id"] == "RID-ABC-1"
        finally:
            router.close()

    def test_resuming_marks_it_running_so_it_cannot_be_collected_twice(
        self, tmp_path: Path
    ):
        # Left as 'waiting' while a worker is already collecting, the interface would offer
        # the button again and a researcher could queue the same work twice.
        router, _run_id, screening_id = self._waiting(tmp_path)
        try:
            router.service.store.resume_screening(screening_id)
            held = router.service.store.screening(screening_id)
            assert held["status"] == "running"
            assert held["error"] is None
            status, _, body = router.dispatch(
                "POST", f"/api/screenings/{screening_id}/resume")
            assert status == 409
            assert "already being collected" in json.loads(body)["detail"]
        finally:
            router.close()

    def test_a_finished_screening_cannot_be_resumed(self, tmp_path: Path):
        service = service_for(tmp_path)
        saved = run_request(service)
        screening_id = service.store.start_screening(saved["app"]["run_id"], "transcript")
        service.store.finish_screening(screening_id, {"specificity": []})
        router = AppRouter(service)
        try:
            status, _, body = router.dispatch(
                "POST", f"/api/screenings/{screening_id}/resume")
            assert status == 409
            assert "already finished" in json.loads(body)["detail"]
        finally:
            router.close()

    def test_restarting_the_app_does_not_interrupt_a_search_held_at_ncbi(
        self, tmp_path: Path
    ):
        # A 'running' screening is this process's work and is interrupted. A 'waiting' one is
        # NCBI's, under an id that outlives the process, so it must survive a restart.
        router, _run_id, screening_id = self._waiting(tmp_path)
        router.close()
        reopened = RunStore(tmp_path / "runs.sqlite3")
        held = reopened.screening(screening_id)
        assert held["status"] == "waiting"
        assert held["request_id"] == "RID-ABC-1"
