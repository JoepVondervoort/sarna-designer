"""Exercise real local jobs and exports using recorded sources, without network or sockets."""

import csv
import io
import json
import sqlite3
import zlib
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from sarna.app.models import DesignRequest
from sarna.app.service import DesignService
from sarna.app.store import RunStore
from sarna.app.web import AppRouter
from sarna.genome.ensembl import EnsemblCacheMissError, EnsemblClient
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

# Known 19-nt target from the cached SCN1A run; deliberately CpG-depleted flanks.
DNA = "GCTCTACTCTCAGTTCCAAACTATGACCTACTATGATACCATAC"


def service_for(tmp_path: Path, routes: dict[str, Any] | None = None) -> DesignService:
    client = client_for(routes or {
        GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
        "/sequence/region/": {"seq": AUTO_SEQUENCE}, "/overlap/region/": [],
    })
    return DesignService(tmp_path, client_factory=lambda _: client)


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
        assert json.loads(data) == saved
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
        assert connection.execute("PRAGMA user_version").fetchone()[0] == 3


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
