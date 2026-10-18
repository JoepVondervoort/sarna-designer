"""Independent local-search oracle, native integration, and resumable reference setup."""

import gzip
import hashlib
import io
import json
import os
import platform
import shutil
import sqlite3
import subprocess
import zipfile
from collections.abc import Callable
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace

import httpx
import pytest

from sarna.offtarget import references
from sarna.offtarget.local import parse_alignment, reverse_complement, run_local_screening
from sarna.offtarget.references import ReferenceManager, installation_lock, reference_bases

QUERY = "GTGAAGCCATAGACTCGTA"


def test_footprints_validate_actual_bases_on_both_strands():
    subject = QUERY[:7] + "A" + QUERY[8:]
    for strand, reference in (("+", subject), ("-", reverse_complement(subject))):
        hit = parse_alignment(
            f"q0\t{strand}\tNC_000001.11\t0\tunused\tI\t0\t\n",
            QUERY,
            "title",
            0,
            io.BytesIO(reference.encode()),
        )
        assert hit.query_sequence == QUERY
        assert hit.subject_sequence == subject
        assert hit.mismatch_positions == (8,)
        assert hit.evalue is None
        assert hit.bit_score is None
        assert (hit.subject_start, hit.subject_end) == ((1, 19) if strand == "+" else (19, 1))
        json.dumps(hit.to_dict(), allow_nan=False)


def test_reference_bases_cross_line_boundary_without_crossing_record():
    sequence = "ACGT" * 40
    fasta = b"".join(sequence[i : i + 60].encode() + b"\n" for i in range(0, len(sequence), 60))
    for start in (0, 42, 58, 60, 119, 141):
        assert reference_bases(io.BytesIO(fasta), 0, start, 19) == sequence[start : start + 19]
    with pytest.raises(ValueError):
        reference_bases(io.BytesIO(fasta + b">next\n"), 0, 158, 19)


def test_status_is_read_only_and_never_downloads(tmp_path: Path):
    root = tmp_path / "not-created"
    status = ReferenceManager(root).status()
    assert not status["ready"]
    assert status["catalog_release"] == references.RELEASE
    assert not root.exists()
    assert sum(a["bytes"] for a in status["download_plan"]) == status["download_bytes"]
    assert next(a for a in status["download_plan"] if a["key"] == "genomic")["bytes"] == 972898531


def test_only_one_process_can_prepare_references(tmp_path: Path):
    with (
        installation_lock(tmp_path),
        pytest.raises(ValueError, match="already running"),
        installation_lock(tmp_path),
    ):
        pass
    with installation_lock(tmp_path):
        pass


def download_client(
    monkeypatch: pytest.MonkeyPatch, handler: Callable[[httpx.Request], httpx.Response]
) -> None:
    real_client = httpx.Client
    monkeypatch.setattr(
        httpx,
        "Client",
        lambda **kwargs: real_client(**kwargs, transport=httpx.MockTransport(handler)),
    )


@pytest.mark.parametrize("supports_range", [True, False])
def test_interrupted_download_resumes_or_restarts_safely(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    supports_range: bool,
):
    content = b"known reference data" * 100
    expected = hashlib.sha256(content).hexdigest()
    manager = ReferenceManager(tmp_path)
    cache = tmp_path / "downloads"
    cache.mkdir()
    (cache / (expected + ".part")).write_bytes(content[:121])
    requests = []

    def reply(request: httpx.Request) -> httpx.Response:
        requests.append(request)
        assert request.headers["range"] == "bytes=121-"
        return httpx.Response(
            206 if supports_range else 200,
            content=content[121:] if supports_range else content,
            headers={"Content-Range": f"bytes 121-{len(content) - 1}/{len(content)}"},
        )

    download_client(monkeypatch, reply)
    downloaded = manager._download("https://example.test/data", expected, len(content), "sha256")
    assert downloaded.read_bytes() == content
    assert len(requests) == 1
    assert (
        manager._download("https://example.test/data", expected, len(content), "sha256")
        == downloaded
    )
    assert len(requests) == 1  # Cached, verified asset never downloaded again.


def test_corrupt_download_cannot_be_activated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    manager = ReferenceManager(tmp_path)
    download_client(monkeypatch, lambda _: httpx.Response(200, content=b"wrong"))
    expected = hashlib.sha256(b"right").hexdigest()
    with pytest.raises(ValueError, match="checksum failed"):
        manager._download("https://example.test/data", expected, 5, "sha256")
    assert not (tmp_path / "downloads" / expected).exists()
    assert not manager.status()["ready"]


def test_wrong_range_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    manager = ReferenceManager(tmp_path)
    (tmp_path / "downloads").mkdir()
    (tmp_path / "downloads" / "expected.part").write_bytes(b"ab")
    download_client(
        monkeypatch,
        lambda _: httpx.Response(206, content=b"cd", headers={"Content-Range": "bytes 0-1/4"}),
    )
    with pytest.raises(ValueError, match="unexpected byte range"):
        manager._download("https://example.test/data", "expected", 4, "sha256")


@pytest.fixture
def native_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[ReferenceManager, str]:
    """Exercise the complete installer with actual Bowtie and tiny, source-shaped references."""
    if platform.system() not in references.TOOLS:
        pytest.skip("Automatic provisioning is unsupported on this platform")
    aligner, builder = shutil.which("bowtie-align-s"), shutil.which("bowtie-build-s")
    if not aligner or not builder:
        pytest.skip("Native integration requires Bowtie 1 binaries on PATH")
    variants = [QUERY]
    for count in (1, 2, 3, 4):
        sequence = list(QUERY)
        for position in (2, 7, 12, 17)[:count]:
            sequence[position] = next(base for base in "ACGT" if base != sequence[position])
        variants.append("".join(sequence))
    # Deliberately distributed mismatches; a seven-base exact-seed heuristic can miss these.
    genome = ("N" * 30).join([*variants, reverse_complement(QUERY), QUERY])
    source = tmp_path / "source"
    source.mkdir()
    contents = {
        "genomic": gzip.compress(f">NC_000001.11 chromosome 1\n{genome}\n".encode()),
        "transcript": gzip.compress(
            f">NM_000001.1 Homo sapiens example (GENE1), mRNA\n{genome}\n".encode()
        ),
        "assembly": (
            b"1\tassembled-molecule\t1\tChromosome\tCM1\t=\tNC_000001.11\t"
            b"Primary Assembly\t1000\tchr1\n"
        ),
        "annotation": gzip.compress(
            b"NC_000001.11\tRefSeq\tmRNA\t1\t1000\t.\t+\t.\t"
            b"ID=rna-NM_000001.1;transcript_id=NM_000001.1;gene=GENE1\n"
            b"NC_000001.11\tRefSeq\tsnoRNA\t500\t600\t.\t-\t.\t"
            b"ID=rna-NR_000002.1;transcript_id=NR_000002.1;gene=SMALLRNA\n"
            b"NC_000001.11\tRefSeq\tmRNA\t300\t550\t.\t-\t.\t"
            b"ID=rna-NM_000003.1;transcript_id=NM_000003.1;gene=MINUSGENE\n"
            b"NC_000001.11\tRefSeq\tmRNA\t300\t600\t.\t+\t.\t"
            b"ID=rna-NM_000004.1;transcript_id=NM_000004.1;gene=PARTIAL;partial=true\n"
        ),
    }
    assets = {}
    downloads = {}
    for key, content in contents.items():
        suffix = references.ASSETS[key][0]
        assets[key] = (suffix, hashlib.md5(content).hexdigest(), len(content))
        path = source / key
        path.write_bytes(content)
        downloads[references.SOURCE + references.PREFIX + suffix] = path
    archive = source / "tools.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipped:
        extension = ".exe" if os.name == "nt" else ""
        zipped.write(aligner, "bowtie/bin/bowtie-align-s" + extension)
        zipped.write(builder, "bowtie/bin/bowtie-build-s" + extension)
    tool_url = "https://example.test/bowtie.zip"
    downloads[tool_url] = archive
    monkeypatch.setattr(references, "ASSETS", assets)
    monkeypatch.setitem(
        references.TOOLS,
        platform.system(),
        (
            references.TOOLS[platform.system()][0],
            tool_url,
            hashlib.sha256(archive.read_bytes()).hexdigest(),
            archive.stat().st_size,
        ),
    )
    manager = ReferenceManager(tmp_path / "data with spaces")
    manager.root.mkdir()
    monkeypatch.setattr(manager, "_download", lambda url, *_: downloads[url])
    manager._prepare()
    assert manager.status()["ready"]
    folder, manifest = manager.active()
    assert manifest["release"] == references.RELEASE
    with closing(sqlite3.connect(folder / "annotation.sqlite3")) as db:
        assert db.execute(
            "SELECT tss,is_tss_proxy FROM transcripts WHERE accession='NR_000002.1'"
        ).fetchone() == (600, 0)
    assert all(
        (folder / (scope + suffix)).exists()
        for scope in ("genomic", "transcript")
        for suffix in references.INDEX_SUFFIXES
    )
    return manager, genome


def oracle(genome: str) -> set[tuple[int, str, str]]:
    """Exhaustively check every callable position independently of Bowtie and its output."""
    matches = set()
    for position in range(len(genome) - len(QUERY) + 1):
        sequence = genome[position : position + len(QUERY)]
        if set(sequence) - set("ACGT"):
            continue
        for strand, oriented in (("plus", sequence), ("minus", reverse_complement(sequence))):
            if sum(a != b for a, b in zip(QUERY, oriented, strict=True)) <= 3:
                matches.add((position, strand, oriented))
    return matches


@pytest.mark.strand
@pytest.mark.parametrize("scope", ["transcript", "genomic"])
def test_native_search_equals_independent_exhaustive_oracle(
    native_reference: tuple[ReferenceManager, str],
    scope: str,
    monkeypatch: pytest.MonkeyPatch,
):
    manager, genome = native_reference
    candidates = [
        {"id": "candidate-one", "rank": 1, "target_site_dna": QUERY},
        {"id": "duplicate-sequence", "rank": 2, "target_site_dna": QUERY},
    ]
    result = run_local_screening(
        manager,
        "run",
        candidates,
        scope,
        19,
        lambda *_: None,
        only=[str(c["id"]) for c in candidates],
        limit=1,
        window=(-100, 100),
    )
    expected = oracle(genome)
    for identifier in ("candidate-one", "duplicate-sequence"):
        actual = {
            (
                min(h["subject_start"], h["subject_end"]) - 1,
                h["subject_strand"],
                h["subject_sequence"],
            )
            for h in result["hits"]
            if h["candidate"] == identifier
        }
        assert actual == expected
    assert len(expected) == 6
    assert result["specificity"][0]["perfect_hits"] == 3
    assert result["specificity"][0]["hit_count"] == 6
    assert result["specificity"][0]["enumeration_complete"]
    assert not result["capped_candidates"]
    assert result["specificity"][0]["best_evalue"] is None
    json.dumps(result, allow_nan=False)
    if scope == "transcript":
        assert result["gene_annotations"]["NM_000001.1"]["chromosome"] == "1"
    else:
        assert len(result["proximities"]) == 12
        assert result["proximities"][0]["nearest_gene"] == "GENE1"
        minus = [p for p in result["proximities"] if p["nearest_gene"] == "MINUSGENE"]
        assert minus
        assert all(p["distance_to_tss"] == 550 - p["position"] for p in minus)
        assert not any(p["nearest_gene"] in {"SMALLRNA", "PARTIAL"} for p in result["proximities"])
    # Repeated batches must use cached native output, without invoking the native tool.
    monkeypatch.setattr("sarna.offtarget.local.software_identity", lambda: {})
    monkeypatch.setattr(
        subprocess, "Popen", lambda *_a, **_k: pytest.fail("Cached query ran again")
    )
    again = run_local_screening(manager, "run-again", candidates[:1], scope, 19, lambda *_: None)
    assert again["specificity"][0]["hit_count"] == 6
    assert again["specificity"][0]["alignment_cache_hit"]


def test_own_locus_removed_and_missing_reference_never_looks_clean(
    native_reference: tuple[ReferenceManager, str],
):
    manager, _ = native_reference
    candidate = {
        "id": "own",
        "rank": 1,
        "target_site_dna": QUERY,
        "genomic": {"contig": "1", "start": 0, "end": 19},
    }
    result = run_local_screening(manager, "run", [candidate], "genomic", 19, lambda *_: None)
    assert result["specificity"][0]["hit_count"] == 5
    assert result["specificity"][0]["on_target_hits"] == 1
    assert len(result["on_target_hits"]) == 1
    folder, _ = manager.active()
    (folder / "genomic.1.ebwt").unlink()
    with pytest.raises(ValueError, match="missing or changed"):
        run_local_screening(manager, "run", [candidate], "genomic", 19, lambda *_: None)


def test_display_limit_does_not_limit_search_counts(
    native_reference: tuple[ReferenceManager, str],
    monkeypatch: pytest.MonkeyPatch,
):
    from sarna.offtarget import local

    manager, _ = native_reference
    monkeypatch.setattr(local, "DETAIL_LIMIT", 2)
    result = run_local_screening(
        manager,
        "run",
        [{"id": "one", "rank": 1, "target_site_dna": QUERY}],
        "transcript",
        19,
        lambda *_: None,
    )
    assert result["specificity"][0]["hit_count"] == 6
    assert result["specificity"][0]["details_retained"] == 2
    row = result["specificity"][0]
    assert sum(row["identity_counts"].values()) == 6
    assert row["identity_counts"]["19"] == 3
    assert row["selection_assessment"]["status"] == "excluded"
    assert len(result["hits"]) == 2
    assert not result["specificity"][0]["alignments_capped"]


def test_resume_ready_installation_does_not_redownload(
    native_reference: tuple[ReferenceManager, str],
    monkeypatch: pytest.MonkeyPatch,
):
    manager, _ = native_reference
    monkeypatch.setattr(manager, "_download", lambda *_: pytest.fail("Ready references downloaded"))
    manager.start()
    assert manager.status()["ready"]
    assert manager._thread is None


def test_reference_routes_are_explicit_and_local_failure_has_no_remote_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    from sarna.app.web import AppRouter
    from tests.unit.test_app_screening import completed_run

    service, run_id, _ = completed_run(tmp_path)
    router = AppRouter(service)
    monkeypatch.setattr(service, "blast_factory", lambda: pytest.fail("Implicit remote request"))
    started = []
    monkeypatch.setattr(service.references, "start", lambda: started.append(True))
    try:
        status, _, body = router.dispatch("GET", "/api/references")
        assert status == 200 and not json.loads(body)["ready"]
        assert not started
        status, _, _ = router.dispatch("POST", "/api/references/setup", b"{}")
        assert status == 202 and started == [True]
        status, _, body = router.dispatch(
            "POST",
            f"/api/runs/{run_id}/screenings",
            json.dumps({"scope": "genomic", "backend": "local"}).encode(),
        )
        assert status == 422
        assert "Set up" in json.loads(body)["detail"]
        assert not service.store.screenings_for(run_id)
        remote, _ = service.store.reserve_screening(run_id, "genomic", ["one"], backend="ncbi")
        local, reused = service.store.reserve_screening(run_id, "genomic", ["one"], backend="local")
        assert not reused and local != remote
        first, _ = service.store.reserve_screening(
            run_id, "genomic", ["one"], backend="local", reference_release="R1"
        )
        second, reused = service.store.reserve_screening(
            run_id, "genomic", ["one"], backend="local", reference_release="R2"
        )
        assert not reused and first != second
    finally:
        router.close()


def test_pause_leaves_partial_download_and_resume_available(tmp_path: Path):
    manager = ReferenceManager(tmp_path)
    manager.pause()
    with pytest.raises(InterruptedError, match="paused"):
        manager._check()
    assert not manager.status()["ready"]


def test_native_service_persists_local_provenance_and_exports_without_remote_calls(
    native_reference: tuple[ReferenceManager, str],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    from sarna.app.web import AppRouter
    from tests.unit.test_app_screening import completed_run

    manager, _ = native_reference
    service, run_id, before = completed_run(tmp_path / "app")
    service.references = manager
    monkeypatch.setattr(service, "blast_factory", lambda: pytest.fail("Local search used BLAST"))
    candidate = next(c for c in before["candidates"] if c["rank"] is not None)
    screening = service.store.start_screening(run_id, "transcript")
    service.screen(screening, run_id, "transcript", only=[candidate["id"]], backend="local")
    saved = service.store.screening(screening)
    assert saved["status"] == "complete", saved
    assert saved["result"]["backend"] == "local"
    assert saved["result"]["reference"]["release"] == references.RELEASE
    assert saved["result"]["specificity"][0]["enumeration_complete"]
    assert service.store.snapshot(run_id) == before
    router = AppRouter(service)
    try:
        for extension in ("json", "html", "xlsx", "order.csv"):
            status, _, body = router.dispatch("GET", f"/api/runs/{run_id}/export.{extension}")
            assert status == 200, (extension, body[:200])
        annotations = service.screening_annotations(
            screening, [h["accession"] for h in saved["result"]["hits"]][:1]
        )
        assert isinstance(annotations, dict)
    finally:
        router.close()


def test_explicit_api_selection_is_not_cut_to_default_fifteen(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    import threading

    from sarna.app.web import AppRouter
    from tests.unit.test_app_screening import completed_run

    service, run_id, snapshot = completed_run(tmp_path)
    candidates = [
        {"id": f"candidate-{i}", "rank": i + 1, "target_site_dna": QUERY} for i in range(25)
    ]
    monkeypatch.setattr(service.store, "snapshot", lambda _: {**snapshot, "candidates": candidates})
    called = threading.Event()
    sizes = []

    def remember(*args: object, **_kwargs: object) -> None:
        sizes.append(args[4])
        called.set()

    monkeypatch.setattr(service, "screen", remember)
    router = AppRouter(service)
    try:
        status, _, body = router.dispatch(
            "POST",
            f"/api/runs/{run_id}/screenings",
            json.dumps({"backend": "ncbi", "candidates": [c["id"] for c in candidates]}).encode(),
        )
        assert status == 202
        assert called.wait(2)
        assert sizes == [25]
        record = service.store.screening(json.loads(body)["id"])
        assert len(record["result"]["candidates_searched"]) == 25
    finally:
        router.close()


def test_tool_self_check_fails_before_large_downloads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    if platform.system() not in references.TOOLS:
        pytest.skip("Automatic provisioning is unsupported on this platform")
    extension = ".exe" if os.name == "nt" else ""
    archive = tmp_path / "tools.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("tools/bowtie-align-s" + extension, b"test")
        zipped.writestr("tools/bowtie-build-s" + extension, b"test")
    manager = ReferenceManager(tmp_path / "references")
    manager.root.mkdir()
    requested = []

    def download(url: str, *_args: object) -> Path:
        requested.append(url)
        return archive

    def refuse(*_args: object) -> None:
        raise ValueError("Native self-check failed")

    monkeypatch.setattr(shutil, "disk_usage",
                        lambda _path: SimpleNamespace(free=25 * 1024**3))
    monkeypatch.setattr(manager, "_download", download)
    monkeypatch.setattr(manager, "command", lambda *_args: None)
    monkeypatch.setattr(manager, "_tool_check", refuse)
    with pytest.raises(ValueError, match="self-check failed"):
        manager._prepare()
    assert requested == [references.TOOLS[platform.system()][1]]
    assert not manager.status()["ready"]


def test_repetitive_search_fails_without_a_partial_cache_or_false_zero(
    native_reference: tuple[ReferenceManager, str],
    monkeypatch: pytest.MonkeyPatch,
):
    from sarna.offtarget import local

    manager, _ = native_reference
    monkeypatch.setattr(local, "SAFETY_HIT_LIMIT", 2)
    with pytest.raises(ValueError, match="No truncated result"):
        run_local_screening(
            manager,
            "run",
            [{"id": "one", "rank": 1, "target_site_dna": QUERY}],
            "transcript",
            19,
            lambda *_: None,
        )
    assert not list((manager.root / "search-cache").glob("*.gz"))


def test_missing_index_can_be_repaired_without_deleting_old_data(
    native_reference: tuple[ReferenceManager, str],
):
    manager, _ = native_reference
    folder, _ = manager.active()
    (folder / "genomic.1.ebwt").unlink()
    assert not manager.status()["ready"]
    manager._prepare()
    assert manager.status()["ready"]
    assert len(list(manager.root.glob(references.RELEASE + ".replaced-*"))) == 1


def test_queued_local_search_stays_bound_to_accepted_reference(
    native_reference: tuple[ReferenceManager, str],
):
    manager, _ = native_reference
    _, manifest = manager.active()
    release = manifest["release"]
    (manager.root / "active.json").write_text(
        json.dumps({"release": "future-release"}), encoding="utf-8"
    )
    result = run_local_screening(
        manager,
        "run",
        [{"id": "one", "rank": 1, "target_site_dna": QUERY}],
        "transcript",
        19,
        lambda *_: None,
        reference_release=release,
    )
    assert result["reference"]["release"] == release
    assert result["specificity"][0]["hit_count"] == 6


@pytest.mark.parametrize("compatible", [True, False])
def test_staged_index_checkpoint_reused_only_for_matching_inputs_and_options(
    native_reference: tuple[ReferenceManager, str],
    monkeypatch: pytest.MonkeyPatch,
    compatible: bool,
):
    manager, _ = native_reference
    folder, _ = manager.active()
    stage = manager.root / (references.RELEASE + ".preparing")
    folder.rename(stage)
    for suffix in references.INDEX_SUFFIXES:
        (stage / ("transcript" + suffix)).unlink()
    (stage / "transcript.index-ready").unlink()
    if not compatible:
        marker = stage / "genomic.index-ready"
        saved = json.loads(marker.read_text(encoding="utf-8"))
        saved["source_md5"] = "another-reference"
        marker.write_text(json.dumps(saved), encoding="utf-8")
    builds = {"genomic": 0, "transcript": 0}
    original = manager.command

    def track(args: list[str], log: Path) -> None:
        for scope in builds:
            if Path(args[-1]) == stage / scope and "bowtie-build-s" in Path(args[0]).name:
                builds[scope] += 1
        original(args, log)

    monkeypatch.setattr(manager, "command", track)
    manager._prepare()
    assert manager.status()["ready"]
    assert builds == {"genomic": 0 if compatible else 1, "transcript": 1}
