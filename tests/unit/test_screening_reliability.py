"""Regressions for incomplete searches, unknown annotation, pacing and saved-job recovery."""

import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import httpx
import pytest

from sarna.app.store import RunStore
from sarna.offtarget.models import GENOMIC_SCOPE, TRANSCRIPT_SCOPE
from sarna.offtarget.ncbi import BlastError, RequestScheduler, parse_report
from sarna.offtarget.proximity import place_hit
from sarna.offtarget.screening import run_screening, specificity_table
from tests.unit.test_app_screening import completed_run, screen_now
from tests.unit.test_offtarget_ncbi import PUT_REPLY, READY_WITH_HITS, REPORT, client_for
from tests.unit.test_offtarget_screening import CANDIDATES, GENOMIC_CANDIDATES, GENOMIC_REPORT


def test_scheduler_spaces_concurrent_clients_and_repeated_rid() -> None:
    current = [0.0]
    def advance(delay: float) -> None:
        current[0] += delay

    scheduler = RequestScheduler(lambda: current[0], advance)
    scheduler.wait({"CMD": "Put"})
    scheduler.wait({"CMD": "Get", "RID": "A", "FORMAT_OBJECT": "SearchInfo"})
    assert current[0] == 10
    scheduler.wait({"CMD": "Get", "RID": "A", "FORMAT_TYPE": "XML"})
    assert current[0] == 20
    scheduler.wait({"CMD": "Get", "RID": "B", "FORMAT_OBJECT": "SearchInfo"})
    assert current[0] == 30
    scheduler.wait({"CMD": "Get", "RID": "A", "FORMAT_OBJECT": "SearchInfo"})
    assert current[0] == 70


@pytest.mark.parametrize("raw", [b"<unexpected/>", b"<BlastOutput/>",
                                 b"<BlastOutput><BlastOutput_iterations/></BlastOutput>"])
def test_unexpected_valid_xml_cannot_become_zero_hits(raw: bytes) -> None:
    with pytest.raises(BlastError):
        parse_report(raw)


def test_missing_query_cannot_become_a_zero_hit_candidate() -> None:
    incomplete = REPORT.replace(b"saRNA-0002 some trailing description", b"different-query")
    client, _ = client_for([PUT_REPLY, READY_WITH_HITS, incomplete])
    with pytest.raises(BlastError, match="coverage mismatch"):
        run_screening("run", CANDIDATES, TRANSCRIPT_SCOPE, 19, client, sleep=lambda _: None)


def test_server_sequence_limit_marks_counts_as_lower_bounds() -> None:
    parsed = parse_report(REPORT, max_sequences=1)
    assert "saRNA-1549" in parsed.capped
    assert "saRNA-0002" not in parsed.capped


def test_footprint_retains_subject_bases_and_leading_mismatch() -> None:
    detailed = REPORT.replace(b"<Hsp_identity>19</Hsp_identity>",
        b"<Hsp_qseq>ACGTACGTACGTACGTACG</Hsp_qseq>"
        b"<Hsp_hseq>TCGTACGTACGTACGTACG</Hsp_hseq>"
        b"<Hsp_midline> ||||||||||||||||||</Hsp_midline>"
        b"<Hsp_identity>18</Hsp_identity>", 1)
    hit = parse_report(detailed).hits[0]
    assert hit.mismatch_positions == (1,)
    assert hit.to_dict()["subject_sequence"].startswith("T")
    assert hit.to_dict()["query_sequence"].startswith("A")


def test_lookup_failure_is_unknown_and_does_not_report_zero_window_hits() -> None:
    client, _ = client_for([PUT_REPLY, READY_WITH_HITS, GENOMIC_REPORT])

    def unavailable(_interval: Any) -> list[dict[str, Any]]:  # noqa: ANN401
        raise RuntimeError("Ensembl unavailable")

    record = run_screening("run", GENOMIC_CANDIDATES, GENOMIC_SCOPE, 19, client,
                           transcripts=unavailable, sleep=lambda _: None)
    assert record.proximities[0].annotation_status == "unavailable"
    row = specificity_table(record, 19)[0]
    assert row["in_run_window_hits"] is None
    assert row["annotation_unknown_hits"] == 1
    assert any("unknown, not absent" in note for note in record.notes)


def test_long_overlapping_gene_does_not_claim_a_tss_within_ten_kb() -> None:
    hit = parse_report(GENOMIC_REPORT).hits[-1]
    placement = place_hit(hit, lambda _: [{"id": "long", "start": 80000,
                                           "end": 100000, "strand": -1}], None)
    assert placement.distance_to_tss is None


def test_submitted_request_and_candidates_survive_restart(tmp_path: Path) -> None:
    service, run_id, _ = completed_run(tmp_path)
    screening_id, reused = service.store.reserve_screening(run_id, TRANSCRIPT_SCOPE, ["saRNA-001"])
    assert not reused
    service.store.screening_progress(screening_id, "Accepted", progress={
        "phase": "waiting", "request_id": "RID-PRESERVED", "candidates_searched": ["saRNA-001"],
    })
    reopened = RunStore(tmp_path / "runs.sqlite3")
    # service_for uses tmp_path directly; reopening must keep the remote handle.
    held = reopened.screening(screening_id)
    assert held["status"] == "waiting"
    assert held["request_id"] == "RID-PRESERVED"
    assert held["result"]["candidates_searched"] == ["saRNA-001"]


def test_duplicate_active_submission_is_reserved_atomically(tmp_path: Path) -> None:
    service, run_id, _ = completed_run(tmp_path)
    with ThreadPoolExecutor(max_workers=4) as workers:
        reservations = list(workers.map(lambda _: service.store.reserve_screening(
            run_id, TRANSCRIPT_SCOPE, ["saRNA-001", "saRNA-002"]), range(8)))
    assert len({screening_id for screening_id, _ in reservations}) == 1
    assert sum(not reused for _, reused in reservations) == 1


def test_explicit_ready_without_hit_status_is_not_a_negative_result() -> None:
    client, _ = client_for([b"Status=READY"])
    with pytest.raises(BlastError, match="explicit hit status"):
        client.status("RID")


def test_transcript_chromosome_annotation_is_saved_with_its_source(tmp_path: Path) -> None:
    from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
    from sarna.genome.ensembl import Gene

    service, run_id, _ = completed_run(tmp_path)
    screening_id = screen_now(service, run_id)
    record = service.store.screening(screening_id)["result"]
    for hit in record["hits"]:
        hit["title"] = "Homo sapiens pipecolic acid oxidase (PIPOX), mRNA"
    service.store.finish_screening(screening_id, record)

    class Genes:
        def gene(self, symbol: str) -> Gene:
            assert symbol == "PIPOX"
            return Gene("ENSG00000000001", "PIPOX",
                        GenomicInterval(Assembly.GRCH38, "17", 100, 1000, Strand.PLUS), "oxidase")

    service.client_factory = lambda _: Genes()  # type: ignore[assignment,return-value]
    annotation = service.screening_annotations(screening_id, ["NM_006920"])["NM_006920"]
    assert annotation["chromosome"] == "17"
    assert annotation["source"] == "Ensembl gene-symbol lookup"
    saved = service.store.screening(screening_id)["result"]["gene_annotations"]
    assert saved["NM_006920"] == annotation


def test_temporary_retrieval_error_retries_but_submission_does_not() -> None:
    from sarna.offtarget.ncbi import BlastClient

    requests = []

    def handler(request: httpx.Request) -> httpx.Response:
        requests.append(request)
        return httpx.Response(503 if len(requests) == 1 else 200, content=b"Status=WAITING")

    client = BlastClient(client=httpx.Client(transport=httpx.MockTransport(handler)))
    assert client.status("RID") == ("WAITING", False)
    assert len(requests) == 2
    requests.clear()
    with pytest.raises(BlastError):
        client.submit({"one": "ACGT"}, TRANSCRIPT_SCOPE)
    assert len(requests) == 1


def test_browser_waiting_stops_batches_and_footprint_displays_real_bases() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is needed for frontend state regressions")
    source = Path("src/sarna/app/static/offtarget.js").resolve()
    script = r'''
const fs=require('fs'), vm=require('vm'), assert=require('assert');
const elements=new Map();
const context={console, setTimeout, clearTimeout, frozen:null,
  state:{checked:new Set(),result:null}, escapeHTML:s=>String(s),
  document:{addEventListener:()=>{}},
  $:id=>{
    if(!elements.has(id)) elements.set(id,{addEventListener:()=>{},textContent:'',innerHTML:''});
    return elements.get(id);
  },
  error:()=>{}, signed:String, fmt:String};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`
(async()=>{
  state.result={app:{run_id:'one'},candidates:Array.from({length:26},(_,i)=>(
    {id:'c'+i,rank:i+1,target_site_dna:'ACGT'}))};
  let submissions=0;
  api=async(path,options)=>{
    if(options){submissions++;return {id:'job'};}
    return {id:'job',scope:'transcript',status:'waiting',error:'Still at NCBI'};
  };
  await startScreening('transcript','all');
  if(submissions!==1) throw Error('Waiting job did not stop later batches');
  if(state.screeningJobs.transcript.status!=='waiting') throw Error('Waiting was overwritten');
  state.result.candidates=[{id:'c',rank:1,target_site_dna:'ACGT'}];
  const html=footprint({candidate:'c',query_start:1,query_end:4,query_sequence:'AC-GT',
    subject_sequence:'ATCGT',subject_strand:'plus'},4);
  if(!html.includes('fp-base mismatch')||!html.includes('fp-base gap')
    ||!html.includes('>T</span>')) throw Error('Actual bases or gaps missing');
  const old=footprint({candidate:'c',query_start:1,query_end:4,
    identities:3,subject_strand:'plus'},4);
  if(!old.includes('unavailable')||!old.includes('fp-base unknown'))
    throw Error('Legacy data falsely coloured as a match');
  if(phaseFraction({phase:'waiting',elapsed_seconds:10,estimate_seconds:100}).fraction!==null)
    throw Error('Estimate drawn as progress');
})()` ,context).then(()=>console.log('PASS')).catch(e=>{console.error(e);process.exitCode=1});
'''
    result = subprocess.run([node, "-e", script, str(source)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "PASS" in result.stdout


def test_documented_xml2_and_classic_xml_produce_the_same_alignment() -> None:
    import xml.etree.ElementTree as ET

    old = ET.fromstring(REPORT)
    root = ET.Element("BlastXML2", xmlns="http://www.ncbi.nlm.nih.gov")
    for iteration in old.iter("Iteration"):
        search = ET.SubElement(root, "Search")
        ET.SubElement(search, "query-title").text = iteration.findtext("Iteration_query-def")
        hits = ET.SubElement(search, "hits")
        for source in iteration.iter("Hit"):
            hit = ET.SubElement(hits, "Hit")
            description = ET.SubElement(ET.SubElement(hit, "description"), "HitDescr")
            ET.SubElement(description, "accession").text = source.findtext("Hit_accession")
            ET.SubElement(description, "title").text = source.findtext("Hit_def")
            hsps = ET.SubElement(hit, "hsps")
            for original in source.iter("Hsp"):
                hsp = ET.SubElement(hsps, "Hsp")
                for child in original:
                    ET.SubElement(hsp, child.tag.removeprefix("Hsp_")).text = child.text
    assert parse_report(ET.tostring(root)).hits == parse_report(REPORT).hits


def test_annotation_budget_does_not_let_the_first_query_starve_the_next(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from sarna.offtarget import screening
    from tests.unit.test_offtarget_ncbi import multi_query_report

    monkeypatch.setattr(screening, "MAX_PLACEMENTS", 2)
    report = multi_query_report({"saRNA-1549": 5, "saRNA-0002": 1})
    client, _ = client_for([PUT_REPLY, READY_WITH_HITS, report])
    record = run_screening("run", CANDIDATES, GENOMIC_SCOPE, 19, client,
                           transcripts=lambda _: [], sleep=lambda _: None)
    assert {p.hit.candidate_id for p in record.proximities} == {"saRNA-1549", "saRNA-0002"}
    assert len(record.hits) == 6


def test_failed_annotation_block_is_not_repeated_for_every_hit() -> None:
    from sarna.offtarget.proximity import BlockedTranscripts

    calls = []

    def unavailable(interval: Any) -> list[dict[str, Any]]:  # noqa: ANN401
        calls.append(interval)
        raise RuntimeError("Service unavailable")

    lookup = BlockedTranscripts(unavailable)
    hits = parse_report(GENOMIC_REPORT).hits
    first = place_hit(hits[0], lookup, None)
    second = place_hit(hits[0], lookup, None)
    assert first.annotation_status == second.annotation_status == "unavailable"
    assert len(calls) == 1


def test_report_download_requests_single_file_xml2_not_the_zip_format() -> None:
    from tests.unit.test_offtarget_ncbi import form

    client, requests = client_for([READY_WITH_HITS, REPORT])
    from sarna.offtarget.ncbi import parameters_for

    report = client.fetch("EXISTING-RID", parameters_for(TRANSCRIPT_SCOPE))
    assert len(report.hits) == 2
    assert form(requests[-1])["FORMAT_TYPE"] == "XML2_S"
    assert all(form(request)["CMD"] == "Get" for request in requests)


def test_unreadable_report_download_is_retried_without_resubmission() -> None:
    from sarna.offtarget.ncbi import parameters_for
    from tests.unit.test_offtarget_ncbi import form

    client, requests = client_for([READY_WITH_HITS, b"<broken", REPORT])
    report = client.fetch("EXISTING-RID", parameters_for(TRANSCRIPT_SCOPE))
    assert len(report.hits) == 2
    assert len(requests) == 3
    assert all(form(request)["RID"] == "EXISTING-RID" for request in requests)
    assert all(form(request)["CMD"] == "Get" for request in requests)


def test_zip_archive_is_explained_instead_of_reporting_an_xml_syntax_error() -> None:
    from sarna.offtarget.ncbi import BlastReportError

    with pytest.raises(BlastReportError, match=r"ZIP archive.*single-file XML"):
        parse_report(b"PK\x03\x04\x14\x00\x08\x08")


def test_unreadable_completed_search_stays_collectable_and_can_be_resumed(tmp_path: Path) -> None:
    from tests.unit.test_app_screening import REPORT as APP_REPORT
    from tests.unit.test_offtarget_ncbi import form

    service, run_id, saved = completed_run(tmp_path)
    ids = [c["id"] for c in saved["candidates"] if c["rank"] is not None]
    screening_id, _ = service.store.reserve_screening(run_id, TRANSCRIPT_SCOPE, ids)
    client, requests = client_for([PUT_REPLY, READY_WITH_HITS, b"<broken"])
    service.screen(screening_id, run_id, TRANSCRIPT_SCOPE,
                   blast_client=client, sleep=lambda _: None)
    held = service.store.screening(screening_id)
    assert held["status"] == "waiting"
    assert held["stage"] == "Report retrieval needs retry"
    assert held["request_id"] == "ABC123XYZ014"
    assert held["result"]["candidates_searched"] == ids
    assert "specificity" not in held["result"]
    assert sum(form(request)["CMD"] == "Put" for request in requests) == 1
    assert service.store.resume_screening(screening_id)
    client, requests = client_for([READY_WITH_HITS, APP_REPORT])
    service.resume(screening_id, blast_client=client, sleep=lambda _: None)
    assert service.store.screening(screening_id)["status"] == "complete"
    assert all(form(request)["CMD"] == "Get" for request in requests)


@pytest.mark.parametrize("reason, expected", [
    ("could not parse NCBI's XML reply: "
     "not well-formed (invalid token): line 1, column 2", "waiting"),
    ("BLAST query coverage mismatch", "error"),
])
def test_restart_recovers_old_xml_download_failures_only(
    tmp_path: Path, reason: str, expected: str,
) -> None:
    service, run_id, _ = completed_run(tmp_path)
    screening_id, _ = service.store.reserve_screening(run_id, TRANSCRIPT_SCOPE, ["saRNA-001"])
    service.store.screening_progress(screening_id, "Accepted", progress={
        "request_id": "EXISTING-RID", "candidates_searched": ["saRNA-001"],
    })
    service.store.fail_screening(screening_id, reason)
    reopened = RunStore(tmp_path / "runs.sqlite3")
    held = reopened.screening(screening_id)
    assert held["status"] == expected
    assert held["request_id"] == "EXISTING-RID"
    assert held["result"]["candidates_searched"] == ["saRNA-001"]


def test_xml2_hit_strand_is_respected_without_a_legacy_hit_frame() -> None:
    raw = b'''<BlastXML2><Search><query-title>candidate-1</query-title><hits><Hit>
    <description><HitDescr><accession>NM_1</accession><title>Example</title></HitDescr></description>
    <hsps><Hsp><identity>4</identity><align-len>4</align-len><query-from>1</query-from>
    <query-to>4</query-to><hit-from>20</hit-from><hit-to>17</hit-to><hit-strand>Minus</hit-strand>
    <qseq>ACGT</qseq><hseq>ACGT</hseq><midline>||||</midline></Hsp></hsps>
    </Hit></hits></Search></BlastXML2>'''
    assert parse_report(raw).hits[0].subject_strand == "minus"
