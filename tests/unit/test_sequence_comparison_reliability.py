"""Comparison reads are frozen, writes are compact and conflicts cannot erase inputs."""

import io
import json
import zlib
from base64 import b64encode
from pathlib import Path
from typing import Any
from unittest.mock import patch
from zipfile import ZipFile

from sarna.app.sequence_import import compare, file_rows, input_identity
from sarna.app.store import RunStore
from sarna.app.web import AppRouter
from tests.unit.test_app import run_request, service_for


def inputs(name: str = "Existing") -> dict[str, str]:
    return {"name": name, "sequence": "AAAAACCCCCGGGGGUUUUU", "role": "sense"}


def test_reopening_and_exporting_never_recalculate_saved_analysis(tmp_path: Path):
    service = service_for(tmp_path)
    run = run_request(service)
    run_id = run["app"]["run_id"]
    router = AppRouter(service)
    try:
        status, _, body = router.dispatch(
            "POST",
            f"/api/runs/{run_id}/sequence-comparisons",
            json.dumps({"revision": 0, "records": [inputs()]}).encode(),
        )
        assert status == 200
        saved = json.loads(body)["records"][0]
        assert saved["analysis_provenance"]["package_sha256"]
        with patch(
            "sarna.app.sequence_import.compare_asymmetry",
            side_effect=AssertionError(
                "A read or ordinary save must not calculate historical predictions",
            ),
        ):
            assert RunStore(service.store.path).snapshot(run_id)["sequence_comparisons"] == [saved]
            status, _, body = router.dispatch("GET", f"/api/runs/{run_id}/export.json")
            assert status == 200
            assert json.loads(body)["sequence_comparisons"] == [saved]
            status, _, body = router.dispatch(
                "POST",
                f"/api/runs/{run_id}/sequence-comparisons",
                json.dumps({"revision": 1, "records": [inputs()]}).encode(),
            )
            assert status == 200
            assert json.loads(body)["records"] == [saved]
        before = service.store.snapshot(run_id)
        with patch("sarna.app.sequence_import.compare", wraps=compare) as calculation:
            status, _, body = router.dispatch(
                "POST",
                f"/api/runs/{run_id}/sequence-comparisons/recalculate",
                b'{"revision":2}',
            )
            assert status == 200
            assert calculation.called
        after = service.store.snapshot(run_id)
        assert after["sequence_comparison_revision"] == 3
        assert (
            before["sequence_comparisons"][0]["analysis_provenance"]["analyzed_at"]
            != (after["sequence_comparisons"][0]["analysis_provenance"]["analyzed_at"])
        )
        with service.store.connect() as connection:
            history = connection.execute(
                "SELECT records FROM sequence_comparison_history WHERE run_id=? AND revision=2",
                (run_id,),
            ).fetchone()
        assert json.loads(history["records"]) == before["sequence_comparisons"]
    finally:
        router.close()


def test_two_windows_cannot_replace_each_others_records(tmp_path: Path):
    service = service_for(tmp_path)
    run_id = run_request(service)["app"]["run_id"]
    router = AppRouter(service)
    path = f"/api/runs/{run_id}/sequence-comparisons"
    try:
        for name, expected in [("window A", 200), ("window B", 409)]:
            status, _, body = router.dispatch(
                "POST",
                path,
                json.dumps(
                    {
                        "revision": 0,
                        "records": [inputs(name)],
                    }
                ).encode(),
            )
            assert status == expected
        assert json.loads(body)["code"] == "sequence-comparison-conflict"
        status, _, body = router.dispatch("GET", path)
        current = json.loads(body)
        assert status == 200 and current["revision"] == 1
        assert [row["name"] for row in current["records"]] == ["window A"]
        status, _, _ = router.dispatch(
            "POST",
            path,
            json.dumps(
                {
                    "revision": 1,
                    "records": [inputs("window A"), inputs("window B")],
                }
            ).encode(),
        )
        assert status == 200
        assert [r["name"] for r in service.store.snapshot(run_id)["sequence_comparisons"]] == [
            "window A",
            "window B",
        ]
    finally:
        router.close()


def test_repeated_additions_use_inputs_and_keep_previous_analysis(tmp_path: Path):
    service = service_for(tmp_path)
    run = run_request(service)
    router = AppRouter(service)
    path = f"/api/runs/{run['app']['run_id']}/sequence-comparisons"
    rows = [
        {
            "name": f"sequence-{i}",
            "sequence": run["promoter"]["sequence"][i : i + 19],
            "role": "sense",
        }
        for i in range(20)
    ]
    try:
        status, _, body = router.dispatch(
            "POST",
            path,
            json.dumps(
                {
                    "revision": 0,
                    "records": rows,
                }
            ).encode(),
        )
        saved = json.loads(body)["records"]
        assert status == 200 and len(json.dumps(saved)) > 200_000
        # Browser input projection excludes all nested candidate and analysis results.
        rows.append(inputs("new input"))
        compact = json.dumps({"revision": 1, "records": rows}).encode()
        assert len(compact) < 10_000
        status, _, body = router.dispatch("POST", path, compact)
        assert status == 200
        assert json.loads(body)["records"][:20] == saved
        assert len(json.loads(body)["records"]) == 21
        assert input_identity(saved[0]) == input_identity(rows[0])
    finally:
        router.close()


def test_legacy_records_are_preserved_and_clearly_unversioned(tmp_path: Path):
    service = service_for(tmp_path)
    run = run_request(service)
    run_id = run["app"]["run_id"]
    legacy = compare(run, [inputs()])
    legacy[0].pop("analysis_provenance")
    with service.store.connect() as connection:
        raw = connection.execute("SELECT snapshot FROM runs WHERE id=?", (run_id,)).fetchone()
        original = json.loads(zlib.decompress(raw["snapshot"]))
        original.pop("sequence_comparisons")
        connection.execute(
            "UPDATE runs SET snapshot=? WHERE id=?",
            (zlib.compress(json.dumps(original).encode()), run_id),
        )
        connection.execute(
            "INSERT INTO sequence_comparisons(run_id,records) VALUES (?,?)",
            (run_id, json.dumps(legacy)),
        )
    with patch(
        "sarna.app.sequence_import.compare", side_effect=AssertionError("no read-time analysis")
    ):
        result = RunStore(service.store.path).snapshot(run_id)
    assert result["sequence_comparisons"] == legacy
    assert "Legacy comparison" in result["sequence_comparison_notice"]
    assert result["sequence_comparison_revision"] == 0


def test_excel_first_sheet_uses_workbook_order_and_relationship_target():
    archive = io.BytesIO()
    with ZipFile(archive, "w") as workbook:
        workbook.writestr(
            "xl/workbook.xml",
            """<workbook
xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets><sheet name="My sequences" sheetId="7" r:id="rId7"/></sheets></workbook>""",
        )
        workbook.writestr(
            "xl/_rels/workbook.xml.rels",
            """<Relationships
xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId7" Target="/xl/worksheets/sheet7.xml"/></Relationships>""",
        )
        workbook.writestr(
            "xl/worksheets/sheet7.xml",
            """<worksheet
xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row>
<c r="A1" t="inlineStr"><is><t>AUGC</t></is></c></row></sheetData></worksheet>""",
        )
        workbook.writestr("xl/worksheets/sheet1.xml", "This is not the selected worksheet")
    assert file_rows("sequences.xlsx", archive.getvalue()) == [{"sequence": "AUGC"}]


def test_large_csv_import_does_not_relax_other_api_body_limits(tmp_path: Path):
    """Exercise the real HTTP handler, rather than bypassing its transport limits."""
    from threading import Thread
    from urllib.error import HTTPError
    from urllib.request import Request, urlopen

    import pytest

    from sarna.app.web import make_server

    router = AppRouter(service_for(tmp_path))
    server = make_server(router, 0)
    worker = Thread(target=server.serve_forever, daemon=True)
    worker.start()
    try:
        csv = (
            "name,sequence,source\n" + "".join(f"row-{i},AUGC,{'x' * 1000}\n" for i in range(180))
        ).encode()
        body = json.dumps(
            {"filename": "sequences.csv", "content": b64encode(csv).decode()}
        ).encode()
        assert 200_000 < len(body) < 2_800_000
        request = Request(
            f"http://127.0.0.1:{server.server_port}/api/sequences/import",
            data=body,
            headers={"Content-Type": "application/json"},
        )
        with urlopen(request) as response:
            assert response.status == 200
            result = json.load(response)
        assert len(result["records"]) == 180 and not result["errors"]
        prefix = (
            "name,sequence,source\n"
            + "".join(f"row-{i},AUGC,{'x' * 9900}\n" for i in range(199))
            + "last,AUGC,"
        ).encode()
        exact_limit = prefix + b"x" * (2_000_000 - len(prefix))
        boundary_request = Request(
            f"http://127.0.0.1:{server.server_port}/api/sequences/import",
            data=json.dumps(
                {"filename": "limit.csv", "content": b64encode(exact_limit).decode()}
            ).encode(),
            headers={"Content-Type": "application/json"},
        )
        with urlopen(boundary_request) as response:
            assert response.status == 200
            assert len(json.load(response)["records"]) == 200
        request = Request(
            f"http://127.0.0.1:{server.server_port}/api/runs",
            data=body,
            headers={"Content-Type": "application/json"},
        )
        with pytest.raises(HTTPError) as error:
            urlopen(request)
        assert error.value.code == 413
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)
        router.close()


def test_selecting_reference_preserves_draft_and_does_not_submit():
    import shutil
    import subprocess

    import pytest

    if not shutil.which("node"):
        pytest.skip("Node needed for reference-selection behavior")
    source = Path("src/sarna/app/static/app.js").read_text()
    function = (
        "async function designAroundSaRNA("
        + source.split(
            "async function designAroundSaRNA(",
            1,
        )[1].split("document.addEventListener('click'", 1)[0]
    )
    program = (
        """
const frozen=null;
const reference={name:'SCN1A-10',known_location:{anchor_id:'scn1a-long-2008'}};
const state={busy:false,gene:{symbol:'SCN1A',sarna_references:[reference]},
 result:{gene:{symbol:'SCN1A'},app:{run_id:'unchanged'}},
 draft:{gc:{minimum:47}},customReferences:{mine:{name:'keep me'}}};
const nodes={reference:{options:[{value:'anchor:scn1a-long-2008'}],value:''},
 upstream:{value:'3500'},downstream:{value:'750'},'sequence-status':{},'job-status':{}};
const $=id=>nodes[id];let submitted=0,restored=0;
const design=async()=>{submitted++;}, restoreDraft=()=>{restored++;};
const chooseOverviewReference=key=>{nodes.reference.value=key;};
const setWorkspaceTab=()=>{},error=message=>{throw Error(message);};
"""
        + function
        + """
(async()=>{
 await designAroundSaRNA('SCN1A-10');
 if(submitted||restored)process.exit(1);
 if(nodes.upstream.value!=='3500'||nodes.downstream.value!=='750')process.exit(2);
 if(state.draft.gc.minimum!==47||!state.customReferences.mine)process.exit(3);
 if(nodes.reference.value!=='anchor:scn1a-long-2008')process.exit(4);
 if(state.result.app.run_id!=='unchanged')process.exit(5);
})().catch(error=>{console.error(error);process.exit(6);});
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)


def test_curated_control_results_are_frozen_inside_new_gene_runs(tmp_path: Path):
    service = service_for(tmp_path)
    run = run_request(
        service, source="gene", symbol="SCN1A", sequence="", upstream=25, downstream=5
    )
    controls = run["sequence_comparisons"]
    assert [record["name"] for record in controls] == ["SCN1A-10", "SCN1A-7"]
    assert all(record["analysis_provenance"]["package_sha256"] for record in controls)
    with patch(
        "sarna.app.sequence_import.compare", side_effect=AssertionError("no read-time analysis")
    ):
        assert service.store.snapshot(run["app"]["run_id"])["sequence_comparisons"] == controls


def test_simultaneous_saves_are_checked_inside_the_transaction(tmp_path: Path):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Barrier

    service = service_for(tmp_path)
    run_id = run_request(service)["app"]["run_id"]
    router = AppRouter(service)
    barrier = Barrier(2)

    def simultaneous_analysis(
        snapshot: dict[str, Any],
        rows: list[dict[str, Any]],
        *,
        refresh_models: bool = False,
    ) -> list[dict[str, Any]]:
        barrier.wait(timeout=10)
        return compare(snapshot, rows, refresh_models=refresh_models)

    def save(name: str) -> int:
        return router.dispatch(
            "POST",
            f"/api/runs/{run_id}/sequence-comparisons",
            json.dumps(
                {
                    "revision": 0,
                    "records": [inputs(name)],
                }
            ).encode(),
        )[0]

    try:
        with (
            patch("sarna.app.sequence_import.compare", side_effect=simultaneous_analysis),
            ThreadPoolExecutor(max_workers=2) as writers,
        ):
            statuses = list(writers.map(save, ["first writer", "second writer"]))
        assert sorted(statuses) == [200, 409]
        assert service.store.snapshot(run_id)["sequence_comparison_revision"] == 1
        assert len(service.store.snapshot(run_id)["sequence_comparisons"]) == 1
    finally:
        router.close()
