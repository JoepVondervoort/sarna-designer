"""Sequence identities must survive imports, complementation and repeated placements."""

import io
from pathlib import Path
from zipfile import ZipFile

import pytest

from sarna.app.sequence_import import compare, core, file_rows, normalize, references


def compact_input(record: dict[str, object]) -> dict[str, object]:
    return {
        key: record[key]
        for key in (
            "name",
            "sequence",
            "antisense",
            "role",
            "direction",
            "source",
        )
        if key in record
    }


def test_rna_conversion_is_not_complementation():
    assert core("AUGC") == ("ATGC", "")
    snapshot = {"promoter": {"sequence": "GCATTTTGCAT"}, "candidates": []}
    record = compare(snapshot, [{"sequence": "AUGC", "role": "unknown"}])[0]
    assert [m["offset"] for m in record["matches"]] == [0, 7]
    assert all(m["input_orientation"] == "reverse complement" for m in record["matches"])
    assert not compare(snapshot, [{"sequence": "AUGC", "role": "sense"}])[0]["matches"]


def test_duplex_chemistry_and_direction():
    record = normalize({"sense": "AUGC[dT][dT]", "antisense": "GCAU[UU]"})
    assert record["overhang"] == "[DT][DT]"
    assert record["partner_overhang"] == "[UU]"
    assert not record["partner_inferred"]
    assert normalize({"sequence": "CGUA", "direction": "3to5"})["dna"] == "ATGC"
    with pytest.raises(ValueError, match="reverse complements"):
        normalize({"sense": "AUGC", "antisense": "AAAA"})
    with pytest.raises(ValueError, match="Mixed"):
        core("AUTG")


def test_csv_txt_and_sparse_excel():
    assert file_rows("s.txt", b"AUGC\nGCAU")[1]["sequence"] == "GCAU"
    assert file_rows("s.csv", b"name,sequence,role\na,AUGC,antisense")[0]["role"] == "antisense"
    output = io.BytesIO()
    with ZipFile(output, "w") as archive:
        archive.writestr(
            "xl/workbook.xml",
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            '<sheets><sheet name="Sequences" sheetId="1" r:id="rId1"/></sheets></workbook>',
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            '<Relationships><Relationship Id="rId1" '
            'Target="worksheets/sheet1.xml"/></Relationships>',
        )
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            """<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row><c r="A1" t="inlineStr"><is><t>name</t></is></c><c r="C1" t="inlineStr"><is><t>sequence</t></is></c></row><row><c r="A2" t="inlineStr"><is><t>test</t></is></c><c r="C2" t="inlineStr"><is><t>AUGC</t></is></c></row></sheetData></worksheet>""",  # noqa: E501
        )
    assert file_rows("s.xlsx", output.getvalue())[0]["sequence"] == "AUGC"


def test_curated_references_preserve_researcher_evidence():
    records = references("SCN1A")
    assert [r["name"] for r in records] == ["SCN1A-10", "SCN1A-7"]
    assert "RNA-seq" in records[1]["source"]
    assert "not activating" in records[1]["source"]
    assert references("CDKN1A")[0]["sequence"] == "CCAACUCAUUCUCCAAGUA"
    assert not references("OTHER")


def test_import_save_reopen_and_exports(tmp_path: Path):
    import json
    from base64 import b64encode

    from sarna.app.web import AppRouter
    from tests.unit.test_app import run_request, service_for

    service = service_for(tmp_path)
    snapshot = run_request(service)
    run_id = snapshot["app"]["run_id"]
    router = AppRouter(service)
    try:
        status, _, data = router.dispatch(
            "POST",
            "/api/sequences/import",
            json.dumps(
                {
                    "filename": "test.csv",
                    "content": b64encode(
                        b"name,sequence,role\nKnown,GCTCTACTCTCAGTTCCAA,sense\nBad,XXX,sense"
                    ).decode(),
                }
            ).encode(),
        )
        imported = json.loads(data)
        assert status == 200 and len(imported["records"]) == 1 and imported["errors"]
        status, _, data = router.dispatch(
            "POST",
            f"/api/runs/{run_id}/sequence-comparisons",
            json.dumps(
                {"revision": 0, "records": [compact_input(r) for r in imported["records"]]}
            ).encode(),
        )
        assert status == 200
        assert json.loads(data)["records"][0]["matches"][0]["offset"] == 0
        assert service.store.snapshot(run_id)["sequence_comparisons"][0]["name"] == "Known"
        status, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.json")
        assert status == 200
        assert json.loads(data)["sequence_comparisons"][0]["matches"]
        status, _, data = router.dispatch("GET", f"/api/runs/{run_id}/export.html")
        assert status == 200 and b"sequence_comparisons" in data and b"Known" in data
    finally:
        router.close()


def test_scn1a_references_outside_mane_window_have_known_sites():
    snapshot = {
        "promoter": {
            "sequence": "A" * 2500,
            "assembly": "GRCh38",
            "contig": "2",
            "start": 166127520,
            "end": 166130020,
        },
        "candidates": [],
    }
    records = compare(snapshot, references("SCN1A"))
    assert all(r["status"] == "Outside current promoter window" for r in records)
    assert all(not r["matches"] for r in records)
    assert records[0]["known_location"]["end"] == 166149338
    assert records[1]["known_location"]["end"] == 166149958
    assert all(r["known_location"]["anchor_id"] == "scn1a-long-2008" for r in records)


def test_scn1a_reference_cores_match_long_promoter_context():
    from sarna.validation.dataset import load_chapter4_locations

    locations = load_chapter4_locations()
    for record in references("SCN1A"):
        context = locations[record["name"]].promoter()
        snapshot = {"promoter": {"sequence": context.sequence.bases}, "candidates": []}
        result = compare(snapshot, [record])[0]
        assert result["status"] == "matched"
        assert result["matches"][0]["offset"] == 250


def test_gene_lookup_exposes_known_sarnas_before_design(tmp_path: Path):
    from tests.unit.test_app import service_for

    gene = service_for(tmp_path).lookup("SCN1A")
    assert [r["name"] for r in gene["sarna_references"]] == ["SCN1A-10", "SCN1A-7"]
    assert all(
        r["known_location"]["anchor_id"] == "scn1a-long-2008" for r in gene["sarna_references"]
    )


def test_setup_window_coverage_and_strand_boundaries():
    import json
    import shutil
    import subprocess

    if not shutil.which("node"):
        pytest.skip("Node needed to verify browser coverage calculations")
    script = Path("src/sarna/app/static/app.js").read_text()
    function = script.split("function promoterSaRNACoverage(", 1)[1].split(
        "function renderKnownSaRNAPromoters()",
        1,
    )[0]
    function = "function promoterSaRNACoverage(" + function
    gene = {
        "assembly": "GRCh38",
        "contig": "2",
        "strand": "-",
        "transcripts": [{"id": "MANE", "tss": 166128019}],
        "anchors": [{"id": "scn1a-long-2008", "reference_point": 166149213}],
        "sarna_references": references("SCN1A"),
    }
    program = (
        function
        + "\nconst gene="
        + json.dumps(gene)
        + """;
const full=promoterSaRNACoverage(gene,{},2000,500);
if(full[0].contains.length!==0 || full[1].contains.length!==2) process.exit(1);
const small=promoterSaRNACoverage(gene,{},200,500);
if(small[1].contains.map(r=>r.name).join(',')!=='SCN1A-10') process.exit(2);
const test={assembly:'GRCh38',contig:'1',strand:'+',transcripts:[{id:'x',tss:100}],
  sarna_references:[{name:'edge',known_location:{assembly:'GRCh38',contig:'1',strand:'+',start:90,end:110}}]};
if(promoterSaRNACoverage(test,{},10,10)[0].contains.length!==1) process.exit(3);
if(promoterSaRNACoverage(test,{},9,10)[0].contains.length!==0) process.exit(4);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)


def test_known_sarna_setup_panel_exists_and_renders():
    import json
    import shutil
    import subprocess
    from html.parser import HTMLParser

    class Elements(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.ids: set[str] = set()

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            self.ids.update(value for key, value in attrs if key == "id" and value)

    html = Path("src/sarna/app/static/index.html").read_text()
    parser = Elements()
    parser.feed(html)
    assert {
        "known-sarna-promoters",
        "known-sarna-selection",
        "reference",
        "upstream",
        "downstream",
    } <= parser.ids
    if not shutil.which("node"):
        pytest.skip("Node needed for setup render smoke check")
    script = Path("src/sarna/app/static/app.js").read_text()
    functions = (
        "function promoterSaRNACoverage("
        + script.split(
            "function promoterSaRNACoverage(",
            1,
        )[1].split("function renderGeneOverview()", 1)[0]
    )
    gene = {
        "assembly": "GRCh38",
        "contig": "2",
        "strand": "-",
        "transcripts": [],
        "anchors": [{"id": "scn1a-long-2008", "reference_point": 166149213}],
        "sarna_references": references("SCN1A"),
    }
    program = """const nodes=Object.fromEntries(IDS.map(id=>[id,
  {innerHTML:'',textContent:'',value:'',options:[]} ]));
const $=id=>nodes[id] || null;
nodes.upstream.value='2000';nodes.downstream.value='500';
nodes.reference.value='anchor:scn1a-long-2008';
const escapeHTML=String, fmt=String, frozen=null;
const state={gene:GENE,customReferences:{}};
""".replace("IDS", json.dumps(sorted(parser.ids))).replace("GENE", json.dumps(gene))
    program += (
        functions
        + """
renderKnownSaRNAPromoters();
if(!nodes['known-sarna-promoters'].innerHTML.includes('SCN1A-7'))process.exit(1);
if(!nodes['known-sarna-selection'].textContent.includes('SCN1A-10'))process.exit(2);
state.gene=null;
renderKnownSaRNAPromoters();
if(nodes['known-sarna-promoters'].innerHTML!=='')process.exit(3);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)


@pytest.mark.parametrize(
    "symbol,name",
    [
        ("VEGFA", "VEGF"),
        ("VEGF", "VEGF"),
        ("TP53", "P53"),
        ("P53", "P53"),
        ("VHL", "VHL-821"),
        ("CEBPA", "CEBPA"),
        ("KLF4", "KLF4"),
    ],
)
def test_thesis_positive_controls_use_recorded_cores_and_locations(symbol: str, name: str):
    from sarna.validation.dataset import load_chapter4_dataset, load_chapter4_locations

    record = references(symbol)[0]
    original = load_chapter4_dataset().by_name(name)
    assert record["sequence"] == original.target.bases
    assert record["source"] == original.evidence
    assert original.is_activating
    location = load_chapter4_locations()[name]
    snapshot = {"promoter": {"sequence": location.context.bases}, "candidates": []}
    matched = compare(snapshot, [record])[0]
    assert matched["matches"][0]["offset"] == location.context_offset
    assert record["known_location"]["start"] == location.target.start


def test_finished_design_keeps_user_tab_until_explicit_review():
    import shutil
    import subprocess

    if not shutil.which("node"):
        pytest.skip("Node needed to verify navigation behavior")
    source = Path("src/sarna/app/static/workflow.js").read_text()
    run_function = (
        "function workflowRun("
        + source.split("function workflowRun(", 1)[1].split(
            "function workflowStatus()",
            1,
        )[0]
    )
    batch_function = (
        "function workflowBatch("
        + source.split("function workflowBatch(", 1)[1].split(
            "function renderExportSelection()",
            1,
        )[0]
    )
    program = (
        """
const state={workspaceTab:'design',selected:null};
const workspaces=new Map(), frozen=null, calls=[];
const nodes={'design-ready':{hidden:true},'design-ready-message':{textContent:''}};
const $=id=>nodes[id];
const document={body:{classList:{add(){}}}};
const renderShortlist=()=>{}, workflowStatus=()=>{};
const setWorkspaceTab=(tab,options)=>{calls.push([tab,options]);state.workspaceTab=tab;};
"""
        + run_function
        + batch_function
        + """
const result={app:{run_id:'one'},counts:{passing:5,generated:50},candidates:[]};
workflowRun(result);
if(state.workspaceTab!=='design'||nodes['design-ready'].hidden)process.exit(1);
if(!nodes['design-ready-message'].textContent.includes('5 of 50'))process.exit(2);
state.workspaceTab='screening';result.app.run_id='two';workflowRun(result);
if(state.workspaceTab!=='screening')process.exit(3);
if(calls.some(call=>call[1].scroll!==false))process.exit(4);
state.workspaceTab='design';state.result=null;workflowBatch({id:'group',status:'complete'});
if(state.workspaceTab!=='design')process.exit(5);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)
    html = Path("src/sarna/app/static/index.html").read_text()
    assert 'id="design-ready"' in html
    assert 'data-workspace-link="candidates">Review candidates when ready' in html


def test_reference_track_separates_crowded_labels_and_marks_downregulation():
    import shutil
    import subprocess

    if not shutil.which("node"):
        pytest.skip("Node needed to verify reference track layout")
    source = Path("src/sarna/app/static/app.js").read_text()
    functions = (
        "function sequenceAppearance("
        + source.split(
            "function sequenceAppearance(",
            1,
        )[1].split("function renderMap()", 1)[0]
    )
    program = (
        """
const escapeHTML=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;');
"""
        + functions
        + """
const match={offset:40,end:59};
const positive={name:'SCN1A-10',reference:true,source:'source',matches:[match]};
const down={name:'SCN1A-7',reference:true,source:'HeLa / AF22',matches:[match]};
if(sequenceAppearance(down).color===sequenceAppearance(positive).color)process.exit(1);
if(!sequenceEffectBadge(down).includes('Downregulating'))process.exit(2);
const single=existingSaRNATrack([positive],0,100,160,800,null);
const crowded=existingSaRNATrack([positive,down],0,100,160,800,null);
if(crowded.height<=single.height)process.exit(3);
if(!crowded.markup.includes('↓ Downregulating'))process.exit(4);
if(!crowded.markup.includes('data-sequence="1" data-match="0"'))process.exit(5);
if(!crowded.markup.includes('role="button" tabindex="0"'))process.exit(6);
const empty=existingSaRNATrack([positive],60,100,160,800,null);
if(empty.markup.includes('data-sequence='))process.exit(7);
const imported={...down,reference:false};
if(sequenceAppearance(imported).kind!=='imported')process.exit(8);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)


@pytest.mark.parametrize("strand,sense_pair,antisense_pair", [("+", "-", "+"), ("-", "+", "-")])
def test_strand_analysis_distinguishes_loading_from_genomic_pairing(
    strand: str,
    sense_pair: str,
    antisense_pair: str,
):
    from sarna.app.sequence_import import strand_analysis

    models = [
        {"predicted_guide": "antisense"},
        {"predicted_guide": "sense"},
        {"predicted_guide": "ambiguous"},
    ]
    analysis = strand_analysis("ATGC", strand, models)
    sense, antisense = analysis["strands"]
    assert sense["rna_5to3"] == "AUGC"
    assert antisense["rna_5to3"] == "GCAU"
    assert sense["counts"] == {"A": 1, "C": 1, "G": 1, "U": 1}
    assert sense["gc_percent"] == 50
    assert sense["complementary_dna_aligned_3to5"] == "TACG"
    assert sense["complementary_genomic_strand"] == sense_pair
    assert antisense["complementary_genomic_strand"] == antisense_pair
    assert [m["predicted_passenger"] for m in analysis["models"]] == ["sense", "antisense", None]


def test_outside_reference_has_analysis_without_claiming_a_loaded_match():
    snapshot = {"promoter": {"sequence": "A" * 40}, "candidates": []}
    record = compare(snapshot, references("SCN1A"))[0]
    assert not record["matches"]
    assert record["strand_analysis"]["strands"][1]["complementary_genomic_strand"] == "-"
    assert record["strand_analysis"]["models"]


def test_candidate_control_identity_uses_exact_core_not_overlap():
    import shutil
    import subprocess

    if not shutil.which("node"):
        pytest.skip("Node needed for candidate control identity")
    source = Path("src/sarna/app/static/app.js").read_text()
    functions = (
        "function candidateControlMatches("
        + source.split(
            "function candidateControlMatches(",
            1,
        )[1].split("function renderMap()", 1)[0]
    )
    program = (
        """const escapeHTML=String;
const state={result:{sequence_comparisons:[
{name:'SCN1A-7',reference:true,dna:'ATGC',matches:[]},
{name:'Imported',reference:false,dna:'ATGC',matches:[]}]}};
"""
        + functions
        + """
const exact={target_site_dna:'ATGC'};
if(candidateControlMatches(exact).length!==1)process.exit(1);
if(!candidateControlBadges(exact).includes('SCN1A-7'))process.exit(2);
if(!candidateControlBadges(exact).includes('Downregulating'))process.exit(3);
if(candidateControlMatches({target_site_dna:'TGCA'}).length)process.exit(4);
if(candidateControlMatches({target_site_dna:'GCAT'}).length!==1)process.exit(5);
if(candidateControlMatches({target_site_dna:'ATGCA'}).length)process.exit(6);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)


def test_existing_reference_displays_candidate_style_roles_without_expanding():
    import json
    import shutil
    import subprocess

    if not shutil.which("node"):
        pytest.skip("Node needed for automatic reference duplex rendering")
    source = Path("src/sarna/app/static/app.js").read_text()
    helpers = (
        "function candidateControlMatches("
        + source.split(
            "function candidateControlMatches(",
            1,
        )[1].split("function renderMap()", 1)[0]
    )
    ladder = (
        "const RNA_PAIR ="
        + source.split("const RNA_PAIR =", 1)[1].split(
            "function duplexPanel(",
            1,
        )[0]
    )
    renderer = (
        "function renderSequenceRecords("
        + source.split(
            "function renderSequenceRecords(",
            1,
        )[1].split("function renderSequenceDetail()", 1)[0]
    )
    records = compare({"promoter": {"sequence": "A" * 40}, "candidates": []}, references("SCN1A"))
    program = (
        """
const nodes={'sequence-import-controls':{},
 'sequence-records':{innerHTML:'',querySelectorAll:()=>[]},
 'sequence-analysis-notice':{},'sequence-recalculate':{}};
const $=id=>nodes[id], escapeHTML=String,fmt=String,frozen=null;
const state={result:{configuration:{rules:[]}},sequenceRecords:RECORDS};
const renderSequenceDetail=()=>{};
""".replace("RECORDS", json.dumps(records))
        + helpers
        + ladder
        + renderer
        + """
renderSequenceRecords();
const html=nodes['sequence-records'].innerHTML;
if(!html.includes('Predicted guide')||!html.includes('Predicted passenger'))process.exit(1);
if(!html.includes('ladder-row guide')||!html.includes('ladder-row passenger'))process.exit(2);
if(!html.includes('Models disagree'))process.exit(3);
if(!html.includes('SCN1A-7 · guide &amp; passenger prediction'))process.exit(4);
const main=html.indexOf('reference-duplex');
const more=html.indexOf('<details class="reference-strand-disclosure"');
if(main<0||more<main)process.exit(5);
if(html.includes('3 prime dTdT overhang'))process.exit(6);
"""
    )
    subprocess.run(["node", "-e", program], check=True, capture_output=True, text=True)
    html = Path("src/sarna/app/static/index.html").read_text()
    assert 'id="candidate-controls"' in html
    off_target = Path("src/sarna/app/static/offtarget.js").read_text()
    assert (
        "<strong>${escapeHTML(d.candidate.id)}</strong>${typeof candidateControlBadges"
        in off_target
    )
