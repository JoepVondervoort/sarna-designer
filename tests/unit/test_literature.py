"""Literature identity, observation provenance and conservative sequence placement."""

from copy import deepcopy
from pathlib import Path
from typing import Any

from sarna.app.literature import for_gene, locate
from sarna.app.models import DesignRequest
from sarna.app.service import DesignService
from sarna.app.web import AppRouter, frozen_page
from tests.unit.test_ensembl import MANE, client_for


def snapshot(sequence: str) -> dict[str, Any]:
    """Supply only fields needed by the independent literature placement operation."""
    return {
        "promoter": {"sequence": sequence},
        "candidates": [{"id": "candidate-1", "offset": 3, "target_site_dna": sequence[3:22]}],
    }


def test_release_links_observations_to_original_sources_and_preserves_unknowns():
    evidence = for_gene("ENSG00000124762")
    assert evidence["version"] == "1.0.0"
    assert evidence["record_count"] == 700
    assert len(evidence["records_sha256"]) == 64
    group = next(group for group in evidence["sequences"] if group["sequence_id"] == "SEQ000024")
    assert len(group["observations"]) > 1
    assert len({record["cell_type_or_model_reported"] for record in group["observations"]}) > 1
    for record in group["observations"]:
        assert record["publication_id"] == record["source"]["publication_id"]
        assert record["source"]["url"].startswith((
            "https://doi.org/10.", "https://pubmed.ncbi.nlm.nih.gov/",
            "https://patents.google.com/patent/",
        )) or not record["source"]["url"]
    assert any(record["effect_reported"] == "not_reported"
               for item in evidence["sequences"] for record in item["observations"])
    # Callers cannot change the cached source data or another run's evidence.
    group["observations"][0]["effect_reported"] = "changed"
    assert all(record["effect_reported"] != "changed"
               for item in for_gene("ENSG00000124762")["sequences"]
               for record in item["observations"])
    assert for_gene("ENSG00000144285")["sequences"] == []


def test_overhang_matching_requires_explicit_annotation_and_complementary_cores():
    evidence = for_gene("ENSG00000124762")
    group = next(group for group in evidence["sequences"] if group["sequence_id"] == "SEQ000024")
    evidence["sequences"] = [group]
    target = "CCAACTCATTCTCCAAGTA"
    run = snapshot("GGG" + target + "CCC")
    original = deepcopy(run)
    placed = locate(evidence, run)
    matches = placed["sequences"][0]["matches"]
    assert len(matches) == 1
    assert matches[0]["candidate_id"] == "candidate-1"
    assert matches[0]["offset"] == 3
    assert matches[0]["end"] == 22
    assert "dTdT" in matches[0]["basis"]
    records = {record["record_id"]: record for record in group["observations"]}
    assert all("dtdt overhang" in records[key]["modifications_reported"].lower()
               for key in matches[0]["record_ids"])
    exported = AppRouter._literature_sheets({"literature": placed})[0]
    assert all(bool(row[-1]) == (row[0] in matches[0]["record_ids"])
               for row in exported.rows)
    assert run == original
    assert "matches" not in evidence["sequences"][0]
    for record in group["observations"]:
        record["modifications_reported"] = "not_reported"
    assert locate(evidence, run)["sequences"][0]["matches"] == []
    group["observations"][0]["modifications_reported"] = "dTdT Overhang antisense strand"
    assert locate(evidence, run)["sequences"][0]["matches"] == []
    group["observations"][0]["modifications_reported"] = "dTdT Overhang"
    group["antisense"] = "A" * 19 + "TT"
    assert locate(evidence, run)["sequences"][0]["matches"] == []


def test_full_strand_matches_both_orientations_without_activity_inference():
    evidence = {"sequences": [{
        "sequence_id": "test", "sense": "CCAACTCATTCTCCAAGTA", "antisense": "no_sequence_reported",
        "observations": [{"record_id": "observation", "modifications_reported": "not_reported"}],
    }]}
    result = locate(evidence, snapshot("GGGTACTTGGAGAATGAGTTGGCCC"))
    match = result["sequences"][0]["matches"][0]
    assert match["offset"] == 3
    assert match["basis"] == "Full reported strand"
    assert "activity" not in match


def test_exports_preserve_observations_and_release_identity():
    evidence = locate(for_gene("ENSG00000124762"), snapshot("G" * 25))
    sheets = AppRouter._literature_sheets({"literature": evidence})
    assert len(sheets[0].rows) == 700
    assert "effect_reported" in sheets[0].headers
    assert ("version", "1.0.0") in sheets[1].rows
    assert AppRouter._literature_sheets({}) == []
    # Frozen reports carry their embedded observations and UI without a live database lookup.
    html = frozen_page({"literature": evidence}, None).decode()
    assert evidence["records_sha256"] in html
    assert 'id="literature-result"' in html
    assert 'renderLiteraturePanel' in html


def literature_service(path: Path) -> DesignService:
    """Use a synthetic annotation/window containing the published p21 core, offline."""
    gene = {
        "id": "ENSG00000124762", "display_name": "CDKN1A", "seq_region_name": "6",
        "start": 36678001, "end": 36680000, "strand": 1, "assembly_name": "GRCh38",
        "description": "Synthetic annotation fixture for literature integration",
    }
    transcript = {
        "id": "ENST_TEST", "biotype": "protein_coding", "seq_region_name": "6",
        "start": 36678426, "end": 36680000, "strand": 1, "is_canonical": 1, "MANE": MANE,
    }
    client = client_for({
        "/lookup/symbol/homo_sapiens/CDKN1A": gene,
        "/lookup/id/ENSG00000124762": {
            **gene, "canonical_transcript": "ENST_TEST", "Transcript": [transcript],
        },
        "/sequence/region/": {"seq": "GGGCCAACTCATTCTCCAAGTACCCCCCCC"},
    })
    return DesignService(path, lambda _: client)


def test_gene_design_freezes_literature_and_marks_excluded_candidates(tmp_path: Path):
    service = literature_service(tmp_path)
    overview = service.lookup("CDKN1A")
    assert overview["literature"]["record_count"] == 700
    request = DesignRequest(symbol="CDKN1A", transcript_id="ENST_TEST", upstream=25,
                            downstream=5, regulatory=False, motifs=False)
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    assert service.store.status(run_id)["status"] == "complete"
    saved = service.store.snapshot(run_id)
    assert saved["app"]["gene_overview"]["literature"] == overview["literature"]
    matches = [match for group in saved["literature"]["sequences"]
               for match in group["matches"] if match["candidate_id"]]
    assert matches
    candidate = next(candidate for candidate in saved["candidates"]
                     if candidate["id"] == matches[0]["candidate_id"])
    assert candidate["rank"] is None
    assert saved["configuration"]["preset"]["id"] == "wang-2015"
    reopened = DesignService(tmp_path)
    assert reopened.store.snapshot(run_id)["literature"] == saved["literature"]
