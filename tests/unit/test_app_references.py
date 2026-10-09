"""General gene/reference inputs, with independently computed coordinate expectations."""

from copy import deepcopy
from pathlib import Path
from typing import Any, Never

import pytest
from pydantic import ValidationError

from sarna.app.models import BatchRequest, CustomReference, DesignRequest
from sarna.app.references import tss_groups
from sarna.app.service import DesignService
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import AssemblyMismatchResponseError, EnsemblUnavailableError, Transcript
from tests.unit.test_ensembl import (
    AUTO_SEQUENCE,
    GENE_PAYLOAD,
    GENE_ROUTE,
    TRANSCRIPTS_PAYLOAD,
    TRANSCRIPTS_ROUTE,
    client_for,
    requests_made,
)

CUSTOM = {
    "id": "my-promoter", "name": "Lab-defined promoter", "tss_position": 200,
    "source": "Researcher notebook, experiment 12",
    "note": "Candidate alternative start; association with SCN1A needs verification.",
}


def recorded_service(tmp_path: Path, strand: int = -1) -> DesignService:
    """Reuse recorded transcript annotation, with a precisely sized genomic sequence."""
    gene = {**GENE_PAYLOAD, "strand": strand}
    transcripts = deepcopy(TRANSCRIPTS_PAYLOAD)
    for transcript in transcripts["Transcript"]:
        transcript["strand"] = strand
    client = client_for({
        GENE_ROUTE: gene, TRANSCRIPTS_ROUTE: transcripts,
        "/sequence/region/": {"seq": AUTO_SEQUENCE},
    })
    return DesignService(tmp_path, lambda _: client)


def complete(service: DesignService, request: DesignRequest) -> dict[str, Any]:
    """Exercise the complete scientific job and inspect its saved record."""
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    status = service.store.status(run_id)
    assert status["status"] == "complete", status["error"]
    return service.store.snapshot(run_id)


@pytest.mark.strand
@pytest.mark.parametrize("strand,start,end,first_start,first_end,sequence", [
    (1, 174, 204, 174, 193, "AAACAAACAAACAAACAAACAAACAAACAA"),
    (-1, 195, 225, 206, 225, "TTGTTTGTTTGTTTGTTTGTTTGTTTGTTT"),
])
def test_custom_reference_orients_each_strand_and_preserves_supplied_source(
    tmp_path: Path, strand: int, start: int, end: int, first_start: int, first_end: int,
    sequence: str,
):
    # Input base 200 is genomic coordinate 199. On minus strand, 25 upstream bases
    # occupy [200,225) and five downstream bases [195,200). Plus mirrors that.
    service = recorded_service(tmp_path, strand)
    snapshot = complete(service, DesignRequest(
        symbol="SCN1A", custom_reference=CustomReference.model_validate(CUSTOM),
        upstream=25, downstream=5, regulatory=False, motifs=False,
    ))
    promoter = snapshot["promoter"]
    assert (promoter["start"], promoter["end"]) == (start, end)
    assert promoter["sequence"] == sequence
    first = next(candidate for candidate in snapshot["candidates"] if candidate["offset"] == 0)
    assert first["tss_span"] == [-25, -7]
    assert (first["genomic"]["start"], first["genomic"]["end"]) == (first_start, first_end)
    anchor = snapshot["app"]["anchor"]
    assert anchor["reference_point"] == 199
    assert anchor["tss_position"] == 200
    assert anchor["reference_kind"] == "inferred"
    assert anchor["researcher_supplied"] is True
    assert anchor["source"] == CUSTOM["source"]
    assert anchor["note"] == CUSTOM["note"]
    assert anchor["name"] == CUSTOM["name"]
    assert "not been independently verified" in anchor["caveat"]
    # It is deliberately outside the annotated gene span; the researcher links it.
    assert anchor["end"] < snapshot["app"]["gene_overview"]["start"]
    assert snapshot["transcript"] is None


def test_custom_reference_outside_chromosome_fails_without_sequence_retrieval(tmp_path: Path):
    service = recorded_service(tmp_path)
    client = service.client_factory(False)
    request = DesignRequest(
        symbol="SCN1A", custom_reference=CustomReference.model_validate({
            **CUSTOM, "tss_position": 242_193_530,
        }), regulatory=False, motifs=False,
    )
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    assert "beyond the end of contig" in service.store.status(run_id)["error"]
    assert not any("/sequence/" in url for url in requests_made(client))


def test_custom_reference_window_cannot_extend_past_chromosome_start(tmp_path: Path):
    service = recorded_service(tmp_path, strand=1)
    request = DesignRequest(
        symbol="SCN1A", custom_reference=CustomReference.model_validate({
            **CUSTOM, "tss_position": 1,
        }), upstream=25, downstream=5, regulatory=False, motifs=False,
    )
    run_id = service.store.create(request.model_dump(mode="json"))
    service.run(run_id, request)
    assert "before the beginning of contig" in service.store.status(run_id)["error"]


@pytest.mark.parametrize("change", [
    {"id": ""}, {"id": "../unsafe"}, {"name": "  "}, {"source": " "}, {"note": ""},
    {"tss_position": 0}, {"tss_position": -1}, {"tss_position": 1.5},
    {"tss_position": True}, {"tss_position": "200"},
])
def test_invalid_custom_reference_is_rejected(change: dict[str, Any]):
    with pytest.raises(ValidationError):
        CustomReference.model_validate({**CUSTOM, **change})


@pytest.mark.parametrize("values", [
    {"symbol": "SCN1A", "transcript_id": "ENST1", "anchor_id": "long"},
    {"symbol": "SCN1A", "transcript_id": "ENST1", "custom_reference": CUSTOM},
    {"symbol": "SCN1A", "anchor_id": "long", "custom_reference": CUSTOM},
    {"source": "sequence", "sequence": "ACGT" * 10, "custom_reference": CUSTOM},
    {"symbol": "SCN1A", "transcript_id": ""},
    {"symbol": "SCN1A", "anchor_id": " "},
])
def test_reference_choices_are_exclusive_and_custom_requires_gene(values: dict[str, Any]):
    with pytest.raises(ValidationError):
        DesignRequest.model_validate(values)


def test_lookup_groups_identical_tss_but_keeps_nearby_sites_separate(tmp_path: Path):
    transcripts = deepcopy(TRANSCRIPTS_PAYLOAD)
    shared = {**transcripts["Transcript"][0], "id": "ENST_SHARED"}
    shared.pop("MANE")
    shared["start"] += 100
    nearby = {**shared, "id": "ENST_NEARBY", "end": shared["end"] + 1}
    transcripts["Transcript"].extend([shared, nearby])
    client = client_for({GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: transcripts})
    overview = DesignService(tmp_path, lambda _: client).lookup("scn1a")
    assert len(overview["transcripts"]) == 5
    assert len(overview["tss_groups"]) == 4
    shared_group = next(group for group in overview["tss_groups"]
                        if group["id"] == "tss:2:-:166128019")
    assert shared_group["transcript_ids"] == ["ENST00000674923", "ENST_SHARED"]
    assert shared_group["mane"] is True
    nearby_group = next(group for group in overview["tss_groups"]
                        if group["tss"] == 166128020)
    assert nearby_group["transcript_ids"] == ["ENST_NEARBY"]
    assert nearby_group["mane"] is False
    assert overview["contig"] == "2"
    assert (overview["start"], overview["end"]) == (165984640, 166182806)


def test_tss_grouping_includes_contig_and_strand_not_only_coordinate():
    transcripts = tuple(Transcript(
        identifier, "protein_coding", GenomicInterval(Assembly.GRCH38, contig, start, end,
                                                       strand), False, False,
    ) for identifier, contig, start, end, strand in (
        ("same-plus", "1", 100, 150, Strand.PLUS),
        ("same-minus", "1", 50, 101, Strand.MINUS),
        ("same-other-contig", "2", 100, 150, Strand.PLUS),
    ))
    groups = tss_groups(transcripts)
    assert [group["id"] for group in groups] == ["tss:1:+:100", "tss:1:-:100", "tss:2:+:100"]


def test_saved_gene_run_freezes_all_reference_choices_without_relookup(tmp_path: Path):
    service = recorded_service(tmp_path)
    lookup = service.lookup("SCN1A")
    snapshot = complete(service, DesignRequest(
        symbol="SCN1A", transcript_id="ENST00000409050", upstream=25, downstream=5,
        regulatory=False, motifs=False,
    ))
    assert snapshot["app"]["gene_overview"] == lookup
    # Opening the existing SQLite snapshot requires no new client/source calls.
    def unavailable(_: bool) -> Never:
        raise AssertionError("saved reference choices must not need a network request")
    reopened = DesignService(tmp_path, unavailable)
    assert reopened.store.snapshot(snapshot["app"]["run_id"])["app"]["gene_overview"] == lookup


@pytest.mark.parametrize("identifier", ["ENSG00000144285", "ENSG00000144285.7"])
def test_human_gene_id_lookup_uses_stable_id_endpoint_and_current_assembly(
    tmp_path: Path, identifier: str,
):
    gene = {**GENE_PAYLOAD, "species": "homo_sapiens", "object_type": "Gene"}
    client = client_for({
        "/lookup/id/ENSG00000144285?expand=0": gene,
        TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
        "/sequence/region/": {"seq": AUTO_SEQUENCE},
    }, cache_dir=tmp_path / "cache")
    service = DesignService(tmp_path / "runs", lambda _: client)
    lookup = service.lookup(identifier.lower())
    assert lookup["symbol"] == "SCN1A"
    assert lookup["id"] == "ENSG00000144285"
    assert requests_made(client)[0].endswith("/lookup/id/ENSG00000144285?expand=0")
    snapshot = complete(service, DesignRequest(
        symbol=identifier.lower(), transcript_id="ENST00000409050", upstream=25,
        downstream=5, regulatory=False, motifs=False,
    ))
    assert snapshot["app"]["gene_overview"] == lookup
    assert snapshot["app"]["request"]["symbol"] == identifier
    cached = client_for({}, cache_dir=tmp_path / "cache", offline=True)
    assert cached.gene(identifier).symbol == "SCN1A"
    assert requests_made(cached) == []


@pytest.mark.parametrize("change,error", [
    ({"species": "mus_musculus"}, ValueError),
    ({"object_type": "Transcript"}, ValueError),
    ({"species": None}, ValueError),
    ({"id": "ENSG00000000001"}, EnsemblUnavailableError),
    ({"assembly_name": None}, EnsemblUnavailableError),
    ({"assembly_name": "GRCh37"}, AssemblyMismatchResponseError),
])
def test_gene_id_lookup_refuses_wrong_species_type_identity_or_assembly(
    change: dict[str, Any], error: type[Exception],
):
    payload = {**GENE_PAYLOAD, "species": "homo_sapiens", "object_type": "Gene", **change}
    with pytest.raises(error):
        client_for({TRANSCRIPTS_ROUTE: payload}).gene("ENSG00000144285")


@pytest.mark.parametrize("identifier", ["ENST00000674923", "ENSMUSG00000000001", "../bad"])
def test_gene_id_method_rejects_nonhuman_gene_ids_without_fetching(identifier: str):
    client = client_for({})
    with pytest.raises(ValueError):
        client.gene_by_id(identifier)
    assert requests_made(client) == []


def test_batch_requires_explicit_references_for_one_normalised_gene():
    batch = BatchRequest.model_validate({"requests": [
        {"symbol": "scn1a", "transcript_id": "ENST00000674923"},
        {"symbol": " SCN1A ", "anchor_id": "scn1a-long-2008"},
        {"symbol": "SCN1A", "custom_reference": CUSTOM},
    ]})
    assert batch.label == ""
    assert len(batch.requests) == 3
    assert {request.symbol for request in batch.requests} == {"SCN1A"}
    assert all(request.resolve().preset_id == "wang-2015" for request in batch.requests)


@pytest.mark.parametrize("requests", [
    [], [{"symbol": "SCN1A"}],
    [{"source": "sequence", "sequence": "ACGT" * 10}],
    [{"symbol": "SCN1A", "transcript_id": "ENST1"},
     {"symbol": "CHD2", "transcript_id": "ENST2"}],
    [{"symbol": "SCN1A", "transcript_id": "ENST1"}] * 33,
])
def test_invalid_batches_are_rejected_before_lookup(requests: list[dict[str, Any]]):
    with pytest.raises(ValidationError):
        BatchRequest.model_validate({"requests": requests})


def test_the_default_reference_is_offered_first_whatever_order_ensembl_returns(tmp_path: Path):
    """The reference list leads with MANE Select, then canonical, then coding, then the rest.

    Built from a payload whose recorded order is the reverse of the wanted one, so the
    assertion is about the sort and not about the order the fixture happens to carry.
    """
    transcripts = deepcopy(TRANSCRIPTS_PAYLOAD)
    recorded = transcripts["Transcript"]
    mane = recorded[0]
    # A non-coding transcript and a second coding one, each on their own start coordinate, so
    # every tier of the ordering is represented and the groups cannot collapse into one.
    extra = [
        {**mane, "id": "ENST_LNC", "biotype": "lncRNA", "end": mane["end"] + 300,
         "is_canonical": 0, "MANE": None},
        {**mane, "id": "ENST_CODING", "end": mane["end"] + 600, "is_canonical": 0,
         "MANE": None},
    ]
    for entry in extra:
        entry.pop("MANE")
    transcripts["Transcript"] = [*extra, *recorded[1:], mane]
    client = client_for({GENE_ROUTE: GENE_PAYLOAD, TRANSCRIPTS_ROUTE: transcripts})
    overview = DesignService(tmp_path, lambda _: client).lookup("SCN1A")

    assert overview["transcripts"][0]["id"] == "ENST00000674923"
    assert overview["transcripts"][0]["mane"] is True
    assert overview["default_transcript"] == "ENST00000674923"
    # The MANE group leads the grouped list the picker renders, not only the flat list.
    assert overview["tss_groups"][0]["mane"] is True
    assert "ENST00000674923" in overview["tss_groups"][0]["transcript_ids"]
    # Non-coding transcripts sort behind every coding one, and ties break on the identifier
    # so two lookups of the same gene offer the same order.
    identifiers = [record["id"] for record in overview["transcripts"]]
    biotypes = {record["id"]: record["biotype"] for record in overview["transcripts"]}
    coding = [i for i, name in enumerate(identifiers) if biotypes[name] == "protein_coding"]
    non_coding = [i for i, name in enumerate(identifiers) if biotypes[name] != "protein_coding"]
    assert max(coding) < min(non_coding)
    # Within the plain-coding tier nothing distinguishes the members but their identifier, so
    # the order is alphabetical and therefore the same on every lookup.
    plain = [name for name in identifiers
             if biotypes[name] == "protein_coding" and name != "ENST00000674923"]
    assert plain == sorted(plain)
