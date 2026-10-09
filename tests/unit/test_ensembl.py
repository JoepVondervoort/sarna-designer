"""Tests for the Ensembl client, run entirely offline against recorded payloads.

The client is the boundary where three silent assumptions become explicit decisions (which
assembly, which transcript, which coordinate convention) so the tests concentrate there rather
than on HTTP mechanics. The payloads are trimmed copies of real Ensembl responses for SCN1A, a
**minus-strand** gene with five distinct protein-coding TSSs, which is the hardest realistic case
and the one the project's primary target happens to be.

A live-network counterpart lives in `tests/network/test_ensembl_live.py`, excluded by default.
"""

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar, cast

import httpx
import pytest

from sarna import settings
from sarna.genome.anchors import PromoterAnchor, ReferenceKind
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import (
    AssemblyMismatchResponseError,
    EnsemblCacheMissError,
    EnsemblClient,
    EnsemblUnavailableError,
    GeneNotFoundError,
    NoManeSelectTranscriptError,
    Transcript,
    promoter_for_symbol,
)
from sarna.rules.provenance import Provenance

SCN1A_TSS = 166128019
"""0-based TSS of SCN1A's MANE Select transcript, the highest coordinate of its span."""

GENE_PAYLOAD: dict[str, Any] = {
    "id": "ENSG00000144285",
    "display_name": "SCN1A",
    "seq_region_name": "2",
    "start": 165984641,
    "end": 166182806,
    "strand": -1,
    "assembly_name": "GRCh38",
    "description": "sodium voltage-gated channel alpha subunit 1",
}

MANE = [{"type": "MANE_Select", "refseq_match": "NM_001165963.4"}]

TRANSCRIPTS_PAYLOAD: dict[str, Any] = {
    "id": "ENSG00000144285",
    "assembly_name": "GRCh38",
    "canonical_transcript": "ENST00000674923.1",
    "Transcript": [
        {
            "id": "ENST00000674923",
            "biotype": "protein_coding",
            "seq_region_name": "2",
            "start": 165985813,
            "end": 166128020,
            "strand": -1,
            "is_canonical": 1,
            "MANE": MANE,
        },
        {
            "id": "ENST00000409050",
            "biotype": "protein_coding",
            "seq_region_name": "2",
            "start": 165985813,
            "end": 166149214,
            "strand": -1,
        },
        {
            "id": "ENST00000650461",
            "biotype": "retained_intron",
            "seq_region_name": "2",
            "start": 166050000,
            "end": 166060000,
            "strand": -1,
        },
    ],
}


AUTO_SEQUENCE = "<auto>"
"""Sentinel for a sequence route that should size its reply to the region requested.

A fixed-length stub reply makes ``Promoter.from_reference_sequence`` reject any window of another
length, correctly, since an interval whose length disagrees with its sequence would shift every
candidate's coordinates. Sizing the reply keeps the tests about the thing under test.
"""


def _sized_sequence(url: str) -> str:
    """Return a sequence of exactly the length the region in ``url`` asks for."""
    region = url.split("/sequence/region/homo_sapiens/")[1].split("?")[0]
    _, span = region.split(":")
    start, end = span.split("..")
    length = int(end) - int(start) + 1
    return ("AAAC" * (length // 4 + 1))[:length]


class StubTransport:
    """Serves recorded payloads by URL fragment, and records what was asked for.

    An ``int`` value in ``routes`` stands for a status code with an HTML body, which is how
    Ensembl's error pages arrive.
    """

    def __init__(self, routes: Mapping[str, object]) -> None:
        """Record the routes to serve."""
        self.routes = dict(routes)
        self.requested: list[str] = []

    def get(self, url: str, headers: dict[str, str], timeout: float) -> httpx.Response:
        """Serve the first route whose fragment appears in ``url``."""
        self.requested.append(url)
        for fragment, payload in self.routes.items():
            if fragment in url:
                if isinstance(payload, int):
                    return httpx.Response(payload, text="<html>error</html>")
                if isinstance(payload, dict) and payload.get("seq") == AUTO_SEQUENCE:
                    return httpx.Response(200, json={"seq": _sized_sequence(url)})
                return httpx.Response(200, json=payload)
        return httpx.Response(404, json={"error": "not found"})


class FixedTransport:
    """Returns one fixed response, or raises one fixed transport error, every time."""

    def __init__(
        self, status: int = 200, body: str = "", error: Exception | None = None
    ) -> None:
        """Record what to return or raise."""
        self.status = status
        self.body = body
        self.error = error
        self.requested: list[str] = []

    def get(self, url: str, headers: dict[str, str], timeout: float) -> httpx.Response:
        """Return the fixed response, or raise the fixed error."""
        self.requested.append(url)
        if self.error is not None:
            raise self.error
        return httpx.Response(self.status, text=self.body)


def client_for(
    routes: Mapping[str, object],
    cache_dir: Path | None = None,
    offline: bool | None = None,
) -> EnsemblClient:
    """Build a client backed by a stub, with retries made instant."""
    return EnsemblClient(
        client=StubTransport(routes),  # type: ignore[arg-type]
        backoff_seconds=0.0,
        cache_dir=cache_dir,
        offline=offline,
    )


def client_returning(
    status: int = 200, body: str = "", error: Exception | None = None
) -> EnsemblClient:
    """Build a client whose transport always answers the same way."""
    return EnsemblClient(
        client=FixedTransport(status, body, error),  # type: ignore[arg-type]
        backoff_seconds=0.0,
    )


def requests_made(client: EnsemblClient) -> list[str]:
    """Return the URLs a stub-backed client asked for."""
    transport = cast("StubTransport | FixedTransport", client._client)
    return transport.requested


GENE_ROUTE = "/lookup/symbol/homo_sapiens/SCN1A"
TRANSCRIPTS_ROUTE = "/lookup/id/ENSG00000144285"

DEFAULT_ROUTES: dict[str, Any] = {
    GENE_ROUTE: GENE_PAYLOAD,
    TRANSCRIPTS_ROUTE: TRANSCRIPTS_PAYLOAD,
}


class TestAssemblyIsNamedNotInherited:
    def test_each_assembly_has_its_own_host(self) -> None:
        assert EnsemblClient(Assembly.GRCH38).base_url == "https://rest.ensembl.org"
        assert EnsemblClient(Assembly.GRCH37).base_url == "https://grch37.rest.ensembl.org"

    def test_the_default_is_grch38(self) -> None:
        assert EnsemblClient().assembly is Assembly.GRCH38

    def test_a_response_from_the_wrong_assembly_is_refused(self) -> None:
        # Defect B15 was mixing assemblies without noticing. A server whose default changed
        # must fail here rather than shift every coordinate silently.
        wrong = {**GENE_PAYLOAD, "assembly_name": "GRCh37"}
        client = client_for({"/lookup/symbol": wrong})
        with pytest.raises(AssemblyMismatchResponseError, match="B15"):
            client.gene("SCN1A")

    def test_a_response_without_an_assembly_field_is_accepted(self) -> None:
        # Some endpoints omit it. Absence is not a mismatch.
        payload = {k: v for k, v in GENE_PAYLOAD.items() if k != "assembly_name"}
        assert client_for({"/lookup/symbol": payload}).gene("SCN1A").symbol == "SCN1A"

    def test_a_sequence_request_for_another_assembly_is_refused(self) -> None:
        client = client_for({})
        interval = GenomicInterval(Assembly.GRCH37, "2", 1000, 2000)
        with pytest.raises(ValueError, match="this client serves GRCh38"):
            client.sequence(interval)


class TestGeneLookup:
    def test_a_gene_is_returned_with_its_span(self) -> None:
        gene = client_for(DEFAULT_ROUTES).gene("SCN1A")
        assert gene.ensembl_id == "ENSG00000144285"
        assert gene.contig == "2"
        assert gene.strand is Strand.MINUS

    def test_coordinates_convert_from_1_based_inclusive_to_0_based_half_open(self) -> None:
        # Ensembl says start=165984641; the model holds 165984640. The conversion happens once,
        # here, rather than two conventions coexisting downstream.
        gene = client_for(DEFAULT_ROUTES).gene("SCN1A")
        assert gene.interval.start == 165984640
        assert gene.interval.end == 166182806
        assert gene.interval.length == 166182806 - 165984640

    def test_an_unknown_symbol_raises_and_mentions_case(self) -> None:
        with pytest.raises(GeneNotFoundError, match="case"):
            client_for({}).gene("NOTAGENE")


class TestTranscriptSelection:
    def test_the_mane_select_transcript_is_chosen(self) -> None:
        client = client_for(DEFAULT_ROUTES)
        transcript = client.select_transcript(client.gene("SCN1A"))
        assert transcript.ensembl_id == "ENST00000674923"
        assert transcript.is_mane_select

    def test_the_refseq_accession_behind_mane_select_is_recorded(self) -> None:
        # MANE Select's value is that two annotation projects agreed; the accession is what
        # lets a reader check the agreement.
        client = client_for(DEFAULT_ROUTES)
        assert client.select_transcript(client.gene("SCN1A")).refseq_match == "NM_001165963.4"

    def test_the_mane_request_asks_for_mane_annotations(self) -> None:
        # Without mane=1 Ensembl omits them entirely and every gene looks as though it has no
        # MANE Select transcript. This was a real bug, found against the live API.
        client = client_for(DEFAULT_ROUTES)
        client.transcripts(client.gene("SCN1A"))
        assert any("mane=1" in url for url in requests_made(client))

    def test_the_canonical_transcript_is_identified_despite_its_version_suffix(self) -> None:
        client = client_for(DEFAULT_ROUTES)
        transcripts = client.transcripts(client.gene("SCN1A"))
        canonical = [t for t in transcripts if t.is_canonical]
        assert [t.ensembl_id for t in canonical] == ["ENST00000674923"]

    def test_with_no_mane_select_the_choice_is_refused_and_alternatives_listed(self) -> None:
        # A promoter is a property of a transcript. Picking the first, longest or
        # lowest-numbered is how a tool designs against a promoter nobody chose.
        no_mane = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {k: v for k, v in t.items() if k != "MANE"}
                for t in TRANSCRIPTS_PAYLOAD["Transcript"]
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: no_mane})
        with pytest.raises(NoManeSelectTranscriptError) as caught:
            client.select_transcript(client.gene("SCN1A"))
        assert "will not make for you" in str(caught.value)
        assert len(caught.value.alternatives) == 2  # the protein-coding ones only
        assert "ENST00000674923" in str(caught.value)

    def test_the_refusal_prefers_protein_coding_alternatives(self) -> None:
        no_mane = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {k: v for k, v in t.items() if k != "MANE"}
                for t in TRANSCRIPTS_PAYLOAD["Transcript"]
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: no_mane})
        with pytest.raises(NoManeSelectTranscriptError) as caught:
            client.select_transcript(client.gene("SCN1A"))
        assert all(t.is_protein_coding for t in caught.value.alternatives)

    def test_more_than_one_mane_select_is_also_refused(self) -> None:
        # Should not happen, and must not be resolved by taking the first.
        two = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {**TRANSCRIPTS_PAYLOAD["Transcript"][0]},
                {**TRANSCRIPTS_PAYLOAD["Transcript"][1], "MANE": MANE},
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: two})
        with pytest.raises(NoManeSelectTranscriptError):
            client.select_transcript(client.gene("SCN1A"))

    def test_a_gene_with_only_noncoding_transcripts_lists_those(self) -> None:
        only_noncoding = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [TRANSCRIPTS_PAYLOAD["Transcript"][2]],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: only_noncoding})
        with pytest.raises(NoManeSelectTranscriptError) as caught:
            client.select_transcript(client.gene("SCN1A"))
        assert len(caught.value.alternatives) == 1


@pytest.mark.strand
class TestTheTssIsOrientedNotAssumed:
    def test_a_minus_strand_tss_is_the_highest_coordinate(self) -> None:
        # Defect B1. Using `start` on both strands puts a minus-strand gene's promoter at the
        # far end of the transcript, tens of kilobases away.
        client = client_for(DEFAULT_ROUTES)
        transcript = client.select_transcript(client.gene("SCN1A"))
        assert transcript.strand is Strand.MINUS
        assert transcript.tss_coordinate == SCN1A_TSS
        assert transcript.tss_coordinate == transcript.interval.end - 1

    def test_a_plus_strand_tss_is_the_lowest_coordinate(self) -> None:
        plus = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {**TRANSCRIPTS_PAYLOAD["Transcript"][0], "strand": 1},
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: plus})
        transcript = client.select_transcript(client.gene("SCN1A"))
        assert transcript.tss_coordinate == transcript.interval.start


@pytest.mark.strand
class TestPromoterWindows:
    SEQUENCE: ClassVar[dict[str, str]] = {"seq": "ACGT" * 625}  # 2500 bp, the default window

    def client(self) -> EnsemblClient:
        return client_for({**DEFAULT_ROUTES, "/sequence/region/": self.SEQUENCE})

    def transcript(self) -> Transcript:
        client = self.client()
        return client.select_transcript(client.gene("SCN1A"))

    def test_a_minus_strand_promoter_extends_to_higher_coordinates(self) -> None:
        # Upstream is toward higher coordinates on the minus strand, so the window's `end`
        # carries the upstream bases and its `start` the downstream ones.
        promoter = self.client().promoter_for(self.transcript())
        assert promoter.interval is not None
        assert promoter.interval.end == SCN1A_TSS + 2000 + 1
        assert promoter.interval.start == SCN1A_TSS - 500 + 1
        assert promoter.length == 2500

    def test_the_promoter_reproduces_the_transcripts_tss_independently(self) -> None:
        # Cross-check: the promoter derives the TSS from its own interval and tss_start, by
        # separate arithmetic from the transcript's five_prime. They must agree.
        transcript = self.transcript()
        promoter = self.client().promoter_for(transcript)
        assert promoter.tss_coordinate == transcript.tss_coordinate

    def test_the_first_base_is_at_minus_the_upstream_width(self) -> None:
        promoter = self.client().promoter_for(self.transcript(), upstream_bp=2000)
        assert promoter.tss_start is not None
        assert promoter.tss_start.value == -2000

    def test_a_minus_strand_promoter_sequence_is_reverse_complemented(self) -> None:
        # Ensembl serves the plus strand; a minus-strand promoter's own 5'->3' sequence is its
        # reverse complement. ACGT repeated reverse-complements to ACGT repeated, so use a
        # sequence that distinguishes them.
        routes = {**DEFAULT_ROUTES, "/sequence/region/": {"seq": "AAAC" * 625}}
        client = client_for(routes)
        promoter = client.promoter_for(client.select_transcript(client.gene("SCN1A")))
        assert promoter.sequence.bases.startswith("GTTT")

    def test_a_window_of_zero_length_is_refused(self) -> None:
        with pytest.raises(ValueError, match="zero length"):
            self.client().promoter_for(self.transcript(), 0, 0)

    def test_negative_window_sizes_are_refused(self) -> None:
        # They are magnitudes, not offsets; a negative would silently invert the window.
        with pytest.raises(ValueError, match="must not be negative"):
            self.client().promoter_for(self.transcript(), -100, 500)

    def test_a_window_running_off_the_contig_is_refused(self) -> None:
        near_start = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {
                    **TRANSCRIPTS_PAYLOAD["Transcript"][0],
                    "start": 10,
                    "end": 50,
                    "strand": 1,
                }
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: near_start})
        transcript = client.select_transcript(client.gene("SCN1A"))
        with pytest.raises(ValueError, match="before the beginning of contig"):
            client.promoter_for(transcript, upstream_bp=2000)

    def test_an_oversized_region_is_refused_before_the_request(self) -> None:
        client = self.client()
        interval = GenomicInterval(Assembly.GRCH38, "2", 0, 20_000_000)
        with pytest.raises(ValueError, match="exceeds Ensembl"):
            client.sequence(interval)


class TestRetriesAndFailures:
    def test_an_html_500_is_retried_and_then_reported(self) -> None:
        # Ensembl's public endpoint serves HTML error pages with 500 statuses under load.
        client = client_for({"/lookup/symbol": 500})
        with pytest.raises(EnsemblUnavailableError, match="after 6 attempts"):
            client.gene("SCN1A")

    def test_a_404_is_not_retried(self) -> None:
        client = client_for({})
        with pytest.raises(GeneNotFoundError):
            client.gene("NOTAGENE")
        assert len(requests_made(client)) == 1

    def test_a_4xx_other_than_404_stops_without_exhausting_retries(self) -> None:
        client = client_for({"/lookup/symbol": 400})
        with pytest.raises(EnsemblUnavailableError, match="HTTP 400"):
            client.gene("SCN1A")
        assert len(requests_made(client)) == 1

    def test_a_transport_error_is_reported_as_unavailable(self) -> None:
        client = client_returning(error=httpx.ConnectError("no route to host"))
        with pytest.raises(EnsemblUnavailableError, match="works offline"):
            client.gene("SCN1A")

    def test_a_200_carrying_html_is_treated_as_a_failure(self) -> None:
        client = client_returning(200, "<html>not json</html>")
        with pytest.raises(EnsemblUnavailableError, match="not JSON"):
            client.gene("SCN1A")


class TestCaching:
    def test_a_response_is_cached_and_reused(self, tmp_path: Path) -> None:
        # What lets a promoter fetched once be reproducible with no network.
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        client.gene("SCN1A")
        first = len(requests_made(client))
        client.gene("SCN1A")
        assert len(requests_made(client)) == first

    def test_a_truncated_cache_file_is_refetched_rather_than_poisoning_the_gene(
        self, tmp_path: Path
    ) -> None:
        """A half-written cache file must not make a gene permanently unusable.

        The cache was written with a plain ``write_text`` and read with a bare ``json.loads``.
        Interrupt the app mid-write -- Ctrl-C, a closed laptop, a full disk -- and the file is
        left truncated. Every later lookup of that gene then raised ``JSONDecodeError`` out of
        the client, on every run, until someone found and deleted the file by hand. There is no
        message anywhere telling a non-computational colleague to do that.
        """
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        client.gene("SCN1A")
        cached = next(f for f in tmp_path.iterdir() if f.suffix == ".json")
        whole = cached.read_text(encoding="utf-8")
        cached.write_text(whole[: len(whole) // 2], encoding="utf-8")

        recovered = client_for(DEFAULT_ROUTES, cache_dir=tmp_path).gene("SCN1A")
        assert recovered.symbol == "SCN1A"
        # And the damaged file has been replaced, so the next call is a cache hit again.
        assert json.loads(cached.read_text(encoding="utf-8"))

    def test_an_interrupted_write_never_leaves_a_partial_cache_file(
        self, tmp_path: Path
    ) -> None:
        # Written to a temporary name and renamed, so a reader sees the old file or the new
        # one and never half of either. This is also what makes a cache directory shared
        # between colleagues safe to read while somebody else is writing it.
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        client.gene("SCN1A")
        for path in tmp_path.iterdir():
            assert path.suffix == ".json", f"{path.name} was left behind by the cache write"
            assert json.loads(path.read_text(encoding="utf-8"))

    def test_the_cache_file_names_the_assembly(self, tmp_path: Path) -> None:
        # So a GRCh37 and a GRCh38 response for the same query cannot collide.
        client_for(DEFAULT_ROUTES, cache_dir=tmp_path).gene("SCN1A")
        assert any(f.name.startswith("GRCh38") for f in tmp_path.iterdir())

    def test_an_unusable_cache_directory_degrades_to_no_cache(self, tmp_path: Path) -> None:
        """A shared cache folder that is not mounted must not stop the tool working.

        INSTALL.md offers ``SARNA_CACHE_DIR`` pointed at a network share so a group fetches each
        gene once. A laptop off the VPN then has an unreachable path, and the constructor's
        ``mkdir`` raised ``FileNotFoundError`` before a single lookup was attempted -- turning a
        missing convenience into a tool that would not start.
        """
        unusable = tmp_path / "not-a-directory" / "cache"
        unusable.parent.write_text("this is a file, not a folder", encoding="utf-8")
        client = client_for(DEFAULT_ROUTES, cache_dir=unusable)
        assert client.gene("SCN1A").symbol == "SCN1A"
        # Caching is simply off; nothing was written and nothing failed.
        assert not unusable.exists()

    def test_the_cache_directory_is_created(self, tmp_path: Path) -> None:
        target = tmp_path / "nested" / "cache"
        client_for(DEFAULT_ROUTES, cache_dir=target).gene("SCN1A")
        assert target.is_dir()

    def test_no_cache_directory_means_no_caching(self) -> None:
        client = client_for(DEFAULT_ROUTES)
        client.gene("SCN1A")
        client.gene("SCN1A")
        assert len(requests_made(client)) == 2


class TestOfflineMode:
    """``SARNA_OFFLINE`` makes a recorded run reproducible by refusing to reach the network.

    The point is not convenience. A design run that silently re-fetches can silently get a
    different answer when Ensembl's annotation changes, and nothing in the output would say so.
    Offline turns that into a loud failure naming what was missing.
    """

    def test_a_cached_response_is_served_without_a_request(self, tmp_path: Path) -> None:
        warm = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        warm.gene("SCN1A")
        offline = client_for(DEFAULT_ROUTES, cache_dir=tmp_path, offline=True)
        assert offline.gene("SCN1A").symbol == "SCN1A"
        assert requests_made(offline) == []

    def test_a_cache_miss_refuses_rather_than_fetching(self, tmp_path: Path) -> None:
        offline = client_for(DEFAULT_ROUTES, cache_dir=tmp_path, offline=True)
        with pytest.raises(EnsemblCacheMissError, match="Cache-only mode") as failure:
            offline.gene("SCN1A")
        assert "SARNA_OFFLINE" not in str(failure.value)
        assert "offline=False" in str(failure.value)
        assert requests_made(offline) == []

    def test_the_refusal_names_what_was_missing(self, tmp_path: Path) -> None:
        offline = client_for(DEFAULT_ROUTES, cache_dir=tmp_path, offline=True)
        with pytest.raises(EnsemblUnavailableError, match="lookup/symbol"):
            offline.gene("SCN1A")

    def test_offline_without_a_cache_refuses_everything(self) -> None:
        # No cache and no network is not a usable configuration, and saying so beats an
        # empty result that looks like a gene with no transcripts.
        offline = client_for(DEFAULT_ROUTES, offline=True)
        with pytest.raises(EnsemblUnavailableError, match="no cache"):
            offline.gene("SCN1A")

    def test_the_environment_supplies_the_default(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, "1")
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        with pytest.raises(EnsemblUnavailableError, match="SARNA_OFFLINE"):
            client.gene("SCN1A")

    def test_an_explicit_argument_overrides_the_environment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # A caller who says False means it, whatever the ambient environment says.
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, "1")
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path, offline=False)
        assert client.gene("SCN1A").symbol == "SCN1A"

    def test_the_mode_is_fixed_when_the_client_is_built(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Re-reading the variable mid-run would let one run be half offline.
        client = client_for(DEFAULT_ROUTES, cache_dir=tmp_path)
        monkeypatch.setenv(settings.OFFLINE_VARIABLE, "1")
        assert client.gene("SCN1A").symbol == "SCN1A"

    def test_array_endpoints_refuse_too(self, tmp_path: Path) -> None:
        # Regulatory features and motifs come back as arrays through a separate fetch path.
        # Guarding only the object path would leave annotation quietly reaching the network.
        interval = GenomicInterval(Assembly.GRCH38, "2", 166_127_500, 166_129_000)
        offline = client_for(DEFAULT_ROUTES, cache_dir=tmp_path, offline=True)
        with pytest.raises(EnsemblCacheMissError, match="Cache-only mode"):
            offline.overlap(interval, "regulatory")
        assert requests_made(offline) == []


class TestMissingRecordsVersusOutages:
    """Ensembl answers an unknown symbol with 400, not 404 (found by the live tests).

    Conflating the two would tell a user the service is down when they mistyped a gene, and
    would retry five more times before saying so.
    """


    def test_a_400_saying_no_valid_lookup_is_a_missing_gene(self) -> None:
        client = client_returning(400, '{"error":"No valid lookup found for symbol NOPE"}')
        with pytest.raises(GeneNotFoundError, match="spelling and case"):
            client.gene("NOPE")

    def test_a_missing_gene_is_not_retried(self) -> None:
        client = client_returning(400, '{"error":"No valid lookup found for symbol NOPE"}')
        with pytest.raises(GeneNotFoundError):
            client.gene("NOPE")
        assert len(requests_made(client)) == 1

    def test_a_404_is_still_a_missing_gene(self) -> None:
        with pytest.raises(GeneNotFoundError):
            client_returning(404, "not found").gene("NOPE")

    def test_a_400_for_another_reason_remains_an_outage(self) -> None:
        # A malformed request is the caller's bug, not a missing gene, and must not be dressed
        # up as one.
        client = client_returning(400, '{"error":"Invalid coordinate system"}')
        with pytest.raises(EnsemblUnavailableError, match="HTTP 400"):
            client.gene("SCN1A")

    def test_the_failure_message_reports_the_attempts_actually_made(self) -> None:
        # A 4xx stops immediately, so claiming six attempts would be false.
        client = client_returning(400, '{"error":"Invalid coordinate system"}')
        with pytest.raises(EnsemblUnavailableError, match="after 1 attempt;"):
            client.gene("SCN1A")

    def test_a_server_error_does_exhaust_the_retries(self) -> None:
        client = client_returning(503, "<html>unavailable</html>")
        with pytest.raises(EnsemblUnavailableError, match="after 6 attempts"):
            client.gene("SCN1A")


class TestManeDetectionShapes:
    """Ensembl has reported MANE membership in more than one shape; both are accepted."""

    def test_a_top_level_is_mane_flag_is_honoured(self) -> None:
        flagged = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {
                    **{k: v for k, v in TRANSCRIPTS_PAYLOAD["Transcript"][0].items()
                       if k != "MANE"},
                    "is_mane": 1,
                }
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: flagged})
        assert client.select_transcript(client.gene("SCN1A")).is_mane_select

    def test_a_mane_annotation_without_a_refseq_match_is_tolerated(self) -> None:
        # MANE Select without the accession is still MANE Select; the accession is provenance,
        # not the identity.
        no_accession = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {
                    **TRANSCRIPTS_PAYLOAD["Transcript"][0],
                    "MANE": [{"type": "MANE_Select"}],
                }
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: no_accession})
        transcript = client.select_transcript(client.gene("SCN1A"))
        assert transcript.is_mane_select
        assert transcript.refseq_match is None

    def test_a_mane_plus_clinical_annotation_is_not_mane_select(self) -> None:
        # MANE Plus Clinical is a different designation and must not be mistaken for Select.
        plus_clinical = {
            **TRANSCRIPTS_PAYLOAD,
            "Transcript": [
                {
                    **TRANSCRIPTS_PAYLOAD["Transcript"][0],
                    "MANE": [{"type": "MANE_Plus_Clinical", "refseq_match": "NM_000333.1"}],
                }
            ],
        }
        client = client_for({**DEFAULT_ROUTES, TRANSCRIPTS_ROUTE: plus_clinical})
        with pytest.raises(NoManeSelectTranscriptError):
            client.select_transcript(client.gene("SCN1A"))


class TestTheWorkflowAHelper:
    def test_it_returns_the_gene_transcript_and_promoter_together(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # A promoter without the transcript it came from cannot be audited, so all three are
        # returned rather than just the promoter.
        stub = StubTransport({**DEFAULT_ROUTES, "/sequence/region/": {"seq": "ACGT" * 625}})
        monkeypatch.setattr(
            httpx, "get", lambda url, **_kwargs: stub.get(url, headers={}, timeout=0)
        )
        gene, transcript, promoter = promoter_for_symbol("SCN1A")
        assert gene.symbol == "SCN1A"
        assert transcript.is_mane_select
        assert promoter.length == 2500


class TestTheRealTransportPath:
    def test_without_an_injected_client_httpx_is_used(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The production path. Worth one test, because every other test bypasses it.
        calls: list[str] = []

        def fake_get(url: str, **kwargs: object) -> httpx.Response:
            calls.append(url)
            return httpx.Response(200, json=GENE_PAYLOAD)

        monkeypatch.setattr(httpx, "get", fake_get)
        gene = EnsemblClient(backoff_seconds=0.0).gene("SCN1A")
        assert gene.symbol == "SCN1A"
        assert calls and calls[0].startswith("https://rest.ensembl.org/")

    def test_a_positive_backoff_actually_sleeps(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        slept: list[float] = []
        monkeypatch.setattr("sarna.genome.ensembl.time.sleep", slept.append)
        with pytest.raises(EnsemblUnavailableError):
            client_returning(503, "down").gene("SCN1A")  # backoff 0, so no sleep
        assert not slept
        client = EnsemblClient(
            client=FixedTransport(503, "down"),  # type: ignore[arg-type]
            backoff_seconds=1.0,
        )
        with pytest.raises(EnsemblUnavailableError):
            client.gene("SCN1A")
        assert slept == [1.0, 2.0, 4.0, 8.0, 16.0]


class TestOverlapQueries:
    """The endpoint the regulatory layer is built on."""

    REGULATORY: ClassVar[list[dict[str, object]]] = [
        {"id": "ENSR1", "description": "promoter", "seq_region_name": "2",
         "start": 166128005, "end": 166128231},
    ]

    def test_records_are_returned_unparsed(self) -> None:
        # Interpreting a regulatory feature is the regulatory layer's job; this client's job is
        # to be the only place that knows how to talk to Ensembl.
        client = client_for({"/overlap/region/": self.REGULATORY})
        interval = GenomicInterval(Assembly.GRCH38, "2", 166127500, 166129000)
        records = client.overlap(interval, "regulatory")
        assert records == self.REGULATORY

    def test_the_query_converts_to_1_based_inclusive(self) -> None:
        client = client_for({"/overlap/region/": self.REGULATORY})
        interval = GenomicInterval(Assembly.GRCH38, "2", 1000, 2000)
        client.overlap(interval, "regulatory")
        assert any("2:1001-2000" in url for url in requests_made(client))

    def test_an_empty_result_is_returned_rather_than_an_error(self) -> None:
        # "Nothing annotated" is an answer, and the regulatory layer attaches the caveat.
        client = client_for({"/overlap/region/": []})
        interval = GenomicInterval(Assembly.GRCH38, "2", 1000, 2000)
        assert client.overlap(interval, "regulatory") == []

    def test_an_interval_on_another_assembly_is_refused(self) -> None:
        client = client_for({"/overlap/region/": []})
        interval = GenomicInterval(Assembly.GRCH37, "2", 1000, 2000)
        with pytest.raises(ValueError, match="this client serves GRCh38"):
            client.overlap(interval, "regulatory")

    def test_an_object_where_an_array_was_expected_is_an_error(self) -> None:
        client = client_for({"/overlap/region/": {"error": "unexpected"}})
        interval = GenomicInterval(Assembly.GRCH38, "2", 1000, 2000)
        with pytest.raises(EnsemblUnavailableError, match="expected a JSON array"):
            client.overlap(interval, "regulatory")

    def test_an_array_where_an_object_was_expected_is_an_error(self) -> None:
        client = client_for({"/lookup/symbol": [{"unexpected": True}]})
        with pytest.raises(EnsemblUnavailableError, match="expected a JSON object"):
            client.gene("SCN1A")

    def test_overlap_results_are_cached(self, tmp_path: Path) -> None:
        client = client_for({"/overlap/region/": self.REGULATORY}, cache_dir=tmp_path)
        interval = GenomicInterval(Assembly.GRCH38, "2", 1000, 2000)
        client.overlap(interval, "regulatory")
        before = len(requests_made(client))
        assert client.overlap(interval, "regulatory") == self.REGULATORY
        assert len(requests_made(client)) == before


@pytest.mark.strand
class TestDesigningAgainstADeclaredPromoter:
    """Anchors are why chapter 2's SCN1A designs are reproducible at all."""

    SEQUENCE: ClassVar[dict[str, str]] = {"seq": AUTO_SEQUENCE}

    def client(self) -> EnsemblClient:
        return client_for({**DEFAULT_ROUTES, "/sequence/region/": self.SEQUENCE})

    @staticmethod
    def declared(reference: int | None, kind: ReferenceKind) -> PromoterAnchor:
        return PromoterAnchor(
            id="scn1a-test",
            gene="SCN1A",
            interval=GenomicInterval(
                Assembly.GRCH38, "2", 166_078_692, 166_080_192, Strand.MINUS
            ),
            reference_point=reference,
            reference_kind=kind,
            provenance=Provenance.citation("Dong et al. 2014", "Mol Neurobiol 50(2):438-448"),
            note="A declared promoter with no annotated transcript nearby.",
        )

    def test_a_declared_promoter_yields_a_tss_anchored_window(self) -> None:
        anchor = self.declared(166_079_092, ReferenceKind.INFERRED)
        promoter = self.client().promoter_for_anchor(anchor, 600, 600)
        assert promoter.tss_start is not None
        assert promoter.tss_start.value == -600
        assert promoter.interval is not None
        assert promoter.interval.end == 166_079_092 + 600 + 1

    def test_an_anchor_with_no_reference_point_yields_a_promoter_with_no_tss(self) -> None:
        # The honest degradation: TSS rules then report not-applicable rather than measuring
        # from an origin nobody established.
        anchor = self.declared(None, ReferenceKind.NONE)
        promoter = self.client().promoter_for_anchor(anchor, region_padding_bp=0)
        assert promoter.tss_start is None
        assert promoter.length == 1500

    def test_a_minus_strand_declared_promoter_is_reverse_complemented(self) -> None:
        anchor = self.declared(166_079_092, ReferenceKind.INFERRED)
        promoter = self.client().promoter_for_anchor(anchor, 600, 600)
        assert promoter.sequence.bases.startswith("GTTT")

    def test_an_anchor_on_another_assembly_is_refused(self) -> None:
        anchor = PromoterAnchor(
            id="wrong-assembly",
            gene="SCN1A",
            interval=GenomicInterval(Assembly.GRCH37, "2", 1000, 2000, Strand.MINUS),
            reference_point=1500,
            reference_kind=ReferenceKind.ANNOTATED_TSS,
            provenance=Provenance.citation("x", "y"),
            note="n",
        )
        with pytest.raises(ValueError, match="this client serves GRCh38"):
            self.client().promoter_for_anchor(anchor)
