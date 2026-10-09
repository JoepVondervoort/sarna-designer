"""Tests for the GTEx client, run entirely offline against recorded payloads.

The payloads are trimmed copies of real ``gtex_v10`` responses for SCN1A, retrieved on
8 October 2026. SCN1A is the right fixture for the same reason it is everywhere else in this
suite, plus one specific to this module: its MANE Select transcript has a median of 0.00 TPM
in every GTEx tissue while two of its other transcripts carry the brain signal, so a client
that quietly confused "the preferred transcript" with "the expressed transcript" would be
caught by a test rather than by a researcher.

A live counterpart lives in ``tests/network/test_gtex_live.py``, excluded by default.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import httpx
import pytest

from sarna.expression.gtex import (
    DATASET_ANNOTATION,
    GeneNotInGtexError,
    GtexClient,
    GtexResponseError,
    UnknownGtexDatasetError,
)
from sarna.expression.models import ExpressionStatus
from sarna.sources.json_http import SourceCacheMissError

MANE = "ENST00000674923"
"""SCN1A's MANE Select transcript, confirmed against Ensembl's own record for the gene."""

CORTEX_LEADER = "ENST00000303395"
"""The transcript that actually carries SCN1A's cortex signal: 1.63 TPM against the MANE
transcript's 0.00."""

GENE_PAYLOAD: dict[str, Any] = {
    "data": [
        {
            "chromosome": "chr2",
            "description": "sodium voltage-gated channel alpha subunit 1",
            "gencodeId": "ENSG00000144285.23",
            "gencodeVersion": "v39",
            "geneSymbol": "SCN1A",
            "genomeBuild": "GRCh38/hg38",
        }
    ],
    "paging_info": {"numberOfPages": 1, "page": 0, "totalNumberOfItems": 1},
}

TISSUE_PAYLOAD: dict[str, Any] = {
    "data": [
        {
            "tissueSiteDetailId": "Brain_Cortex",
            "tissueSiteDetail": "Brain - Cortex",
            "ontologyId": "UBERON:0001870",
            "datasetId": "gtex_v10",
            "rnaSeqSampleSummary": {"totalCount": 270},
        },
        {
            "tissueSiteDetailId": "Brain_Substantia_nigra",
            "tissueSiteDetail": "Brain - Substantia nigra",
            "ontologyId": "UBERON:0002038",
            "datasetId": "gtex_v10",
            "rnaSeqSampleSummary": {"totalCount": 147},
        },
        {
            "tissueSiteDetailId": "Whole_Blood",
            "tissueSiteDetail": "Whole Blood",
            "ontologyId": "UBERON:0013756",
            "datasetId": "gtex_v10",
        },
    ],
    "paging_info": {"numberOfPages": 1, "page": 0, "totalNumberOfItems": 3},
}


def median_row(transcript: str, tissue: str, median: float) -> dict[str, Any]:
    """Build one median row in the shape GTEx returns."""
    return {
        "median": median,
        "transcriptId": transcript,
        "tissueSiteDetailId": tissue,
        "ontologyId": "UBERON:0001870",
        "datasetId": "gtex_v10",
        "gencodeId": "ENSG00000144285.23",
        "geneSymbol": "SCN1A",
        "unit": "TPM",
    }


MEDIAN_ROWS: list[dict[str, Any]] = [
    median_row(f"{MANE}.1", "Brain_Cortex", 0.0),
    median_row(f"{MANE}.1", "Brain_Substantia_nigra", 0.0),
    median_row(f"{MANE}.1", "Whole_Blood", 0.0),
    median_row(f"{CORTEX_LEADER}.9", "Brain_Cortex", 1.63),
    median_row(f"{CORTEX_LEADER}.9", "Brain_Substantia_nigra", 0.52),
    median_row(f"{CORTEX_LEADER}.9", "Whole_Blood", 0.0),
    median_row("ENST00000375405.7", "Brain_Cortex", 0.02),
    median_row("ENST00000375405.7", "Brain_Substantia_nigra", 0.0),
    median_row("ENST00000375405.7", "Whole_Blood", 0.0),
]
"""Nine rows: three transcripts across three tissues, values as the portal reports them."""

MEDIAN_PAYLOAD: dict[str, Any] = {
    "data": MEDIAN_ROWS,
    "paging_info": {
        "numberOfPages": 1, "page": 0, "maxItemsPerPage": 250, "totalNumberOfItems": 9,
    },
}


class StubTransport:
    """Serves recorded payloads by URL fragment, and records what was asked for.

    An ``int`` value stands for a status code with an HTML body, which is how an overloaded
    origin reports a fault.
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
                return httpx.Response(200, json=payload)
        return httpx.Response(404, json={"detail": "not found"})


DEFAULT_ROUTES: dict[str, object] = {
    "reference/gene": GENE_PAYLOAD,
    "tissueSiteDetail": TISSUE_PAYLOAD,
    "medianTranscriptExpression": MEDIAN_PAYLOAD,
}


def client_for(
    routes: Mapping[str, object] | None = None,
    cache_dir: Path | None = None,
    offline: bool | None = False,
    dataset_id: str = "gtex_v10",
) -> tuple[GtexClient, StubTransport]:
    """Build a client backed by a stub, with retries made instant."""
    transport = StubTransport(DEFAULT_ROUTES if routes is None else routes)
    client = GtexClient(
        dataset_id=dataset_id,
        client=transport,  # type: ignore[arg-type]
        backoff_seconds=0.0,
        cache_dir=cache_dir,
        offline=offline,
    )
    return client, transport


class TestTheDatasetDecidesTheAnnotationRelease:
    def test_each_known_dataset_names_the_gencode_release_it_was_built_on(self) -> None:
        # Not cosmetic: querying gtex_v10 with a v26 identifier returns an empty result and
        # no error, which reads exactly like "this gene is not expressed".
        assert DATASET_ANNOTATION["gtex_v8"] == "v26"
        assert DATASET_ANNOTATION["gtex_v10"] == "v39"

    def test_an_unknown_dataset_is_refused_rather_than_guessed(self) -> None:
        with pytest.raises(UnknownGtexDatasetError, match="annotation release"):
            GtexClient(dataset_id="gtex_v11")

    def test_the_release_is_sent_with_the_gene_lookup(self) -> None:
        client, transport = client_for()
        client.release_gene_id("ENSG00000144285")
        assert "gencodeVersion=v39" in transport.requested[0]

    def test_a_version_the_caller_supplies_is_discarded(self) -> None:
        # Ensembl says .26 today; gtex_v10 holds .23. Sending the caller's version would
        # return nothing.
        client, transport = client_for()
        assert client.release_gene_id("ENSG00000144285.26") == "ENSG00000144285.23"
        assert "ENSG00000144285&" in transport.requested[0]

    def test_an_answer_about_a_different_release_is_rejected(self) -> None:
        wrong = {"data": [{"gencodeId": "ENSG00000144285.17", "gencodeVersion": "v26"}]}
        client, _ = client_for({**DEFAULT_ROUTES, "reference/gene": wrong})
        with pytest.raises(GtexResponseError, match="two releases must not be mixed"):
            client.release_gene_id("ENSG00000144285")


class TestAnEmptyResultIsNotAZero:
    def test_a_gene_absent_from_the_release_raises_rather_than_reporting_zeros(self) -> None:
        # GTEx answers HTTP 200 with an empty data list for a gene it has never heard of.
        empty: dict[str, Any] = {"data": [], "paging_info": {"totalNumberOfItems": 0}}
        client, _ = client_for({**DEFAULT_ROUTES, "reference/gene": empty})
        with pytest.raises(GeneNotInGtexError, match="not a measurement of zero"):
            client.expression_for("ENSG00000000000")

    def test_a_gene_the_release_does_not_match_raises(self) -> None:
        other = {"data": [{"gencodeId": "ENSG00000999999.1", "gencodeVersion": "v39"}]}
        client, _ = client_for({**DEFAULT_ROUTES, "reference/gene": other})
        with pytest.raises(GeneNotInGtexError):
            client.release_gene_id("ENSG00000144285")


class TestTheProfileItBuilds:
    def test_every_transcript_gtex_quantified_comes_back_measured(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285")
        assert len(profile.transcripts) == 3
        assert all(entry.status is ExpressionStatus.MEASURED for entry in profile.transcripts)

    def test_it_carries_the_dataset_and_the_release_with_the_numbers(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285")
        assert profile.dataset_id == "gtex_v10"
        assert profile.annotation_release == "v39"
        assert profile.release_gene_id == "ENSG00000144285.23"
        assert profile.gene_id == "ENSG00000144285"

    def test_the_versioned_identifier_is_kept_beside_the_unversioned_one(self) -> None:
        client, _ = client_for()
        entry = client.expression_for("ENSG00000144285").for_transcript(MANE)
        assert entry.transcript_id == MANE
        assert entry.release_id == f"{MANE}.1"

    def test_readable_tissue_names_and_donor_counts_are_joined_on(self) -> None:
        # The expression endpoint reports only Brain_Cortex; the name and the 270 donors
        # behind the median come from the tissue endpoint.
        client, _ = client_for()
        cortex = client.expression_for("ENSG00000144285").for_transcript(
            CORTEX_LEADER
        ).median_in("Brain_Cortex")
        assert cortex is not None
        assert cortex.tissue_name == "Brain - Cortex"
        assert cortex.ontology_id == "UBERON:0001870"
        assert cortex.sample_count == 270

    def test_a_tissue_that_reports_no_donor_count_says_zero_rather_than_inventing_one(
        self,
    ) -> None:
        client, _ = client_for()
        blood = client.expression_for("ENSG00000144285").for_transcript(MANE).median_in(
            "Whole_Blood"
        )
        assert blood is not None
        assert blood.sample_count == 0

    def test_the_medians_are_the_values_gtex_reported(self) -> None:
        # Hand-checked against the payload above rather than recomputed by the same code.
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285")
        leader = profile.for_transcript(CORTEX_LEADER)
        assert [entry.median for entry in leader.medians] == [1.63, 0.52, 0.0]
        assert {entry.unit for entry in leader.medians} == {"TPM"}

    def test_the_highest_tissue_is_the_largest_median_not_the_first_row(self) -> None:
        client, _ = client_for()
        highest = client.expression_for("ENSG00000144285").for_transcript(
            CORTEX_LEADER
        ).highest
        assert highest is not None
        assert highest.tissue_id == "Brain_Cortex"
        assert highest.median == 1.63


class TestTheManeTranscriptIsNotTheExpressedOne:
    def test_scn1a_mane_select_is_measured_at_zero_in_every_tissue(self) -> None:
        # The finding this module exists for. A measured zero, not a gap: the status says
        # MEASURED, so an interface must not render it as "no data".
        client, _ = client_for()
        entry = client.expression_for("ENSG00000144285").for_transcript(MANE)
        assert entry.status is ExpressionStatus.MEASURED
        assert [median.median for median in entry.medians] == [0.0, 0.0, 0.0]

    def test_another_transcript_carries_the_cortex_signal(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285")
        in_cortex = {
            entry.transcript_id: median.median
            for entry in profile.measured_transcripts
            if (median := entry.median_in("Brain_Cortex")) is not None
        }
        assert in_cortex == {MANE: 0.0, CORTEX_LEADER: 1.63, "ENST00000375405": 0.02}
        assert max(in_cortex, key=lambda key: in_cortex[key]) == CORTEX_LEADER


class TestTranscriptsTheReleaseNeverSaw:
    def test_a_requested_transcript_absent_from_the_release_is_recorded_as_such(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285", ["ENST00000713692", MANE])
        entry = profile.for_transcript("ENST00000713692")
        assert entry.status is ExpressionStatus.NOT_IN_RELEASE
        assert entry.medians == ()
        assert "never quantified" in entry.status.note

    def test_the_gap_appears_in_the_profile_rather_than_only_on_request(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285", ["ENST00000713692"])
        absent = [
            entry for entry in profile.transcripts
            if entry.status is ExpressionStatus.NOT_IN_RELEASE
        ]
        assert [entry.transcript_id for entry in absent] == ["ENST00000713692"]

    def test_the_caveat_counts_them_and_refuses_the_word_zero(self) -> None:
        client, _ = client_for()
        profile = client.expression_for(
            "ENSG00000144285", ["ENST00000713692", "ENST00000714447"]
        )
        gap = [caveat for caveat in profile.caveats if "GENCODE v39" in caveat]
        assert len(gap) == 1
        assert gap[0].startswith("2 requested transcripts are absent")
        assert "not a measurement of zero" in gap[0]

    def test_a_requested_transcript_is_not_counted_twice_or_versioned(self) -> None:
        client, _ = client_for()
        profile = client.expression_for(
            "ENSG00000144285", ["ENST00000713692.1", "ENST00000713692"]
        )
        absent = [
            entry for entry in profile.transcripts
            if entry.status is ExpressionStatus.NOT_IN_RELEASE
        ]
        assert [entry.transcript_id for entry in absent] == ["ENST00000713692"]

    def test_a_transcript_gtex_holds_but_ensembl_has_retired_is_still_reported(self) -> None:
        # The disagreement runs both ways: ENST00000375405 is quantified in GENCODE v39 and
        # no longer served by Ensembl. Dropping it would hide a measurement.
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285", [MANE, CORTEX_LEADER])
        assert "ENST00000375405" in {
            entry.transcript_id for entry in profile.measured_transcripts
        }

    def test_no_gap_caveat_when_every_requested_transcript_was_quantified(self) -> None:
        client, _ = client_for()
        profile = client.expression_for("ENSG00000144285", [MANE, CORTEX_LEADER])
        assert not [caveat for caveat in profile.caveats if "GENCODE v39" in caveat]


class TestTheCaveatsItAlwaysCarries:
    def test_it_says_gtex_has_no_cell_lines(self) -> None:
        # The researcher's question was about AF22. The answer travels with the profile rather
        # than living only in the page, so a methods export carries it too.
        client, _ = client_for()
        caveats = client.expression_for("ENSG00000144285").caveats
        assert "AF22" in caveats[0]
        assert "no cell lines" in caveats[0]

    def test_it_says_bulk_tissue_averages_over_cell_types(self) -> None:
        # Distinct from the line above, and the more easily missed of the two: a transcript
        # used by few cells reads as near zero in the tissue those cells sit in.
        client, _ = client_for()
        caveats = " ".join(client.expression_for("ENSG00000144285").caveats)
        assert "average of many cell types" in caveats
        assert "can read as near zero" in caveats

    def test_it_says_a_median_hides_the_spread_and_the_donor_counts_differ(self) -> None:
        client, _ = client_for()
        caveats = " ".join(client.expression_for("ENSG00000144285").caveats)
        assert "median across donors" in caveats
        assert "not equally certain" in caveats
        # The measured spread, not a vague "several fold": 11 donors to 818 in gtex_v10.
        assert "74-fold" in caveats

    def test_it_says_isoform_values_are_apportioned_estimates(self) -> None:
        client, _ = client_for()
        caveats = " ".join(client.expression_for("ENSG00000144285").caveats)
        assert "apportioned between overlapping isoforms" in caveats

    def test_the_provenance_names_the_dataset_and_the_release(self) -> None:
        client, _ = client_for()
        provenance = str(client.expression_for("ENSG00000144285").provenance)
        assert "gtex_v10" in provenance
        assert "GENCODE v39" in provenance


class TestPaging:
    def test_every_page_is_fetched_and_concatenated(self) -> None:
        first = {
            "data": MEDIAN_ROWS[:6],
            "paging_info": {"numberOfPages": 2, "page": 0, "totalNumberOfItems": 9},
        }
        second = {
            "data": MEDIAN_ROWS[6:],
            "paging_info": {"numberOfPages": 2, "page": 1, "totalNumberOfItems": 9},
        }
        pages = iter([first, second])

        class PagingTransport(StubTransport):
            def get(
                self, url: str, headers: dict[str, str], timeout: float
            ) -> httpx.Response:
                if "medianTranscriptExpression" in url:
                    self.requested.append(url)
                    return httpx.Response(200, json=next(pages))
                return super().get(url, headers, timeout)

        transport = PagingTransport(DEFAULT_ROUTES)
        client = GtexClient(
            client=transport,  # type: ignore[arg-type]
            backoff_seconds=0.0, offline=False,
        )
        profile = client.expression_for("ENSG00000144285")
        asked = [url for url in transport.requested if "medianTranscriptExpression" in url]
        assert len(asked) == 2
        assert "page=0" in asked[0]
        assert "page=1" in asked[1]
        assert len(profile.transcripts) == 3
        assert len(profile.for_transcript("ENST00000375405").medians) == 3

    def test_a_single_page_is_not_followed_by_a_second_request(self) -> None:
        client, transport = client_for()
        client.expression_for("ENSG00000144285")
        asked = [url for url in transport.requested if "medianTranscriptExpression" in url]
        assert len(asked) == 1


class TestMalformedResponses:
    def test_a_response_without_a_data_list_is_refused(self) -> None:
        client, _ = client_for({**DEFAULT_ROUTES, "tissueSiteDetail": {"oops": True}})
        with pytest.raises(GtexResponseError, match="no 'data' list"):
            client.expression_for("ENSG00000144285")

    def test_a_row_without_a_median_is_refused_rather_than_read_as_zero(self) -> None:
        broken = dict(median_row(f"{MANE}.1", "Brain_Cortex", 0.0))
        del broken["median"]
        client, _ = client_for({
            **DEFAULT_ROUTES,
            "medianTranscriptExpression": {"data": [broken], "paging_info": {}},
        })
        with pytest.raises(GtexResponseError, match="missing its transcript, tissue or value"):
            client.expression_for("ENSG00000144285")

    def test_a_row_about_another_gene_is_refused(self) -> None:
        stray = {**median_row(f"{MANE}.1", "Brain_Cortex", 9.9), "gencodeId": "ENSG00000111.1"}
        client, _ = client_for({
            **DEFAULT_ROUTES,
            "medianTranscriptExpression": {"data": [stray], "paging_info": {}},
        })
        with pytest.raises(GtexResponseError, match="but a row reports"):
            client.expression_for("ENSG00000144285")

    def test_a_dataset_reporting_no_tissues_is_refused(self) -> None:
        client, _ = client_for({
            **DEFAULT_ROUTES, "tissueSiteDetail": {"data": [], "paging_info": {}},
        })
        with pytest.raises(GtexResponseError, match="reported no tissues"):
            client.expression_for("ENSG00000144285")


class TestCachingAndOffline:
    def test_a_second_lookup_is_served_from_disk(self, tmp_path: Path) -> None:
        client, transport = client_for(cache_dir=tmp_path)
        client.expression_for("ENSG00000144285")
        before = len(transport.requested)
        client.expression_for("ENSG00000144285")
        assert len(transport.requested) == before

    def test_an_offline_client_serves_what_was_cached(self, tmp_path: Path) -> None:
        warm, _ = client_for(cache_dir=tmp_path)
        expected = warm.expression_for("ENSG00000144285")
        offline, transport = client_for(cache_dir=tmp_path, offline=True)
        assert offline.expression_for("ENSG00000144285").to_dict() == expected.to_dict()
        assert transport.requested == []

    def test_an_offline_client_refuses_rather_than_returning_nothing(
        self, tmp_path: Path
    ) -> None:
        client, _ = client_for(cache_dir=tmp_path, offline=True)
        with pytest.raises(SourceCacheMissError, match="Cache-only mode"):
            client.expression_for("ENSG00000144285")

    def test_the_offline_flag_is_readable(self) -> None:
        client, _ = client_for(offline=True)
        assert client.offline is True
