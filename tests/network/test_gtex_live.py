"""Live GTEx tests, excluded from the default run.

Recorded payloads cannot catch the failure that matters most here: GTEx publishing a new
dataset, retiring an old one, or changing which GENCODE release a dataset is keyed on. That
last one is not hypothetical, it is why :data:`sarna.expression.gtex.DATASET_ANNOTATION`
exists, and the symptom of getting it wrong is an empty result with no error, which reads as
"this gene is not expressed".

Run with::

    .venv/bin/python -m pytest -m network

Failures here mean "check GTEx", not "the engine is broken".
"""

import pytest

from sarna.expression.gtex import (
    DATASET_ANNOTATION,
    GeneNotInGtexError,
    GtexClient,
)
from sarna.expression.models import ExpressionStatus

pytestmark = pytest.mark.network

SCN1A = "ENSG00000144285"

SCN1A_MANE = "ENST00000674923"
"""SCN1A's MANE Select transcript, verified against Ensembl's record for the gene."""

SCN1A_V10_ID = "ENSG00000144285.23"
"""What gtex_v10 calls SCN1A. Verified against the portal 2026-10-08."""


@pytest.fixture(scope="module")
def client() -> GtexClient:
    return GtexClient()


class TestTheReleaseMappingIsStillTrue:
    def test_v10_still_resolves_scn1a_through_gencode_v39(self, client: GtexClient) -> None:
        # The regression test for the versioning trap. If GTEx rebuilds v10 on another
        # GENCODE release, this fails here rather than silently returning no expression.
        assert client.annotation_release == "v39"
        assert client.release_gene_id(SCN1A) == SCN1A_V10_ID

    @pytest.mark.parametrize("dataset_id", sorted(DATASET_ANNOTATION))
    def test_every_dataset_this_module_claims_to_know_still_resolves_the_gene(
        self, dataset_id: str
    ) -> None:
        resolved = GtexClient(dataset_id=dataset_id).release_gene_id(SCN1A)
        assert resolved.startswith(f"{SCN1A}.")

    def test_the_two_datasets_disagree_about_the_version(self) -> None:
        # The whole reason the mapping is explicit. If these ever matched, the trap would be
        # gone and this test should be reconsidered rather than deleted.
        v8 = GtexClient(dataset_id="gtex_v8").release_gene_id(SCN1A)
        v10 = GtexClient(dataset_id="gtex_v10").release_gene_id(SCN1A)
        assert v8 != v10


class TestTheLiveProfile:
    def test_it_covers_the_thirteen_brain_regions(self, client: GtexClient) -> None:
        profile = client.expression_for(SCN1A)
        brain = {
            median.tissue_id
            for entry in profile.measured_transcripts
            for median in entry.medians
            if median.tissue_id.startswith("Brain_")
        }
        assert len(brain) == 13

    def test_every_median_carries_a_readable_name_and_a_unit(
        self, client: GtexClient
    ) -> None:
        profile = client.expression_for(SCN1A)
        medians = [
            median for entry in profile.measured_transcripts for median in entry.medians
        ]
        assert medians
        assert all(median.tissue_name and median.unit for median in medians)

    def test_most_tissues_report_how_many_donors_are_behind_the_median(
        self, client: GtexClient
    ) -> None:
        # Measured rather than asserted at 100%: GTEx does not publish a donor count for
        # every group, and a missing count is recorded as zero rather than invented.
        profile = client.expression_for(SCN1A)
        counts = [
            median.sample_count
            for entry in profile.measured_transcripts
            for median in entry.medians
        ]
        assert sum(1 for count in counts if count > 0) / len(counts) > 0.9

    def test_the_mane_transcript_is_measured_and_is_zero_in_the_brain(
        self, client: GtexClient
    ) -> None:
        # The finding that motivates showing expression beside the transcript picker at all.
        # If GTEx ever revises this, the app's wording should be revisited, so it fails here.
        entry = client.expression_for(SCN1A).for_transcript(SCN1A_MANE)
        assert entry.status is ExpressionStatus.MEASURED
        cortex = entry.median_in("Brain_Cortex")
        assert cortex is not None
        assert cortex.median == 0.0

    def test_some_other_transcript_is_expressed_in_the_cortex(
        self, client: GtexClient
    ) -> None:
        profile = client.expression_for(SCN1A)
        expressed = [
            entry.transcript_id
            for entry in profile.measured_transcripts
            if (median := entry.median_in("Brain_Cortex")) is not None and median.median > 1
        ]
        assert expressed
        assert SCN1A_MANE not in expressed


class TestLiveFailureModes:
    def test_a_gene_gtex_has_never_heard_of_raises_rather_than_reporting_zeros(self) -> None:
        # Live because the behaviour under test is GTEx's: it answers HTTP 200 with an empty
        # data list, which a recorded payload can only imitate.
        client = GtexClient()
        with pytest.raises(GeneNotInGtexError):
            client.release_gene_id("ENSG00000000000")
