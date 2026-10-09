"""Tests for the expression records, concentrating on the three silences they keep apart.

The whole point of these types is that "measured zero", "absent from this release" and "this
resource does not cover that sample" are different statements. A record that lets the first two
collapse into one number is the defect, so most of these tests are about the status rather than
about the value.
"""

import pytest

from sarna.expression.models import (
    ExpressionStatus,
    GeneExpression,
    TissueMedian,
    TranscriptExpression,
)
from sarna.rules.provenance import Provenance

CORTEX = TissueMedian(
    tissue_id="Brain_Cortex", tissue_name="Brain - Cortex", median=1.63, unit="TPM",
    ontology_id="UBERON:0001870", sample_count=270,
)

NIGRA = TissueMedian(
    tissue_id="Brain_Substantia_nigra", tissue_name="Brain - Substantia nigra",
    median=0.52, unit="TPM", ontology_id="UBERON:0002038", sample_count=147,
)

BLOOD = TissueMedian(
    tissue_id="Whole_Blood", tissue_name="Whole Blood", median=0.0, unit="TPM",
)


def profile(*transcripts: TranscriptExpression) -> GeneExpression:
    """Build a profile around the given transcript records."""
    return GeneExpression(
        gene_id="ENSG00000144285", release_gene_id="ENSG00000144285.23",
        dataset_id="gtex_v10", annotation_release="v39", transcripts=transcripts,
        provenance=Provenance.citation("GTEx v10", "test fixture"),
    )


class TestAStatusQualifiesEveryNumber:
    def test_a_measured_zero_says_it_was_measured(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000674923", status=ExpressionStatus.MEASURED,
            medians=(BLOOD,),
        )
        assert entry.status.is_measured
        assert entry.medians[0].median == 0.0
        assert "median of zero is a measurement" in entry.status.note

    def test_absence_from_the_release_is_not_measured_and_says_why(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000713692", status=ExpressionStatus.NOT_IN_RELEASE,
        )
        assert not entry.status.is_measured
        assert entry.medians == ()
        assert "never quantified" in entry.status.note
        assert "not a measurement of zero" in entry.status.note

    def test_the_two_notes_are_different_text(self) -> None:
        # If they were ever made the same string, an interface would render ignorance and a
        # finding identically, which is the error this module exists to prevent.
        assert ExpressionStatus.MEASURED.note != ExpressionStatus.NOT_IN_RELEASE.note

    def test_a_record_cannot_be_unmeasured_and_carry_a_number(self) -> None:
        with pytest.raises(ValueError, match="measurement that was never made"):
            TranscriptExpression(
                transcript_id="ENST00000713692",
                status=ExpressionStatus.NOT_IN_RELEASE,
                medians=(CORTEX,),
            )

    def test_the_status_renders_as_its_label(self) -> None:
        assert str(ExpressionStatus.MEASURED) == "measured"
        assert str(ExpressionStatus.NOT_IN_RELEASE) == "not-in-release"


class TestReadingOneTranscript:
    def test_highest_returns_the_largest_median_not_the_first(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
            medians=(BLOOD, CORTEX, NIGRA),
        )
        highest = entry.highest
        assert highest is not None
        assert highest.tissue_id == "Brain_Cortex"

    def test_highest_is_none_when_nothing_was_measured(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000713692", status=ExpressionStatus.NOT_IN_RELEASE,
        )
        assert entry.highest is None

    def test_a_tissue_that_was_measured_comes_back(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
            medians=(CORTEX, NIGRA),
        )
        found = entry.median_in("Brain_Substantia_nigra")
        assert found is not None
        assert found.median == 0.52

    def test_a_tissue_that_was_not_measured_is_none_rather_than_zero(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
            medians=(CORTEX,),
        )
        assert entry.median_in("Liver") is None


class TestAMedianCarriesItsQualifiers:
    def test_the_donor_count_travels_with_the_number(self) -> None:
        # 270 donors for cortex against 147 for substantia nigra: the same printed precision
        # does not mean the same certainty.
        assert CORTEX.sample_count == 270
        assert NIGRA.sample_count == 147

    def test_an_unreported_donor_count_is_zero_rather_than_a_guess(self) -> None:
        assert BLOOD.sample_count == 0

    def test_the_unit_is_whatever_the_resource_said(self) -> None:
        assert CORTEX.unit == "TPM"

    def test_it_serialises_every_field_it_holds(self) -> None:
        assert CORTEX.to_dict() == {
            "tissue_id": "Brain_Cortex", "tissue_name": "Brain - Cortex", "median": 1.63,
            "unit": "TPM", "ontology_id": "UBERON:0001870", "sample_count": 270,
        }


class TestReadingAGene:
    def test_an_unknown_transcript_comes_back_qualified_rather_than_as_none(self) -> None:
        # Returning None would let a caller write `or 0`. Returning a record forces the
        # status to be carried alongside.
        entry = profile().for_transcript("ENST00000713692")
        assert entry.status is ExpressionStatus.NOT_IN_RELEASE
        assert entry.transcript_id == "ENST00000713692"

    def test_a_known_transcript_comes_back_as_stored(self) -> None:
        stored = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
            medians=(CORTEX,), release_id="ENST00000303395.9",
        )
        assert profile(stored).for_transcript("ENST00000303395") is stored

    def test_measured_transcripts_exclude_the_ones_never_quantified(self) -> None:
        measured = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
            medians=(CORTEX,),
        )
        absent = TranscriptExpression(
            transcript_id="ENST00000713692", status=ExpressionStatus.NOT_IN_RELEASE,
        )
        assert profile(measured, absent).measured_transcripts == (measured,)

    def test_it_serialises_its_dataset_release_and_caveats(self) -> None:
        payload = GeneExpression(
            gene_id="ENSG00000144285", release_gene_id="ENSG00000144285.23",
            dataset_id="gtex_v10", annotation_release="v39",
            transcripts=(
                TranscriptExpression(
                    transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
                    medians=(CORTEX,), release_id="ENST00000303395.9",
                ),
            ),
            provenance=Provenance.citation("GTEx v10", "test fixture"),
            caveats=("bulk tissue",),
        ).to_dict()
        assert payload["dataset_id"] == "gtex_v10"
        assert payload["annotation_release"] == "v39"
        assert payload["caveats"] == ["bulk tissue"]
        transcripts = payload["transcripts"]
        assert isinstance(transcripts, list)
        assert transcripts[0]["status"] == "measured"
        assert transcripts[0]["release_id"] == "ENST00000303395.9"
        assert "status_note" in transcripts[0]


class TestTheRecordsAreImmutable:
    def test_a_median_cannot_be_rewritten(self) -> None:
        with pytest.raises(AttributeError):
            CORTEX.median = 99.0  # type: ignore[misc]

    def test_a_transcript_record_cannot_be_rewritten(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED,
        )
        with pytest.raises(AttributeError):
            entry.status = ExpressionStatus.NOT_IN_RELEASE  # type: ignore[misc]
