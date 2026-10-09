"""Transcript-level expression, and the three different silences it has to tell apart.

A tool that prints "0" beside a transcript is making a claim, and there are three ways a number
can be missing that must not be drawn alike:

1. **Measured, and it is zero.** The resource quantified this transcript in this tissue and the
   median was 0.00 TPM. A real finding, and often the interesting one.
2. **Not in this resource's annotation.** The transcript exists in Ensembl today but is absent
   from the GENCODE release the expression resource was built against, so nothing was ever
   quantified for it. Rendering that as zero invents a measurement.
3. **The resource does not cover this sample at all.** GTEx is bulk post-mortem *tissue*. AF22
   and HeLa are not in it and never will be, so "no GTEx value for AF22" says nothing about
   AF22.

The second case is not hypothetical and not rare. GTEx v8 is built on GENCODE v26, from 2017,
and SCN1A's MANE Select transcript (the one this app offers first) does not appear in it at
all; v10 (GENCODE v39) does contain it. Which release a number came from is therefore carried
with the number rather than assumed.

This mirrors :class:`sarna.regulatory.evidence.CellTypeStatus`, for the same reason: absence of
data and a measured absence are different statements, and an interface that renders them alike
turns ignorance into a finding.
"""

from dataclasses import dataclass
from enum import Enum

from sarna.rules.provenance import Provenance


class ExpressionStatus(Enum):
    """Why a transcript does or does not carry an expression value."""

    MEASURED = "measured"
    """The resource quantified this transcript. The value may still be zero."""

    NOT_IN_RELEASE = "not-in-release"
    """Absent from the annotation release the resource was built on. Nothing was measured."""

    def __str__(self) -> str:
        """Return the status label."""
        return self.value

    @property
    def is_measured(self) -> bool:
        """Return whether a number here came from a measurement."""
        return self is ExpressionStatus.MEASURED

    @property
    def note(self) -> str:
        """Return the caveat an interface must show beside this status."""
        if self is ExpressionStatus.MEASURED:
            return "Quantified in this resource. A median of zero is a measurement, not a gap."
        return (
            "This transcript is not in the annotation release the resource was built on, so it "
            "was never quantified. This is not a measurement of zero, and says nothing about "
            "whether the transcript is used."
        )


@dataclass(frozen=True, slots=True)
class TissueMedian:
    """One sample group's median expression of one transcript."""

    tissue_id: str
    """The resource's own identifier, e.g. ``Brain_Cortex``. Stable; used for lookups."""

    tissue_name: str
    """The readable name, e.g. ``Brain - Cortex``."""

    median: float
    """Median across that group's samples, in :attr:`unit`. A median hides the spread."""

    unit: str
    """Units as the resource states them, e.g. ``TPM``. Never assumed."""

    ontology_id: str = ""
    """UBERON term for the sample site, where the resource supplies one."""

    sample_count: int = 0
    """Donors behind the median, or ``0`` where the resource did not say.

    Carried because a median is only as good as what it summarises, and GTEx's groups differ
    74-fold: 11 donors for kidney medulla against 818 for skeletal muscle. Printed to the same
    two decimal places they look equally certain, and are not.
    """

    def to_dict(self) -> dict[str, object]:
        """Return the median as JSON-serialisable data."""
        return {
            "tissue_id": self.tissue_id,
            "tissue_name": self.tissue_name,
            "median": self.median,
            "unit": self.unit,
            "ontology_id": self.ontology_id,
            "sample_count": self.sample_count,
        }


@dataclass(frozen=True, slots=True)
class TranscriptExpression:
    """What one expression resource holds for one transcript.

    ``transcript_id`` is the **unversioned** Ensembl identifier, which is what the rest of this
    package uses and what a reader sees in the picker. ``release_id`` is the versioned
    identifier the resource itself used, kept because the two annotations can disagree about a
    transcript's structure even while sharing a stem.
    """

    transcript_id: str
    status: ExpressionStatus
    medians: tuple[TissueMedian, ...] = ()
    release_id: str | None = None

    def __post_init__(self) -> None:
        """Enforce that the status and the values agree.

        Raises:
            ValueError: If a transcript absent from the release carries medians anyway, which
                would mean a number was attached to something nobody quantified.
        """
        if not self.status.is_measured and self.medians:
            raise ValueError(
                f"{self.transcript_id} is recorded as absent from the resource's annotation "
                f"release but carries {len(self.medians)} median(s); a value here would imply "
                f"a measurement that was never made"
            )

    @property
    def highest(self) -> TissueMedian | None:
        """Return the sample group with the largest median, or ``None`` if nothing was measured.

        Ties keep the first in the order the resource returned, which is its own ordering and
        carries no ranking of its own.
        """
        if not self.medians:
            return None
        return max(self.medians, key=lambda entry: entry.median)

    def median_in(self, tissue_id: str) -> TissueMedian | None:
        """Return this transcript's median in one sample group, if it was measured there."""
        return next((entry for entry in self.medians if entry.tissue_id == tissue_id), None)

    def to_dict(self) -> dict[str, object]:
        """Return the record as JSON-serialisable data."""
        return {
            "transcript_id": self.transcript_id,
            "status": str(self.status),
            "status_note": self.status.note,
            "release_id": self.release_id,
            "medians": [entry.to_dict() for entry in self.medians],
        }


@dataclass(frozen=True, slots=True)
class GeneExpression:
    """Every transcript of one gene, as one expression resource holds them.

    Carries the resource's identity and release with the numbers, because "5.35 TPM" means
    nothing without knowing which dataset, which annotation and which sample group produced it.
    """

    gene_id: str
    """Unversioned Ensembl gene identifier, as the rest of this package uses it."""

    release_gene_id: str
    """The versioned identifier the resource used, e.g. ``ENSG00000144285.23``."""

    dataset_id: str
    """The resource's dataset, e.g. ``gtex_v10``."""

    annotation_release: str
    """The annotation the dataset was quantified against, e.g. ``v39`` for GENCODE v39."""

    transcripts: tuple[TranscriptExpression, ...]
    provenance: Provenance
    caveats: tuple[str, ...] = ()

    def for_transcript(self, transcript_id: str) -> TranscriptExpression:
        """Return what is held for one transcript, by unversioned identifier.

        A transcript the resource never saw comes back as
        :attr:`ExpressionStatus.NOT_IN_RELEASE` rather than as ``None``, so a caller cannot
        reach a number without passing the status that qualifies it.
        """
        found = next(
            (entry for entry in self.transcripts if entry.transcript_id == transcript_id), None
        )
        if found is not None:
            return found
        return TranscriptExpression(
            transcript_id=transcript_id, status=ExpressionStatus.NOT_IN_RELEASE
        )

    @property
    def measured_transcripts(self) -> tuple[TranscriptExpression, ...]:
        """Return only the transcripts the resource actually quantified."""
        return tuple(entry for entry in self.transcripts if entry.status.is_measured)

    def to_dict(self) -> dict[str, object]:
        """Return the profile as JSON-serialisable data."""
        return {
            "gene_id": self.gene_id,
            "release_gene_id": self.release_gene_id,
            "dataset_id": self.dataset_id,
            "annotation_release": self.annotation_release,
            "provenance": str(self.provenance),
            "caveats": list(self.caveats),
            "transcripts": [entry.to_dict() for entry in self.transcripts],
        }


__all__ = [
    "ExpressionStatus",
    "GeneExpression",
    "TissueMedian",
    "TranscriptExpression",
]
