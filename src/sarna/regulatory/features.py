"""Ensembl Regulatory Build features, and what their absence does and does not mean.

The Regulatory Build classifies genomic segments as promoters, enhancers, CTCF binding sites, open
chromatin and TF binding sites, from epigenomic data **pooled across a panel of cell types**. That
pooling is the whole character of the resource and governs how it may be read:

- A feature present means "regulatory in at least one of the cell types surveyed". It does not say
  which, and a UI must not imply it says anything about the cell being transfected.
- A feature absent means "not detected in that panel". The panel does not include AF22 or any
  neural progenitor, so a promoter active only in neurons has every reason to be missing.

Both halves are enforced here rather than left to the caller. Every feature carries
:data:`~sarna.regulatory.evidence.EvidenceKind.AGGREGATED_ANNOTATION`, which the evidence model
forbids attaching a cell type to, and a query that finds nothing returns a
:class:`RegulatoryAnnotation` that says so in words — `no_features_note` — rather than an empty list
a caller can quietly render as a clean promoter.

Measured consequence, from `docs/scientific/scn1a-design-coordinates.md`: the Build annotates a
promoter at SCN1A's MANE Select TSS and nothing at either of the two TSSs the thesis designed
against. That is suggestive and it is not evidence that those are not promoters.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any

from sarna.genome.coordinates import GenomicInterval, Strand
from sarna.genome.ensembl import EnsemblClient
from sarna.regulatory.evidence import Evidence, aggregated_annotation
from sarna.rules.provenance import Provenance

ENSEMBL_FEATURE = "regulatory"
"""The Ensembl overlap feature type that returns Regulatory Build segments."""

REGULATORY_BUILD_PROVENANCE = Provenance.citation(
    "Ensembl Regulatory Build",
    "Zerbino et al. 2015, Genome Biol 16:56; segments pooled across cell types, "
    "queried via Ensembl REST /overlap/region",
)


class RegulatoryFeatureKind(Enum):
    """What the Regulatory Build calls a segment."""

    PROMOTER = "promoter"
    ENHANCER = "enhancer"
    CTCF_BINDING_SITE = "CTCF_binding_site"
    TF_BINDING_SITE = "TF_binding_site"
    OPEN_CHROMATIN = "open_chromatin_region"
    PROMOTER_FLANKING = "promoter_flanking_region"
    OTHER = "other"

    @classmethod
    def parse(cls, description: str | None) -> "RegulatoryFeatureKind":
        """Return the matching kind, or :attr:`OTHER` for anything unrecognised.

        Ensembl reports the kind in the record's ``description`` field, not in ``feature_type``,
        which is always the literal string ``regulatory``.
        """
        label = (description or "").strip().lower()
        for member in cls:
            if member.value.lower() == label:
                return member
        return cls.OTHER

    def __str__(self) -> str:
        """Return the kind label."""
        return self.value

    @property
    def is_promoter_associated(self) -> bool:
        """Return whether this kind marks promoter activity rather than something else.

        The distinction that matters for saRNA design: a target inside an annotated promoter has
        independent support for being in a promoter, whereas one inside an enhancer or a CTCF
        site does not, and may be doing something other than promoter-targeted activation.
        """
        return self in (
            RegulatoryFeatureKind.PROMOTER,
            RegulatoryFeatureKind.PROMOTER_FLANKING,
        )


@dataclass(frozen=True, slots=True)
class RegulatoryFeature:
    """One Regulatory Build segment."""

    feature_id: str
    """Ensembl regulatory feature stable id, e.g. ``ENSR2_86TTNX``."""

    kind: RegulatoryFeatureKind
    interval: GenomicInterval

    @property
    def evidence(self) -> Evidence:
        """Return this feature's evidential standing: pooled across cell types, never one."""
        return aggregated_annotation(REGULATORY_BUILD_PROVENANCE)

    def __str__(self) -> str:
        """Return a one-line rendering."""
        return (
            f"{self.kind} {self.interval.contig}:{self.interval.start}-{self.interval.end} "
            f"({self.feature_id})"
        )


@dataclass(frozen=True, slots=True)
class RegulatoryAnnotation:
    """Every Regulatory Build feature overlapping one region, and the caveats on reading it."""

    interval: GenomicInterval
    features: tuple[RegulatoryFeature, ...]

    @property
    def promoters(self) -> tuple[RegulatoryFeature, ...]:
        """Return the promoter and promoter-flanking features."""
        return tuple(f for f in self.features if f.kind.is_promoter_associated)

    def of_kind(self, kind: RegulatoryFeatureKind) -> tuple[RegulatoryFeature, ...]:
        """Return the features of one kind."""
        return tuple(f for f in self.features if f.kind is kind)

    def overlapping(self, interval: GenomicInterval) -> tuple[RegulatoryFeature, ...]:
        """Return the features overlapping ``interval``, for annotating one candidate."""
        return tuple(f for f in self.features if f.interval.overlaps(interval))

    @property
    def is_empty(self) -> bool:
        """Return whether nothing was annotated in this region."""
        return not self.features

    @property
    def no_features_note(self) -> str | None:
        """Return the caveat to display when nothing was found, or ``None`` if something was.

        A separate property rather than a comment, so an interface that renders the annotation
        has to deal with the empty case explicitly instead of drawing a clean promoter.
        """
        if not self.is_empty:
            return None
        return (
            "The Ensembl Regulatory Build annotates nothing in this region. It is built from "
            "epigenomic data pooled over a panel of cell types that includes neither AF22 nor any "
            "neural progenitor, so absence here is weak evidence of absence: a promoter active "
            "only in neurons would be missing from it."
        )

    @property
    def summary(self) -> str:
        """Return a plain-language statement of what is annotated here."""
        if self.is_empty:
            return "No regulatory features annotated in this region."
        counts: dict[str, int] = {}
        for feature in self.features:
            counts[str(feature.kind)] = counts.get(str(feature.kind), 0) + 1
        listing = ", ".join(f"{count} {kind}" for kind, count in sorted(counts.items()))
        return f"Annotated across pooled cell types: {listing}."

    def to_dict(self) -> dict[str, Any]:
        """Return the annotation as JSON-serialisable data, caveat included."""
        return {
            "region": {
                "contig": self.interval.contig,
                "start": self.interval.start,
                "end": self.interval.end,
                "assembly": str(self.interval.assembly),
            },
            "evidence_kind": str(aggregated_annotation(REGULATORY_BUILD_PROVENANCE).kind),
            "provenance": str(REGULATORY_BUILD_PROVENANCE),
            "summary": self.summary,
            "caveat": self.no_features_note,
            "features": [
                {
                    "id": f.feature_id,
                    "kind": str(f.kind),
                    "start": f.interval.start,
                    "end": f.interval.end,
                    "promoter_associated": f.kind.is_promoter_associated,
                }
                for f in self.features
            ],
        }


def fetch_regulatory_features(
    client: EnsemblClient, interval: GenomicInterval
) -> RegulatoryAnnotation:
    """Fetch Regulatory Build features overlapping ``interval``.

    Args:
        client: An Ensembl client, pinned to the assembly ``interval`` is on.
        interval: The region to annotate.

    Returns:
        The features, in coordinate order, with the caveats on reading them.

    Raises:
        ValueError: If the interval is on another assembly than the client serves.
        EnsemblUnavailableError: If Ensembl cannot be reached.
    """
    records = client.overlap(interval, ENSEMBL_FEATURE)
    return RegulatoryAnnotation(
        interval=interval,
        features=tuple(
            sorted(
                (_feature_from(record, interval) for record in records),
                key=lambda f: (f.interval.start, f.interval.end, f.feature_id),
            )
        ),
    )


def _feature_from(record: dict[str, Any], region: GenomicInterval) -> RegulatoryFeature:
    """Build a feature from an Ensembl record, converting to 0-based half-open."""
    return RegulatoryFeature(
        feature_id=str(record.get("id", "unknown")),
        # Ensembl reports the kind in `description`; `feature_type` is always "regulatory".
        kind=RegulatoryFeatureKind.parse(record.get("description")),
        interval=GenomicInterval(
            assembly=region.assembly,
            contig=str(record.get("seq_region_name", region.contig)),
            start=int(record["start"]) - 1,
            end=int(record["end"]),
            # Regulatory features are not strand-specific; Ensembl reports strand 0 for most.
            strand=Strand.PLUS,
        ),
    )


def promoter_support(
    annotation: RegulatoryAnnotation, interval: GenomicInterval
) -> tuple[bool, str]:
    """Report whether an annotated promoter overlaps ``interval``, and say what that means.

    Returns:
        Whether a promoter-associated feature overlaps, and a sentence stating the strength of
        the claim either way. The sentence exists because both answers are easy to over-read.
    """
    overlapping = tuple(
        f for f in annotation.overlapping(interval) if f.kind.is_promoter_associated
    )
    if overlapping:
        ids = ", ".join(f.feature_id for f in overlapping)
        return True, (
            f"Overlaps an annotated promoter ({ids}). Independent support for this being a "
            f"promoter, from pooled cell types — not evidence that it is active in yours."
        )
    if annotation.is_empty:
        return False, (
            "Nothing is annotated in the surrounding region at all, so this says nothing either "
            "way. See the region's caveat."
        )
    return False, (
        "Other regulatory features are annotated nearby but no promoter overlaps this target. "
        "Worth checking what it does overlap: an enhancer or CTCF site is a different "
        "mechanistic setting from a promoter."
    )


def features_by_kind(
    features: Sequence[RegulatoryFeature],
) -> dict[RegulatoryFeatureKind, tuple[RegulatoryFeature, ...]]:
    """Group features by kind, for rendering one track per kind."""
    grouped: dict[RegulatoryFeatureKind, list[RegulatoryFeature]] = {}
    for feature in features:
        grouped.setdefault(feature.kind, []).append(feature)
    return {kind: tuple(members) for kind, members in grouped.items()}


__all__ = [
    "ENSEMBL_FEATURE",
    "REGULATORY_BUILD_PROVENANCE",
    "RegulatoryAnnotation",
    "RegulatoryFeature",
    "RegulatoryFeatureKind",
    "features_by_kind",
    "fetch_regulatory_features",
    "promoter_support",
]
