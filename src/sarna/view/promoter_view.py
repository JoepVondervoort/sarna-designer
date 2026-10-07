"""Everything an interface needs to draw one promoter, assembled once and serialisable.

The engine already decides everything worth showing; this module's job is to put it into one object
in **promoter coordinates**, so that a renderer never does biology. Nothing here recomputes a
verdict, and nothing here does strand arithmetic — genomic features are placed through
:meth:`Promoter.span_of`, which is the inverse of the mapping the design pipeline used.

Three things the view insists on, because they are what a promoter diagram usually gets wrong.

**Every TSS, not just the chosen one.** A gene has many; SCN1A has sixteen. A diagram showing only
the transcript the run selected hides the fact that a different choice would have moved the whole
window. The others are drawn as what they are — alternatives the run did not use (ADR 0012).

**Rejected candidates stay on the drawing.** They carry the rules that excluded them. A promoter
view showing only survivors is the behaviour this project replaces.

**Annotations carry their evidence kind.** A motif match and an observed ChIP-seq peak must not be
drawn the same way, and a region with nothing annotated carries the caveat explaining why that is
weak evidence. The view passes both through from :mod:`sarna.regulatory` rather than flattening
them into "features".
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from sarna.design.candidate import Candidate
from sarna.design.pipeline import DesignResult
from sarna.genome.ensembl import Gene, Transcript
from sarna.genome.promoter import Promoter
from sarna.regulatory.chromatin import Af22ChromatinReader, ChromatinProfile
from sarna.regulatory.evidence import (
    CELL_TYPES,
    CellType,
    CellTypeStatus,
    DataAxis,
)
from sarna.regulatory.features import RegulatoryAnnotation
from sarna.regulatory.motifs import MotifScan


def _compact(candidate: Candidate) -> dict[str, Any]:
    """Return a rejected candidate's identity, position and reasons, without every measurement.

    Keeps the project's rule that a rejection is never silently dropped while keeping the document
    loadable: the candidate is present, positioned, and names what excluded it.
    """
    return {
        "id": candidate.id,
        "rank": None,
        "offset": candidate.offset,
        "target_site_dna": candidate.target.bases,
        "gc_percent": round(candidate.target.gc_percent, 1),
        "tss_span": (
            [candidate.tss_span[0].value, candidate.tss_span[1].value]
            if candidate.tss_span
            else None
        ),
        "preference_points": candidate.preference_points,
        "max_preference_points": candidate.max_preference_points,
        "excluded_by": [result.rule_id for result in candidate.excluded_by],
        "rejection_summary": candidate.rejection_summary,
        "abbreviated": True,
    }


@dataclass(frozen=True, slots=True)
class PlacedSpan:
    """A feature placed in promoter coordinates, ready to draw."""

    start: int
    """0-based offset into the promoter sequence."""

    end: int
    """0-based exclusive end."""

    label: str
    kind: str
    detail: str
    """What a reader should be told when they hover or click."""

    clipped: bool = False
    """Whether the feature extends past the promoter edge and is drawn truncated."""

    def to_dict(self) -> dict[str, Any]:
        """Return the span as JSON-serialisable data."""
        return {
            "start": self.start,
            "end": self.end,
            "label": self.label,
            "kind": self.kind,
            "detail": self.detail,
            "clipped": self.clipped,
        }


@dataclass(frozen=True, slots=True)
class CellTypePanel:
    """What is known about one cell type at this promoter, including when that is nothing.

    Exists so an interface cannot render an empty column as a negative result. A cell type with no
    data is ``NOT_ASSAYED`` with an ``absence_note`` saying so in words; one with chromatin tracks
    carries a real :class:`~sarna.regulatory.chromatin.ChromatinProfile` and its own verdict.

    Note what each half covers. The chromatin profile says whether the promoter is **accessible**
    in this cell type. It says nothing about which transcription factor is bound there, and the
    motif predictions elsewhere in the view say nothing about this cell type. Both gaps are stated
    rather than papered over.
    """

    cell_type: CellType
    status: CellTypeStatus
    observations: tuple[PlacedSpan, ...] = ()
    absence_note: str | None = None
    chromatin: ChromatinProfile | None = None
    signal_tracks: Mapping[str, tuple[float, ...]] = field(default_factory=dict)
    """Per-mark binned signal across the promoter, in transcript 5'-to-3' order."""

    @property
    def has_chromatin(self) -> bool:
        """Return whether real chromatin data was read for this cell type."""
        return self.chromatin is not None and self.chromatin.is_readable

    def to_dict(self) -> dict[str, Any]:
        """Return the panel as JSON-serialisable data."""
        return {
            "id": self.cell_type.id,
            "name": self.cell_type.name,
            "description": self.cell_type.description,
            "proliferating": self.cell_type.is_proliferating,
            "status": str(self.status),
            "known": self.status.is_known,
            "observations": [span.to_dict() for span in self.observations],
            "absence_note": self.absence_note,
            "data_axes": [
                self.cell_type.availability(axis).to_dict() for axis in DataAxis
            ],
            "chromatin": self.chromatin.to_dict() if self.chromatin else None,
            "signal_tracks": {
                mark: [round(value, 3) for value in values]
                for mark, values in self.signal_tracks.items()
            },
        }


@dataclass(frozen=True, slots=True)
class PromoterView:
    """One promoter, its candidates and its annotations, all in promoter coordinates."""

    result: DesignResult
    gene: Gene | None = None
    transcript: Transcript | None = None
    alternative_tss: tuple[PlacedSpan, ...] = ()
    regulatory: RegulatoryAnnotation | None = None
    motif_scans: tuple[MotifScan, ...] = ()
    cell_panels: tuple[CellTypePanel, ...] = field(default_factory=tuple)

    @property
    def promoter(self) -> Promoter:
        """Return the promoter being drawn."""
        return self.result.promoter

    @property
    def tss_offset(self) -> int | None:
        """Return the TSS's offset into the promoter, or ``None`` if it falls outside."""
        coordinate = self.promoter.tss_coordinate
        if coordinate is None:
            return None
        return self.promoter.offset_of(coordinate)

    def island_spans(self) -> tuple[PlacedSpan, ...]:
        """Return the CpG islands as placed spans.

        Islands are detected on the promoter sequence itself, so they are already in promoter
        coordinates and need no mapping.
        """
        if self.result.cpg_islands is None:
            return ()
        return tuple(
            PlacedSpan(
                start=island.start,
                end=island.end,
                label="CpG island",
                kind="cpg-island",
                detail=(
                    f"{island.end - island.start} bp, GC {island.gc_fraction:.0%}, "
                    f"observed/expected CpG {island.observed_expected_ratio:.2f}. Targets "
                    f"overlapping an island are excluded by default."
                ),
            )
            for island in self.result.cpg_islands
        )

    def regulatory_spans(self) -> tuple[PlacedSpan, ...]:
        """Return the Regulatory Build features placed on the promoter, clipped to it."""
        if self.regulatory is None:
            return ()
        placed = []
        for feature in self.regulatory.features:
            span = self.promoter.span_of(feature.interval)
            if span is None:
                continue
            start, end = span
            placed.append(
                PlacedSpan(
                    start=start,
                    end=end,
                    label=str(feature.kind),
                    kind=f"regulatory-{feature.kind.value}",
                    detail=(
                        f"{feature.feature_id}. {feature.evidence.kind.confidence_note}"
                    ),
                    clipped=(end - start) != feature.interval.length,
                )
            )
        return tuple(placed)

    def motif_spans(self, minimum_enrichment: float = 1.0) -> tuple[PlacedSpan, ...]:
        """Return motif matches as placed spans, from scans that beat chance.

        Args:
            minimum_enrichment: Scans whose enrichment over a composition-matched shuffle is
                below this are dropped entirely. The default of 1.0 keeps only motifs found more
                often than base composition alone explains, which in a GC-rich promoter removes
                most of them.

        Returns:
            The matches, in offset order. Motif scans are performed on the promoter sequence, so
            their offsets need no mapping.
        """
        placed = []
        for scan in self.motif_scans:
            if (scan.enrichment or 0.0) < minimum_enrichment:
                continue
            for site in scan.sites:
                placed.append(
                    PlacedSpan(
                        start=site.offset,
                        end=site.end,
                        label=scan.motif.name,
                        kind="motif",
                        detail=(
                            f"{scan.motif.name} ({scan.motif.matrix_id}), relative score "
                            f"{site.relative_score:.2f} on the {site.strand} strand. "
                            f"{scan.summary} {scan.motif.data_type} matrix. "
                            f"{site.evidence.kind.confidence_note}"
                        ),
                    )
                )
        return tuple(sorted(placed, key=lambda s: (s.start, s.label)))

    def ranked_candidates(self) -> tuple[Candidate, ...]:
        """Return every generated candidate, in position order, carrying its rank.

        ``DesignResult.ranked`` holds *copies* that carry a rank, because ranking must not
        retroactively alter candidates a caller already holds. ``DesignResult.generated`` therefore
        has ``rank is None`` throughout, and exporting it directly produces a list a renderer
        cannot order — which it does silently, by falling back to position. This merges the two.
        """
        ranks = {candidate.id: candidate for candidate in self.result.ranked}
        return tuple(
            ranks.get(candidate.id, candidate) for candidate in self.result.generated
        )

    def rejection_breakdown(self) -> tuple[dict[str, Any], ...]:
        """Return how many candidates each mandatory rule excluded, most first.

        What a researcher actually wants from 2 000 rejections: not the list, but which rule is
        doing the work. A rule rejecting almost everything is either the real constraint on this
        promoter or a setting worth revisiting, and either way it should be the first thing read.
        """
        counts: dict[str, int] = {}
        for candidate in self.result.rejected:
            for failing in candidate.excluded_by:
                counts[failing.rule_id] = counts.get(failing.rule_id, 0) + 1
        total = len(self.result.generated)
        return tuple(
            {
                "rule": rule_id,
                "rejected": count,
                "percent_of_generated": round(100.0 * count / total, 1) if total else 0.0,
            }
            for rule_id, count in sorted(counts.items(), key=lambda item: -item[1])
        )

    def to_dict(
        self,
        motif_minimum_enrichment: float = 1.0,
        full_rules_for_rejected: bool = True,
    ) -> dict[str, Any]:
        """Return the whole view as JSON-serialisable data for a renderer.

        Args:
            motif_minimum_enrichment: Motif scans below this enrichment are left out of the
                drawable track, though they stay in ``motif_scans`` so a reader can see what was
                looked for and not found.
            full_rules_for_rejected: When false, rejected candidates carry only their identity,
                position and the rules that excluded them, rather than every rule's measurements.
                Every candidate is still present — nothing is dropped — but a 2 500 bp promoter
                goes from about 21 MB to under 2 MB, which is the difference between a document a
                browser can load and one it cannot. Use ``rejection_breakdown`` for the summary.
        """
        promoter = self.promoter
        interval = promoter.interval
        return {
            "gene": (
                {
                    "symbol": self.gene.symbol,
                    "ensembl_id": self.gene.ensembl_id,
                    "description": self.gene.description,
                }
                if self.gene
                else None
            ),
            "transcript": (
                {
                    "ensembl_id": self.transcript.ensembl_id,
                    "biotype": self.transcript.biotype,
                    "refseq_match": self.transcript.refseq_match,
                    "is_mane_select": self.transcript.is_mane_select,
                    "tss_coordinate": self.transcript.tss_coordinate,
                }
                if self.transcript
                else None
            ),
            "promoter": {
                "length": promoter.length,
                "sequence": promoter.sequence.bases,
                "gc_percent": round(promoter.sequence.gc_percent, 1),
                "tss_start": promoter.tss_start.value if promoter.tss_start else None,
                "tss_offset": self.tss_offset,
                "assembly": str(interval.assembly) if interval else None,
                "contig": interval.contig if interval else None,
                "strand": str(interval.strand) if interval else None,
                "start": interval.start if interval else None,
                "end": interval.end if interval else None,
            },
            "configuration": self.result.config.to_dict(),
            "counts": {
                "generated": len(self.result.generated),
                "passing": len(self.result.passing),
                "rejected": len(self.result.rejected),
            },
            "candidates": [
                candidate.to_dict()
                if full_rules_for_rejected or not candidate.is_rejected
                else _compact(candidate)
                for candidate in self.ranked_candidates()
            ],
            "rejection_breakdown": list(self.rejection_breakdown()),
            "tracks": {
                "cpg_islands": [s.to_dict() for s in self.island_spans()],
                "alternative_tss": [s.to_dict() for s in self.alternative_tss],
                "regulatory": [s.to_dict() for s in self.regulatory_spans()],
                "motifs": [
                    s.to_dict() for s in self.motif_spans(motif_minimum_enrichment)
                ],
            },
            "regulatory_summary": (
                self.regulatory.to_dict() if self.regulatory else None
            ),
            "motif_scans": [
                {
                    "motif": scan.motif.name,
                    "matrix_id": scan.motif.matrix_id,
                    "relevance": scan.motif.relevance,
                    "data_type": str(scan.motif.data_type),
                    "in_vivo": scan.motif.data_type.is_in_vivo,
                    "sites": len(scan.sites),
                    "expected": round(scan.expected_hits, 1),
                    "enrichment": (
                        round(scan.enrichment, 2) if scan.enrichment is not None else None
                    ),
                    "exceeds_chance": scan.exceeds_chance,
                    "summary": scan.summary,
                }
                for scan in self.motif_scans
            ],
            "cell_types": [panel.to_dict() for panel in self.cell_panels],
            "unevaluated_rules": list(self.result.unevaluated_rule_ids),
        }


def alternative_tss_spans(
    promoter: Promoter,
    transcripts: Sequence[Transcript],
    selected: Transcript | None,
) -> tuple[PlacedSpan, ...]:
    """Place every other transcript's TSS on the promoter, as one-base spans.

    The transcripts whose TSS falls outside the promoter window are omitted — not hidden, since
    the view's caller reports the total — because there is nowhere on the drawing to put them.

    Args:
        promoter: The promoter being drawn.
        transcripts: Every transcript of the gene.
        selected: The transcript the run used, excluded from the result.

    Returns:
        One span per distinct alternative TSS inside the window, in offset order.
    """
    by_offset: dict[int, list[Transcript]] = {}
    for transcript in transcripts:
        if selected is not None and transcript.ensembl_id == selected.ensembl_id:
            continue
        offset = promoter.offset_of(transcript.tss_coordinate)
        if offset is None:
            continue
        by_offset.setdefault(offset, []).append(transcript)
    return tuple(
        PlacedSpan(
            start=offset,
            end=offset + 1,
            label="alternative TSS",
            kind="alternative-tss",
            detail=(
                f"Transcription start site of "
                f"{', '.join(t.ensembl_id for t in sorted(members, key=lambda x: x.ensembl_id))}"
                f". This run designed against a different transcript, so positions on this "
                f"drawing are not relative to this TSS."
            ),
        )
        for offset, members in sorted(by_offset.items())
    )


def cell_panels_for(
    cell_types: Sequence[CellType] | None = None,
    promoter: Promoter | None = None,
    chromatin_reader: Af22ChromatinReader | None = None,
    bins: int = 200,
) -> tuple[CellTypePanel, ...]:
    """Build cell-type panels, filling in chromatin state where tracks are available.

    Args:
        cell_types: Which cell types to build panels for. Defaults to those this project
            characterises.
        promoter: The promoter being viewed, needed to read signal across it.
        chromatin_reader: An open track reader, or ``None``. It states which cell line its
            tracks belong to, and only that line's panel gets chromatin.
        bins: How many slices to read the signal in, for drawing.

    Returns:
        One panel per cell type. A panel with no data is ``NOT_ASSAYED`` and says why, which is
        what stops an interface from rendering an empty panel as "nothing binds here".
    """
    selected = list(CELL_TYPES.values()) if cell_types is None else list(cell_types)
    panels = []
    for cell in selected:
        profile = None
        tracks: dict[str, tuple[float, ...]] = {}
        if (
            chromatin_reader is not None
            and promoter is not None
            and promoter.interval is not None
            # The reader states which cell line its tracks belong to, so ask it rather than
            # testing identity against AF22. One reader still means one cell line can carry
            # chromatin; a mapping of cell type to reader is what lifts that.
            and cell is chromatin_reader.cell_type
        ):
            profile = chromatin_reader.profile(promoter.interval)
            for mark in chromatin_reader.marks:
                values = chromatin_reader.signal_track(promoter.interval, mark, bins)
                if values is not None:
                    tracks[str(mark)] = values
        panels.append(
            CellTypePanel(
                cell_type=cell,
                status=(
                    CellTypeStatus.OBSERVED_BOUND
                    if profile is not None and profile.activating_marks
                    else CellTypeStatus.ASSAYED_NOT_BOUND
                    if profile is not None and profile.is_readable
                    else CellTypeStatus.NOT_ASSAYED
                ),
                absence_note=_absence_note(cell, profile),
                chromatin=profile,
                signal_tracks=tracks,
            )
        )
    return tuple(panels)


def _absence_note(cell: CellType, profile: ChromatinProfile | None) -> str:
    """Return what this panel cannot tell the reader, which is never nothing.

    Even with chromatin tracks loaded, accessibility is not occupancy: nothing here says which
    factor is bound. Saying so is the point of the note.
    """
    if profile is not None and profile.is_readable:
        return (
            f"Chromatin accessibility in {cell.name} is measured here; factor occupancy is not. "
            f"No ChIP-seq for individual transcription factors exists for this line, so the motif "
            f"matches shown elsewhere remain sequence predictions and say nothing about "
            f"{cell.name}."
        )
    base = (
        f"No factor-occupancy data is loaded for {cell.name}. Motif matches shown elsewhere are "
        f"sequence predictions and say nothing about this cell type. "
    )
    if cell.is_proliferating:
        return base + (
            "ENCODE holds extensive ChIP-seq for HeLa-S3 and wiring it in would populate this "
            "panel properly."
        )
    return base + (
        "No public ChIP-seq exists for this line. Its histone tracks give chromatin state rather "
        "than factor occupancy, and are not loaded in this run."
    )


__all__ = [
    "CellTypePanel",
    "PlacedSpan",
    "PromoterView",
    "alternative_tss_spans",
    "cell_panels_for",
]
