"""Placing a genomic off-target hit relative to the nearest annotated start site.

An saRNA duplex carries a risk ordinary siRNA screening does not look for: it may target
*another promoter* and change the expression of a gene nobody intended. A transcript search
cannot see that, so genomic hits are placed against annotated TSSs here.

**No threshold decides whether a hit "is" in a promoter.** This project has not defended a
definition of promoter proximity, and inventing one — ±2 kb, say — would turn a convenience into
a scientific claim. What is reported instead is the distance to the nearest annotated TSS and
which gene it belongs to, which invents nothing. One flag is set,
:attr:`~sarna.offtarget.models.PromoterProximity.inside_run_window`, and it says only that the
hit falls inside *the same window this run designed against*: a statement about the run's own
configuration, not a new definition (ADR 0019).

The Ensembl Regulatory Build was considered for this and rejected. `docs/scientific/tf-binding.md`
and `docs/scientific/af22-chromatin.md` already record that it pools observations across cell
types and misses promoters active only in neurons, so its silence would be weak evidence
presented as reassurance.
"""

from collections.abc import Callable, Sequence
from typing import Any

from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.offtarget.models import OffTargetHit, PromoterProximity

REFSEQ_CHROMOSOMES = {
    f"NC_{n:06d}": name
    for n, name in enumerate(
        [*(str(i) for i in range(1, 23)), "X", "Y"], start=1
    )
}
"""RefSeq chromosome accessions to Ensembl contig names, for GRCh38.

``NC_000001`` is chromosome 1 and so on to ``NC_000023`` (X) and ``NC_000024`` (Y). A hit on a
scaffold or an unplaced contig is not in this map and is reported without a placement rather
than guessed at.
"""

TranscriptLookup = Callable[[GenomicInterval], Sequence[dict[str, Any]]]
"""Returns annotated transcripts overlapping an interval, as Ensembl's overlap endpoint does."""


def contig_for(accession: str) -> str | None:
    """Return the Ensembl contig name for a RefSeq accession, or ``None`` if unmapped."""
    return REFSEQ_CHROMOSOMES.get(accession.split(".")[0].upper())


def tss_of(transcript: dict[str, Any]) -> int | None:
    """Return a transcript's start site, which is its high coordinate on the minus strand."""
    start, end, strand = transcript.get("start"), transcript.get("end"), transcript.get("strand")
    if not isinstance(start, int) or not isinstance(end, int):
        return None
    return end if strand == -1 else start


def signed_distance(position: int, tss: int, strand: int) -> int:
    """Return the distance from ``tss`` to ``position``, negative upstream.

    Upstream and downstream are defined by the transcript's own direction, which is why the sign
    flips on the minus strand. Reporting an unsigned distance would lose the distinction between
    a hit in a promoter and one inside the gene body.
    """
    offset = position - tss
    return -offset if strand == -1 else offset


def place_hit(
    hit: OffTargetHit,
    lookup: TranscriptLookup,
    window: tuple[int, int] | None,
    search_radius: int = 10_000,
) -> PromoterProximity:
    """Place one genomic hit against the nearest annotated start site.

    Args:
        hit: A hit whose accession is a RefSeq chromosome.
        lookup: Returns transcripts overlapping an interval.
        window: The run's own design window as ``(start, end)`` in TSS coordinates, or ``None``
            when the run used no window. Only used to set ``inside_run_window``.
        search_radius: How far either side to look for a transcript. A hit with nothing within
            this distance is reported with no nearest gene rather than with a distant one.

    Returns:
        The placement, with ``None`` distances where nothing could be established.
    """
    contig = contig_for(hit.accession)
    position = min(hit.subject_start, hit.subject_end)
    if contig is None:
        return PromoterProximity(
            hit=hit, contig=hit.accession, position=position, nearest_gene=None,
            nearest_transcript=None, distance_to_tss=None, inside_run_window=False,
        )
    interval = GenomicInterval(
        Assembly.GRCH38, contig,
        max(0, position - search_radius), position + search_radius, Strand.PLUS,
    )
    nearest: dict[str, Any] | None = None
    nearest_distance: int | None = None
    for transcript in lookup(interval):
        tss = tss_of(transcript)
        if tss is None:
            continue
        distance = signed_distance(position, tss, int(transcript.get("strand", 1)))
        if nearest_distance is None or abs(distance) < abs(nearest_distance):
            nearest, nearest_distance = transcript, distance
    inside = (
        nearest_distance is not None
        and window is not None
        and window[0] <= nearest_distance <= window[1]
    )
    return PromoterProximity(
        hit=hit,
        contig=contig,
        position=position,
        nearest_gene=(nearest or {}).get("external_name") or (nearest or {}).get("Parent"),
        nearest_transcript=(nearest or {}).get("id"),
        distance_to_tss=nearest_distance,
        inside_run_window=inside,
    )


def describe(placement: PromoterProximity) -> str:
    """Return one sentence stating where the hit is, and claiming nothing more."""
    if placement.distance_to_tss is None:
        return (
            f"On {placement.contig} at {placement.position:,}, with no annotated transcript "
            f"nearby. Distance to a start site could not be established."
        )
    direction = "upstream of" if placement.distance_to_tss < 0 else "downstream of"
    gene = placement.nearest_gene or placement.nearest_transcript or "an unnamed transcript"
    sentence = (
        f"On {placement.contig} at {placement.position:,}, "
        f"{abs(placement.distance_to_tss):,} bp {direction} the start site of {gene}."
    )
    if placement.inside_run_window:
        sentence += (
            " That falls inside the same window this run designed against, so this sequence "
            "would have been a legal target for that gene too. Proximity to an annotated start "
            "site does not establish that the promoter is active or that the duplex binds it."
        )
    return sentence
