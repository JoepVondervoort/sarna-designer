"""Placing a genomic off-target hit relative to the nearest annotated start site.

An saRNA duplex carries a risk ordinary siRNA screening does not look for: it may target
*another promoter* and change the expression of a gene nobody intended. A transcript search
cannot see that, so genomic hits are placed against annotated TSSs here.

**No threshold decides whether a hit "is" in a promoter.** This project has not defended a
definition of promoter proximity, and inventing one (±2 kb, say) would turn a convenience into
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


BLOCK_SIZE = 1_000_000
"""How much genome one cached transcript query covers.

Chosen by measurement over simulated screenings, counting requests against placing each hit
directly. Smaller blocks lose: at 100 kb a hit's own 20 kb interval straddles a boundary often
enough that scattered hits cost **more** requests than no cache at all.

    block     scattered hits      clustered hits
    100 kb        -6%                 +61%
    250 kb       +18%                 +67%
    500 kb       +37%                 +68%
      1 Mb       +57%                 +70%
      2 Mb       +76%                 +74%

1 Mb is where the scattered case is comfortably positive while one response stays a modest
number of transcripts; Ensembl's overlap endpoint accepts up to 5 Mb, so there is headroom if
a later measurement wants it. Real screenings sit nearer the clustered column, because the
candidates of one run share a promoter and draw repeated alignments to the same loci.

Snapping to a grid is also what makes the client's on-disk response cache work across runs:
per-hit intervals are all different and never hit it twice, while block boundaries repeat.
"""


def _overlaps(transcript: dict[str, Any], interval: GenomicInterval) -> bool:
    """Return whether an Ensembl transcript record overlaps a half-open interval.

    Ensembl reports 1-based inclusive coordinates and :class:`GenomicInterval` is 0-based
    half-open, which is the conversion the whole adapter layer exists to keep in one place.
    """
    start, end = transcript.get("start"), transcript.get("end")
    if not isinstance(start, int) or not isinstance(end, int):
        return False
    return start - 1 < interval.end and end > interval.start


class BlockedTranscripts:
    """Answer many small overlap queries out of a few large ones.

    Placing genomic hits was an N+1: :func:`place_hit` asks for the transcripts around one hit,
    and a screening returning the full five hundred alignments therefore made five hundred
    sequential Ensembl requests. Ensembl rate-limits, so that is both slow and close to rude.

    Queries are snapped to a fixed grid and each block fetched once. **The result is filtered
    back down to exactly the interval that was asked for**, so this is a caching layer and not
    a change to which transcripts are considered; a test asserts the two agree hit for hit.
    """

    __slots__ = ("_block", "_blocks", "_lookup", "requests")

    def __init__(self, lookup: TranscriptLookup, block: int = BLOCK_SIZE) -> None:
        """Wrap ``lookup``, caching whole blocks of ``block`` bases."""
        self._lookup = lookup
        self._block = block
        self._blocks: dict[tuple[str, int], list[dict[str, Any]]] = {}
        self.requests = 0
        """Underlying lookups performed, so the saving can be measured rather than claimed."""

    def __call__(self, interval: GenomicInterval) -> list[dict[str, Any]]:
        """Return the transcripts overlapping ``interval``, from cached blocks."""
        first = max(0, interval.start) // self._block
        last = max(first, (max(0, interval.end) - 1) // self._block)
        seen: dict[str, dict[str, Any]] = {}
        for index in range(first, last + 1):
            for transcript in self._block_at(interval, index):
                # A transcript spanning a block boundary comes back from both blocks.
                key = str(transcript.get("id") or (transcript.get("start"),
                                                   transcript.get("end")))
                seen[key] = transcript
        return [t for t in seen.values() if _overlaps(t, interval)]

    def _block_at(self, interval: GenomicInterval, index: int) -> list[dict[str, Any]]:
        """Return one cached block of transcripts, fetching it if this is the first ask.

        The block is clipped to the contig. A fixed grid's last block always ends past the end
        of the chromosome, and ``GenomicInterval`` refuses that, so a hit near a telomere --
        ordinary in BLAST output -- would otherwise have failed the whole screening.
        """
        key = (interval.contig, index)
        cached = self._blocks.get(key)
        if cached is None:
            limit = interval.assembly.contig_length(interval.contig)
            start = min(index * self._block, max(0, limit - 1))
            end = min((index + 1) * self._block, limit)
            self.requests += 1
            cached = list(self._lookup(GenomicInterval(
                interval.assembly, interval.contig, start, end, Strand.PLUS,
            )))
            self._blocks[key] = cached
        return cached


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
    # Both ends are clipped to the contig. The low end always was; the high end was not, so a
    # hit within search_radius of a telomere built an interval past the end of the chromosome
    # and GenomicInterval refused it -- failing the entire screening over one alignment.
    limit = Assembly.GRCH38.contig_length(contig)
    interval = GenomicInterval(
        Assembly.GRCH38, contig,
        max(0, min(position - search_radius, limit - 1)),
        min(position + search_radius, limit),
        Strand.PLUS,
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
