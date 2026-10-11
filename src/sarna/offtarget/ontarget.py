"""Telling a candidate's own locus apart from somewhere else it could bind.

Every candidate is a 19 nt stretch of a real human promoter, so a genomic search finds that
promoter. Left alone, the column a researcher reads first says "1 perfect match elsewhere" for
every candidate in the run and means nothing at all.

**The on-target hit is identified by coordinate, never by sequence or by name.** The run already
records each candidate's genomic interval, and a hit carries the chromosome and the position it
aligned to. Overlap between those two is a fact. Matching on a gene name inside a BLAST
description would be a guess, and this project does not guess silently (ADR 0023).

Two limits, both of which the caller is expected to report rather than paper over:

- A run designed from pasted sequence has no coordinates, so nothing here can identify anything
  and every hit stays where it is.
- The transcript scope gets no removal. A promoter is usually not transcribed, so there is
  normally nothing to remove, and the only available handle would be the gene name.
"""

from collections.abc import Mapping, Sequence
from typing import Any

from sarna.offtarget.models import OffTargetHit
from sarna.offtarget.proximity import contig_for

CandidateLoci = Mapping[str, tuple[str, int, int]]
"""Candidate id to its own ``(contig, start, end)``, 0-based half-open, as the run records it."""


def loci_from_candidates(records: Sequence[Mapping[str, Any]]) -> dict[str, tuple[str, int, int]]:
    """Return each candidate's own genomic interval, skipping any that has none.

    A candidate from a pasted sequence carries ``genomic: null``, and is absent from the result
    rather than present with invented coordinates.
    """
    loci: dict[str, tuple[str, int, int]] = {}
    for record in records:
        genomic = record.get("genomic")
        identifier = record.get("id")
        if not isinstance(genomic, Mapping) or not isinstance(identifier, str):
            continue
        contig, start, end = genomic.get("contig"), genomic.get("start"), genomic.get("end")
        if isinstance(contig, str) and isinstance(start, int) and isinstance(end, int):
            loci[identifier] = (contig, start, end)
    return loci


def is_on_target(hit: OffTargetHit, loci: CandidateLoci) -> bool:
    """Return whether this hit is the candidate's own locus rather than another place.

    The hit's accession is a RefSeq chromosome and its subject coordinates are 1-based
    inclusive, both of which are converted here, because this is the boundary where BLAST's
    conventions meet this project's: intervals are 0-based half-open everywhere inside.
    """
    locus = loci.get(hit.candidate_id)
    if locus is None:
        return False
    contig, start, end = locus
    if contig_for(hit.accession) != contig:
        return False
    low = min(hit.subject_start, hit.subject_end) - 1
    high = max(hit.subject_start, hit.subject_end)
    return low < end and high > start


def split_on_target(
    hits: Sequence[OffTargetHit], loci: CandidateLoci
) -> tuple[tuple[OffTargetHit, ...], tuple[OffTargetHit, ...]]:
    """Return ``(off_target, on_target)``, preserving the order within each.

    Nothing is deleted. The on-target alignments travel on the record so that a count which had
    something taken out of it can say so, and so a reader can check the removal was right.
    """
    off: list[OffTargetHit] = []
    own: list[OffTargetHit] = []
    for hit in hits:
        (own if is_on_target(hit, loci) else off).append(hit)
    return tuple(off), tuple(own)
