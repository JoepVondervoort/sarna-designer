"""What an off-target search found, and under exactly what conditions it found it.

A screening result is only meaningful with its own provenance. NCBI's databases change week to
week, so "no off-target hits" is a statement about one database on one day with one set of
parameters, and re-running it in six months can legitimately give a different answer. Every type
here therefore carries when, where and how it was searched, which is also why a screening is a
separate dated record rather than something written back into the frozen design record
(ADR 0019).
"""

from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any

SearchScope = str
"""Which question a search asked. See :data:`TRANSCRIPT_SCOPE` and :data:`GENOMIC_SCOPE`."""

TRANSCRIPT_SCOPE = "transcript"
"""Could the duplex silence an mRNA? The classic siRNA off-target risk, searched in RefSeq RNA."""

GENOMIC_SCOPE = "genomic"
"""Could the duplex target another promoter? The risk specific to saRNA, searched in the genome."""


@dataclass(frozen=True, slots=True)
class OffTargetHit:
    """One alignment between a candidate and some other sequence.

    Field names follow BLAST's own vocabulary so a reader can check them against the raw report.
    ``identities`` is matching bases, not a percentage: a percentage of a short alignment hides
    how much of the duplex actually matched.
    """

    candidate_id: str
    accession: str
    title: str
    identities: int
    alignment_length: int
    gaps: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    subject_strand: str
    evalue: float
    bit_score: float

    def is_perfect(self, core_length: int) -> bool:
        """Return whether the whole core matched with no gaps and no mismatches."""
        return (
            self.gaps == 0
            and self.identities == core_length
            and self.alignment_length == core_length
        )

    def to_dict(self) -> dict[str, Any]:
        """Export the hit."""
        return {
            "candidate": self.candidate_id,
            "accession": self.accession,
            "title": self.title,
            "identities": self.identities,
            "alignment_length": self.alignment_length,
            "gaps": self.gaps,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "subject_start": self.subject_start,
            "subject_end": self.subject_end,
            "subject_strand": self.subject_strand,
            "evalue": self.evalue,
            "bit_score": self.bit_score,
        }


@dataclass(frozen=True, slots=True)
class PromoterProximity:
    """A genomic hit placed relative to the nearest annotated transcription start site.

    **No threshold decides whether this "is" a promoter hit.** Reporting the distance and the
    gene invents nothing; choosing a cutoff would invent a scientific claim this project has not
    defended. The one flag that is set, :attr:`inside_run_window`, says only that the hit falls
    inside *the same window this run designed against*, a statement about the run's own
    configuration rather than a new definition of a promoter (ADR 0019).
    """

    hit: OffTargetHit
    contig: str
    position: int
    nearest_gene: str | None
    nearest_transcript: str | None
    distance_to_tss: int | None
    """Signed distance in bp; negative is upstream of the TSS. ``None`` when nothing is near."""

    inside_run_window: bool

    def to_dict(self) -> dict[str, Any]:
        """Export the placement."""
        return {
            "hit": self.hit.to_dict(),
            "contig": self.contig,
            "position": self.position,
            "nearest_gene": self.nearest_gene,
            "nearest_transcript": self.nearest_transcript,
            "distance_to_tss": self.distance_to_tss,
            "inside_run_window": self.inside_run_window,
        }


@dataclass(frozen=True, slots=True)
class SearchParameters:
    """Exactly how a search was run, so the same question can be asked again.

    ``blastn-short`` with a small word size is not a refinement here, it is a requirement: at
    19 nt the default ``blastn`` word size and E-value threshold return nothing at all, which
    would look like a clean result.
    """

    scope: SearchScope
    database: str
    organism: str = "txid9606[ORGN]"
    """Restricts the search to human.

    Not a refinement: ``refseq_rna`` and ``refseq_genomic`` span every species RefSeq holds, and
    an unrestricted search of a 19-mer returns hundreds of alignments in plants, fish and frogs.
    They are real alignments and entirely irrelevant, they swamp the measures and would make a
    specific duplex look promiscuous. Measured on a live search before this was added: 500
    alignments returned, 230 of them "near-perfect", none of them human.
    """

    program: str = "blastn"
    task: str = "blastn-short"
    word_size: int = 7
    expect: float = 1000.0
    filter_low_complexity: bool = False
    max_hits: int = 500

    def to_dict(self) -> dict[str, Any]:
        """Export the parameters."""
        return {
            "scope": self.scope,
            "database": self.database,
            "organism": self.organism,
            "program": self.program,
            "task": self.task,
            "word_size": self.word_size,
            "expect": self.expect,
            "filter_low_complexity": self.filter_low_complexity,
            "max_hits": self.max_hits,
        }


@dataclass(frozen=True, slots=True)
class ScreeningRecord:
    """One dated screening of one run: the separate record ADR 0019 requires.

    Not written into the design snapshot, which stays frozen. Re-screening the same run produces
    a second record that can be compared with this one rather than silently replacing it.
    """

    run_id: str
    scope: SearchScope
    parameters: SearchParameters
    hits: tuple[OffTargetHit, ...]
    proximities: tuple[PromoterProximity, ...] = ()
    searched_at: datetime = field(default_factory=lambda: datetime.now(UTC))
    request_id: str | None = None
    """NCBI's RID for this job, so the original report can be retrieved while it lives."""

    candidates_searched: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def hits_for(self, candidate_id: str) -> tuple[OffTargetHit, ...]:
        """Return the hits belonging to one candidate."""
        return tuple(h for h in self.hits if h.candidate_id == candidate_id)

    def to_dict(self) -> dict[str, Any]:
        """Export the whole screening."""
        return {
            "run_id": self.run_id,
            "scope": self.scope,
            "parameters": self.parameters.to_dict(),
            "searched_at": self.searched_at.isoformat(),
            "request_id": self.request_id,
            "candidates_searched": list(self.candidates_searched),
            "hit_count": len(self.hits),
            "hits": [h.to_dict() for h in self.hits],
            "proximities": [p.to_dict() for p in self.proximities],
            "notes": list(self.notes),
        }
