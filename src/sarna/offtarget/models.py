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
    mismatch_positions: tuple[int, ...] = ()
    """Which bases of the *query* did not match, 1-based, within the aligned span.

    Taken from BLAST's own alignment strings rather than inferred from the identity count,
    because "17 of 19" does not say *which* two, and that is the whole question when a reader
    is deciding whether a partial match matters. Empty when the report carried no alignment
    strings, which is not the same as a perfect match: see :attr:`has_alignment_detail`.
    """

    has_alignment_detail: bool = False
    """Whether the report carried the aligned sequences these positions were read from.

    Without this an empty :attr:`mismatch_positions` would be drawn as "every base matched",
    turning a missing field into a finding.
    """

    query_sequence: str = ""
    subject_sequence: str = ""
    """Aligned strings from BLAST, including gaps; DNA alphabet, oriented as reported."""

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
            "mismatch_positions": list(self.mismatch_positions),
            "has_alignment_detail": self.has_alignment_detail,
            "query_sequence": self.query_sequence,
            "subject_sequence": self.subject_sequence,
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
    on_chromosome: bool = True
    """Whether the hit's accession is an assembled chromosome this project can place.

    Two different ``None`` distances otherwise read alike: a hit on a scaffold that cannot be
    placed at all, and a hit placed on a chromosome with no annotated transcript near it. The
    first is a limit of the lookup, the second is a finding about the locus.
    """

    annotation_status: str = "available"
    """available, unavailable, or unmapped; an unavailable lookup is never a negative finding."""

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
            "on_chromosome": self.on_chromosome,
            "annotation_status": self.annotation_status,
        }


BY_DATABASE = "database"
"""Human because the database holds one species. The stronger of the two guarantees."""

BY_ENTREZ_FILTER = "entrez-filter"
"""Human because an Entrez query restricted a multi-species database.

Weaker, and worth naming as such: ``ENTREZ_QUERY`` is not a documented URL API parameter, and
NCBI's own page says parameters it does not list "should not be relied upon". A scope resting
on it is human by request rather than by construction.
"""


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

    species_restriction: str = BY_ENTREZ_FILTER
    """How this search was confined to human. See :data:`BY_DATABASE`, :data:`BY_ENTREZ_FILTER`.

    Recorded because the two are not equally strong and the interface must not word them
    alike: one is a property of the database searched, the other a request that an
    undocumented parameter be honoured.
    """

    program: str = "blastn"
    task: str = "explicit-short-query"
    reward: int | None = 1
    penalty: int | None = -3
    gap_costs: str | None = "5 2"
    word_size: int = 7
    expect: float = 1000.0
    filter_low_complexity: bool = False
    max_hits: int = 100
    """``HITLIST_SIZE``: database sequences kept **per query**, which is NCBI's own unit.

    Was 500, which at fifteen queries is up to 7,500 alignments of which the overwhelming
    majority match 12 or 13 of 19 bases and mean nothing. BLAST keeps the best by score, so
    lowering this drops the weakest and never a strong match.
    """

    max_alignments_per_candidate: int = 200
    """Alignments this client keeps per candidate when parsing.

    A separate number from :attr:`max_hits` because the units differ: that one counts database
    sequences, this one counts HSPs, and one sequence can align more than once. Headroom over
    ``max_hits`` on purpose.
    """

    def to_dict(self) -> dict[str, Any]:
        """Export the parameters."""
        return {
            "scope": self.scope,
            "database": self.database,
            "organism": self.organism,
            "species_restriction": self.species_restriction,
            "program": self.program,
            "task": self.task,
            "reward": self.reward,
            "penalty": self.penalty,
            "gap_costs": self.gap_costs,
            "word_size": self.word_size,
            "expect": self.expect,
            "filter_low_complexity": self.filter_low_complexity,
            "max_hits": self.max_hits,
            "max_alignments_per_candidate": self.max_alignments_per_candidate,
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
    on_target_hits: tuple[OffTargetHit, ...] = ()
    """Alignments removed as a candidate's own locus (ADR 0023).

    Removed from :attr:`hits` and from every count, kept here. "Dropped" must not mean
    "deleted": a reader has to be able to check that the removal was right, and a count that
    had something taken out of it has to be able to say so.
    """

    capped_candidates: tuple[str, ...] = ()
    """Candidates whose alignment list reached the per-candidate limit, so their counts are
    floors rather than totals."""

    window: tuple[int, int] | None = None
    """The run's own design window, which is what ``inside_run_window`` was measured against.

    On the record rather than only in the run, because the flag is meaningless without it: a
    reader of this JSON, or of the plot drawn from it, has to be able to see which window was
    used. ``None`` when the run had no window rule, in which case nothing was flagged.
    """

    def hits_for(self, candidate_id: str) -> tuple[OffTargetHit, ...]:
        """Return the hits belonging to one candidate."""
        return tuple(h for h in self.hits if h.candidate_id == candidate_id)

    def on_target_for(self, candidate_id: str) -> tuple[OffTargetHit, ...]:
        """Return the removed own-locus alignments belonging to one candidate."""
        return tuple(h for h in self.on_target_hits if h.candidate_id == candidate_id)

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
            "on_target_hits": [h.to_dict() for h in self.on_target_hits],
            "capped_candidates": list(self.capped_candidates),
            "window": list(self.window) if self.window else None,
        }
