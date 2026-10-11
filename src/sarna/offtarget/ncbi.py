"""A client for NCBI's BLAST URL API, shaped for 19 nt queries.

Three things about this interface drive the design.

**Short queries need explicit parameters.** At 19 nt the default `blastn` word size (11 or 28)
and E-value threshold return *nothing*, which is indistinguishable from a clean result and is
exactly the misreading this project exists to prevent. The web interface calls the right
configuration ``blastn-short``; ``TASK`` is not a documented URL API parameter, so the
configuration is reproduced explicitly, ``MEGABLAST=off``, ``WORD_SIZE=7``, a generous
``EXPECT`` and low-complexity filtering off. :class:`~sarna.offtarget.models.SearchParameters`
records what was actually sent, so a reader never has to guess.

**It is asynchronous by nature.** ``CMD=Put`` returns a request id; results arrive later. That
suits the accepted design (the researcher reads their candidates while the search runs) and it
is why nothing here blocks: :meth:`submit` and :meth:`poll` are separate calls.

**NCBI's usage policy is a shared resource, not a rate limit to tune.** Exceeding it can get the
caller's whole institution blocked, so every candidate in a run goes in **one** multi-FASTA job
rather than one job each, polling never runs faster than
:data:`MINIMUM_POLL_SECONDS`, and the client refuses to poll sooner.

Responses are untrusted external input. The XML is size-capped and rejected outright if it
declares a DTD or entities, because `xml.etree` expands entities and a crafted response could
otherwise exhaust memory.
"""

import re
import xml.etree.ElementTree as ElementTree
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import httpx

from sarna import settings
from sarna.offtarget.models import (
    GENOMIC_SCOPE,
    TRANSCRIPT_SCOPE,
    OffTargetHit,
    SearchParameters,
)

BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

TRANSCRIPT_DATABASE = "refseq_rna"
"""RefSeq RNA: could this duplex silence an mRNA? The classic siRNA off-target question."""

GENOMIC_DATABASE = "refseq_genomic"
"""RefSeq genomic: could this duplex target another promoter? The saRNA-specific question."""

MINIMUM_POLL_SECONDS = 15.0
"""NCBI asks callers not to poll more often than this. Polling faster risks an IP block."""

MAX_RESPONSE_BYTES = 64 * 1024 * 1024
"""Refuse a response larger than this rather than parse it."""

MAX_ALIGNMENTS_PER_CANDIDATE = 200
"""Alignments kept **per candidate**, not per job.

This was one shared budget of 500 for the whole multi-query report, and the bug it caused is
the reason ADR 0023 exists: the first query filled it, the parser returned, and the other 21
candidates of a real run were reported with zero alignments each. They had not been searched
clean, they had not been read at all, and the table printed their zeros in the same ink as a
genuine result.

A per-candidate budget cannot starve a later query, and a candidate that reaches it is recorded
as capped so "nothing found" and "more than we kept" stay different statements. BLAST orders a
query's alignments best first, so what a cap drops is always the weakest.
"""

DEFAULT_TIMEOUT_SECONDS = 60.0

_RID = re.compile(r"^\s*RID\s*=\s*(\S+)", re.MULTILINE)
_RTOE = re.compile(r"^\s*RTOE\s*=\s*(\d+)", re.MULTILINE)
_STATUS = re.compile(r"^\s*Status\s*=\s*(\S+)", re.MULTILINE)
_THERE_ARE_HITS = re.compile(r"^\s*ThereAreHits\s*=\s*(\S+)", re.MULTILINE)
_ENTITY_DECLARATION = re.compile(rb"<!ENTITY", re.IGNORECASE)
"""The actual attack vector, and the only one worth refusing.

NCBI's real replies open with a ``<!DOCTYPE BlastOutput PUBLIC ... NCBI_BlastOutput.dtd>``, so
refusing every DOCTYPE (as an earlier version of this guard did) rejects every genuine report.
What makes "billion laughs" and quadratic blowup possible is an *internal entity declaration*;
``xml.etree`` never fetches an external DTD, so a bare DOCTYPE expands nothing."""


class BlastError(RuntimeError):
    """A BLAST search could not be run or its result could not be read."""


class BlastStillRunning(Exception):
    """The search has not finished. Not an error: poll again later."""


@dataclass(frozen=True, slots=True)
class ParsedReport:
    """Everything one BLAST report said, including where it was cut short."""

    hits: tuple[OffTargetHit, ...]
    capped: frozenset[str]
    """Candidates whose alignment list reached the per-candidate limit."""


@dataclass(frozen=True, slots=True)
class SubmittedSearch:
    """A job NCBI has accepted but not yet finished."""

    request_id: str
    estimated_seconds: int
    parameters: SearchParameters

    def to_dict(self) -> dict[str, Any]:
        """Export the submission."""
        return {
            "request_id": self.request_id,
            "estimated_seconds": self.estimated_seconds,
            "parameters": self.parameters.to_dict(),
        }


def parameters_for(scope: str) -> SearchParameters:
    """Return the search parameters for one scope.

    Raises:
        ValueError: If the scope is not one this project searches.
    """
    if scope == TRANSCRIPT_SCOPE:
        return SearchParameters(scope=scope, database=TRANSCRIPT_DATABASE)
    if scope == GENOMIC_SCOPE:
        return SearchParameters(scope=scope, database=GENOMIC_DATABASE)
    raise ValueError(
        f"unknown search scope {scope!r}; expected {TRANSCRIPT_SCOPE!r} or {GENOMIC_SCOPE!r}"
    )


def as_fasta(sequences: Mapping[str, str]) -> str:
    """Return one multi-FASTA document, so a whole run is a single BLAST job.

    Raises:
        ValueError: If there is nothing to search, or an identifier contains whitespace that
            would split the FASTA header and silently rename the query.
    """
    if not sequences:
        raise ValueError("there are no sequences to search")
    for name in sequences:
        if not name or any(c.isspace() for c in name):
            raise ValueError(
                f"candidate identifier {name!r} contains whitespace; BLAST would truncate it "
                f"and the hits could not be attributed back to the candidate"
            )
    return "\n".join(f">{name}\n{bases}" for name, bases in sequences.items())


class BlastClient:
    """Submits and retrieves NCBI BLAST searches over the URL API."""

    __slots__ = ("_client", "_owned", "base_url", "timeout")

    def __init__(
        self,
        base_url: str = BASE_URL,
        timeout: float = DEFAULT_TIMEOUT_SECONDS,
        client: httpx.Client | None = None,
    ) -> None:
        """Build a client, optionally over a supplied transport for testing."""
        self.base_url = base_url
        self.timeout = timeout
        self._owned = client is None
        self._client = client or httpx.Client(
            timeout=timeout, headers={"User-Agent": settings.user_agent()}
        )

    def close(self) -> None:
        """Close the transport if this client created it."""
        if self._owned:
            self._client.close()

    def submit(self, sequences: Mapping[str, str], scope: str) -> SubmittedSearch:
        """Submit every sequence as one job and return its request id.

        Raises:
            BlastError: If NCBI did not accept the job or returned no request id.
            ValueError: If the scope or the sequences are unusable.
        """
        parameters = parameters_for(scope)
        payload: dict[str, str] = {
            "CMD": "Put",
            "PROGRAM": parameters.program,
            "DATABASE": parameters.database,
            # Explicitly off: megablast is tuned for long, highly similar sequences and finds
            # nothing at 19 nt.
            "MEGABLAST": "off",
            "WORD_SIZE": str(parameters.word_size),
            "EXPECT": str(parameters.expect),
            "FILTER": "F" if not parameters.filter_low_complexity else "T",
            "HITLIST_SIZE": str(parameters.max_hits),
            "QUERY": as_fasta(sequences),
        }
        if parameters.organism:
            # This project designs against the human genome only, so a hit in another species
            # is noise that would swamp every measure.
            payload["ENTREZ_QUERY"] = parameters.organism
        text = self._request(payload)
        rid = _RID.search(text)
        if rid is None:
            raise BlastError(
                "NCBI accepted the request but returned no RID, so the search cannot be "
                "retrieved. This usually means the query or database was rejected. The first "
                f"200 characters of the reply were: {text.strip()[:200]!r}"
            )
        estimate = _RTOE.search(text)
        return SubmittedSearch(
            request_id=rid.group(1),
            estimated_seconds=int(estimate.group(1)) if estimate else 0,
            parameters=parameters,
        )

    def status(self, request_id: str) -> tuple[str, bool]:
        """Return ``(status, there_are_hits)`` for a submitted search.

        Status is NCBI's own word: ``READY``, ``WAITING`` or ``UNKNOWN``. ``UNKNOWN`` means the
        request id has expired or never existed, which is a failure rather than a wait.
        """
        text = self._request(
            {"CMD": "Get", "RID": request_id, "FORMAT_OBJECT": "SearchInfo"}
        )
        found = _STATUS.search(text)
        if found is None:
            raise BlastError(
                f"NCBI's status reply for {request_id} contained no Status line: "
                f"{text.strip()[:200]!r}"
            )
        hits = _THERE_ARE_HITS.search(text)
        return found.group(1).upper(), bool(hits and hits.group(1).lower() == "yes")

    def fetch(self, request_id: str, parameters: SearchParameters) -> ParsedReport:
        """Retrieve a finished search's alignments, and which candidates were cut short.

        Raises:
            BlastStillRunning: If the search has not finished. Poll again, no sooner than
                :data:`MINIMUM_POLL_SECONDS`.
            BlastError: If the search failed or its request id has expired.
        """
        state, there_are_hits = self.status(request_id)
        if state == "WAITING":
            raise BlastStillRunning(request_id)
        if state == "UNKNOWN":
            raise BlastError(
                f"NCBI no longer recognises search {request_id}. Request ids expire after about "
                f"24 hours; run the search again."
            )
        if state != "READY":
            raise BlastError(f"NCBI reported search {request_id} as {state}.")
        if not there_are_hits:
            return ParsedReport(hits=(), capped=frozenset())
        raw = self._request_bytes(
            {"CMD": "Get", "RID": request_id, "FORMAT_TYPE": "XML"}
        )
        return parse_report(raw, per_candidate=parameters.max_alignments_per_candidate)

    def _request(self, payload: Mapping[str, str]) -> str:
        """Return a reply as text."""
        return self._request_bytes(payload).decode("utf-8", "replace")

    def _request_bytes(self, payload: Mapping[str, str]) -> bytes:
        """Post to the URL API and return the raw body.

        Raises:
            BlastError: On a transport failure, an error status, or an oversized body.
        """
        try:
            response = self._client.post(self.base_url, data=dict(payload), timeout=self.timeout)
        except httpx.HTTPError as error:
            raise BlastError(
                f"could not reach NCBI BLAST ({type(error).__name__}: {error}). Off-target "
                f"screening needs internet access; the designs themselves are unaffected."
            ) from error
        if response.status_code >= 400:
            raise BlastError(
                f"NCBI BLAST returned HTTP {response.status_code}. If this persists, the usage "
                f"policy may have been exceeded: searches are limited in rate and large jobs "
                f"are asked to run outside US business hours."
            )
        body = response.content
        if len(body) > MAX_RESPONSE_BYTES:
            raise BlastError(
                f"NCBI returned {len(body):,} bytes, over the {MAX_RESPONSE_BYTES:,} byte limit "
                f"this client will parse."
            )
        return body


def parse_hits(
    raw: bytes, max_hits: int = MAX_ALIGNMENTS_PER_CANDIDATE
) -> tuple[OffTargetHit, ...]:
    """Return just the alignments from a report. See :func:`parse_report`.

    ``max_hits`` is **per candidate**. A caller that needs to know which candidates were cut
    short wants :func:`parse_report` instead.
    """
    return parse_report(raw, per_candidate=max_hits).hits


def parse_report(
    raw: bytes, per_candidate: int = MAX_ALIGNMENTS_PER_CANDIDATE
) -> ParsedReport:
    """Parse BLAST XML into hits, attributing each to the candidate that produced it.

    Args:
        raw: The XML body, treated as untrusted input. A DOCTYPE is accepted (every
            genuine NCBI report has one) but an entity declaration is not.
        per_candidate: Alignments kept for each query. A budget shared across the whole report
            lets one promiscuous query starve every later one, which is the defect ADR 0023
            was written for; see :data:`MAX_ALIGNMENTS_PER_CANDIDATE`.

    Returns:
        Every HSP as its own hit, because two alignments to the same accession are two separate
        places the duplex could bind, plus the candidates whose list reached the limit.

    Raises:
        BlastError: If the document declares entities, or is not parseable. `xml.etree`
            expands internal entities, so a crafted reply could otherwise exhaust memory.
        ValueError: If the per-candidate budget is not positive.
    """
    if per_candidate <= 0:
        raise ValueError(
            f"the per-candidate alignment budget must be positive, got {per_candidate}; "
            f"zero would report every candidate as having no alignments"
        )
    if _ENTITY_DECLARATION.search(raw):
        raise BlastError(
            "the BLAST reply declares XML entities, which this client refuses to parse because "
            "expanding them can exhaust memory"
        )
    try:
        root = ElementTree.fromstring(raw)
    except ElementTree.ParseError as error:
        raise BlastError(f"could not parse NCBI's XML reply: {error}") from error

    hits: list[OffTargetHit] = []
    kept: dict[str, int] = {}
    capped: set[str] = set()
    for iteration in root.iter("Iteration"):
        candidate = _text(iteration, "Iteration_query-def") or ""
        # BLAST keeps only the first token of a FASTA header as the id in some views.
        candidate = candidate.split()[0] if candidate.split() else candidate
        for hit in iteration.iter("Hit"):
            # Per candidate, so a promiscuous query cannot consume a later query's budget.
            # Moving on to the next iteration rather than returning is the whole fix.
            if kept.get(candidate, 0) >= per_candidate:
                capped.add(candidate)
                break
            accession = _text(hit, "Hit_accession") or _text(hit, "Hit_id") or "unknown"
            title = _text(hit, "Hit_def") or ""
            for hsp in hit.iter("Hsp"):
                if kept.get(candidate, 0) >= per_candidate:
                    capped.add(candidate)
                    break
                kept[candidate] = kept.get(candidate, 0) + 1
                frame = _int(hsp, "Hsp_hit-frame", 1)
                query_from = _int(hsp, "Hsp_query-from", 0)
                aligned = _text(hsp, "Hsp_qseq")
                midline = _text(hsp, "Hsp_midline")
                detailed = bool(aligned and midline)
                hits.append(
                    OffTargetHit(
                        candidate_id=candidate,
                        accession=accession,
                        title=title,
                        identities=_int(hsp, "Hsp_identity", 0),
                        alignment_length=_int(hsp, "Hsp_align-len", 0),
                        gaps=_int(hsp, "Hsp_gaps", 0),
                        query_start=query_from,
                        query_end=_int(hsp, "Hsp_query-to", 0),
                        subject_start=_int(hsp, "Hsp_hit-from", 0),
                        subject_end=_int(hsp, "Hsp_hit-to", 0),
                        subject_strand="minus" if frame < 0 else "plus",
                        evalue=_float(hsp, "Hsp_evalue", float("inf")),
                        bit_score=_float(hsp, "Hsp_bit-score", 0.0),
                        mismatch_positions=(
                            mismatch_positions(aligned or "", midline or "", query_from)
                            if detailed else ()
                        ),
                        has_alignment_detail=detailed,
                    )
                )
    return ParsedReport(hits=tuple(hits), capped=frozenset(capped))


def mismatch_positions(query: str, midline: str, query_start: int) -> tuple[int, ...]:
    """Return the 1-based query positions that did not match, from BLAST's alignment strings.

    ``Hsp_midline`` carries a bar under every matching column and a space under every other
    one, and ``Hsp_qseq`` carries the query with gap characters inserted. Walking them together
    converts column positions into query positions, which is what a reader needs: an identity
    count says seventeen of nineteen matched and never says *which* two did not.

    A column where the query itself is gapped belongs to no query position and is attributed to
    the position before it, so an insertion in the subject is still marked somewhere rather
    than dropped.
    """
    found: list[int] = []
    position = query_start - 1
    for index, base in enumerate(query):
        if base != "-":
            position += 1
        matched = index < len(midline) and midline[index] == "|"
        if not matched and position >= query_start and (
            not found or found[-1] != position
        ):
            found.append(position)
    return tuple(found)


def _text(element: ElementTree.Element, tag: str) -> str | None:
    """Return a child's text, or ``None``."""
    found = element.find(tag)
    return None if found is None or found.text is None else found.text.strip()


def _int(element: ElementTree.Element, tag: str, default: int) -> int:
    """Return a child's text as an integer, or ``default``."""
    raw = _text(element, tag)
    try:
        return default if raw is None else int(raw)
    except ValueError:
        return default


def _float(element: ElementTree.Element, tag: str, default: float) -> float:
    """Return a child's text as a float, or ``default``."""
    raw = _text(element, tag)
    try:
        return default if raw is None else float(raw)
    except ValueError:
        return default


def candidate_sequences(
    records: Sequence[Mapping[str, Any]],
    limit: int | None = None,
    only: Collection[str] | None = None,
) -> dict[str, str]:
    """Return the target DNA of each passing candidate, keyed by candidate id, in rank order.

    The **target site DNA** is searched, not the RNA core: an off-target search asks where else
    this sequence occurs in the genome or transcriptome, which is a question about DNA. Only
    candidates that passed the gates are included, because a duplex nobody would order does not
    need screening.

    Args:
        records: The run's candidates.
        limit: Keep at most this many, the least promising dropped first.
        only: Screen exactly these candidates. Every batch the interface offers, a count, the
            next fifteen, the ticked selection, resolves to one of these lists, so there is one
            selection path rather than one per button (ADR 0023). Ids that are unknown or did
            not pass the gates are ignored here and rejected at the boundary, where a researcher
            can be told which ones.
    """
    passing = [r for r in records if r.get("rank") is not None]
    if only is not None:
        wanted = set(only)
        passing = [r for r in passing if r.get("id") in wanted]
    passing.sort(key=lambda r: r["rank"])
    if limit is not None:
        passing = passing[:limit]
    return {r["id"]: r["target_site_dna"] for r in passing}
