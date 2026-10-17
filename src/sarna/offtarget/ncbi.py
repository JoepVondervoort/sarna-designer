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

import heapq
import os
import re
import threading
import time
import xml.etree.ElementTree as ElementTree
from collections.abc import Callable, Collection, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import httpx

from sarna import settings
from sarna.offtarget.models import (
    BY_DATABASE,
    GENOMIC_SCOPE,
    TRANSCRIPT_SCOPE,
    OffTargetHit,
    SearchParameters,
)

BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

TRANSCRIPT_DATABASE = "refseq_rna"
"""RefSeq RNA: could this duplex silence an mRNA? The classic siRNA off-target question."""

GENOMIC_DATABASE = "GPIPE/9606/current/ref_top_level"
"""The human GRCh38 reference assembly: could this duplex target another human promoter?

**Was `refseq_genomic`, which is not a database NCBI offers any more.** It appears nowhere in
NCBI's BLAST database directory and nowhere in the database list on NCBI's own blastn page.
The URL API still accepted it and still returned a request id, so the search ran, against
something that is certainly not 3.3 Gbp of human genome: a real 22 candidate run was still
unfinished after the full half-hour budget and returned nothing at all.

Measured sizes from NCBI's own metadata, which is the whole argument: this assembly is
3,272,089,205 bases in 639 sequences, where `refseq_rna` is 251,911,205,044 and `core_nt` is
1,086,004,251,965. Searching a 19-mer at word size 7 costs roughly in proportion, so pointing
the genomic scope at the one genome it is asking about is worth two orders of magnitude.

This is the name NCBI's own blastn page submits for "Human genomic plus transcript", minus the
transcript half, which this project searches separately as its other scope (ADR 0024).
"""

MINIMUM_POLL_SECONDS = 60.0
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


class BlastReportError(BlastError):
    """A submitted search exists but its downloaded report could not be interpreted."""


class BlastStillRunning(Exception):
    """The search has not finished. Not an error: poll again later."""


@dataclass(frozen=True, slots=True)
class ParsedReport:
    """Everything one BLAST report said, including where it was cut short."""

    hits: tuple[OffTargetHit, ...]
    capped: frozenset[str]
    queries: frozenset[str] | None = None
    """Query IDs actually present; None only for an explicit READY/no-hits reply."""
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
        # RefSeq RNA spans every species RefSeq holds, so this scope is human *by filter*.
        return SearchParameters(scope=scope, database=TRANSCRIPT_DATABASE)
    if scope == GENOMIC_SCOPE:
        # The assembly contains one species, so no filter is sent and none is needed. That
        # also drops a reliance on ENTREZ_QUERY, which the URL API does not document and which
        # its own page says should not be relied upon.
        return SearchParameters(
            scope=scope, database=GENOMIC_DATABASE, organism="",
            species_restriction=BY_DATABASE,
        )
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



class RequestScheduler:
    """Pace all clients in this process, including resumed RIDs and concurrent scopes."""

    def __init__(
        self, now: Callable[[], float] = time.monotonic,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        """Allow a virtual clock in network-free tests."""
        self.now, self.sleep = now, sleep
        self.lock = threading.Lock()
        self.last: float | None = None
        self.rids: dict[str, float] = {}

    def wait(self, payload: Mapping[str, str]) -> None:
        """Reserve one request slot; status polls also respect the per-RID minute."""
        with self.lock:
            current = self.now()
            due = current if self.last is None else self.last + 10.0
            rid = payload.get("RID", "")
            polling = payload.get("FORMAT_OBJECT") == "SearchInfo"
            if polling and rid in self.rids:
                due = max(due, self.rids[rid] + MINIMUM_POLL_SECONDS)
            if due > current:
                self.sleep(due - current)
            self.last = self.now()
            if polling:
                self.rids[rid] = self.last
            # RIDs expire; bound process memory independently of server lifetime.
            if len(self.rids) > 10000:
                self.rids = {key: at for key, at in self.rids.items()
                             if self.last - at < 86400}


_SCHEDULER = RequestScheduler()

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
            "NUCL_REWARD": str(parameters.reward),
            "NUCL_PENALTY": str(parameters.penalty),
            "GAPCOSTS": str(parameters.gap_costs),
            "tool": "sarna-designer",
        }
        email = os.environ.get("SARNA_NCBI_EMAIL", "").strip()
        if email:
            payload["email"] = email
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
        if found.group(1).upper() == "READY" and (
            hits is None or hits.group(1).lower() not in ("yes", "no")
        ):
            raise BlastError(
                "NCBI reported READY without an explicit hit status; no findings recorded."
            )
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
        # XML2 is the multiple-file format: NCBI returns a ZIP for multi-query jobs.
        # XML2_S returns every query in one XML document. Retrying downloads uses the
        # existing RID, never Put, and the shared scheduler paces both attempts.
        for attempt in range(2):
            raw = self._request_bytes(
                {"CMD": "Get", "RID": request_id, "FORMAT_TYPE": "XML2_S"}
            )
            try:
                return parse_report(raw, per_candidate=parameters.max_alignments_per_candidate,
                                    max_sequences=parameters.max_hits)
            except BlastReportError:
                if attempt:
                    raise
        raise BlastReportError("NCBI's report could not be retrieved.")

    def _request(self, payload: Mapping[str, str]) -> str:
        """Return a reply as text."""
        return self._request_bytes(payload).decode("utf-8", "replace")

    def _request_bytes(self, payload: Mapping[str, str]) -> bytes:
        """Post to the URL API and return the raw body.

        Raises:
            BlastError: On a transport failure, an error status, or an oversized body.
        """
        # Retry only retrievals: an interrupted submission may already have been accepted,
        # and retrying Put would create a duplicate search with an unknown RID.
        attempts = 3 if payload.get("CMD") == "Get" else 1
        for attempt in range(attempts):
            _SCHEDULER.wait(payload)
            try:
                with self._client.stream("POST", self.base_url, data=dict(payload),
                                         timeout=self.timeout) as response:
                    if response.status_code in (429, 502, 503, 504) and attempt < attempts - 1:
                        # Scheduling spaces retries as well. Honour a bounded Retry-After.
                        try:
                            retry_after = float(response.headers.get("Retry-After", "10"))
                        except ValueError:
                            retry_after = 10.0
                        _SCHEDULER.sleep(min(max(retry_after, 10), 120))
                        continue
                    if response.status_code >= 400:
                        raise BlastError(
                            f"NCBI BLAST returned HTTP {response.status_code}. The usage policy "
                            "may have been exceeded or the service may be unavailable. "
                            "The designs themselves are unaffected."
                        )
                    body = bytearray()
                    for chunk in response.iter_bytes():
                        if len(body) + len(chunk) > MAX_RESPONSE_BYTES:
                            raise BlastError("NCBI returned a report over the response size limit.")
                        body.extend(chunk)
                    return bytes(body)
            except httpx.HTTPError as error:
                if attempt < attempts - 1:
                    continue
                raise BlastError(
                    f"could not reach NCBI BLAST ({type(error).__name__}: {error}). Off-target "
                    "screening needs internet access; the designs themselves are unaffected."
                ) from error
        raise BlastError("NCBI retrieval could not be completed.")


def parse_hits(
    raw: bytes, max_hits: int = MAX_ALIGNMENTS_PER_CANDIDATE
) -> tuple[OffTargetHit, ...]:
    """Return just the alignments from a report. See :func:`parse_report`.

    ``max_hits`` is **per candidate**. A caller that needs to know which candidates were cut
    short wants :func:`parse_report` instead.
    """
    return parse_report(raw, per_candidate=max_hits).hits


def _normalise_xml2(root: ElementTree.Element) -> ElementTree.Element:
    """Read the documented XML2 format while retaining support for recorded classic XML."""
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]
    if root.tag not in ("BlastXML2", "BlastOutput2"):
        return root
    converted = ElementTree.Element("BlastOutput")
    iterations = ElementTree.SubElement(converted, "BlastOutput_iterations")
    for search in root.iter("Search"):
        iteration = ElementTree.SubElement(iterations, "Iteration")
        ElementTree.SubElement(iteration, "Iteration_query-def").text = _text(search, "query-title")
        hit_list = ElementTree.SubElement(iteration, "Iteration_hits")
        for source in search.findall("./hits/Hit"):
            hit = ElementTree.SubElement(hit_list, "Hit")
            description = source.find("./description/HitDescr")
            if description is None:
                raise BlastReportError("XML2 hit omitted its sequence description.")
            for name, tag in (("accession", "Hit_accession"), ("title", "Hit_def"),
                              ("id", "Hit_id")):
                ElementTree.SubElement(hit, tag).text = _text(description, name)
            hsps = ElementTree.SubElement(hit, "Hit_hsps")
            for source_hsp in source.findall("./hsps/Hsp"):
                hsp = ElementTree.SubElement(hsps, "Hsp")
                for field in source_hsp:
                    # Preserve whitespace in midline: a leading space is a mismatch.
                    ElementTree.SubElement(hsp, "Hsp_" + field.tag).text = field.text
    return converted


def parse_report(
    raw: bytes, per_candidate: int = MAX_ALIGNMENTS_PER_CANDIDATE,
    max_sequences: int | None = None,
) -> ParsedReport:
    """Parse BLAST XML into hits, attributing each to the candidate that produced it.

    Args:
        raw: The XML body, treated as untrusted input. A DOCTYPE is accepted (every
            genuine NCBI report has one) but an entity declaration is not.
        max_sequences: Server sequence limit; reaching it marks counts as potentially capped.
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
    if raw.startswith(b"PK"):
        raise BlastReportError(
            "NCBI returned a ZIP archive instead of single-file XML. "
            "The existing search can be retrieved again using XML2_S."
        )
    if per_candidate <= 0:
        raise ValueError(
            f"the per-candidate alignment budget must be positive, got {per_candidate}; "
            f"zero would report every candidate as having no alignments"
        )
    if _ENTITY_DECLARATION.search(raw):
        raise BlastReportError(
            "the BLAST reply declares XML entities, which this client refuses to parse because "
            "expanding them can exhaust memory"
        )
    try:
        root = ElementTree.fromstring(raw)
    except ElementTree.ParseError as error:
        raise BlastReportError(f"could not parse NCBI's XML reply: {error}") from error

    root = _normalise_xml2(root)
    if root.tag != "BlastOutput" or root.find("BlastOutput_iterations") is None:
        raise BlastReportError(
            "NCBI returned an unexpected BLAST report format; no findings recorded."
        )
    iterations = list(root.iter("Iteration"))
    if not iterations:
        raise BlastReportError("NCBI's report contains no query results; no findings recorded.")
    queries: set[str] = set()
    # A bounded best-of per candidate rather than the first N encountered. BLAST orders a
    # report by database sequence, and on a chromosome-level database one sequence carries
    # thousands of alignments, so taking them in document order fills the budget with
    # chromosome 1's weak partial matches and never reaches chromosome 7's strong one. The
    # heap keeps the highest bit scores and stays bounded, which is what the cap is for.
    best: dict[str, list[tuple[float, int, OffTargetHit]]] = {}
    order: list[str] = []
    capped: set[str] = set()
    seen = 0
    for iteration in iterations:
        candidate = _text(iteration, "Iteration_query-def") or ""
        # BLAST keeps only the first token of a FASTA header as the id in some views.
        candidate = candidate.split()[0] if candidate.split() else candidate
        if not candidate or candidate in queries:
            raise BlastReportError("NCBI returned a missing or duplicate query identifier.")
        queries.add(candidate)
        if candidate not in best:
            best[candidate] = []
            order.append(candidate)
        sequence_hits = list(iteration.iter("Hit"))
        if max_sequences is not None and len(sequence_hits) >= max_sequences:
            capped.add(candidate)
        for hit in sequence_hits:
            accession = _text(hit, "Hit_accession") or _text(hit, "Hit_id") or "unknown"
            title = _text(hit, "Hit_def") or ""
            for hsp in hit.iter("Hsp"):
                frame = _int(hsp, "Hsp_hit-frame", 1)
                query_from = _int(hsp, "Hsp_query-from", 0)
                aligned = _text(hsp, "Hsp_qseq")
                midline_node = hsp.find("Hsp_midline")
                midline = midline_node.text if midline_node is not None else None
                detailed = bool(aligned and midline)
                parsed = (
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
                        subject_strand=(
                            "minus" if frame < 0
                            or (_text(hsp, "Hsp_hit-strand") or "").lower() == "minus"
                            else "plus"
                        ),
                        evalue=_float(hsp, "Hsp_evalue", float("inf")),
                        bit_score=_float(hsp, "Hsp_bit-score", 0.0),
                        mismatch_positions=(
                            mismatch_positions(aligned or "", midline or "", query_from)
                            if detailed else ()
                        ),
                        has_alignment_detail=detailed,
                        query_sequence=aligned or "",
                        subject_sequence=_text(hsp, "Hsp_hseq") or "",
                    )
                )
                # Per candidate, so a promiscuous query cannot consume a later query's budget.
                seen += 1
                held = best[candidate]
                heapq.heappush(held, (parsed.bit_score, -seen, parsed))
                if len(held) > per_candidate:
                    heapq.heappop(held)
                    capped.add(candidate)
    hits: list[OffTargetHit] = []
    for candidate in order:
        # Strongest first within each candidate, and candidates in the order BLAST reported
        # them, so the output is stable rather than heap order.
        hits.extend(
            parsed for _score, _tie, parsed in sorted(best[candidate], reverse=True)
        )
    return ParsedReport(hits=tuple(hits), capped=frozenset(capped), queries=frozenset(queries))


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
