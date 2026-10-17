"""Running a screening end to end: submit, wait, place, measure.

This is the piece that turns a finished design run into a dated screening record. It is written
to be called from a background worker, and to be interruptible in the sense that matters: every
stage it reaches is reported, so the interface can say *what* it is waiting for instead of
spinning.

Sleeping is injected rather than imported so tests do not wait fifteen real seconds per poll.
"""

import re
import time
from collections.abc import Callable, Collection, Mapping, Sequence
from typing import Any

from sarna.offtarget.models import (
    GENOMIC_SCOPE,
    PromoterProximity,
    ScreeningRecord,
    SearchParameters,
)
from sarna.offtarget.ncbi import (
    MINIMUM_POLL_SECONDS,
    BlastClient,
    BlastError,
    BlastStillRunning,
    SubmittedSearch,
    candidate_sequences,
    parameters_for,
)
from sarna.offtarget.ontarget import loci_from_candidates, split_on_target
from sarna.offtarget.proximity import TranscriptLookup, place_hit
from sarna.offtarget.specificity import measure_specificity

MAX_CANDIDATES = 50
"""Most candidates to screen in **one job**.

A courtesy to NCBI rather than a technical limit: a run can generate thousands of candidates,
and only the ones a researcher might order are worth a shared public service's time. Ranked
order, so the cut falls on the least promising.

Since ADR 0023 this is a per-job ceiling rather than the whole allowance. A researcher asking
for more gets several sequential jobs, each its own dated record, never parallel submissions.
"""

DEFAULT_CANDIDATES = 15
"""How many a screening covers when nobody says otherwise.

The old fixed 50 is what made a genomic search time out after half an hour and return nothing.
Fifteen finishes, and the next fifteen are one button away.
"""

MAX_WAIT_SECONDS = 1800.0
"""Stop *waiting* after half an hour.

Not "give up": NCBI holds a request id for about a day, so a search that overruns this is
still there to be collected. :exc:`SearchStillAtNcbi` carries the id out so the interface can
offer to come back for it, which is the difference between a slow answer and no answer.
"""

MAX_ANNOTATION_SECONDS = 60.0
"""Soft budget for optional TSS annotation, checked between remote lookups."""

MAX_PLACEMENTS = 300
"""Most genomic hits placed against annotated start sites in one screening.

Placement is one Ensembl request per uncached megabase block, measured at 1,022 requests for
3,000 scattered hits, so an unbounded screening spends minutes on alignments nobody will read.
The strongest are placed first and the rest are reported without a distance, which the record
states rather than leaving the difference invisible.
"""

ProgressReporter = Callable[[str, "Mapping[str, Any] | None"], None]
"""Called with a sentence and, where there is one, the countable state behind it.

A sentence cannot drive a progress bar and a number cannot be read aloud, so both travel.
"""


def _silent(_stage: str, _detail: "Mapping[str, Any] | None" = None) -> None:
    """Swallow progress when nobody asked for it."""


class SearchStillAtNcbi(Exception):
    """The wait budget ran out while the search was still running.

    Not an error. The search exists, NCBI is still working on it, and the request id is the
    handle to collect it later.
    """

    def __init__(
        self,
        request_id: str,
        waited: float,
        searched: Sequence[str],
        estimate_seconds: int = 0,
    ) -> None:
        """Carry everything a resume needs, and NCBI's own estimate if it gave one."""
        super().__init__(
            
                f"NCBI estimates about {_plainly(estimate_seconds)} for search {request_id}, "
                f"which is longer than this app will sit and wait. The search is queued at "
                f"NCBI and can be collected when it is done; the designs are unaffected."
                if estimate_seconds > waited
                else
                f"NCBI search {request_id} has run for {_plainly(waited)} without finishing. "
                f"It is still queued there and can be collected when it is done; the designs "
                f"are unaffected."
            
        )
        self.request_id = request_id
        self.waited = waited
        self.searched = tuple(searched)
        self.estimate_seconds = estimate_seconds


def _plainly(seconds: float) -> str:
    """Return a duration a reader can judge at a glance rather than count zeros in."""
    seconds = int(seconds)
    if seconds < 90:
        return f"{seconds} seconds"
    if seconds < 5400:
        return f"{round(seconds / 60)} minutes"
    return f"{seconds / 3600:.1f} hours"


def run_screening(
    run_id: str,
    candidates: Sequence[Mapping[str, Any]],
    scope: str,
    core_length: int,
    client: BlastClient,
    progress: ProgressReporter | None = None,
    transcripts: TranscriptLookup | None = None,
    window: tuple[int, int] | None = None,
    sleep: Callable[[float], None] = time.sleep,
    now: Callable[[], float] = time.monotonic,
    poll_seconds: float = MINIMUM_POLL_SECONDS,
    max_wait: float = MAX_WAIT_SECONDS,
    only: Collection[str] | None = None,
    limit: int = MAX_CANDIDATES,
) -> ScreeningRecord:
    """Screen some of one run's passing candidates and return the dated record.

    Args:
        run_id: The design run being screened.
        candidates: The run's candidate records.
        scope: ``transcript`` or ``genomic``.
        core_length: Duplex core length, for judging a match complete.
        client: The BLAST client.
        progress: Called with each stage, for the interface.
        transcripts: Looks up annotated transcripts near a genomic hit. Required for the
            genomic scope; ignored otherwise.
        window: The run's own design window, used only to flag hits that fall inside it.
        sleep: Injected so tests need not wait.
        now: Monotonic clock, injected for the same reason.
        poll_seconds: Never below NCBI's requested minimum.
        max_wait: Stop local collection after this long; the RID can be resumed.
        only: Screen exactly these candidates rather than the best ``limit`` of them.
        limit: Most candidates in this one job, never above :data:`MAX_CANDIDATES`.

    Returns:
        The screening, including the parameters and the moment it was made.

    Raises:
        BlastError: If the search could not be run or did not finish in time.
        ValueError: If there is nothing to screen, or the genomic scope has no lookup.
    """
    say = progress or _silent
    if scope == GENOMIC_SCOPE and transcripts is None:
        raise ValueError(
            "the genomic scope places each hit against annotated start sites, so it needs a "
            "transcript lookup; without one a hit would be a coordinate with no meaning"
        )
    sequences = candidate_sequences(
        candidates, limit=min(limit, MAX_CANDIDATES), only=only
    )
    if not sequences:
        raise ValueError(
            "no candidate passed the gates in this run, so there is nothing to screen"
            if only is None
            else "none of the chosen candidates passed the gates in this run"
        )

    say(f"Submitting {len(sequences)} sequences to NCBI",
        {"phase": "submitting", "candidates": len(sequences)})
    submitted = client.submit(sequences, scope)
    say(f"NCBI accepted the search (id {submitted.request_id})",
        {"phase": "waiting", "request_id": submitted.request_id,
         "estimate_seconds": submitted.estimated_seconds, "elapsed_seconds": 0,
         "candidates": len(sequences), "candidates_searched": list(sequences),
         "parameters": submitted.parameters.to_dict()})
    return collect_screening(
        run_id=run_id, candidates=candidates, scope=scope, core_length=core_length,
        client=client, request_id=submitted.request_id, searched=tuple(sequences),
        estimate_seconds=submitted.estimated_seconds, progress=progress,
        transcripts=transcripts, window=window, sleep=sleep, now=now,
        poll_seconds=poll_seconds, max_wait=max_wait,
    )


def collect_screening(
    run_id: str,
    candidates: Sequence[Mapping[str, Any]],
    scope: str,
    core_length: int,
    client: BlastClient,
    request_id: str,
    searched: Sequence[str],
    estimate_seconds: int = 0,
    progress: ProgressReporter | None = None,
    transcripts: TranscriptLookup | None = None,
    window: tuple[int, int] | None = None,
    sleep: Callable[[float], None] = time.sleep,
    now: Callable[[], float] = time.monotonic,
    poll_seconds: float = MINIMUM_POLL_SECONDS,
    max_wait: float = MAX_WAIT_SECONDS,
    waited_already: float = 0.0,
    parameters: SearchParameters | None = None,
) -> ScreeningRecord:
    """Wait for a submitted search, read it, place it, and build the dated record.

    Separate from submission so a search that outlived its wait budget can be collected later
    from its request id instead of being thrown away and run again (ADR 0024).

    Args:
        run_id: The design run being screened.
        candidates: The run's candidate records, for the on-target split.
        scope: ``transcript`` or ``genomic``.
        core_length: Duplex core length, for judging a match complete.
        client: The BLAST client.
        request_id: NCBI's id for a search already submitted.
        searched: The candidates that went into that search, in the order they were sent.
        estimate_seconds: NCBI's own RTOE for this search, for the first wait and the bar.
        progress: Called with each stage and its countable state, for the interface.
        transcripts: Looks up annotated transcripts near a genomic hit.
        window: The run's own design window, used only to flag hits that fall inside it.
        sleep: Injected so tests need not wait.
        now: Monotonic clock, injected for the same reason.
        poll_seconds: Never below NCBI's requested minimum.
        max_wait: Stop waiting after this long.
        parameters: Original submitted parameters when resuming; defaults to current settings.
        waited_already: Seconds an earlier attempt spent waiting, so the interface can show
            the total rather than restarting the clock.

    Raises:
        SearchStillAtNcbi: If the budget runs out while NCBI is still working.
        BlastError: If the search failed or its request id has expired.
    """
    say = progress or _silent
    parameters = parameters or parameters_for(scope)
    sequences = {name: "" for name in searched}

    deadline = now() + max_wait
    interval = max(poll_seconds, MINIMUM_POLL_SECONDS)
    submitted = SubmittedSearch(
        request_id=request_id, estimated_seconds=estimate_seconds, parameters=parameters,
    )
    # RTOE is a forecast, not a lower bound. Check at the policy interval even when
    # the estimate is large, so an early result does not wait fifteen unnecessary minutes.
    beyond_budget = submitted.estimated_seconds > max_wait
    wait_next = interval
    started = now()
    waited = waited_already
    while True:
        sleep(wait_next)
        waited = waited_already + max(now() - started, waited - waited_already + wait_next)
        # Continue at the policy interval regardless of the forecast.
        wait_next = interval
        try:
            report = client.fetch(submitted.request_id, parameters)
            break
        except BlastStillRunning:
            say(
                f"NCBI is still searching ({int(waited)}s elapsed)",
                {"phase": "waiting", "request_id": submitted.request_id,
                 "estimate_seconds": submitted.estimated_seconds,
                 "elapsed_seconds": int(waited), "candidates": len(sequences)},
            )
            if beyond_budget or now() > deadline:
                raise SearchStillAtNcbi(
                    submitted.request_id, waited, tuple(sequences),
                    estimate_seconds=submitted.estimated_seconds,
                ) from None

    if report.queries is not None and set(searched) != set(report.queries):
        missing = set(searched) - set(report.queries)
        unexpected = set(report.queries) - set(searched)
        raise BlastError(
            f"BLAST query coverage mismatch: missing {sorted(missing)}, unexpected "
            f"{sorted(unexpected)}. No zero-hit findings were inferred."
        )
    say(f"Reading {len(report.hits)} alignments",
        {"phase": "reading", "alignments": len(report.hits),
         "candidates": len(sequences)})
    # The candidate's own promoter is a perfect match to itself, so without this every genomic
    # candidate reports at least one perfect hit and the column means nothing (ADR 0023).
    # Only the genomic scope: a promoter is usually not transcribed, and the only handle in a
    # transcript report would be a gene name, which is a guess.
    loci = loci_from_candidates(candidates) if scope == GENOMIC_SCOPE else {}
    hits, own = split_on_target(report.hits, loci)

    proximities: tuple[PromoterProximity, ...] = ()
    unplaced = 0
    if scope == GENOMIC_SCOPE and transcripts is not None:
        # Strongest first, and bounded: placement costs an Ensembl request per uncached block,
        # and the weakest alignments are the ones nobody will open.
        # Round-robin strongest hits per candidate so one repetitive query cannot consume
        # every annotation slot. Raw alignments are preserved independently of this budget.
        groups = {name: sorted((h for h in hits if h.candidate_id == name),
                              key=lambda h: h.bit_score, reverse=True) for name in searched}
        ordered = [group[index] for index in range(max(map(len, groups.values()), default=0))
                   for group in groups.values() if index < len(group)]
        to_place, rest = ordered[:MAX_PLACEMENTS], ordered[MAX_PLACEMENTS:]
        unplaced = len(rest)
        placed: list[PromoterProximity] = []
        annotation_deadline = now() + MAX_ANNOTATION_SECONDS
        for index, hit in enumerate(to_place, start=1):
            if now() >= annotation_deadline:
                unplaced += len(to_place) - len(placed)
                break
            if index == 1 or index % 20 == 0 or index == len(to_place):
                say(
                    f"Placing hits against annotated start sites ({index} of "
                    f"{len(to_place)})",
                    {"phase": "placing", "placed": index, "to_place": len(to_place),
                     "candidates": len(sequences)},
                )
            placed.append(place_hit(hit, transcripts, window))
        proximities = tuple(placed)

    passing = sum(1 for c in candidates if c.get("rank") is not None)
    notes = []
    if len(sequences) < passing:
        notes.append(
            f"{len(sequences)} of the {passing} candidates that passed the gates were screened "
            f"in this search. The rest have not been searched, which is not the same as having "
            f"been searched clean."
        )
    if report.capped:
        notes.append(
            f"{len(report.capped)} candidate"
            f"{'s' if len(report.capped) != 1 else ''} reached a reporting limit "
            f"({parameters.max_hits} database sequences or "
            f"{parameters.max_alignments_per_candidate} retained alignments each); "
            "counts may be floors rather than totals."
        )
    if own:
        notes.append(
            f"{len(own)} alignment{'s' if len(own) != 1 else ''} at the candidates' own target "
            f"loci {'were' if len(own) != 1 else 'was'} removed from these findings. A sequence "
            f"matching the promoter it was designed against is the intended target."
        )
    if scope == GENOMIC_SCOPE and not loci:
        notes.append(
            "This run has no genomic coordinates, so a hit at a candidate's own target site "
            "could not be told apart from one elsewhere and none were removed."
        )
    unavailable = sum(p.annotation_status == "unavailable" for p in proximities)
    if unavailable:
        notes.append(f"Transcript annotation was unavailable for {unavailable} alignments; "
                     "their TSS context is unknown, not absent.")
    if unplaced:
        notes.append(
            f"Up to {MAX_PLACEMENTS} alignments, strongest first per candidate, "
            "were placed against annotated start "
            f"sites; {unplaced} weaker ones were not, so they carry no distance. They are "
            f"still counted and still listed. Annotation has a {MAX_ANNOTATION_SECONDS:g}s "
            "soft time budget, checked between lookups, as well as a count limit."
        )
    return ScreeningRecord(
        run_id=run_id,
        scope=scope,
        parameters=parameters,
        hits=hits,
        proximities=proximities,
        request_id=submitted.request_id,
        candidates_searched=tuple(sequences),
        notes=tuple(notes),
        on_target_hits=own,
        capped_candidates=tuple(sorted(report.capped)),
        window=window if scope == GENOMIC_SCOPE else None,
    )


def specificity_table(
    record: ScreeningRecord, core_length: int
) -> list[dict[str, Any]]:
    """Return one row per screened candidate, each measure separately sortable.

    Candidates with no hits are present with zeroes rather than absent: a missing row would read
    as "not screened", which is a different statement from "screened, nothing found".
    """
    capped = set(record.capped_candidates)
    rows = []
    for candidate_id in record.candidates_searched:
        measured = measure_specificity(
            candidate_id,
            record.hits_for(candidate_id),
            core_length,
            capped=candidate_id in capped,
            on_target_hits=len(record.on_target_for(candidate_id)),
        )
        row = measured.to_dict()
        row["in_run_window_hits"] = sum(
            1
            for p in record.proximities
            if p.hit.candidate_id == candidate_id and p.inside_run_window
        )
        hits = record.hits_for(candidate_id)
        placements = [p for p in record.proximities if p.hit.candidate_id == candidate_id]
        complete = len(placements) == len(hits) and all(
            p.annotation_status == "available" and p.on_chromosome for p in placements
        )
        row["annotation_complete"] = complete if record.scope == GENOMIC_SCOPE else None
        row["annotation_unknown_hits"] = len(hits) - sum(
            p.annotation_status == "available" and p.on_chromosome for p in placements
        ) if record.scope == GENOMIC_SCOPE else None
        if record.scope == GENOMIC_SCOPE and not complete:
            row["in_run_window_hits"] = None
        names = [re.search(r"\(([A-Za-z0-9_.-]{1,80})\),", hit.title) for hit in hits]
        row["named_genes"] = len({match.group(1) for match in names if match})
        row["unnamed_alignments"] = sum(match is None for match in names)
        row["distinct_sequences"] = len({h.accession.split(".")[0] for h in hits})
        row["distinct_loci"] = (len({(h.accession, min(h.subject_start, h.subject_end),
                                      max(h.subject_start, h.subject_end), h.subject_strand)
                                     for h in hits}) if record.scope == GENOMIC_SCOPE else None)
        rows.append(row)
    return rows
