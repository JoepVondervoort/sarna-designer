"""Running a screening end to end: submit, wait, place, measure.

This is the piece that turns a finished design run into a dated screening record. It is written
to be called from a background worker, and to be interruptible in the sense that matters: every
stage it reaches is reported, so the interface can say *what* it is waiting for instead of
spinning.

Sleeping is injected rather than imported so tests do not wait fifteen real seconds per poll.
"""

import time
from collections.abc import Callable, Collection, Mapping, Sequence
from typing import Any

from sarna.offtarget.models import (
    GENOMIC_SCOPE,
    PromoterProximity,
    ScreeningRecord,
)
from sarna.offtarget.ncbi import (
    MINIMUM_POLL_SECONDS,
    BlastClient,
    BlastError,
    BlastStillRunning,
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
"""Give up after half an hour. A search that has not finished by then will not be waited for."""

ProgressReporter = Callable[[str], None]


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
        max_wait: Give up after this long.
        only: Screen exactly these candidates rather than the best ``limit`` of them.
        limit: Most candidates in this one job, never above :data:`MAX_CANDIDATES`.

    Returns:
        The screening, including the parameters and the moment it was made.

    Raises:
        BlastError: If the search could not be run or did not finish in time.
        ValueError: If there is nothing to screen, or the genomic scope has no lookup.
    """
    say = progress or (lambda _stage: None)
    parameters = parameters_for(scope)
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

    say(f"Submitting {len(sequences)} sequences to NCBI")
    submitted = client.submit(sequences, scope)
    say(f"NCBI accepted the search (id {submitted.request_id})")

    deadline = now() + max_wait
    interval = max(poll_seconds, MINIMUM_POLL_SECONDS)
    # NCBI answers every submission with RTOE, its own estimate of how long the search will
    # take, and the search cannot be ready before then. Waiting it out costs no latency and
    # removes every poll that could only have returned WAITING: at the fifteen-second floor a
    # seven-minute search made twenty-eight requests against a shared service ADR 0019 notes
    # can block the whole institution's address. Capped at half the budget so a wild estimate
    # cannot turn an early finish into a wait for the deadline, and floored at the rate limit.
    wait_next = min(max(float(submitted.estimated_seconds), interval), max_wait / 2)
    waited = 0.0
    while True:
        sleep(wait_next)
        waited += wait_next
        # Only the first wait follows the estimate. Afterwards return to the floor, so a job
        # that overruns is still collected within one poll of finishing rather than backing off.
        wait_next = interval
        try:
            report = client.fetch(submitted.request_id, parameters)
            break
        except BlastStillRunning:
            say(f"NCBI is still searching ({int(waited)}s elapsed)")
            if now() > deadline:
                raise BlastError(
                    f"NCBI search {submitted.request_id} had not finished after "
                    f"{int(max_wait)}s. It may still complete; the designs are unaffected."
                ) from None

    say(f"Reading {len(report.hits)} alignments")
    # The candidate's own promoter is a perfect match to itself, so without this every genomic
    # candidate reports at least one perfect hit and the column means nothing (ADR 0023).
    # Only the genomic scope: a promoter is usually not transcribed, and the only handle in a
    # transcript report would be a gene name, which is a guess.
    loci = loci_from_candidates(candidates) if scope == GENOMIC_SCOPE else {}
    hits, own = split_on_target(report.hits, loci)

    proximities: tuple[PromoterProximity, ...] = ()
    if scope == GENOMIC_SCOPE and transcripts is not None:
        say(f"Placing {len(hits)} hits against annotated start sites")
        proximities = tuple(place_hit(hit, transcripts, window) for hit in hits)

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
            f"{'s' if len(report.capped) != 1 else ''} returned more alignments than were kept "
            f"({parameters.max_alignments_per_candidate} each), so those counts are floors."
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
        rows.append(row)
    return rows
