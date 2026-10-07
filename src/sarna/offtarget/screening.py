"""Running a screening end to end: submit, wait, place, measure.

This is the piece that turns a finished design run into a dated screening record. It is written
to be called from a background worker, and to be interruptible in the sense that matters: every
stage it reaches is reported, so the interface can say *what* it is waiting for instead of
spinning.

Sleeping is injected rather than imported so tests do not wait fifteen real seconds per poll.
"""

import time
from collections.abc import Callable, Mapping, Sequence
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
from sarna.offtarget.proximity import TranscriptLookup, place_hit
from sarna.offtarget.specificity import measure_specificity

MAX_CANDIDATES = 50
"""Most candidates to screen in one job.

A courtesy to NCBI rather than a technical limit: a run can generate thousands of candidates,
and only the ones a researcher might order are worth a shared public service's time. Ranked
order, so the cut falls on the least promising.
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
) -> ScreeningRecord:
    """Screen one run's passing candidates and return the dated record.

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
    sequences = candidate_sequences(candidates, limit=MAX_CANDIDATES)
    if not sequences:
        raise ValueError(
            "no candidate passed the gates in this run, so there is nothing to screen"
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
            hits = client.fetch(submitted.request_id, parameters)
            break
        except BlastStillRunning:
            say(f"NCBI is still searching ({int(waited)}s elapsed)")
            if now() > deadline:
                raise BlastError(
                    f"NCBI search {submitted.request_id} had not finished after "
                    f"{int(max_wait)}s. It may still complete; the designs are unaffected."
                ) from None

    say(f"Reading {len(hits)} alignments")
    proximities: tuple[PromoterProximity, ...] = ()
    if scope == GENOMIC_SCOPE and transcripts is not None:
        say("Placing hits against annotated start sites")
        proximities = tuple(place_hit(hit, transcripts, window) for hit in hits)

    notes = []
    if len(candidates) and len(sequences) < sum(1 for c in candidates if c.get("rank")):
        notes.append(
            f"Only the top {len(sequences)} ranked candidates were screened, out of "
            f"{sum(1 for c in candidates if c.get('rank'))} that passed the gates."
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
    )


def specificity_table(
    record: ScreeningRecord, core_length: int
) -> list[dict[str, Any]]:
    """Return one row per screened candidate, each measure separately sortable.

    Candidates with no hits are present with zeroes rather than absent: a missing row would read
    as "not screened", which is a different statement from "screened, nothing found".
    """
    rows = []
    for candidate_id in record.candidates_searched:
        measured = measure_specificity(
            candidate_id, record.hits_for(candidate_id), core_length
        )
        row = measured.to_dict()
        row["in_run_window_hits"] = sum(
            1
            for p in record.proximities
            if p.hit.candidate_id == candidate_id and p.inside_run_window
        )
        rows.append(row)
    return rows
