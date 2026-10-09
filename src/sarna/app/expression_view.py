"""Turn a GTEx profile into the records the transcript picker shows.

The picker's order is annotation's job and stays that way (ADR 0012): this module never sorts
transcripts by expression, because "most expressed in bulk cortex" is not the same claim as
"the reference transcript" and quietly promoting one over the other would make a scientific
choice on the researcher's behalf. What it does is put the evidence next to each option.

Two wordings here are load-bearing. A transcript measured at 0.00 everywhere reads as
*measured at zero*, never as "no data"; and GTEx's nervous-system group is described as what it
is (13 sample groups, one of which is cervical spinal cord rather than brain) instead of
being summed into a single "brain" number that implies a tissue nobody sampled.
"""

from typing import Any

from sarna.expression.gtex import CELL_LINE_CAVEAT
from sarna.expression.models import (
    ExpressionStatus,
    GeneExpression,
    TissueMedian,
    TranscriptExpression,
)
from sarna.expression.organs import BY_ID, organ_for

BRAIN_PREFIX = "Brain_"
"""GTEx's own prefix for its nervous-system groups."""

BRAIN_LABEL = "Brain (GTEx's 13 nervous-system groups)"

BRAIN_NOTE = (
    "GTEx files 13 sample groups under Brain, one of which is cervical spinal cord rather "
    "than brain. The figure shown is the highest single group, not an average over them: "
    "averaging regions nobody pooled would invent a measurement."
)

ORGAN_NOTE = (
    "An organ symbol means GTEx measured a median above zero in at least one of that organ's "
    "sample groups, and it is shown with the group that carried the value. Its absence means "
    "the median was zero there, not that the organ was not looked at. Organs without an "
    "honest symbol are named in words rather than given a lookalike one."
)

CELL_LINE_NOTE = (
    "Not assayed here means exactly that: it is a statement about GTEx and not about those "
    "cells. Nothing on this panel supports or contradicts expression in AF22 or HeLa."
)
"""The interpretation, shown beside
:data:`sarna.expression.gtex.CELL_LINE_CAVEAT`, which states the fact. Kept apart so the two
read as one thought instead of saying the same thing twice, which is how the panel first
rendered."""

TOP_TISSUES = 5
"""How many sample groups each transcript lists before the full table is opened."""

TOP_ORGANS = 4
"""How many organs a transcript's summary shows before it says how many more there are."""


def organ_summary(entry: TranscriptExpression) -> list[dict[str, Any]]:
    """Return the organs this transcript was measured in, strongest first.

    An organ is included only where at least one of its sample groups has a median **above
    zero**. A symbol beside an organ therefore means "measured here", and its absence means
    "not measured above zero here", which is a statement about the median and not about the
    organ. Organs appear with the group that carried the value, so a reader can see that
    "brain" came from frontal cortex rather than from an average nobody computed.
    """
    if not entry.status.is_measured:
        return []
    best: dict[str, TissueMedian] = {}
    for median in entry.medians:
        organ = organ_for(median.tissue_id)
        if organ is None or median.median <= 0:
            continue
        current = best.get(organ.id)
        if current is None or median.median > current.median:
            best[organ.id] = median
    ordered = sorted(best.items(), key=lambda pair: -pair[1].median)
    return [
        {
            "organ_id": organ_id,
            "name": BY_ID[organ_id].name,
            "chip": BY_ID[organ_id].chip,
            "symbol": BY_ID[organ_id].symbol,
            "is_body_site": BY_ID[organ_id].is_body_site,
            "median": median.median,
            "unit": median.unit,
            "tissue_name": median.tissue_name,
            "sample_count": median.sample_count,
        }
        for organ_id, median in ordered
    ]


def _summarise(entry: TranscriptExpression) -> dict[str, Any]:
    """Describe one transcript's expression in the words the picker shows.

    Returns a record whose ``label`` is safe to print beside an option without further
    qualification, including when nothing was measured.
    """
    if entry.status is not ExpressionStatus.MEASURED:
        return {
            "transcript_id": entry.transcript_id,
            "status": str(entry.status),
            "status_note": entry.status.note,
            "silent": False,
            "label": "not in this release, so nothing was quantified",
            "short_label": "not in GTEx",
            "measured_groups": 0,
            "highest": None,
            "brain": None,
            "brain_label": "not in this release, so nothing was quantified",
            "top": [],
            "organs": [],
            "all_organs": [],
            "other_organs": 0,
        }

    ranked = sorted(entry.medians, key=lambda median: -median.median)
    organs = organ_summary(entry)
    brain = [median for median in entry.medians if median.tissue_id.startswith(BRAIN_PREFIX)]
    best_brain = max(brain, key=lambda median: median.median) if brain else None
    highest = ranked[0] if ranked else None
    unit = highest.unit if highest else ""
    silent = highest is not None and highest.median == 0.0
    return {
        "transcript_id": entry.transcript_id,
        "release_id": entry.release_id,
        "status": str(entry.status),
        "status_note": entry.status.note,
        "silent": silent,
        "label": (
            f"measured at 0.00 {unit} in all {len(entry.medians)} sample groups" if silent
            else f"up to {highest.median:.2f} {unit} in {highest.tissue_name}"
            if highest else "no sample groups reported"
        ),
        "short_label": (
            f"brain 0.00 {unit}" if silent
            else f"brain {best_brain.median:.2f} {best_brain.unit}"
            if best_brain and best_brain.median > 0
            else f"brain 0.00 {unit}" if best_brain
            else "no brain group"
        ),
        "measured_groups": len(entry.medians),
        "highest": highest.to_dict() if highest else None,
        "brain": best_brain.to_dict() if best_brain else None,
        "brain_label": _brain_label(best_brain, len(brain)),
        # Nothing to rank when every group is zero: the first five would be an alphabetical
        # slice of zeros, and printing them reads as though those tissues were notable.
        "top": [] if silent else [median.to_dict() for median in ranked[:TOP_TISSUES]],
        "all": [median.to_dict() for median in ranked],
        "organs": organs[:TOP_ORGANS],
        "all_organs": organs,
        "other_organs": max(0, len(organs) - TOP_ORGANS),
    }


def _brain_label(best: TissueMedian | None, groups: int) -> str:
    """Describe the nervous-system groups without naming a region arbitrarily.

    When every group is zero the highest one is whichever happened to sort first, and printing
    its name reads as though that region were distinguished. It is not, so the count is given
    instead.
    """
    if best is None:
        return "no nervous-system group reported for this transcript"
    if best.median == 0.0:
        return f"0.00 {best.unit} in every one of the {groups}"
    return f"{best.median:.2f} {best.unit} in {best.tissue_name}"


def _caveat_list(caveats: tuple[str, ...]) -> list[str]:
    """Return the caveats in reading order, with the interpretation after what it interprets.

    Rendered as two separate fields, the panel printed "not assayed here means exactly that"
    above the line naming the cells GTEx does not hold, which is backwards. Matched on the
    caveat itself rather than on its position, so a source that does not carry it still gets
    the note, at the end.
    """
    ordered: list[str] = []
    for caveat in caveats:
        ordered.append(caveat)
        if caveat == CELL_LINE_CAVEAT:
            ordered.append(CELL_LINE_NOTE)
    if CELL_LINE_CAVEAT not in caveats:
        ordered.append(CELL_LINE_NOTE)
    return ordered


def expression_overview(profile: GeneExpression) -> dict[str, Any]:
    """Freeze a whole gene's expression evidence beside the picker's own ordering.

    ``transcripts`` is keyed by unversioned transcript identifier so the page can annotate an
    option it already has, rather than rebuilding its list from this payload, which is what
    keeps expression out of the ordering decision.
    """
    summaries = {
        entry.transcript_id: _summarise(entry) for entry in profile.transcripts
    }
    measured = profile.measured_transcripts
    expressed = [
        entry.transcript_id for entry in measured
        if any(median.median > 0 for median in entry.medians)
    ]
    return {
        "source": "GTEx Portal",
        "dataset_id": profile.dataset_id,
        "annotation_release": profile.annotation_release,
        "release_gene_id": profile.release_gene_id,
        "provenance": str(profile.provenance),
        "caveats": list(profile.caveats),
        "caveat_list": _caveat_list(profile.caveats),
        "cell_line_note": CELL_LINE_NOTE,
        "brain_label": BRAIN_LABEL,
        "brain_note": BRAIN_NOTE,
        "organ_note": ORGAN_NOTE,
        "measured_transcripts": len(measured),
        "silent_transcripts": len(measured) - len(expressed),
        "transcripts": summaries,
    }


def unavailable(reason: str, code: str) -> dict[str, Any]:
    """Describe a profile that could not be retrieved, without implying absence of expression.

    A failed lookup and a measured zero must not render alike, so this carries no transcript
    records at all rather than a set of empty ones.
    """
    return {
        "source": "GTEx Portal",
        "available": False,
        "code": code,
        "detail": reason,
        "note": "No expression evidence was retrieved. This says nothing about whether these "
                "transcripts are expressed.",
        "transcripts": {},
    }


__all__ = [
    "BRAIN_LABEL",
    "BRAIN_NOTE",
    "CELL_LINE_NOTE",
    "TOP_TISSUES",
    "expression_overview",
    "unavailable",
]
