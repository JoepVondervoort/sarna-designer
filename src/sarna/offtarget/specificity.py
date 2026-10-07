"""Turning hits into the several numbers a researcher can sort on.

There is deliberately **no specificity score**. A single number would have to weigh "one perfect
match in another transcript" against "forty weak partial matches", and those are different
experimental problems: the first can silence a specific unintended gene, the second is mostly
noise. Any weighting ranks them identically for some choice of weights, and the choice would be
ours rather than the researcher's. ADR 0006 already refuses weighted sums for ranking candidates;
ADR 0019 applies the same reasoning here.

So each measure is reported and sortable on its own, and the interface offers them as separate
columns.

**Nothing here filters, and a clean result is not proof of specificity.** BLAST finding nothing
is one search of one database on one day, with parameters that can miss things. The summary text
says so rather than congratulating the candidate.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from sarna.offtarget.models import OffTargetHit

NEAR_PERFECT_IDENTITIES = 17
"""Matching bases at which an alignment counts as a near-perfect match of a 19 nt core.

Two mismatches. Chosen as a reporting threshold, not a biological one: it is loose enough to
surface the matches a researcher would want to look at and tight enough that the count stays
meaningful. It changes what is *counted*, never what is kept, and the raw hits travel beside it
so a different threshold can be applied by eye.
"""

NO_EVALUE_SORTS_LAST = float("inf")
"""Sort key for a candidate with no alignments at all.

``None`` cannot be compared with a float, and mapping it to 0.0 would sort a candidate with no
hits as though it had the strongest possible one.
"""


@dataclass(frozen=True, slots=True)
class CandidateSpecificity:
    """Several independent measures of what else one candidate could bind."""

    candidate_id: str
    hit_count: int
    best_evalue: float | None
    near_perfect_hits: int
    perfect_hits: int
    longest_complementary_stretch: int

    @property
    def sort_key_best_evalue(self) -> float:
        """Return a sortable best E-value, with "no hits at all" sorting last."""
        return NO_EVALUE_SORTS_LAST if self.best_evalue is None else self.best_evalue

    def to_dict(self) -> dict[str, Any]:
        """Export every measure separately. There is no combined figure to export."""
        return {
            "candidate": self.candidate_id,
            "hit_count": self.hit_count,
            "best_evalue": self.best_evalue,
            "near_perfect_hits": self.near_perfect_hits,
            "perfect_hits": self.perfect_hits,
            "longest_complementary_stretch": self.longest_complementary_stretch,
            "summary": summarise(self),
        }


def measure_specificity(
    candidate_id: str, hits: Sequence[OffTargetHit], core_length: int
) -> CandidateSpecificity:
    """Measure what ``hits`` say about one candidate.

    Args:
        candidate_id: The candidate these hits belong to.
        hits: Its alignments, which may be empty.
        core_length: The duplex core length, against which a match is judged complete.

    Returns:
        The separate measures, with no combined score.

    Raises:
        ValueError: If a hit belongs to a different candidate, or ``core_length`` is not
            positive. Measuring another candidate's hits would attribute a risk to the wrong
            duplex, silently.
    """
    if core_length <= 0:
        raise ValueError(f"core length must be positive, got {core_length}")
    wrong = {h.candidate_id for h in hits} - {candidate_id}
    if wrong:
        raise ValueError(
            f"hits for {', '.join(sorted(wrong))} were passed while measuring {candidate_id!r}; "
            f"each candidate's specificity must be measured from its own alignments"
        )
    return CandidateSpecificity(
        candidate_id=candidate_id,
        hit_count=len(hits),
        best_evalue=min((h.evalue for h in hits), default=None),
        near_perfect_hits=sum(1 for h in hits if h.identities >= NEAR_PERFECT_IDENTITIES),
        perfect_hits=sum(1 for h in hits if h.is_perfect(core_length)),
        longest_complementary_stretch=max((h.identities for h in hits), default=0),
    )


def summarise(measured: CandidateSpecificity) -> str:
    """Return one sentence a researcher can read without opening the alignment table.

    States findings. Never recommends discarding a candidate: that is the researcher's judgement,
    and accepted constraint 15 makes screening report rather than exclude.
    """
    if measured.hit_count == 0:
        return (
            "No alignments were returned for this candidate. That does not establish "
            "specificity: it is one search, of one database, on one day, with parameters that "
            "can miss short or partial matches."
        )
    parts = [f"{measured.hit_count} alignment{'s' if measured.hit_count != 1 else ''} returned"]
    if measured.perfect_hits:
        parts.append(
            f"{measured.perfect_hits} of them a perfect full-length match elsewhere"
        )
    elif measured.near_perfect_hits:
        parts.append(
            f"{measured.near_perfect_hits} matching at least "
            f"{NEAR_PERFECT_IDENTITIES} bases"
        )
    else:
        parts.append(
            f"the closest matching {measured.longest_complementary_stretch} bases"
        )
    if measured.best_evalue is not None:
        parts.append(f"best E-value {measured.best_evalue:.2g}")
    return (
        ", ".join(parts)
        + ". Sequence similarity is not measured binding, and these findings do not change "
        "the candidate's rank."
    )
