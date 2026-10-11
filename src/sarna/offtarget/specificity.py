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

import math
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

**Locally chosen, with no source.** Kept because saved runs and the xlsx export already carry
this column and silently changing what a stored number means is worse than carrying two, but
:data:`WANG_HOMOLOGY_FRACTION` is what the interface leads with (ADR 0023).
"""

WANG_HOMOLOGY_FRACTION = 0.80
"""The one homology threshold Wang 2015 actually states.

The paper rejects a **control** duplex showing "significant homology (>= 80%)" to anything.
Applying it to candidate off-targets is an extension by analogy and ADR 0023 records that cost;
what it buys is that the number a researcher reads first has a citation behind it, which
:data:`NEAR_PERFECT_IDENTITIES` never had. Nothing is rejected here either way: ADR 0019 makes
screening report rather than filter.

Source: Wang J, Place RF, Portnoy V, Huang V, Kang MR, Kosaka M, Ho MKC, Li L-C. *Inducing gene
expression by targeting promoter sequences using small activating RNAs.* J Biol Methods, 2015
(PMC4379447).
"""

PROMISCUOUS_STRONG_HITS = 5
"""Alignments at Wang's threshold after which a candidate is summarised rather than listed.

A **display** rule and nothing else. It changes no rank, excludes no candidate and alters no
stored measure; the alignments stay in the record and the list expands on demand. The number is
the user's and no literature supports it, which is exactly why it is declared here rather than
buried in the interface (ADR 0023).
"""


def wang_homology_identities(core_length: int) -> int:
    """Return the matching bases that reach Wang's 80% for a core of this length.

    Rounded up, so the threshold is the smallest count that genuinely reaches 80%: at 19 nt
    that is 16, because 15 of 19 is 78.9% and would fall short of the published figure.

    Raises:
        ValueError: If the core length is not positive.
    """
    if core_length <= 0:
        raise ValueError(f"core length must be positive, got {core_length}")
    return math.ceil(core_length * WANG_HOMOLOGY_FRACTION)

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
    wang_homology_hits: int = 0
    """Alignments reaching Wang 2015's 80%, the one threshold the paper states."""

    strong_threshold: int = 0
    """The matching bases ``wang_homology_hits`` was counted at, so the column can label itself."""

    alignments_capped: bool = False
    """Whether this candidate's alignment list was cut at the per-candidate limit.

    Three outcomes again: no alignments, some alignments, and more than were kept. A capped list
    read as a complete one would understate exactly the candidates that deserve attention.
    """

    on_target_hits: int = 0
    """Alignments removed as this candidate's own locus (ADR 0023).

    Reported so a zero that had something taken out of it cannot be read as a zero that never
    did. The hits themselves stay on the record.
    """

    @property
    def sort_key_best_evalue(self) -> float:
        """Return a sortable best E-value, with "no hits at all" sorting last."""
        return NO_EVALUE_SORTS_LAST if self.best_evalue is None else self.best_evalue

    @property
    def is_promiscuous(self) -> bool:
        """Whether this candidate reached the count at which the list is summarised, not drawn.

        A display state. Nothing is excluded, reordered or rescored by it.
        """
        return self.wang_homology_hits >= PROMISCUOUS_STRONG_HITS

    def to_dict(self) -> dict[str, Any]:
        """Export every measure separately. There is no combined figure to export."""
        return {
            "candidate": self.candidate_id,
            "hit_count": self.hit_count,
            "best_evalue": self.best_evalue,
            "near_perfect_hits": self.near_perfect_hits,
            "perfect_hits": self.perfect_hits,
            "longest_complementary_stretch": self.longest_complementary_stretch,
            "wang_homology_hits": self.wang_homology_hits,
            "strong_threshold": self.strong_threshold,
            "alignments_capped": self.alignments_capped,
            "on_target_hits": self.on_target_hits,
            "promiscuous": self.is_promiscuous,
            "summary": summarise(self),
        }


def measure_specificity(
    candidate_id: str,
    hits: Sequence[OffTargetHit],
    core_length: int,
    capped: bool = False,
    on_target_hits: int = 0,
) -> CandidateSpecificity:
    """Measure what ``hits`` say about one candidate.

    Args:
        candidate_id: The candidate these hits belong to.
        hits: Its alignments, which may be empty. On-target hits are expected to have been
            removed already, which is why their count is passed separately.
        core_length: The duplex core length, against which a match is judged complete.
        capped: Whether the alignment list was cut at the per-candidate limit.
        on_target_hits: How many alignments were removed as this candidate's own locus.

    Returns:
        The separate measures, with no combined score.

    Raises:
        ValueError: If a hit belongs to a different candidate, or ``core_length`` is not
            positive. Measuring another candidate's hits would attribute a risk to the wrong
            duplex, silently.
    """
    strong = wang_homology_identities(core_length)
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
        wang_homology_hits=sum(1 for h in hits if h.identities >= strong),
        strong_threshold=strong,
        alignments_capped=capped,
        on_target_hits=on_target_hits,
    )


def summarise(measured: CandidateSpecificity) -> str:
    """Return one sentence a researcher can read without opening the alignment table.

    States findings. Never recommends discarding a candidate: that is the researcher's judgement,
    and accepted constraint 15 makes screening report rather than exclude.
    """
    removed = (
        f" {measured.on_target_hits} alignment"
        f"{'s' if measured.on_target_hits != 1 else ''} at this candidate's own locus "
        f"{'were' if measured.on_target_hits != 1 else 'was'} removed before counting."
        if measured.on_target_hits
        else ""
    )
    if measured.hit_count == 0:
        return (
            "No alignments were returned for this candidate. That does not establish "
            "specificity: it is one search, of one database, on one day, with parameters that "
            "can miss short or partial matches." + removed
        )
    counted = "kept" if measured.alignments_capped else "returned"
    parts = [f"{measured.hit_count} alignment{'s' if measured.hit_count != 1 else ''} {counted}"]
    if measured.is_promiscuous:
        parts.append(
            f"{measured.wang_homology_hits} of them matching at least "
            f"{measured.strong_threshold} bases, so this sequence occurs in many places"
        )
    elif measured.perfect_hits:
        parts.append(
            f"{measured.perfect_hits} of them a perfect full-length match elsewhere"
        )
    elif measured.wang_homology_hits:
        parts.append(
            f"{measured.wang_homology_hits} matching at least "
            f"{measured.strong_threshold} bases, Wang 2015's 80%"
        )
    else:
        parts.append(
            f"the closest matching {measured.longest_complementary_stretch} bases"
        )
    if measured.best_evalue is not None:
        parts.append(f"best E-value {measured.best_evalue:.2g}")
    capped = (
        f" The list was cut at {measured.hit_count} alignments, so there are more than are "
        f"shown and every count here is a floor."
        if measured.alignments_capped
        else ""
    )
    return (
        ", ".join(parts)
        + ". Sequence similarity is not measured binding, and these findings do not change "
        "the candidate's rank."
        + capped
        + removed
    )
