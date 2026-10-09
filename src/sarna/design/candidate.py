"""One designed candidate, with every verdict that was reached about it.

A candidate is the unit a researcher reads, so it carries the whole evaluation rather than a
summary of it: every rule's outcome, reason and measurements, the per-table guide predictions, and
the arithmetic behind its rank. Nothing here is computed again later for display, if the
interface shows it, this object already said it.

The design follows from two commitments made earlier:

**A rejected candidate is still a candidate.** It keeps its id, its position, its sequences and
its full rule results, and names the rules that excluded it. A sequence that disappears from the
output without explanation is the behaviour this project exists to replace.

**The rank is an itemised checklist, not a number.** :attr:`Candidate.preference_points` out of
:attr:`Candidate.max_preference_points`, with the met and unmet preferences listed (ADR 0006).
Preferences that could not be evaluated are excluded from the *evaluated* denominator, because a
preference nobody could check is not one the candidate failed to meet.
"""

from dataclasses import dataclass
from typing import Any

from sarna.genome.coordinates import GenomicInterval
from sarna.genome.tss import TssPosition
from sarna.rules.results import RuleOutcome, RuleResult
from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.guide import AsymmetryComparison


@dataclass(frozen=True, slots=True)
class Candidate:
    """A generated candidate and everything concluded about it."""

    id: str
    """Stable identifier, numbered by position in the promoter rather than by rank."""

    offset: int
    """0-based offset into the promoter sequence."""

    target: NucleicSequence
    """The target site, as DNA. A genomic locus, not a delivered molecule."""

    sense_strand: NucleicSequence
    """The duplex strand matching the target site, as RNA."""

    antisense_strand: NucleicSequence
    """The duplex strand complementary to the target site, as RNA. The intended guide."""

    tss_span: tuple[TssPosition, TssPosition] | None
    interval: GenomicInterval | None
    results: tuple[RuleResult, ...]
    guide_comparison: AsymmetryComparison
    max_preference_points: int
    rank: int | None = None
    """Position in the ranked list, 1-based. ``None`` for a rejected candidate."""

    def with_rank(self, rank: int) -> "Candidate":
        """Return a copy carrying ``rank``.

        A copy rather than a mutation, so ranking cannot retroactively alter the candidates a
        caller already holds.
        """
        return Candidate(
            id=self.id,
            offset=self.offset,
            target=self.target,
            sense_strand=self.sense_strand,
            antisense_strand=self.antisense_strand,
            tss_span=self.tss_span,
            interval=self.interval,
            results=self.results,
            guide_comparison=self.guide_comparison,
            max_preference_points=self.max_preference_points,
            rank=rank,
        )

    def result(self, rule_id: str) -> RuleResult:
        """Return this candidate's result for one rule.

        Raises:
            KeyError: If the rule was not part of the configuration that produced it.
        """
        for result in self.results:
            if result.rule_id == rule_id:
                return result
        raise KeyError(
            f"candidate {self.id} has no result for rule {rule_id!r}; it was not selected by "
            f"the configuration this candidate was designed with"
        )

    @property
    def excluded_by(self) -> tuple[RuleResult, ...]:
        """Return the mandatory rules this candidate failed. Empty if it passed.

        Every failing gate is reported, not just the first, so a researcher sees whether a
        candidate missed by one criterion or by five.
        """
        return tuple(result for result in self.results if result.excludes)

    @property
    def is_rejected(self) -> bool:
        """Return whether a mandatory rule excluded this candidate."""
        return bool(self.excluded_by)

    @property
    def rejection_summary(self) -> str:
        """Return why this candidate was excluded, in plain language. Empty if it passed."""
        if not self.is_rejected:
            return ""
        reasons = " ".join(result.detail for result in self.excluded_by)
        names = ", ".join(result.rule_id for result in self.excluded_by)
        return f"Excluded by {names}. {reasons}"

    @property
    def preference_points(self) -> int:
        """Return how many preferences this candidate meets."""
        return sum(result.points for result in self.results)

    @property
    def met_preferences(self) -> tuple[RuleResult, ...]:
        """Return the preferences this candidate meets."""
        return tuple(result for result in self.results if result.points)

    @property
    def unmet_preferences(self) -> tuple[RuleResult, ...]:
        """Return the preferences evaluated and not met. Excludes unevaluable ones."""
        return tuple(
            result
            for result in self.results
            if result.kind.value == "preference" and result.outcome is RuleOutcome.FAIL
        )

    @property
    def unevaluated(self) -> tuple[RuleResult, ...]:
        """Return the rules that could not be evaluated for this candidate."""
        return tuple(result for result in self.results if result.is_unevaluated)

    @property
    def evaluated_preference_points(self) -> int:
        """Return how many preferences could be evaluated at all.

        The honest denominator where some preference could not be checked, a candidate at the
        promoter's 3' edge has no flank to read, and scoring those positions zero would silently
        penalise the last candidates of every promoter.
        """
        return len(self.met_preferences) + len(self.unmet_preferences)

    @property
    def tiebreak_value(self) -> float | None:
        """Return the first tie-break's ordering value, or ``None`` if none was computed."""
        for result in self.results:
            if result.ordering_value is not None:
                return result.ordering_value
        return None

    @property
    def score_summary(self) -> str:
        """Return the rank as an itemised statement, e.g. ``meets 5 of 7 preferences``."""
        summary = f"Meets {self.preference_points} of {self.max_preference_points} preferences"
        if self.evaluated_preference_points != self.max_preference_points:
            unchecked = self.max_preference_points - self.evaluated_preference_points
            summary += f" ({unchecked} could not be evaluated)"
        return summary

    def to_dict(self) -> dict[str, Any]:
        """Return the candidate as JSON-serialisable data, including every rule verdict."""
        return {
            "id": self.id,
            "rank": self.rank,
            "offset": self.offset,
            "target_site_dna": self.target.bases,
            "sense_strand_rna": self.sense_strand.bases,
            "antisense_strand_rna": self.antisense_strand.bases,
            "gc_percent": round(self.target.gc_percent, 1),
            "tss_span": (
                [self.tss_span[0].value, self.tss_span[1].value] if self.tss_span else None
            ),
            "genomic": (
                {
                    "assembly": str(self.interval.assembly),
                    "contig": self.interval.contig,
                    "start": self.interval.start,
                    "end": self.interval.end,
                    "strand": str(self.interval.strand),
                }
                if self.interval
                else None
            ),
            "preference_points": self.preference_points,
            "max_preference_points": self.max_preference_points,
            "score_summary": self.score_summary,
            "excluded_by": [result.rule_id for result in self.excluded_by],
            "rejection_summary": self.rejection_summary,
            "guide_prediction": {
                "tables_agree": self.guide_comparison.tables_agree,
                "summary": self.guide_comparison.summary,
                "assessments": [
                    assessment.to_dict() for assessment in self.guide_comparison.assessments
                ],
            },
            "rules": [
                {
                    "id": result.rule_id,
                    "kind": str(result.kind),
                    "outcome": str(result.outcome),
                    "detail": result.detail,
                    "measurements": dict(result.measurements),
                    "points": result.points,
                }
                for result in self.results
            ],
        }
