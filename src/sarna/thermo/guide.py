"""Which strand of the duplex is predicted to become the guide, and why.

RISC loads whichever strand has the more weakly paired 5' end (Schwarz et al. 2003; Khvorova et
al. 2003). The antisense strand's 5' end pairs with the **3' end of the target**, so comparing the
target's two ends answers the question:

- target 3' end less stable than its 5' end → the antisense 5' end is the weaker → **antisense**
  is loaded. This is the desirable case for a promoter-targeting saRNA, and the legacy script's
  positive ``stabilitydiff``.
- target 5' end less stable → the **sense** strand would be loaded instead, diluting activity and
  risking off-target effects from the wrong strand.
- the two within tolerance → **ambiguous**. Not a strand picked by coin flip.

The sign convention is the legacy script's, preserved: ``difference = E(3') - E(5')``, positive
being the desirable condition.

ADR 0009 requires that a prediction never arrive as a bare verdict. Each assessment therefore
carries the table it came from, both end energies with their per-stack breakdown, the difference,
and an explanation in words naming the strand, the mechanism and its sources. Where several tables
are assessed, :func:`compare_asymmetry` reports all of them and says whether they agree; it
deliberately offers **no** single combined prediction, because the parameter sets are not
independent measurements of one quantity and averaging them would manufacture a number with no
source.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Any

from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.ends import (
    DEFAULT_END_LENGTH,
    DuplexEnd,
    EndStability,
    end_stability,
)
from sarna.thermo.model import NearestNeighbourTable

LEGACY_TOLERANCE = 1e-8
"""The legacy script's ``fudge`` constant, retained so boundary behaviour is reproduced."""

_MECHANISM = (
    "RISC loads whichever strand has the more weakly paired 5' end (Schwarz et al. 2003; "
    "Khvorova et al. 2003)"
)


class GuideStrand(Enum):
    """Which strand of the duplex is predicted to be loaded as the guide."""

    ANTISENSE = "antisense"
    """Complementary to the target site. The desirable outcome for a promoter-targeting saRNA."""

    SENSE = "sense"
    """Same sequence as the target site. Loading this strand is the undesirable outcome."""

    AMBIGUOUS = "ambiguous"
    """The two ends differ by less than the tolerance, so neither strand is favoured."""

    def __str__(self) -> str:
        """Return the strand label."""
        return self.value


@dataclass(frozen=True, slots=True)
class AsymmetryAssessment:
    """One nearest-neighbour table's verdict on which strand becomes the guide."""

    table_id: str
    table_name: str
    five_prime: EndStability
    three_prime: EndStability
    difference: float
    """``E(3') - E(5')``. Positive means the target's 5' end is the more stable."""

    tolerance: float
    """Below this magnitude the difference is treated as no difference."""

    predicted_guide: GuideStrand

    @property
    def unit(self) -> str:
        """Return the energy unit the difference is expressed in."""
        return self.five_prime.unit

    @property
    def explanation(self) -> str:
        """Return the reasoning in full, for display beside the verdict.

        Names the measured energies, the strand, and the mechanism the inference rests on, so a
        researcher can answer "why does it think that?" without reading the source (ADR 0009).
        """
        measured = (
            f"Using {self.table_name}, the target's 5' end scores "
            f"{self.five_prime.energy:.2f} {self.unit} "
            f"({self.five_prime.breakdown}) and its 3' end "
            f"{self.three_prime.energy:.2f} {self.unit} "
            f"({self.three_prime.breakdown}), a difference of "
            f"{self.difference:+.2f} {self.unit}."
        )
        if self.predicted_guide is GuideStrand.AMBIGUOUS:
            return (
                f"{measured} The two ends are equally stable to within "
                f"{self.tolerance:g} {self.unit}, so neither strand is favoured for loading and "
                f"no guide is predicted. {_MECHANISM}, and here that difference is absent."
            )
        if self.predicted_guide is GuideStrand.ANTISENSE:
            return (
                f"{measured} The 3' end is the less stable of the two, and the antisense "
                f"strand's 5' end pairs with it, so the antisense 5' end is the weaker. "
                f"{_MECHANISM}, so the antisense strand is predicted to be the guide."
            )
        return (
            f"{measured} The 5' end is the less stable of the two, which is the 5' end of the "
            f"sense strand itself. {_MECHANISM}, so the sense strand is predicted to be the "
            f"guide, the undesirable outcome, since the antisense strand is the one intended to "
            f"engage the promoter."
        )

    def to_dict(self) -> dict[str, Any]:
        """Export the verdict together with the two measurements it was derived from.

        ADR 0009 requires that a prediction never arrive as a bare verdict. A reader who
        receives only ``difference`` has to trust the sign; one who receives both ends and
        their stacks can redo the subtraction. :attr:`explanation` says the same thing in
        words, which a person can read but an interface cannot draw.
        """
        return {
            "table": self.table_id,
            "table_name": self.table_name,
            "predicted_guide": str(self.predicted_guide),
            "difference": self.difference,
            "tolerance": self.tolerance,
            "unit": self.unit,
            "five_prime": self.five_prime.to_dict(),
            "three_prime": self.three_prime.to_dict(),
            "explanation": self.explanation,
        }


def assess_asymmetry(
    target: NucleicSequence,
    table: NearestNeighbourTable,
    end_length: int = DEFAULT_END_LENGTH,
    tolerance: float = LEGACY_TOLERANCE,
) -> AsymmetryAssessment:
    """Assess which strand becomes the guide, according to one parameter set.

    Args:
        target: The target site, normally DNA.
        table: The stacking parameters to use.
        end_length: Nucleotides defining each end.
        tolerance: Differences smaller than this are reported as ambiguous.

    Returns:
        The assessment, including its explanation.

    Raises:
        ValueError: If the tolerance is negative, or the end length does not fit the target.
    """
    if tolerance < 0:
        raise ValueError(f"tolerance must not be negative, got {tolerance}")
    five_prime = end_stability(target, table, DuplexEnd.FIVE_PRIME, end_length)
    three_prime = end_stability(target, table, DuplexEnd.THREE_PRIME, end_length)
    difference = three_prime.energy - five_prime.energy
    if abs(difference) <= tolerance:
        predicted = GuideStrand.AMBIGUOUS
    elif difference > 0:
        predicted = GuideStrand.ANTISENSE
    else:
        predicted = GuideStrand.SENSE
    return AsymmetryAssessment(
        table_id=table.id,
        table_name=table.name,
        five_prime=five_prime,
        three_prime=three_prime,
        difference=difference,
        tolerance=tolerance,
        predicted_guide=predicted,
    )


@dataclass(frozen=True, slots=True)
class AsymmetryComparison:
    """Several tables' verdicts on one candidate, side by side.

    There is deliberately no combined prediction. Where the tables disagree, that disagreement is
    the finding: it means the asymmetry is small enough for the parameter set to decide the
    answer, which is itself a reason to trust the prediction less (ADR 0009).
    """

    assessments: tuple[AsymmetryAssessment, ...]

    @property
    def predictions(self) -> frozenset[GuideStrand]:
        """Return the distinct guide strands predicted across the tables."""
        return frozenset(assessment.predicted_guide for assessment in self.assessments)

    @property
    def tables_agree(self) -> bool:
        """Return whether every table predicts the same guide strand."""
        return len(self.predictions) == 1

    @property
    def summary(self) -> str:
        """Return a one-line statement of agreement or disagreement, naming the tables."""
        if self.tables_agree:
            agreed = next(iter(self.predictions))
            tables = len(self.assessments)
            if tables == 1:
                only = self.assessments[0]
                return f"{only.table_name} predicts the {agreed} strand as the guide."
            return (
                f"All {tables} parameter sets agree: the {agreed} strand is predicted to be "
                f"the guide."
            )
        verdicts = "; ".join(
            f"{assessment.table_name} predicts {assessment.predicted_guide}"
            for assessment in self.assessments
        )
        return (
            f"The parameter sets disagree about which strand becomes the guide, {verdicts}. "
            f"The asymmetry is small enough that the choice of parameters decides the answer, "
            f"so treat this candidate's strand selection as uncertain."
        )


def compare_asymmetry(
    target: NucleicSequence,
    tables: tuple[NearestNeighbourTable, ...],
    end_length: int = DEFAULT_END_LENGTH,
    tolerance: float = LEGACY_TOLERANCE,
) -> AsymmetryComparison:
    """Assess ``target`` against several parameter sets and report them together.

    Args:
        target: The target site, normally DNA.
        tables: The parameter sets to assess, in the order they should be displayed.
        end_length: Nucleotides defining each end.
        tolerance: Differences smaller than this are reported as ambiguous.

    Returns:
        Every table's assessment, with whether they agree.

    Raises:
        ValueError: If no tables are given, or one appears more than once.
    """
    if not tables:
        raise ValueError("at least one nearest-neighbour table is required")
    identifiers = [table.id for table in tables]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(
            f"a table appears more than once in {identifiers}; each parameter set yields one "
            f"assessment, so a repeat would double-report the same verdict"
        )
    return AsymmetryComparison(
        assessments=tuple(
            assess_asymmetry(target, table, end_length, tolerance) for table in tables
        )
    )
