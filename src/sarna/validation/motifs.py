"""Positional motif enrichment, and whether an enrichment means anything.

Chapter 4 reports that functional saRNAs carry G at position 10, U at position 15 and so on,
with differences around 20-25 percentage points against the non-functional set. Those numbers
reproduce exactly from the published sequences. The question this module exists to answer is
the one the chapter leaves open, having set out to be exploratory: how large a difference would
appear if the labels meant nothing at all?

The answer governs how the project is allowed to use the finding, so it is computed here rather
than asserted. Two properties of the data make the question unavoidable:

**The grid is large.** Nineteen positions times four bases is 76 single-nucleotide cells, 288
dinucleotide cells and 1088 trinucleotide cells. Scanning that many cells for the largest
difference and then reporting that difference is a maximum, not a sample: its null distribution
is far to the right of any single cell's.

**One sequence is worth 7.1 points.** With 14 activators, a frequency can only take the values
0, 7.1, 14.3 ... percent. A 21-point difference is three sequences. Nothing finer than a third
of a sequence is measurable, so apparent precision in these numbers is an artefact of the
percentage scale.

The test is a max-statistic permutation test, which handles both at once: shuffle the labels,
recompute every cell, keep the largest absolute difference, and repeat. The fraction of shuffles
whose maximum reaches the observed maximum is a family-wise error rate over the whole grid, with
no independence assumption — appropriate here, because neighbouring positions and overlapping
k-mers are strongly dependent.
"""

import random
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from sarna.sequence.nucleic import NucleicSequence

DEFAULT_PERMUTATIONS = 20_000
"""Shuffles per test. Enough to resolve a family-wise p down to about 0.0001."""

DEFAULT_SEED = 20261006
"""Fixed so a reported p-value is reproducible. Stated in the output, never left implicit."""


@dataclass(frozen=True, slots=True)
class MotifDifference:
    """One cell of the comparison: a motif at a position, and how the two groups differ there."""

    position: int
    """1-based position in the sequence, where the motif starts."""

    motif: str
    """The motif itself, e.g. ``G`` or ``UG``."""

    group_percent: float
    """Percentage of the labelled group carrying this motif here."""

    reference_percent: float
    """Percentage of the reference group carrying it."""

    @property
    def difference(self) -> float:
        """Return the difference in percentage points, group minus reference."""
        return self.group_percent - self.reference_percent

    def __str__(self) -> str:
        """Return a one-line rendering, e.g. ``pos 10 G: +24.8 pp``."""
        return f"pos {self.position} {self.motif}: {self.difference:+.1f} pp"


@dataclass(frozen=True, slots=True)
class PermutationResult:
    """What a max-statistic permutation test concluded about a whole grid of differences."""

    observed_maximum: float
    """The largest absolute difference actually seen, in percentage points."""

    null_median: float
    """The median largest-absolute-difference across shuffles."""

    null_95th_percentile: float
    family_wise_p: float
    """Fraction of shuffles whose maximum reached the observed maximum, plus-one corrected."""

    permutations: int
    seed: int
    cells: int
    """How many position-by-motif cells the maximum was taken over."""

    @property
    def is_significant(self) -> bool:
        """Return whether the observed maximum survives a 0.05 family-wise threshold."""
        return self.family_wise_p < 0.05

    @property
    def exceeds_chance(self) -> bool:
        """Return whether the observed maximum is even above the median of the null.

        Worth reporting separately from significance. A maximum *below* the null median is not
        merely non-significant: the labelled split produced less positional structure than a
        typical random split of the same sequences, which is a stronger statement than "we could
        not detect an effect".
        """
        return self.observed_maximum > self.null_median

    @property
    def summary(self) -> str:
        """Return a plain-language verdict, including the direction against chance."""
        verdict = (
            f"family-wise p = {self.family_wise_p:.3f} over {self.cells} cells"
            f" ({self.permutations} shuffles, seed {self.seed})"
        )
        if self.is_significant:
            return (
                f"The largest difference, {self.observed_maximum:.1f} pp, is larger than chance "
                f"would usually produce: {verdict}."
            )
        if self.exceeds_chance:
            return (
                f"The largest difference, {self.observed_maximum:.1f} pp, is within what chance "
                f"produces: a random split of these sequences reaches "
                f"{self.null_median:.1f} pp on average; {verdict}."
            )
        return (
            f"The largest difference, {self.observed_maximum:.1f} pp, is smaller than chance "
            f"usually produces: a random split of these same sequences reaches "
            f"{self.null_median:.1f} pp on average, so the labels carry no positional signal "
            f"this test can detect; {verdict}."
        )


def motifs_at(sequence: NucleicSequence, position: int, length: int) -> str:
    """Return the ``length``-mer starting at 1-based ``position``.

    Raises:
        ValueError: If the motif would run past the end of the sequence.
    """
    if position < 1 or position + length - 1 > len(sequence):
        raise ValueError(
            f"a {length}-mer at position {position} does not fit in a "
            f"{len(sequence)} nt sequence"
        )
    return sequence.subsequence(position, position + length - 1).bases


def positional_frequencies(
    sequences: Sequence[NucleicSequence], motif_length: int = 1
) -> Mapping[tuple[int, str], float]:
    """Return the percentage of ``sequences`` carrying each motif at each position.

    Args:
        sequences: Equal-length sequences to tabulate. All must be the same length.
        motif_length: 1 for single nucleotides, 2 for dinucleotides, and so on.

    Returns:
        Percentages keyed by (1-based start position, motif). Only observed motifs appear.

    Raises:
        ValueError: If ``sequences`` is empty or the sequences differ in length.
    """
    length = _common_length(sequences)
    counts: Counter[tuple[int, str]] = Counter()
    for sequence in sequences:
        for position in range(1, length - motif_length + 2):
            counts[(position, motifs_at(sequence, position, motif_length))] += 1
    total = len(sequences)
    return {cell: 100.0 * count / total for cell, count in counts.items()}


def frequency_differences(
    group: Sequence[NucleicSequence],
    reference: Sequence[NucleicSequence],
    motif_length: int = 1,
) -> tuple[MotifDifference, ...]:
    """Return every cell's difference between ``group`` and ``reference``, largest first.

    A motif absent from one group counts as 0% there rather than being omitted, so a depletion
    is reported as readily as an enrichment.
    """
    in_group = positional_frequencies(group, motif_length)
    in_reference = positional_frequencies(reference, motif_length)
    differences = tuple(
        MotifDifference(
            position=position,
            motif=motif,
            group_percent=in_group.get((position, motif), 0.0),
            reference_percent=in_reference.get((position, motif), 0.0),
        )
        for position, motif in sorted(set(in_group) | set(in_reference))
    )
    return tuple(sorted(differences, key=lambda d: -abs(d.difference)))


def permutation_test(
    group: Sequence[NucleicSequence],
    reference: Sequence[NucleicSequence],
    motif_length: int = 1,
    permutations: int = DEFAULT_PERMUTATIONS,
    seed: int = DEFAULT_SEED,
) -> PermutationResult:
    """Test whether the largest positional difference exceeds what label shuffling produces.

    The statistic is the maximum absolute difference over the whole position-by-motif grid, so
    the resulting p-value is family-wise across every cell examined. Group sizes are preserved
    under shuffling, which is what makes the null the right comparison for a 14-against-111
    split rather than a balanced one.

    Args:
        group: The labelled group, e.g. the activating sequences.
        reference: The group to compare against.
        motif_length: Motif width to scan.
        permutations: Number of label shuffles.
        seed: Random seed, reported in the result so the number can be reproduced.

    Returns:
        The observed maximum, the null distribution's shape, and the family-wise p-value.

    Raises:
        ValueError: If either group is empty or the sequences differ in length.
    """
    observed = frequency_differences(group, reference, motif_length)
    observed_maximum = abs(observed[0].difference)

    pooled = [*group, *reference]
    _common_length(pooled)
    split = len(group)
    rng = random.Random(seed)
    null: list[float] = []
    for _ in range(permutations):
        shuffled = pooled[:]
        rng.shuffle(shuffled)
        differences = frequency_differences(
            shuffled[:split], shuffled[split:], motif_length
        )
        null.append(abs(differences[0].difference))
    null.sort()

    at_least_as_large = sum(1 for maximum in null if maximum >= observed_maximum)
    return PermutationResult(
        observed_maximum=observed_maximum,
        null_median=_percentile(null, 50),
        null_95th_percentile=_percentile(null, 95),
        # Plus-one correction: with no shuffle reaching the observed value the honest statement
        # is "below 1/(n+1)", not "exactly zero".
        family_wise_p=(at_least_as_large + 1) / (permutations + 1),
        permutations=permutations,
        seed=seed,
        cells=len(observed),
    )


def _common_length(sequences: Iterable[NucleicSequence]) -> int:
    """Return the shared length of ``sequences``.

    Raises:
        ValueError: If empty, or if the sequences are not all the same length. Positional
            analysis of ragged sequences would silently compare different molecules' positions.
    """
    lengths = {len(sequence) for sequence in sequences}
    if not lengths:
        raise ValueError("no sequences to analyse")
    if len(lengths) > 1:
        raise ValueError(
            f"positional analysis needs equal-length sequences; found lengths "
            f"{sorted(lengths)}. Aligning them is a decision for the caller, not this function."
        )
    return lengths.pop()


def _percentile(sorted_values: Sequence[float], percentile: float) -> float:
    """Return a percentile of an already-sorted sequence by nearest rank.

    Callers pass a non-empty list of permutation maxima; ``permutations`` is at least one.
    """
    rank = round(percentile / 100.0 * (len(sorted_values) - 1))
    index = min(len(sorted_values) - 1, max(0, rank))
    return sorted_values[index]


__all__ = [
    "DEFAULT_PERMUTATIONS",
    "DEFAULT_SEED",
    "MotifDifference",
    "PermutationResult",
    "frequency_differences",
    "motifs_at",
    "permutation_test",
    "positional_frequencies",
]
