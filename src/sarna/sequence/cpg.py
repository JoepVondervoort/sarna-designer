"""CpG-island detection from sequence.

A CpG island is a stretch of DNA unusually rich in CpG dinucleotides relative to what its base
composition predicts. Two definitions are in circulation and this module ships both, because the
choice materially changes which candidates a design run excludes:

``TAKAI_JONES_2002``
    The default. At least 500 bp, GC at least 55 %, observed/expected CpG at least 0.65. Takai &
    Jones proposed these precisely because the older thresholds over-call: applying them to human
    chromosomes 21 and 22 produced 14 062 islands of which more than half were Alu repeats, where
    the stricter criteria gave 1 101.

``GARDINER_GARDEN_FROMMER_1987``
    The original and still the most widely cited definition: at least 200 bp, GC above 50 %,
    observed/expected CpG above 0.60. Also the thresholds the legacy script's own detector used.

The distinction matters more here than it would in an annotation pipeline, because this feature
**excludes** candidates. A mis-called island costs a usable candidate for no biological reason, so
the project prefers the definition that calls fewer, better-supported islands.

Note the convention difference against the published criteria: both definitions are implemented
with inclusive comparisons, where Gardiner-Garden & Frommer state theirs as strict inequalities.
The two differ only for a window sitting exactly on a threshold.

**Locally computed, not annotated.** These islands are derived from sequence by this module. They
are not the same thing as an island record from a genome annotation database, and the application
keeps the two conceptually separate.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence

DEFAULT_MINIMUM_CPG_COUNT = 7
"""Observed CpGs a window must contain regardless of its ratio.

A GC-rich window whose Cs and Gs are segregated (``CCCC…GGGG``) has a very low *expected* CpG
count, so the observed-to-expected ratio can clear any threshold on a handful of real CpGs. The
guard rejects that arithmetic artefact. The legacy detector carries the same guard and states the
same reason, though it tests ``> 7`` where its own comment says "at least 7".
"""


@dataclass(frozen=True, slots=True)
class CpGIslandCriteria:
    """Thresholds defining what counts as a CpG island."""

    minimum_length: int
    """Shortest qualifying window, which is also the sliding-window size."""

    minimum_gc_fraction: float
    minimum_observed_expected: float
    minimum_cpg_count: int = DEFAULT_MINIMUM_CPG_COUNT

    def __post_init__(self) -> None:
        """Validate the thresholds.

        Raises:
            ValueError: If the length is not positive, the GC fraction is outside 0 to 1, or
                either the ratio or the count threshold is negative.
        """
        if self.minimum_length <= 0:
            raise ValueError(
                f"minimum island length must be positive, got {self.minimum_length}"
            )
        if not 0.0 <= self.minimum_gc_fraction <= 1.0:
            raise ValueError(
                f"minimum GC fraction must lie between 0 and 1, got "
                f"{self.minimum_gc_fraction}"
            )
        if self.minimum_observed_expected < 0:
            raise ValueError(
                f"minimum observed/expected CpG ratio must not be negative, got "
                f"{self.minimum_observed_expected}"
            )
        if self.minimum_cpg_count < 0:
            raise ValueError(
                f"minimum CpG count must not be negative, got {self.minimum_cpg_count}"
            )


TAKAI_JONES_2002 = CpGIslandCriteria(
    minimum_length=500,
    minimum_gc_fraction=0.55,
    minimum_observed_expected=0.65,
)
"""Takai & Jones 2002, PNAS 99(6):3740-3745. The project default."""

GARDINER_GARDEN_FROMMER_1987 = CpGIslandCriteria(
    minimum_length=200,
    minimum_gc_fraction=0.50,
    minimum_observed_expected=0.60,
)
"""Gardiner-Garden & Frommer 1987, J Mol Biol 196(2):261-282. Also the legacy script's values."""


@dataclass(frozen=True, slots=True)
class CpGIsland:
    """A detected island, as a half-open ``[start, end)`` span of sequence offsets.

    The statistics describe the **merged** span, not the first window that qualified. That window
    is by construction the least GC-rich part of a long island, so reporting its numbers would
    understate the island.
    """

    start: int
    end: int
    gc_fraction: float
    observed_cpg: int
    expected_cpg: float

    @property
    def length(self) -> int:
        """Return the island length in base pairs."""
        return self.end - self.start

    @property
    def observed_expected_ratio(self) -> float:
        """Return observed CpGs over expected, or 0 where none are expected."""
        return self.observed_cpg / self.expected_cpg if self.expected_cpg else 0.0

    def overlaps(self, start: int, end: int) -> bool:
        """Return whether this island shares at least one base with ``[start, end)``."""
        return start < self.end and self.start < end


def _qualifies(
    window_gc: int,
    window_cpg: int,
    window_c: int,
    window_g: int,
    window_length: int,
    criteria: CpGIslandCriteria,
) -> bool:
    """Return whether one window satisfies every criterion."""
    if window_cpg < criteria.minimum_cpg_count:
        return False
    if window_gc / window_length < criteria.minimum_gc_fraction:
        return False
    expected = (window_c * window_g) / window_length
    if expected == 0:
        return False
    return window_cpg / expected >= criteria.minimum_observed_expected


def find_cpg_islands(
    sequence: NucleicSequence, criteria: CpGIslandCriteria = TAKAI_JONES_2002
) -> tuple[CpGIsland, ...]:
    """Return the CpG islands in ``sequence``, in sequence order.

    A window of ``criteria.minimum_length`` is slid one base at a time; qualifying windows are
    merged where they overlap, and each merged span's statistics are then computed over the span
    as a whole. All counts are maintained incrementally, so this is linear in the sequence length
    the legacy detector rescans the window on every step, making it O(n · window) (defect B7).

    Args:
        sequence: The promoter or region to scan. Must be DNA.
        criteria: The island definition to apply.

    Returns:
        The islands found, possibly empty. A sequence shorter than the window yields none: that
        is a fact about the sequence, so no warning is emitted and the caller decides what to
        report.

    Raises:
        ValueError: If the sequence is RNA, where CpG content is not defined.
    """
    if sequence.alphabet is not Nucleic.DNA:
        raise ValueError(
            "CpG islands are a property of DNA; this sequence is RNA. Call .as_dna() on the "
            "target region if that is what you meant"
        )
    bases = sequence.bases
    window = criteria.minimum_length
    if len(bases) < window:
        return ()

    gc = sum(1 for base in bases[:window] if base in "GC")
    c_count = bases.count("C", 0, window)
    g_count = bases.count("G", 0, window)
    cpg = bases.count("CG", 0, window)

    spans: list[list[int]] = []
    for start in range(len(bases) - window + 1):
        if start > 0:
            leaving, entering = bases[start - 1], bases[start + window - 1]
            gc += (entering in "GC") - (leaving in "GC")
            c_count += (entering == "C") - (leaving == "C")
            g_count += (entering == "G") - (leaving == "G")
            # A CpG leaves when the departing base was the C of a pair, and one arrives when the
            # entering base completes a pair with its predecessor.
            cpg -= leaving == "C" and bases[start] == "G"
            cpg += bases[start + window - 2] == "C" and entering == "G"
        if not _qualifies(gc, cpg, c_count, g_count, window, criteria):
            continue
        end = start + window
        if spans and start <= spans[-1][1]:
            spans[-1][1] = end
        else:
            spans.append([start, end])

    return tuple(_describe(bases, start, end) for start, end in spans)


def _describe(bases: str, start: int, end: int) -> CpGIsland:
    """Compute a merged span's statistics over the whole span."""
    span = bases[start:end]
    length = end - start
    c_count = span.count("C")
    g_count = span.count("G")
    return CpGIsland(
        start=start,
        end=end,
        gc_fraction=(c_count + g_count) / length,
        observed_cpg=span.count("CG"),
        expected_cpg=(c_count * g_count) / length,
    )


def overlaps_any_island(
    islands: Sequence[CpGIsland], start: int, end: int
) -> bool:
    """Return whether ``[start, end)`` shares a base with any of ``islands``.

    **Any** overlap counts, not containment. The asymmetry against the TSS window rule, which
    requires total containment (ADR 0008), is deliberate: both choices are the cautious one in
    their own direction. A window the design is *permitted* to use must be entirely inside it; a
    region the design must *avoid* is avoided entirely. A methylatable CpG anywhere in the target
    site can block the duplex from engaging, so a partial overlap is still a reason for concern.

    Args:
        islands: Islands to test against, typically from :func:`find_cpg_islands`.
        start: 0-based inclusive start of the span.
        end: 0-based exclusive end of the span.

    Returns:
        Whether the span touches any island.

    Raises:
        ValueError: If the span is empty or inverted.
    """
    if end <= start:
        raise ValueError(f"span [{start}, {end}) is empty or inverted")
    return any(island.overlaps(start, end) for island in islands)
