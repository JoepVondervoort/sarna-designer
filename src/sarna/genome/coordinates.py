"""Assembly-tagged genomic coordinates.

Defect B15 in the legacy pipeline was the silent combination of GRCh38 transcript coordinates
with GRCh37 signal tracks. ADR 0004 responds by making a coordinate without an assembly
*unrepresentable*: :class:`GenomicInterval` requires one, refuses to compare or measure across
assemblies, and validates itself against that assembly's contig lengths.

Intervals are **0-based, half-open**, ``[start, end)``, so ``(0, 19)`` is the first 19 bases.
This matches bigWig and BED. Ensembl REST is 1-based inclusive, so its client converts at the
boundary rather than letting two conventions coexist in the model.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

from sarna.genome.contigs import GRCH37_CONTIG_LENGTHS, GRCH38_CONTIG_LENGTHS


class AssemblyMismatchError(ValueError):
    """Raised when an operation would combine coordinates from two different assemblies.

    This is deliberately an error rather than an implicit conversion. Converting requires a
    chain file and may fail or be ambiguous, so it is the explicit job of the liftover adapter,
    never a silent side effect of arithmetic.
    """


class Assembly(Enum):
    """A human genome reference assembly."""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"

    def __str__(self) -> str:
        """Return the conventional assembly name, e.g. ``GRCh38``."""
        return self.value

    @property
    def contig_lengths(self) -> Mapping[str, int]:
        """Return the primary-assembly contig lengths for this assembly."""
        return (
            GRCH37_CONTIG_LENGTHS if self is Assembly.GRCH37 else GRCH38_CONTIG_LENGTHS
        )

    def contig_length(self, contig: str) -> int:
        """Return the length of ``contig`` in this assembly.

        Args:
            contig: Ensembl-style contig name, e.g. ``"1"``, ``"X"``, ``"MT"``.

        Returns:
            The contig length in base pairs.

        Raises:
            ValueError: If the contig is not part of this assembly's primary set.
        """
        try:
            return self.contig_lengths[contig]
        except KeyError:
            raise ValueError(
                f"unknown contig {contig!r} for {self}; expected an Ensembl-style primary "
                f"contig name without the 'chr' prefix"
            ) from None


class Strand(Enum):
    """The DNA strand a feature lies on."""

    PLUS = "+"
    MINUS = "-"

    def __str__(self) -> str:
        """Return ``+`` or ``-``."""
        return self.value

    @property
    def sign(self) -> int:
        """Return ``+1`` on the plus strand and ``-1`` on the minus strand.

        Used to orient distances so that "upstream" is negative on both strands.
        """
        return 1 if self is Strand.PLUS else -1


@dataclass(frozen=True, slots=True)
class GenomicInterval:
    """A half-open ``[start, end)`` interval on a named assembly and strand.

    Frozen, so an interval can be passed around without any risk of a caller mutating it.
    Operations that would combine assemblies raise :class:`AssemblyMismatchError`.
    """

    assembly: Assembly
    contig: str
    start: int
    end: int
    strand: Strand = Strand.PLUS

    def __post_init__(self) -> None:
        """Validate the interval against its assembly.

        Raises:
            ValueError: If the contig is unknown or ``chr``-prefixed, the coordinates are
                negative or empty, or the interval runs past the end of the contig.
        """
        if self.contig.startswith("chr"):
            raise ValueError(
                f"contig {self.contig!r} is UCSC-style; the domain model uses Ensembl-style "
                f"names without the 'chr' prefix (adapters normalise this)"
            )
        contig_length = self.assembly.contig_length(self.contig)
        if self.start < 0:
            raise ValueError(f"start must not be negative, got {self.start}")
        if self.end <= self.start:
            raise ValueError(
                f"end must be greater than start, got start={self.start}, end={self.end}"
            )
        if self.end > contig_length:
            raise ValueError(
                f"interval ends at {self.end}, beyond the end of contig {self.contig!r} "
                f"on {self.assembly} (length {contig_length})"
            )

    @property
    def length(self) -> int:
        """Return the interval length in base pairs."""
        return self.end - self.start

    @property
    def five_prime(self) -> int:
        """Return the coordinate of the 5' end *in the feature's own orientation*.

        On the plus strand this is ``start``; on the minus strand it is the **highest**
        coordinate, ``end - 1``. Assuming otherwise is defect B1.
        """
        return self.start if self.strand is Strand.PLUS else self.end - 1

    def _require_same_assembly(self, other: "GenomicInterval") -> None:
        if self.assembly is not other.assembly:
            raise AssemblyMismatchError(
                f"cannot combine coordinates on {self.assembly} with coordinates on "
                f"{other.assembly}; convert explicitly via the liftover adapter first"
            )

    def overlaps(self, other: "GenomicInterval") -> bool:
        """Return whether this interval shares at least one base with ``other``.

        Strand is ignored: overlap is a statement about position, not orientation.

        Raises:
            AssemblyMismatchError: If the two intervals are on different assemblies.
        """
        self._require_same_assembly(other)
        if self.contig != other.contig:
            return False
        return self.start < other.end and other.start < self.end

    def signed_distance_to(self, reference: "GenomicInterval") -> int:
        """Return the distance from ``reference`` to this interval, oriented by strand.

        Measured between the two 5' ends in the feature's own orientation, so that
        **upstream is negative and downstream positive on both strands**. A target 500 bp
        upstream of a TSS reports ``-500`` whether the gene is on the plus or minus strand.

        Args:
            reference: The anchor to measure from, typically a transcription start site.

        Returns:
            Signed distance in base pairs.

        Raises:
            AssemblyMismatchError: If the two intervals are on different assemblies.
            ValueError: If they are on different contigs or different strands, where a
                signed distance has no meaning.
        """
        self._require_same_assembly(reference)
        if self.contig != reference.contig:
            raise ValueError(
                f"cannot measure a distance between contigs {self.contig!r} and "
                f"{reference.contig!r}"
            )
        if self.strand is not reference.strand:
            raise ValueError(
                f"cannot measure an oriented distance between features on opposite strands "
                f"({self.strand} and {reference.strand})"
            )
        return (self.five_prime - reference.five_prime) * self.strand.sign
