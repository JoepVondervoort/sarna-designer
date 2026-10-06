"""The promoter region a design run scans, and the mapping from a window to coordinates.

A :class:`Promoter` holds its sequence **5' to 3' in the transcript's own orientation**. For a
minus-strand gene that is the reverse complement of the reference sequence, and
:meth:`Promoter.from_reference_sequence` performs that flip at the boundary. The consequence is
worth stating plainly: past this point, **nothing in the engine needs to know about strand**. A
candidate at offset 37 is at offset 37, its TSS-relative position is the same number on either
strand, and only the mapping back to contig coordinates differs.

That mapping is the one piece of strand arithmetic in the engine, and it is defect B1's home. The
legacy EpiMarker requested a minus-strand region from Ensembl, received the reverse complement,
then computed ``abspos = start + idx`` as though the index counted forward from the low
coordinate. It does not: on the minus strand the transcript's 5' end is the **highest**
coordinate, so offsets run downwards. Every minus-strand position and every TSS-relative distance
that pipeline reported was wrong.

Both anchors are optional and independent:

- **No interval** — a pasted sequence with no genomic context (Workflow B). Candidates are still
  designed; coordinate-dependent annotations are unavailable.
- **No TSS** — the promoter's position relative to a start site is unknown, so the TSS window
  rules report ``not-applicable`` rather than passing (ADR 0007).
"""

from dataclasses import dataclass

from sarna.genome.coordinates import GenomicInterval
from sarna.genome.tss import TssPosition
from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence


@dataclass(frozen=True, slots=True)
class Promoter:
    """A promoter region to scan, optionally anchored to a genome and to a TSS.

    Attributes:
        sequence: The promoter as DNA, 5' to 3' **in the transcript's orientation**.
        tss_start: TSS-relative position of the sequence's first base, e.g. ``-1000``.
        interval: Genomic span of the whole promoter, on the transcript's strand.
    """

    sequence: NucleicSequence
    tss_start: TssPosition | None = None
    interval: GenomicInterval | None = None

    def __post_init__(self) -> None:
        """Validate the sequence alphabet and that any interval matches the sequence length.

        Raises:
            ValueError: If the sequence is not DNA, or if an interval is given whose length
                differs from the sequence length — a mismatch that would silently shift every
                candidate's coordinates.
        """
        if self.sequence.alphabet is not Nucleic.DNA:
            raise ValueError(
                "a promoter is DNA, but the sequence given is RNA; call .as_dna() if the "
                "sequence really is the promoter and not the duplex"
            )
        if self.interval is not None and self.interval.length != len(self.sequence):
            raise ValueError(
                f"the interval spans {self.interval.length} bp but the sequence is "
                f"{len(self.sequence)} nt; these must agree, or every candidate's coordinates "
                f"would be shifted"
            )

    @classmethod
    def from_reference_sequence(
        cls,
        reference_sequence: NucleicSequence,
        interval: GenomicInterval,
        tss_start: TssPosition | None = None,
    ) -> "Promoter":
        """Build a promoter from a sequence read off the **plus strand** of the reference.

        This is the boundary where strand is resolved. A reference genome serves the plus strand;
        for a minus-strand transcript the promoter's own 5'-to-3' sequence is its reverse
        complement, and performing that flip here is what keeps the rest of the engine
        strand-agnostic.

        Args:
            reference_sequence: The span's sequence as the reference stores it, plus strand.
            interval: The span, carrying the transcript's strand.
            tss_start: TSS-relative position of the promoter's first base, if known.

        Returns:
            A promoter whose sequence runs 5' to 3' in the transcript's orientation.
        """
        oriented = (
            reference_sequence
            if interval.strand.sign > 0
            else reference_sequence.reverse_complement()
        )
        return cls(sequence=oriented, tss_start=tss_start, interval=interval)

    @property
    def length(self) -> int:
        """Return the promoter length in base pairs."""
        return len(self.sequence)

    @property
    def is_anchored(self) -> bool:
        """Return whether this promoter's position relative to a TSS is known."""
        return self.tss_start is not None

    def _check_window(self, offset: int, length: int) -> None:
        """Validate that a window of ``length`` at ``offset`` lies inside this promoter.

        Raises:
            ValueError: If the length is not positive or the window runs outside the promoter.
        """
        if length <= 0:
            raise ValueError(f"window length must be positive, got {length}")
        if offset < 0 or offset + length > self.length:
            raise ValueError(
                f"a {length} nt window at offset {offset} does not fit inside a "
                f"{self.length} bp promoter"
            )

    def tss_span(self, offset: int, length: int) -> tuple[TssPosition, TssPosition] | None:
        """Return the TSS-relative positions of a window's first and last base.

        The span is closed: a 19 nt window starting at ``-347`` ends at ``-329``, not ``-328``.
        The result is strand-independent, because TSS-relative positions are measured in the
        transcript's own frame (ADR 0008).

        Args:
            offset: 0-based offset of the window into the promoter sequence.
            length: Window length in nucleotides.

        Returns:
            ``(first, last)`` positions, or ``None`` if this promoter has no known TSS.

        Raises:
            ValueError: If the window does not fit inside the promoter.
        """
        self._check_window(offset, length)
        if self.tss_start is None:
            return None
        first = self.tss_start.shifted_by(offset)
        return first, first.shifted_by(length - 1)

    def genomic_interval(self, offset: int, length: int) -> GenomicInterval | None:
        """Return the genomic span of a window, on the promoter's strand.

        On the plus strand an offset moves the span to higher coordinates; on the minus strand it
        moves it to **lower** ones, because the transcript's 5' end is the highest coordinate of
        the promoter. This is the whole of the engine's strand arithmetic.

        Args:
            offset: 0-based offset of the window into the promoter sequence.
            length: Window length in nucleotides.

        Returns:
            The window's half-open genomic interval, or ``None`` for an unanchored promoter.

        Raises:
            ValueError: If the window does not fit inside the promoter.
        """
        self._check_window(offset, length)
        if self.interval is None:
            return None
        if self.interval.strand.sign > 0:
            start = self.interval.start + offset
        else:
            start = self.interval.end - offset - length
        return GenomicInterval(
            assembly=self.interval.assembly,
            contig=self.interval.contig,
            start=start,
            end=start + length,
            strand=self.interval.strand,
        )

    @property
    def tss_coordinate(self) -> int | None:
        """Return the genomic coordinate of the first transcribed base, if determinable.

        Requires both anchors: an interval to measure from and a TSS-relative start to measure.
        On the minus strand the TSS lies at a *lower* coordinate than the promoter, since the
        promoter is upstream and upstream means higher there.

        Returns:
            The 0-based coordinate of position ``+1``, or ``None`` if either anchor is missing.
        """
        if self.interval is None or self.tss_start is None:
            return None
        # Distance in base pairs from the promoter's first base to the TSS, in the transcript's
        # own direction of travel; then stepped off along the strand.
        to_tss = self.tss_start.distance_to(TssPosition(1))
        return self.interval.five_prime + to_tss * self.interval.strand.sign
