"""Positions measured relative to a transcription start site.

ADR 0008 fixes the convention: ``+1`` is the first transcribed base, ``-1`` is the base
immediately upstream, and **there is no position 0**. That is what Wang et al. 2015 and the
promoter literature mean by "-500 bp", and what the interface's ruler shows, so the engine
computes in the same units it displays.

The convention has a discontinuity at the TSS, which is a genuine source of off-by-one error.
This module is the only place that discontinuity exists. Every operation converts to a
**continuous internal offset** first:

===============  ===============
TSS position     internal offset
===============  ===============
``-2``           ``-2``
``-1``           ``-1``
``+1``           ``0``
``+2``           ``1``
===============  ===============

Arithmetic on a raw ``int`` TSS position anywhere else in the engine is a defect.

This type says nothing about strand or assembly. It is a position in a transcript's own frame of
reference, which is what makes it strand-free: a position 500 bp upstream of the TSS is ``-500``
whether the gene is on the plus or the minus strand. Mapping it onto a contig is the job of
:mod:`sarna.genome.coordinates`, and that is where strand matters.
"""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True, order=True)
class TssPosition:
    """A position relative to a transcription start site, in the biological convention.

    Ordered, so positions sort naturally from furthest upstream to furthest downstream: the
    integer values happen to sort correctly despite the missing zero, because the gap is between
    ``-1`` and ``+1`` and no value lies in it.
    """

    value: int
    """The position as a researcher would write it: negative upstream, positive from the TSS on."""

    def __post_init__(self) -> None:
        """Reject position 0.

        Raises:
            ValueError: If ``value`` is 0, which is not a position in this convention.
        """
        if self.value == 0:
            raise ValueError(
                "there is no position 0 relative to a TSS: +1 is the first transcribed base "
                "and -1 is the base immediately upstream of it (ADR 0008)"
            )

    def __str__(self) -> str:
        """Return the position with an explicit sign, e.g. ``-347`` or ``+1``.

        The sign is the information (it says which side of the TSS the position is on) so it is
        never omitted for downstream positions.
        """
        return f"{self.value:+d}"

    @property
    def offset(self) -> int:
        """Return the continuous internal offset for this position.

        Upstream positions are unchanged; downstream positions are one lower, closing the gap
        left by the missing zero. Use this for any arithmetic, then convert back with
        :meth:`from_offset`.
        """
        return self.value - 1 if self.value > 0 else self.value

    @classmethod
    def from_offset(cls, offset: int) -> "TssPosition":
        """Return the TSS position for a continuous internal ``offset``.

        Args:
            offset: A continuous offset, as produced by :attr:`offset`.

        Returns:
            The corresponding position in the biological convention.
        """
        return cls(offset + 1 if offset >= 0 else offset)

    def distance_to(self, other: "TssPosition") -> int:
        """Return the signed distance in base pairs from this position to ``other``.

        Positive downstream, negative upstream. The distance from ``-1`` to ``+1`` is 1, because
        they are adjacent bases.
        """
        return other.offset - self.offset

    def shifted_by(self, base_pairs: int) -> "TssPosition":
        """Return the position ``base_pairs`` downstream of this one.

        Negative values shift upstream. Shifting across the TSS skips the non-existent zero, so
        ``TssPosition(-1).shifted_by(1)`` is ``+1``.
        """
        return TssPosition.from_offset(self.offset + base_pairs)

    @property
    def is_upstream(self) -> bool:
        """Return whether this position lies upstream of the transcription start site."""
        return self.value < 0
