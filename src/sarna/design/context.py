"""Everything a rule evaluator needs to judge one candidate.

ADR 0007 settles the evaluator contract: every rule is a pure function of a
:class:`CandidateContext` and its resolved parameters. The context exists so that one uniform
signature serves rules with very different appetites — ``gc-content`` reads only the target,
``flank-20-au`` needs the promoter beyond it, ``tss-window`` needs a strand-aware distance to a
start site — and so that promoter-level work happens once per promoter rather than once per rule
per window.

Two conventions are enforced here rather than in each evaluator:

**Positions are Table 1's positions.** :meth:`CandidateContext.base_at` numbers the target 1 to
19, and :meth:`CandidateContext.flank_base_at` continues into the promoter at 20 to 23, which is
how Wang et al. 2015 numbers them. No evaluator does index arithmetic.

**Absent is not false.** Every field that depends on genomic anchoring is ``None`` when the
promoter has none, and :meth:`CandidateContext.tss_window_contains` returns ``None`` rather than
``False``. A rule reading ``None`` must report ``not-applicable`` (ADR 0007); returning ``False``
here would let a sequence-only run reject every candidate for a check that never ran.

CpG islands are passed **in** rather than computed here. They are a property of the whole
promoter, so detecting them once per promoter and sharing the result across every candidate keeps
the scan linear; computing them per candidate would make it quadratic. Passing ``None`` means
they were not computed, which is distinct from an empty tuple meaning none were found — a rule
reading ``None`` reports not-applicable rather than passing.

Not yet present, deliberately: the thermodynamic duplex features. Those are assessed per
nearest-neighbour table (ADR 0009) rather than being a single value, so they belong with the
candidate rather than its context.
"""

from dataclasses import dataclass

from sarna.genome.coordinates import GenomicInterval
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.sequence.cpg import CpGIsland, overlaps_any_island
from sarna.sequence.nucleic import NucleicSequence

FLANK_LENGTH = 4
"""Bases 3' of the target that Table 1's flank preferences cover: positions 20 to 23."""


@dataclass(frozen=True, slots=True)
class CandidateContext:
    """One candidate window, with its surroundings, ready to be evaluated.

    Build with :meth:`at` rather than directly, so the derived fields cannot disagree with the
    promoter and offset they came from.
    """

    promoter: Promoter
    """The promoter this window was taken from, 5' to 3' in the transcript's orientation."""

    offset: int
    """0-based offset of the target into the promoter sequence."""

    target: NucleicSequence
    """The target site itself, as DNA. Table 1 positions 1 to 19."""

    flank: NucleicSequence | None
    """Up to four bases 3' of the target, or ``None`` if the promoter ends at the target."""

    tss_span: tuple[TssPosition, TssPosition] | None
    """TSS-relative positions of the target's first and last base, inclusive, if anchored."""

    interval: GenomicInterval | None
    """The target's genomic span on the transcript's strand, if the promoter is anchored."""

    cpg_islands: tuple[CpGIsland, ...] | None
    """Islands in the whole promoter, or ``None`` if island detection was not run."""

    @classmethod
    def at(
        cls,
        promoter: Promoter,
        offset: int,
        length: int,
        flank_length: int = FLANK_LENGTH,
        cpg_islands: tuple[CpGIsland, ...] | None = None,
    ) -> "CandidateContext":
        """Build the context for the window of ``length`` at ``offset`` in ``promoter``.

        Args:
            promoter: The promoter being scanned.
            offset: 0-based offset of the window into the promoter sequence.
            length: Target length in nucleotides, e.g. 19.
            flank_length: How many bases 3' of the target to carry. Truncated silently where the
                promoter is too short, because a short flank is a fact about the promoter rather
                than an error; the affected preferences then report not-applicable.
            cpg_islands: Islands detected across the whole promoter, computed once by the
                caller. ``None`` means detection was not run, which the island rule reports as
                not-applicable rather than as a pass.

        Returns:
            The context for that window.

        Raises:
            ValueError: If the window does not fit inside the promoter.
        """
        promoter.validate_window(offset, length)
        target = promoter.sequence.subsequence(offset + 1, offset + length)
        flank_start = offset + length + 1
        flank_end = min(flank_start + flank_length - 1, promoter.length)
        flank = (
            promoter.sequence.subsequence(flank_start, flank_end)
            if flank_start <= flank_end
            else None
        )
        return cls(
            promoter=promoter,
            offset=offset,
            target=target,
            flank=flank,
            tss_span=promoter.tss_span(offset=offset, length=length),
            interval=promoter.genomic_interval(offset=offset, length=length),
            cpg_islands=cpg_islands,
        )

    @property
    def length(self) -> int:
        """Return the target length in nucleotides."""
        return len(self.target)

    def base_at(self, position: int) -> str:
        """Return the target base at Table 1 ``position``, counted from 1.

        Args:
            position: 1-based position within the target, 1 to :attr:`length`.

        Returns:
            A single-character base.

        Raises:
            IndexError: If the position lies outside the target. Positions beyond the target are
                the flank's business — see :meth:`flank_base_at`.
        """
        return self.target.base_at(position)

    def flank_base_at(self, position: int) -> str | None:
        """Return the flank base at Table 1 ``position``, or ``None`` if the promoter ends first.

        Table 1 numbers the flank continuously with the target, so for a 19 nt target the flank
        occupies positions 20 to 23.

        Args:
            position: 1-based position, greater than :attr:`length`.

        Returns:
            The base, or ``None`` if the promoter does not extend that far. ``None`` means the
            preference cannot be evaluated, not that it was missed.

        Raises:
            ValueError: If the position lies inside the target, which would make the caller's
                intent ambiguous.
        """
        if position <= self.length:
            raise ValueError(
                f"position {position} is inside the {self.length} nt target, not its flank; "
                f"use base_at for target positions"
            )
        index = position - self.length
        if self.flank is None or index > len(self.flank):
            return None
        return self.flank.base_at(index)

    @property
    def has_tss_context(self) -> bool:
        """Return whether TSS-relative rules can be evaluated for this candidate."""
        return self.tss_span is not None

    @property
    def has_coordinates(self) -> bool:
        """Return whether this candidate has genomic coordinates, for annotation and export."""
        return self.interval is not None

    def tss_window_contains(self, start_bp: int, end_bp: int) -> bool | None:
        """Return whether the **whole** target lies within a TSS-relative window.

        Containment is total, per ADR 0008: a target whose 5' end is inside the window but whose
        3' end runs past it does not lie in the window. Both the mandatory ``tss-window`` rule and
        the ``tss-preferred-window`` preference use this one implementation.

        Args:
            start_bp: Furthest upstream edge of the window, e.g. ``-1000``.
            end_bp: Nearest edge of the window, e.g. ``-100``.

        Returns:
            ``True`` or ``False``, or ``None`` if this promoter has no known TSS and the question
            cannot be answered.

        Raises:
            ValueError: If the bounds are the wrong way round.
        """
        if start_bp > end_bp:
            raise ValueError(
                f"window start {start_bp} lies downstream of its end {end_bp}; bounds run "
                f"from the furthest upstream edge to the nearest"
            )
        if self.tss_span is None:
            return None
        first, last = self.tss_span
        window_start = TssPosition(start_bp)
        window_end = TssPosition(end_bp)
        return window_start <= first and last <= window_end

    @property
    def overlaps_cpg_island(self) -> bool | None:
        """Return whether the target touches a CpG island, or ``None`` if none were detected.

        Any overlap counts, not containment: see
        :func:`sarna.sequence.cpg.overlaps_any_island` for why the two directions differ.
        """
        if self.cpg_islands is None:
            return None
        return overlaps_any_island(
            self.cpg_islands, start=self.offset, end=self.offset + self.length
        )
