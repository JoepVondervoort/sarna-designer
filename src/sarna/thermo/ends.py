"""Stacking stability of a duplex's terminal ends.

The quantity that decides which strand RISC loads is not the stability of the whole duplex but of
its two **ends**. Schwarz et al. 2003 report that the unwinding helicase reads roughly four bases
at a time, which is why the legacy script defines an end as four nucleotides — three
nearest-neighbour stacks — and the default here is the same.

Energies are sums of stacking free energies from a chosen :class:`~sarna.thermo.model.
NearestNeighbourTable` and are therefore negative: a **more negative** value is a more stable end.
No initiation or terminal penalty is applied, because the quantity of interest is the difference
between two ends of the same duplex, where those terms cancel. This also matches the legacy
calculation, which omits them entirely.

The per-stack breakdown travels with the result. It is what lets the application answer "why is
this end less stable?" with the actual steps it summed rather than a bare number (ADR 0009).
"""

from dataclasses import dataclass
from enum import Enum

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.model import NearestNeighbourTable

DEFAULT_END_LENGTH = 4
"""Nucleotides defining one end, after Schwarz et al. 2003's helicase footprint."""


class DuplexEnd(Enum):
    """Which end of the target sequence to measure."""

    FIVE_PRIME = "5'"
    THREE_PRIME = "3'"

    def __str__(self) -> str:
        """Return ``5'`` or ``3'``."""
        return self.value


@dataclass(frozen=True, slots=True)
class EndStability:
    """The stacking stability of one terminal end, with the steps it was summed from."""

    end: DuplexEnd
    energy: float
    """Summed stacking free energy. Negative; more negative means more stable."""

    stacks: tuple[tuple[str, float], ...]
    """The dinucleotide steps summed, in 5'-to-3' order, each with its energy."""

    unit: str
    """The table's energy unit, e.g. ``kcal/mol``. Carried so the two are never separated."""

    @property
    def breakdown(self) -> str:
        """Return the summed steps as text, e.g. ``GC -9.36 + CG -9.07 + GC -9.36``."""
        return " + ".join(f"{step} {energy:.2f}" for step, energy in self.stacks)


def end_stability(
    target: NucleicSequence,
    table: NearestNeighbourTable,
    end: DuplexEnd,
    end_length: int = DEFAULT_END_LENGTH,
) -> EndStability:
    """Return the stacking stability of one end of ``target``.

    The target is converted to the table's alphabet first: an RNA parameter set describes the
    duplex that is actually loaded, while the target site itself is DNA. T-to-U is a one-to-one
    relabelling, so the steps correspond exactly.

    Args:
        target: The target site, normally DNA.
        table: The stacking parameters to use.
        end: Which end to measure.
        end_length: Nucleotides defining the end. ``end_length - 1`` stacks are summed.

    Returns:
        The end's stability, with its per-stack breakdown.

    Raises:
        ValueError: If ``end_length`` is below 2 or longer than the target.
    """
    if end_length < 2:
        raise ValueError(
            f"an end must span at least 2 nucleotides to contain one stack, got {end_length}"
        )
    if end_length > len(target):
        raise ValueError(
            f"an end of {end_length} nt does not fit in a {len(target)} nt target"
        )
    sequence = target.as_rna() if table.alphabet is Nucleic.RNA else target.as_dna()
    bases = sequence.bases
    terminal = bases[:end_length] if end is DuplexEnd.FIVE_PRIME else bases[-end_length:]
    stacks = tuple(
        (terminal[index : index + 2], table.stack(terminal[index : index + 2]))
        for index in range(end_length - 1)
    )
    return EndStability(
        end=end,
        energy=sum(energy for _, energy in stacks),
        stacks=stacks,
        unit=table.unit,
    )
