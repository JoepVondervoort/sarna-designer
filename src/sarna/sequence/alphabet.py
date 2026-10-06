"""Which nucleic acid a sequence or a thermodynamic table is written in.

This is the most primitive concept in the engine, so it lives at the bottom: the sequence
primitives, the nearest-neighbour tables and the rule evaluators all depend on it, and it
depends on nothing. It was originally defined inside :mod:`sarna.thermo.model`, which inverted
that relationship.
"""

from collections.abc import Mapping
from enum import Enum


class Nucleic(Enum):
    """Whether a sequence or parameter set is DNA or RNA.

    The distinction is not cosmetic. A promoter target site is DNA and its CpG content is a
    meaningful property of it; the duplex delivered to a cell is RNA and has no CpGs. The legacy
    Scanner converted everything to RNA on input (``checksequence`` replaces T with U), after
    which "is this inside a CpG island" is a question about a sequence that no longer records
    the information needed to answer it.
    """

    RNA = "RNA"
    DNA = "DNA"

    def __str__(self) -> str:
        """Return ``RNA`` or ``DNA``."""
        return self.value

    @property
    def thymine_analogue(self) -> str:
        """Return ``U`` for RNA and ``T`` for DNA."""
        return "U" if self is Nucleic.RNA else "T"

    @property
    def bases(self) -> tuple[str, ...]:
        """Return the four bases of this alphabet, in alphabetical order."""
        return ("A", "C", "G", "U") if self is Nucleic.RNA else ("A", "C", "G", "T")

    @property
    def complement(self) -> Mapping[str, str]:
        """Return the base-pairing map for this alphabet."""
        partner = self.thymine_analogue
        return {"A": partner, partner: "A", "G": "C", "C": "G"}

    @property
    def other(self) -> "Nucleic":
        """Return the opposite alphabet."""
        return Nucleic.DNA if self is Nucleic.RNA else Nucleic.RNA
