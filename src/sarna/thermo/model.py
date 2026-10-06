"""The nearest-neighbour table type.

A stacking table is reference data with a citation, not a pile of module constants. Modelling
it this way lets the application offer more than one table, show the user where each came from,
and — importantly — *detect* a table that is internally inconsistent instead of trusting it.

In the nearest-neighbour model a stack and its reverse complement are the same physical
interaction, so ``AC`` must carry the same energy as ``GU``. A table where they differ has its
parameters mis-assigned. That is exactly defect B2 in the legacy script, which is why
:attr:`NearestNeighbourTable.is_self_consistent` exists.
"""

from collections.abc import Mapping
from itertools import product
from typing import Self

from pydantic import BaseModel, ConfigDict, model_validator

from sarna.rules.provenance import Provenance
from sarna.sequence.alphabet import Nucleic


class NearestNeighbourTable(BaseModel):
    """A set of nearest-neighbour stacking free energies, with its source."""

    model_config = ConfigDict(frozen=True)

    id: str
    name: str
    alphabet: Nucleic
    unit: str
    provenance: Provenance
    stacks: Mapping[str, float]
    initiation: float
    terminal_au_penalty: float
    symmetry_correction: float
    is_recommended: bool = False
    note: str | None = None

    @model_validator(mode="after")
    def _check_stacks_are_complete_and_in_alphabet(self) -> Self:
        """Require all sixteen stacks, and only bases from this alphabet."""
        expected = {a + b for a, b in product(self.alphabet.bases, repeat=2)}
        unknown = set(self.stacks) - expected
        if unknown:
            raise ValueError(
                f"stacks {sorted(unknown)} use bases outside the {self.alphabet} alphabet"
            )
        missing = expected - set(self.stacks)
        if missing:
            raise ValueError(f"table {self.id!r} is missing stacks {sorted(missing)}")
        return self

    def stack(self, step: str) -> float:
        """Return the stacking free energy for a two-base step.

        Raises:
            KeyError: If the step is not in this table.
        """
        try:
            return self.stacks[step]
        except KeyError:
            raise KeyError(f"no stack {step!r} in table {self.id!r}") from None

    def equivalent_stack(self, step: str) -> str:
        """Return the reverse-complement step, which describes the same interaction."""
        complement = self.alphabet.complement
        return complement[step[1]] + complement[step[0]]

    @property
    def inconsistent_stacks(self) -> tuple[tuple[str, str], ...]:
        """Return stack pairs that should be equal but are not.

        Each pair is reported once, in sorted order. An empty result means the table assigns
        every parameter to the interaction it belongs to.
        """
        seen: set[frozenset[str]] = set()
        found: list[tuple[str, str]] = []
        for step in sorted(self.stacks):
            twin = self.equivalent_stack(step)
            key = frozenset((step, twin))
            if step == twin or key in seen:
                continue
            seen.add(key)
            if self.stacks[step] != self.stacks[twin]:
                found.append((step, twin))
        return tuple(found)

    @property
    def is_self_consistent(self) -> bool:
        """Return whether every stack equals its reverse complement, as physics requires."""
        return not self.inconsistent_stacks


__all__ = ["NearestNeighbourTable", "Nucleic"]
