"""Sequence primitives: the validated, alphabet-tagged building blocks of every candidate.

This layer knows nothing about rules, presets, genomes or scoring. It only knows what a
nucleotide sequence is and what can be measured about one.
"""

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.cpg import (
    GARDINER_GARDEN_FROMMER_1987,
    TAKAI_JONES_2002,
    CpGIsland,
    CpGIslandCriteria,
    find_cpg_islands,
    overlaps_any_island,
)
from sarna.sequence.nucleic import InvalidSequenceError, NucleicSequence

__all__ = [
    "GARDINER_GARDEN_FROMMER_1987",
    "TAKAI_JONES_2002",
    "CpGIsland",
    "CpGIslandCriteria",
    "InvalidSequenceError",
    "Nucleic",
    "NucleicSequence",
    "find_cpg_islands",
    "overlaps_any_island",
]
