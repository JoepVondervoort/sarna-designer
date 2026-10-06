"""Sequence primitives: the validated, alphabet-tagged building blocks of every candidate.

This layer knows nothing about rules, presets, genomes or scoring. It only knows what a
nucleotide sequence is and what can be measured about one.
"""

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import InvalidSequenceError, NucleicSequence

__all__ = ["InvalidSequenceError", "Nucleic", "NucleicSequence"]
