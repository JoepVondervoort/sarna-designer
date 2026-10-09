"""The validated nucleic-acid sequence value object.

Everything the engine computes about a candidate starts here, so this module is deliberately
conservative:

**It fails fast.** The legacy ``checksequence`` returns the *string* ``'failure'`` on bad input,
which is then happily sliced into 19-mers by the caller (defect B9 is the same failure mode in
``endbase``). Here a malformed sequence raises, naming the position and the character, because
the researcher pasting a promoter needs to know *where* the problem is.

**It is immutable and alphabet-tagged.** A ``NucleicSequence`` cannot be mutated in place and
always knows whether it is DNA or RNA, so a DNA-only question such as CpG content cannot be
asked of an RNA duplex.

**Positions are 1-based.** Wang et al. 2015 Table 1 speaks of "nucleotide 1", "nucleotide 18",
"nucleotide 19" and the "20-23" flank. :meth:`NucleicSequence.base_at` and
:meth:`NucleicSequence.subsequence` use that same convention so the position rules read like the
paper. Genomic *intervals* remain 0-based half-open (see :mod:`sarna.genome.coordinates`); the
two conventions never meet, because one indexes inside a sequence and the other locates a
sequence on a contig.
"""

from collections.abc import Iterator
from dataclasses import dataclass

from sarna.sequence.alphabet import Nucleic

_WHITESPACE = str.maketrans("", "", " \t\r\n\v\f")


class InvalidSequenceError(ValueError):
    """Raised when text cannot be read as a sequence in the requested alphabet.

    The message names the offending characters and their 1-based positions, so the error can be
    shown to a researcher as-is rather than translated by the UI layer.
    """


@dataclass(frozen=True, slots=True)
class NucleicSequence:
    """An immutable, validated run of nucleotides in a known alphabet.

    Construct with :meth:`parse` rather than directly: the constructor performs no normalisation,
    so ``bases`` is trusted to be already uppercase and already in ``alphabet``.
    """

    bases: str
    alphabet: Nucleic = Nucleic.DNA

    def __str__(self) -> str:
        """Return the bases, so a sequence can be written straight into a FASTA record."""
        return self.bases

    def __len__(self) -> int:
        """Return the number of nucleotides."""
        return len(self.bases)

    @classmethod
    def parse(cls, text: str, alphabet: Nucleic = Nucleic.DNA) -> "NucleicSequence":
        """Read ``text`` as a sequence, normalising case and discarding whitespace.

        Whitespace is stripped everywhere rather than only at the ends, because pasted and
        FASTA-derived sequences arrive wrapped across lines.

        Args:
            text: Raw sequence text. May be lower case and may contain line breaks.
            alphabet: The alphabet to validate against. DNA by default.

        Returns:
            The parsed sequence.

        Raises:
            InvalidSequenceError: If the text is empty once whitespace is removed, or contains
                any character outside ``alphabet``. ``U`` in DNA and ``T`` in RNA are errors
                rather than silent conversions, use :meth:`as_rna` or :meth:`as_dna` to
                convert deliberately.
        """
        bases = text.translate(_WHITESPACE).upper()
        if not bases:
            raise InvalidSequenceError(
                "the sequence is empty; paste at least one nucleotide"
            )
        permitted = set(alphabet.bases)
        offenders = [
            (index + 1, character)
            for index, character in enumerate(bases)
            if character not in permitted
        ]
        if offenders:
            raise InvalidSequenceError(cls._offender_message(offenders, alphabet))
        return cls(bases=bases, alphabet=alphabet)

    @staticmethod
    def _offender_message(
        offenders: list[tuple[int, str]], alphabet: Nucleic
    ) -> str:
        """Build the parse-error message, reporting each bad character once."""
        reported: dict[str, int] = {}
        for position, character in offenders:
            reported.setdefault(character, position)
        detail = ", ".join(
            f"{character!r} at position {position}"
            for character, position in reported.items()
        )
        expected = "".join(alphabet.bases)
        hint = ""
        wrong_alphabet = alphabet.other.thymine_analogue
        if wrong_alphabet in reported:
            hint = (
                f" {wrong_alphabet!r} belongs to the {alphabet.other} alphabet; convert "
                f"explicitly rather than relying on this parser to do it."
            )
        return (
            f"not a valid {alphabet} sequence: {detail}. Expected only {expected} "
            f"(case and whitespace are ignored).{hint}"
        )

    def as_rna(self) -> "NucleicSequence":
        """Return this sequence as RNA, replacing T with U. A no-op if already RNA."""
        if self.alphabet is Nucleic.RNA:
            return self
        return NucleicSequence(bases=self.bases.replace("T", "U"), alphabet=Nucleic.RNA)

    def as_dna(self) -> "NucleicSequence":
        """Return this sequence as DNA, replacing U with T. A no-op if already DNA."""
        if self.alphabet is Nucleic.DNA:
            return self
        return NucleicSequence(bases=self.bases.replace("U", "T"), alphabet=Nucleic.DNA)

    def complement(self) -> "NucleicSequence":
        """Return the base-by-base complement, **without** reversing.

        Rarely what a caller wants on its own; :meth:`reverse_complement` is the biologically
        meaningful operation. Both exist separately because the legacy script conflated them in
        its naming (defect B13).
        """
        mapping = self.alphabet.complement
        return NucleicSequence(
            bases="".join(mapping[base] for base in self.bases),
            alphabet=self.alphabet,
        )

    def reverse_complement(self) -> "NucleicSequence":
        """Return the reverse complement: the opposing strand read 5' to 3'.

        This is the operation that turns a plus-strand target site into the sequence a
        minus-strand reader sees, and the antisense strand of a duplex into its sense strand.
        """
        mapping = self.alphabet.complement
        return NucleicSequence(
            bases="".join(mapping[base] for base in reversed(self.bases)),
            alphabet=self.alphabet,
        )

    @property
    def gc_fraction(self) -> float:
        """Return the proportion of G or C nucleotides, between 0 and 1.

        Reported as a fraction rather than a rounded percentage so that a rule comparing against
        a 40-60 % window is not decided by display rounding: no 19-mer has a GC content that
        lands exactly on either bound.
        """
        return (self.bases.count("G") + self.bases.count("C")) / len(self.bases)

    @property
    def gc_percent(self) -> float:
        """Return the GC content as a percentage, for display."""
        return self.gc_fraction * 100.0

    def base_at(self, position: int) -> str:
        """Return the nucleotide at a **1-based** ``position``.

        Args:
            position: 1-based position, counted from the 5' end.

        Returns:
            A single-character base.

        Raises:
            IndexError: If the position is outside the sequence. Note that ``base_at(-1)`` is an
                error rather than the last base: negative indexing would silently return a
                plausible answer for an off-by-one bug in a rule.
        """
        if not 1 <= position <= len(self.bases):
            raise IndexError(
                f"position {position} is outside this {len(self.bases)} nt sequence; "
                f"positions are 1-based, so valid values are 1 to {len(self.bases)}"
            )
        return self.bases[position - 1]

    def subsequence(self, start: int, end: int) -> "NucleicSequence":
        """Return the subsequence from ``start`` to ``end``, 1-based and **inclusive**.

        Args:
            start: 1-based first position to include.
            end: 1-based last position to include.

        Returns:
            The subsequence, in the same alphabet.

        Raises:
            IndexError: If the range is empty or runs outside the sequence.
        """
        length = len(self.bases)
        if not 1 <= start <= end <= length:
            raise IndexError(
                f"cannot take positions {start}-{end} of a {length} nt sequence; the range "
                f"must satisfy 1 <= start <= end <= {length}"
            )
        return NucleicSequence(bases=self.bases[start - 1 : end], alphabet=self.alphabet)

    def windows(self, length: int) -> Iterator[tuple[int, "NucleicSequence"]]:
        """Yield every window of ``length`` nucleotides, with its 0-based offset.

        The offset is 0-based because it is an *offset into* this sequence, which is what
        genomic-coordinate arithmetic adds to a start position. A sequence shorter than
        ``length`` yields nothing rather than raising, since "no candidates fit in this promoter"
        is a legitimate outcome the caller reports to the user.

        Args:
            length: Window length in nucleotides, e.g. 19 for a Wang 2015 target.

        Yields:
            ``(offset, window)`` pairs, 5' to 3'.

        Raises:
            ValueError: If ``length`` is not positive.
        """
        if length <= 0:
            raise ValueError(f"window length must be positive, got {length}")
        for offset in range(len(self.bases) - length + 1):
            yield (
                offset,
                NucleicSequence(
                    bases=self.bases[offset : offset + length], alphabet=self.alphabet
                ),
            )

    @property
    def longest_homopolymer_run(self) -> int:
        """Return the length of the longest run of a single repeated nucleotide.

        The *threshold* is not applied here. Table 1 excludes runs of five or more, but which
        number counts as too long is a configurable rule parameter, so this returns the measured
        value and the rule layer compares it (defect B5 was precisely a hardcoded threshold
        disagreeing with its own documentation).
        """
        longest = 1
        current = 1
        for previous, base in zip(self.bases, self.bases[1:], strict=False):
            current = current + 1 if base == previous else 1
            longest = max(longest, current)
        return longest

    @property
    def cpg_count(self) -> int:
        """Return the number of CpG dinucleotides: an ordered C followed by a G.

        Overlapping occurrences each count, so ``CGCGCG`` has three. ``GC`` is a different
        dinucleotide and is not counted.

        Raises:
            ValueError: If this sequence is RNA. CpG methylation is a property of the DNA target
                site, not of the delivered duplex, so answering for RNA would answer a different
                question than the caller asked.
        """
        if self.alphabet is not Nucleic.DNA:
            raise ValueError(
                "CpG content is defined on the DNA target site, but this sequence is RNA; "
                "call .as_dna() on the target site if that is what you meant"
            )
        return sum(
            1
            for first, second in zip(self.bases, self.bases[1:], strict=False)
            if first == "C" and second == "G"
        )

    def longest_tandem_repeat(self, unit_length: int) -> int:
        """Return the most consecutive copies of any repeating unit of ``unit_length``.

        Used for the Table 1 simple-repeat criterion, which excludes di- and tri-nucleotide
        repeats such as ``ATATAT`` or ``CAGCAGCAG``. A unit of identical bases is **not**
        counted, so ``AAAAAA`` reports one copy rather than three: that is a homopolymer run,
        already measured by :attr:`longest_homopolymer_run`, and counting it twice would
        penalise one candidate against two separate criteria for the same feature.

        Args:
            unit_length: Length of the repeating unit, e.g. 2 for dinucleotide repeats.

        Returns:
            The largest number of consecutive identical units found, at least 1.

        Raises:
            ValueError: If ``unit_length`` is not positive.
        """
        if unit_length <= 0:
            raise ValueError(f"unit length must be positive, got {unit_length}")
        longest = 1
        for start in range(len(self.bases) - unit_length + 1):
            unit = self.bases[start : start + unit_length]
            if len(set(unit)) == 1:
                continue
            copies = 1
            while self.bases.startswith(unit, start + copies * unit_length):
                copies += 1
            longest = max(longest, copies)
        return longest
