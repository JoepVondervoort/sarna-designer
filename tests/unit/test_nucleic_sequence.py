"""Tests for the validated nucleic-acid sequence primitive.

These are the lowest-level sequence operations every design rule is built on, so they are
tested harder than their simplicity suggests: a silent error here would be invisible in the
rule layer and would corrupt every candidate.
"""

import pytest

from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import InvalidSequenceError, NucleicSequence


class TestParsing:
    def test_parses_a_plain_dna_sequence(self) -> None:
        assert NucleicSequence.parse("ACGT").bases == "ACGT"

    def test_uppercases_lowercase_input(self) -> None:
        assert NucleicSequence.parse("acgt").bases == "ACGT"

    def test_strips_internal_and_surrounding_whitespace(self) -> None:
        # Pasted promoter sequences routinely arrive wrapped across lines.
        assert NucleicSequence.parse("  AC\n GT\t\r\n").bases == "ACGT"

    def test_rejects_an_unknown_base(self) -> None:
        with pytest.raises(InvalidSequenceError, match="position 3"):
            NucleicSequence.parse("ACNT")

    def test_names_every_offending_character_once(self) -> None:
        with pytest.raises(InvalidSequenceError) as excinfo:
            NucleicSequence.parse("ACXTZ")
        message = str(excinfo.value)
        assert "'X'" in message
        assert "'Z'" in message

    def test_rejects_u_in_a_dna_sequence(self) -> None:
        with pytest.raises(InvalidSequenceError, match="DNA"):
            NucleicSequence.parse("ACGU")

    def test_rejects_t_in_an_rna_sequence(self) -> None:
        with pytest.raises(InvalidSequenceError, match="RNA"):
            NucleicSequence.parse("ACGT", alphabet=Nucleic.RNA)

    def test_rejects_an_empty_sequence(self) -> None:
        # Length zero has no GC content and no 5' end; rejecting it removes a whole class of
        # division-by-zero and index errors from every rule downstream.
        with pytest.raises(InvalidSequenceError, match="empty"):
            NucleicSequence.parse("   ")

    def test_parses_rna_when_told_to(self) -> None:
        sequence = NucleicSequence.parse("ACGU", alphabet=Nucleic.RNA)
        assert sequence.bases == "ACGU"
        assert sequence.alphabet is Nucleic.RNA


class TestComplement:
    def test_complements_dna_without_reversing(self) -> None:
        # Defect B13: the legacy `complement` is documented as a reverse complement but is
        # a plain complement. Both exist here, separately named.
        assert NucleicSequence.parse("ACGT").complement().bases == "TGCA"

    def test_reverse_complements_dna(self) -> None:
        assert NucleicSequence.parse("ACGT").reverse_complement().bases == "ACGT"

    def test_reverse_complement_of_an_asymmetric_sequence(self) -> None:
        assert NucleicSequence.parse("AAGCTC").reverse_complement().bases == "GAGCTT"

    def test_reverse_complement_is_an_involution(self) -> None:
        sequence = NucleicSequence.parse("GATTACAGGCTAC")
        assert sequence.reverse_complement().reverse_complement() == sequence

    def test_reverse_complements_rna_using_u(self) -> None:
        sequence = NucleicSequence.parse("ACGU", alphabet=Nucleic.RNA)
        assert sequence.reverse_complement().bases == "ACGU"
        assert sequence.reverse_complement().alphabet is Nucleic.RNA

    def test_complement_preserves_the_alphabet(self) -> None:
        sequence = NucleicSequence.parse("AAAA", alphabet=Nucleic.RNA)
        assert sequence.complement().bases == "UUUU"


class TestTranscription:
    def test_dna_to_rna_replaces_t_with_u(self) -> None:
        rna = NucleicSequence.parse("ACGTTT").as_rna()
        assert rna.bases == "ACGUUU"
        assert rna.alphabet is Nucleic.RNA

    def test_rna_to_dna_replaces_u_with_t(self) -> None:
        dna = NucleicSequence.parse("ACGUUU", alphabet=Nucleic.RNA).as_dna()
        assert dna.bases == "ACGTTT"
        assert dna.alphabet is Nucleic.DNA

    def test_conversion_of_dna_to_dna_is_a_no_op(self) -> None:
        sequence = NucleicSequence.parse("ACGT")
        assert sequence.as_dna() is sequence

    def test_conversion_of_rna_to_rna_is_a_no_op(self) -> None:
        sequence = NucleicSequence.parse("ACGU", alphabet=Nucleic.RNA)
        assert sequence.as_rna() is sequence


class TestGcContent:
    @pytest.mark.parametrize(
        ("bases", "fraction"),
        [
            ("AAAA", 0.0),
            ("GGGG", 1.0),
            ("ACGT", 0.5),
            ("GCGCGCGCAA", 0.8),
        ],
    )
    def test_gc_fraction(self, bases: str, fraction: float) -> None:
        assert NucleicSequence.parse(bases).gc_fraction == pytest.approx(fraction)

    def test_gc_percent_is_the_fraction_scaled(self) -> None:
        assert NucleicSequence.parse("ACGT").gc_percent == pytest.approx(50.0)

    def test_gc_fraction_of_rna_counts_the_same_bases(self) -> None:
        sequence = NucleicSequence.parse("ACGU", alphabet=Nucleic.RNA)
        assert sequence.gc_fraction == pytest.approx(0.5)

    def test_a_nineteen_mer_at_the_wang_boundary(self) -> None:
        # Table 1 requires 40-60 %. A 19-mer cannot hit either bound exactly, which is why the
        # rule layer must compare fractions and not rounded percentages.
        sequence = NucleicSequence.parse("GGGGGGGGAAAAAAAAAAA")
        assert sequence.gc_fraction == pytest.approx(8 / 19)
        assert sequence.gc_percent == pytest.approx(42.105, abs=1e-3)


class TestPositions:
    def test_base_at_is_one_based(self) -> None:
        # Wang Table 1 speaks of "nucleotide 1" and "nucleotide 19". Indexing in the same
        # convention as the paper keeps the off-by-one out of every position rule.
        sequence = NucleicSequence.parse("ACGT")
        assert sequence.base_at(1) == "A"
        assert sequence.base_at(4) == "T"

    @pytest.mark.parametrize("position", [0, -1, 5])
    def test_base_at_rejects_a_position_outside_the_sequence(self, position: int) -> None:
        with pytest.raises(IndexError):
            NucleicSequence.parse("ACGT").base_at(position)

    def test_subsequence_is_one_based_and_inclusive(self) -> None:
        sequence = NucleicSequence.parse("ACGTACGT")
        assert sequence.subsequence(2, 4).bases == "CGT"

    def test_subsequence_rejects_a_range_running_past_the_end(self) -> None:
        with pytest.raises(IndexError, match="8"):
            NucleicSequence.parse("ACGTACGT").subsequence(6, 9)

    def test_windows_yields_every_window_of_a_given_length(self) -> None:
        windows = list(NucleicSequence.parse("ACGTA").windows(3))
        assert [offset for offset, _ in windows] == [0, 1, 2]
        assert [window.bases for _, window in windows] == ["ACG", "CGT", "GTA"]

    def test_windows_yields_nothing_when_the_sequence_is_too_short(self) -> None:
        assert list(NucleicSequence.parse("AC").windows(3)) == []

    def test_windows_rejects_a_non_positive_length(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            list(NucleicSequence.parse("ACGT").windows(0))


class TestHomopolymers:
    @pytest.mark.parametrize(
        ("bases", "longest"),
        [
            ("ACGT", 1),
            ("AACGT", 2),
            ("ACGTTTTT", 5),
            ("AAAAACGT", 5),
            ("GGGG", 4),
        ],
    )
    def test_longest_homopolymer_run(self, bases: str, longest: int) -> None:
        assert NucleicSequence.parse(bases).longest_homopolymer_run == longest

    def test_detects_a_run_of_five_as_the_legacy_scanner_does(self) -> None:
        # Defect B5: the legacy `checkconsecutive` docstring claims "four or more" but the loop
        # detects runs of five. Table 1 also says five, so behaviour is preserved and the
        # threshold lives in the rule layer, not here.
        assert NucleicSequence.parse("ACGAAAAACGT").longest_homopolymer_run == 5
        assert NucleicSequence.parse("ACGAAAACGT").longest_homopolymer_run == 4


class TestCpG:
    def test_counts_cpg_dinucleotides(self) -> None:
        assert NucleicSequence.parse("ACGCGT").cpg_count == 2

    def test_does_not_count_a_gc_step(self) -> None:
        # CpG is the ordered pair C then G. GC is a different dinucleotide entirely.
        assert NucleicSequence.parse("AGCT").cpg_count == 0

    def test_counts_overlapping_occurrences_once_each(self) -> None:
        assert NucleicSequence.parse("CGCGCG").cpg_count == 3

    def test_cpg_count_on_rna_is_rejected(self) -> None:
        # CpG is a property of the DNA target site. Counting it on the RNA duplex would
        # silently answer a different question.
        with pytest.raises(ValueError, match="DNA"):
            _ = NucleicSequence.parse("ACGU", alphabet=Nucleic.RNA).cpg_count


class TestTandemRepeats:
    @pytest.mark.parametrize(
        ("bases", "unit_length", "copies"),
        [
            ("ATATATAT", 2, 4),
            ("ACGACGACG", 3, 3),
            ("ACGT", 2, 1),
            ("CAGCAGCAGCAG", 3, 4),
        ],
    )
    def test_longest_tandem_repeat(self, bases: str, unit_length: int, copies: int) -> None:
        sequence = NucleicSequence.parse(bases)
        assert sequence.longest_tandem_repeat(unit_length) == copies

    def test_a_homopolymer_is_not_reported_as_a_dinucleotide_repeat(self) -> None:
        # AAAA is a run, not an AT-class simple repeat. Conflating them would double-count
        # one candidate against two different Table 1 criteria.
        assert NucleicSequence.parse("AAAAAAAA").longest_tandem_repeat(2) == 1

    def test_rejects_a_non_positive_unit_length(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            NucleicSequence.parse("ACGT").longest_tandem_repeat(0)


class TestValueSemantics:
    def test_equal_sequences_compare_equal(self) -> None:
        assert NucleicSequence.parse("ACGT") == NucleicSequence.parse("acgt")

    def test_the_same_bases_in_different_alphabets_are_not_equal(self) -> None:
        assert NucleicSequence.parse("ACGA") != NucleicSequence.parse(
            "ACGA", alphabet=Nucleic.RNA
        )

    def test_is_hashable(self) -> None:
        assert len({NucleicSequence.parse("ACGT"), NucleicSequence.parse("ACGT")}) == 1

    def test_is_immutable(self) -> None:
        sequence = NucleicSequence.parse("ACGT")
        with pytest.raises((AttributeError, TypeError)):
            sequence.bases = "TTTT"  # type: ignore[misc]

    def test_length(self) -> None:
        assert len(NucleicSequence.parse("ACGTA")) == 5

    def test_repr_shows_the_bases_and_alphabet(self) -> None:
        assert "ACGT" in repr(NucleicSequence.parse("ACGT"))

    def test_str_is_the_bases(self) -> None:
        assert str(NucleicSequence.parse("ACGT")) == "ACGT"
