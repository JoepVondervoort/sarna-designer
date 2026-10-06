"""Coordinates carry their assembly, and refuse to be mixed across assemblies.

These tests exist because of defect B15: the legacy pipeline combined GRCh38 TSS coordinates
with GRCh37 signal tracks, silently. ADR 0004 makes that unrepresentable rather than merely
discouraged, so the guarantee is asserted here first.
"""

import pytest

from sarna.genome.coordinates import (
    Assembly,
    AssemblyMismatchError,
    GenomicInterval,
    Strand,
)


class TestAssembly:
    def test_renders_as_the_conventional_assembly_name(self):
        assert str(Assembly.GRCH38) == "GRCh38"
        assert str(Assembly.GRCH37) == "GRCh37"

    def test_known_contig_lengths_distinguish_them(self):
        # The exact values that identified the AF22 bigWigs as GRCh37 (inventory section 12).
        assert Assembly.GRCH37.contig_length("1") == 249_250_621
        assert Assembly.GRCH38.contig_length("1") == 248_956_422
        assert Assembly.GRCH37.contig_length("X") == 155_270_560
        assert Assembly.GRCH38.contig_length("X") == 156_040_895

    def test_unknown_contig_is_rejected(self):
        with pytest.raises(ValueError, match="unknown contig"):
            Assembly.GRCH38.contig_length("chr1")


class TestGenomicIntervalConstruction:
    def test_builds_a_plain_interval(self):
        iv = GenomicInterval(Assembly.GRCH38, "7", 1_000, 1_019, Strand.PLUS)
        assert iv.length == 19
        assert iv.assembly is Assembly.GRCH38

    def test_strand_defaults_to_plus(self):
        assert GenomicInterval(Assembly.GRCH38, "7", 10, 29).strand is Strand.PLUS

    def test_coordinates_are_zero_based_half_open(self):
        # [0, 19) is the first 19 bases, not [1, 19].
        assert GenomicInterval(Assembly.GRCH38, "1", 0, 19).length == 19

    def test_rejects_empty_interval(self):
        with pytest.raises(ValueError, match="end must be greater than start"):
            GenomicInterval(Assembly.GRCH38, "1", 500, 500)

    def test_rejects_reversed_interval(self):
        with pytest.raises(ValueError, match="end must be greater than start"):
            GenomicInterval(Assembly.GRCH38, "1", 500, 499)

    def test_rejects_negative_start(self):
        with pytest.raises(ValueError, match="start must not be negative"):
            GenomicInterval(Assembly.GRCH38, "1", -1, 19)

    def test_rejects_interval_running_past_the_end_of_the_contig(self):
        length = Assembly.GRCH38.contig_length("1")
        with pytest.raises(ValueError, match="beyond the end of contig"):
            GenomicInterval(Assembly.GRCH38, "1", length - 5, length + 5)

    def test_rejects_ucsc_style_contig_names(self):
        # chr-prefixing is a bigWig file-format detail and must not reach the domain model.
        with pytest.raises(ValueError, match="without the 'chr' prefix"):
            GenomicInterval(Assembly.GRCH38, "chr1", 0, 19)

    def test_rejects_unknown_contig(self):
        with pytest.raises(ValueError, match="unknown contig"):
            GenomicInterval(Assembly.GRCH38, "banana", 0, 19)

    def test_is_immutable(self):
        iv = GenomicInterval(Assembly.GRCH38, "1", 0, 19)
        with pytest.raises((AttributeError, TypeError)):
            iv.start = 5  # type: ignore[misc]


class TestAssemblySafety:
    def test_intervals_on_different_assemblies_are_never_equal(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 100, 119)
        b = GenomicInterval(Assembly.GRCH37, "1", 100, 119)
        assert a != b

    def test_overlap_across_assemblies_raises_rather_than_answering(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 100, 119)
        b = GenomicInterval(Assembly.GRCH37, "1", 100, 119)
        with pytest.raises(AssemblyMismatchError):
            a.overlaps(b)

    def test_overlap_within_one_assembly_is_computed(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 100, 120)
        assert a.overlaps(GenomicInterval(Assembly.GRCH38, "1", 110, 130))
        assert not a.overlaps(GenomicInterval(Assembly.GRCH38, "1", 120, 130))
        assert not a.overlaps(GenomicInterval(Assembly.GRCH38, "2", 110, 130))

    def test_assembly_mismatch_error_names_both_assemblies(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 100, 119)
        b = GenomicInterval(Assembly.GRCH37, "1", 100, 119)
        with pytest.raises(AssemblyMismatchError) as excinfo:
            a.overlaps(b)
        message = str(excinfo.value)
        assert "GRCh38" in message
        assert "GRCh37" in message


class TestStrandAwareOperations:
    @pytest.mark.strand
    def test_five_prime_position_depends_on_strand(self):
        plus = GenomicInterval(Assembly.GRCH38, "1", 100, 119, Strand.PLUS)
        minus = GenomicInterval(Assembly.GRCH38, "1", 100, 119, Strand.MINUS)
        assert plus.five_prime == 100
        # The 5' end of a minus-strand feature is its highest coordinate. B1 is exactly the
        # class of bug that comes from assuming otherwise.
        assert minus.five_prime == 118

    @pytest.mark.strand
    def test_signed_distance_from_a_reference_respects_strand(self):
        # Upstream of a plus-strand TSS is negative; upstream of a minus-strand TSS is at a
        # *higher* coordinate and must still read as negative.
        tss_plus = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001, Strand.PLUS)
        upstream_plus = GenomicInterval(Assembly.GRCH38, "1", 500, 519, Strand.PLUS)
        assert upstream_plus.signed_distance_to(tss_plus) == -500

        # The exact mirror image: the minus-strand target's 5' end (its *highest* coordinate,
        # 1500) sits 500 bp upstream of the TSS, so it must report the same -500.
        tss_minus = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001, Strand.MINUS)
        upstream_minus = GenomicInterval(Assembly.GRCH38, "1", 1_482, 1_501, Strand.MINUS)
        assert upstream_minus.five_prime == 1_500
        assert upstream_minus.signed_distance_to(tss_minus) == -500

    @pytest.mark.strand
    def test_downstream_is_positive_on_both_strands(self):
        tss_plus = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001, Strand.PLUS)
        assert (
            GenomicInterval(Assembly.GRCH38, "1", 1_200, 1_219, Strand.PLUS).signed_distance_to(
                tss_plus
            )
            == 200
        )
        tss_minus = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001, Strand.MINUS)
        assert (
            GenomicInterval(Assembly.GRCH38, "1", 782, 801, Strand.MINUS).signed_distance_to(
                tss_minus
            )
            == 200
        )

    @pytest.mark.strand
    def test_signed_distance_across_assemblies_raises(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 500, 519)
        tss = GenomicInterval(Assembly.GRCH37, "1", 1_000, 1_001)
        with pytest.raises(AssemblyMismatchError):
            a.signed_distance_to(tss)

    @pytest.mark.strand
    def test_signed_distance_across_contigs_raises(self):
        a = GenomicInterval(Assembly.GRCH38, "2", 500, 519)
        tss = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001)
        with pytest.raises(ValueError, match="between contigs"):
            a.signed_distance_to(tss)

    @pytest.mark.strand
    def test_signed_distance_requires_matching_strand(self):
        a = GenomicInterval(Assembly.GRCH38, "1", 500, 519, Strand.PLUS)
        tss = GenomicInterval(Assembly.GRCH38, "1", 1_000, 1_001, Strand.MINUS)
        with pytest.raises(ValueError, match="strand"):
            a.signed_distance_to(tss)
