"""Tests for the chapter 4 empirical dataset.

These pin the dataset's *contents* as well as its mechanics, because the file was transcribed
from a PDF appendix. A silent transcription error would propagate into every validation claim
made from it, and would look like a scientific finding rather than a parsing bug.
"""

import pytest

from sarna.genome.coordinates import Assembly, Strand
from sarna.sequence.alphabet import Nucleic
from sarna.validation.dataset import (
    CONTEXT_CAVEATS,
    Activity,
    ValidationDataset,
    load_chapter4_dataset,
    load_chapter4_locations,
    target_dna_sequences,
)

DATASET = load_chapter4_dataset()


class TestDatasetContents:
    def test_the_chapter_reports_125_sequences(self) -> None:
        assert len(DATASET) == 125

    def test_fourteen_activated_their_target(self) -> None:
        # The chapter's own count. If transcription dropped or duplicated a row this breaks.
        assert len(DATASET.activating) == 14
        assert len(DATASET.not_activating) == 111

    def test_every_sequence_is_19_nucleotides(self) -> None:
        # Wang's framework is written for 19-mers; a 18 or 20 here means a parse error.
        assert all(len(s.target) == 19 for s in DATASET)

    def test_names_are_unique(self) -> None:
        names = [s.name for s in DATASET]
        assert len(set(names)) == len(names)

    def test_seventeen_genes_are_represented(self) -> None:
        assert len(DATASET.genes) == 17
        assert "SCN1A" in DATASET.genes
        assert "CHD2" in DATASET.genes

    def test_sequences_are_held_as_rna(self) -> None:
        assert all(s.target.alphabet is Nucleic.RNA for s in DATASET)

    def test_a_known_published_activator_is_present_verbatim(self) -> None:
        # Li et al. 2006 dsP21-322, the field's canonical positive control.
        p21 = DATASET.by_name("P21-322")
        assert p21.target.bases == "CCAACUCAUUCUCCAAGUA"
        assert p21.is_activating

    def test_the_published_control_reverse_complements_to_its_known_antisense(self) -> None:
        # dsP21-322's antisense strand as published. This is the check that the stored strand is
        # the sense strand: reading the stored sequence as a guide would give a different duplex.
        assert DATASET.by_name("P21-322").guide.bases == "UACUUGGAGAAUGAGUUGG"

    def test_the_scn1a_activators_from_earlier_chapters_are_present(self) -> None:
        for name in ("SCN1A-10", "SCN1A-12", "SCN1A-14"):
            assert DATASET.by_name(name).is_activating

    def test_a_negative_is_labelled_not_activating(self) -> None:
        assert DATASET.by_name("SCN1A-1").activity is Activity.NOT_ACTIVATING


class TestStrandConvention:
    """The stored strand decides what every position-based finding means (ADR 0011)."""

    def test_most_sequences_carry_a_at_position_19(self) -> None:
        # A Table 1 criterion written about the target. 113 of 125; the evidence that these are
        # target strands rather than guides.
        with_a = sum(1 for s in DATASET if s.target.base_at(19) == "A")
        assert with_a == 113

    def test_nearly_all_carry_a_or_u_at_position_18(self) -> None:
        weak = sum(1 for s in DATASET if s.target.base_at(18) in "AU")
        assert weak == 123

    def test_read_as_guides_the_same_criteria_would_not_hold(self) -> None:
        # If these were guide strands, Table 1's target rules would appear reverse-complemented,
        # putting U at position 1. Only 15 of 125 do, so the guide reading is not tenable.
        starting_u = sum(1 for s in DATASET if s.target.base_at(1) == "U")
        assert starting_u == 15


class TestSequenceAccessors:
    def test_the_target_is_available_as_dna_for_the_design_rules(self) -> None:
        p21 = DATASET.by_name("P21-322")
        assert p21.target_dna.bases == "CCAACTCATTCTCCAAGTA"
        assert p21.target_dna.alphabet is Nucleic.DNA

    def test_the_guide_is_the_reverse_complement_of_the_target(self) -> None:
        for sequence in DATASET:
            assert sequence.guide.bases == sequence.target.reverse_complement().bases

    def test_target_dna_sequences_converts_a_whole_selection(self) -> None:
        converted = target_dna_sequences(DATASET.activating)
        assert len(converted) == 14
        assert all(s.alphabet is Nucleic.DNA for s in converted)

    def test_an_unknown_name_raises_and_names_the_dataset(self) -> None:
        with pytest.raises(KeyError, match="vondervoort-2026-chapter4"):
            DATASET.by_name("NO-SUCH-1")


class TestProvenanceAndCaveats:
    def test_the_dataset_cites_its_source(self) -> None:
        assert "van de Vondervoort 2026" in str(DATASET.provenance)
        assert "chapter 4" in str(DATASET.provenance)

    def test_the_range_restriction_caveat_is_carried(self) -> None:
        # The most important limitation: the set cannot measure the mandatory gate, because
        # every sequence in it was built to pass that gate. An export must say so.
        joined = " ".join(DATASET.caveats).lower()
        assert "range-restricted" in joined
        assert "cannot measure" in joined

    def test_the_cell_type_caveat_is_carried(self) -> None:
        joined = " ".join(DATASET.caveats).lower()
        assert "af22" in joined

    def test_the_duplicate_name_and_sequence_issues_are_recorded(self) -> None:
        joined = " ".join(DATASET.caveats)
        assert "GAD1-1" in joined
        assert "GAD1-5b" in joined

    def test_the_duplicate_negative_sequence_is_really_there(self) -> None:
        # Asserted so the caveat and the data cannot drift apart.
        assert DATASET.by_name("GAD1-5b").target == DATASET.by_name("GAD1-6b").target

    def test_the_negatives_contain_110_distinct_sequences(self) -> None:
        distinct = {s.target.bases for s in DATASET.not_activating}
        assert len(distinct) == 110

    def test_the_export_states_its_own_limits(self) -> None:
        exported = DATASET.to_dict()
        assert exported["counts"] == {
            "total": 125,
            "activating": 14,
            "not_activating": 111,
        }
        assert exported["caveats"]
        assert exported["provenance"]

    def test_the_export_is_json_serialisable(self) -> None:
        import json

        json.dumps(DATASET.to_dict())


class TestImmutability:
    def test_a_sequence_cannot_be_mutated(self) -> None:
        with pytest.raises((AttributeError, TypeError)):
            DATASET.sequences[0].name = "other"  # type: ignore[misc]

    def test_a_dataset_cannot_be_mutated(self) -> None:
        with pytest.raises((AttributeError, TypeError)):
            DATASET.id = "other"  # type: ignore[misc]

    def test_loading_twice_yields_equal_datasets(self) -> None:
        assert load_chapter4_dataset().sequences == DATASET.sequences

    def test_the_loader_returns_the_documented_type(self) -> None:
        assert isinstance(DATASET, ValidationDataset)


LOCATIONS = load_chapter4_locations()


class TestRecoveredCoordinates:
    """The addition that makes the location rules testable against measured activity."""

    def test_123_of_125_designs_are_located(self) -> None:
        # CHD8-9 and CHD8-10 each match 12 places in their locus, so they are absent rather
        # than present with a guess.
        assert len(LOCATIONS) == 123
        assert "CHD8-9" not in LOCATIONS
        assert "CHD8-10" not in LOCATIONS

    def test_every_located_design_is_in_the_dataset(self) -> None:
        for name in LOCATIONS:
            assert DATASET.by_name(name)

    def test_coordinates_are_on_grch38(self) -> None:
        assert all(loc.target.assembly is Assembly.GRCH38 for loc in LOCATIONS.values())

    def test_each_target_interval_is_19_bases(self) -> None:
        assert all(loc.target.length == 19 for loc in LOCATIONS.values())

    def test_the_stored_context_contains_its_target_at_the_stated_offset(self) -> None:
        # The check that the minus-strand reverse-complement and the offset agree. Without it a
        # mirrored context would look plausible and place every flank on the wrong side.
        for name, loc in LOCATIONS.items():
            expected = DATASET.by_name(name).target_dna.bases
            start = loc.context_offset + 1
            assert loc.context.subsequence(start, start + 18).bases == expected

    def test_every_design_has_the_flank_the_preferences_need(self) -> None:
        # Positions 20-23 lie in the promoter beyond the target, so the context must reach them.
        for loc in LOCATIONS.values():
            assert len(loc.context) >= loc.context_offset + 19 + 4


@pytest.mark.strand
class TestTheRecoveredPromoters:
    def test_every_promoter_reproduces_its_target_interval(self) -> None:
        # The round trip: the promoter's own strand arithmetic, run forwards, must land back on
        # the interval the coordinates record. This is the test that catches a sign error.
        for loc in LOCATIONS.values():
            promoter = loc.promoter()
            assert promoter.genomic_interval(loc.context_offset, 19) == loc.target

    def test_both_strands_are_represented(self) -> None:
        strands = {loc.target.strand for loc in LOCATIONS.values()}
        assert strands == {Strand.PLUS, Strand.MINUS}

    def test_the_tss_span_matches_the_recorded_offset(self) -> None:
        for loc in LOCATIONS.values():
            span = loc.promoter().tss_span(loc.context_offset, 19)
            assert span is not None
            assert span[0] == loc.tss_position

    def test_a_minus_strand_design_is_placed_correctly(self) -> None:
        # SCN1A is on the minus strand, where upstream means higher coordinates.
        loc = LOCATIONS["SCN1A-11"]
        assert loc.target.strand is Strand.MINUS
        assert loc.nearest_offset == -6010


class TestExternalValidationOfTheCoordinates:
    """Two published saRNAs are named for their TSS-relative position. That is the control.

    Nothing in the chain — Ensembl lookup, MANE Select choice, strand handling, locus fetch,
    sequence search, offset arithmetic — can be meaningfully wrong and still land within a few
    bases of a position published twenty years ago.
    """

    def test_dsp21_322_recovers_near_minus_322(self) -> None:
        location = LOCATIONS["P21-322"]
        assert location.mane_offset is not None
        assert abs(location.mane_offset - (-322)) < 20

    def test_vhl_821_recovers_near_minus_821(self) -> None:
        location = LOCATIONS["VHL-821"]
        assert location.mane_offset is not None
        assert abs(location.mane_offset - (-821)) < 10

    def test_the_caveats_record_this_validation(self) -> None:
        joined = " ".join(CONTEXT_CAVEATS)
        assert "VHL-821" in joined
        assert "annotation drift" in joined

    def test_the_caveats_record_what_the_context_cannot_support(self) -> None:
        joined = " ".join(CONTEXT_CAVEATS)
        assert "island" in joined
        assert "67 annotated TSSs" in joined


class TestTheNearestTssIsAWeakAnchor:
    def test_genes_with_many_tss_are_flagged_by_the_count(self) -> None:
        # TSC2 has 67 annotated TSSs, so "nearest" is close to meaningless there.
        crowded = max(LOCATIONS.values(), key=lambda loc: loc.distinct_tss)
        assert crowded.distinct_tss > 40

    def test_both_anchors_are_available(self) -> None:
        location = LOCATIONS["TSC2-8"]
        assert location.mane_offset is not None
        assert location.mane_offset != location.nearest_offset

    def test_the_mane_offset_is_preferred_when_present(self) -> None:
        location = LOCATIONS["TSC2-8"]
        assert location.offset == location.mane_offset
