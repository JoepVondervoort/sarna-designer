"""Tests for the chapter 4 empirical dataset.

These pin the dataset's *contents* as well as its mechanics, because the file was transcribed
from a PDF appendix. A silent transcription error would propagate into every validation claim
made from it, and would look like a scientific finding rather than a parsing bug.
"""

import pytest

from sarna.sequence.alphabet import Nucleic
from sarna.validation.dataset import (
    Activity,
    ValidationDataset,
    load_chapter4_dataset,
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
