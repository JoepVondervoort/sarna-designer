"""Grouping GTEx's sample groups into organs, and the symbols that stand for them.

The risk this module carries is not a crash. It is a symbol that means more than the data
says: an organ shown for a transcript that was measured at zero there, a sample group quietly
dropped because nobody mapped it, or a cultured cell group filed under an organ as though it
had been sampled from a body.
"""

import pytest

from sarna.app.expression_view import TOP_ORGANS, organ_summary
from sarna.expression.models import (
    ExpressionStatus,
    TissueMedian,
    TranscriptExpression,
)
from sarna.expression.organs import BY_ID, ORGANS, TISSUE_ORGAN, organ_for

CULTURED = {"Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes"}


def median(tissue_id: str, value: float, name: str = "") -> TissueMedian:
    """Build one median in a named sample group."""
    return TissueMedian(
        tissue_id=tissue_id, tissue_name=name or tissue_id.replace("_", " "),
        median=value, unit="TPM", sample_count=200,
    )


def measured(*medians: TissueMedian) -> TranscriptExpression:
    """Build a measured transcript record over the given medians."""
    return TranscriptExpression(
        transcript_id="ENST00000303395", status=ExpressionStatus.MEASURED, medians=medians,
    )


class TestEverySampleGroupIsAccountedFor:
    def test_the_map_covers_the_whole_gtex_v10_tissue_set(self) -> None:
        # 54 in gtex_v10, listed against the portal on 8 October 2026. A group nobody mapped
        # would be silently absent from every organ summary.
        assert len(TISSUE_ORGAN) == 54

    def test_every_mapped_organ_exists(self) -> None:
        assert set(TISSUE_ORGAN.values()) <= set(BY_ID)

    def test_every_declared_organ_is_used(self) -> None:
        # An organ nothing maps to is dead weight that will drift out of date unnoticed.
        assert {organ.id for organ in ORGANS} == set(TISSUE_ORGAN.values())

    def test_organ_identifiers_are_unique(self) -> None:
        assert len({organ.id for organ in ORGANS}) == len(ORGANS)

    def test_an_unmapped_group_returns_none_rather_than_guessing(self) -> None:
        assert organ_for("Pineal_Gland") is None


class TestTheGroupingItself:
    def test_all_thirteen_brain_groups_are_one_organ(self) -> None:
        brain = [t for t in TISSUE_ORGAN if t.startswith("Brain_")]
        assert len(brain) == 13
        assert {TISSUE_ORGAN[t] for t in brain} == {"nervous"}

    def test_the_tibial_nerve_is_nervous_system_despite_its_name(self) -> None:
        # A prefix rule on "Brain_" would miss it, which is why the map is written out.
        nerve = organ_for("Nerve_Tibial")
        assert nerve is not None
        assert nerve.id == "nervous"

    def test_the_salivary_gland_is_not_filed_as_muscle(self) -> None:
        # "Minor_Salivary_Gland" contains none of the muscle words, but a loose substring
        # rule over organ names is the kind of thing that would put it somewhere odd.
        salivary = organ_for("Minor_Salivary_Gland")
        assert salivary is not None
        assert salivary.id == "salivary"

    @pytest.mark.parametrize("tissue_id", sorted(CULTURED))
    def test_cultured_cells_are_not_a_body_site(self, tissue_id: str) -> None:
        organ = organ_for(tissue_id)
        assert organ is not None
        assert organ.id == "cultured"
        assert organ.is_body_site is False

    def test_every_other_organ_is_a_body_site(self) -> None:
        assert [o.id for o in ORGANS if not o.is_body_site] == ["cultured"]


class TestSymbolsOnlyWhereTheyAreHonest:
    def test_the_organs_with_symbols_are_the_ones_a_symbol_can_depict(self) -> None:
        # Deliberately a short list. Most organs have no Unicode emoji, and a lookalike object
        # would be ambiguous and undignified in a tool that reports measurements.
        assert {o.id for o in ORGANS if o.symbol} == {
            "nervous", "heart", "lung", "blood", "muscle", "thyroid", "kidney", "cultured",
        }

    def test_an_organ_without_a_symbol_carries_an_empty_string_not_a_placeholder(self) -> None:
        liver = BY_ID["liver"]
        assert liver.symbol == ""
        assert liver.name == "Liver"

    def test_no_two_organs_share_a_symbol(self) -> None:
        symbols = [o.symbol for o in ORGANS if o.symbol]
        assert len(set(symbols)) == len(symbols)

    def test_cultured_cells_are_marked_as_a_dish_not_an_organ(self) -> None:
        assert BY_ID["cultured"].symbol == "\N{PETRI DISH}"


class TestWhatASymbolClaims:
    def test_an_organ_appears_only_when_a_median_is_above_zero(self) -> None:
        entry = measured(median("Brain_Cortex", 1.6), median("Lung", 0.0))
        assert [o["organ_id"] for o in organ_summary(entry)] == ["nervous"]

    def test_a_measured_zero_produces_no_symbol_at_all(self) -> None:
        # The alternative, a greyed-out symbol, reads as "looked at and unclear". Zero is a
        # measurement and the number says so; the strip is for where it is made.
        entry = measured(median("Brain_Cortex", 0.0), median("Liver", 0.0))
        assert organ_summary(entry) == []

    def test_a_transcript_absent_from_the_release_has_no_organs(self) -> None:
        entry = TranscriptExpression(
            transcript_id="ENST00000713692", status=ExpressionStatus.NOT_IN_RELEASE,
        )
        assert organ_summary(entry) == []

    def test_an_organ_carries_the_group_the_value_came_from(self) -> None:
        # So a reader sees that "brain" came from frontal cortex, not from an average over
        # thirteen regions that nobody pooled.
        entry = measured(
            median("Brain_Cortex", 1.6, "Brain - Cortex"),
            median("Brain_Frontal_Cortex_BA9", 2.5, "Brain - Frontal Cortex (BA9)"),
        )
        [brain] = organ_summary(entry)
        assert brain["median"] == 2.5
        assert brain["tissue_name"] == "Brain - Frontal Cortex (BA9)"

    def test_an_organ_takes_its_highest_group_not_their_sum(self) -> None:
        entry = measured(median("Heart_Atrial_Appendage", 1.0),
                         median("Heart_Left_Ventricle", 3.0))
        [heart] = organ_summary(entry)
        assert heart["median"] == 3.0

    def test_organs_come_back_strongest_first(self) -> None:
        entry = measured(median("Lung", 0.2), median("Brain_Cortex", 2.5),
                         median("Liver", 1.0))
        assert [o["organ_id"] for o in organ_summary(entry)] == ["nervous", "liver", "lung"]

    def test_the_cultured_group_is_flagged_so_the_page_can_say_so(self) -> None:
        entry = measured(median("Cells_Cultured_fibroblasts", 4.0))
        [cultured] = organ_summary(entry)
        assert cultured["is_body_site"] is False
        assert cultured["organ_id"] == "cultured"


class TestWhatThePageIsGiven:
    def test_the_summary_is_capped_and_says_how_many_it_left_out(self) -> None:
        from sarna.app.expression_view import _summarise

        groups = ["Brain_Cortex", "Heart_Left_Ventricle", "Lung", "Liver", "Kidney_Cortex",
                  "Spleen", "Thyroid"]
        entry = measured(*[median(g, 10.0 - i) for i, g in enumerate(groups)])
        summary = _summarise(entry)
        assert len(summary["organs"]) == TOP_ORGANS
        assert summary["other_organs"] == len(groups) - TOP_ORGANS

    def test_nothing_is_left_out_when_there_are_few_organs(self) -> None:
        from sarna.app.expression_view import _summarise

        summary = _summarise(measured(median("Brain_Cortex", 2.5), median("Lung", 0.1)))
        assert len(summary["organs"]) == 2
        assert summary["other_organs"] == 0

    def test_each_organ_carries_its_name_for_the_ones_with_no_symbol(self) -> None:
        from sarna.app.expression_view import _summarise

        summary = _summarise(measured(median("Liver", 3.0)))
        assert summary["organs"][0]["symbol"] == ""
        assert summary["organs"][0]["name"] == "Liver"

    def test_every_organ_has_a_chip_label_short_enough_for_a_table_cell(self) -> None:
        # The text chip is the fallback on a machine with no colour emoji font, which is the
        # normal state of WSL and of most Linux containers. It has to fit a narrow cell, and
        # "Oesophagus, stomach and intestine" does not.
        for organ in ORGANS:
            assert organ.chip, organ.id
            assert len(organ.chip) <= 12, f"{organ.id}: {organ.chip!r} is too long for a chip"

    def test_the_chip_falls_back_to_the_full_name_when_it_is_short(self) -> None:
        assert BY_ID["liver"].chip == "Liver"
        assert BY_ID["nervous"].chip == "brain"
        assert BY_ID["digestive"].chip == "digestive"
