"""Tests for the order-ready duplex: what a synthesis vendor is actually asked to make.

A candidate carries a 19 nt RNA core. What gets ordered is two 21-mers, each the core plus a
2 nt dTdT 3' overhang. Chapter 2's supplementary material settles that the overhang is
deoxythymidine written ``[dT][dT]``, not a ribo-``UU``. `docs/scientific/` records the finding in
`scn1a-design-coordinates.md` and `GLOSSARY.md` the vocabulary. These tests pin the chemistry
against that source, so a change to it breaks a test that names the paper.

The expected strands below are hand-derived, not taken from the implementation:

    target DNA  GGGTAATTTTGACAGGTCA
    sense RNA   GGGUAAUUUUGACAGGUCA   (transcribed, T -> U)
    complement  CCCAUUAAAACUGUCCAGU
    reversed    UGACCUGUCAAAAUUACCC   (the antisense core)
"""

import pytest

from sarna.design.order import (
    CSV_COLUMNS,
    ORDERED_LENGTH_NT,
    OVERHANG_ANNOTATED,
    OVERHANG_LENGTH_NT,
    SUPPORTED_CORE_LENGTH_NT,
    OrderOutcome,
    csv_row,
    order_for,
    order_from_record,
)
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.promoter import Promoter
from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence

LEFT = "GCTTAACGGATTCAAGTTACCGATTACAGGCTAACGTTAGCCATTAAGCGTTACA"
RIGHT = "ATCGGTTACAGCTTAAGGCTATTACGGCATTAACGGTTACAGCTTAAGGCATTAA"

TARGET = "GGGTAATTTTGACAGGTCA"
"""All three tables agree this loads the designed antisense strand."""

SENSE_CORE = "GGGUAAUUUUGACAGGUCA"
ANTISENSE_CORE = "UGACCUGUCAAAAUUACCC"
BASES = LEFT + TARGET + RIGHT

DIVIDED = "TTTCCTCATGCAATTCAAA"
"""Divides the tables at the default end length: RNA says antisense, both DNA tables say sense."""

DIVIDED_BASES = LEFT + DIVIDED + RIGHT

SENSE_PREDICTED = "TTCGCTCTATTGACTACGA"
"""All three tables agree this loads the *sense* strand, against the design intent."""

SENSE_BASES = LEFT + SENSE_PREDICTED + RIGHT


def candidate_for(bases: str = BASES, target: str = TARGET):
    """Return the generated candidate whose target is ``target``."""
    result = design_sarnas(
        Promoter(sequence=NucleicSequence.parse(bases)), WANG_2015.resolve()
    )
    return next(c for c in result.generated if c.target.bases == target)


def check(order, check_id: str):  # noqa: ANN001
    """Return the named sanity check from ``order``."""
    return next(c for c in order.checks if c.id == check_id)


class TestTheChemistryMatchesChapterTwo:
    """The overhang is deoxythymidine on both strands, and the core is untouched."""

    def test_both_strands_carry_the_overhang(self) -> None:
        # Chapter 2's supplementary material: dTdT on both strands, not one. The expected
        # strings are written out rather than built from OVERHANG_ANNOTATED, so that changing
        # that constant breaks this test instead of being carried into its expectation.
        order = order_for(candidate_for())
        assert order.sense.annotated == "GGGUAAUUUUGACAGGUCA[dT][dT]"
        assert order.antisense.annotated == "UGACCUGUCAAAAUUACCC[dT][dT]"

    def test_the_overhang_is_deoxythymidine_rather_than_ribo_uracil(self) -> None:
        # The distinction chapter 2 settles. A ribo-UU overhang would read "UU" here, and the
        # core's own uracils must not be mistaken for it.
        assert OVERHANG_ANNOTATED == "[dT][dT]"
        assert "U" not in OVERHANG_ANNOTATED
        order = order_for(candidate_for())
        assert order.sense.overhang_bases == "TT"
        assert order.sense.overhang_alphabet is Nucleic.DNA

    def test_what_is_ordered_is_two_twenty_one_mers(self) -> None:
        order = order_for(candidate_for())
        assert SUPPORTED_CORE_LENGTH_NT == 19
        assert OVERHANG_LENGTH_NT == 2
        assert ORDERED_LENGTH_NT == 21
        assert order.sense.length_nt == 21
        assert order.antisense.length_nt == 21

    def test_the_core_stays_rna_and_the_overhang_is_not_part_of_it(self) -> None:
        # The glossary: the overhang "is not part of the Target and is never scanned or scored".
        order = order_for(candidate_for())
        assert order.sense.core.bases == SENSE_CORE
        assert order.sense.core.alphabet is Nucleic.RNA
        assert "T" not in order.sense.core.bases
        assert len(order.sense.core) == 19

    def test_the_target_site_is_reported_as_dna_and_unchanged(self) -> None:
        # Invariant: target sites are DNA, duplex strands are RNA. Ordering changes neither.
        order = order_for(candidate_for())
        assert order.target_site.bases == TARGET
        assert order.target_site.alphabet is Nucleic.DNA

    def test_the_two_cores_anneal(self) -> None:
        order = order_for(candidate_for())
        assert order.antisense.core.bases == order.sense.core.reverse_complement().bases


class TestTheSanityReport:
    """Constraint 15: the export carries a sanity report, whose findings never exclude."""

    def test_a_clean_duplex_reports_no_concern_on_its_chemistry(self) -> None:
        order = order_for(candidate_for())
        assert check(order, "core-length").outcome is OrderOutcome.OK
        assert check(order, "strand-complementarity").outcome is OrderOutcome.OK

    def test_off_target_screening_is_always_reported_as_not_performed(self) -> None:
        # Missing data is not a negative assay. No off-target analysis exists in this project,
        # so the report must say so rather than omit the row and imply a clean screen.
        order = order_for(candidate_for())
        screening = check(order, "off-target-screening")
        assert screening.outcome is OrderOutcome.NOT_PERFORMED
        assert "not" in screening.detail.lower()

    def test_a_concern_never_withholds_the_sequences(self) -> None:
        # Findings report rather than exclude: the researcher still gets something to order.
        order = order_for(candidate_for())
        assert order.sense.annotated
        assert order.antisense.annotated

    def test_the_report_names_the_model_behind_the_strand_prediction(self) -> None:
        order = order_for(candidate_for())
        assert "rna-xia-1998" in check(order, "guide-matches-intent").detail

    def test_disagreeing_tables_are_a_concern_rather_than_silence(self) -> None:
        # 13.72% of inputs divide the tables (docs/scientific/guide-prediction-divergence.md).
        # A researcher about to spend money should see that, not an averaged verdict.
        order = order_for(candidate_for(bases=DIVIDED_BASES, target=DIVIDED))
        agreement = check(order, "model-agreement")
        assert agreement.outcome is OrderOutcome.CONCERN
        assert "disagree" in agreement.detail.lower()
        assert "13.72%" in agreement.detail

    def test_the_same_duplex_can_match_intent_under_one_model_and_not_another(self) -> None:
        # The point of reporting every table rather than averaging them (ADR 0009). On this
        # target the RNA table expects the designed antisense strand and both DNA tables expect
        # the sense strand, so the intent verdict legitimately depends on which is named.
        candidate = candidate_for(bases=DIVIDED_BASES, target=DIVIDED)
        assert check(order_for(candidate, model="rna-xia-1998"), "guide-matches-intent").outcome \
            is OrderOutcome.OK
        for dna_table in ("dna-santalucia-1998", "legacy-dna-kj"):
            concern = check(order_for(candidate, model=dna_table), "guide-matches-intent")
            assert concern.outcome is OrderOutcome.CONCERN
            assert "not the antisense strand" in concern.detail

    def test_a_sense_guide_prediction_is_flagged_against_the_design_intent(self) -> None:
        # The duplex is designed so the antisense strand is the guide. All three tables agree
        # this target loads the sense strand, so ordering it risks the strand nobody intended.
        order = order_for(candidate_for(bases=SENSE_BASES, target=SENSE_PREDICTED))
        intent = check(order, "guide-matches-intent")
        assert intent.outcome is OrderOutcome.CONCERN
        assert "may act through the strand you did not intend" in intent.detail
        # The disagreement check is a separate question and must not fire here.
        assert check(order, "model-agreement").outcome is OrderOutcome.OK

    def test_an_antisense_prediction_matching_the_design_passes(self) -> None:
        order = order_for(candidate_for())
        assert check(order, "guide-matches-intent").outcome is OrderOutcome.OK
        assert order.concerns == ()


class TestUnsupportedLengths:
    def test_a_core_that_is_not_nineteen_is_a_concern_and_names_the_supported_length(self) -> None:
        # Order chemistry is fixed at 19 nt + dTdT (ADR 0013, accepted constraint 4). An 18 nt
        # core is designable but not orderable under that contract, and must say so.
        result = design_sarnas(
            Promoter(sequence=NucleicSequence.parse(BASES)),
            WANG_2015.resolve({"target-length-19": {"target-length": 18}}),
        )
        order = order_for(result.generated[0])
        length = check(order, "core-length")
        assert length.outcome is OrderOutcome.CONCERN
        assert "19" in length.detail
        assert order.sense.length_nt == 20


class TestExport:
    def test_an_order_exports_every_field_a_vendor_form_needs(self) -> None:
        exported = order_for(candidate_for()).to_dict()
        assert exported["sense"]["annotated"] == "GGGUAAUUUUGACAGGUCA[dT][dT]"
        assert exported["antisense"]["annotated"] == "UGACCUGUCAAAAUUACCC[dT][dT]"
        assert exported["sense"]["length_nt"] == 21
        assert exported["target_site_dna"] == TARGET
        assert [c["id"] for c in exported["checks"]]

    def test_the_export_separates_core_from_overhang_for_rewriting(self) -> None:
        # Vendor notation for a chimeric oligo differs between suppliers. Emitting the core and
        # the overhang separately lets a researcher reformat without re-deriving the chemistry.
        exported = order_for(candidate_for()).to_dict()
        assert exported["sense"]["core_rna_5to3"] == SENSE_CORE
        assert exported["sense"]["overhang"] == OVERHANG_ANNOTATED

    def test_a_saved_record_produces_the_identical_duplex_and_report(self) -> None:
        # A reopened run must order the same molecule as the live run did. Two code paths that
        # could disagree about one duplex is the failure this guards: the record path is what
        # the app's export actually uses.
        for target, bases in ((TARGET, BASES), (DIVIDED, DIVIDED_BASES),
                              (SENSE_PREDICTED, SENSE_BASES)):
            candidate = candidate_for(bases=bases, target=target)
            assert order_from_record(candidate.to_dict()).to_dict() \
                == order_for(candidate).to_dict()

    def test_the_csv_row_carries_the_sanity_report_beside_the_sequences(self) -> None:
        # The caveats must not be separable from the sequences by pasting one row elsewhere.
        order = order_for(candidate_for(bases=SENSE_BASES, target=SENSE_PREDICTED))
        row = dict(zip(CSV_COLUMNS, csv_row(order, "run-1"), strict=True))
        assert row["sense_order_5to3"].endswith("[dT][dT]")
        assert row["antisense_order_5to3"].endswith("[dT][dT]")
        assert row["ordered_length_nt"] == "21"
        assert "off-target-screening [not-performed]" in row["sanity_report"]
        assert "guide-matches-intent [concern]" in row["sanity_report"]
        assert row["concerns"] == "1"

    def test_an_unknown_model_is_rejected_rather_than_silently_unreported(self) -> None:
        with pytest.raises(KeyError, match="no-such-table"):
            order_for(candidate_for(), model="no-such-table")
