"""Specificity measures: several explicit numbers, never one composite score.

[ADR 0019](../../docs/adr/0019-off-target-screening.md) decides that a candidate's specificity is
reported as separate, separately sortable measures rather than a single score. The reason is in
these tests: a candidate with one perfect 19/19 match elsewhere is a different experimental
problem from one with forty weak partial matches, and any weighted sum can rank them the same.
ADR 0006 already forbids weighted sums for ranking; this applies the same rule to screening.

Nothing here filters. A finding reports.
"""

import pytest

from sarna.offtarget.models import OffTargetHit
from sarna.offtarget.specificity import (
    NEAR_PERFECT_IDENTITIES,
    WANG_HOMOLOGY_FRACTION,
    measure_specificity,
    summarise,
    wang_homology_identities,
)

CORE = 19


def hit(
    candidate: str = "a",
    identities: int = 19,
    length: int = 19,
    evalue: float = 1e-5,
    accession: str = "NM_000001.1",
    gaps: int = 0,
) -> OffTargetHit:
    """Build one hit, defaulting to a perfect full-length match."""
    return OffTargetHit(
        candidate_id=candidate,
        accession=accession,
        title=f"Homo sapiens something ({accession})",
        identities=identities,
        alignment_length=length,
        gaps=gaps,
        query_start=1,
        query_end=length,
        subject_start=100,
        subject_end=100 + length,
        subject_strand="plus",
        evalue=evalue,
        bit_score=38.2,
    )


class TestTheMeasuresAreSeparate:
    def test_a_candidate_with_no_hits_reports_zero_rather_than_a_missing_value(self) -> None:
        measured = measure_specificity("saRNA-001", (), core_length=CORE)
        assert measured.hit_count == 0
        assert measured.perfect_hits == 0
        assert measured.near_perfect_hits == 0
        assert measured.longest_complementary_stretch == 0
        # No alignment means no E-value exists. Reporting 0.0 would read as a perfect match.
        assert measured.best_evalue is None

    def test_one_perfect_match_and_forty_weak_ones_are_not_the_same_finding(self) -> None:
        # The distinction a single composite score would destroy.
        one_perfect = measure_specificity(
            "a", (hit(identities=19, evalue=1e-8),), core_length=CORE
        )
        many_weak = measure_specificity(
            "b",
            tuple(hit(candidate="b", identities=11, length=11, evalue=5.0,
                      accession=f"NM_{i:06d}.1")
                  for i in range(40)),
            core_length=CORE,
        )
        assert one_perfect.perfect_hits == 1 and many_weak.perfect_hits == 0
        assert one_perfect.hit_count == 1 and many_weak.hit_count == 40
        assert one_perfect.longest_complementary_stretch == 19
        assert many_weak.longest_complementary_stretch == 11

    def test_a_perfect_hit_needs_full_length_and_no_gaps(self) -> None:
        assert measure_specificity("a", (hit(identities=19, length=19),), CORE).perfect_hits == 1
        # Same identity count but shorter alignment: not the whole duplex.
        assert measure_specificity("a", (hit(identities=19, length=25),), CORE).perfect_hits == 0
        # A gapped alignment is not a perfect complement even at full length.
        assert measure_specificity("a", (hit(identities=19, gaps=1),), CORE).perfect_hits == 0

    def test_near_perfect_uses_the_declared_identity_threshold(self) -> None:
        assert NEAR_PERFECT_IDENTITIES == 17
        at = measure_specificity("a", (hit(identities=17, length=17),), CORE)
        below = measure_specificity("a", (hit(identities=16, length=16),), CORE)
        assert at.near_perfect_hits == 1
        assert below.near_perfect_hits == 0

    def test_a_perfect_hit_also_counts_as_near_perfect(self) -> None:
        # Nested, not exclusive: "how many close matches" should not drop the closest one.
        measured = measure_specificity("a", (hit(identities=19),), CORE)
        assert measured.perfect_hits == 1
        assert measured.near_perfect_hits == 1

    def test_the_best_evalue_is_the_smallest_one(self) -> None:
        measured = measure_specificity(
            "a", (hit(evalue=3.0), hit(evalue=1e-9), hit(evalue=0.2)), CORE
        )
        assert measured.best_evalue == 1e-9


class TestSortingWithoutScoring:
    def test_every_measure_is_separately_available_for_sorting(self) -> None:
        exported = measure_specificity("a", (hit(),), CORE).to_dict()
        for column in (
            "hit_count", "best_evalue", "near_perfect_hits",
            "perfect_hits", "longest_complementary_stretch",
        ):
            assert column in exported, column

    def test_there_is_no_single_specificity_score(self) -> None:
        # ADR 0006 and ADR 0019: a composite would hide which risk a candidate carries.
        exported = measure_specificity("a", (hit(),), CORE).to_dict()
        assert not any(
            key in exported for key in ("score", "specificity_score", "rank", "weighted")
        ), exported.keys()

    def test_candidates_with_no_hits_sort_cleanly_against_those_with_some(self) -> None:
        # None for best_evalue must not crash a sort or silently sort as the best value.
        clean = measure_specificity("clean", (), CORE)
        dirty = measure_specificity("dirty", (hit(candidate="dirty", evalue=1e-9),), CORE)
        ordered = sorted([dirty, clean], key=lambda m: m.sort_key_best_evalue)
        assert [m.candidate_id for m in ordered] == ["dirty", "clean"]


class TestTheSummaryIsHonestAboutWhatItMeans:
    def test_a_clean_candidate_does_not_claim_to_be_specific(self) -> None:
        # BLAST finding nothing is not proof of specificity: it is one search, of one database,
        # at one moment, with parameters that can miss things.
        text = summarise(measure_specificity("a", (), CORE))
        assert "no" in text.lower()
        assert "does not establish" in text.lower() or "not proof" in text.lower()

    def test_a_perfect_off_target_is_stated_plainly(self) -> None:
        text = summarise(measure_specificity("a", (hit(identities=19),), CORE))
        assert "perfect" in text.lower()
        assert "1" in text

    def test_the_summary_never_tells_the_researcher_to_discard_a_candidate(self) -> None:
        # Findings report rather than exclude (accepted constraint 15).
        text = summarise(measure_specificity("a", (hit(identities=19),), CORE)).lower()
        for word in ("discard", "reject", "excluded", "do not use", "unusable"):
            assert word not in text, f"the summary should not instruct: {word!r}"


class TestInputValidation:
    def test_a_hit_belonging_to_another_candidate_is_refused(self) -> None:
        # Silently measuring another candidate's hits would attribute a risk to the wrong duplex.
        with pytest.raises(ValueError, match="saRNA-002"):
            measure_specificity("saRNA-001", (hit(candidate="saRNA-002"),), CORE)

    def test_an_impossible_core_length_is_refused(self) -> None:
        with pytest.raises(ValueError, match="core length"):
            measure_specificity("a", (), core_length=0)


class TestWangsOwnHomologyThreshold:
    """80% is the one homology figure Wang 2015 states, and ADR 0023 reports it beside 19/19.

    The paper uses it to reject a **control** duplex showing "significant homology", not to
    assess a candidate's off-targets, so carrying it across is an extension by analogy. The
    cost is in the ADR; what these tests pin is that the number means what the paper says.
    """

    def test_the_threshold_is_the_smallest_count_that_reaches_eighty_percent(self) -> None:
        # Derived from the definition rather than from the implementation: whatever the
        # function returns must clear 80%, and one base fewer must fall short of it. A plain
        # assertion that 19 gives 16 would pass just as well against round() or int().
        for core in range(5, 31):
            threshold = wang_homology_identities(core)
            assert threshold / core >= WANG_HOMOLOGY_FRACTION, core
            assert (threshold - 1) / core < WANG_HOMOLOGY_FRACTION, core

    def test_a_nineteen_base_core_asks_for_sixteen(self) -> None:
        # 15 of 19 is 78.9% and does not reach the published figure; 16 of 19 is 84.2%.
        assert wang_homology_identities(19) == 16

    def test_sixteen_of_nineteen_counts_and_fifteen_does_not(self) -> None:
        measured = measure_specificity(
            "a", (hit(identities=16), hit(identities=15)), CORE
        )
        assert measured.wang_homology_hits == 1
        assert measured.strong_threshold == 16

    def test_the_locally_chosen_seventeen_is_still_reported_beside_it(self) -> None:
        # Saved runs and the xlsx export carry near_perfect_hits. Changing what a stored
        # column means is worse than carrying two, so both travel.
        measured = measure_specificity("a", (hit(identities=16),), CORE)
        assert measured.near_perfect_hits == 0, "16 does not reach the local 17"
        assert measured.wang_homology_hits == 1, "but it does reach Wang's 80%"
        assert measured.to_dict()["near_perfect_hits"] == 0

    def test_the_threshold_travels_in_the_export_so_a_column_can_label_itself(self) -> None:
        exported = measure_specificity("a", (), CORE).to_dict()
        assert exported["strong_threshold"] == 16


class TestPromiscuityIsADisplayStateAndSaysSo:
    def test_five_matches_at_wangs_threshold_is_called_promiscuous(self) -> None:
        measured = measure_specificity("a", tuple(hit(identities=16) for _ in range(5)), CORE)
        assert measured.is_promiscuous
        assert "reported alignments" in summarise(measured)

    def test_four_is_not(self) -> None:
        measured = measure_specificity("a", tuple(hit(identities=16) for _ in range(4)), CORE)
        assert not measured.is_promiscuous

    def test_weak_alignments_do_not_count_towards_it(self) -> None:
        # Forty 13-base partial matches are noise, and the whole point of counting at a stated
        # threshold is that noise cannot make a candidate look promiscuous.
        measured = measure_specificity("a", tuple(hit(identities=13) for _ in range(40)), CORE)
        assert not measured.is_promiscuous
        assert measured.hit_count == 40

    def test_it_changes_no_measure_and_no_rank(self) -> None:
        # A display state must leave every number alone, or it is a filter wearing a disguise.
        hits = tuple(hit(identities=16) for _ in range(6))
        measured = measure_specificity("a", hits, CORE)
        assert measured.hit_count == 6
        assert measured.wang_homology_hits == 6
        assert "rank" not in measured.to_dict()


class TestACappedListIsNotACleanOne:
    def test_a_capped_candidate_says_its_counts_are_floors(self) -> None:
        measured = measure_specificity("a", (hit(identities=16),), CORE, capped=True)
        assert measured.alignments_capped
        assert "floor" in summarise(measured)

    def test_an_uncapped_candidate_makes_no_such_claim(self) -> None:
        assert "floor" not in summarise(measure_specificity("a", (hit(),), CORE))


class TestRemovedOwnLocusHitsAreCounted:
    def test_a_zero_that_had_something_taken_out_of_it_says_so(self) -> None:
        # Otherwise it reads exactly like a zero that never did, and ADR 0023 removes the hit
        # from the findings precisely so that the remaining number can be trusted.
        measured = measure_specificity("a", (), CORE, on_target_hits=1)
        assert measured.on_target_hits == 1
        assert "own locus" in summarise(measured)

    def test_nothing_is_said_when_nothing_was_removed(self) -> None:
        assert "own locus" not in summarise(measure_specificity("a", (), CORE))
