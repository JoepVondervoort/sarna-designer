"""Tests for the end-to-end design run: generate, gate, rank.

These are the first tests of the engine as a whole. They care about four things:

- that the four candidate populations stay distinguishable (generated, rejected, passing, ranked);
- that a rejected candidate says which rule rejected it, because a sequence that simply vanishes
  is the behaviour this project exists to replace;
- that ranking follows ADR 0006 and is deterministic;
- that a run works, with reduced information, on a pasted promoter with no coordinates.
"""

from itertools import pairwise

import pytest

from sarna.design.config import DesignConfig
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import LEGACY_SCANNER, WANG_2015
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.rules.catalog import CATALOGUE
from sarna.rules.results import RuleOutcome
from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.guide import GuideStrand

# A 300 bp synthetic promoter with varied composition, so a scan yields a mixed candidate set
# rather than all-pass or all-fail.
PROMOTER_BASES = (
    "GCTACGGTATCATGTCATAGGCTTAACGGATTCAAGTTACCGATTACAGGCTAACGTTAGCCATTAAGCGTTACA"
    "ATCGGTTACAGCTTAAGGCTATTACGGCATTAACGGTTACAGCTTAAGGCATTAACGGTTACAGCTTAAGGCTAT"
    "TACGGCATTAACGGTTACAGCTTAAGGCATTAACGTTACAGGCTTAAGGCATTAACGGTTACAGCTTAAGGCATT"
    "AACGGTTACAGCTTAAGGCATTAACGGTTACAGCTTAAGGCATTAACGGTTACAGCTTAAGGCATTAACGGTTAC"
)


def sequence_only() -> Promoter:
    """A pasted promoter: no assembly, no coordinates, no TSS (Workflow B)."""
    return Promoter(sequence=NucleicSequence.parse(PROMOTER_BASES))


def anchored(strand: Strand = Strand.PLUS) -> Promoter:
    """A promoter anchored to a genome and a TSS, 1000 bp upstream of it."""
    length = len(PROMOTER_BASES)
    return Promoter(
        sequence=NucleicSequence.parse(PROMOTER_BASES),
        tss_start=TssPosition(-1000),
        interval=GenomicInterval(Assembly.GRCH38, "1", 100_000, 100_000 + length, strand),
    )


class TestGeneration:
    def test_every_window_of_the_target_length_becomes_a_candidate(self) -> None:
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert len(result.generated) == len(PROMOTER_BASES) - 19 + 1

    def test_the_scan_starts_at_the_first_base(self) -> None:
        # Defect B6: the original script starts at index 3 and silently discards the first three
        # windows with no stated reason.
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert result.generated[0].offset == 0
        assert result.generated[0].target.bases == PROMOTER_BASES[:19]

    def test_candidate_ids_are_stable_and_follow_position(self) -> None:
        # Numbered by position, not by rank, so an id refers to the same sequence however the
        # ranking changes.
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert result.generated[0].id == "saRNA-001"
        assert result.generated[1].id == "saRNA-002"
        assert [c.offset for c in result.generated] == sorted(
            c.offset for c in result.generated
        )

    def test_a_promoter_shorter_than_the_target_yields_nothing(self) -> None:
        # A legitimate outcome reported as such, not an error.
        short = Promoter(sequence=NucleicSequence.parse("GCTACGGTA"))
        result = design_sarnas(short, WANG_2015.resolve())
        assert result.generated == ()
        assert result.ranked == ()

    def test_the_result_records_the_promoter_and_configuration(self) -> None:
        promoter = sequence_only()
        config = WANG_2015.resolve()
        result = design_sarnas(promoter, config)
        assert result.promoter is promoter
        assert result.config is config


class TestPartitioning:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_the_populations_partition_the_generated_set(self) -> None:
        assert len(self.RESULT.rejected) + len(self.RESULT.passing) == len(
            self.RESULT.generated
        )

    def test_ranked_contains_exactly_the_passing_candidates(self) -> None:
        assert {c.id for c in self.RESULT.ranked} == {c.id for c in self.RESULT.passing}

    def test_some_candidates_are_rejected_and_some_pass(self) -> None:
        # Guards the guard: a fixture where everything passes or everything fails would make the
        # partitioning tests vacuous.
        assert self.RESULT.rejected
        assert self.RESULT.passing

    def test_a_rejected_candidate_is_not_silently_dropped(self) -> None:
        assert all(candidate in self.RESULT.generated for candidate in self.RESULT.rejected)


class TestRejectionIsExplained:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_every_rejected_candidate_names_the_rules_that_excluded_it(self) -> None:
        for candidate in self.RESULT.rejected:
            assert candidate.excluded_by, f"{candidate.id} was rejected with no reason"

    def test_the_named_rule_actually_failed(self) -> None:
        for candidate in self.RESULT.rejected:
            for result in candidate.excluded_by:
                assert result.outcome is RuleOutcome.FAIL
                assert result.excludes

    def test_the_rejection_reason_states_what_was_measured(self) -> None:
        candidate = self.RESULT.rejected[0]
        assert candidate.rejection_summary
        assert candidate.excluded_by[0].detail in candidate.rejection_summary

    def test_a_passing_candidate_has_nothing_excluding_it(self) -> None:
        for candidate in self.RESULT.passing:
            assert candidate.excluded_by == ()

    def test_every_candidate_carries_a_result_for_every_selected_rule(self) -> None:
        config = WANG_2015.resolve()
        for candidate in self.RESULT.generated:
            assert {r.rule_id for r in candidate.results} == {
                resolved.rule.id for resolved in config
            }


class TestRanking:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_ranks_are_consecutive_from_one(self) -> None:
        assert [c.rank for c in self.RESULT.ranked] == list(
            range(1, len(self.RESULT.ranked) + 1)
        )

    def test_preference_points_never_increase_down_the_ranking(self) -> None:
        points = [c.preference_points for c in self.RESULT.ranked]
        assert points == sorted(points, reverse=True)

    def test_a_tie_on_points_is_broken_by_asymmetry_descending(self) -> None:
        ranked = self.RESULT.ranked
        for earlier, later in pairwise(ranked):
            if earlier.preference_points != later.preference_points:
                continue
            if earlier.tiebreak_value is None or later.tiebreak_value is None:
                continue
            assert earlier.tiebreak_value >= later.tiebreak_value

    def test_the_denominator_is_the_number_of_selected_preferences(self) -> None:
        # ADR 0006 presents a rank as "meets 5 of 7 preferences".
        for candidate in self.RESULT.ranked:
            assert candidate.max_preference_points == 7
            assert 0 <= candidate.preference_points <= 7

    def test_an_unevaluated_preference_is_left_out_of_the_denominator(self) -> None:
        # A preference nobody could check is not one the candidate failed to meet, so it must
        # not count against it.
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        last = result.generated[-1]
        unevaluated = [
            r for r in last.results if r.outcome is RuleOutcome.NOT_APPLICABLE and r.points == 0
        ]
        assert unevaluated, "the last candidate should have an unevaluable flank preference"
        assert last.evaluated_preference_points < last.max_preference_points

    def test_the_ranking_is_deterministic(self) -> None:
        again = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert [c.id for c in again.ranked] == [c.id for c in self.RESULT.ranked]

    def test_ranking_explains_itself_per_candidate(self) -> None:
        top = self.RESULT.ranked[0]
        assert "of 7" in top.score_summary
        assert top.met_preferences or top.unmet_preferences


class TestStrandHandling:
    def test_a_minus_strand_run_reports_coordinates_descending_with_offset(self) -> None:
        result = design_sarnas(anchored(Strand.MINUS), WANG_2015.resolve())
        first, second = result.generated[0], result.generated[1]
        assert first.interval is not None
        assert second.interval is not None
        assert second.interval.start < first.interval.start

    def test_tss_positions_are_identical_on_both_strands(self) -> None:
        # Measured in the transcript's own frame, so strand cannot affect them.
        plus = design_sarnas(anchored(Strand.PLUS), WANG_2015.resolve())
        minus = design_sarnas(anchored(Strand.MINUS), WANG_2015.resolve())
        assert [c.tss_span for c in plus.generated] == [
            c.tss_span for c in minus.generated
        ]

    def test_the_same_candidates_pass_on_both_strands(self) -> None:
        # The promoter sequence is identical in the transcript's frame, so the verdicts must be.
        plus = design_sarnas(anchored(Strand.PLUS), WANG_2015.resolve())
        minus = design_sarnas(anchored(Strand.MINUS), WANG_2015.resolve())
        assert [c.id for c in plus.ranked] == [c.id for c in minus.ranked]


class TestDuplexSequences:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_the_sense_strand_is_the_target_transcribed(self) -> None:
        candidate = self.RESULT.generated[0]
        assert candidate.sense_strand.bases == candidate.target.as_rna().bases

    def test_the_antisense_strand_is_the_reverse_complement(self) -> None:
        candidate = self.RESULT.generated[0]
        expected = candidate.target.as_rna().reverse_complement()
        assert candidate.antisense_strand == expected

    def test_both_strands_are_rna(self) -> None:
        candidate = self.RESULT.generated[0]
        assert "T" not in candidate.sense_strand.bases
        assert "T" not in candidate.antisense_strand.bases

    def test_the_target_site_stays_dna(self) -> None:
        # It is a genomic locus, not a delivered molecule.
        assert "U" not in self.RESULT.generated[0].target.bases


class TestGuideComparison:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_every_candidate_carries_an_assessment_per_comparison_table(self) -> None:
        # ADR 0009: all selected tables reported side by side, never averaged.
        candidate = self.RESULT.generated[0]
        assert len(candidate.guide_comparison.assessments) == 3

    def test_each_assessment_explains_itself(self) -> None:
        for assessment in self.RESULT.generated[0].guide_comparison.assessments:
            assert "Schwarz" in assessment.explanation
            assert assessment.table_name

    def test_the_gate_uses_only_the_configured_table(self) -> None:
        # The comparison is context; it must not change which candidates pass.
        candidate = self.RESULT.generated[0]
        gate = next(
            r for r in candidate.results if r.rule_id == "thermodynamic-asymmetry"
        )
        assert gate.measurements["nearest-neighbour-table"] == "rna-xia-1998"

    def test_disagreement_is_visible_where_it_occurs(self) -> None:
        disagreeing = [
            candidate
            for candidate in self.RESULT.generated
            if not candidate.guide_comparison.tables_agree
        ]
        assert disagreeing, "the fixture should contain at least one disagreeing candidate"
        assert "disagree" in disagreeing[0].guide_comparison.summary


class TestWorkflowBDegradesGracefully:
    def test_a_pasted_promoter_still_yields_ranked_candidates(self) -> None:
        # The failure this guards: if the TSS rules failed rather than reporting
        # not-applicable, every candidate of every pasted promoter would be rejected.
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert result.ranked

    def test_the_tss_rules_report_not_applicable_rather_than_failing(self) -> None:
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        candidate = result.generated[0]
        window = next(r for r in candidate.results if r.rule_id == "tss-window")
        assert window.outcome is RuleOutcome.NOT_APPLICABLE
        assert not window.excludes

    def test_the_run_lists_which_rules_could_not_be_evaluated(self) -> None:
        # Surfaced at the run level so a user is told once, rather than discovering it per
        # candidate or not at all.
        result = design_sarnas(sequence_only(), WANG_2015.resolve())
        assert "tss-window" in result.unevaluated_rule_ids
        assert "tss-preferred-window" in result.unevaluated_rule_ids

    def test_an_anchored_run_evaluates_the_tss_rules(self) -> None:
        result = design_sarnas(anchored(), WANG_2015.resolve())
        assert "tss-window" not in result.unevaluated_rule_ids


class TestIslandDetectionRunsOncePerPromoter:
    def test_islands_are_detected_when_the_rule_is_selected(self) -> None:
        promoter = Promoter(sequence=NucleicSequence.parse("CG" * 400 + PROMOTER_BASES))
        result = design_sarnas(promoter, WANG_2015.resolve())
        assert result.cpg_islands
        island_rule = next(
            r for r in result.generated[0].results if r.rule_id == "not-in-cpg-island"
        )
        assert island_rule.outcome is RuleOutcome.FAIL

    def test_detection_is_skipped_when_the_rule_is_not_selected(self) -> None:
        # legacy-scanner excludes it (defect B8), so running detection would be wasted work
        # and would imply a definition the configuration never chose.
        result = design_sarnas(sequence_only(), LEGACY_SCANNER.resolve())
        assert result.cpg_islands is None


class TestPresetsDiffer:
    def test_the_legacy_preset_admits_a_different_candidate_set(self) -> None:
        promoter = sequence_only()
        wang = design_sarnas(promoter, WANG_2015.resolve())
        legacy = design_sarnas(promoter, LEGACY_SCANNER.resolve())
        assert {c.id for c in wang.passing} != {c.id for c in legacy.passing}

    def test_both_presets_score_out_of_seven_but_not_the_same_seven(self) -> None:
        # Legacy drops the preferred-window preference, having no TSS concept, and adds the
        # unsourced pos7-u, which the original script weighted as a dominant term. The totals
        # coincide at 7; the composition does not, so the two scores are not comparable.
        wang_config = WANG_2015.resolve()
        legacy_config = LEGACY_SCANNER.resolve()
        wang_preferences = {r.rule.id for r in wang_config.preferences}
        legacy_preferences = {r.rule.id for r in legacy_config.preferences}

        assert len(wang_preferences) == len(legacy_preferences) == 7
        assert wang_preferences != legacy_preferences
        assert wang_preferences - legacy_preferences == {"tss-preferred-window"}
        assert legacy_preferences - wang_preferences == {"pos7-u"}

    def test_the_legacy_preset_scores_the_unsourced_position_seven_rule(self) -> None:
        legacy = design_sarnas(sequence_only(), LEGACY_SCANNER.resolve())
        candidate = legacy.generated[0]
        assert candidate.max_preference_points == 7
        assert any(result.rule_id == "pos7-u" for result in candidate.results)


class TestExport:
    def test_the_result_exports_with_its_configuration(self) -> None:
        import json

        result = design_sarnas(anchored(), WANG_2015.resolve())
        exported = result.to_dict()
        json.dumps(exported)
        assert exported["configuration"]["preset"]["id"] == "wang-2015"
        assert exported["counts"]["generated"] == len(result.generated)

    def test_each_exported_candidate_carries_its_rule_results(self) -> None:
        result = design_sarnas(anchored(), WANG_2015.resolve())
        first = result.to_dict()["candidates"][0]
        assert first["id"] == "saRNA-001"
        assert len(first["rules"]) == len(WANG_2015.resolve())
        assert first["rules"][0]["detail"]

    def test_rejected_candidates_appear_in_the_export(self) -> None:
        # Never simply make filtered sequences disappear.
        result = design_sarnas(anchored(), WANG_2015.resolve())
        exported = result.to_dict()
        assert any(entry["excluded_by"] for entry in exported["candidates"])


class TestValidation:
    def test_a_configuration_without_a_target_length_is_rejected(self) -> None:
        from sarna.design.config import DesignConfig

        with pytest.raises(KeyError, match="target-length"):
            design_sarnas(sequence_only(), DesignConfig.of(["gc-content"]))


class TestCandidateAccessors:
    RESULT = design_sarnas(sequence_only(), WANG_2015.resolve())

    def test_a_result_can_be_fetched_by_rule_id(self) -> None:
        candidate = self.RESULT.generated[0]
        assert candidate.result("gc-content").rule_id == "gc-content"

    def test_asking_for_an_unselected_rule_explains_why_it_is_missing(self) -> None:
        # "not selected by the configuration" rather than "no such rule": a user reading this
        # needs to know their preset omitted it, not that they mistyped.
        candidate = self.RESULT.generated[0]
        with pytest.raises(KeyError, match="not selected by"):
            candidate.result("pos7-u")

    def test_a_candidate_without_a_tiebreak_rule_has_no_tiebreak_value(self) -> None:
        from sarna.design.config import DesignConfig

        config = DesignConfig.of(["target-length-19", "gc-content"])
        result = design_sarnas(sequence_only(), config)
        assert result.generated[0].tiebreak_value is None

    def test_candidates_without_tiebreaks_still_rank_deterministically(self) -> None:
        # An uncomputed asymmetry is not a neutral asymmetry, so these sort below any candidate
        # that has one; with none of them having one, position decides.
        from sarna.design.config import DesignConfig

        config = DesignConfig.of(["target-length-19", "gc-content", "start-1-gc"])
        result = design_sarnas(sequence_only(), config)
        assert [c.rank for c in result.ranked] == list(range(1, len(result.ranked) + 1))


class TestDeclaredParameterRangesAreReachable:
    """A declared range is a promise to the user, so every value in it must run.

    ``target-length`` advertised 18-21 while 20 and 21 raised ``ValueError`` out of
    ``flank_base_at``: the flank preferences are registered at literal positions 20-23 and
    collide with the target as it lengthens. Nothing caught it because no test ran the engine
    at a non-default length. A range nobody exercises is a range nobody has checked.
    """

    def test_every_declared_target_length_produces_candidates(self) -> None:
        parameter = CATALOGUE["target-length-19"].parameter("target-length")
        assert parameter.minimum is not None and parameter.maximum is not None
        for length in range(int(parameter.minimum), int(parameter.maximum) + 1):
            config = WANG_2015.resolve({"target-length-19": {"target-length": length}})
            result = design_sarnas(sequence_only(), config)
            assert result.generated, f"target-length {length} produced no candidates"

    def test_the_default_lies_inside_the_declared_range(self) -> None:
        parameter = CATALOGUE["target-length-19"].parameter("target-length")
        assert parameter.minimum is not None
        assert parameter.maximum is not None
        assert float(parameter.minimum) <= float(parameter.default) <= float(parameter.maximum)


class TestTheStrandPredictionHonoursTheConfiguredEndLength:
    """The reported guide strand must be computed at the end length the run actually evaluates.

    ``_compare`` once counted four terminal bases regardless of configuration. A run that
    deliberately evaluated a different end length would then be shown the *opposite* strand from
    the one its own asymmetry rule admitted the candidate on, a wrong duplex to order, reported
    with full confidence and agreed on by every table.

    ``ASYMMETRY_FLIP`` is the fixture that makes the difference observable: a 19-mer found by
    enumerating random targets until one had all three tables agreeing on antisense at the
    default four-base end and all three agreeing on sense at six. Both verdicts are unanimous, so
    neither assertion rests on a near-tolerance tie.
    """

    # All three tables agree: antisense at end-length 4, sense at end-length 6.
    ASYMMETRY_FLIP = "GGGTAATTTTGACAGGTCA"
    BASES = (
        "GCTTAACGGATTCAAGTTACCGATTACAGGCTAACGTTAGCCATTAAGCGTTACA"
        + ASYMMETRY_FLIP
        + "ATCGGTTACAGCTTAAGGCTATTACGGCATTAACGGTTACAGCTTAAGGCATTAA"
    )
    OFFSET = BASES.index(ASYMMETRY_FLIP)

    def predictions(self, config: DesignConfig) -> frozenset[GuideStrand]:
        """Return the strands predicted for the flipping candidate under ``config``."""
        result = design_sarnas(
            Promoter(sequence=NucleicSequence.parse(self.BASES)), config
        )
        candidate = next(c for c in result.generated if c.offset == self.OFFSET)
        assert candidate.target.bases == self.ASYMMETRY_FLIP
        return candidate.guide_comparison.predictions

    def test_the_default_four_base_end_predicts_the_antisense_strand(self) -> None:
        assert self.predictions(WANG_2015.resolve()) == {GuideStrand.ANTISENSE}

    def test_a_configured_six_base_end_predicts_the_other_strand(self) -> None:
        # Fails if the comparison hard-codes four: it would still answer antisense here.
        config = WANG_2015.resolve({"thermodynamic-asymmetry": {"end-length-nt": 6}})
        assert self.predictions(config) == {GuideStrand.SENSE}

    def test_the_gate_decides_the_end_length_rather_than_the_tiebreak(self) -> None:
        # Both rules declare end-length-nt. The gate admitted the candidate, so the strand
        # reported beside it is the one that gate's end length implies.
        config = WANG_2015.resolve({
            "thermodynamic-asymmetry": {"end-length-nt": 4},
            "asymmetry-margin": {"end-length-nt": 6},
        })
        assert self.predictions(config) == {GuideStrand.ANTISENSE}

    def test_the_tiebreak_supplies_the_end_length_when_the_gate_is_disabled(self) -> None:
        config = WANG_2015.resolve(
            {"asymmetry-margin": {"end-length-nt": 6}},
            role_overrides={"thermodynamic-asymmetry": None},
        )
        assert "thermodynamic-asymmetry" not in config
        assert self.predictions(config) == {GuideStrand.SENSE}

    def test_a_run_with_neither_asymmetry_rule_still_predicts_a_strand(self) -> None:
        # Nothing configured an end length, so the documented default of four applies rather
        # than the comparison failing or reporting nothing.
        config = WANG_2015.resolve(role_overrides={
            "thermodynamic-asymmetry": None, "asymmetry-margin": None,
        })
        assert "asymmetry-margin" not in config
        assert self.predictions(config) == {GuideStrand.ANTISENSE}


class TestEmptyRun:
    def test_a_run_with_no_candidates_reports_no_unevaluated_rules(self) -> None:
        # Nothing was evaluated, so claiming a rule "could not be evaluated" would be a
        # statement about zero candidates.
        short = Promoter(sequence=NucleicSequence.parse("GCTACGGTA"))
        result = design_sarnas(short, WANG_2015.resolve())
        assert result.unevaluated_rule_ids == ()
        assert result.to_dict()["counts"]["generated"] == 0
