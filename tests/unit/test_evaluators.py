"""Tests for the rule evaluators: sequence, position and location rules.

Each evaluator is a pure function of a context and a resolved rule (ADR 0007), so every case
here is built from a hand-made context with no promoter file, no network and no fixtures.

Three things are checked for every rule: the verdict, that the stated reason names the measured
value, and — where the rule depends on something a sequence-only run does not have — that it
reports not-applicable rather than passing or failing.
"""

import pytest

from sarna.design.context import CandidateContext
from sarna.design.evaluators import EVALUATORS, evaluate, evaluator_for
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.rules.catalog import CATALOGUE
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleOutcome
from sarna.rules.rule import RuleKind
from sarna.sequence.cpg import (
    GARDINER_GARDEN_FROMMER_1987,
    CpGIsland,
    find_cpg_islands,
)
from sarna.sequence.nucleic import NucleicSequence


def context_for(
    target: str,
    flank: str = "TTTT",
    tss_start: int | None = None,
    cpg_islands: tuple[CpGIsland, ...] | None = None,
) -> CandidateContext:
    """Build a context whose target is exactly ``target``, followed by ``flank``."""
    promoter = Promoter(
        sequence=NucleicSequence.parse(target + flank),
        tss_start=TssPosition(tss_start) if tss_start is not None else None,
    )
    return CandidateContext.at(
        promoter, offset=0, length=len(target), cpg_islands=cpg_islands
    )


def check(rule_id: str, context: CandidateContext, **overrides: object):
    """Evaluate one catalogued rule against ``context``."""
    rule = CATALOGUE[rule_id]
    return evaluate(context, ResolvedRule(rule, overrides))  # type: ignore[arg-type]


# A 19-mer satisfying every mandatory sequence rule: GC 42.1 %, one CpG, longest run 2, no
# simple repeats, position 18 T, position 19 A, positions 1-2 G/C, position 7 G.
GOOD = "GCTACGGTATCATGTCATA"


class TestEveryCatalogueRuleHasAnEvaluator:
    def test_no_rule_is_silently_unevaluated(self) -> None:
        # A rule with no evaluator would never fire, and a rule that never fires is worse than
        # a missing rule: it looks like it is protecting you.
        missing = [rule.id for rule in CATALOGUE if rule.id not in EVALUATORS]
        assert missing == [], f"rules with no evaluator: {missing}"

    def test_no_evaluator_is_registered_for_an_unknown_rule(self) -> None:
        unknown = [rule_id for rule_id in EVALUATORS if rule_id not in CATALOGUE]
        assert unknown == []


class TestTargetLength:
    def test_a_nineteen_mer_passes(self) -> None:
        assert check("target-length-19", context_for(GOOD)).outcome is RuleOutcome.PASS

    def test_a_shorter_target_fails(self) -> None:
        result = check("target-length-19", context_for(GOOD[:18]))
        assert result.outcome is RuleOutcome.FAIL
        assert "18" in result.detail

    def test_the_required_length_is_configurable(self) -> None:
        result = check("target-length-19", context_for(GOOD[:18]), **{"target-length": 18})
        assert result.outcome is RuleOutcome.PASS

    def test_the_measurement_carries_the_length(self) -> None:
        result = check("target-length-19", context_for(GOOD))
        assert result.measurements == {"length": 19}


class TestGcContent:
    def test_a_target_inside_the_band_passes(self) -> None:
        result = check("gc-content", context_for(GOOD))
        assert result.outcome is RuleOutcome.PASS
        assert "42.1" in result.detail

    def test_a_gc_rich_target_fails_and_names_the_bound_it_crossed(self) -> None:
        result = check("gc-content", context_for("GCGCGCGCGCGCGCGCGCG"))
        assert result.outcome is RuleOutcome.FAIL
        assert "60" in result.detail

    def test_a_gc_poor_target_fails(self) -> None:
        result = check("gc-content", context_for("ATATATATATATATATATA"))
        assert result.outcome is RuleOutcome.FAIL
        assert "40" in result.detail

    def test_the_legacy_band_admits_a_target_the_paper_rejects(self) -> None:
        # The single parameter difference that accounts for much of the gap between presets.
        target = "GCGCGCATGCGTCAGTCAT"  # 63.2 % GC: above the paper's 60 %, below the legacy 65 %
        strict = check("gc-content", context_for(target))
        legacy = check("gc-content", context_for(target), **{"maximum-gc-fraction": 0.65})
        assert strict.outcome is RuleOutcome.FAIL
        assert legacy.outcome is RuleOutcome.PASS

    def test_the_measurement_is_a_fraction_not_a_rounded_percentage(self) -> None:
        result = check("gc-content", context_for(GOOD))
        assert result.measurements["gc-fraction"] == pytest.approx(8 / 19)


class TestHomopolymerRun:
    def test_a_target_with_no_long_run_passes(self) -> None:
        assert check("no-homopolymer-run", context_for(GOOD)).outcome is RuleOutcome.PASS

    def test_a_run_of_five_fails(self) -> None:
        result = check("no-homopolymer-run", context_for("GCTTAAAAAGGATTCAAGT"))
        assert result.outcome is RuleOutcome.FAIL
        assert "5" in result.detail

    def test_a_run_of_four_passes_at_the_default(self) -> None:
        # The default permits runs up to 4, so 5 is the first rejected length.
        result = check("no-homopolymer-run", context_for("GCTTAAAAGGATTCAAGTT"))
        assert result.outcome is RuleOutcome.PASS


class TestCpgCount:
    def test_a_target_with_one_cpg_passes_at_the_default(self) -> None:
        result = check("cpg-count", context_for(GOOD))
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["cpg-count"] == 1

    def test_the_legacy_setting_rejects_the_same_target(self) -> None:
        # The original script tolerated no CpG at all, which is stricter than the paper.
        result = check("cpg-count", context_for(GOOD), **{"maximum-cpg-count": 0})
        assert result.outcome is RuleOutcome.FAIL

    def test_two_cpgs_fail_at_the_default(self) -> None:
        result = check("cpg-count", context_for("GCTTACGGATTCAACGTTA"))
        assert result.outcome is RuleOutcome.FAIL
        assert "2" in result.detail


class TestSimpleRepeats:
    def test_a_target_with_no_simple_repeat_passes(self) -> None:
        assert check("no-simple-repeats", context_for(GOOD)).outcome is RuleOutcome.PASS

    def test_a_dinucleotide_repeat_fails(self) -> None:
        result = check("no-simple-repeats", context_for("GCATATATATATTCAAGTT"))
        assert result.outcome is RuleOutcome.FAIL

    def test_a_trinucleotide_repeat_fails(self) -> None:
        result = check("no-simple-repeats", context_for("GCCAGCAGCAGCAGCAAGT"))
        assert result.outcome is RuleOutcome.FAIL

    def test_a_homopolymer_is_not_counted_as_a_simple_repeat(self) -> None:
        # Measured by the homopolymer rule instead; counting it here would penalise one feature
        # against two separate criteria.
        result = check("no-simple-repeats", context_for("GCTTAAAAAGGATTCAAGT"))
        assert result.outcome is RuleOutcome.PASS

    def test_the_detail_names_the_repeat_unit_length(self) -> None:
        result = check("no-simple-repeats", context_for("GCATATATATATTCAAGTT"))
        assert "dinucleotide" in result.detail


class TestPositionRules:
    def test_position_eighteen_a_or_t_passes(self) -> None:
        assert check("pos18-at", context_for(GOOD)).outcome is RuleOutcome.PASS

    def test_position_eighteen_g_fails(self) -> None:
        result = check("pos18-at", context_for("GCTTAACGGATTCAAGTGA"))
        assert result.outcome is RuleOutcome.FAIL
        assert "G" in result.detail

    def test_position_nineteen_a_passes(self) -> None:
        assert check("pos19-a", context_for("GCTTAACGGATTCAAGTTA")).outcome is RuleOutcome.PASS

    def test_position_nineteen_not_a_fails(self) -> None:
        result = check("pos19-a", context_for("GCTTAACGGATTCAAGTTG"))
        assert result.outcome is RuleOutcome.FAIL

    def test_a_target_too_short_to_have_the_position_is_not_applicable(self) -> None:
        # Not a failure: a 15-mer has no position 18 to measure, so the rule cannot be
        # evaluated and must say so rather than rejecting the candidate.
        result = check("pos18-at", context_for(GOOD[:15]))
        assert result.outcome is RuleOutcome.NOT_APPLICABLE
        assert "18" in result.detail

    def test_start_position_one_gc_is_a_preference_worth_a_point(self) -> None:
        result = check("start-1-gc", context_for(GOOD))
        assert result.outcome is RuleOutcome.PASS
        assert result.kind is RuleKind.PREFERENCE
        assert result.points == 1

    def test_start_position_one_not_gc_scores_nothing_but_does_not_exclude(self) -> None:
        result = check("start-1-gc", context_for("ACTTAACGGATTCAAGTTA"))
        assert result.outcome is RuleOutcome.FAIL
        assert result.points == 0
        assert result.excludes is False

    def test_start_position_two_gc(self) -> None:
        assert check("start-2-gc", context_for(GOOD)).outcome is RuleOutcome.PASS

    def test_the_exploratory_position_seven_rule_reads_position_seven(self) -> None:
        # GOOD position 7 is G, so this does not fire; the rule is off by default anyway.
        assert check("pos7-u", context_for(GOOD)).outcome is RuleOutcome.FAIL
        assert check("pos7-u", context_for("GCTTAATGGATTCAAGTTA")).outcome is RuleOutcome.PASS


class TestFlankPreferences:
    def test_an_a_or_t_flank_position_passes(self) -> None:
        result = check("flank-20-au", context_for(GOOD, flank="TTTT"))
        assert result.outcome is RuleOutcome.PASS

    def test_a_g_or_c_flank_position_does_not_score(self) -> None:
        result = check("flank-20-au", context_for(GOOD, flank="GTTT"))
        assert result.outcome is RuleOutcome.FAIL
        assert result.points == 0

    @pytest.mark.parametrize(
        ("rule_id", "flank"),
        [
            ("flank-21-au", "GTGG"),
            ("flank-22-au", "GGTG"),
            ("flank-23-au", "GGGT"),
        ],
    )
    def test_each_flank_position_reads_its_own_base(self, rule_id: str, flank: str) -> None:
        assert check(rule_id, context_for(GOOD, flank=flank)).outcome is RuleOutcome.PASS

    def test_a_promoter_ending_at_the_target_makes_the_flank_not_applicable(self) -> None:
        # The flank lies in the promoter, not the duplex, so a target at the promoter's 3' edge
        # has nothing to measure. Scoring it zero would silently penalise the last candidates
        # of every promoter.
        result = check("flank-20-au", context_for(GOOD, flank=""))
        assert result.outcome is RuleOutcome.NOT_APPLICABLE
        assert result.points == 0

    def test_a_truncated_flank_is_partly_evaluable(self) -> None:
        context = context_for(GOOD, flank="TT")
        assert check("flank-20-au", context).outcome is RuleOutcome.PASS
        assert check("flank-22-au", context).outcome is RuleOutcome.NOT_APPLICABLE


class TestTssWindow:
    def test_a_target_inside_the_window_passes(self) -> None:
        result = check("tss-window", context_for(GOOD, tss_start=-500))
        assert result.outcome is RuleOutcome.PASS
        assert "-500" in result.detail

    def test_a_target_outside_the_window_fails(self) -> None:
        result = check("tss-window", context_for(GOOD, tss_start=-50))
        assert result.outcome is RuleOutcome.FAIL

    def test_a_target_extending_past_the_window_fails(self) -> None:
        # ADR 0008: total containment. A 19-mer from -105 runs to -87.
        result = check("tss-window", context_for(GOOD, tss_start=-105))
        assert result.outcome is RuleOutcome.FAIL
        assert "-87" in result.detail

    def test_a_promoter_with_no_tss_is_not_applicable_and_does_not_exclude(self) -> None:
        # The Workflow B case. Failing here would empty the results of every pasted-sequence run.
        result = check("tss-window", context_for(GOOD))
        assert result.outcome is RuleOutcome.NOT_APPLICABLE
        assert result.excludes is False
        assert "TSS" in result.detail

    def test_the_preferred_window_is_a_preference_inside_the_mandatory_one(self) -> None:
        inside = check("tss-preferred-window", context_for(GOOD, tss_start=-400))
        outside = check("tss-preferred-window", context_for(GOOD, tss_start=-900))
        assert inside.outcome is RuleOutcome.PASS
        assert inside.points == 1
        assert outside.outcome is RuleOutcome.FAIL
        assert outside.excludes is False

    def test_the_window_bounds_are_configurable(self) -> None:
        result = check(
            "tss-window",
            context_for(GOOD, tss_start=-50),
            **{"window-start-bp": -1000, "window-end-bp": -10},
        )
        assert result.outcome is RuleOutcome.PASS


class TestCpgIsland:
    ISLAND_PROMOTER = "CG" * 300

    def _island_context(self, offset: int) -> CandidateContext:
        sequence = NucleicSequence.parse(self.ISLAND_PROMOTER + "AT" * 300)
        islands = find_cpg_islands(sequence, GARDINER_GARDEN_FROMMER_1987)
        promoter = Promoter(sequence=sequence)
        return CandidateContext.at(promoter, offset=offset, length=19, cpg_islands=islands)

    def test_a_target_inside_an_island_fails(self) -> None:
        result = check("not-in-cpg-island", self._island_context(offset=0))
        assert result.outcome is RuleOutcome.FAIL
        assert result.excludes is True

    def test_a_target_outside_every_island_passes(self) -> None:
        result = check("not-in-cpg-island", self._island_context(offset=1000))
        assert result.outcome is RuleOutcome.PASS

    def test_not_running_detection_is_not_applicable(self) -> None:
        # Distinct from finding no islands. Passing here would claim a check nobody performed.
        result = check("not-in-cpg-island", context_for(GOOD))
        assert result.outcome is RuleOutcome.NOT_APPLICABLE
        assert result.excludes is False

    def test_finding_no_islands_is_a_pass(self) -> None:
        promoter = Promoter(sequence=NucleicSequence.parse(GOOD + "TTTT"))
        context = CandidateContext.at(promoter, offset=0, length=19, cpg_islands=())
        result = check("not-in-cpg-island", context)
        assert result.outcome is RuleOutcome.PASS


class TestThermodynamicAsymmetry:
    # 5' end GCGC is strongly stacked, 3' end TTAA weakly, so the antisense 5' end is weaker.
    ANTISENSE_FAVOURED = "GCGCACGGATCCAAGTTAA"
    SENSE_FAVOURED = "TTAAACGGATCCAAGGCGC"
    # 3' end UUAG stacks to -4.34 against the same 5' GCGC at -9.20, so the margin is only
    # +4.86 kcal/mol — enough to pass the default, not enough for a raised threshold.
    WEAKLY_ASYMMETRIC = "GCGCACGGATCCAAGTTAG"

    def test_an_asymmetric_target_favouring_the_antisense_strand_passes(self) -> None:
        result = check("thermodynamic-asymmetry", context_for(self.ANTISENSE_FAVOURED))
        assert result.outcome is RuleOutcome.PASS
        assert result.measurements["predicted-guide"] == "antisense"

    def test_a_target_favouring_the_sense_strand_fails(self) -> None:
        result = check("thermodynamic-asymmetry", context_for(self.SENSE_FAVOURED))
        assert result.outcome is RuleOutcome.FAIL
        assert result.measurements["predicted-guide"] == "sense"

    def test_the_reason_is_the_full_explanation_not_a_bare_number(self) -> None:
        result = check("thermodynamic-asymmetry", context_for(self.ANTISENSE_FAVOURED))
        assert "Schwarz" in result.detail
        assert "antisense" in result.detail

    def test_the_result_records_which_table_was_used(self) -> None:
        # Required for reproducibility: the verdict is meaningless without it (ADR 0009).
        result = check("thermodynamic-asymmetry", context_for(self.ANTISENSE_FAVOURED))
        assert result.measurements["nearest-neighbour-table"] == "rna-xia-1998"

    def test_the_legacy_table_can_be_selected(self) -> None:
        result = check(
            "thermodynamic-asymmetry",
            context_for(self.ANTISENSE_FAVOURED),
            **{"nearest-neighbour-table": "legacy-dna-kj"},
        )
        assert result.measurements["nearest-neighbour-table"] == "legacy-dna-kj"
        assert result.measurements["asymmetry-unit"] == "kJ/mol"

    def test_a_raised_threshold_rejects_a_weakly_asymmetric_target(self) -> None:
        context = context_for(self.WEAKLY_ASYMMETRIC)
        assert check("thermodynamic-asymmetry", context).outcome is RuleOutcome.PASS
        raised = check(
            "thermodynamic-asymmetry", context, **{"minimum-asymmetry-kcal": 5.0}
        )
        assert raised.outcome is RuleOutcome.FAIL
        assert "below the" in raised.detail

    def test_the_threshold_cannot_be_set_outside_the_parameters_own_range(self) -> None:
        # The bound belongs to the parameter, so an absurd threshold is rejected before any
        # evaluation happens rather than producing an empty candidate set.
        with pytest.raises(ValueError, match="above the maximum"):
            check(
                "thermodynamic-asymmetry",
                context_for(self.ANTISENSE_FAVOURED),
                **{"minimum-asymmetry-kcal": 7.0},
            )

    def test_the_threshold_is_converted_into_the_tables_unit(self) -> None:
        """A kcal/mol threshold must not be compared against a kJ/mol asymmetry.

        The parameter is stated in kcal/mol; the legacy table reports kJ/mol. The error is
        invisible at the default threshold of 0.0, where only the sign matters, and silently
        wrong at any other value.

        This target's legacy asymmetry is +16.77 kJ/mol, which is +4.01 kcal/mol. A threshold of
        5.0 kcal/mol is 20.92 kJ/mol, so the candidate must be rejected. Compared unconverted,
        5.0 against 16.77 would wrongly pass — so this threshold distinguishes the two.
        """
        context = context_for(self.ANTISENSE_FAVOURED)
        legacy = check(
            "thermodynamic-asymmetry",
            context,
            **{"nearest-neighbour-table": "legacy-dna-kj"},
        )
        assert legacy.measurements["asymmetry"] == pytest.approx(16.77)

        result = check(
            "thermodynamic-asymmetry",
            context,
            **{"nearest-neighbour-table": "legacy-dna-kj", "minimum-asymmetry-kcal": 5.0},
        )
        assert result.outcome is RuleOutcome.FAIL
        # The reported threshold is in the table's unit, so the message is self-consistent.
        assert "20.92 kJ/mol required" in result.detail

    def test_a_threshold_just_below_the_converted_margin_still_passes(self) -> None:
        # The other side of the same boundary: 3.0 kcal/mol is 12.55 kJ/mol, below the 16.77
        # the target actually has.
        result = check(
            "thermodynamic-asymmetry",
            context_for(self.ANTISENSE_FAVOURED),
            **{"nearest-neighbour-table": "legacy-dna-kj", "minimum-asymmetry-kcal": 3.0},
        )
        assert result.outcome is RuleOutcome.PASS

    def test_the_tiebreak_reports_an_ordering_value_and_scores_nothing(self) -> None:
        result = check("asymmetry-margin", context_for(self.ANTISENSE_FAVOURED))
        assert result.kind is RuleKind.TIEBREAK
        assert result.ordering_value is not None
        assert result.points == 0
        assert result.excludes is False

    def test_the_tiebreak_records_its_table(self) -> None:
        result = check("asymmetry-margin", context_for(self.ANTISENSE_FAVOURED))
        assert result.measurements["nearest-neighbour-table"] == "rna-xia-1998"


class TestEveryResultExplainsItself:
    def test_no_evaluator_returns_a_blank_reason(self) -> None:
        # Enforced by RuleResult's own validation, but asserted across every rule so that a new
        # evaluator cannot be added without a reason.
        context = context_for(GOOD, tss_start=-400, cpg_islands=())
        for rule in CATALOGUE:
            result = evaluate(context, ResolvedRule(rule))
            assert result.detail.strip(), f"{rule.id} gave no reason"
            assert result.rule_id == rule.id


class TestEvaluatorLookupFailsLoudly:
    def test_asking_for_an_unregistered_rule_raises(self) -> None:
        # A rule with no evaluator would never fire while still appearing in the settings UI,
        # so this must not return a no-op.
        with pytest.raises(KeyError, match="no evaluator is registered"):
            evaluator_for("a-rule-that-does-not-exist")

    def test_the_error_says_what_to_do_about_it(self) -> None:
        with pytest.raises(KeyError, match="rather than leaving the rule inert"):
            evaluator_for("another-missing-rule")
