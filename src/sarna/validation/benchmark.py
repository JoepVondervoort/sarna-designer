"""Run the design engine against measured activity, and report how well it discriminates.

This is the harness that can embarrass a rule. Everything else in the project asks whether a rule
was implemented as published; this asks whether candidates the engine prefers were more likely to
work.

Two structural limits shape what any answer means, and both are properties of the data rather
than of the engine.

**The gate cannot be assessed.** Every sequence in the dataset was designed to pass Wang's
structural criteria, so the set contains almost no candidate the gate would reject. A gate that
admits everything presented to it has not been shown to be permissive; it has been shown nothing.

**Assayed sequences arrive without their promoters.** The dataset gives target sites, not genomic
coordinates, so the rules that need promoter context — the TSS windows, CpG-island overlap, the
3' flank positions — report not-applicable rather than contributing. The measured discrimination
is therefore the discrimination of the *sequence-intrinsic* rules alone, which is a subset of
what a real run applies. Recovering the coordinates would fix this and is worth doing; it needs
the promoter each duplex was designed against, which the chapter does not tabulate.

The headline statistic is the probability that a randomly chosen activator outranks a randomly
chosen non-activator, with ties counted as half. That is the area under the ROC curve, and it is
reported as a probability because "0.5 means the ranking is uninformative" is a sentence a
reader can check against their intuition, which "p = 0.3" is not.
"""

import random
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from sarna.design.candidate import Candidate
from sarna.design.config import DesignConfig
from sarna.design.context import CandidateContext
from sarna.design.evaluators import evaluate
from sarna.design.pipeline import DEFAULT_COMPARISON_TABLES
from sarna.genome.promoter import Promoter
from sarna.rules.results import RuleOutcome
from sarna.thermo.guide import compare_asymmetry
from sarna.validation.dataset import ValidationDataset, ValidationSequence

DEFAULT_PERMUTATIONS = 10_000
DEFAULT_SEED = 20261006

UNINFORMATIVE_AUC = 0.5
"""The value a ranking with no information takes. Stated so comparisons read unambiguously."""


@dataclass(frozen=True, slots=True)
class ScoredSequence:
    """One assayed sequence and the verdict the engine reached about it."""

    sequence: ValidationSequence
    candidate: Candidate

    @property
    def is_activating(self) -> bool:
        """Return whether the assay saw activation."""
        return self.sequence.is_activating

    @property
    def points(self) -> int:
        """Return the preference points the engine awarded."""
        return self.candidate.preference_points

    @property
    def was_rejected(self) -> bool:
        """Return whether a mandatory rule excluded it."""
        return self.candidate.is_rejected


@dataclass(frozen=True, slots=True)
class RuleDiscrimination:
    """How one rule's verdicts divided along the activity labels."""

    rule_id: str
    activating_pass_percent: float
    not_activating_pass_percent: float
    unevaluated_count: int
    """How many sequences the rule could not be evaluated for at all."""

    @property
    def difference(self) -> float:
        """Return the pass-rate difference in percentage points, activators minus the rest."""
        return self.activating_pass_percent - self.not_activating_pass_percent

    @property
    def was_never_evaluated(self) -> bool:
        """Return whether the rule contributed nothing, for want of promoter context."""
        return self.activating_pass_percent == 0.0 and self.not_activating_pass_percent == 0.0

    def __str__(self) -> str:
        """Return a one-line rendering of the rule's separation."""
        if self.was_never_evaluated:
            return f"{self.rule_id}: not evaluable without promoter context"
        return (
            f"{self.rule_id}: {self.activating_pass_percent:.0f}% of activators pass vs "
            f"{self.not_activating_pass_percent:.0f}% of the rest ({self.difference:+.0f} pp)"
        )


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """What one preset achieved against one dataset."""

    preset_id: str | None
    scored: tuple[ScoredSequence, ...]
    auc: float
    """Probability an activator outranks a non-activator on preference points, ties as half."""

    auc_p_value: float
    """One-sided permutation p-value for the AUC exceeding chance."""

    permutations: int
    seed: int
    rule_discrimination: tuple[RuleDiscrimination, ...]
    """Per-rule pass rates in each class, largest separation first."""

    max_preference_points: int

    @property
    def activating(self) -> tuple[ScoredSequence, ...]:
        """Return the scored activators."""
        return tuple(s for s in self.scored if s.is_activating)

    @property
    def not_activating(self) -> tuple[ScoredSequence, ...]:
        """Return the scored non-activators."""
        return tuple(s for s in self.scored if not s.is_activating)

    @property
    def rejected_count(self) -> int:
        """Return how many sequences a mandatory rule excluded."""
        return sum(1 for s in self.scored if s.was_rejected)

    @property
    def rejected_activators(self) -> tuple[ScoredSequence, ...]:
        """Return activators the gate excluded.

        The most diagnostic output here. A duplex that demonstrably worked and that the engine
        would have discarded is direct evidence against a mandatory rule, and it is the one
        finding this range-restricted dataset *can* deliver.
        """
        return tuple(s for s in self.activating if s.was_rejected)

    @property
    def mean_points_activating(self) -> float:
        """Return mean preference points among activators."""
        return _mean([s.points for s in self.activating])

    @property
    def mean_points_not_activating(self) -> float:
        """Return mean preference points among non-activators."""
        return _mean([s.points for s in self.not_activating])

    @property
    def discriminates(self) -> bool:
        """Return whether the ranking beat chance at the 0.05 level."""
        return self.auc_p_value < 0.05

    @property
    def summary(self) -> str:
        """Return a plain-language verdict on the ranking's usefulness."""
        direction = "above" if self.auc > UNINFORMATIVE_AUC else "below"
        if self.discriminates:
            return (
                f"Ranking carries signal: an activator outranks a non-activator "
                f"{self.auc:.0%} of the time (p = {self.auc_p_value:.3f})."
            )
        return (
            f"Ranking is uninformative on this set: an activator outranks a non-activator "
            f"{self.auc:.0%} of the time, {direction} the {UNINFORMATIVE_AUC:.0%} expected from "
            f"a coin flip, p = {self.auc_p_value:.3f}."
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the result as JSON-serialisable data."""
        return {
            "preset": self.preset_id,
            "counts": {
                "activating": len(self.activating),
                "not_activating": len(self.not_activating),
                "rejected_by_a_gate": self.rejected_count,
                "activators_rejected_by_a_gate": len(self.rejected_activators),
            },
            "preference_points": {
                "maximum": self.max_preference_points,
                "mean_activating": round(self.mean_points_activating, 2),
                "mean_not_activating": round(self.mean_points_not_activating, 2),
            },
            "auc": round(self.auc, 4),
            "auc_p_value": round(self.auc_p_value, 4),
            "permutations": self.permutations,
            "seed": self.seed,
            "discriminates": self.discriminates,
            "summary": self.summary,
            "rules": [
                {
                    "id": rule.rule_id,
                    "activating_pass_percent": round(rule.activating_pass_percent, 1),
                    "not_activating_pass_percent": round(rule.not_activating_pass_percent, 1),
                    "difference_pp": round(rule.difference, 1),
                    "unevaluated": rule.unevaluated_count,
                }
                for rule in self.rule_discrimination
            ],
        }


def score_sequence(sequence: ValidationSequence, config: DesignConfig) -> ScoredSequence:
    """Evaluate one assayed sequence as though it were a candidate.

    The sequence is wrapped in a promoter consisting of nothing but itself, with no TSS and no
    coordinates. Rules needing that context report not-applicable, which is the correct outcome
    rather than a failure: the dataset does not record where these targets came from.
    """
    promoter = Promoter(sequence=sequence.target_dna)
    context = CandidateContext.at(promoter, offset=0, length=len(sequence.target))
    duplex = context.target.as_rna()
    candidate = Candidate(
        id=sequence.name,
        offset=0,
        target=context.target,
        sense_strand=duplex,
        antisense_strand=duplex.reverse_complement(),
        tss_span=None,
        interval=None,
        results=tuple(evaluate(context, resolved) for resolved in config),
        guide_comparison=compare_asymmetry(context.target, DEFAULT_COMPARISON_TABLES),
        max_preference_points=config.max_preference_points,
    )
    return ScoredSequence(sequence=sequence, candidate=candidate)


def benchmark(
    dataset: ValidationDataset,
    config: DesignConfig,
    permutations: int = DEFAULT_PERMUTATIONS,
    seed: int = DEFAULT_SEED,
) -> BenchmarkResult:
    """Score every sequence in ``dataset`` and measure how well the ranking separates activity.

    Args:
        dataset: The labelled sequences.
        config: The resolved rules to apply, normally from a preset.
        permutations: Label shuffles for the AUC's null distribution.
        seed: Random seed, carried into the result so the p-value is reproducible.

    Returns:
        The scored sequences, the ranking's AUC with its p-value, and per-rule pass rates.

    Raises:
        ValueError: If the dataset has no activators or no non-activators, leaving nothing to
            discriminate between.
    """
    scored = tuple(score_sequence(sequence, config) for sequence in dataset)
    activating = [s.points for s in scored if s.is_activating]
    inactive = [s.points for s in scored if not s.is_activating]
    if not activating or not inactive:
        raise ValueError(
            "both classes must be present to measure discrimination; this dataset has "
            f"{len(activating)} activating and {len(inactive)} non-activating sequences"
        )

    auc = _auc(activating, inactive)
    pooled = activating + inactive
    rng = random.Random(seed)
    at_least_as_high = 0
    for _ in range(permutations):
        shuffled = pooled[:]
        rng.shuffle(shuffled)
        if _auc(shuffled[: len(activating)], shuffled[len(activating) :]) >= auc:
            at_least_as_high += 1

    return BenchmarkResult(
        preset_id=config.preset_id,
        scored=scored,
        auc=auc,
        auc_p_value=(at_least_as_high + 1) / (permutations + 1),
        permutations=permutations,
        seed=seed,
        rule_discrimination=_rule_discrimination(scored, config),
        max_preference_points=config.max_preference_points,
    )


def _auc(group: Sequence[float], reference: Sequence[float]) -> float:
    """Return the probability a ``group`` value exceeds a ``reference`` value, ties as half.

    Computed by direct enumeration rather than through a rank formula. The dataset is small
    enough that the clarity is worth more than the speed, and a rank-based implementation is
    where tie handling usually goes quietly wrong.
    """
    wins = 0.0
    for high in group:
        for low in reference:
            if high > low:
                wins += 1.0
            elif high == low:
                wins += 0.5
    return wins / (len(group) * len(reference))


def _rule_discrimination(
    scored: Sequence[ScoredSequence], config: DesignConfig
) -> tuple[RuleDiscrimination, ...]:
    """Return each rule's pass rate in both classes, largest separation first."""
    activating = [s for s in scored if s.is_activating]
    inactive = [s for s in scored if not s.is_activating]
    discriminations = []
    for resolved in config:
        rule_id = resolved.rule.id
        unevaluated = sum(
            1 for s in scored if s.candidate.result(rule_id).is_unevaluated
        )
        discriminations.append(
            RuleDiscrimination(
                rule_id=rule_id,
                activating_pass_percent=_pass_percent(activating, rule_id),
                not_activating_pass_percent=_pass_percent(inactive, rule_id),
                unevaluated_count=unevaluated,
            )
        )
    return tuple(sorted(discriminations, key=lambda d: -abs(d.difference)))


def _pass_percent(scored: Sequence[ScoredSequence], rule_id: str) -> float:
    """Return the percentage of ``scored`` whose result for ``rule_id`` passed.

    Never called with an empty group: :func:`benchmark` refuses a dataset missing either class.
    """
    passed = sum(
        1 for s in scored if s.candidate.result(rule_id).outcome is RuleOutcome.PASS
    )
    return 100.0 * passed / len(scored)


def _mean(values: Sequence[float]) -> float:
    """Return the arithmetic mean of a non-empty sequence."""
    return sum(values) / len(values)


__all__ = [
    "DEFAULT_PERMUTATIONS",
    "DEFAULT_SEED",
    "UNINFORMATIVE_AUC",
    "BenchmarkResult",
    "RuleDiscrimination",
    "ScoredSequence",
    "benchmark",
    "score_sequence",
]
