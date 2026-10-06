"""The evaluator lookup: which function judges which rule.

Evaluation is kept out of the rule catalogue (ADR 0003, ADR 0007), so the two are joined here by
rule id. The arrangement has one property worth the indirection: importing
:mod:`sarna.rules.catalog` to render provenance in a UI pulls in no evaluation code at all.

It has one hazard to match, which :func:`evaluate` closes: a rule with no evaluator would never
fire, and a rule that never fires is worse than a missing rule, because it looks like it is
protecting you. Asking for one raises, and a test asserts the catalogue and this table agree.
"""

from collections.abc import Callable, Mapping

from sarna.design.context import CandidateContext
from sarna.design.evaluators import (
    location_rules,
    position_rules,
    sequence_rules,
    thermodynamic_rules,
)
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult

Evaluator = Callable[[CandidateContext, ResolvedRule], RuleResult]
"""The evaluator contract fixed by ADR 0007: a pure function, context first."""


EVALUATORS: Mapping[str, Evaluator] = {
    # Sequence composition.
    "target-length-19": sequence_rules.target_length,
    "gc-content": sequence_rules.gc_content,
    "no-homopolymer-run": sequence_rules.no_homopolymer_run,
    "cpg-count": sequence_rules.cpg_count,
    "no-simple-repeats": sequence_rules.no_simple_repeats,
    # Position identity.
    "pos18-at": position_rules.pos18_at,
    "pos19-a": position_rules.pos19_a,
    "pos7-u": position_rules.pos7_u,
    "start-1-gc": position_rules.start_gc(1),
    "start-2-gc": position_rules.start_gc(2),
    "flank-20-au": position_rules.flank_au(20),
    "flank-21-au": position_rules.flank_au(21),
    "flank-22-au": position_rules.flank_au(22),
    "flank-23-au": position_rules.flank_au(23),
    # Location, which a sequence-only run cannot answer.
    "tss-window": location_rules.tss_window,
    "tss-preferred-window": location_rules.tss_preferred_window,
    "not-in-cpg-island": location_rules.not_in_cpg_island,
    # Duplex thermodynamics.
    "thermodynamic-asymmetry": thermodynamic_rules.thermodynamic_asymmetry,
    "asymmetry-margin": thermodynamic_rules.asymmetry_margin,
}


def evaluator_for(rule_id: str) -> Evaluator:
    """Return the evaluator for ``rule_id``.

    Raises:
        KeyError: If no evaluator is registered. This is a programming error rather than a
            configuration one: a catalogued rule without an evaluator would silently never fire.
    """
    try:
        return EVALUATORS[rule_id]
    except KeyError:
        raise KeyError(
            f"no evaluator is registered for rule {rule_id!r}, so it would never be applied; "
            f"register one in sarna.design.evaluators rather than leaving the rule inert"
        ) from None


def evaluate(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Evaluate one resolved rule against one candidate.

    Args:
        context: The candidate being judged.
        resolved: The rule and the parameter values this run uses.

    Returns:
        The rule's verdict, with its reason and measurements.

    Raises:
        KeyError: If the rule has no registered evaluator.
    """
    return evaluator_for(resolved.rule.id)(context, resolved)


__all__ = ["EVALUATORS", "Evaluator", "evaluate", "evaluator_for"]
