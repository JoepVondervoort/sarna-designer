"""Evaluators for rules about the base at a particular position.

Positions are Table 1's positions, counted from 1 at the target's 5' end, and continuing past
the target into the promoter for the 20-23 flank. The context owns that numbering, so nothing
here does index arithmetic.

**A known ambiguity, stated rather than hidden.** Table 1's position rules are written for a
19 nt target: position 18 is the second-to-last base and position 19 the last. Because target
length is configurable (18-21 nt), "position 18" could mean the eighteenth base or the
second-to-last. These evaluators use the **literal** position, which is what Table 1 says, and
report not-applicable where the target is too short to have it. For a 21-mer that means position
19 is no longer the 3' terminal base, so the thermodynamic reasoning behind ``pos19-a`` — an A
there puts a U at the guide's 5' end — no longer holds. Changing target length away from 19
therefore changes what these rules mean, and the project has not decided which reading is
intended. Until it does, leave the length at 19.
"""

from collections.abc import Callable

from sarna.design.context import CandidateContext
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult

Evaluator = Callable[[CandidateContext, ResolvedRule], RuleResult]

PURINE_FREE_BASES = "AT"
"""A or T on a DNA target, equivalently A or U on the RNA duplex."""

STRONG_BASES = "GC"


def _target_base(
    context: CandidateContext, resolved: ResolvedRule, position: int
) -> tuple[str | None, RuleResult | None]:
    """Return the base at ``position``, or a not-applicable result if the target is too short."""
    if position > context.length:
        return None, RuleResult.not_applicable(
            resolved.rule,
            f"The target is only {context.length} nt long, so it has no position {position} "
            f"to measure.",
        )
    return context.base_at(position), None


def _at_target_position(
    context: CandidateContext,
    resolved: ResolvedRule,
    position: int,
    permitted: str,
    description: str,
) -> RuleResult:
    """Evaluate a rule requiring one of ``permitted`` at a target position."""
    base, unevaluable = _target_base(context, resolved, position)
    if unevaluable is not None:
        return unevaluable
    assert base is not None
    measurements = {f"position-{position}": base}
    if base in permitted:
        return RuleResult.passed(
            resolved.rule, f"Position {position} is {base}, which is {description}.", measurements
        )
    return RuleResult.failed(
        resolved.rule,
        f"Position {position} is {base}, not {description}.",
        measurements,
    )


def pos18_at(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check position 18 is A or T, weakening the guide strand's 5' end."""
    return _at_target_position(context, resolved, 18, PURINE_FREE_BASES, "A or T")


def pos19_a(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check position 19 is A, placing a U at the guide strand's 5' end."""
    return _at_target_position(context, resolved, 19, "A", "A")


def pos7_u(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check position 7 is T, the exploratory rule carried from the original script.

    No source supports this rule; it is off by default and its provenance is marked unsourced.
    The target is DNA, so T here is the U the rule's name refers to.
    """
    return _at_target_position(context, resolved, 7, "T", "T (U on the duplex)")


def start_gc(position: int) -> Evaluator:
    """Build an evaluator checking for G or C at a target start position."""

    def evaluate(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
        return _at_target_position(context, resolved, position, STRONG_BASES, "G or C")

    evaluate.__doc__ = f"Check position {position} is G or C."
    return evaluate


def flank_au(position: int) -> Evaluator:
    """Build an evaluator checking for A or T at a flank position in the promoter."""

    def evaluate(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
        base = context.flank_base_at(position)
        if base is None:
            return RuleResult.not_applicable(
                resolved.rule,
                f"The promoter ends before position {position}, so the 3' flank cannot be "
                f"read there. The flank lies in the promoter rather than in the duplex.",
            )
        measurements = {f"position-{position}": base}
        if base in PURINE_FREE_BASES:
            return RuleResult.passed(
                resolved.rule,
                f"Flank position {position} is {base}, which is A or T.",
                measurements,
            )
        return RuleResult.failed(
            resolved.rule,
            f"Flank position {position} is {base}, not A or T.",
            measurements,
        )

    evaluate.__doc__ = f"Check flank position {position} is A or T."
    return evaluate


def requires_base(position: int, permitted: str, description: str) -> Evaluator:
    """Build an evaluator checking the base at ``position`` is one of ``permitted``.

    The general form behind the named evaluators above, used by the optional rules from
    Reynolds et al. 2004 and thesis chapter 4. ``description`` is the phrase that appears in the
    verdict, so it reads as prose: "Position 13 is G, not A, C or T".

    Args:
        position: 1-based position in the target.
        permitted: The bases that satisfy the rule, as DNA.
        description: How to describe the requirement to a reader.
    """

    def evaluate(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
        return _at_target_position(context, resolved, position, permitted, description)

    evaluate.__doc__ = f"Check position {position} is {description}."
    return evaluate
