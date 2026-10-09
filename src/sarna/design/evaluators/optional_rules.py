"""Evaluators for the optional siRNA-literature and chapter 4 rules.

All of these read the target site as DNA, so ``T`` is the base tested where a paper writes ``U``.
The papers state their criteria on the sense strand, which is the target, so positions need no
translation, but the *reasoning* often concerns the guide, and the two run in opposite
directions. Target position 19 is guide position 1. Comments here name whichever strand the
mechanism is about, because that is where mistakes of this kind hide.
"""

from sarna.design.context import CandidateContext
from sarna.rules.parameters import ParameterValue
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult

WEAK_BASES = "AT"
"""A or T on a DNA target; A or U on the RNA duplex. 'Weak' as in two hydrogen bonds."""

STRONG_BASES = "GC"


def au_rich_guide_5_prime(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target's 3' terminal window is A/T-rich, weakening the guide's 5' end.

    The window is counted from the target's 3' end, so it tracks the guide's 5' end as target
    length changes rather than drifting away from the region the criterion is about.
    """
    window = resolved.integer("window-length-nt")
    required = resolved.integer("minimum-au-count")
    if window > context.length:
        return RuleResult.not_applicable(
            resolved,
            f"The target is {context.length} nt long, shorter than the {window} nt window this "
            f"rule examines, so the window cannot be read.",
        )
    first = context.length - window + 1
    bases = [context.base_at(position) for position in range(first, context.length + 1)]
    weak = sum(1 for base in bases if base in WEAK_BASES)
    measurements: dict[str, ParameterValue] = {
        "window-start": first,
        "window-end": context.length,
        "au-count": weak,
        "window-bases": "".join(bases),
    }
    region = f"positions {first}-{context.length}"
    if weak >= required:
        return RuleResult.passed(
            resolved,
            f"{weak} of the {window} bases at {region} are A or T, meeting the {required} "
            f"required. A weak 3' end on the target is a weak 5' end on the guide, which is "
            f"what biases AGO2 to retain the intended strand.",
            measurements,
        )
    return RuleResult.failed(
        resolved,
        f"Only {weak} of the {window} bases at {region} are A or T, below the {required} "
        f"required, so the guide's 5' end is more stable than asymmetric loading prefers.",
        measurements,
    )


def no_gc_stretch(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target contains no G/C run longer than the permitted length.

    Counts G and C together, in any order, which is what distinguishes this from the mandatory
    homopolymer rule: ``GCGCGCGCGC`` is ten bases of G/C and no repeat of a single base.
    """
    permitted = resolved.integer("maximum-gc-stretch-nt")
    longest, running = 0, 0
    for position in range(1, context.length + 1):
        running = running + 1 if context.base_at(position) in STRONG_BASES else 0
        longest = max(longest, running)
    measurements: dict[str, ParameterValue] = {
        "longest-gc-stretch": longest,
        "maximum-permitted": permitted,
    }
    if longest <= permitted:
        return RuleResult.passed(
            resolved,
            f"The longest run of G or C is {longest} nt, within the {permitted} nt permitted.",
            measurements,
        )
    return RuleResult.failed(
        resolved,
        f"The longest run of G or C is {longest} nt, above the {permitted} nt permitted. A "
        f"tract that stable resists unwinding and favours intramolecular structure.",
        measurements,
    )


def central_stability(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target's central region carries enough G/C bases.

    Chapter 4's stability finding, carried as a hypothesis. The region is given as literal
    positions rather than fractions of the length, because that is how the chapter reports it.
    """
    start = resolved.integer("central-region-start")
    end = resolved.integer("central-region-end")
    required = resolved.integer("minimum-strong-count")
    if start > end:
        raise ValueError(
            f"central-region-start ({start}) is after central-region-end ({end}); the region "
            f"would be empty, so no candidate could ever satisfy the rule"
        )
    if end > context.length:
        return RuleResult.not_applicable(
            resolved,
            f"The target is only {context.length} nt long, so it has no position {end} and the "
            f"central region {start}-{end} cannot be read.",
        )
    bases = [context.base_at(position) for position in range(start, end + 1)]
    strong = sum(1 for base in bases if base in STRONG_BASES)
    measurements: dict[str, ParameterValue] = {
        "region-start": start,
        "region-end": end,
        "strong-count": strong,
        "region-bases": "".join(bases),
    }
    if strong >= required:
        return RuleResult.passed(
            resolved,
            f"{strong} of the {len(bases)} bases at positions {start}-{end} are G or C, "
            f"meeting the {required} required.",
            measurements,
        )
    return RuleResult.failed(
        resolved,
        f"Only {strong} of the {len(bases)} bases at positions {start}-{end} are G or C, "
        f"below the {required} required.",
        measurements,
    )
