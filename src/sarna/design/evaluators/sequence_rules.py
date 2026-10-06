"""Evaluators for rules that read only the target sequence.

Each is a pure function of a context and a resolved rule (ADR 0007). None of them touch
coordinates, so none can report not-applicable: a target always has a length, a GC content and a
CpG count.

A note on thresholds: the measurement lives in :mod:`sarna.sequence.nucleic` and the comparison
lives here. That split is what stopped defect B5 recurring — the legacy ``checkconsecutive``
hardcoded a run length that disagreed with its own documentation, because the measuring and the
judging were the same function.
"""

from sarna.design.context import CandidateContext
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult

REPEAT_UNIT_LENGTHS = (2, 3)
"""Unit lengths Table 1's simple-repeat criterion covers: di- and tri-nucleotide repeats."""

_UNIT_NAMES = {2: "dinucleotide", 3: "trinucleotide"}


def target_length(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target is the configured length."""
    required = resolved.integer("target-length")
    measurements = {"length": context.length}
    if context.length == required:
        return RuleResult.passed(
            resolved.rule,
            f"The target is {context.length} nt, as required.",
            measurements,
        )
    return RuleResult.failed(
        resolved.rule,
        f"The target is {context.length} nt, but {required} nt is required.",
        measurements,
    )


def gc_content(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target's GC content lies inside the configured band."""
    minimum = resolved.number("minimum-gc-fraction")
    maximum = resolved.number("maximum-gc-fraction")
    fraction = context.target.gc_fraction
    measurements = {"gc-fraction": fraction}
    band = f"{minimum:.0%} to {maximum:.0%} band"
    if fraction < minimum:
        return RuleResult.failed(
            resolved.rule,
            f"GC content is {fraction:.1%}, below the {minimum:.0%} minimum.",
            measurements,
        )
    if fraction > maximum:
        return RuleResult.failed(
            resolved.rule,
            f"GC content is {fraction:.1%}, above the {maximum:.0%} maximum.",
            measurements,
        )
    return RuleResult.passed(
        resolved.rule, f"GC content is {fraction:.1%}, inside the {band}.", measurements
    )


def no_homopolymer_run(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target has no run of identical bases longer than permitted."""
    permitted = resolved.integer("maximum-run-length")
    longest = context.target.longest_homopolymer_run
    measurements = {"longest-run": longest}
    if longest > permitted:
        return RuleResult.failed(
            resolved.rule,
            f"The target contains a run of {longest} identical bases; at most {permitted} "
            f"are permitted.",
            measurements,
        )
    return RuleResult.passed(
        resolved.rule,
        f"The longest run of identical bases is {longest}, within the limit of {permitted}.",
        measurements,
    )


def cpg_count(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target contains no more CpG dinucleotides than permitted."""
    permitted = resolved.integer("maximum-cpg-count")
    found = context.target.cpg_count
    measurements = {"cpg-count": found}
    plural = "" if found == 1 else "s"
    if found > permitted:
        return RuleResult.failed(
            resolved.rule,
            f"The target contains {found} CpG{plural}; at most {permitted} is permitted.",
            measurements,
        )
    return RuleResult.passed(
        resolved.rule,
        f"The target contains {found} CpG{plural}, within the limit of {permitted}.",
        measurements,
    )


def no_simple_repeats(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target is not a di- or tri-nucleotide simple repeat.

    A homopolymer is deliberately not counted: :meth:`NucleicSequence.longest_tandem_repeat`
    ignores units of identical bases, so ``AAAAAA`` is measured by the homopolymer rule alone and
    not penalised twice for one feature.
    """
    permitted = resolved.integer("maximum-tandem-repeat-units")
    copies = {
        unit: context.target.longest_tandem_repeat(unit) for unit in REPEAT_UNIT_LENGTHS
    }
    measurements = {f"tandem-copies-{unit}nt": count for unit, count in copies.items()}
    worst_unit, worst = max(copies.items(), key=lambda item: item[1])
    if worst > permitted:
        return RuleResult.failed(
            resolved.rule,
            f"The target contains {worst} tandem copies of a {_UNIT_NAMES[worst_unit]}; at "
            f"most {permitted} are permitted.",
            measurements,
        )
    return RuleResult.passed(
        resolved.rule,
        f"The longest simple repeat is {worst} tandem copies, within the limit of {permitted}.",
        measurements,
    )
