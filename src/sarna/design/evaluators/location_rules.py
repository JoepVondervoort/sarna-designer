"""Evaluators for rules about where the target sits, rather than what it says.

These are the rules a sequence-only run cannot answer. Each one checks for the context it needs
and reports ``not-applicable`` when it is absent, which under ADR 0007 neither excludes the
candidate nor counts as a pass. Without that, Workflow B, a pasted promoter with no coordinates
would either lose every candidate to the TSS window or quietly claim checks that never ran.
"""

from sarna.design.context import CandidateContext
from sarna.genome.tss import TssPosition
from sarna.rules.resolved import ResolvedRule
from sarna.rules.results import RuleResult

_NO_TSS = (
    "This promoter has no known transcription start site, so the target's distance from the "
    "TSS cannot be measured. Supply a TSS, or design from a gene rather than a pasted sequence."
)


def _tss_window(
    context: CandidateContext,
    resolved: ResolvedRule,
    start_parameter: str,
    end_parameter: str,
) -> RuleResult:
    """Evaluate total containment within a TSS-relative window (ADR 0008)."""
    start_bp = resolved.integer(start_parameter)
    end_bp = resolved.integer(end_parameter)
    contained = context.tss_window_contains(start_bp, end_bp)
    if contained is None:
        return RuleResult.not_applicable(resolved, _NO_TSS)
    assert context.tss_span is not None
    first, last = context.tss_span
    window = f"{TssPosition(start_bp)} to {TssPosition(end_bp)}"
    measurements = {"tss-start": first.value, "tss-end": last.value}
    span = f"{first} to {last}"
    if contained:
        return RuleResult.passed(
            resolved,
            f"The target spans {span} relative to the TSS, inside the {window} window.",
            measurements,
        )
    return RuleResult.failed(
        resolved,
        f"The target spans {span} relative to the TSS, which is not wholly inside the "
        f"{window} window.",
        measurements,
    )


def tss_window(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the whole target lies inside the mandatory TSS-relative window."""
    return _tss_window(context, resolved, "window-start-bp", "window-end-bp")


def tss_preferred_window(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the whole target lies inside the preferred sub-window."""
    return _tss_window(
        context, resolved, "preferred-window-start-bp", "preferred-window-end-bp"
    )


def not_in_cpg_island(context: CandidateContext, resolved: ResolvedRule) -> RuleResult:
    """Check the target does not overlap a CpG island.

    The island criteria are the rule's own parameters, but detection runs once per promoter
    rather than per candidate, so the islands arrive in the context. ``None`` there means
    detection was not run, which is reported rather than treated as a clean result.
    """
    minimum = resolved.integer("island-minimum-length")
    overlaps = context.overlaps_cpg_island
    if context.promoter.length < minimum:
        return RuleResult.not_applicable(
            resolved,
            f"Only {context.promoter.length} bp of sequence are available; the selected CpG "
            f"island definition needs at least {minimum} bp. Island overlap is unknown.",
        )
    if overlaps is None:
        return RuleResult.not_applicable(
            resolved,
            "CpG island detection was not run for this promoter, so whether the target lies "
            "in an island is unknown.",
        )
    islands = context.cpg_islands or ()
    measurements = {"islands-in-promoter": len(islands)}
    if not overlaps:
        if (
            context.offset < minimum - 1
            or context.promoter.length - context.offset - len(context.target) < minimum - 1
        ):
            return RuleResult.not_applicable(
                resolved,
                "No overlapping island was detected in the supplied sequence, but the target "
                "is too close to a sequence boundary to exclude an island extending beyond it. "
                "Island overlap is unknown; supply more flanking sequence.",
            )
        found = (
            f"{len(islands)} island(s) were found in the promoter, none overlapping the target"
            if islands
            else "no CpG islands were found in the promoter"
        )
        return RuleResult.passed(
            resolved, f"The target is outside any CpG island: {found}.", measurements
        )
    touched = next(
        island
        for island in islands
        if island.overlaps(context.offset, context.offset + context.length)
    )
    return RuleResult.failed(
        resolved,
        f"The target overlaps a CpG island spanning promoter positions "
        f"{touched.start + 1} to {touched.end} "
        f"({touched.length} bp, {touched.gc_fraction:.0%} GC, observed/expected CpG "
        f"{touched.observed_expected_ratio:.2f}).",
        measurements,
    )
