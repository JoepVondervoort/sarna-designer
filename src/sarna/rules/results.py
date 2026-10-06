"""What a rule reports once it has been evaluated against a candidate.

A rule is data (:mod:`sarna.rules.rule`); evaluating one produces a :class:`RuleResult`. Every
claim the UI, the ranking and the exports make about a candidate is one of these, which is why
the type is strict about three things in particular.

**A verdict always states its reason.** ``detail`` is required and may not be blank, so no rule
can report a failure the researcher cannot interpret (ADR 0003). The measured values travel
alongside the prose in ``measurements``, so a table column or a plot never has to parse English.

**"Could not be evaluated" is its own outcome.** ADR 0007 gives coordinate-dependent rules a
third answer for the Workflow B case, where a promoter is pasted with no genomic coordinates and
a TSS-relative rule has nothing to measure. :attr:`RuleResult.excludes` is deliberately false
for a not-applicable mandatory rule: treating it as a failure would reject every candidate in a
sequence-only run, and treating it as a pass would claim a check that never happened.

**A tie-break never scores.** Only a met preference is worth a point, and only a tie-break may
carry an ordering value (ADR 0006). The legacy ``0.01·Δstab + 1000·pos18 + …`` score cannot be
expressed with this type at all, which is the point.
"""

from collections.abc import Mapping
from enum import Enum
from typing import Self

from pydantic import BaseModel, ConfigDict, model_validator

from sarna.rules.parameters import ParameterValue
from sarna.rules.rule import Rule, RuleKind


class RuleOutcome(Enum):
    """The three answers a rule can give about a candidate."""

    PASS = "pass"
    """The candidate satisfies the rule. For a preference, this is what earns the point."""

    FAIL = "fail"
    """The candidate does not satisfy the rule. For a mandatory rule, this excludes it."""

    NOT_APPLICABLE = "not-applicable"
    """The rule could not be evaluated, and says why. Neither a pass nor a failure."""

    def __str__(self) -> str:
        """Return the outcome label as it appears in exports."""
        return self.value


class RuleResult(BaseModel):
    """One rule's verdict on one candidate, with the reason and the numbers behind it.

    Prefer the :meth:`passed`, :meth:`failed`, :meth:`not_applicable` and :meth:`tiebreak`
    constructors over building one directly: they take the :class:`~sarna.rules.rule.Rule` and
    copy its kind across, so an evaluator cannot accidentally report a preference as a gate.
    """

    model_config = ConfigDict(frozen=True)

    rule_id: str
    """The ``id`` of the rule that produced this result. Appears in exports."""

    kind: RuleKind
    """Copied from the rule, so an exported row says what the verdict *means* on its own."""

    outcome: RuleOutcome

    detail: str
    """Why, in plain language, including the measured value. Required, never blank."""

    measurements: Mapping[str, ParameterValue] = {}
    """The values behind the verdict, keyed by name, for tables, plots and exports."""

    ordering_value: float | None = None
    """A tie-break's magnitude. Only a tie-break may set it; it never affects the score."""

    @classmethod
    def passed(
        cls,
        rule: Rule,
        detail: str,
        measurements: Mapping[str, ParameterValue] | None = None,
    ) -> "RuleResult":
        """Report that ``rule`` is satisfied.

        Args:
            rule: The rule evaluated.
            detail: What was measured and why it satisfies the rule.
            measurements: The values behind the verdict.
        """
        return cls(
            rule_id=rule.id,
            kind=rule.kind,
            outcome=RuleOutcome.PASS,
            detail=detail,
            measurements=measurements or {},
        )

    @classmethod
    def failed(
        cls,
        rule: Rule,
        detail: str,
        measurements: Mapping[str, ParameterValue] | None = None,
    ) -> "RuleResult":
        """Report that ``rule`` is not satisfied.

        Args:
            rule: The rule evaluated.
            detail: What was measured and how it misses the rule. For a mandatory rule this is
                the text shown to explain why the candidate was excluded, so it names both the
                measured value and the limit it crossed.
            measurements: The values behind the verdict.
        """
        return cls(
            rule_id=rule.id,
            kind=rule.kind,
            outcome=RuleOutcome.FAIL,
            detail=detail,
            measurements=measurements or {},
        )

    @classmethod
    def not_applicable(cls, rule: Rule, reason: str) -> "RuleResult":
        """Report that ``rule`` could not be evaluated, and why.

        Used when the inputs the rule needs are absent rather than unsatisfactory — most often a
        TSS-relative rule on a pasted promoter with no coordinates (ADR 0007).

        Args:
            rule: The rule that could not be evaluated.
            reason: What was missing. Shown to the user as-is, so it names the missing input
                rather than the internal field that was ``None``.
        """
        return cls(
            rule_id=rule.id,
            kind=rule.kind,
            outcome=RuleOutcome.NOT_APPLICABLE,
            detail=reason,
        )

    @classmethod
    def tiebreak(
        cls,
        rule: Rule,
        value: float,
        detail: str,
        measurements: Mapping[str, ParameterValue] | None = None,
    ) -> "RuleResult":
        """Report a tie-break's ordering value.

        Args:
            rule: The tie-break rule.
            value: The magnitude candidates are ordered by, larger first.
            detail: The value in plain language, with its unit.
            measurements: The values behind it.
        """
        return cls(
            rule_id=rule.id,
            kind=rule.kind,
            outcome=RuleOutcome.PASS,
            detail=detail,
            measurements=measurements or {},
            ordering_value=value,
        )

    @property
    def is_unevaluated(self) -> bool:
        """Return whether the rule could not be evaluated at all."""
        return self.outcome is RuleOutcome.NOT_APPLICABLE

    @property
    def excludes(self) -> bool:
        """Return whether this result excludes the candidate from the results.

        Only a mandatory rule that was actually evaluated and actually failed. A mandatory rule
        reporting ``not-applicable`` does **not** exclude: see the module docstring and ADR 0007.
        """
        return self.kind is RuleKind.MANDATORY and self.outcome is RuleOutcome.FAIL

    @property
    def points(self) -> int:
        """Return the points this result contributes to the rank: one, or none.

        One only for a preference that was met. Mandatory rules and tie-breaks score nothing,
        and an unevaluated preference scores nothing — it should also be left out of the
        denominator when a rank is presented as "meets 5 of 7", since a preference nobody could
        check is not one the candidate failed to meet.
        """
        if self.kind is RuleKind.PREFERENCE and self.outcome is RuleOutcome.PASS:
            return 1
        return 0

    @model_validator(mode="after")
    def _check_detail_is_present(self) -> Self:
        """Require a non-blank reason for every verdict."""
        if not self.detail.strip():
            raise ValueError(
                f"rule {self.rule_id!r} reported {self.outcome} with no stated reason; every "
                f"verdict must say what was measured, or why it could not be"
            )
        return self

    @model_validator(mode="after")
    def _check_ordering_value_belongs_to_a_tiebreak(self) -> Self:
        """Confine ordering values to tie-breaks, and require one unless unevaluated."""
        is_tiebreak = self.kind is RuleKind.TIEBREAK
        if self.ordering_value is not None and not is_tiebreak:
            raise ValueError(
                f"rule {self.rule_id!r} is a {self.kind} rule but carries an ordering value; "
                f"only a tie-break may, because anything else is a weighted score (ADR 0006)"
            )
        if is_tiebreak and self.ordering_value is None and not self.is_unevaluated:
            raise ValueError(
                f"tie-break {self.rule_id!r} reported {self.outcome} but no ordering value; a "
                f"tie-break exists only to supply one"
            )
        return self
