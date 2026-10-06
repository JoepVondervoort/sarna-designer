"""The design-rule record.

A rule is data, not a function. Separating *what a rule is* from *how it is evaluated* is what
lets the UI render provenance next to every setting, lets exports stay interpretable away from
this codebase (ADR 0003), and lets a ranking be reported as an itemised checklist rather than
one opaque number (ADR 0006).
"""

import re
from enum import Enum
from typing import Annotated, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from sarna.rules.provenance import EvidenceClass, Provenance

_KEBAB_CASE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")

NonBlank = Annotated[str, Field(min_length=1)]


class RuleKind(Enum):
    """What a rule does to a candidate."""

    MANDATORY = "mandatory"
    """An admission requirement. Failing it excludes the candidate, naming this rule."""

    PREFERENCE = "preference"
    """Worth exactly one point if met. All preferences are weighted equally (ADR 0006)."""

    TIEBREAK = "tiebreak"
    """Orders candidates that are otherwise tied. Never folded into the score."""

    def __str__(self) -> str:
        """Return the kind label."""
        return self.value


class Rule(BaseModel):
    """A single design rule, with everything a reader needs to judge it.

    The explanatory fields are required rather than optional because a rule nobody can evaluate
    the merit of is exactly what this project exists to stop shipping.
    """

    model_config = ConfigDict(frozen=True)

    id: str
    """Stable kebab-case identifier. Referenced by presets and appears in exports."""

    name: NonBlank
    """Short human-readable label, as shown in the UI."""

    kind: RuleKind
    evidence_class: EvidenceClass
    provenance: Provenance

    rationale: NonBlank
    """Why this rule is believed to matter, biologically."""

    consequence: NonBlank
    """What changing it does to the candidate set, in plain language."""

    enabled_by_default: bool = True
    """Whether a preset enables this rule unless it says otherwise."""

    group: str | None = None
    """Optional UI grouping, so related single-point preferences can collapse into one line."""

    @model_validator(mode="after")
    def _check_id_and_integrity(self) -> Self:
        """Enforce the identifier format and the unsourced/established contradiction."""
        if not _KEBAB_CASE.match(self.id):
            raise ValueError(
                f"rule id {self.id!r} must be kebab-case, e.g. 'gc-content' or 'flank-20-au'"
            )
        if self.provenance.is_warning and self.evidence_class is EvidenceClass.ESTABLISHED:
            raise ValueError(
                f"rule {self.id!r} claims evidence class 'established' but its provenance is "
                f"unsourced; an unsourced rule cannot be established biology"
            )
        return self

    @model_validator(mode="after")
    def _check_explanatory_text_is_not_whitespace(self) -> Self:
        """Reject explanatory fields that are present but blank."""
        for field in ("name", "rationale", "consequence"):
            if not str(getattr(self, field)).strip():
                raise ValueError(f"{field} must not be blank for rule {self.id!r}")
        return self

    @property
    def points(self) -> int:
        """Return the points this rule contributes when met.

        One for a preference, zero for anything else. There is deliberately no other value:
        unequal weights are what made the legacy score uninterpretable (ADR 0006).
        """
        return 1 if self.kind is RuleKind.PREFERENCE else 0

    @property
    def is_gate(self) -> bool:
        """Return whether failing this rule excludes the candidate outright."""
        return self.kind is RuleKind.MANDATORY

    @property
    def needs_warning(self) -> bool:
        """Return whether the UI must show a provenance warning beside this rule."""
        return self.provenance.is_warning
