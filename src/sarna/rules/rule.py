"""The design-rule record.

A rule is data, not a function. Separating *what a rule is* from *how it is evaluated* is what
lets the UI render provenance next to every setting, lets exports stay interpretable away from
this codebase (ADR 0003), and lets a ranking be reported as an itemised checklist rather than
one opaque number (ADR 0006).
"""

import re
from collections.abc import Mapping
from enum import Enum
from typing import Annotated, Any, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from sarna.rules.options import RuleOption
from sarna.rules.parameters import ParameterValue, RuleParameter
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

    parameters: tuple[RuleParameter, ...] = ()
    """Tunable values this rule exposes, with their defaults, bounds and units (ADR 0003)."""

    options: tuple[RuleOption, ...] = ()
    """Named alternative settings the user may choose, each with its own provenance."""

    @property
    def default_option(self) -> RuleOption | None:
        """Return the option marked default, or ``None`` if this rule offers no options."""
        for option in self.options:
            if option.is_default:
                return option
        return None

    def option(self, option_id: str) -> RuleOption:
        """Return the option called ``option_id``.

        Raises:
            KeyError: If this rule has no such option.
        """
        for option in self.options:
            if option.id == option_id:
                return option
        raise KeyError(f"rule {self.id!r} has no option {option_id!r}")

    def parameter(self, name: str) -> RuleParameter:
        """Return the parameter called ``name``.

        Raises:
            KeyError: If this rule has no such parameter.
        """
        for parameter in self.parameters:
            if parameter.name == name:
                return parameter
        raise KeyError(f"rule {self.id!r} has no parameter {name!r}")

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

    @model_validator(mode="after")
    def _check_options_are_coherent(self) -> Self:
        """Require options to be unique, defaulted once, and to set real, valid parameters."""
        if not self.options:
            return self
        ids = [option.id for option in self.options]
        if len(set(ids)) != len(ids):
            raise ValueError(f"rule {self.id!r} has duplicate option ids")
        defaults = [option for option in self.options if option.is_default]
        if len(defaults) != 1:
            raise ValueError(
                f"rule {self.id!r} offers options, so exactly one must be marked default; "
                f"found {len(defaults)}"
            )
        by_name = {parameter.name: parameter for parameter in self.parameters}
        for option in self.options:
            for name, value in option.values.items():
                if name not in by_name:
                    raise ValueError(
                        f"option {option.id!r} of rule {self.id!r} sets unknown "
                        f"parameter {name!r}"
                    )
                by_name[name].validated(value)
        for name, value in defaults[0].values.items():
            if by_name[name].default != value:
                raise ValueError(
                    f"rule {self.id!r}: the default option {defaults[0].id!r} sets "
                    f"{name}={value!r} but the parameter default is "
                    f"{by_name[name].default!r}; these must agree"
                )
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

    @property
    def evidence_note(self) -> str:
        """State the scope of support without implying validated activation prediction."""
        source = self.provenance.source or ""
        if self.provenance.is_warning:
            return "No traceable supporting source. This is an unsupported design hypothesis."
        if self.provenance.kind.value == "legacy-script":
            return "Inherited implementation choice; legacy behavior is not biological evidence."
        if source in ("Reynolds et al. 2004", "Ui-Tei et al. 2004"):
            return (
                "Evidence from siRNA knockdown experiments. Transfer to promoter-targeted "
                "saRNA activation is an extrapolation, not demonstrated efficacy."
            )
        if source == "van de Vondervoort 2026":
            return (
                "Project thesis evidence or hypothesis, not independent validation. "
                "Rules derived from these data cannot be validated on the same data."
            )
        if self.evidence_class is EvidenceClass.EXPLORATORY:
            return "Exploratory hypothesis; improved saRNA activation has not been demonstrated."
        if self.kind is RuleKind.TIEBREAK:
            return (
                "RNAi strand-selection research motivates this calculation. Using terminal "
                "energy as a ranking tiebreak is a project choice without demonstrated "
                "predictive benefit for saRNA activation."
            )
        if self.id == "not-in-cpg-island":
            return (
                "Wang recommends avoiding CpG islands. The detection definition and "
                "sliding-window implementation are project choices; island sequence does "
                "not establish methylation, accessibility or saRNA response."
            )
        if source == "Wang et al. 2015":
            return (
                "Published saRNA protocol guidance. A citation supports the recommendation, "
                "not a universal biological requirement or this application's ranking."
            )
        return (
            "Published support applies to the phenomenon studied in the cited work. "
            "It does not establish that this rule predicts saRNA activation."
        )

    def to_dict(
        self, values: Mapping[str, ParameterValue] | None = None
    ) -> dict[str, Any]:
        """Return the rule record as JSON-serialisable data.

        ADR 0003 commits to exports that carry the rule records rather than only the values a run
        resolved. ``rationale`` and ``consequence`` are the largest body of "why this choice" text
        in the project, and without them an export states a threshold while withholding the reason
        for it.

        Args:
            values: The parameter values a run resolved, keyed by parameter name. Supplied, each
                parameter is exported with its value and whether that value is the default, so a
                preset's divergence is visible without a second lookup.

        Returns:
            The record, with ``parameters`` and ``options`` as lists of objects. A list rather
            than a mapping because each entry is self-describing, and a renderer that wants a
            lookup can build one in a line.
        """
        return {
            "id": self.id,
            "name": self.name,
            "kind": str(self.kind),
            "evidence_class": str(self.evidence_class),
            "evidence_note": self.evidence_note,
            "provenance": self.provenance.to_dict(),
            "rationale": self.rationale,
            "consequence": self.consequence,
            "enabled_by_default": self.enabled_by_default,
            "group": self.group,
            "parameters": [
                parameter.to_dict(None if values is None else values.get(parameter.name))
                for parameter in self.parameters
            ],
            "options": [option.to_dict() for option in self.options],
        }
