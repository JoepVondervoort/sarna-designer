"""A rule paired with the parameter values a particular run will use.

This is the second argument of the evaluator signature fixed by ADR 0007. Its job is to make an
evaluator's life boring: by the time one runs, every parameter has a value, every value has been
checked against the parameter's own type and bounds, and asking for a value by the wrong type
raises instead of silently coercing.

The alternative (handing evaluators the raw override mapping) would scatter "what if this is
missing" and "what if this is out of range" across ten rules, which is where a preset typo turns
into a default quietly left in place.

Resolution is per rule rather than global. A preset is a selection of rules each with its own
overrides, so the unit that gets resolved is one rule's parameters, and the resolved values stay
inspectable as data for export (ADR 0003).
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

from sarna.rules.parameters import ParameterValue
from sarna.rules.rule import Rule, RuleKind


@dataclass(frozen=True, slots=True)
class ResolvedRule:
    """A rule and the validated, fully-defaulted parameter values to evaluate it with.

    Args:
        rule: The rule being evaluated.
        overrides: Values replacing the parameter defaults, typically from a preset or the UI.
            Copied, so a later change to the caller's mapping cannot alter this object.
        kind_override: Mandatory or preference for this run, leaving the catalogue unchanged.
        option_id: A deliberately selected sourced option. A matching number is not provenance.
    """

    rule: Rule
    overrides: Mapping[str, ParameterValue] = field(default_factory=dict)
    kind_override: RuleKind | None = None
    option_id: str | None = None
    values: Mapping[str, ParameterValue] = field(init=False)

    def __post_init__(self) -> None:
        """Validate the overrides and merge them over the parameter defaults.

        Raises:
            KeyError: If an override names a parameter this rule does not have, which is almost
                always a typo in a preset and must not be ignored.
            ValueError: If a value is the wrong type for its parameter or outside its bounds.
        """
        if self.kind_override is not None:
            if not isinstance(self.kind_override, RuleKind):
                raise ValueError("kind_override must be a RuleKind")
            if self.kind_override is not self.rule.kind:
                if RuleKind.TIEBREAK in (self.rule.kind, self.kind_override):
                    raise ValueError(
                        f"rule {self.rule.id!r}: a continuous tie-break cannot become a "
                        "gate or preference, or vice versa"
                    )
                if self.rule.id == "target-length-19":
                    raise ValueError(
                        "target-length-19 defines candidate generation and cannot become "
                        "a scored preference"
                    )
        resolved: dict[str, ParameterValue] = {
            parameter.name: parameter.default for parameter in self.rule.parameters
        }
        if self.option_id is not None:
            resolved.update(self.rule.option(self.option_id).values)
        for name, value in self.overrides.items():
            if name not in resolved:
                known = ", ".join(sorted(resolved)) or "none"
                raise KeyError(
                    f"rule {self.rule.id!r} has no parameter {name!r}; it has {known}"
                )
            resolved[name] = self.rule.parameter(name).validated(value)
        self._check_relationships(resolved)
        object.__setattr__(self, "overrides", MappingProxyType(dict(self.overrides)))
        object.__setattr__(self, "values", MappingProxyType(resolved))

    @property
    def id(self) -> str:
        """Return the unchanged catalogue identifier."""
        return self.rule.id

    @property
    def kind(self) -> RuleKind:
        """Return the effective role used by evaluation, gating, scoring and exports."""
        return self.kind_override if self.kind_override is not None else self.rule.kind

    def _check_relationships(self, values: Mapping[str, ParameterValue]) -> None:
        """Reject contradictory settings before they reach an evaluator."""
        ordered_pairs = (
            ("minimum-gc-fraction", "maximum-gc-fraction"),
            ("window-start-bp", "window-end-bp"),
            ("preferred-window-start-bp", "preferred-window-end-bp"),
            ("central-region-start", "central-region-end"),
        )
        for first, last in ordered_pairs:
            left, right = values.get(first), values.get(last)
            if isinstance(left, int | float) and isinstance(right, int | float) and left > right:
                raise ValueError(f"rule {self.id!r}: {first} must not exceed {last}")

    def to_dict(self) -> dict[str, Any]:
        """Export effective behavior alongside the catalogue role and choice origin."""
        exported = self.rule.to_dict(self.values)
        exported["kind"] = str(self.kind)
        exported["default_kind"] = str(self.rule.kind)
        exported["kind_overridden"] = self.kind is not self.rule.kind
        exported["application_note"] = (
            "Gate, preference and off describe software behavior, not evidence strength. "
            "Custom values and roles are researcher choices; the cited source does not "
            "automatically support those changes."
        )
        exported["selected_option"] = (
            self.rule.option(self.option_id).to_dict() if self.option_id is not None else None
        )
        for parameter in exported["parameters"]:
            name = parameter["name"]
            parameter["origin"] = (
                "override"
                if name in self.overrides
                else "named-option"
                if self.option_id is not None and name in self.rule.option(self.option_id).values
                else "catalogue-default"
            )
        return exported

    @classmethod
    def from_option(cls, rule: Rule, option_id: str) -> "ResolvedRule":
        """Resolve ``rule`` using one of its published options.

        Choosing a named option is how a documented alternative should be selected, so that the
        resulting configuration records *which* alternative was chosen and can show its
        provenance, rather than recording a bare number typed in from nowhere.

        Args:
            rule: The rule being evaluated.
            option_id: The option to apply, as published by the rule.

        Returns:
            The rule resolved with that option's values.

        Raises:
            KeyError: If the rule has no such option.
        """
        return cls(rule, option_id=option_id)

    def _value(self, name: str) -> ParameterValue:
        """Return a resolved value, or raise if the parameter does not exist."""
        try:
            return self.values[name]
        except KeyError:
            known = ", ".join(sorted(self.values)) or "none"
            raise KeyError(
                f"rule {self.rule.id!r} has no parameter {name!r}; it has {known}"
            ) from None

    def number(self, name: str) -> float:
        """Return a numeric parameter's value.

        Raises:
            KeyError: If the parameter does not exist.
            TypeError: If the parameter is textual. Coercing a table name to a number would
                produce a nonsensical comparison rather than a traceback.
        """
        value = self._value(name)
        if isinstance(value, str):
            raise TypeError(
                f"parameter {name!r} of rule {self.rule.id!r} is a choice, not a number"
            )
        return float(value)

    def integer(self, name: str) -> int:
        """Return a whole-number parameter's value.

        Raises:
            KeyError: If the parameter does not exist.
            TypeError: If the parameter is textual, or numeric but not a whole number.
        """
        value = self.number(name)
        if value != int(value):
            raise TypeError(
                f"parameter {name!r} of rule {self.rule.id!r} is {value}, which is not a whole "
                f"number"
            )
        return int(value)

    def text(self, name: str) -> str:
        """Return a choice parameter's value.

        Raises:
            KeyError: If the parameter does not exist.
            TypeError: If the parameter is numeric or boolean.
        """
        value = self._value(name)
        if not isinstance(value, str):
            raise TypeError(
                f"parameter {name!r} of rule {self.rule.id!r} is numeric, not a choice"
            )
        return value
