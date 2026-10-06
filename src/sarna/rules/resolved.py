"""A rule paired with the parameter values a particular run will use.

This is the second argument of the evaluator signature fixed by ADR 0007. Its job is to make an
evaluator's life boring: by the time one runs, every parameter has a value, every value has been
checked against the parameter's own type and bounds, and asking for a value by the wrong type
raises instead of silently coercing.

The alternative — handing evaluators the raw override mapping — would scatter "what if this is
missing" and "what if this is out of range" across ten rules, which is where a preset typo turns
into a default quietly left in place.

Resolution is per rule rather than global. A preset is a selection of rules each with its own
overrides, so the unit that gets resolved is one rule's parameters, and the resolved values stay
inspectable as data for export (ADR 0003).
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from sarna.rules.parameters import ParameterValue
from sarna.rules.rule import Rule


@dataclass(frozen=True, slots=True)
class ResolvedRule:
    """A rule and the validated, fully-defaulted parameter values to evaluate it with.

    Args:
        rule: The rule being evaluated.
        overrides: Values replacing the parameter defaults, typically from a preset or the UI.
            Copied, so a later change to the caller's mapping cannot alter this object.
    """

    rule: Rule
    overrides: Mapping[str, ParameterValue] = field(default_factory=dict)
    values: Mapping[str, ParameterValue] = field(init=False)

    def __post_init__(self) -> None:
        """Validate the overrides and merge them over the parameter defaults.

        Raises:
            KeyError: If an override names a parameter this rule does not have, which is almost
                always a typo in a preset and must not be ignored.
            ValueError: If a value is the wrong type for its parameter or outside its bounds.
        """
        resolved: dict[str, ParameterValue] = {
            parameter.name: parameter.default for parameter in self.rule.parameters
        }
        for name, value in self.overrides.items():
            if name not in resolved:
                known = ", ".join(sorted(resolved)) or "none"
                raise KeyError(
                    f"rule {self.rule.id!r} has no parameter {name!r}; it has {known}"
                )
            resolved[name] = self.rule.parameter(name).validated(value)
        object.__setattr__(self, "overrides", MappingProxyType(dict(self.overrides)))
        object.__setattr__(self, "values", MappingProxyType(resolved))

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
        return cls(rule, rule.option(option_id).values)

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
