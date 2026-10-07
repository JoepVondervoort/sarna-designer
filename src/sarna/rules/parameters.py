"""Tunable values attached to a design rule.

ADR 0003 requires a rule to publish its ``default``, ``valid_range`` and ``unit`` where they
are meaningful. Those three together are what let the UI render a setting the user can reason
about, and what makes the difference between presets inspectable: ``wang-2015`` and
``legacy-scanner`` disagree about the upper GC bound, and that disagreement should read as one
parameter with two values, not as two different programs.
"""

import re
from math import isfinite
from typing import Any, Self

from pydantic import BaseModel, ConfigDict, model_validator

_KEBAB_CASE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")

ParameterValue = bool | int | float | str


class RuleParameter(BaseModel):
    """A single tunable value belonging to a rule."""

    model_config = ConfigDict(frozen=True)

    name: str
    """Stable kebab-case identifier, unique within its rule."""

    default: ParameterValue
    """The value used unless a preset overrides it."""

    minimum: float | None = None
    maximum: float | None = None

    unit: str | None = None
    """Unit of measurement, e.g. ``nt``, ``bp``, ``fraction``. ``None`` for booleans."""

    choices: tuple[str, ...] | None = None
    """For a parameter chosen from a named set, the legal values. ``None`` for numerics."""

    description: str
    """What this value controls, in plain language."""

    @property
    def is_boolean(self) -> bool:
        """Return whether this parameter is on/off rather than numeric."""
        return isinstance(self.default, bool)

    @property
    def is_choice(self) -> bool:
        """Return whether this parameter is picked from a named set rather than numeric."""
        return isinstance(self.default, str)

    @model_validator(mode="after")
    def _check_name_bounds_and_default(self) -> Self:
        """Validate the name format, the bounds, and that the default satisfies them."""
        if not _KEBAB_CASE.match(self.name):
            raise ValueError(
                f"parameter name {self.name!r} must be kebab-case, e.g. 'minimum-gc-fraction'"
            )
        if not self.description.strip():
            raise ValueError(f"description must not be blank for parameter {self.name!r}")
        if self.is_boolean:
            if self.minimum is not None or self.maximum is not None:
                raise ValueError(
                    f"parameter {self.name!r} is boolean and cannot carry numeric bounds"
                )
            return self
        if self.is_choice:
            if not self.choices:
                raise ValueError(
                    f"parameter {self.name!r} is chosen from a set and must list its choices"
                )
            if self.default not in self.choices:
                raise ValueError(
                    f"default {self.default!r} is not among the choices for "
                    f"parameter {self.name!r}"
                )
            return self
        if (
            self.minimum is not None
            and self.maximum is not None
            and self.minimum > self.maximum
        ):
            raise ValueError(
                f"minimum {self.minimum} exceeds maximum {self.maximum} for "
                f"parameter {self.name!r}"
            )
        assert not isinstance(self.default, str)
        self._check_in_bounds(self.default, what="default")
        return self

    def _check_in_bounds(self, value: float, what: str) -> None:
        if not isfinite(value):
            raise ValueError(f"{what} must be finite for parameter {self.name!r}")
        if self.minimum is not None and value < self.minimum:
            raise ValueError(
                f"{what} {value} is below the minimum {self.minimum} for "
                f"parameter {self.name!r}"
            )
        if self.maximum is not None and value > self.maximum:
            raise ValueError(
                f"{what} {value} is above the maximum {self.maximum} for "
                f"parameter {self.name!r}"
            )

    def validated(self, value: ParameterValue) -> ParameterValue:
        """Return ``value`` if it is a legal override for this parameter.

        Args:
            value: The candidate override, typically from a preset or the UI.

        Returns:
            The value unchanged.

        Raises:
            ValueError: If the value is the wrong type or outside the bounds.
        """
        if self.is_boolean:
            if not isinstance(value, bool):
                raise ValueError(
                    f"parameter {self.name!r} is boolean; got {value!r}"
                )
            return value
        if self.is_choice:
            if value not in (self.choices or ()):
                raise ValueError(
                    f"parameter {self.name!r} must be one of "
                    f"{', '.join(self.choices or ())}; got {value!r}"
                )
            return value
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise ValueError(f"parameter {self.name!r} is numeric; got {value!r}")
        self._check_in_bounds(value, what="value")
        if isinstance(self.default, int) and isinstance(value, float) and not value.is_integer():
            raise ValueError(f"parameter {self.name!r} requires a whole number; got {value!r}")
        return value

    def to_dict(self, value: ParameterValue | None = None) -> dict[str, Any]:
        """Return the parameter as JSON-serialisable data.

        Args:
            value: The value a run resolved this parameter to. Supplied, the export carries the
                declaration and the value used in one object, because a reader needs both at once:
                ``0.6`` means nothing without "of a possible 0.0 to 1.0, as a fraction". Omitted,
                only the declaration is exported, with no ``value`` key.

        Returns:
            The parameter. With a value, ``is_default`` records whether the run left the default
            alone, which is how a preset's divergence becomes visible without a second lookup.
        """
        exported: dict[str, Any] = {
            "name": self.name,
            "default": self.default,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "unit": self.unit,
            "choices": list(self.choices) if self.choices is not None else None,
            "description": self.description,
        }
        if value is not None:
            exported["value"] = value
            exported["is_default"] = value == self.default
        return exported
