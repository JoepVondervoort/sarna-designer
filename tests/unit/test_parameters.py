"""Rule parameters carry their default, range and unit.

ADR 0003 requires these alongside the rule itself. They are what presets differ by: the
`wang-2015` GC band is 40-60 %, the `legacy-scanner` one is 40-65 %, and both must be visible
as a value with bounds rather than buried in a module constant the way the original script
buried them.
"""

import pytest
from pydantic import ValidationError

from sarna.rules.parameters import RuleParameter


def gc_minimum(**overrides: object) -> RuleParameter:
    defaults: dict[str, object] = {
        "name": "minimum-gc-fraction",
        "default": 0.40,
        "minimum": 0.0,
        "maximum": 1.0,
        "unit": "fraction",
        "description": "Lowest GC fraction a target may have.",
    }
    return RuleParameter(**{**defaults, **overrides})  # type: ignore[arg-type]


class TestConstruction:
    def test_builds_a_numeric_parameter(self):
        p = gc_minimum()
        assert p.default == 0.40
        assert p.unit == "fraction"

    def test_names_must_be_kebab_case(self):
        with pytest.raises(ValidationError, match="kebab-case"):
            gc_minimum(name="minimum_GC")

    def test_description_is_required(self):
        with pytest.raises(ValidationError, match="description"):
            gc_minimum(description="  ")

    def test_a_boolean_parameter_needs_no_bounds_or_unit(self):
        p = RuleParameter(
            name="enabled",
            default=True,
            description="Whether the rule applies.",
        )
        assert p.default is True
        assert p.minimum is None
        assert p.unit is None


class TestBounds:
    def test_default_must_lie_within_the_bounds(self):
        with pytest.raises(ValidationError, match="default"):
            gc_minimum(default=1.5)

    def test_minimum_must_not_exceed_maximum(self):
        with pytest.raises(ValidationError, match="minimum"):
            gc_minimum(minimum=0.9, maximum=0.1)

    def test_a_boolean_may_not_carry_numeric_bounds(self):
        with pytest.raises(ValidationError, match="boolean"):
            RuleParameter(
                name="enabled",
                default=True,
                minimum=0,
                maximum=1,
                description="Whether the rule applies.",
            )

    def test_bounds_are_optional(self):
        p = RuleParameter(
            name="window-size",
            default=4,
            unit="nt",
            description="How many flanking bases to inspect.",
        )
        assert p.minimum is None and p.maximum is None


class TestValidatingAnOverride:
    def test_accepts_a_value_inside_the_bounds(self):
        assert gc_minimum().validated(0.65) == 0.65

    def test_rejects_a_value_below_the_minimum(self):
        with pytest.raises(ValueError, match="minimum-gc-fraction"):
            gc_minimum().validated(-0.1)

    def test_rejects_a_value_above_the_maximum(self):
        with pytest.raises(ValueError, match="minimum-gc-fraction"):
            gc_minimum().validated(1.1)

    def test_rejects_a_value_of_the_wrong_type(self):
        with pytest.raises(ValueError, match="numeric"):
            gc_minimum().validated(True)

    def test_a_boolean_parameter_rejects_a_number(self):
        boolean = RuleParameter(
            name="enabled", default=True, description="Whether the rule applies."
        )
        with pytest.raises(ValueError, match="boolean"):
            boolean.validated(1)

    def test_rejects_a_non_numeric_value_entirely(self):
        with pytest.raises(ValueError, match="numeric"):
            gc_minimum().validated("0.5")  # type: ignore[arg-type]

    def test_a_boolean_parameter_accepts_a_boolean(self):
        boolean = RuleParameter(
            name="enabled", default=True, description="Whether the rule applies."
        )
        assert boolean.validated(False) is False

    def test_an_unbounded_parameter_accepts_any_number(self):
        unbounded = RuleParameter(
            name="window-size", default=4, unit="nt", description="Flank width."
        )
        assert unbounded.validated(99) == 99


class TestParameterIsData:
    def test_is_immutable(self):
        with pytest.raises(ValidationError):
            gc_minimum().default = 0.5

    def test_round_trips_through_serialisation(self):
        p = gc_minimum()
        assert RuleParameter.model_validate(p.model_dump()) == p
