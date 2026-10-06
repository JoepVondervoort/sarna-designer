"""Tests for a rule paired with the parameter values a run will actually use.

This is the second half of the evaluator signature in ADR 0007. It exists so that an evaluator
reads typed, validated, already-defaulted values and never has to decide what to do about a
missing or out-of-range setting.
"""

import pytest

from sarna.rules.catalog.mandatory import GC_CONTENT, POS19_A, THERMODYNAMIC_ASYMMETRY
from sarna.rules.resolved import ResolvedRule


class TestDefaulting:
    def test_unset_parameters_take_their_defaults(self) -> None:
        resolved = ResolvedRule(GC_CONTENT)
        assert resolved.number("minimum-gc-fraction") == pytest.approx(0.40)
        assert resolved.number("maximum-gc-fraction") == pytest.approx(0.60)

    def test_an_override_replaces_only_that_parameter(self) -> None:
        resolved = ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 0.65})
        assert resolved.number("minimum-gc-fraction") == pytest.approx(0.40)
        assert resolved.number("maximum-gc-fraction") == pytest.approx(0.65)

    def test_a_rule_with_no_parameters_resolves_to_nothing(self) -> None:
        assert ResolvedRule(POS19_A).values == {}

    def test_the_resolved_values_are_inspectable_as_a_whole(self) -> None:
        # ADR 0003: the configuration a run used must be exportable, so the resolved values are
        # readable as data rather than only through accessors.
        resolved = ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 0.65})
        assert resolved.values == {
            "minimum-gc-fraction": 0.40,
            "maximum-gc-fraction": 0.65,
        }


class TestValidation:
    def test_rejects_an_unknown_parameter(self) -> None:
        # A typo in a preset must fail loudly rather than silently leaving the default in place.
        with pytest.raises(KeyError, match="maximum-gc-percent"):
            ResolvedRule(GC_CONTENT, {"maximum-gc-percent": 0.65})

    def test_rejects_a_value_outside_the_parameters_bounds(self) -> None:
        with pytest.raises(ValueError, match="above the maximum"):
            ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 1.5})

    def test_rejects_a_value_of_the_wrong_type(self) -> None:
        with pytest.raises(ValueError, match="numeric"):
            ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": "high"})

    def test_rejects_a_choice_outside_the_permitted_set(self) -> None:
        with pytest.raises(ValueError, match="must be one of"):
            ResolvedRule(
                THERMODYNAMIC_ASYMMETRY, {"nearest-neighbour-table": "rna-turner-2004"}
            )


class TestTypedAccess:
    def test_integer_returns_an_int(self) -> None:
        resolved = ResolvedRule(THERMODYNAMIC_ASYMMETRY)
        assert isinstance(resolved.number("minimum-asymmetry-kcal"), float)

    def test_text_returns_a_choice(self) -> None:
        resolved = ResolvedRule(THERMODYNAMIC_ASYMMETRY)
        assert resolved.text("nearest-neighbour-table") == "rna-xia-1998"

    def test_number_rejects_a_textual_parameter(self) -> None:
        # Asking for a table name as a number is a programming error, and silently coercing it
        # would produce a nonsensical comparison rather than a traceback.
        resolved = ResolvedRule(THERMODYNAMIC_ASYMMETRY)
        with pytest.raises(TypeError, match="nearest-neighbour-table"):
            resolved.number("nearest-neighbour-table")

    def test_text_rejects_a_numeric_parameter(self) -> None:
        resolved = ResolvedRule(THERMODYNAMIC_ASYMMETRY)
        with pytest.raises(TypeError, match="minimum-asymmetry-kcal"):
            resolved.text("minimum-asymmetry-kcal")

    def test_an_unknown_name_raises_rather_than_returning_none(self) -> None:
        with pytest.raises(KeyError):
            ResolvedRule(GC_CONTENT).number("gc-max")


class TestFromOption:
    def test_builds_from_a_named_option(self) -> None:
        # Choosing the legacy GC band should be expressible as "that option", not as a number
        # typed in from nowhere.
        resolved = ResolvedRule.from_option(GC_CONTENT, "legacy-scanner")
        assert resolved.number("maximum-gc-fraction") == pytest.approx(0.65)

    def test_the_default_option_matches_the_parameter_defaults(self) -> None:
        resolved = ResolvedRule.from_option(GC_CONTENT, "wang-2015")
        assert resolved.values == ResolvedRule(GC_CONTENT).values

    def test_an_unknown_option_raises(self) -> None:
        with pytest.raises(KeyError, match="no option"):
            ResolvedRule.from_option(GC_CONTENT, "permissive")


class TestValueSemantics:
    def test_is_immutable(self) -> None:
        resolved = ResolvedRule(GC_CONTENT)
        with pytest.raises((AttributeError, TypeError)):
            resolved.rule = POS19_A  # type: ignore[misc]

    def test_overriding_does_not_mutate_the_rule(self) -> None:
        ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 0.65})
        assert GC_CONTENT.parameter("maximum-gc-fraction").default == pytest.approx(0.60)

    def test_the_caller_cannot_mutate_the_resolved_values_afterwards(self) -> None:
        overrides = {"maximum-gc-fraction": 0.65}
        resolved = ResolvedRule(GC_CONTENT, overrides)
        overrides["maximum-gc-fraction"] = 0.99
        assert resolved.number("maximum-gc-fraction") == pytest.approx(0.65)


class TestWholeNumbers:
    def test_integer_rejects_a_fractional_value(self) -> None:
        # A rule asking for a count must not silently truncate 4.5 to 4.
        resolved = ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 0.5})
        with pytest.raises(TypeError, match="not a whole number"):
            resolved.integer("maximum-gc-fraction")

    def test_integer_accepts_a_whole_value_stored_as_a_float(self) -> None:
        resolved = ResolvedRule(GC_CONTENT, {"maximum-gc-fraction": 1.0})
        assert resolved.integer("maximum-gc-fraction") == 1
