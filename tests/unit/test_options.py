"""Optional alternatives, each carrying its own source.

Where the paper, the original script and the wider literature disagree, the alternatives are
offered as options rather than buried in preset definitions. Every option must say where it
comes from and what choosing it costs.
"""

from collections.abc import Sequence

import pytest
from pydantic import ValidationError

from sarna.rules.catalog import CATALOGUE
from sarna.rules.options import RuleOption
from sarna.rules.provenance import Provenance, ProvenanceKind
from sarna.thermo.tables import TABLES

RULES_WITH_OPTIONS = ["gc-content", "cpg-count", "thermodynamic-asymmetry"]


class TestEveryOptionIsExplained:
    @pytest.mark.parametrize("rule_id", RULES_WITH_OPTIONS)
    def test_options_carry_provenance_and_a_stated_trade_off(self, rule_id: str):
        for option in CATALOGUE[rule_id].options:
            assert str(option.provenance)
            assert option.note.strip()
            assert option.label.strip()

    @pytest.mark.parametrize("rule_id", RULES_WITH_OPTIONS)
    def test_exactly_one_option_is_the_default(self, rule_id: str):
        rule = CATALOGUE[rule_id]
        assert len([o for o in rule.options if o.is_default]) == 1
        assert rule.default_option is not None

    def test_a_rule_without_options_has_no_default_option(self):
        assert CATALOGUE["pos19-a"].default_option is None

    def test_an_unknown_option_id_raises(self):
        with pytest.raises(KeyError, match="no-such-option"):
            CATALOGUE["gc-content"].option("no-such-option")


class TestTheDefaultsFollowThePaper:
    def test_the_default_gc_option_is_the_published_band(self):
        default = CATALOGUE["gc-content"].default_option
        assert default is not None
        assert default.values["maximum-gc-fraction"] == 0.60
        assert default.provenance.kind is ProvenanceKind.CITATION

    def test_the_legacy_gc_option_is_offered_and_attributed_to_the_script(self):
        legacy = CATALOGUE["gc-content"].option("legacy-scanner")
        assert legacy.values["maximum-gc-fraction"] == 0.65
        assert legacy.provenance.kind is ProvenanceKind.LEGACY_SCRIPT
        assert legacy.is_default is False

    def test_the_zero_cpg_option_is_attributed_to_the_script_not_the_paper(self):
        legacy = CATALOGUE["cpg-count"].option("legacy-scanner")
        assert legacy.values["maximum-cpg-count"] == 0
        assert legacy.provenance.kind is ProvenanceKind.LEGACY_SCRIPT

    def test_an_option_beyond_the_literature_is_marked_unsourced(self):
        permissive = CATALOGUE["cpg-count"].option("permissive")
        assert permissive.provenance.kind is ProvenanceKind.UNSOURCED
        assert permissive.provenance.is_warning


class TestNearestNeighbourTableOptions:
    def test_each_option_names_a_table_that_exists(self):
        for option in CATALOGUE["thermodynamic-asymmetry"].options:
            assert option.values["nearest-neighbour-table"] in TABLES

    def test_the_default_is_the_rna_table(self):
        default = CATALOGUE["thermodynamic-asymmetry"].default_option
        assert default is not None
        assert default.values["nearest-neighbour-table"] == "rna-xia-1998"

    def test_the_default_option_agrees_with_the_recommended_table(self):
        # The engine's default must not drift from the table marked recommended.
        default = CATALOGUE["thermodynamic-asymmetry"].default_option
        assert default is not None
        chosen = TABLES[str(default.values["nearest-neighbour-table"])]
        assert chosen.is_recommended is True

    def test_all_three_tables_are_offered(self):
        offered = {
            o.values["nearest-neighbour-table"]
            for o in CATALOGUE["thermodynamic-asymmetry"].options
        }
        assert offered == set(TABLES)


class TestCoherenceIsEnforced:
    def test_the_default_option_must_match_the_parameter_defaults(self):
        rule = CATALOGUE["gc-content"]
        default = rule.default_option
        assert default is not None
        for name, value in default.values.items():
            assert rule.parameter(name).default == value

    def test_an_option_setting_an_unknown_parameter_is_refused(self):
        rule = CATALOGUE["gc-content"]
        bad = RuleOption(
            id="broken",
            label="Broken",
            provenance=Provenance.unsourced("test"),
            values={"no-such-parameter": 1.0},
            note="test",
        )
        with pytest.raises(ValidationError, match="unknown parameter"):
            rule.model_copy(update={"options": (*rule.options, bad)}).model_validate(
                rule.model_dump() | {"options": [*rule.options, bad]}
            )

    def test_an_option_value_outside_the_parameter_bounds_is_refused(self):
        rule = CATALOGUE["cpg-count"]
        bad = RuleOption(
            id="broken",
            label="Broken",
            provenance=Provenance.unsourced("test"),
            values={"maximum-cpg-count": 99},
            note="test",
        )
        with pytest.raises(ValidationError, match="above the maximum"):
            type(rule).model_validate(rule.model_dump() | {"options": [*rule.options, bad]})


class TestSerialisation:
    @pytest.mark.parametrize("rule_id", RULES_WITH_OPTIONS)
    def test_options_survive_a_round_trip(self, rule_id: str):
        rule = CATALOGUE[rule_id]
        assert type(rule).model_validate(rule.model_dump()) == rule

    def test_options_appear_in_json_with_their_provenance(self):
        dumped = CATALOGUE["gc-content"].model_dump(mode="json")
        legacy = next(o for o in dumped["options"] if o["id"] == "legacy-scanner")
        assert legacy["provenance"]["kind"] == "legacy-script"


class TestOptionValidation:
    def _option(self, **overrides: object) -> RuleOption:
        defaults: dict[str, object] = {
            "id": "an-option",
            "label": "An option",
            "provenance": Provenance.unsourced("test"),
            "values": {"maximum-cpg-count": 1},
            "note": "A stated trade-off.",
        }
        return RuleOption(**{**defaults, **overrides})  # type: ignore[arg-type]

    def test_ids_must_be_kebab_case(self):
        with pytest.raises(ValidationError, match="kebab-case"):
            self._option(id="An_Option")

    def test_a_label_is_required(self):
        with pytest.raises(ValidationError, match="label"):
            self._option(label="   ")

    def test_an_option_with_no_stated_trade_off_is_refused(self):
        with pytest.raises(ValidationError, match="trade-off"):
            self._option(note="   ")

    def test_an_option_setting_nothing_is_refused(self):
        with pytest.raises(ValidationError, match="no parameter values"):
            self._option(values={})


class TestRuleLevelOptionCoherence:
    def _with_options(self, rule_id: str, options: Sequence[object]) -> None:
        rule = CATALOGUE[rule_id]
        type(rule).model_validate(rule.model_dump() | {"options": options})

    def test_duplicate_option_ids_are_refused(self):
        rule = CATALOGUE["cpg-count"]
        duplicate = rule.options[0].model_dump()
        with pytest.raises(ValidationError, match="duplicate option ids"):
            self._with_options("cpg-count", [duplicate, duplicate])

    def test_two_defaults_are_refused(self):
        rule = CATALOGUE["cpg-count"]
        first = rule.options[0].model_dump()
        second = rule.options[1].model_dump() | {"is_default": True}
        with pytest.raises(ValidationError, match="exactly one must be marked default"):
            self._with_options("cpg-count", [first, second])

    def test_no_default_is_refused(self):
        rule = CATALOGUE["cpg-count"]
        options = [o.model_dump() | {"is_default": False} for o in rule.options]
        with pytest.raises(ValidationError, match="exactly one must be marked default"):
            self._with_options("cpg-count", options)

    def test_a_default_option_disagreeing_with_the_parameter_default_is_refused(self):
        # Catches the drift where someone edits a parameter default but not the option.
        rule = CATALOGUE["cpg-count"]
        options = [o.model_dump() for o in rule.options]
        options[0]["values"] = {"maximum-cpg-count": 3}
        with pytest.raises(ValidationError, match="these must agree"):
            self._with_options("cpg-count", options)
