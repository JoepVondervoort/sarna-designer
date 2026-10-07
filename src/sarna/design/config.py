"""The resolved design configuration: exactly which rules a run applies, and with what values.

ADR 0010 makes a preset a *selection* of rules rather than a mutation of them, and this is what a
selection resolves to. A `DesignConfig` holds one :class:`~sarna.rules.resolved.ResolvedRule` per
selected rule, so by the time the engine runs there is no further question about what applies.

Two distinctions this type exists to keep straight:

**Not selected is not the same as not applicable.** A rule the configuration omits is never
evaluated and produces no result. A rule that *is* selected but whose inputs are missing reports
``not-applicable`` (ADR 0007). "This preset does not use the TSS window" and "the TSS window could
not be checked because the promoter has no TSS" are different statements about a candidate and the
interface must not merge them.

**The configuration is the reproducible record, not the preset name.** A preset's meaning can
change between versions of this software, so :meth:`DesignConfig.to_dict` exports the resolved
rules with their values, provenance and evidence class. A stored analysis referencing only
``"wang-2015"`` would not be reproducible.
"""

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Any

from sarna.rules.catalog import CATALOGUE
from sarna.rules.parameters import ParameterValue
from sarna.rules.resolved import ResolvedRule
from sarna.rules.rule import RuleKind
from sarna.sequence.cpg import CpGIslandCriteria

ISLAND_RULE = "not-in-cpg-island"
LENGTH_RULE = "target-length-19"


@dataclass(frozen=True, slots=True)
class DesignConfig:
    """Every rule a run applies, resolved, in catalogue order."""

    rules: tuple[ResolvedRule, ...]
    preset_id: str | None = None
    preset_name: str | None = None
    preset_divergences: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Reject a selection naming the same rule twice.

        Raises:
            ValueError: If a rule id appears more than once, which would evaluate it twice and
                double-count its preference point.
        """
        seen = [resolved.rule.id for resolved in self.rules]
        duplicates = sorted({rule_id for rule_id in seen if seen.count(rule_id) > 1})
        if duplicates:
            raise ValueError(
                f"rule(s) {', '.join(duplicates)} appear more than once in this configuration; "
                f"a rule evaluated twice would count its preference point twice"
            )

    @classmethod
    def of(
        cls,
        rule_ids: Iterable[str],
        overrides: Mapping[str, Mapping[str, ParameterValue]] | None = None,
        preset_id: str | None = None,
        preset_name: str | None = None,
        preset_divergences: tuple[str, ...] = (),
    ) -> "DesignConfig":
        """Resolve a selection of catalogued rule ids into a configuration.

        Args:
            rule_ids: The rules to apply, by id.
            overrides: Parameter overrides, keyed by rule id then parameter name.
            preset_id: The preset this came from, if any.
            preset_name: Its human-readable name.
            preset_divergences: Known ways it departs from what it reproduces.

        Returns:
            The resolved configuration.

        Raises:
            KeyError: If a rule id is not in the catalogue.
            ValueError: If an override value is invalid for its parameter.
        """
        by_rule = overrides or {}
        return cls(
            rules=tuple(
                ResolvedRule(CATALOGUE[rule_id], by_rule.get(rule_id, {}))
                for rule_id in rule_ids
            ),
            preset_id=preset_id,
            preset_name=preset_name,
            preset_divergences=preset_divergences,
        )

    def __len__(self) -> int:
        """Return how many rules this configuration applies."""
        return len(self.rules)

    def __iter__(self) -> Iterator[ResolvedRule]:
        """Iterate the resolved rules in catalogue order."""
        return iter(self.rules)

    def __contains__(self, rule_id: object) -> bool:
        """Return whether a rule id is selected by this configuration."""
        return any(resolved.rule.id == rule_id for resolved in self.rules)

    def rule(self, rule_id: str) -> ResolvedRule:
        """Return a selected rule by id.

        Raises:
            KeyError: If the rule is not selected. The message distinguishes "no such rule" from
                "not part of this configuration", because they call for different fixes.
        """
        for resolved in self.rules:
            if resolved.rule.id == rule_id:
                return resolved
        known = "an unknown rule" if rule_id not in CATALOGUE else "not selected"
        raise KeyError(
            f"rule {rule_id!r} is {known} in configuration "
            f"{self.preset_id or '(no preset)'!r}"
        )

    def _of_kind(self, kind: RuleKind) -> tuple[ResolvedRule, ...]:
        return tuple(resolved for resolved in self.rules if resolved.rule.kind is kind)

    @property
    def mandatory(self) -> tuple[ResolvedRule, ...]:
        """Return the gates. Failing any excludes the candidate."""
        return self._of_kind(RuleKind.MANDATORY)

    @property
    def preferences(self) -> tuple[ResolvedRule, ...]:
        """Return the scored preferences, each worth one point."""
        return self._of_kind(RuleKind.PREFERENCE)

    @property
    def tiebreaks(self) -> tuple[ResolvedRule, ...]:
        """Return the tie-break rules, in the order they are applied."""
        return self._of_kind(RuleKind.TIEBREAK)

    @property
    def max_preference_points(self) -> int:
        """Return the denominator in "meets 5 of 7 preferences"."""
        return len(self.preferences)

    @property
    def target_length(self) -> int:
        """Return the candidate length the scan generates.

        Needed before any rule runs, so it is surfaced here rather than dug out of the rule
        list by the generator.

        Raises:
            KeyError: If the configuration does not select the target-length rule, without which
                there is no defined candidate length.
        """
        return self.rule(LENGTH_RULE).integer("target-length")

    @property
    def cpg_island_criteria(self) -> CpGIslandCriteria | None:
        """Return the island definition to detect with, or ``None`` if islands are not used.

        ``None`` means the configuration does not select the island rule, so detection should not
        run at all. It is not a signal to fall back on a default definition: the engine does not
        have an opinion about islands that the configuration has not stated.
        """
        if ISLAND_RULE not in self:
            return None
        resolved = self.rule(ISLAND_RULE)
        return CpGIslandCriteria(
            minimum_length=resolved.integer("island-minimum-length"),
            minimum_gc_fraction=resolved.number("island-minimum-gc-fraction"),
            minimum_observed_expected=resolved.number(
                "island-minimum-observed-expected-cpg"
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the configuration as JSON-serialisable data, for export and provenance.

        Records the resolved parameter values rather than the preset name, because a preset's
        meaning can change between versions. Each entry is the full rule *record* — its rationale,
        consequence, structured provenance, every parameter's unit and legal range, and the options
        not chosen with their own provenance — so a reader can judge the configuration without this
        codebase to hand (ADR 0003).

        Each rule appears once here, and a candidate's results reference it by ``id``. Repeating
        the record on every result would be the same text thousands of times over: a 2.5 kb
        promoter yields upwards of 2,400 candidates against 26 rules, and the page has to parse
        whatever it is given.
        """
        return {
            "preset": {
                "id": self.preset_id,
                "name": self.preset_name,
                "divergences": list(self.preset_divergences),
            },
            "rules": [resolved.rule.to_dict(resolved.values) for resolved in self.rules],
        }
