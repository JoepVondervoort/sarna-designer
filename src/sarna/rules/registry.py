"""An immutable, ordered catalogue of design rules.

The registry is deliberately dumb: it holds rules and partitions them by kind, and knows
nothing about defaults, configuration or evaluation. A preset is a *subset* of a registry
(see :meth:`RuleRegistry.select`), which keeps "which rules exist" and "which rules this run
uses" as separate, inspectable things.
"""

from collections.abc import Iterable, Iterator

from sarna.rules.rule import Rule, RuleKind


class DuplicateRuleError(ValueError):
    """Raised when two rules share an id."""


class RuleNotFoundError(KeyError):
    """Raised when a rule id is not in the registry."""

    def __str__(self) -> str:
        """Return the message without ``KeyError``'s extra quoting."""
        return str(self.args[0]) if self.args else ""


class RuleRegistry:
    """An ordered, immutable collection of :class:`Rule`, keyed by id.

    Iteration order is construction order, because that is the order the UI presents rules in
    and it must not depend on dictionary internals or on how a preset was assembled.
    """

    __slots__ = ("_by_id",)

    def __init__(self, rules: Iterable[Rule]) -> None:
        """Build a registry from ``rules``.

        Args:
            rules: The rules to hold, in presentation order.

        Raises:
            DuplicateRuleError: If two rules share an id.
        """
        by_id: dict[str, Rule] = {}
        for rule in rules:
            if rule.id in by_id:
                raise DuplicateRuleError(
                    f"duplicate rule id {rule.id!r}; rule ids appear in exports and presets "
                    f"and must be unique"
                )
            by_id[rule.id] = rule
        self._by_id = by_id

    def __len__(self) -> int:
        """Return the number of rules held."""
        return len(self._by_id)

    def __iter__(self) -> Iterator[Rule]:
        """Iterate the rules in presentation order."""
        return iter(self._by_id.values())

    def __contains__(self, rule_id: object) -> bool:
        """Return whether ``rule_id`` is in the registry."""
        return rule_id in self._by_id

    def __getitem__(self, rule_id: str) -> Rule:
        """Return the rule with ``rule_id``.

        Raises:
            RuleNotFoundError: If no such rule exists.
        """
        try:
            return self._by_id[rule_id]
        except KeyError:
            raise RuleNotFoundError(f"no rule with id {rule_id!r} in this registry") from None

    def _of_kind(self, kind: RuleKind) -> tuple[Rule, ...]:
        return tuple(rule for rule in self if rule.kind is kind)

    @property
    def mandatory(self) -> tuple[Rule, ...]:
        """Return the rules that gate admission, in order."""
        return self._of_kind(RuleKind.MANDATORY)

    @property
    def preferences(self) -> tuple[Rule, ...]:
        """Return the scored preference rules, in order."""
        return self._of_kind(RuleKind.PREFERENCE)

    @property
    def tiebreaks(self) -> tuple[Rule, ...]:
        """Return the tie-breaking rules, in the order they are applied."""
        return self._of_kind(RuleKind.TIEBREAK)

    @property
    def max_preference_points(self) -> int:
        """Return the highest achievable preference score.

        This is the denominator in "meets 5 of 7 preferences". Because every preference is
        worth exactly one point (ADR 0006), it is simply how many there are.
        """
        return sum(rule.points for rule in self.preferences)

    @property
    def needing_warning(self) -> tuple[Rule, ...]:
        """Return the rules whose provenance the UI must flag as a warning."""
        return tuple(rule for rule in self if rule.needs_warning)

    def with_rule(self, rule: Rule) -> "RuleRegistry":
        """Return a new registry with ``rule`` appended.

        Raises:
            DuplicateRuleError: If the id is already present.
        """
        return RuleRegistry([*self, rule])

    def select(self, rule_ids: Iterable[str]) -> "RuleRegistry":
        """Return a new registry holding only ``rule_ids``, in the order given.

        This is how a preset is expressed: a named selection from the full catalogue.

        Raises:
            RuleNotFoundError: If any id is not in this registry.
        """
        return RuleRegistry([self[rule_id] for rule_id in rule_ids])
