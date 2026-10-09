"""The shipped presets (ADR 0002, ADR 0010).

A preset names the rules it applies and the values it overrides. It never edits a rule, so the
catalogue stays the single immutable statement of what rules exist and a preset is a projection
of it.

``wang-2015``
    The default. Every rule enabled by default in the catalogue, with default parameters.
    ``pos7-u`` is excluded, being unsourced.

``strict``
    Wang et al. 2015 with its softer statements promoted to requirements. Fewer candidates,
    all conformant on every count the paper makes.

``exploratory``
    Wang et al. 2015 plus optional criteria from the siRNA literature and thesis chapter 4.
    Same candidate set as the default, reordered. Its score denominator differs, so its
    numbers are not comparable with any other preset's.

``legacy-scanner``
    Reproduces the original Promoter Scanner's filtering closely, including its defects.
    **Not bit-for-bit**, and the two exceptions are declared on the preset itself (in
    :attr:`Preset.divergences`) rather than only in an ADR, so that they reach an export and the
    interface rather than living where nobody reads them at the bench.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field

from sarna.design.config import DesignConfig
from sarna.rules.catalog import CATALOGUE
from sarna.rules.parameters import ParameterValue
from sarna.rules.rule import RuleKind

ParameterOverrides = Mapping[str, Mapping[str, ParameterValue]]


@dataclass(frozen=True, slots=True)
class Preset:
    """A named selection of rules with their overrides."""

    id: str
    name: str
    description: str
    rule_ids: tuple[str, ...]
    overrides: ParameterOverrides = field(default_factory=dict)
    is_default: bool = False
    divergences: tuple[str, ...] = ()
    """Known ways this preset departs from what it claims to reproduce. Exported and displayed."""

    def resolve(
        self,
        user_overrides: ParameterOverrides | None = None,
        *,
        role_overrides: Mapping[str, RuleKind | None] | None = None,
        option_ids: Mapping[str, str] | None = None,
    ) -> DesignConfig:
        """Resolve this preset, optionally layering user overrides over its own.

        Args:
            user_overrides: Further parameter overrides, keyed by rule id then parameter name.
                Merged per rule, so overriding one parameter leaves the preset's other values
                for that rule in place.
            role_overrides: Custom gate/preference roles; ``None`` disables a rule. An optional
                catalogue rule can be added by explicitly assigning it a supported role.
            option_ids: Deliberately chosen sourced alternatives; manual values take precedence.

        Returns:
            The resolved configuration.

        Raises:
            KeyError: If an override names a rule this preset does not select. Ignoring it
                would let a user believe they had changed something.
            ValueError: If an override value is invalid for its parameter.
        """
        roles = role_overrides or {}
        options = option_ids or {}
        for rule_id in roles:
            CATALOGUE[rule_id]  # Reject unknown ids, including attempts to disable them.
        selected = (*self.rule_ids, *(
            rule.id for rule in CATALOGUE
            if rule.id not in self.rule_ids and rule.id in roles
        ))
        for rule_id in set(user_overrides or {}) | set(options):
            if rule_id not in selected or (rule_id in roles and roles[rule_id] is None):
                raise KeyError(
                    f"cannot override rule {rule_id!r}: it is not selected by preset "
                    f"{self.id!r} or is disabled, so changing its parameters would have no effect"
                )
        merged: dict[str, dict[str, ParameterValue]] = {
            rule_id: dict(values) for rule_id, values in self.overrides.items()
        }
        for rule_id, option_id in options.items():
            for name in CATALOGUE[rule_id].option(option_id).values:
                merged.get(rule_id, {}).pop(name, None)
        for rule_id, values in (user_overrides or {}).items():
            merged.setdefault(rule_id, {}).update(values)
        baseline = DesignConfig.of(self.rule_ids, self.overrides).rules
        return DesignConfig.of(
            selected,
            merged,
            preset_id=self.id,
            preset_name=self.name,
            preset_divergences=self.divergences,
            role_overrides=roles,
            option_ids=options,
            baseline_rules=baseline,
        )


_DEFAULT_RULE_IDS = tuple(
    rule.id for rule in CATALOGUE if rule.enabled_by_default
)

WANG_2015 = Preset(
    id="wang-2015",
    name="Wang et al. 2015 (default)",
    description=(
        "Conformant with Wang et al. 2015 Table 1. Every rule the catalogue enables by "
        "default, at its published value. A literature-based starting point for research; "
        "conformance does not establish efficacy for a new gene or cell type."
    ),
    rule_ids=_DEFAULT_RULE_IDS,
    is_default=True,
)


_LEGACY_EXCLUDED = frozenset(
    {
        # The script scans a pasted string and has no concept of a transcription start site.
        "tss-window",
        "tss-preferred-window",
        # Defect B8: the script's `inisland` call is commented out, so the rule never ran.
        "not-in-cpg-island",
        # Absent from the script entirely.
        "no-simple-repeats",
    }
)

LEGACY_SCANNER = Preset(
    id="legacy-scanner",
    name="Original Promoter Scanner (reproduction)",
    description=(
        "Reproduces the original script's filtering, defects included, so historical designs "
        "can be re-derived and compared. Not for new work: it uses a GC band wider than the "
        "paper's, a nearest-neighbour table with three stacks mis-assigned, and an unsourced "
        "position-7 rule."
    ),
    rule_ids=tuple(
        rule_id
        for rule_id in (
            *(rule.id for rule in CATALOGUE if rule.enabled_by_default),
            "pos7-u",  # Off by default, but one of the script's two dominant ranking terms.
        )
        if rule_id not in _LEGACY_EXCLUDED
    ),
    overrides={
        "gc-content": {"maximum-gc-fraction": 0.65},
        "cpg-count": {"maximum-cpg-count": 0},
        "thermodynamic-asymmetry": {"nearest-neighbour-table": "legacy-dna-kj"},
        # Ordering by a different parameter set than the one that admitted the candidate would
        # be incoherent, so the tie-break follows the gate.
        "asymmetry-margin": {"nearest-neighbour-table": "legacy-dna-kj"},
    },
    divergences=(
        "Position 1 G/C is scored as a preference here, as Wang et al. have it, but the "
        "original script enforces it as a requirement through its `endbase` check (finding "
        "B16). This preset therefore admits candidates the original rejected (every candidate "
        "whose position 1 is A or T) so its candidate set is wider than the original's, not "
        "merely ordered differently.",
        "Candidates are ranked by counting criteria met (ADR 0006), not by the original's "
        "weighted sum of 0.01 x stability difference + 1000 x position 18 + 1000 x position 7 "
        "+ 1 x flank score. Filtering is reproduced closely; ordering is not reproduced at all.",
    ),
)


STRICT = Preset(
    id="strict",
    name="Strict (narrower location criteria)",
    description=(
        "Wang et al. 2015 with a narrower -500 to -200 bp admission window and a more "
        "permissive CpG-island detector, which excludes more candidate regions. Sequence "
        "criteria retain Wang's values. These stricter location choices do not establish "
        "greater confidence in activation."
    ),
    rule_ids=_DEFAULT_RULE_IDS,
    overrides={
        # The paper permits -1000 to -100 but reports activity concentrated in -500 to -200.
        # Strict gates on the reported region instead of merely preferring it.
        "tss-window": {"window-start-bp": -500, "window-end-bp": -200},
        # The looser island definition finds more and shorter islands, so excluding any overlap
        # with one rejects more candidates. Gardiner-Garden & Frommer is the stricter choice
        # here precisely because it is the more permissive island detector.
        "not-in-cpg-island": {
            "island-minimum-length": 200,
            "island-minimum-gc-fraction": 0.50,
            "island-minimum-observed-expected-cpg": 0.60,
        },
    },
    divergences=(
        "The -500 to -200 bp admission window comes from the paper's narrative text, not from "
        "Table 1, which permits -1000 to -100. Candidates Wang et al. would admit are excluded "
        "here, so this preset is stricter than the paper rather than conformant with it.",
        "Only location criteria are tightened. The sequence criteria stay at Wang's values "
        "because the one empirical dataset available supports them exactly as published: "
        "GC 40-60 %, at most one CpG, and no run over four bases each exclude none of the 14 "
        "measured activators (docs/scientific/empirical-validation.md).",
        "Two tightenings that look obvious were measured and rejected. Requiring zero CpG "
        "excludes 4 of the 14 known activators, and activators in fact carry slightly more CpG "
        "than non-activators, not less. Requiring an asymmetry margin of 0.5 kcal/mol excludes "
        "a further 2. Both are available as manual overrides; neither is in this preset.",
    ),
)

_EXPLORATORY_ADDED = (
    # Best supported of the optional rules: strand-selection mechanism, three agreeing sources.
    "au-rich-guide-5-prime",
    "no-gc-stretch",
    # Reynolds' positional set. pos10-u is excluded: it contradicts pos10-g below, and its
    # mechanism is cleavage, which activation does not use.
    "pos3-a",
    "no-g-at-13",
    # Chapter 4's hypotheses.
    "pos10-g",
    "central-stability",
)

EXPLORATORY = Preset(
    id="exploratory",
    name="Exploratory (literature and chapter 4 hypotheses)",
    description=(
        "Wang et al. 2015 plus optional criteria from the siRNA design literature (Reynolds et "
        "al. 2004, Ui-Tei et al. 2004) and from thesis chapter 4's positional analysis. Admits "
        "exactly the same candidates as the default (every added rule is a preference, so none "
        "can exclude anything) and reorders them. Use it to see whether the hypotheses would "
        "have changed your shortlist; do not use it to justify a design."
    ),
    rule_ids=(*_DEFAULT_RULE_IDS, *_EXPLORATORY_ADDED),
    overrides={
        # Reynolds' form of the A/U-rich criterion rather than Ui-Tei's. Ui-Tei's threshold is
        # missed by 7 of the 14 measured activators and Reynolds' by 1, so at Ui-Tei's setting
        # the rule would mostly penalise the sequences it is meant to reward.
        "au-rich-guide-5-prime": {"window-length-nt": 5, "minimum-au-count": 3},
    },
    divergences=(
        "Four of the six added rules are EXPLORATORY: pos3-a and no-g-at-13 have no established "
        "mechanism, and pos10-g and central-stability come from an analysis that does not "
        "survive correction for multiple comparisons "
        "(docs/scientific/positional-motif-analysis.md). The ranking this preset produces is "
        "a hypothesis about what matters, not a better estimate than the default's.",
        "The preference denominator rises from 7 to 13, so scores from this preset are not "
        "comparable with scores from any other. Compare orderings, never numbers.",
        "central-stability and au-rich-guide-5-prime both overlap the mandatory thermodynamic "
        "asymmetry rule, which measures the same physics continuously. A candidate can "
        "therefore gain points for a property it was already admitted on.",
    ),
)


PRESETS: Mapping[str, Preset] = {
    preset.id: preset
    for preset in (WANG_2015, STRICT, EXPLORATORY, LEGACY_SCANNER)
}


def preset(preset_id: str) -> Preset:
    """Return a shipped preset by id.

    Raises:
        KeyError: If no such preset exists. The message lists those that do.
    """
    try:
        return PRESETS[preset_id]
    except KeyError:
        raise KeyError(
            f"unknown preset {preset_id!r}; available presets are "
            f"{', '.join(sorted(PRESETS))}"
        ) from None


__all__ = [
    "EXPLORATORY",
    "LEGACY_SCANNER",
    "PRESETS",
    "STRICT",
    "WANG_2015",
    "Preset",
    "preset",
]
