"""The two shipped presets (ADR 0002, ADR 0010).

A preset names the rules it applies and the values it overrides. It never edits a rule, so the
catalogue stays the single immutable statement of what rules exist and a preset is a projection
of it.

``wang-2015``
    The default. Every rule enabled by default in the catalogue, with default parameters.
    ``pos7-u`` is excluded, being unsourced.

``legacy-scanner``
    Reproduces the original Promoter Scanner's filtering closely, including its defects.
    **Not bit-for-bit**, and the two exceptions are declared on the preset itself — in
    :attr:`Preset.divergences` — rather than only in an ADR, so that they reach an export and the
    interface rather than living where nobody reads them at the bench.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field

from sarna.design.config import DesignConfig
from sarna.rules.catalog import CATALOGUE
from sarna.rules.parameters import ParameterValue

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

    def resolve(self, user_overrides: ParameterOverrides | None = None) -> DesignConfig:
        """Resolve this preset, optionally layering user overrides over its own.

        Args:
            user_overrides: Further parameter overrides, keyed by rule id then parameter name.
                Merged per rule, so overriding one parameter leaves the preset's other values
                for that rule in place.

        Returns:
            The resolved configuration.

        Raises:
            KeyError: If an override names a rule this preset does not select. Ignoring it
                would let a user believe they had changed something.
            ValueError: If an override value is invalid for its parameter.
        """
        merged: dict[str, dict[str, ParameterValue]] = {
            rule_id: dict(values) for rule_id, values in self.overrides.items()
        }
        for rule_id, values in (user_overrides or {}).items():
            if rule_id not in self.rule_ids:
                raise KeyError(
                    f"cannot override rule {rule_id!r}: it is not selected by preset "
                    f"{self.id!r}, so changing its parameters would have no effect"
                )
            merged.setdefault(rule_id, {}).update(values)
        return DesignConfig.of(
            self.rule_ids,
            merged,
            preset_id=self.id,
            preset_name=self.name,
            preset_divergences=self.divergences,
        )


_DEFAULT_RULE_IDS = tuple(
    rule.id for rule in CATALOGUE if rule.enabled_by_default
)

WANG_2015 = Preset(
    id="wang-2015",
    name="Wang et al. 2015 (default)",
    description=(
        "Conformant with Wang et al. 2015 Table 1. Every rule the catalogue enables by "
        "default, at its published value. Recommended for new designs."
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
        "B16). This preset therefore admits candidates the original rejected — every candidate "
        "whose position 1 is A or T — so its candidate set is wider than the original's, not "
        "merely ordered differently.",
        "Candidates are ranked by counting criteria met (ADR 0006), not by the original's "
        "weighted sum of 0.01 x stability difference + 1000 x position 18 + 1000 x position 7 "
        "+ 1 x flank score. Filtering is reproduced closely; ordering is not reproduced at all.",
    ),
)


PRESETS: Mapping[str, Preset] = {
    preset.id: preset for preset in (WANG_2015, LEGACY_SCANNER)
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


__all__ = ["LEGACY_SCANNER", "PRESETS", "WANG_2015", "Preset", "preset"]
