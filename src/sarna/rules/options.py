"""Named alternatives a user may choose for a rule.

Several rules have more than one defensible setting: the paper says one thing, the original
script did another, and a third convention exists in the wider literature. Rather than hide the
alternatives inside preset definitions, each rule publishes them as options the user can pick,
**each carrying its own provenance**.

The point is that choosing the legacy GC band should not feel like editing a number out of
nowhere. It should read as "this is where that value comes from, and this is what picking it
costs you".
"""

import re
from collections.abc import Mapping
from typing import Any, Self

from pydantic import BaseModel, ConfigDict, model_validator

from sarna.rules.parameters import ParameterValue
from sarna.rules.provenance import Provenance

_KEBAB_CASE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")


class RuleOption(BaseModel):
    """One selectable setting for a rule, with where it comes from and what it costs."""

    model_config = ConfigDict(frozen=True)

    id: str
    """Stable kebab-case identifier, unique within its rule."""

    label: str
    """Short human-readable name, as shown in the UI."""

    provenance: Provenance
    """Where this particular setting comes from. Required, like any other rule metadata."""

    values: Mapping[str, ParameterValue]
    """Parameter values this option applies, keyed by parameter name."""

    note: str
    """Why a user might choose it, and what it costs them. Plain language."""

    is_default: bool = False
    """Whether this is the setting used unless the user picks another."""

    @model_validator(mode="after")
    def _check_identifier_and_text(self) -> Self:
        """Validate the identifier format and require the explanatory text."""
        if not _KEBAB_CASE.match(self.id):
            raise ValueError(f"option id {self.id!r} must be kebab-case, e.g. 'legacy-scanner'")
        if not self.label.strip():
            raise ValueError(f"option {self.id!r} must have a label")
        if not self.note.strip():
            raise ValueError(
                f"option {self.id!r} must explain why a user would choose it; an option with "
                f"no stated trade-off is a trap"
            )
        if not self.values:
            raise ValueError(f"option {self.id!r} sets no parameter values")
        return self

    def to_dict(self) -> dict[str, Any]:
        """Return the option as JSON-serialisable data.

        Exported even when it was not chosen: a reader judging a threshold needs to know what else
        was on offer and what the alternative costs, which is what ``note`` carries.
        """
        return {
            "id": self.id,
            "label": self.label,
            "provenance": self.provenance.to_dict(),
            "values": dict(self.values),
            "note": self.note,
            "is_default": self.is_default,
        }
