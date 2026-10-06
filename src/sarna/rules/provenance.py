"""Where a design rule comes from, and how strong the evidence for it is.

ADR 0003 makes this mandatory metadata rather than documentation. The legacy Scanner buried its
science in bare module constants, so a reader could not tell that GC <= 0.65 was a local
deviation or that the position-7 rule had no source at all. Here a rule cannot be constructed
without stating where it came from, and ``unsourced`` is a permitted, *visible* answer — the
one thing that is not permitted is silence.
"""

from enum import Enum
from typing import Self

from pydantic import BaseModel, ConfigDict, model_validator


class EvidenceClass(Enum):
    """How well supported a rule is."""

    ESTABLISHED = "established"
    """Supported by primary literature beyond this project's own source paper."""

    EMPIRICAL = "empirical"
    """Observed or asserted by the source, without a mechanistic account."""

    EXPLORATORY = "exploratory"
    """Speculative. Carried for traceability, not because it is believed."""

    def __str__(self) -> str:
        """Return the class label, e.g. ``exploratory``."""
        return self.value


class ProvenanceKind(Enum):
    """The kind of origin a rule has."""

    CITATION = "citation"
    LEGACY_SCRIPT = "legacy-script"
    UNSOURCED = "unsourced"

    def __str__(self) -> str:
        """Return the kind label."""
        return self.value


class Provenance(BaseModel):
    """The traceable origin of a design rule.

    Construct through :meth:`citation`, :meth:`legacy_script` or :meth:`unsourced` rather than
    directly; each enforces the fields that kind requires.
    """

    model_config = ConfigDict(frozen=True)

    kind: ProvenanceKind
    source: str | None = None
    locator: str | None = None
    note: str | None = None

    @model_validator(mode="after")
    def _check_required_fields_for_kind(self) -> Self:
        """Enforce the fields each kind of provenance must carry."""
        if self.kind is ProvenanceKind.CITATION:
            if not (self.source or "").strip():
                raise ValueError("a citation requires a source, e.g. 'Wang et al. 2015'")
            if not (self.locator or "").strip():
                raise ValueError(
                    "a citation requires a locator, e.g. 'Table 1'; a bare reference cannot "
                    "be checked by a reader, which is the entire point of recording it"
                )
        elif self.kind is ProvenanceKind.LEGACY_SCRIPT:
            if not (self.locator or "").strip():
                raise ValueError(
                    "legacy-script provenance requires a locator naming the function or line "
                    "it came from, e.g. 'sarna_promoter_scanner.py:checkpos19'"
                )
        elif not (self.note or "").strip():
            raise ValueError(
                "unsourced provenance requires a note recording what was searched for and "
                "not found; 'unsourced' is a finding, not a shrug"
            )
        return self

    @classmethod
    def citation(cls, source: str, locator: str) -> "Provenance":
        """Return provenance citing published literature.

        Args:
            source: The work, e.g. ``"Wang et al. 2015"``.
            locator: Where in it, e.g. ``"Table 1"``. Required.
        """
        return cls(kind=ProvenanceKind.CITATION, source=source, locator=locator)

    @classmethod
    def legacy_script(cls, locator: str) -> "Provenance":
        """Return provenance for behaviour inherited from the original scripts.

        Args:
            locator: The function or line it came from.
        """
        return cls(kind=ProvenanceKind.LEGACY_SCRIPT, locator=locator)

    @classmethod
    def unsourced(cls, note: str) -> "Provenance":
        """Return provenance recording that no source could be traced.

        Args:
            note: What was searched for and not found.
        """
        return cls(kind=ProvenanceKind.UNSOURCED, note=note)

    @property
    def is_warning(self) -> bool:
        """Return whether the UI must render this as a warning.

        True only for :attr:`ProvenanceKind.UNSOURCED`. An unsourced rule may still be useful,
        but it must never be presented as established biology.
        """
        return self.kind is ProvenanceKind.UNSOURCED

    def __str__(self) -> str:
        """Return a human-readable, reader-checkable rendering."""
        if self.kind is ProvenanceKind.CITATION:
            return f"{self.source}, {self.locator}"
        if self.kind is ProvenanceKind.LEGACY_SCRIPT:
            return f"inherited from the legacy script ({self.locator})"
        return f"unsourced — {self.note}"
