"""Validated browser requests, resolved by the same configuration path as the library."""

from typing import Any, Literal, Self

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from sarna.design.config import DesignConfig
from sarna.design.presets import preset
from sarna.rules.parameters import ParameterValue
from sarna.rules.rule import RuleKind


class CustomReference(BaseModel):
    """A researcher-supplied genomic reference, without a verified-TSS claim."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9_.-]+$")
    name: str = Field(min_length=1, max_length=150)
    tss_position: int = Field(gt=0, strict=True)
    source: str = Field(min_length=1, max_length=300)
    note: str = Field(min_length=1, max_length=2000)


class ShortlistEntry(BaseModel):
    """A researcher selection with a reason, independent of scientific ranking."""

    model_config = ConfigDict(extra="forbid")
    candidate: str = Field(min_length=1, max_length=80)
    reason: str = Field(default="", max_length=2000)


class ShortlistRequest(BaseModel):
    """Replace one run's shortlist only if its saved revision still matches."""

    model_config = ConfigDict(extra="forbid")
    revision: int = Field(ge=0, strict=True)
    entries: list[ShortlistEntry] = Field(max_length=40_000)

    @model_validator(mode="after")
    def unique_candidates(self) -> Self:
        """Reject duplicate selections rather than silently losing a reason."""
        if len({entry.candidate for entry in self.entries}) != len(self.entries):
            raise ValueError("Each shortlisted candidate must occur once.")
        return self


class DesignRequest(BaseModel):
    """A deliberate input, promoter reference and editable scientific configuration."""

    model_config = ConfigDict(extra="forbid")

    source: Literal["gene", "sequence"] = "gene"
    symbol: str = Field(default="", max_length=80, pattern=r"^[A-Za-z0-9_.-]*$")
    sequence: str = Field(default="", max_length=40_000)
    transcript_id: str | None = Field(default=None, min_length=1, max_length=80)
    anchor_id: str | None = Field(default=None, min_length=1, max_length=80)
    custom_reference: CustomReference | None = None
    upstream: int = Field(default=2000, ge=0, le=20_000)
    downstream: int = Field(default=500, ge=0, le=20_000)
    tss_start: int | None = None
    cell_line: str = Field(default="Not specified", max_length=150)
    preset_id: str = "wang-2015"
    overrides: dict[str, dict[str, ParameterValue]] = Field(default_factory=dict)
    roles: dict[str, RuleKind | Literal["off"]] = Field(default_factory=dict)
    options: dict[str, str] = Field(default_factory=dict)
    offline: bool = False
    regulatory: bool = True
    motifs: bool = True
    motif_shuffles: int = Field(default=5, ge=1, le=100)
    chromatin: bool = False

    @field_validator("symbol", mode="before")
    @classmethod
    def normalise_symbol(cls, value: Any) -> Any:  # noqa: ANN401
        """Save the same human symbol used for lookup, independent of input casing."""
        return value.strip().upper() if isinstance(value, str) else value

    @model_validator(mode="after")
    def check_input(self) -> Self:
        """Reject ambiguous sources and unsupported workload sizes before lookup."""
        if self.source == "gene" and not self.symbol.strip():
            raise ValueError("Enter a human gene symbol.")
        if self.source == "sequence" and not self.sequence.strip():
            raise ValueError("Paste a promoter DNA sequence.")
        if self.source == "gene" and not 19 <= self.upstream + self.downstream <= 20_000:
            raise ValueError("The supported promoter span is 19 to 20,000 bp.")
        references = (self.transcript_id, self.anchor_id, self.custom_reference)
        if any(isinstance(reference, str) and not reference.strip() for reference in references):
            raise ValueError("A selected promoter reference must have a nonempty identifier.")
        if sum(reference is not None for reference in references) > 1:
            raise ValueError(
                "Choose one transcript, one literature anchor or one custom reference."
            )
        if self.source != "gene" and self.custom_reference is not None:
            raise ValueError("A custom genomic reference requires a confirmed human gene.")
        if self.tss_start == 0:
            raise ValueError("TSS positions have no zero: use -1 or +1.")
        return self

    def resolve(self) -> DesignConfig:
        """Validate roles, options and parameter relationships in the Python engine."""
        config = preset(self.preset_id).resolve(
            self.overrides,
            role_overrides={key: None if value == "off" else value
                            for key, value in self.roles.items()},
            option_ids=self.options,
        )
        for rule_id, names in (
            ("tss-window", ("window-start-bp", "window-end-bp")),
            ("tss-preferred-window", ("preferred-window-start-bp", "preferred-window-end-bp")),
        ):
            if rule_id in config:
                for name in names:
                    if config.rule(rule_id).integer(name) == 0:
                        raise ValueError(f"{name}: there is no TSS position 0; use -1 or +1.")
        return config


class BatchRequest(BaseModel):
    """Separate designs for deliberately selected references of the same gene."""

    model_config = ConfigDict(extra="forbid")

    label: str = Field(default="", max_length=150)
    requests: list[DesignRequest] = Field(min_length=1, max_length=32)

    @model_validator(mode="after")
    def check_references(self) -> Self:
        """Keep every batch member tied to an explicit reference and one input gene."""
        symbols = {request.symbol for request in self.requests}
        if any(request.source != "gene" for request in self.requests) or len(symbols) != 1:
            raise ValueError("A promoter batch requires one human gene shared by every request.")
        for request in self.requests:
            if sum(reference is not None for reference in (
                request.transcript_id, request.anchor_id, request.custom_reference,
            )) != 1:
                raise ValueError("Choose an explicit promoter reference for every batch request.")
        return self
