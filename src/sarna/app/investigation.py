"""Validated map annotations, exact sequence searches and interval measurements."""

import re
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, StrictInt, model_validator

from sarna.app.sequence_import import reverse


class Annotation(BaseModel):
    """A researcher region in zero-based, half-open promoter coordinates."""

    model_config = ConfigDict(extra="forbid")
    name: str = Field(min_length=1, max_length=150)
    start: Annotated[StrictInt, Field(ge=0)]
    end: Annotated[StrictInt, Field(gt=0)]
    source: str = Field(default="", max_length=2000)
    orientation: Literal["same", "reverse complement", "region"] = "region"

    @model_validator(mode="after")
    def ordered(self) -> "Annotation":
        """Require a nonempty region."""
        if self.end <= self.start or not self.name.strip():
            raise ValueError("Name the region and put its end after its start.")
        return self


class InvestigationRequest(BaseModel):
    """Revision-protected researcher annotations and motif display choices."""

    model_config = ConfigDict(extra="forbid")
    revision: Annotated[StrictInt, Field(ge=0)]
    annotations: list[Annotation] = Field(default_factory=list, max_length=200)
    motifs: list[str] | None = Field(default=None, max_length=100)


class SequenceSearch(BaseModel):
    """An exact search of the loaded promoter, with explicit orientation."""

    model_config = ConfigDict(extra="forbid")
    sequence: str = Field(min_length=4, max_length=20000)
    orientation: Literal["both", "same", "reverse complement"] = "both"

    def find(self, promoter: str) -> list[dict[str, Any]]:
        """Return every overlapping occurrence, deduplicating palindromes."""
        dna = re.sub(r"\s+", "", self.sequence).upper()
        if not re.fullmatch(r"[ACGTU]{4,10000}", dna) or ("T" in dna and "U" in dna):
            raise ValueError("Use 4-10000 DNA or RNA bases, without mixed T and U.")
        dna = dna.replace("U", "T")
        targets = [(dna, "same"), (reverse(dna), "reverse complement")]
        found: list[dict[str, Any]] = []
        seen = set()
        for target, orientation in targets:
            if self.orientation not in ("both", orientation) or target in seen:
                continue
            seen.add(target)
            offset = promoter.find(target)
            while offset >= 0:
                if len(found) >= 500:
                    raise ValueError("More than 500 matches; use a longer sequence.")
                found.append({"start": offset, "end": offset + len(target),
                              "orientation": orientation})
                offset = promoter.find(target, offset + 1)
        return sorted(found, key=lambda item: (item["start"], item["orientation"]))
