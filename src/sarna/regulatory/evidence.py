"""What is actually known about a regulatory feature, and in which cell type.

A tool that draws a transcription factor on a promoter is making a claim, and there are three very
different claims it might be making:

1. **The sequence could bind this factor.** A motif match. Computed from sequence alone, so it is
   cell-type agnostic by construction and says nothing about any particular cell.
2. **This region is regulatory in some cell type.** An Ensembl Regulatory Build annotation,
   derived from epigenomic data aggregated over a panel of cell types. It does not say which.
3. **This factor was observed bound here, in this cell line.** A ChIP-seq peak. The only one of
   the three that answers "would this bind in my cells?".

Drawing all three the same way would be the most misleading thing this tool could do, because the
first is cheap and abundant and the third is expensive and rare. A GC-rich promoter contains an
SP1 motif match every few dozen bases; whether SP1 is bound at any of them in AF22 neuroepithelial
stem cells is a separate question, usually unanswered.

So :class:`EvidenceKind` is required on every annotation, and the model enforces the one rule that
keeps the kinds from blurring: **a motif prediction may not name a cell type.** Attaching one would
imply a measurement that was never made.

The fourth state matters as much as the three claims. :class:`CellTypeStatus` distinguishes
"observed bound", "assayed and not bound" and **"not assayed"**. Most factor-by-cell-type
combinations are in the third state, and a UI that renders absence of data the same as absence of
binding converts ignorance into a negative finding.
"""

import os
import tomllib
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from importlib import resources
from pathlib import Path
from typing import Any

from sarna.rules.provenance import Provenance


class EvidenceKind(Enum):
    """What kind of claim an annotation makes."""

    MOTIF_PREDICTION = "motif-prediction"
    """The sequence matches a binding motif. Cell-type agnostic; abundant; many false positives."""

    AGGREGATED_ANNOTATION = "aggregated-annotation"
    """Annotated as regulatory from data pooled across cell types. Does not say which."""

    CELL_TYPE_OBSERVATION = "cell-type-observation"
    """Observed experimentally in one named cell type. The only kind that answers "in my cells"."""

    def __str__(self) -> str:
        """Return the kind label."""
        return self.value

    @property
    def is_cell_type_specific(self) -> bool:
        """Return whether this kind of evidence is about a particular cell type."""
        return self is EvidenceKind.CELL_TYPE_OBSERVATION

    @property
    def confidence_note(self) -> str:
        """Return the caveat a UI must display beside this kind of evidence."""
        if self is EvidenceKind.MOTIF_PREDICTION:
            return (
                "Sequence match only. Most motif matches are not bound in any given cell type; "
                "treat this as 'the site exists', not 'the factor is there'."
            )
        if self is EvidenceKind.AGGREGATED_ANNOTATION:
            return (
                "Annotated from data pooled across a panel of cell types, which may not include "
                "yours. Absence here is weak evidence of absence."
            )
        return "Observed in the named cell type by experiment."


class CellTypeStatus(Enum):
    """Whether a factor is known to bind in a particular cell type."""

    OBSERVED_BOUND = "observed-bound"
    ASSAYED_NOT_BOUND = "assayed-not-bound"
    NOT_ASSAYED = "not-assayed"
    """No experiment has looked. The commonest state, and not the same as 'not bound'."""

    def __str__(self) -> str:
        """Return the status label."""
        return self.value

    @property
    def is_known(self) -> bool:
        """Return whether anything was actually measured."""
        return self is not CellTypeStatus.NOT_ASSAYED


class DataAxis(Enum):
    """A kind of evidence a cell line may or may not have, and what it answers."""

    EXPRESSION = "expression"
    """Is the target gene transcribed here at all? The first question to ask."""

    CHROMATIN = "chromatin"
    """Is the promoter accessible here? Histone marks."""

    TF_OCCUPANCY = "tf-occupancy"
    """Which factors are actually bound here? ChIP-seq."""

    def __str__(self) -> str:
        """Return the axis label."""
        return self.value

    @property
    def toml_key(self) -> str:
        """Return the key this axis is declared under in a cell-line file.

        Separate from :attr:`value` because the label uses a hyphen and TOML keys here use an
        underscore. Reading one as the other silently drops the axis, which is how this started:
        ``tf_occupancy`` was declared, looked for as ``tf-occupancy``, and never found.
        """
        return self.value.replace("-", "_")

    @property
    def question(self) -> str:
        """Return the question this axis answers, for a UI to label it with."""
        return {
            DataAxis.EXPRESSION: "Is this gene transcribed in these cells?",
            DataAxis.CHROMATIN: "Is this promoter accessible in these cells?",
            DataAxis.TF_OCCUPANCY: "Which factors are bound here in these cells?",
        }[self]


@dataclass(frozen=True, slots=True)
class DataAvailability:
    """Whether one axis of evidence exists for one cell line, and what to do about it.

    ``available = False`` is a first-class answer that must stay visible. An empty panel means
    nobody looked, which is not the same as nothing being there, and the ``note`` says which.
    """

    axis: DataAxis
    available: bool
    note: str
    source: str | None = None
    """Where the data is or would come from, e.g. ``ENCODE``. ``None`` if there is no route."""

    @property
    def is_reachable(self) -> bool:
        """Return whether a source exists that could supply this, even if unwired."""
        return self.available or self.source is not None

    def to_dict(self) -> dict[str, object]:
        """Return the availability as JSON-serialisable data."""
        return {
            "axis": str(self.axis),
            "question": self.axis.question,
            "available": self.available,
            "reachable": self.is_reachable,
            "source": self.source,
            "note": self.note,
        }


@dataclass(frozen=True, slots=True)
class CellType:
    """A cell line or primary cell type an annotation can be scoped to."""

    id: str
    """Stable identifier, e.g. ``hela`` or ``af22``."""

    name: str
    description: str
    """What the cell type is, and why it is in this project."""

    is_proliferating: bool
    """Whether the cells divide. Bears directly on chromatin accessibility and RNAa competence."""

    provenance: Provenance
    data: Mapping[DataAxis, DataAvailability] = field(default_factory=dict)
    """What evidence exists, per axis. A missing axis is treated as unavailable with no route."""

    def availability(self, axis: DataAxis) -> DataAvailability:
        """Return what is known along one axis, defaulting to "nothing, and no route"."""
        known = self.data.get(axis)
        if known is not None:
            return known
        return DataAvailability(
            axis=axis,
            available=False,
            note=f"No {axis} data is declared for {self.name}.",
        )

    @property
    def available_axes(self) -> tuple[DataAxis, ...]:
        """Return the axes with data loaded, in declaration order."""
        return tuple(axis for axis in DataAxis if self.availability(axis).available)

    def __str__(self) -> str:
        """Return the display name."""
        return self.name

    def to_dict(self) -> dict[str, object]:
        """Return the cell type as JSON-serialisable data."""
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "proliferating": self.is_proliferating,
            "provenance": str(self.provenance),
            "data": [self.availability(axis).to_dict() for axis in DataAxis],
        }


DATA_PACKAGE = "sarna.regulatory.data"
CELL_LINES_FILENAME = "cell_lines.toml"
CELL_LINES_ENVIRONMENT_VARIABLE = "SARNA_CELL_LINES"


def cell_line_paths(extra: Sequence[Path] | None = None) -> tuple[Path, ...]:
    """Return extra TOML files of cell lines to load, from arguments and the environment."""
    paths = [Path(candidate) for candidate in extra or ()]
    override = os.environ.get(CELL_LINES_ENVIRONMENT_VARIABLE)
    if override:
        location = Path(override)
        if location.is_dir():
            paths.extend(sorted(location.glob("*.toml")))
        elif location.exists():
            paths.append(location)
    return tuple(paths)


def load_cell_types(extra: Sequence[Path] | None = None) -> dict[str, CellType]:
    """Load the shipped cell lines plus any the user has added.

    Raises:
        ValueError: If a file is malformed or two lines share an id.
    """
    source = resources.files(DATA_PACKAGE).joinpath(CELL_LINES_FILENAME)
    with resources.as_file(source) as path:
        loaded = list(_read_cell_lines(path))
    for path in cell_line_paths(extra):
        loaded.extend(_read_cell_lines(path))
    registry: dict[str, CellType] = {}
    for cell in loaded:
        if cell.id in registry:
            raise ValueError(
                f"duplicate cell line id {cell.id!r}; ids appear in exports and must be unique"
            )
        registry[cell.id] = cell
    return registry


def _read_cell_lines(path: Path) -> Iterator[CellType]:
    """Read one TOML file of cell lines.

    Raises:
        ValueError: If the file cannot be parsed or an entry is missing a required field.
    """
    try:
        with path.open("rb") as handle:
            payload = tomllib.load(handle)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"{path} is not valid TOML: {error}") from error
    for index, entry in enumerate(payload.get("cell_line", []), start=1):
        try:
            yield _cell_type_from(entry)
        except (KeyError, ValueError) as error:
            raise ValueError(f"{path}, cell_line {index}: {error}") from error


def _cell_type_from(entry: dict[str, Any]) -> CellType:
    """Build one cell type from a parsed TOML table."""
    data = {}
    for axis in DataAxis:
        declared = entry.get(axis.toml_key)
        if declared is None:
            continue
        data[axis] = DataAvailability(
            axis=axis,
            available=bool(declared["available"]),
            note=str(declared.get("note", "")).strip(),
            source=declared.get("source"),
        )
    return CellType(
        id=entry["id"],
        name=entry["name"],
        description=str(entry["description"]).strip(),
        is_proliferating=bool(entry["proliferating"]),
        provenance=Provenance.citation(entry["source"], entry["locator"]),
        data=data,
    )


CELL_TYPES: dict[str, CellType] = load_cell_types()
"""The cell types this project has data or published characterisation for."""

HELA = CELL_TYPES["hela"]
AF22 = CELL_TYPES["af22"]


def cell_type(cell_type_id: str) -> CellType:
    """Return a known cell type by id.

    Raises:
        KeyError: If it is not one this project characterises. The message lists those that are.
    """
    try:
        return CELL_TYPES[cell_type_id]
    except KeyError:
        raise KeyError(
            f"unknown cell type {cell_type_id!r}; this project characterises "
            f"{', '.join(sorted(CELL_TYPES))}. Add one by pointing "
            f"{CELL_LINES_ENVIRONMENT_VARIABLE} at a TOML file, so its properties and provenance "
            f"travel with it."
        ) from None


@dataclass(frozen=True, slots=True)
class Evidence:
    """The evidential standing of one annotation."""

    kind: EvidenceKind
    provenance: Provenance
    cell_type: CellType | None = None
    """The cell type observed, for a cell-type observation. Must be ``None`` otherwise."""

    status: CellTypeStatus = CellTypeStatus.NOT_ASSAYED
    """What was found in that cell type. Defaults to not-assayed, which is usually the truth."""

    def __post_init__(self) -> None:
        """Enforce that the evidence kind and the cell-type claim agree.

        Raises:
            ValueError: If a motif prediction or aggregated annotation names a cell type, or if a
                cell-type observation does not.
        """
        if self.kind is EvidenceKind.CELL_TYPE_OBSERVATION:
            if self.cell_type is None:
                raise ValueError(
                    "a cell-type observation must name the cell type it was observed in; "
                    "without it the observation cannot be interpreted or reproduced"
                )
            if not self.status.is_known:
                raise ValueError(
                    f"a cell-type observation in {self.cell_type} must record what was found; "
                    f"'not-assayed' contradicts having observed it"
                )
            return
        if self.cell_type is not None:
            raise ValueError(
                f"evidence of kind {self.kind} is computed without reference to any cell type, "
                f"so it must not name one ({self.cell_type} given). Attaching a cell type here "
                f"would imply a measurement that was never made."
            )

    @property
    def is_cell_type_specific(self) -> bool:
        """Return whether this evidence speaks to a particular cell type."""
        return self.kind.is_cell_type_specific

    @property
    def summary(self) -> str:
        """Return a one-line statement of what is and is not claimed."""
        if self.kind is EvidenceKind.CELL_TYPE_OBSERVATION:
            assert self.cell_type is not None
            verb = (
                "observed bound"
                if self.status is CellTypeStatus.OBSERVED_BOUND
                else "assayed and not found bound"
            )
            return f"{verb} in {self.cell_type} ({self.provenance})"
        if self.kind is EvidenceKind.MOTIF_PREDICTION:
            return f"sequence matches the motif; no cell type implied ({self.provenance})"
        return f"annotated as regulatory across pooled cell types ({self.provenance})"


def motif_prediction(provenance: Provenance) -> Evidence:
    """Return evidence for a sequence-only motif match."""
    return Evidence(kind=EvidenceKind.MOTIF_PREDICTION, provenance=provenance)


def aggregated_annotation(provenance: Provenance) -> Evidence:
    """Return evidence for an annotation pooled across cell types."""
    return Evidence(kind=EvidenceKind.AGGREGATED_ANNOTATION, provenance=provenance)


def observed_in(
    cell: CellType,
    provenance: Provenance,
    status: CellTypeStatus = CellTypeStatus.OBSERVED_BOUND,
) -> Evidence:
    """Return evidence for an experimental observation in one cell type."""
    return Evidence(
        kind=EvidenceKind.CELL_TYPE_OBSERVATION,
        provenance=provenance,
        cell_type=cell,
        status=status,
    )


__all__ = [
    "AF22",
    "CELL_LINES_ENVIRONMENT_VARIABLE",
    "CELL_TYPES",
    "HELA",
    "CellType",
    "CellTypeStatus",
    "DataAvailability",
    "DataAxis",
    "Evidence",
    "EvidenceKind",
    "aggregated_annotation",
    "cell_line_paths",
    "cell_type",
    "load_cell_types",
    "motif_prediction",
    "observed_in",
]
