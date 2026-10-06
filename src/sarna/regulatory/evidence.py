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

from dataclasses import dataclass
from enum import Enum

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

    def __str__(self) -> str:
        """Return the display name."""
        return self.name


HELA = CellType(
    id="hela",
    name="HeLa",
    description=(
        "Cervical adenocarcinoma line, continuously proliferating, with globally open and "
        "transcriptionally active chromatin. The permissive comparator for RNAa: chapter 4's "
        "P21 positive control activates reliably here, up to 14-fold on some plates."
    ),
    is_proliferating=True,
    provenance=Provenance.citation(
        "van de Vondervoort 2026", "chapter 4, Material and methods — cell culture"
    ),
)

AF22 = CellType(
    id="af22",
    name="AF22 (lt-NES)",
    description=(
        "Long-term neuroepithelial-like stem cells derived from human iPSCs (Falk et al. 2012). "
        "Divide slowly, retain tighter epigenetic control over lineage-specific genes, and carry "
        "fewer active histone marks at many promoters. Chapter 4 found no saRNA activated any "
        "locus here, including the P21 control, which is the largest single effect in that "
        "screen and is a property of the cell rather than of any sequence."
    ),
    is_proliferating=False,
    provenance=Provenance.citation(
        "Falk et al. 2012", "PLOS ONE 7(1):e29597, capture of neuroepithelial-like stem cells"
    ),
)

CELL_TYPES: dict[str, CellType] = {cell.id: cell for cell in (HELA, AF22)}
"""The cell types this project has data or published characterisation for."""


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
            f"{', '.join(sorted(CELL_TYPES))}. Add one in sarna.regulatory.evidence rather than "
            f"passing a bare string, so its properties and provenance travel with it."
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
    "CELL_TYPES",
    "HELA",
    "CellType",
    "CellTypeStatus",
    "Evidence",
    "EvidenceKind",
    "aggregated_annotation",
    "cell_type",
    "motif_prediction",
    "observed_in",
]
