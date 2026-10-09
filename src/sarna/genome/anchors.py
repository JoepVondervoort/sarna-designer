"""Promoter anchors: what a promoter window is measured from, when it is not a transcript's TSS.

The glossary says a TSS belongs to a transcript and a promoter window is meaningless without one.
Chapter 2 of the thesis breaks that, for a good reason. Its SCN1A duplexes target promoters mapped
by **Dong et al. 2014** and **Long et al. 2008**, and the Dong distal region (which produced two of
the project's three SCN1A activators) has no Ensembl transcript within 5.5 kb. A tool that can only
measure from an annotated TSS could not have produced those two designs.

So the thing a window is measured from is generalised here into a **promoter anchor**: a named,
cited region with an optional reference point. A transcript's TSS is one kind of anchor and no
longer the only kind.

**The reference point carries its own evidence tier**, for the same reason the regulatory layer
does. "-500 bp" means something different when the zero is a GENCODE-annotated TSS, a TSS mapped in
a publication, or a point inferred from where somebody's designs happened to land.
:class:`ReferenceKind` makes a caller say which, and an anchor with no reference point produces a
promoter with no TSS at all, so the TSS rules report not-applicable rather than measuring from an
invented origin.

Anchors are **data, not code**: shipped ones live in ``data/promoter_anchors.toml`` with their
citations, and a user adds more by pointing ``SARNA_ANCHORS`` at a file or directory. A promoter
definition should arrive as a reviewable diff citing a paper, which is what this buys.
"""

import os
import tomllib
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from enum import Enum
from importlib import resources
from pathlib import Path
from typing import Any

from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.tss import TssPosition
from sarna.rules.provenance import Provenance

DATA_PACKAGE = "sarna.genome.data"
ANCHORS_FILENAME = "promoter_anchors.toml"
ANCHORS_ENVIRONMENT_VARIABLE = "SARNA_ANCHORS"

DEFAULT_UPSTREAM_BP = 2000
DEFAULT_DOWNSTREAM_BP = 500
DEFAULT_REGION_PADDING_BP = 500
"""Padding either side of a declared region when there is no reference point to measure from."""


class ReferenceKind(Enum):
    """How well established an anchor's zero point is."""

    ANNOTATED_TSS = "annotated-tss"
    """A TSS from the genome annotation. The strongest case, and the default path."""

    MAPPED_TSS = "mapped-tss"
    """A start site mapped experimentally in the cited publication."""

    INFERRED = "inferred"
    """Deduced rather than reported; the anchor's note must say how, and how uncertain."""

    NONE = "none"
    """The source describes a region but no start site. Offsets are not computable."""

    def __str__(self) -> str:
        """Return the kind label."""
        return self.value

    @property
    def is_measurable(self) -> bool:
        """Return whether offsets from this anchor mean anything."""
        return self is not ReferenceKind.NONE

    @property
    def caveat(self) -> str | None:
        """Return the warning a UI must show beside offsets from this anchor."""
        if self is ReferenceKind.ANNOTATED_TSS:
            return None
        if self is ReferenceKind.MAPPED_TSS:
            return (
                "Offsets are measured from a start site mapped in the cited publication, not from "
                "the genome annotation. They will not match a browser's TSS."
            )
        if self is ReferenceKind.INFERRED:
            return (
                "The zero point is inferred, not reported. Every offset from this anchor carries "
                "that uncertainty; see the anchor's note for how it was derived."
            )
        return (
            "This region has no start site, so no offset can be computed from it. Rules about "
            "position relative to the TSS will report not-applicable."
        )


@dataclass(frozen=True, slots=True)
class AnchorWindow:
    """A region to fetch, and what the first base's TSS position is if that is knowable."""

    interval: GenomicInterval
    tss_start: TssPosition | None
    """TSS position of the window's first base, or ``None`` for an unmeasurable anchor."""


@dataclass(frozen=True, slots=True)
class PromoterAnchor:
    """A named promoter region to design against, with where it comes from."""

    id: str
    gene: str
    interval: GenomicInterval
    """The region the source describes. Strand is the gene's own."""

    reference_kind: ReferenceKind
    provenance: Provenance
    note: str
    """Why this region, what the source says, and how a reference point was arrived at."""

    reference_point: int | None = None
    """0-based genomic coordinate of position ``+1``, or ``None`` if the source maps none."""

    def __post_init__(self) -> None:
        """Enforce that the reference point and its declared kind agree.

        Raises:
            ValueError: If a measurable kind has no reference point or vice versa, or if the
                reference point lies outside the region it belongs to.
        """
        if self.reference_kind.is_measurable and self.reference_point is None:
            raise ValueError(
                f"anchor {self.id!r} declares reference kind {self.reference_kind} but gives no "
                f"reference point; use 'none' if the source maps no start site"
            )
        if not self.reference_kind.is_measurable and self.reference_point is not None:
            raise ValueError(
                f"anchor {self.id!r} gives a reference point but declares kind 'none'; state "
                f"which kind it is so a reader knows how far to trust an offset"
            )
        if self.reference_point is not None and not (
            self.interval.start <= self.reference_point < self.interval.end
        ):
            raise ValueError(
                f"anchor {self.id!r} puts its reference point at {self.reference_point}, outside "
                f"its own region {self.interval.start}-{self.interval.end}. One of the two is "
                f"wrong, and offsets computed from them would be meaningless either way."
            )
        if not self.note.strip():
            raise ValueError(
                f"anchor {self.id!r} must explain itself; an uncited region with no note is "
                f"indistinguishable from a guess"
            )

    @property
    def strand(self) -> Strand:
        """Return the strand the gene is transcribed from."""
        return self.interval.strand

    @property
    def is_measurable(self) -> bool:
        """Return whether offsets from this anchor can be computed."""
        return self.reference_kind.is_measurable

    @classmethod
    def from_transcript(
        cls, gene: str, transcript_id: str, tss: int, interval: GenomicInterval
    ) -> "PromoterAnchor":
        """Build an anchor from an annotated transcript, the ordinary case.

        Args:
            gene: The gene symbol.
            transcript_id: The transcript whose TSS this is.
            tss: 0-based genomic coordinate of the TSS.
            interval: The transcript's span, which supplies contig, assembly and strand.
        """
        return cls(
            id=f"{gene.lower()}-{transcript_id.lower()}",
            gene=gene,
            interval=interval,
            reference_point=tss,
            reference_kind=ReferenceKind.ANNOTATED_TSS,
            provenance=Provenance.citation(
                "Ensembl annotation", f"transcript {transcript_id}, TSS at {tss}"
            ),
            note=f"Transcription start site of {transcript_id}, from the genome annotation.",
        )

    def window(
        self,
        upstream_bp: int = DEFAULT_UPSTREAM_BP,
        downstream_bp: int = DEFAULT_DOWNSTREAM_BP,
        region_padding_bp: int = DEFAULT_REGION_PADDING_BP,
    ) -> AnchorWindow:
        """Return the region to fetch and the TSS position of its first base.

        With a reference point, the window is ``upstream_bp`` and ``downstream_bp`` around it, and
        the first base's TSS position follows. Without one, the window is the declared region plus
        ``region_padding_bp`` either side and there is no TSS position, which is the honest answer,
        not a degraded one.

        Raises:
            ValueError: If a width is negative, if a measurable window would be empty, or if the
                window runs off the start of the contig.
        """
        if min(upstream_bp, downstream_bp, region_padding_bp) < 0:
            raise ValueError(
                "window widths are magnitudes and must not be negative; got "
                f"upstream={upstream_bp}, downstream={downstream_bp}, "
                f"padding={region_padding_bp}"
            )
        if self.reference_point is None:
            start = self.interval.start - region_padding_bp
            end = self.interval.end + region_padding_bp
            tss_start = None
        else:
            if upstream_bp + downstream_bp == 0:
                raise ValueError("a promoter window of zero length contains no candidates")
            if self.strand.sign > 0:
                start = self.reference_point - upstream_bp
                end = self.reference_point + downstream_bp
            else:
                # Upstream is toward higher coordinates on the minus strand, mirroring
                # Promoter.genomic_interval and EnsemblClient.promoter_for.
                start = self.reference_point - downstream_bp + 1
                end = self.reference_point + upstream_bp + 1
            tss_start = TssPosition.from_offset(-upstream_bp)
        if start < 0:
            raise ValueError(
                f"anchor {self.id!r} would start at {start}, before the beginning of contig "
                f"{self.interval.contig}; reduce the upstream width"
            )
        return AnchorWindow(
            interval=GenomicInterval(
                assembly=self.interval.assembly,
                contig=self.interval.contig,
                start=start,
                end=end,
                strand=self.strand,
            ),
            tss_start=tss_start,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the anchor as JSON-serialisable data, caveat included."""
        return {
            "id": self.id,
            "gene": self.gene,
            "assembly": str(self.interval.assembly),
            "contig": self.interval.contig,
            "start": self.interval.start,
            "end": self.interval.end,
            "strand": str(self.strand),
            "reference_point": self.reference_point,
            "reference_kind": str(self.reference_kind),
            "measurable": self.is_measurable,
            "provenance": str(self.provenance),
            "note": self.note,
            "caveat": self.reference_kind.caveat,
        }


class AnchorRegistry:
    """The promoter anchors available, keyed by id."""

    __slots__ = ("_by_id",)

    def __init__(self, anchors: Sequence[PromoterAnchor]) -> None:
        """Build a registry.

        Raises:
            ValueError: If two anchors share an id, which would make a run's record ambiguous.
        """
        by_id: dict[str, PromoterAnchor] = {}
        for anchor in anchors:
            if anchor.id in by_id:
                raise ValueError(
                    f"duplicate anchor id {anchor.id!r}; ids appear in exports and must be unique"
                )
            by_id[anchor.id] = anchor
        self._by_id = by_id

    def __len__(self) -> int:
        """Return how many anchors are held."""
        return len(self._by_id)

    def __iter__(self) -> Iterator[PromoterAnchor]:
        """Iterate the anchors in load order."""
        return iter(self._by_id.values())

    def __contains__(self, anchor_id: object) -> bool:
        """Return whether an id is present."""
        return anchor_id in self._by_id

    def __getitem__(self, anchor_id: str) -> PromoterAnchor:
        """Return the anchor with ``anchor_id``.

        Raises:
            KeyError: If no such anchor exists. The message lists what the gene does have.
        """
        try:
            return self._by_id[anchor_id]
        except KeyError:
            raise KeyError(
                f"no promoter anchor {anchor_id!r}; known anchors are "
                f"{', '.join(sorted(self._by_id)) or 'none'}"
            ) from None

    def for_gene(self, gene: str) -> tuple[PromoterAnchor, ...]:
        """Return the anchors declared for one gene, in load order.

        Empty is a normal answer: most genes have no literature-mapped promoter on file, and the
        caller falls back to the annotated transcript path.
        """
        wanted = gene.upper()
        return tuple(a for a in self if a.gene.upper() == wanted)

    @property
    def genes(self) -> tuple[str, ...]:
        """Return the genes with declared anchors, alphabetically."""
        return tuple(sorted({a.gene for a in self}))


def anchor_paths(extra: Sequence[Path] | None = None) -> tuple[Path, ...]:
    """Return the TOML files to load, shipped last so user files cannot be shadowed by them.

    ``SARNA_ANCHORS`` may name a file or a directory; a directory contributes its ``*.toml`` in
    sorted order.
    """
    paths: list[Path] = []
    for candidate in extra or ():
        paths.append(Path(candidate))
    override = os.environ.get(ANCHORS_ENVIRONMENT_VARIABLE)
    if override:
        location = Path(override)
        if location.is_dir():
            paths.extend(sorted(location.glob("*.toml")))
        elif location.exists():
            paths.append(location)
    return tuple(paths)


def load_anchors(extra: Sequence[Path] | None = None) -> AnchorRegistry:
    """Load the shipped anchors plus any the user has added.

    Args:
        extra: Further TOML files to load, beyond the shipped set and ``SARNA_ANCHORS``.

    Returns:
        The registry.

    Raises:
        ValueError: If a file is malformed, an anchor is internally inconsistent, or two anchors
            share an id.
    """
    source = resources.files(DATA_PACKAGE).joinpath(ANCHORS_FILENAME)
    with resources.as_file(source) as path:
        anchors = list(_read(path))
    for path in anchor_paths(extra):
        anchors.extend(_read(path))
    return AnchorRegistry(anchors)


def _read(path: Path) -> Iterator[PromoterAnchor]:
    """Read one TOML file of anchors.

    Raises:
        ValueError: If the file cannot be parsed or an entry is missing a required field.
    """
    try:
        with path.open("rb") as handle:
            payload = tomllib.load(handle)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"{path} is not valid TOML: {error}") from error
    for index, entry in enumerate(payload.get("anchor", []), start=1):
        try:
            yield _anchor_from(entry)
        except (KeyError, ValueError) as error:
            raise ValueError(f"{path}, anchor {index}: {error}") from error


def _anchor_from(entry: dict[str, Any]) -> PromoterAnchor:
    """Build one anchor from a parsed TOML table."""
    kind = ReferenceKind(entry.get("reference_kind", "none"))
    return PromoterAnchor(
        id=entry["id"],
        gene=entry["gene"],
        interval=GenomicInterval(
            assembly=Assembly(entry.get("assembly", "GRCh38")),
            contig=str(entry["contig"]),
            start=int(entry["start"]),
            end=int(entry["end"]),
            strand=Strand.PLUS if entry.get("strand", "+") == "+" else Strand.MINUS,
        ),
        reference_point=(
            int(entry["reference_point"]) if entry.get("reference_point") is not None else None
        ),
        reference_kind=kind,
        provenance=Provenance.citation(entry["source"], entry["locator"]),
        note=entry["note"],
    )


__all__ = [
    "ANCHORS_ENVIRONMENT_VARIABLE",
    "DEFAULT_DOWNSTREAM_BP",
    "DEFAULT_REGION_PADDING_BP",
    "DEFAULT_UPSTREAM_BP",
    "AnchorRegistry",
    "AnchorWindow",
    "PromoterAnchor",
    "ReferenceKind",
    "anchor_paths",
    "load_anchors",
]
