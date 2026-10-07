"""AF22 chromatin state, read through the liftover boundary ADR 0004 put it behind.

This is the only cell-type-specific evidence the project has, and it answers the question chapter 4
left open. That chapter's largest effect was not a sequence effect: **no saRNA activated any locus
in AF22**, including the P21 control that works reliably in HeLa, and the stated hypothesis was
more restrictive chromatin. Six histone marks for AF22 are available (GSE108506), so that
hypothesis is testable per locus rather than argued in general.

**Everything here crosses an assembly boundary, which is defect B15.** The tracks are
GRCh37/hg19; the engine is GRCh38. Measured on SCN1A's MANE Select TSS: GRCh38 chr2:166,128,019 is
hg19 chr2:166,984,529 — **856,510 bp apart**. Reading a GRCh38 coordinate straight out of these
files does not return noisy data, it returns a different gene's data, and the legacy pipeline did
exactly that. So conversion happens here, explicitly, and four ways it can go wrong are checked
rather than assumed:

- **The interval may not map at all.** Some regions are absent from the other assembly.
- **It may split.** The two ends can land on non-contiguous coordinates across an indel, in which
  case the span between them is not the same span.
- **It may change contig or strand.** Rare, and silently catastrophic.
- **It may be the mitochondrion**, where hg19 (16,571 bp, NC_001807) and GRCh37 (16,569 bp, rCRS)
  are *different sequences* that no liftover reconciles. Refused outright.

Each failure returns a stated reason, never a number. A signal value that might be from the wrong
locus is worse than no signal value, because it looks like data.
"""

import os
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Protocol

from sarna.genome.coordinates import Assembly, GenomicInterval
from sarna.regulatory.evidence import AF22, CellType, CellTypeStatus, Evidence, observed_in
from sarna.rules.provenance import Provenance

DEFAULT_TRACK_DIR = Path.home() / "af22_data"
TRACK_DIR_ENVIRONMENT_VARIABLE = "SARNA_AF22_DIR"

GSE108506 = Provenance.citation(
    "GSE108506 (AF22 lt-NES ChIP-seq)",
    "Falk et al. AF22 neuroepithelial stem cells, 50 bp bins, aligned to GRCh37/hg19",
)

MITOCHONDRIAL_CONTIGS = frozenset({"MT", "M", "chrM", "chrMT"})
"""Refused outright: hg19 and GRCh37 hold different mitochondrial sequences (inventory §12)."""

BACKGROUND_WINDOWS = 250
"""Random same-contig windows sampled to establish what "no signal" looks like.

Needed because these are **read-depth tracks with no input or IgG normalisation**. Every mark
carries counts everywhere, so asking whether a mark is "present" at a locus is meaningless without
a baseline: the SCN1A promoter reads H3K4me3 around 1.0, and so does random intergenic DNA.
"""

BACKGROUND_SEED = 11
"""Fixed, so a reported fold-enrichment is reproducible."""

ENRICHMENT_THRESHOLD = 5.0
"""Fold over the background median at which a mark is called present.

A convention, and nothing examined sits near it. Measured in AF22 at +/-500 bp of the MANE Select
TSS: GAPDH 87.8x, ACTB 144.7x, TSC1 148.1x, DEPDC5 159.7x, CDKN1A 41.6x, CHD2 26.6x, NPRL3 17.0x
against SCN1A at 1.15x. The housekeeping promoters are the positive control that shows the
measurement works at all; see docs/scientific/af22-chromatin.md.
"""

MAXIMUM_LENGTH_DRIFT_FRACTION = 0.02
"""How much an interval's length may change under liftover before it is treated as unmappable.

A few bases of drift across a small indel still describes the same locus. A span that changes by
more than this is a different span, and averaging signal across it answers a different question
from the one asked.
"""


class HistoneMark(Enum):
    """A histone modification, and what its presence at a promoter means."""

    H3K4ME3 = "H3K4me3"
    H3K27AC = "H3K27ac"
    H3K4ME1 = "H3K4me1"
    H3K27ME3 = "H3K27me3"
    H3K9ME3 = "H3K9me3"
    H3K36ME3 = "H3K36me3"

    def __str__(self) -> str:
        """Return the mark's conventional name."""
        return self.value

    @property
    def interpretation(self) -> str:
        """Return what this mark indicates at a promoter, in one sentence."""
        return {
            HistoneMark.H3K4ME3: (
                "Active or poised promoter. The single most informative mark for whether a "
                "promoter is available to be activated at all."
            ),
            HistoneMark.H3K27AC: "Active promoter or enhancer; distinguishes active from poised.",
            HistoneMark.H3K4ME1: (
                "Primed or enhancer-like chromatin rather than an active promoter."
            ),
            HistoneMark.H3K27ME3: (
                "Polycomb repression. Present at a promoter, it is the most likely chromatin "
                "explanation for a saRNA failing to activate a gene that is otherwise expressed."
            ),
            HistoneMark.H3K9ME3: (
                "Constitutive heterochromatin; a promoter here is effectively shut."
            ),
            HistoneMark.H3K36ME3: (
                "Active transcription through a gene body, not a promoter mark."
            ),
        }[self]

    @property
    def is_activating(self) -> bool:
        """Return whether the mark indicates an accessible, active promoter."""
        return self in (HistoneMark.H3K4ME3, HistoneMark.H3K27AC)

    @property
    def is_repressive(self) -> bool:
        """Return whether the mark indicates a repressed locus."""
        return self in (HistoneMark.H3K27ME3, HistoneMark.H3K9ME3)


class LiftoverProblem(Enum):
    """Why a GRCh38 interval could not be read from a GRCh37 track."""

    MITOCHONDRIAL = "mitochondrial"
    UNMAPPED = "unmapped"
    SPLIT = "split-across-an-indel"
    CONTIG_CHANGED = "contig-changed"
    STRAND_FLIPPED = "strand-flipped"

    def __str__(self) -> str:
        """Return the problem label."""
        return self.value


class LiftoverError(ValueError):
    """Raised when an interval cannot be converted between assemblies.

    Carries the problem so a caller can report it rather than treating every failure alike.
    """

    def __init__(self, problem: LiftoverProblem, detail: str) -> None:
        """Record which way the conversion failed, and why it matters."""
        self.problem = problem
        super().__init__(f"{problem}: {detail}")


class CoordinateConverter(Protocol):
    """The slice of ``pyliftover.LiftOver`` this module uses, so it can be stubbed in tests."""

    def convert_coordinate(
        self, chromosome: str, position: int
    ) -> Sequence[tuple[str, int, str, int]] | None:
        """Return candidate conversions of one coordinate, or ``None`` for an unknown contig."""
        ...


class SignalReader(Protocol):
    """The slice of ``pyBigWig`` this module uses."""

    def chroms(self) -> Mapping[str, int]:
        """Return contig lengths as the file declares them."""
        ...

    def stats(
        self, chromosome: str, start: int, end: int, type: str = ..., nBins: int = ...
    ) -> Sequence[float | None]:
        """Return a summary statistic over a span."""
        ...


@dataclass(frozen=True, slots=True)
class BackgroundModel:
    """What one mark reads over random windows of a given width on one contig.

    The baseline that makes a signal value mean anything. Without it, "H3K4me3 is 1.0 here" is
    unreadable; against it, 1.0 is background and 88 is a promoter.
    """

    contig: str
    mark: HistoneMark
    window_length: int
    median: float
    percentile_95: float
    windows: int
    seed: int

    def fold(self, value: float) -> float:
        """Return ``value`` as a multiple of the background median."""
        return value / self.median if self.median > 0 else float("inf")


@dataclass(frozen=True, slots=True)
class ChromatinSignal:
    """One mark's signal over one interval, in one cell type, against its background."""

    mark: HistoneMark
    cell_type: CellType
    mean: float
    maximum: float
    background: BackgroundModel
    source_interval: GenomicInterval
    """The interval actually read, on the **track's** assembly, after conversion."""

    provenance: Provenance
    """Where this track came from. Carried per signal rather than read from a module constant,
    so a reader over another cell line's tracks cannot cite AF22's accession."""

    @property
    def fold_over_background(self) -> float:
        """Return the mean signal as a multiple of random same-contig windows."""
        return self.background.fold(self.mean)

    @property
    def is_present(self) -> bool:
        """Return whether the mark is enriched above background at this locus.

        Enrichment, not raw signal: see :data:`BACKGROUND_WINDOWS` for why a raw value cannot
        answer this, and :data:`ENRICHMENT_THRESHOLD` for where the line sits and what validates it.
        """
        return self.fold_over_background >= ENRICHMENT_THRESHOLD

    @property
    def evidence(self) -> Evidence:
        """Return this signal's evidential standing: measured, in a named cell type.

        ``OBSERVED_BOUND`` reads oddly for a histone mark and is literally right: modified
        nucleosomes are bound at that locus. Absence is ``ASSAYED_NOT_BOUND`` rather than
        not-assayed, because the track was read and reported no enrichment.
        """
        return observed_in(
            self.cell_type,
            self.provenance,
            CellTypeStatus.OBSERVED_BOUND
            if self.is_present
            else CellTypeStatus.ASSAYED_NOT_BOUND,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the signal as JSON-serialisable data."""
        return {
            "mark": str(self.mark),
            "cell_type": self.cell_type.id,
            "mean": round(self.mean, 3),
            "maximum": round(self.maximum, 3),
            "background_median": round(self.background.median, 3),
            "fold_over_background": round(self.fold_over_background, 2),
            "present": self.is_present,
            "threshold": ENRICHMENT_THRESHOLD,
            "interpretation": self.mark.interpretation,
            "read_from": {
                "assembly": str(self.source_interval.assembly),
                "contig": self.source_interval.contig,
                "start": self.source_interval.start,
                "end": self.source_interval.end,
            },
        }


@dataclass(frozen=True, slots=True)
class ChromatinProfile:
    """Every available mark over one interval, with a verdict about promoter accessibility."""

    interval: GenomicInterval
    """The interval as asked for, on the engine's assembly."""

    cell_type: CellType
    signals: tuple[ChromatinSignal, ...]
    provenance: Provenance
    """Where these tracks came from, stated by the reader that produced the profile."""

    problem: LiftoverProblem | None = None
    problem_detail: str = ""

    @property
    def is_readable(self) -> bool:
        """Return whether the interval could be read at all."""
        return self.problem is None

    def signal(self, mark: HistoneMark) -> ChromatinSignal | None:
        """Return one mark's signal, or ``None`` if that track was not available."""
        for signal in self.signals:
            if signal.mark is mark:
                return signal
        return None

    @property
    def activating_marks(self) -> tuple[HistoneMark, ...]:
        """Return the activating marks present here."""
        return tuple(s.mark for s in self.signals if s.mark.is_activating and s.is_present)

    @property
    def repressive_marks(self) -> tuple[HistoneMark, ...]:
        """Return the repressive marks present here."""
        return tuple(s.mark for s in self.signals if s.mark.is_repressive and s.is_present)

    @property
    def summary(self) -> str:
        """Return a plain-language reading, with its limits stated.

        Deliberately does not produce a score. These are read-depth tracks with no peak calls and
        no input normalisation available here, so the honest output is which marks carry signal and
        what that suggests, not a number that invites comparison between loci.
        """
        if not self.is_readable:
            return (
                f"Chromatin state in {self.cell_type} could not be read: "
                f"{self.problem} — {self.problem_detail}"
            )
        if not self.signals:
            return f"No {self.cell_type} tracks were available to read."
        active = ", ".join(
            f"{s.mark} ({s.fold_over_background:.0f}x)"
            for s in self.signals
            if s.mark.is_activating and s.is_present
        ) or "none"
        repressed = ", ".join(
            f"{s.mark} ({s.fold_over_background:.0f}x)"
            for s in self.signals
            if s.mark.is_repressive and s.is_present
        ) or "none"
        reading = (
            f"Activating marks present: {active}. Repressive marks present: {repressed}."
        )
        if self.activating_marks and not self.repressive_marks:
            return f"{reading} Consistent with an accessible promoter in {self.cell_type}."
        if self.repressive_marks and not self.activating_marks:
            return (
                f"{reading} Consistent with a repressed promoter in {self.cell_type}, which "
                f"would explain a saRNA failing to activate here."
            )
        if self.activating_marks and self.repressive_marks:
            return f"{reading} Mixed, which at a promoter can indicate a bivalent state."
        return (
            f"{reading} No mark is enriched above background, so this promoter looks closed in "
            f"{self.cell_type}. By the mechanism RNAa is understood to use, a saRNA cannot "
            f"activate a locus that is not already transcribed."
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the profile as JSON-serialisable data."""
        return {
            "cell_type": self.cell_type.id,
            "readable": self.is_readable,
            "problem": str(self.problem) if self.problem else None,
            "problem_detail": self.problem_detail or None,
            "summary": self.summary,
            "provenance": str(self.provenance),
            "assembly_note": (
                "Tracks are GRCh37/hg19 and queries are GRCh38; coordinates are converted per "
                "request and the conversion is checked for splits, contig changes and strand "
                "flips (defect B15, ADR 0004)."
            ),
            "signals": [signal.to_dict() for signal in self.signals],
        }


def convert_interval(
    interval: GenomicInterval, converter: CoordinateConverter, target: Assembly
) -> GenomicInterval:
    """Convert ``interval`` to ``target``, refusing anything the conversion cannot guarantee.

    Args:
        interval: The interval to convert, on the engine's assembly.
        converter: A liftover chain, GRCh38 to GRCh37.
        target: The assembly being converted to, used to tag the result.

    Returns:
        The converted interval, on ``target``, on the same strand.

    Raises:
        LiftoverError: If the interval is mitochondrial, does not map, splits across an indel,
            changes contig, or flips strand. Each carries its :class:`LiftoverProblem`.
    """
    if interval.contig in MITOCHONDRIAL_CONTIGS:
        raise LiftoverError(
            LiftoverProblem.MITOCHONDRIAL,
            "hg19 holds NC_001807 (16,571 bp) and GRCh37 holds rCRS (16,569 bp). These are "
            "different sequences, so no liftover reconciles them. Histone signal is not "
            "interpretable on mitochondrial DNA anyway, which is not nucleosome-packaged.",
        )
    ucsc = f"chr{interval.contig}"
    first = _single(converter, ucsc, interval.start)
    last = _single(converter, ucsc, interval.end - 1)
    if first[0] != last[0]:
        raise LiftoverError(
            LiftoverProblem.CONTIG_CHANGED,
            f"the interval's ends map to {first[0]} and {last[0]}, so it spans a rearrangement",
        )
    if first[2] != last[2]:
        raise LiftoverError(
            LiftoverProblem.STRAND_FLIPPED,
            f"the interval's ends map to opposite strands ({first[2]} and {last[2]})",
        )
    low, high = min(first[1], last[1]), max(first[1], last[1]) + 1
    converted_length = high - low
    drift = abs(converted_length - interval.length) / interval.length
    if drift > MAXIMUM_LENGTH_DRIFT_FRACTION:
        raise LiftoverError(
            LiftoverProblem.SPLIT,
            f"the interval is {interval.length} bp on {interval.assembly} but "
            f"{converted_length} bp after conversion, a {drift:.0%} change. The two spans do not "
            f"describe the same locus, so signal averaged over the second would answer a "
            f"different question",
        )
    return GenomicInterval(
        assembly=target,
        contig=first[0].removeprefix("chr"),
        start=low,
        end=high,
        strand=interval.strand,
    )


def _single(
    converter: CoordinateConverter, contig: str, position: int
) -> tuple[str, int, str, int]:
    """Return the one unambiguous conversion of a coordinate.

    Raises:
        LiftoverError: If the position does not map, or maps to more than one place. An ambiguous
            position is not a usable one: picking the first candidate is how a pipeline ends up
            reading a paralogue.
    """
    candidates = converter.convert_coordinate(contig, position) or []
    if not candidates:
        raise LiftoverError(
            LiftoverProblem.UNMAPPED,
            f"{contig}:{position} has no equivalent in the target assembly",
        )
    if len(candidates) > 1:
        places = ", ".join(f"{c[0]}:{c[1]}" for c in candidates[:4])
        raise LiftoverError(
            LiftoverProblem.UNMAPPED,
            f"{contig}:{position} maps to {len(candidates)} places ({places}), so no single "
            f"locus can be read",
        )
    return tuple(candidates[0])  # type: ignore[return-value]


class Af22ChromatinReader:
    """Reads one cell line's histone signal for GRCh38 intervals, converting at the boundary.

    The cell line and its source are constructor state, not per-call arguments: a set of tracks
    belongs to one cell line, and letting a caller name a different one at the point of use is
    how a HeLa profile came to cite AF22's accession.
    """

    __slots__ = (
        "_backgrounds",
        "_converter",
        "_readers",
        "cell_type",
        "provenance",
        "track_assembly",
    )

    def __init__(
        self,
        readers: Mapping[HistoneMark, SignalReader],
        converter: CoordinateConverter,
        *,
        cell_type: CellType,
        provenance: Provenance,
        track_assembly: Assembly = Assembly.GRCH37,
    ) -> None:
        """Build a reader over already-opened tracks.

        Args:
            readers: One signal reader per available mark.
            converter: A liftover chain from the engine's assembly to ``track_assembly``.
            cell_type: The cell line these tracks were measured in. Required, because every
                signal this reader produces is labelled with it.
            provenance: Where the tracks came from. Required for the same reason: it is cited
                on every signal and in every export.
            track_assembly: The assembly the tracks are on. Stated, never inferred.
        """
        self._readers = dict(readers)
        self._converter = converter
        self.cell_type = cell_type
        self.provenance = provenance
        self.track_assembly = track_assembly
        self._backgrounds: dict[tuple[str, HistoneMark, int], BackgroundModel] = {}

    @property
    def marks(self) -> tuple[HistoneMark, ...]:
        """Return the marks this reader has tracks for, in enum order."""
        return tuple(mark for mark in HistoneMark if mark in self._readers)

    def profile(self, interval: GenomicInterval) -> ChromatinProfile:
        """Return every available mark's signal over ``interval``.

        The cell type is this reader's own, set when it was built, so a profile cannot be
        attributed to a cell line whose tracks it does not hold.

        Args:
            interval: The region to read, on the engine's assembly.

        Returns:
            The profile. If conversion failed, it carries the reason and no signals rather than
            signal read from an unverified locus.

        Raises:
            ValueError: If the interval is already on the tracks' assembly, which would mean the
                caller has mixed conventions somewhere upstream.
        """
        if interval.assembly is self.track_assembly:
            raise ValueError(
                f"interval is already on {self.track_assembly}; this reader converts *to* it, so "
                f"being handed one suggests an assembly was mixed up upstream"
            )
        try:
            converted = convert_interval(interval, self._converter, self.track_assembly)
        except LiftoverError as error:
            return ChromatinProfile(
                interval=interval,
                cell_type=self.cell_type,
                signals=(),
                provenance=self.provenance,
                problem=error.problem,
                problem_detail=str(error).split(": ", 1)[-1],
            )
        signals = []
        for mark, reader in self._readers.items():
            signal = self._read(reader, mark, converted)
            if signal is not None:
                signals.append(signal)
        return ChromatinProfile(
            interval=interval,
            cell_type=self.cell_type,
            signals=tuple(sorted(signals, key=lambda s: list(HistoneMark).index(s.mark))),
            provenance=self.provenance,
        )

    def _read(
        self,
        reader: SignalReader,
        mark: HistoneMark,
        interval: GenomicInterval,
    ) -> ChromatinSignal | None:
        """Read one mark over a converted interval, or ``None`` if the contig is absent."""
        contig = f"chr{interval.contig}"
        available = reader.chroms()
        if contig not in available:
            return None
        end = min(interval.end, available[contig])
        if end <= interval.start:
            return None
        mean = reader.stats(contig, interval.start, end, type="mean", nBins=1)[0] or 0.0
        maximum = reader.stats(contig, interval.start, end, type="max", nBins=1)[0] or 0.0
        return ChromatinSignal(
            mark=mark,
            cell_type=self.cell_type,
            provenance=self.provenance,
            mean=float(mean),
            maximum=float(maximum),
            background=self.background(reader, mark, contig, end - interval.start),
            source_interval=interval,
        )

    def signal_track(
        self, interval: GenomicInterval, mark: HistoneMark, bins: int = 200
    ) -> tuple[float, ...] | None:
        """Return ``mark``'s mean signal in ``bins`` equal slices across ``interval``.

        For drawing a profile rather than judging a locus, so it returns raw means and leaves
        scaling to the caller. Bins are in the interval's **own 5'-to-3' order**: on the minus
        strand the track is reversed after reading, because the promoter is drawn in transcript
        orientation and a profile read in reference order would be mirrored.

        Returns:
            One mean per bin, or ``None`` if the mark has no track or the interval cannot be
            converted. ``None`` rather than zeros: no data is not an absence of signal.
        """
        reader = self._readers.get(mark)
        if reader is None:
            return None
        try:
            converted = convert_interval(interval, self._converter, self.track_assembly)
        except LiftoverError:
            return None
        contig = f"chr{converted.contig}"
        available = reader.chroms()
        if contig not in available:
            return None
        end = min(converted.end, available[contig])
        if end <= converted.start:
            return None
        values = reader.stats(contig, converted.start, end, type="mean", nBins=bins)
        read = tuple(float(v or 0.0) for v in values)
        return tuple(reversed(read)) if interval.strand.sign < 0 else read

    def background(
        self, reader: SignalReader, mark: HistoneMark, contig: str, window_length: int
    ) -> BackgroundModel:
        """Return what ``mark`` reads over random ``window_length`` windows on ``contig``.

        Sampled once per contig, mark and width, then cached: a promoter view asks for the same
        three numbers for every candidate it draws, and resampling for each would make the page
        take minutes to build.

        Args:
            reader: The track to sample.
            mark: Which mark it holds, for the returned model's label.
            contig: Contig name as the **track** spells it, e.g. ``chr2``.
            window_length: Width of the windows to sample, matching the query's width.

        Returns:
            The background model, with its seed and sample count recorded.
        """
        key = (contig, mark, window_length)
        cached = self._backgrounds.get(key)
        if cached is not None:
            return cached
        length = reader.chroms()[contig]
        rng = random.Random(BACKGROUND_SEED)
        # Kept clear of the telomeres, where coverage is unrepresentative.
        margin = min(1_000_000, max(0, (length - window_length) // 4))
        highest = max(margin + 1, length - margin - window_length)
        samples = []
        for _ in range(BACKGROUND_WINDOWS):
            start = rng.randrange(margin, highest)
            value = reader.stats(
                contig, start, start + window_length, type="mean", nBins=1
            )[0]
            samples.append(float(value or 0.0))
        samples.sort()
        model = BackgroundModel(
            contig=contig,
            mark=mark,
            window_length=window_length,
            median=samples[len(samples) // 2],
            percentile_95=samples[min(len(samples) - 1, int(0.95 * len(samples)))],
            windows=BACKGROUND_WINDOWS,
            seed=BACKGROUND_SEED,
        )
        self._backgrounds[key] = model
        return model


def track_directory() -> Path:
    """Return where the AF22 tracks live, from the environment or the default location."""
    override = os.environ.get(TRACK_DIR_ENVIRONMENT_VARIABLE)
    return Path(override) if override else DEFAULT_TRACK_DIR


def find_tracks(directory: Path | None = None) -> dict[HistoneMark, Path]:
    """Return the bigWig file for each mark found in ``directory``.

    Matches on the mark's name appearing in the filename, which is how GSE108506 names its files
    (``GSM2902406_No0028_AF22_H3K4me3_bs50.bw``). Marks with no file are simply absent, so a
    partial download still yields a usable reader over the marks that are present.
    """
    folder = directory or track_directory()
    if not folder.is_dir():
        return {}
    found: dict[HistoneMark, Path] = {}
    for path in sorted(folder.glob("*.bw")) + sorted(folder.glob("*.bigWig")):
        for mark in HistoneMark:
            if mark.value.lower() in path.name.lower() and mark not in found:
                found[mark] = path
    return found


def open_af22_reader(directory: Path | None = None) -> Af22ChromatinReader | None:
    """Open a reader over whatever AF22 tracks are present, or return ``None`` if none are.

    Returns ``None`` rather than raising, because the tracks are 890 MB of data this project does
    not ship: their absence is the normal case and must not stop a design run. The caller reports
    chromatin state as unavailable, which the evidence model already distinguishes from absent.

    Raises:
        ImportError: If the ``genomics`` extra is not installed but tracks are present, which is a
            configuration mistake rather than a missing-data case.
    """
    paths = find_tracks(directory)
    if not paths:
        return None
    try:
        import pyBigWig
        from pyliftover import LiftOver
    except ImportError as error:  # pragma: no cover - depends on the install
        raise ImportError(
            f"AF22 tracks are present in {directory or track_directory()} but the genomics extra "
            f'is not installed. Run: pip install -e ".[genomics]"'
        ) from error
    readers = {mark: pyBigWig.open(str(path)) for mark, path in paths.items()}
    return Af22ChromatinReader(
        readers=readers,
        converter=LiftOver("hg38", "hg19"),
        cell_type=AF22,
        provenance=GSE108506,
    )


__all__ = [
    "BACKGROUND_WINDOWS",
    "ENRICHMENT_THRESHOLD",
    "GSE108506",
    "MAXIMUM_LENGTH_DRIFT_FRACTION",
    "Af22ChromatinReader",
    "BackgroundModel",
    "ChromatinProfile",
    "ChromatinSignal",
    "HistoneMark",
    "LiftoverError",
    "LiftoverProblem",
    "convert_interval",
    "find_tracks",
    "open_af22_reader",
    "track_directory",
]
