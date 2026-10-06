"""Transcription factor binding motifs, and honest scanning for them.

A position weight matrix finds where a factor *could* bind. The gap between that and where a
factor *does* bind is the whole difficulty, and it is large: a 9 bp motif like SP1's matches
roughly once every few hundred bases of random sequence at the conventional threshold, and promoter
sequence is far from random — it is GC-rich, which is exactly what SP1, NRF1, EGR1 and GABPA
recognise. Scanning a 2 kb CpG-island promoter for sixteen factors produces dozens of hits, and
most of them are not bound in any cell.

This module therefore does three things beyond scanning:

**It reports the expected false-positive rate alongside the hits.** :attr:`MotifScan.expected_hits`
is how many matches a shuffled sequence of the same length and base composition would produce at
the same threshold. A factor whose observed count is at or below that number has told you nothing.

**It refuses to claim a cell type.** Every site it produces carries
:data:`~sarna.regulatory.evidence.EvidenceKind.MOTIF_PREDICTION`, and the evidence model forbids
attaching a cell type to that kind.

**It keeps the matrices' provenance.** Each motif carries its JASPAR accession, the experiment type
behind it, and the PubMed id, because "SP1 binds here" is a different claim depending on whether
the matrix came from ChIP-seq or from SELEX.

Scores are log-odds against a uniform background, reported as a **relative score** in [0, 1]:
``(score - minimum) / (maximum - minimum)`` over the matrix's achievable range, which is the
convention JASPAR's own tools use. The default threshold of 0.80 is that convention's default. It
is a convention, not a statistic — a relative score has no p-value — which is why the shuffled
expectation is reported beside it.
"""

import json
import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from importlib import resources

from sarna.genome.coordinates import Strand
from sarna.regulatory.evidence import Evidence, motif_prediction
from sarna.rules.provenance import Provenance
from sarna.sequence.nucleic import NucleicSequence

DATA_PACKAGE = "sarna.regulatory.data"
MATRIX_FILENAME = "jaspar_core_selection.json"

BASES = "ACGT"

DEFAULT_RELATIVE_THRESHOLD = 0.80
"""JASPAR's conventional relative-score cutoff. A convention, not a significance level."""

MINIMUM_RAW_SCORE = 0.0
"""Hard floor on the raw log-odds, applied in addition to the relative threshold.

The relative score rescales a matrix's achievable range onto [0, 1], which makes thresholds
comparable between matrices but silently depends on how negative the *minimum* is. A matrix with
one near-zero count has a very negative minimum, so 0.80 of its range can still be a negative
log-odds — meaning the match is *less* likely under the motif than under background.

Measured over this project's sixteen curated matrices, 4 of them are in that state at a relative
score of 0.80: NRF1 (-10.2), MEF2C (-1.7), YY1 (-1.4) and ETS1 (-1.1). Before this floor, NRF1
matched 141 times in 408 bp of GC-rich sequence.

A log-odds above zero is the weakest claim worth drawing on a diagram: the sequence is at least
more consistent with the motif than with background. The rigorous alternative is a per-matrix
p-value from the exact score distribution, which this module does not implement; see
docs/scientific/tf-binding.md.
"""

PSEUDOCOUNT = 0.25
"""Added to every cell before taking logs, so a zero count does not give negative infinity.

0.25 is the usual square-root-of-N-scaled choice simplified to a flat prior of one observation
spread over four bases.
"""

SHUFFLES_FOR_EXPECTATION = 20
"""Shuffles used to estimate the false-positive rate. Enough for a number to compare against."""

JASPAR_LICENCE_NOTE = (
    "JASPAR matrices are distributed under CC BY-NC-SA 4.0, which is non-commercial. Bundling "
    "them constrains this package's own licensing; see docs/scientific/tf-binding.md."
)


class MotifDataType(Enum):
    """The experiment a matrix was derived from. Bears on how much it is worth."""

    CHIP_SEQ = "ChIP-seq"
    HT_SELEX = "HT-SELEX"
    SELEX = "SELEX"
    OTHER = "other"

    @classmethod
    def parse(cls, label: str | None) -> "MotifDataType":
        """Return the matching type, or :attr:`OTHER` for anything unrecognised."""
        for member in cls:
            if member.value.lower() == (label or "").lower():
                return member
        return cls.OTHER

    def __str__(self) -> str:
        """Return the type label."""
        return self.value

    @property
    def is_in_vivo(self) -> bool:
        """Return whether the matrix came from binding measured in cells.

        ChIP-seq matrices reflect where a factor was actually found, so they carry chromatin and
        cofactor context. SELEX matrices describe the protein's intrinsic sequence preference, in
        vitro, which is cleaner and less predictive of occupancy.
        """
        return self is MotifDataType.CHIP_SEQ


# Deliberately without ``slots``, unlike every other dataclass here: the derived log-odds matrix
# is consulted once per scanned window and must be cached, and ``cached_property`` needs a
# ``__dict__``. Recomputing it per window made scanning a 2 kb promoter about forty times slower.
@dataclass(frozen=True)
class Motif:
    """One position weight matrix, with where it came from and why it is in this project."""

    matrix_id: str
    """JASPAR accession, e.g. ``MA0079.5``. Versioned: matrices change between releases."""

    name: str
    counts: Mapping[str, tuple[int, ...]]
    """The position frequency matrix: observed counts per base per position."""

    tf_class: tuple[str, ...]
    data_type: MotifDataType
    pubmed_ids: tuple[str, ...]
    relevance: str
    """Why this factor is in the curated set, tied to this project rather than generic."""

    @property
    def length(self) -> int:
        """Return the motif width in base pairs."""
        return len(self.counts["A"])

    @property
    def provenance(self) -> Provenance:
        """Return the matrix's citation, naming the accession a reader can look up."""
        return Provenance.citation(
            f"JASPAR CORE, {self.name} ({self.matrix_id})",
            f"{self.data_type} matrix, {self.length} bp"
            + (f", PMID {self.pubmed_ids[0]}" if self.pubmed_ids else ""),
        )

    @cached_property
    def log_odds(self) -> tuple[Mapping[str, float], ...]:
        """Return per-position log-odds scores against a uniform background."""
        matrix = []
        for position in range(self.length):
            column = {base: self.counts[base][position] for base in BASES}
            total = sum(column.values()) + 4 * PSEUDOCOUNT
            matrix.append(
                {
                    base: math.log2(((column[base] + PSEUDOCOUNT) / total) / 0.25)
                    for base in BASES
                }
            )
        return tuple(matrix)

    @cached_property
    def maximum_score(self) -> float:
        """Return the highest score any sequence could achieve."""
        return sum(max(column.values()) for column in self.log_odds)

    @cached_property
    def minimum_score(self) -> float:
        """Return the lowest score any sequence could achieve."""
        return sum(min(column.values()) for column in self.log_odds)

    def score(self, window: str) -> float:
        """Return the raw log-odds score of ``window``.

        Raises:
            ValueError: If the window length does not match the motif, or holds a base outside
                ACGT. An ``N`` in the reference is not scoreable and must not be scored as if it
                were an A.
        """
        if len(window) != self.length:
            raise ValueError(
                f"motif {self.name} is {self.length} bp but the window is {len(window)} bp"
            )
        total = 0.0
        for column, base in zip(self.log_odds, window, strict=True):
            if base not in column:
                raise ValueError(
                    f"cannot score base {base!r} against motif {self.name}; only ACGT are "
                    f"scoreable, and treating an ambiguity code as a definite base would "
                    f"invent evidence"
                )
            total += column[base]
        return total

    def relative_score(self, window: str) -> float:
        """Return the score of ``window`` scaled to [0, 1] over the matrix's achievable range."""
        span = self.maximum_score - self.minimum_score
        return (self.score(window) - self.minimum_score) / span


@dataclass(frozen=True, slots=True)
class MotifSite:
    """One place a motif matches, on one strand."""

    motif: Motif
    offset: int
    """0-based offset into the scanned sequence where the match starts."""

    strand: Strand
    """Which strand the match is on, relative to the sequence as given."""

    relative_score: float
    raw_score: float
    matched_bases: str
    """The sequence as matched, already reverse-complemented for a minus-strand hit."""

    @property
    def end(self) -> int:
        """Return the 0-based exclusive end of the match."""
        return self.offset + self.motif.length

    @property
    def evidence(self) -> Evidence:
        """Return this site's evidential standing: a prediction, with no cell type."""
        return motif_prediction(self.motif.provenance)

    def __str__(self) -> str:
        """Return a one-line rendering of the match."""
        return (
            f"{self.motif.name} at {self.offset}-{self.end} ({self.strand}) "
            f"score {self.relative_score:.2f}"
        )


@dataclass(frozen=True, slots=True)
class MotifScan:
    """Every match for one motif in one sequence, with its false-positive expectation."""

    motif: Motif
    sites: tuple[MotifSite, ...]
    threshold: float
    sequence_length: int
    expected_hits: float
    """Matches a shuffled sequence of the same composition yields at the same threshold.

    The number that makes a hit count interpretable. Promoter sequence is GC-rich, and so are the
    motifs of most promoter-binding factors, so an impressive-looking count can be entirely
    explained by base composition.
    """

    @property
    def enrichment(self) -> float | None:
        """Return observed hits divided by expected, or ``None`` if none were expected."""
        if self.expected_hits == 0:
            return None
        return len(self.sites) / self.expected_hits

    @property
    def exceeds_chance(self) -> bool:
        """Return whether more sites were found than composition alone explains."""
        return len(self.sites) > self.expected_hits

    @property
    def summary(self) -> str:
        """Return a plain-language statement of what the scan found."""
        if not self.sites:
            return f"No {self.motif.name} match above a relative score of {self.threshold:.2f}."
        verdict = (
            f"{len(self.sites)} match(es)"
            if self.exceeds_chance
            else f"{len(self.sites)} match(es), at or below the {self.expected_hits:.1f} expected "
            f"from base composition alone, so this is not evidence of anything"
        )
        return f"{self.motif.name}: {verdict}."


def load_motifs() -> tuple[Motif, ...]:
    """Load the curated JASPAR selection shipped with this package.

    A subset rather than the whole CORE collection: sixteen factors, each with a stated reason
    tied to promoter-targeted saRNA design, to the two cell lines in use, or to the target genes.
    Scanning all 800-odd CORE matrices over a promoter would produce several hundred hits and
    nothing a researcher could act on.
    """
    source = resources.files(DATA_PACKAGE).joinpath(MATRIX_FILENAME)
    with resources.as_file(source) as path:
        payload = json.loads(path.read_text())
    return tuple(
        Motif(
            matrix_id=entry["matrix_id"],
            name=entry["name"],
            counts={base: tuple(entry["pfm"][base]) for base in BASES},
            tf_class=tuple(entry.get("tf_class") or ()),
            data_type=MotifDataType.parse(entry.get("data_type")),
            pubmed_ids=tuple(str(p) for p in entry.get("pubmed_ids") or ()),
            relevance=entry["relevance"],
        )
        for entry in payload["matrices"].values()
    )


def scan(
    sequence: NucleicSequence,
    motif: Motif,
    threshold: float = DEFAULT_RELATIVE_THRESHOLD,
    shuffles: int = SHUFFLES_FOR_EXPECTATION,
    seed: int = 20261006,
) -> MotifScan:
    """Scan ``sequence`` for ``motif`` on both strands, with a false-positive expectation.

    Args:
        sequence: The region to scan, as DNA.
        motif: The matrix to scan with.
        threshold: Minimum relative score, in [0, 1].
        shuffles: Shuffled replicates used to estimate the expected hit count. Zero skips the
            estimate, which is faster but leaves the result uninterpretable.
        seed: Random seed for the shuffles, so the expectation is reproducible.

    Returns:
        The matches, and how many a composition-matched shuffle would have produced.

    Raises:
        ValueError: If the threshold is outside [0, 1].
    """
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(f"threshold must be a relative score in [0, 1], got {threshold}")
    sites = _sites_in(sequence.bases, motif, threshold)
    expected = _expected_hits(sequence.bases, motif, threshold, shuffles, seed)
    return MotifScan(
        motif=motif,
        sites=sites,
        threshold=threshold,
        sequence_length=len(sequence),
        expected_hits=expected,
    )


def scan_all(
    sequence: NucleicSequence,
    motifs: Sequence[Motif] | None = None,
    threshold: float = DEFAULT_RELATIVE_THRESHOLD,
    shuffles: int = SHUFFLES_FOR_EXPECTATION,
    seed: int = 20261006,
) -> tuple[MotifScan, ...]:
    """Scan ``sequence`` for every motif, most-enriched first.

    Ordered by enrichment over the shuffled expectation rather than by raw hit count, because raw
    counts rank the motifs that match GC-rich sequence most readily, which in a promoter is every
    one of them.
    """
    selected = load_motifs() if motifs is None else motifs
    scans = tuple(scan(sequence, motif, threshold, shuffles, seed) for motif in selected)
    return tuple(
        sorted(scans, key=lambda s: (-(s.enrichment or 0.0), -len(s.sites), s.motif.name))
    )


def _sites_in(bases: str, motif: Motif, threshold: float) -> tuple[MotifSite, ...]:
    """Return every match above ``threshold``, on both strands, in offset order."""
    found: list[MotifSite] = []
    reverse = _reverse_complement(bases)
    length = len(bases)
    for offset in range(length - motif.length + 1):
        window = bases[offset : offset + motif.length]
        if not set(window) <= set(BASES):
            # An N in the reference. Skipped rather than scored, on either strand.
            continue
        for strand in (Strand.PLUS, Strand.MINUS):
            if strand is Strand.PLUS:
                candidate = window
            else:
                start = length - offset - motif.length
                candidate = reverse[start : start + motif.length]
            raw = motif.score(candidate)
            if raw <= MINIMUM_RAW_SCORE:
                # Less likely under the motif than under background, whatever the relative
                # score says. See MINIMUM_RAW_SCORE.
                continue
            relative = motif.relative_score(candidate)
            if relative >= threshold:
                found.append(
                    MotifSite(
                        motif=motif,
                        offset=offset,
                        strand=strand,
                        relative_score=relative,
                        raw_score=raw,
                        matched_bases=candidate,
                    )
                )
    return tuple(found)


def _expected_hits(
    bases: str, motif: Motif, threshold: float, shuffles: int, seed: int
) -> float:
    """Return the mean hit count over shuffles of ``bases``, preserving base composition.

    Shuffling rather than sampling from a background model, so the expectation matches the actual
    sequence's GC content. It destroys dinucleotide structure, which overstates the expectation
    slightly for CpG-dependent motifs — stated rather than corrected, because a dinucleotide-
    preserving shuffle is a larger piece of machinery than this estimate warrants.
    """
    if shuffles <= 0:
        return 0.0
    rng = random.Random(seed)
    characters = list(bases)
    total = 0
    for _ in range(shuffles):
        rng.shuffle(characters)
        total += len(_sites_in("".join(characters), motif, threshold))
    return total / shuffles


def _reverse_complement(bases: str) -> str:
    """Return the reverse complement of ``bases``, leaving unknown characters in place."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(bases))


__all__ = [
    "DEFAULT_RELATIVE_THRESHOLD",
    "JASPAR_LICENCE_NOTE",
    "Motif",
    "MotifDataType",
    "MotifScan",
    "MotifSite",
    "load_motifs",
    "scan",
    "scan_all",
]
