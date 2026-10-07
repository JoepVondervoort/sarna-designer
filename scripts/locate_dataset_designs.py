"""Recover genomic coordinates for every assayed design, so the location rules become testable.

`docs/scientific/empirical-validation.md` names this as the single most valuable addition to the
dataset and the only one that needs no new experiments. The chapter 4 appendix tabulates target
sequences without coordinates, so the TSS windows, CpG-island overlap and the 3' flank preferences
report not-applicable for all 125 — five of the seven preferences and three of the gates.

Each 19-mer is searched for in its gene's locus plus padding, in both orientations, and its offset
is computed against **every** annotated TSS rather than an assumed one. A design is only accepted
when it matches exactly once; an ambiguous or absent match is reported, never guessed.

Needs the network. Writes a CSV beside the dataset so the result is reusable offline::

    .venv/bin/python scripts/locate_dataset_designs.py
    .venv/bin/python scripts/locate_dataset_designs.py --gene CHD2 --dry-run
"""

import argparse
import csv
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from sarna import settings
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import (
    EnsemblClient,
    EnsemblError,
    Gene,
    Transcript,
)
from sarna.sequence.nucleic import NucleicSequence
from sarna.validation.dataset import ValidationSequence, load_chapter4_dataset

CACHE = settings.cache_dir("ensembl")
OUTPUT = Path("src/sarna/validation/data/vondervoort_2026_chapter4_coordinates.csv")
CONTEXT_FLANK_BP = 250
"""Genomic context stored either side of each target, in transcript orientation.

Enough to make the 3' flank preferences (positions 20-23) evaluable and to show local GC and CpG
density, at about 64 KB for the whole dataset. Not enough for CpG *island* detection, which needs
500 bp of island before it can call one.
"""

LOCUS_PADDING_BP = 25_000
"""Searched either side of the gene span. Promoter-targeting designs sit upstream of the TSS,
which for a long gene can be well outside the annotated span."""

ALIASES = {
    "P21": "CDKN1A",
    "P53": "TP53",
    "VEGF": "VEGFA",
}
"""Dataset names that are not HGNC symbols. Ensembl matches HGNC exactly."""


@dataclass(frozen=True, slots=True)
class Located:
    """One design placed on the genome, with its position relative to the nearest TSS."""

    design: ValidationSequence
    symbol: str
    contig: str
    five_prime: int
    """Genomic coordinate of the target's 5' base, in the gene's own orientation."""

    strand: Strand
    context: str
    """``CONTEXT_FLANK_BP`` either side of the target, 5'->3' in transcript orientation."""

    context_offset: int
    """0-based offset of the target within :attr:`context`."""

    mane_transcript: str | None
    mane_offset: int | None
    """Offset from the MANE Select TSS, the anchor ADR 0012 makes the tool's default."""

    nearest_transcript: str
    nearest_tss: int
    nearest_offset: int
    """Offset from whichever TSS is closest. Weak for a gene with dozens of annotated TSSs."""

    distinct_tss: int

    @property
    def offset(self) -> int:
        """Return the MANE Select offset where there is one, else the nearest."""
        return self.mane_offset if self.mane_offset is not None else self.nearest_offset


def main(argv: list[str] | None = None) -> int:
    """Locate every design and write the coordinate table. Returns a process exit code."""
    args = _parse(argv)
    client = EnsemblClient(cache_dir=CACHE)
    dataset = load_chapter4_dataset()
    genes = sorted({s.gene for s in dataset})
    if args.gene:
        genes = [g for g in genes if g.upper() == args.gene.upper()]
        if not genes:
            print(f"no gene {args.gene!r} in the dataset", file=sys.stderr)
            return 2

    located: list[Located] = []
    unresolved: list[tuple[str, str]] = []
    for symbol in genes:
        designs = [s for s in dataset if s.gene == symbol]
        resolved = ALIASES.get(symbol, symbol)
        try:
            gene = client.gene(resolved)
            transcripts = client.transcripts(gene)
        except EnsemblError as error:
            print(f"{symbol:9s} lookup failed ({type(error).__name__})")
            unresolved.extend((d.name, "gene lookup failed") for d in designs)
            continue
        reference = _locus_sequence(client, gene)
        found, missing = _place(designs, reference, gene, transcripts)
        located.extend(found)
        unresolved.extend(missing)
        in_window = sum(1 for item in found if -1000 <= item.nearest_offset <= -100)
        print(
            f"{symbol:9s} {resolved:9s} {len(found):2d}/{len(designs):2d} located, "
            f"{in_window:2d} inside -1000..-100, "
            f"{len({t.tss_coordinate for t in transcripts})} TSSs"
        )

    print(f"\n{len(located)} of {len(dataset)} designs located; {len(unresolved)} not")
    for name, why in unresolved:
        print(f"    {name}: {why}")

    _report(located)
    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0
    _write(located, args.out or OUTPUT)
    return 0


def _locus_sequence(client: EnsemblClient, gene: Gene) -> tuple[str, GenomicInterval]:
    """Fetch the gene's locus plus padding, as the plus strand stores it."""
    contig_length = gene.interval.assembly.contig_length(gene.contig)
    interval = GenomicInterval(
        assembly=gene.interval.assembly,
        contig=gene.contig,
        start=max(0, gene.interval.start - LOCUS_PADDING_BP),
        end=min(contig_length, gene.interval.end + LOCUS_PADDING_BP),
    )
    return client.sequence(interval).bases, interval


def _place(
    designs: Sequence[ValidationSequence],
    reference: tuple[str, GenomicInterval],
    gene: Gene,
    transcripts: Sequence[Transcript],
) -> tuple[list[Located], list[tuple[str, str]]]:
    """Locate each design in the locus, refusing ambiguous matches."""
    bases, interval = reference
    found: list[Located] = []
    missing: list[tuple[str, str]] = []
    for design in designs:
        target = design.target_dna.bases
        reverse = design.target_dna.reverse_complement().bases
        forward_hits = _all_indices(bases, target)
        reverse_hits = _all_indices(bases, reverse)
        total = len(forward_hits) + len(reverse_hits)
        if total == 0:
            missing.append((design.name, "not found in the locus"))
            continue
        if total > 1:
            missing.append((design.name, f"matches {total} places; ambiguous"))
            continue
        if forward_hits:
            strand = Strand.PLUS
            index = forward_hits[0]
            five_prime = interval.start + index
        else:
            strand = Strand.MINUS
            index = reverse_hits[0]
            five_prime = interval.start + index + len(target) - 1
        context, context_offset = _context(bases, index, len(target), strand)
        # Upstream is negative on both strands: lower coordinates on the plus strand, higher on
        # the minus. The sign multiplier goes on the *difference*, not on the result of the
        # already-oriented subtraction, which is easy to get backwards.
        offsets = {
            t.ensembl_id: (five_prime - t.tss_coordinate) * strand.sign for t in transcripts
        }
        nearest_id = min(offsets, key=lambda key: abs(offsets[key]))
        nearest = next(t for t in transcripts if t.ensembl_id == nearest_id)
        mane = next((t for t in transcripts if t.is_mane_select), None)
        found.append(
            Located(
                design=design,
                symbol=gene.symbol,
                contig=gene.contig,
                five_prime=five_prime,
                strand=strand,
                context=context,
                context_offset=context_offset,
                mane_transcript=mane.ensembl_id if mane else None,
                mane_offset=offsets[mane.ensembl_id] if mane else None,
                nearest_transcript=nearest_id,
                nearest_tss=nearest.tss_coordinate,
                nearest_offset=offsets[nearest_id],
                distinct_tss=len({t.tss_coordinate for t in transcripts}),
            )
        )
    return found, missing


def _context(bases: str, index: int, length: int, strand: Strand) -> tuple[str, int]:
    """Return genomic context around a match, oriented 5'->3' for the transcript.

    The reference is the plus strand, so a minus-strand target's context is reverse-complemented,
    which also moves the target's offset within it to the other end. Clipped at the locus edges,
    so the returned offset is computed rather than assumed to be ``CONTEXT_FLANK_BP``.
    """
    low = max(0, index - CONTEXT_FLANK_BP)
    high = min(len(bases), index + length + CONTEXT_FLANK_BP)
    window = bases[low:high]
    if strand is Strand.PLUS:
        return window, index - low
    flipped = NucleicSequence.parse(window).reverse_complement().bases
    # After reversing, the bases that were 3' of the target are now 5' of it.
    return flipped, high - (index + length)


def _all_indices(haystack: str, needle: str) -> list[int]:
    """Return every index at which ``needle`` occurs, so a second match is detectable."""
    hits: list[int] = []
    start = haystack.find(needle)
    while start >= 0:
        hits.append(start)
        start = haystack.find(needle, start + 1)
    return hits


def _report(located: Sequence[Located]) -> None:
    """Summarise how the located designs fall relative to the rule window."""
    if not located:
        return
    _window_report(located, "MANE Select TSS", lambda item: item.mane_offset)
    _window_report(located, "nearest annotated TSS", lambda item: item.nearest_offset)

    print("\nFurthest upstream of any TSS:")
    for item in sorted(located, key=lambda i: i.nearest_offset)[:6]:
        label = "ACT" if item.design.is_activating else ""
        print(
            f"    {item.design.name:10s} {label:4s} {item.nearest_offset:+7d} "
            f"from {item.nearest_transcript}"
        )

    downstream = [item for item in located if item.nearest_offset > 0]
    if downstream:
        print(f"\n{len(downstream)} design(s) lie downstream of their nearest TSS:")
        for item in sorted(downstream, key=lambda i: -i.nearest_offset)[:6]:
            label = "ACT" if item.design.is_activating else ""
            print(f"    {item.design.name:10s} {label:4s} {item.nearest_offset:+7d}")

    crowded = sorted({(i.symbol, i.distinct_tss) for i in located}, key=lambda e: -e[1])[:4]
    print(
        "\nThe nearest-TSS figure is weak where a gene has many: "
        + ", ".join(f"{symbol} {count}" for symbol, count in crowded)
        + ". With dozens of annotated TSSs almost any position is near one, which is why the"
        " MANE Select offset is reported alongside it (ADR 0012)."
    )


def _window_report(
    located: Sequence[Located], label: str, offset_of: Callable[[Located], int | None]
) -> None:
    """Report how many designs fall inside Wang's window against one TSS anchor."""
    measurable = [item for item in located if offset_of(item) is not None]
    inside = [item for item in measurable if _in_window(offset_of(item))]
    activators = [item for item in measurable if item.design.is_activating]
    activators_inside = [item for item in inside if item.design.is_activating]
    print(
        f"\nRelative to the {label}: {len(inside)}/{len(measurable)} designs inside Wang's "
        f"-1000..-100 window, {len(activators_inside)}/{len(activators)} activators"
    )


def _in_window(offset: int | None) -> bool:
    """Return whether an offset lies in the window Wang et al. Table 1 admits."""
    return offset is not None and -1000 <= offset <= -100


def _write(located: Sequence[Located], destination: Path) -> None:
    """Write the coordinate table."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "name",
                "gene",
                "assembly",
                "contig",
                "five_prime_coordinate",
                "strand",
                "context_sequence",
                "target_offset_in_context",
                "mane_transcript",
                "mane_tss_offset",
                "nearest_transcript",
                "nearest_tss_coordinate",
                "nearest_tss_offset",
                "distinct_tss_in_gene",
            ]
        )
        for item in sorted(located, key=lambda i: (i.symbol, i.design.name)):
            writer.writerow(
                [
                    item.design.name,
                    item.symbol,
                    str(Assembly.GRCH38),
                    item.contig,
                    item.five_prime,
                    str(item.strand),
                    item.context,
                    item.context_offset,
                    item.mane_transcript or "",
                    "" if item.mane_offset is None else item.mane_offset,
                    item.nearest_transcript,
                    item.nearest_tss,
                    item.nearest_offset,
                    item.distinct_tss,
                ]
            )
    print(f"\nwrote {destination} ({len(located)} rows)")


def _parse(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", default=None, help="limit to one gene from the dataset")
    parser.add_argument("--out", type=Path, default=None, help="output CSV path")
    parser.add_argument("--dry-run", action="store_true", help="report without writing")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
