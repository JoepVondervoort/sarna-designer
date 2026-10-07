"""Read AF22 chromatin state at every chapter 4 target promoter, with positive controls.

Chapter 4's largest effect was not a sequence effect: no saRNA activated any locus in AF22,
including the P21 control that works reliably in HeLa, and the discussion proposed that AF22
chromatin is more restrictive. The AF22 histone tracks (GSE108506) make that testable per locus.

The housekeeping promoters are not decoration. Without GAPDH and ACTB showing strong H3K4me3, a
null result at a target promoter is indistinguishable from the adapter reading the wrong place —
and it would be reading the wrong place by about 856 kb, since the tracks are GRCh37 and the
engine is GRCh38 (defect B15).

Needs the network, and the tracks in ``~/af22_data`` or ``$SARNA_AF22_DIR``::

    .venv/bin/python scripts/af22_chromatin_at_targets.py

Written up in docs/scientific/af22-chromatin.md.
"""

import sys

from sarna import settings
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.ensembl import EnsemblClient, EnsemblError
from sarna.regulatory.chromatin import (
    ENRICHMENT_THRESHOLD,
    Af22ChromatinReader,
    HistoneMark,
    open_af22_reader,
)

CACHE = settings.cache_dir("ensembl")
WINDOW_BP = 500
"""Half-width around the TSS. H3K4me3 peaks sharply there; a wider window dilutes it."""

CONTROLS = ("GAPDH", "ACTB")
"""Transcribed in every cell type. If these are not marked, the measurement is broken."""

TARGETS = ("CDKN1A", "CHD2", "DEPDC5", "NPRL2", "NPRL3", "TSC1", "TSC2", "SCN1A")
"""Chapter 4's six targets, the P21 positive control, and SCN1A from chapters 2-3."""

SCN1A_TSS = {
    "MANE Select": 166_128_019,
    "cluster A (5 designs)": 166_149_213,
    "cluster B (5 designs)": 166_073_669,
    "cluster C (5 designs, no annotated TSS)": 166_079_200,
}
"""The TSSs the thesis's SCN1A designs actually cluster around, per scn1a-design-coordinates.md."""

REPORTED = (HistoneMark.H3K4ME3, HistoneMark.H3K27AC, HistoneMark.H3K27ME3)


def main() -> int:
    """Report chromatin state at the controls, the targets, and SCN1A's several TSSs."""
    reader = open_af22_reader()
    if reader is None:
        print(
            "AF22 tracks not found. Point SARNA_AF22_DIR at them, or place the GSE108506 "
            "bigWigs in ~/af22_data.",
            file=sys.stderr,
        )
        return 1
    client = EnsemblClient(cache_dir=CACHE)

    print(f"AF22 lt-NES chromatin, +/-{WINDOW_BP} bp of each TSS, as fold over a same-contig")
    print(
        f"background of 250 random windows. A mark is called present above "
        f"{ENRICHMENT_THRESHOLD}x."
    )
    print("Tracks are GRCh37/hg19; queries are GRCh38 and converted per request.\n")

    print("Positive controls — these must be strongly marked or nothing below means anything:")
    _table(reader, client, CONTROLS)

    print("\nTarget promoters:")
    _table(reader, client, TARGETS)

    print("\nSCN1A's several TSSs, including the two the thesis designed against:")
    header = f"  {'TSS':40s} " + " ".join(f"{m!s:>10s}" for m in REPORTED) + "   verdict"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for label, tss in SCN1A_TSS.items():
        interval = GenomicInterval(
            Assembly.GRCH38, "2", tss - WINDOW_BP, tss + WINDOW_BP, Strand.MINUS
        )
        profile = reader.profile(interval)
        folds = []
        for mark in REPORTED:
            signal = profile.signal(mark)
            folds.append(f"{signal.fold_over_background:9.2f}x" if signal else "        —")
        verdict = "open" if profile.activating_marks else "closed"
        print(f"  {label:40s} " + " ".join(folds) + f"   {verdict}")

    print(
        "\nRead the controls first. If GAPDH and ACTB are above 20x and a target is near 1x, the\n"
        "target's promoter carries no active mark in AF22. See docs/scientific/af22-chromatin.md."
    )
    return 0


def _table(reader: Af22ChromatinReader, client: EnsemblClient, symbols: tuple[str, ...]) -> None:
    """Print one row per gene, at its MANE Select TSS."""
    header = f"  {'gene':10s} " + " ".join(f"{m!s:>10s}" for m in REPORTED) + "   verdict"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for symbol in symbols:
        try:
            gene = client.gene(symbol)
            transcript = client.select_transcript(gene)
        except EnsemblError as error:
            print(f"  {symbol:10s} lookup failed: {type(error).__name__}")
            continue
        interval = GenomicInterval(
            Assembly.GRCH38,
            gene.contig,
            transcript.tss_coordinate - WINDOW_BP,
            transcript.tss_coordinate + WINDOW_BP,
            Strand.PLUS,
        )
        profile = reader.profile(interval)
        if not profile.is_readable:
            print(f"  {symbol:10s} {profile.problem}: {profile.problem_detail}")
            continue
        folds = []
        for mark in REPORTED:
            signal = profile.signal(mark)
            folds.append(f"{signal.fold_over_background:9.2f}x" if signal else "        —")
        verdict = "open" if profile.activating_marks else "CLOSED"
        print(f"  {symbol:10s} " + " ".join(folds) + f"   {verdict}")


if __name__ == "__main__":
    raise SystemExit(main())
