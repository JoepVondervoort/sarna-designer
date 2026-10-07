"""Map the thesis's SCN1A designs back onto the genome, and report which TSS each sits near.

The chapter 4 dataset records target sequences but no coordinates, which is what stops the
location rules — the TSS windows, CpG-island overlap — from being testable against measured
activity. This script recovers the coordinates for the 15 SCN1A designs by searching the locus,
then reports each one's offset from **every** annotated TSS rather than from an assumed one.

Needs the network (Ensembl). Run with::

    .venv/bin/python scripts/locate_scn1a_designs.py

Written up in docs/scientific/scn1a-design-coordinates.md.
"""


from sarna import settings
from sarna.genome.coordinates import GenomicInterval
from sarna.genome.ensembl import EnsemblClient
from sarna.validation.dataset import ValidationSequence, load_chapter4_dataset

CACHE = settings.cache_dir("ensembl")
LOCUS_PADDING_BP = 25_000
WANG_WINDOW = (-1000, -100)
"""The admission window from Wang et al. 2015 Table 1, as TSS-relative positions."""


def main() -> None:
    """Locate every SCN1A design and report its position relative to the nearest TSS."""
    client = EnsemblClient(cache_dir=CACHE)
    gene = client.gene("SCN1A")
    transcripts = client.transcripts(gene)
    mane = next((t for t in transcripts if t.is_mane_select), None)

    locus = GenomicInterval(
        assembly=client.assembly,
        contig=gene.contig,
        start=gene.interval.start - LOCUS_PADDING_BP,
        end=gene.interval.end + LOCUS_PADDING_BP,
    )
    reference = client.sequence(locus).bases

    designs = [s for s in load_chapter4_dataset() if s.gene == "SCN1A"]
    all_tss = sorted({t.tss_coordinate for t in transcripts})
    print(f"{gene.symbol} {gene.ensembl_id} chr{gene.contig} strand {gene.strand}")
    print(f"  {len(transcripts)} transcripts, {len(all_tss)} distinct TSSs")
    mane_label = f"{mane.ensembl_id} at {mane.tss_coordinate}" if mane else "none"
    print(f"  MANE Select: {mane_label}")
    print(f"  searching {len(designs)} designs in {locus.length} bp\n")

    header = (
        f"{'design':10s} {'act':4s} {'5p coord':>10s} {'nearest TSS':>12s} "
        f"{'offset':>8s} {'in window':>10s}  transcripts"
    )
    print(header)
    print("-" * len(header))
    rows = []
    for design in designs:
        five_prime = _locate(design, reference, locus)
        if five_prime is None:
            print(f"{design.name:10s} {'':4s} {'not found':>10s}")
            continue
        offsets = {t.ensembl_id: t.tss_coordinate - five_prime for t in transcripts}
        nearest_id = min(offsets, key=lambda k: abs(offsets[k]))
        nearest_tss = next(t.tss_coordinate for t in transcripts if t.ensembl_id == nearest_id)
        offset = offsets[nearest_id]
        sharing = sorted(t.ensembl_id for t in transcripts if t.tss_coordinate == nearest_tss)
        inside = WANG_WINDOW[0] <= offset <= WANG_WINDOW[1]
        rows.append((design, five_prime, nearest_tss, offset, inside))
        print(
            f"{design.name:10s} {'ACT' if design.is_activating else '':4s} {five_prime:>10d} "
            f"{nearest_tss:>12d} {offset:>+8d} {'yes' if inside else 'NO':>10s}  "
            f"{sharing[0]}{f' +{len(sharing) - 1}' if len(sharing) > 1 else ''}"
        )

    _summarise(rows, mane.tss_coordinate if mane else None)


def _locate(
    design: ValidationSequence, reference: str, locus: GenomicInterval
) -> int | None:
    """Return the design's 5' genomic coordinate in transcript orientation, or ``None``.

    SCN1A is on the minus strand, so a target reading 5'->3' in transcript orientation appears
    reverse-complemented in the plus-strand reference, and its 5' end is the **highest**
    coordinate of the match.
    """
    plus_strand_form = design.target_dna.reverse_complement().bases
    index = reference.find(plus_strand_form)
    if index < 0:
        if (forward := reference.find(design.target_dna.bases)) >= 0:
            return locus.start + forward
        return None
    return locus.start + index + len(design.target_dna) - 1


def _summarise(
    rows: list[tuple[ValidationSequence, int, int, int, bool]], mane_tss: int | None
) -> None:
    """Group the designs by the TSS they sit nearest and report the clusters."""
    print()
    clusters: dict[int, list[tuple[ValidationSequence, int, bool]]] = {}
    for design, _five_prime, nearest_tss, offset, inside in rows:
        clusters.setdefault(nearest_tss, []).append((design, offset, inside))

    print(f"{len(clusters)} cluster(s), by the TSS each design sits nearest:\n")
    for tss in sorted(clusters, reverse=True):
        members = clusters[tss]
        activators = [d.name for d, _o, _i in members if d.is_activating]
        offsets = sorted(o for _d, o, _i in members)
        in_window = sum(1 for _d, _o, inside in members if inside)
        marker = "  <- MANE Select" if tss == mane_tss else ""
        print(f"  TSS {tss}{marker}")
        print(f"    {len(members)} designs, offsets {offsets[0]:+d} to {offsets[-1]:+d}")
        print(f"    {in_window}/{len(members)} inside Wang's -1000..-100 window")
        print(f"    activators: {', '.join(activators) if activators else 'none'}")
    if mane_tss is not None and mane_tss not in clusters:
        print(
            "\n  No design sits nearest the MANE Select TSS. The tool's default transcript is\n"
            "  therefore not the one this project's own SCN1A work targeted."
        )


if __name__ == "__main__":
    main()
