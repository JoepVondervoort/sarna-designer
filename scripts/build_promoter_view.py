"""Build a complete promoter view for one gene and write it as JSON.

The end-to-end path: gene symbol -> transcript -> promoter -> candidates -> annotations -> one JSON
document a renderer can draw without doing any biology.

Needs the network for the first run; afterwards Ensembl responses come from `.cache/ensembl`.

    .venv/bin/python scripts/build_promoter_view.py SCN1A
    .venv/bin/python scripts/build_promoter_view.py CHD2 --preset exploratory --out view.json

The transcript is MANE Select, or the script reports the alternatives and stops (ADR 0012).
"""

import argparse
import json
import sys
from pathlib import Path

from sarna import settings
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import PRESETS, preset
from sarna.genome.anchors import load_anchors
from sarna.genome.ensembl import (
    DEFAULT_DOWNSTREAM_BP,
    DEFAULT_UPSTREAM_BP,
    EnsemblClient,
    EnsemblError,
    NoManeSelectTranscriptError,
    Transcript,
)
from sarna.regulatory.chromatin import open_af22_reader
from sarna.regulatory.features import fetch_regulatory_features
from sarna.regulatory.motifs import scan_all
from sarna.view.html import render as render_html
from sarna.view.promoter_view import (
    PromoterView,
    alternative_tss_spans,
    cell_panels_for,
)

CACHE = settings.cache_dir("ensembl")


def main(argv: list[str] | None = None) -> int:
    """Build and write a promoter view. Returns a process exit code."""
    args = _parse(argv)
    client = EnsemblClient(cache_dir=CACHE)
    anchors = load_anchors()

    if args.list_anchors:
        declared = anchors.for_gene(args.symbol)
        if not declared:
            print(
                f"No promoter anchor is declared for {args.symbol}; the annotated MANE Select TSS "
                f"will be used. Genes with declared anchors: "
                f"{', '.join(anchors.genes) or 'none'}."
            )
            return 0
        for entry in declared:
            print(f"{entry.id}  [{entry.reference_kind}]  {entry.provenance}")
            if entry.reference_kind.caveat:
                print(f"    {entry.reference_kind.caveat}")
        return 0

    try:
        gene = client.gene(args.symbol)
        transcript: Transcript | None = None
        if args.anchor:
            anchor = anchors[args.anchor]
            if anchor.gene.upper() != gene.symbol.upper():
                print(
                    f"anchor {anchor.id!r} belongs to {anchor.gene}, not {gene.symbol}",
                    file=sys.stderr,
                )
                return 2
            promoter = client.promoter_for_anchor(anchor, args.upstream, args.downstream)
            print(f"  anchored on {anchor.id} [{anchor.reference_kind}]: {anchor.provenance}")
            if anchor.reference_kind.caveat:
                print(f"  CAVEAT {anchor.reference_kind.caveat}")
        else:
            transcript = client.select_transcript(gene)
            promoter = client.promoter_for(transcript, args.upstream, args.downstream)
            others = anchors.for_gene(gene.symbol)
            if others:
                print(
                    f"  note: {len(others)} literature promoter(s) are declared for "
                    f"{gene.symbol} and were NOT used. Run with --list-anchors to see them."
                )
    except NoManeSelectTranscriptError as error:
        # A question, not a failure: the promoter depends on a choice the tool will not make.
        print(error, file=sys.stderr)
        return 2
    except EnsemblError as error:
        print(f"Ensembl: {error}", file=sys.stderr)
        return 1

    config = preset(args.preset).resolve()
    result = design_sarnas(promoter, config)

    print(f"{gene.symbol} {gene.ensembl_id}  chr{gene.contig} strand {gene.strand}")
    if transcript is not None:
        print(f"  transcript {transcript.ensembl_id} ({transcript.refseq_match}), MANE Select")
    print(f"  promoter {promoter.length} bp, {promoter.sequence.gc_percent:.1f}% GC")
    print(
        f"  {len(result.generated)} candidates generated, {len(result.passing)} passing, "
        f"{len(result.rejected)} rejected"
    )

    transcripts = client.transcripts(gene)
    alternatives = alternative_tss_spans(promoter, transcripts, transcript)
    print(
        f"  {len({t.tss_coordinate for t in transcripts})} distinct TSSs in this gene, "
        f"{len(alternatives)} of them inside this window"
    )

    assert promoter.interval is not None
    regulatory = fetch_regulatory_features(client, promoter.interval)
    print(f"  {regulatory.summary}")

    motif_scans = scan_all(promoter.sequence, shuffles=args.shuffles)
    above_chance = [s for s in motif_scans if s.exceeds_chance]
    print(
        f"  {len(motif_scans)} motifs scanned, {len(above_chance)} above what base composition "
        f"alone explains"
    )
    for scan in above_chance[:5]:
        print(f"      {scan.motif.name:7s} {len(scan.sites):4d} sites, {scan.summary}")

    chromatin = None if args.no_chromatin else open_af22_reader()
    if chromatin is not None:
        marks = ", ".join(str(m) for m in chromatin.marks)
        print(f"  AF22 chromatin tracks loaded: {marks}")
    else:
        print("  AF22 chromatin tracks not loaded (absent, or --no-chromatin)")

    view = PromoterView(
        result=result,
        gene=gene,
        transcript=transcript,
        alternative_tss=alternatives,
        regulatory=regulatory,
        motif_scans=motif_scans,
        cell_panels=cell_panels_for(
            promoter=promoter, chromatin_reader=chromatin
        ),
    )
    payload = view.to_dict(full_rules_for_rejected=args.full_rules)
    destination = args.out or Path(f"{gene.symbol.lower()}_promoter_view.json")
    destination.write_text(json.dumps(payload, indent=1))
    size_kb = destination.stat().st_size / 1024
    print(f"\nwrote {destination} ({size_kb:.0f} KB)")
    if args.html is not None:
        # --full-rules reached the JSON and not the page, so asking for full detail silently got
        # an abbreviated one. It is off by default here for size, not because it is unwanted.
        args.html.write_text(
            render_html(view, full_rules_for_rejected=args.full_rules), encoding="utf-8"
        )
        print(f"wrote {args.html} ({args.html.stat().st_size / 1024:.0f} KB, self-contained)")
    print("  top rejection reasons:")
    for entry in payload["rejection_breakdown"][:4]:
        print(
            f"      {entry['rule']:28s} {entry['rejected']:5d} "
            f"({entry['percent_of_generated']}% of generated)"
        )
    print(
        f"  tracks: {len(payload['tracks']['cpg_islands'])} CpG islands, "
        f"{len(payload['tracks']['alternative_tss'])} alternative TSSs, "
        f"{len(payload['tracks']['regulatory'])} regulatory features, "
        f"{len(payload['tracks']['motifs'])} motif sites"
    )
    return 0


def _parse(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("symbol", help="HGNC gene symbol, e.g. SCN1A")
    parser.add_argument(
        "--preset",
        default="wang-2015",
        choices=sorted(PRESETS),
        help="rule preset to design with (default: wang-2015)",
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=DEFAULT_UPSTREAM_BP,
        help=f"bases upstream of the TSS (default: {DEFAULT_UPSTREAM_BP})",
    )
    parser.add_argument(
        "--downstream",
        type=int,
        default=DEFAULT_DOWNSTREAM_BP,
        help=f"bases downstream of the TSS (default: {DEFAULT_DOWNSTREAM_BP})",
    )
    parser.add_argument(
        "--shuffles",
        type=int,
        default=10,
        help="shuffles for each motif's false-positive expectation (default: 10)",
    )
    parser.add_argument(
        "--full-rules",
        action="store_true",
        help=(
            "include every rule's measurements for rejected candidates too, in both the JSON and "
            "the HTML page. Much larger: about 21 MB for a 2.5 kb promoter against under 2 MB. "
            "Rejected candidates are always present either way, with the rules that excluded "
            "them; what this adds is the rules they *passed* and every measurement behind both."
        ),
    )
    parser.add_argument(
        "--no-chromatin",
        action="store_true",
        help="skip AF22 histone tracks even if they are present",
    )
    parser.add_argument(
        "--anchor",
        default=None,
        help=(
            "design against a declared promoter anchor by id instead of the MANE Select TSS. "
            "Use --list-anchors to see what is declared for a gene."
        ),
    )
    parser.add_argument(
        "--list-anchors",
        action="store_true",
        help="list the promoter anchors declared for this gene and exit",
    )
    parser.add_argument("--out", type=Path, default=None, help="output JSON path")
    parser.add_argument(
        "--html",
        type=Path,
        default=None,
        help="also write a self-contained interactive HTML page to this path",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
