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

from sarna.design.pipeline import design_sarnas
from sarna.design.presets import PRESETS, preset
from sarna.genome.ensembl import (
    DEFAULT_DOWNSTREAM_BP,
    DEFAULT_UPSTREAM_BP,
    EnsemblClient,
    EnsemblError,
    NoManeSelectTranscriptError,
)
from sarna.regulatory.features import fetch_regulatory_features
from sarna.regulatory.motifs import scan_all
from sarna.view.promoter_view import (
    PromoterView,
    alternative_tss_spans,
    cell_panels_for,
)

CACHE = Path(".cache/ensembl")


def main(argv: list[str] | None = None) -> int:
    """Build and write a promoter view. Returns a process exit code."""
    args = _parse(argv)
    client = EnsemblClient(cache_dir=CACHE)
    try:
        gene = client.gene(args.symbol)
        transcript = client.select_transcript(gene)
    except NoManeSelectTranscriptError as error:
        # A question, not a failure: the promoter depends on a choice the tool will not make.
        print(error, file=sys.stderr)
        return 2
    except EnsemblError as error:
        print(f"Ensembl: {error}", file=sys.stderr)
        return 1

    promoter = client.promoter_for(transcript, args.upstream, args.downstream)
    config = preset(args.preset).resolve()
    result = design_sarnas(promoter, config)

    print(f"{gene.symbol} {gene.ensembl_id}  chr{gene.contig} strand {gene.strand}")
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

    view = PromoterView(
        result=result,
        gene=gene,
        transcript=transcript,
        alternative_tss=alternatives,
        regulatory=regulatory,
        motif_scans=motif_scans,
        cell_panels=cell_panels_for(),
    )
    payload = view.to_dict(full_rules_for_rejected=args.full_rules)
    destination = args.out or Path(f"{gene.symbol.lower()}_promoter_view.json")
    destination.write_text(json.dumps(payload, indent=1))
    size_kb = destination.stat().st_size / 1024
    print(f"\nwrote {destination} ({size_kb:.0f} KB)")
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
            "include every rule's measurements for rejected candidates too. Much larger: about "
            "21 MB for a 2.5 kb promoter against under 2 MB. Rejected candidates are always "
            "present either way, with the rules that excluded them."
        ),
    )
    parser.add_argument("--out", type=Path, default=None, help="output JSON path")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
