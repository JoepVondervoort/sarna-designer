"""Measure local genomic-report processing without submitting work to NCBI.

Run: python scripts/benchmark_offtarget.py --queries 25 --alignments 120
This measures parsing and cached annotation overhead, not NCBI queue/search latency.
"""

import argparse
import json
from time import perf_counter

from sarna.offtarget.ncbi import parse_report
from sarna.offtarget.proximity import BlockedTranscripts, place_hit


def main() -> None:
    """Build a deterministic chromosome-level report and print measured wall times."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", type=int, default=25, choices=range(1, 51), metavar="1..50")
    parser.add_argument("--alignments", type=int, default=120,
                        choices=range(1, 1001), metavar="1..1000")
    args = parser.parse_args()
    iterations = []
    for query in range(args.queries):
        hsps = "".join(
            f"<Hsp><Hsp_identity>17</Hsp_identity><Hsp_align-len>19</Hsp_align-len>"
            f"<Hsp_gaps>0</Hsp_gaps><Hsp_query-from>1</Hsp_query-from>"
            f"<Hsp_query-to>19</Hsp_query-to><Hsp_hit-from>{1000 + n * 5000}</Hsp_hit-from>"
            f"<Hsp_hit-to>{1018 + n * 5000}</Hsp_hit-to><Hsp_hit-frame>1</Hsp_hit-frame>"
            f"<Hsp_bit-score>30</Hsp_bit-score><Hsp_evalue>1</Hsp_evalue></Hsp>"
            for n in range(args.alignments)
        )
        iterations.append(
            f"<Iteration><Iteration_query-def>candidate-{query}</Iteration_query-def>"
            f"<Iteration_hits><Hit><Hit_accession>NC_000017.11</Hit_accession>"
            f"<Hit_def>Human chromosome 17</Hit_def><Hit_hsps>{hsps}</Hit_hsps></Hit>"
            f"</Iteration_hits></Iteration>"
        )
    raw = ("<BlastOutput><BlastOutput_iterations>" + "".join(iterations)
           + "</BlastOutput_iterations></BlastOutput>").encode()
    start = perf_counter()
    report = parse_report(raw)
    parsed_at = perf_counter()
    lookup = BlockedTranscripts(lambda _: [])
    for hit in report.hits[:300]:
        place_hit(hit, lookup, (-1000, -100))
    end = perf_counter()
    print(json.dumps({
        "queries": len(report.queries or ()),
        "alignments_in_report": args.queries * args.alignments,
        "retained_alignments": len(report.hits), "report_bytes": len(raw),
        "parse_seconds": round(parsed_at - start, 4),
        "cached_placement_seconds": round(end - parsed_at, 4),
        "underlying_annotation_lookups": lookup.requests,
        "remote_requests": 0, "includes_ncbi_queue_time": False,
    }, indent=2))


if __name__ == "__main__":
    main()
