"""Full-length, ungapped local screening against a pinned human reference.

Bowtie's ``-v 3 -a --best`` enumerates all strata on both strands. No seed-only policy,
``--strata``, ``-m`` or ``-k`` is used. Completeness is limited to this explicit mismatch
policy, callable reference bases, and the recorded reference release; it is not biological
specificity, a bulge/gap search, a variant-genome search, or an RNA seed-effect prediction.
"""

import gzip
import hashlib
import json
import sqlite3
import subprocess
import tempfile
import time
from collections import defaultdict
from collections.abc import Collection, Mapping, Sequence
from contextlib import ExitStack, closing
from dataclasses import replace
from pathlib import Path
from typing import Any, BinaryIO

from sarna import __version__
from sarna.offtarget.models import (
    OffTargetHit,
    PromoterProximity,
    ScreeningRecord,
    SearchParameters,
)
from sarna.offtarget.ncbi import candidate_sequences
from sarna.offtarget.ontarget import loci_from_candidates, split_on_target
from sarna.offtarget.references import TSS_PROXY_FEATURES, ReferenceManager, reference_bases
from sarna.offtarget.screening import ProgressReporter, specificity_table
from sarna.reproducibility import software_identity

CHROMOSOMES = frozenset(str(i) for i in range(1, 23)) | {"X", "Y", "MT"}
DETAIL_LIMIT = 200
SAFETY_HIT_LIMIT = 250_000
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    """Return the DNA reverse complement without changing coordinate conventions."""
    return sequence.translate(_COMPLEMENT)[::-1]


def parse_alignment(
    line: str, query: str, title: str, offset: int, fasta: BinaryIO
) -> OffTargetHit:
    """Validate a native alignment against actual reference bases, including reverse strands."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 8 or fields[1] not in ("+", "-"):
        raise ValueError("Malformed local aligner output; no screening result was saved.")
    start = int(fields[3])
    subject = reference_bases(fasta, offset, start, len(query))
    if fields[1] == "-":
        subject = reverse_complement(subject)
    mismatch = tuple(i + 1 for i, (a, b) in enumerate(zip(query, subject, strict=True)) if a != b)
    if len(mismatch) > 3:
        raise ValueError("A local alignment failed reference validation.")
    return OffTargetHit(
        "",
        fields[2],
        title,
        len(query) - len(mismatch),
        len(query),
        0,
        1,
        len(query),
        start + 1 if fields[1] == "+" else start + len(query),
        start + len(query) if fields[1] == "+" else start + 1,
        "plus" if fields[1] == "+" else "minus",
        None,
        None,
        mismatch,
        True,
        query,
        subject,
    )


def _align(
    manager: ReferenceManager,
    folder: Path,
    manifest: dict[str, Any],
    scope: str,
    queries: Mapping[str, str],
    say: ProgressReporter,
) -> tuple[dict[str, Path], frozenset[str]]:
    if manager.shutdown_event.is_set():
        raise ValueError("The app is closing; run this local screening again after restarting.")
    cache = manager.root / "search-cache"
    cache.mkdir(exist_ok=True)
    files = {
        sequence: cache
        / (
            hashlib.sha256(
                json.dumps(
                    {
                        "sequence": sequence,
                        "scope": scope,
                        "reference": manifest["sources"],
                        "tool": manifest["tool_sha256"],
                        "policy": "v3-all-both-strands-v1",
                    },
                    sort_keys=True,
                ).encode()
            ).hexdigest()
            + ".gz"
        )
        for sequence in set(queries.values())
    }
    for path in files.values():
        if path.exists() and path.stat().st_size < 20:
            path.unlink()
    reused = frozenset(sequence for sequence, path in files.items() if path.exists())
    missing = [sequence for sequence, path in files.items() if not path.exists()]
    if not missing:
        say("Reusing completed local alignments", {"phase": "reading", "backend": "local"})
        return files, reused
    with tempfile.TemporaryDirectory(prefix="search-", dir=cache) as temporary:
        work = Path(temporary)
        reads = work / "queries.fa"
        reads.write_text(
            "".join(f">q{i}\n{sequence}\n" for i, sequence in enumerate(missing)), encoding="utf-8"
        )
        say(
            "Searching the local human reference on both strands",
            {"phase": "reading", "backend": "local"},
        )
        # One native index at a time avoids concurrent genome/transcript RAM peaks.
        args = [
            str(folder / manifest["aligner"]),
            "-f",
            "-v",
            "3",
            "-a",
            "--best",
            "-p",
            "2",
            str(folder / scope),
            str(reads),
        ]
        with (
            (work / "alignments.tsv").open("wb") as output,
            (work / "stderr.log").open("wb") as errors,
            subprocess.Popen(args, stdout=output, stderr=errors) as process,
        ):
            started = time.monotonic()
            while process.poll() is None:
                if (
                    manager.shutdown_event.is_set()
                    or time.monotonic() - started > 1800
                    or output.tell() > 2 * 1024**3
                ):
                    process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    if manager.shutdown_event.is_set():
                        raise ValueError(
                            "Local search stopped because the app is closing. "
                            "Run it again after restarting."
                        )
                    raise ValueError(
                        "Local search exceeded its resource budget. No partial result was "
                        "reported; screen fewer candidates or review low-complexity "
                        "sequences."
                    )
                time.sleep(0.1)
            if process.returncode:
                raise ValueError(
                    "Local search failed: "
                    + (work / "stderr.log").read_text(encoding="utf-8", errors="replace")[-1200:]
                )
        with ExitStack() as stack:
            streams = [
                stack.enter_context(gzip.open(work / f"q{i}.gz", "wt", encoding="utf-8"))
                for i in range(len(missing))
            ]
            count = 0
            with (work / "alignments.tsv").open(encoding="utf-8") as alignments:
                for line in alignments:
                    fields = line.split("\t")
                    if not fields[0].startswith("q"):
                        raise ValueError("Unexpected query identifier in local output.")
                    index = int(fields[0][1:])
                    if not 0 <= index < len(missing):
                        raise ValueError("Unexpected query identifier in local output.")
                    streams[index].write(line)
                    count += 1
                    if count > SAFETY_HIT_LIMIT:
                        raise ValueError(
                            "This batch exceeds 250,000 alignments. No truncated result was "
                            "saved. Screen fewer candidates or review repetitive sequences."
                        )
        for i, sequence in enumerate(missing):
            (work / f"q{i}.gz").replace(files[sequence])
    return files, reused


def run_local_screening(
    manager: ReferenceManager,
    run_id: str,
    candidates: Sequence[Mapping[str, Any]],
    scope: str,
    core_length: int,
    progress: ProgressReporter,
    only: Collection[str] | None = None,
    limit: int = 15,
    window: tuple[int, int] | None = None,
    reference_release: str | None = None,
) -> dict[str, Any]:
    """Search locally, count all eligible alignments, and retain bounded UI details."""
    if scope not in ("genomic", "transcript"):
        raise ValueError("Unknown local search scope.")
    queries = candidate_sequences(candidates, limit=None if only is not None else limit, only=only)
    if (
        core_length != 19
        or not queries
        or any(len(sequence) != 19 or set(sequence) - set("ACGT") for sequence in queries.values())
    ):
        raise ValueError("Local screening currently requires unambiguous 19-base DNA cores.")
    progress("Queued for local reference access", {"phase": "reading", "backend": "local"})
    with manager.search_lock:
        folder, manifest = manager.active(reference_release)
        files, reused = _align(manager, folder, manifest, scope, queries, progress)
        progress(
            "Validating alignments and placing hits using local annotations",
            {"phase": "placing", "backend": "local"},
        )
        hits: list[OffTargetHit] = []
        annotations: dict[str, Any] = {}
        with (
            closing(sqlite3.connect(folder / "annotation.sqlite3")) as db,
            (folder / f"{scope}.fa").open("rb") as fasta,
        ):
            metadata: dict[str, tuple[str, int, int]] = {}
            for candidate, query in queries.items():
                with gzip.open(files[query], "rt", encoding="utf-8") as cached:
                    for line in cached:
                        fields = line.rstrip("\r\n").split("\t")
                        accession = fields[2]
                        if accession not in metadata:
                            row = db.execute(
                                "SELECT title,offset,length FROM sequences WHERE scope=? AND "
                                "accession=?",
                                (scope, accession),
                            ).fetchone()
                            if row is None:
                                raise ValueError(
                                    "The aligner returned an unknown reference sequence."
                                )
                            metadata[accession] = row
                        title, offset, length = metadata[accession]
                        if int(fields[3]) < 0 or int(fields[3]) + len(query) > length:
                            raise ValueError(
                                "The aligner returned an out-of-bounds reference position."
                            )
                        validated = parse_alignment(line, query, title, offset, fasta)
                        if set(validated.subject_sequence) - set("ACGT"):
                            continue  # Uncalled bases are outside the search policy.
                        hits.append(replace(validated, candidate_id=candidate))
                        if len(hits) > SAFETY_HIT_LIMIT:
                            raise ValueError(
                                "This batch exceeds the local alignment resource budget. No "
                                "partial result was saved."
                            )
            loci = loci_from_candidates(candidates)
            off, own = split_on_target(hits, loci) if scope == "genomic" else (tuple(hits), ())
            placements: list[PromoterProximity] = []
            if scope == "genomic":
                for hit in off:
                    position = (hit.subject_start + hit.subject_end) // 2
                    nearby = db.execute(
                        "SELECT accession,tss,strand,gene FROM transcripts WHERE contig=? "
                        "AND is_tss_proxy=1 AND tss BETWEEN ? AND ? "
                        "ORDER BY ABS(tss-?),accession LIMIT 1",
                        (hit.accession, position - 10000, position + 10000, position),
                    ).fetchone()
                    label = db.execute(
                        "SELECT label FROM contigs WHERE accession=?", (hit.accession,)
                    ).fetchone()
                    distance = (
                        (position - nearby[1]) * (1 if nearby[2] == "+" else -1) if nearby else None
                    )
                    inside = bool(
                        window and distance is not None and window[0] <= distance <= window[1]
                    )
                    placements.append(
                        PromoterProximity(
                            hit,
                            label[0] if label else hit.accession,
                            position,
                            nearby[3] if nearby else None,
                            nearby[0] if nearby else None,
                            distance,
                            inside,
                            bool(label and label[0] in CHROMOSOMES),
                            "available",
                        )
                    )
            else:
                for accession in {h.accession for h in off}:
                    rows = db.execute(
                        "SELECT DISTINCT c.label,t.gene FROM transcripts t LEFT JOIN contigs c "
                        "ON c.accession=t.contig WHERE t.accession=?",
                        (accession,),
                    ).fetchall()
                    labels = sorted({r[0] for r in rows if r[0] in CHROMOSOMES})
                    annotations[accession] = {
                        "status": "available" if labels else "unavailable",
                        "chromosome": ", ".join(labels) if labels else None,
                        "source": "NCBI RefSeq 110 GFF transcript mapping",
                        "reference_contigs": sorted({r[0] for r in rows if r[0]}),
                        "reason": None
                        if labels
                        else "No primary chromosome mapping in this pinned GFF.",
                        "gene": rows[0][1] if rows else None,
                    }
        parameters = SearchParameters(
            scope,
            manifest["release"],
            species_restriction="database",
            program="bowtie",
            task="full-length-ungapped-v3-all-both-strands",
            reward=None,
            penalty=None,
            gap_costs=None,
            word_size=None,
            expect=None,
            max_hits=0,
            max_alignments_per_candidate=DETAIL_LIMIT,
        )
        record = ScreeningRecord(
            run_id,
            scope,
            parameters,
            off,
            tuple(placements),
            candidates_searched=tuple(queries),
            on_target_hits=own,
            window=window,
        )
        rows = specificity_table(record, core_length)
        retained: list[OffTargetHit] = []
        counts: dict[str, int] = defaultdict(int)
        for hit in sorted(
            off,
            key=lambda h: (
                h.candidate_id,
                -h.identities,
                h.accession,
                h.subject_start,
                h.subject_end,
                h.subject_strand,
            ),
        ):
            if counts[hit.candidate_id] < DETAIL_LIMIT:
                retained.append(hit)
                counts[hit.candidate_id] += 1
        retained_set = set(retained)
        notes: tuple[str, ...] = (
            "Complete enumeration of full-length 19-base alignments with 0-3 mismatches on both"
            " strands of the recorded linear reference records. Gaps, bulges, shorter seed matches,"
            " individual"
            " variation and biological binding are not assessed. Alignments containing uncalled"
            " reference bases are excluded.",
            "Counts include all enumerated alignments after coordinate-based own-locus removal."
            " At most 200 alignment details per candidate are retained; this is a "
            "display/storage limit, not a search-count limit.",
            "Nearest annotated TSS is searched within 10 kb using the same pinned RefSeq GFF. "
            "This annotation window is a reporting choice, not a validated definition of a "
            "promoter. Alternate and patch sequences can represent overlapping biological loci.",
        )
        if scope == "genomic":
            notes += (
                "TSS proxies use the annotated 5-prime boundaries of non-partial mRNA, "
                "lnc_RNA, transcript and primary_transcript models. Mature small-RNA feature "
                "starts are not automatically treated as promoters. These are annotation "
                "proxies, not measured cell-specific transcription starts.",
            )
            unlocated = sum(candidate not in loci for candidate in queries)
            if unlocated:
                notes += (
                    f"Own-locus removal unavailable for {unlocated} candidate(s) without "
                    "genomic coordinates; their counts may include the intended target.",
                )
        for row in rows:
            row["details_retained"] = counts[row["candidate"]]
            row["enumeration_complete"] = True
            row["alignment_cache_hit"] = queries[row["candidate"]] in reused
            if scope == "genomic":
                row["annotation_complete"] = True
                row["annotation_unknown_hits"] = 0
                row["in_run_window_hits"] = sum(
                    p.inside_run_window
                    for p in placements
                    if p.hit.candidate_id == row["candidate"]
                )
        compact = replace(
            record,
            hits=tuple(retained),
            proximities=tuple(p for p in placements if p.hit in retained_set),
            notes=notes,
        )
        return {
            **compact.to_dict(),
            "hit_count": len(off),
            "specificity": rows,
            "backend": "local",
            "software": {"version": __version__, **software_identity()},
            "reference": manifest,
            "gene_annotations": annotations,
            "search_policy": {
                "mismatches": 3,
                "full_length": 19,
                "gaps": False,
                "reference_records": "linear",
                "circular_origin_spanning": False,
                "both_strands": True,
                "enumeration_complete": True,
                "detail_limit": DETAIL_LIMIT,
                "tss_proxy_features": sorted(TSS_PROXY_FEATURES),
                "partial_tss_models_excluded": True,
            },
        }
