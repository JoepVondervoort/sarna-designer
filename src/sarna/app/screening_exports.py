"""Export recorded sequence evidence without merging references or inventing zero counts."""

import json
from collections.abc import Mapping, Sequence
from typing import Any

from sarna.app.workbook import Sheet

SUMMARY_COLUMNS = (
    "run_id", "candidate", "design_rank", "preferences_met", "preferences_evaluated",
    "scope", "status", "search_id", "searched_at", "backend", "reference",
    "reported_alignments", "matches_19_of_19", "matches_18_of_19", "matches_17_of_19",
    "matches_16_of_19", "details_retained", "enumeration_complete", "alignments_capped",
    "own_locus_alignments_removed", "count_unit", "limitations",
)
EXPORT_NOTES = (
    "Summary: latest recorded evidence per candidate, scope, backend and reference. "
    "Separate searches and references are never pooled.",
    "Counts describe reference alignments, not independent genes, binding or biological effects.",
    "Missing, running or failed searches have blank counts, not zero counts. Historical "
    "completed evidence remains in the search records and retained alignment sheets.",
    "Local counts enumerate the declared search policy; retained details are a subset "
    "(normally up to 200 per candidate). Gaps, variants and seed effects are not assessed locally.",
    "Remote BLAST reporting limits may make counts lower bounds. Read search methods and notes.",
)


def record_reference(record: Mapping[str, Any]) -> tuple[str, str]:
    """Identify the method and reference exactly as recorded, including pending jobs."""
    result = record.get("result") or {}
    backend = result.get("backend") or (record.get("progress") or {}).get("backend") or "ncbi"
    reference = ((result.get("reference") or {}).get("release")
                 or result.get("reference_release")
                 or (result.get("parameters") or {}).get("database") or "Not recorded")
    return str(backend), str(reference)


def screening_summary(
    run_id: str, candidates: Sequence[Mapping[str, Any]], records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Return latest evidence per recorded reference, plus explicit unscreened scopes."""
    latest: dict[tuple[str, str, str, str], tuple[Mapping[str, Any], Mapping[str, Any]]] = {}
    for record in sorted(records, key=lambda r: str(r.get("created", ""))):
        result = record.get("result") or {}
        rows = {r["candidate"]: r for r in result.get("specificity", [])}
        backend, reference = record_reference(record)
        for candidate in set(result.get("candidates_searched", [])) | rows.keys():
            key = (candidate, str(record.get("scope")), backend, reference)
            latest[key] = record, rows.get(candidate, {})
    exported = []
    for candidate in sorted(candidates, key=lambda c: (c.get("rank") is None,
                                                     c.get("rank") or 0, c["id"])):
        rating = candidate.get("rating") or {}
        for scope in ("transcript", "genomic"):
            matches = [(key, pair) for key, pair in latest.items()
                       if key[:2] == (candidate["id"], scope)]
            if not matches:
                matches = [((candidate["id"], scope, "", ""), ({}, {}))]
            for key, (record, row) in matches:
                complete = record.get("status") == "complete" and bool(row)
                counts = row.get("identity_counts") or {} if complete else {}
                result = record.get("result") or {}
                exported.append(dict(zip(SUMMARY_COLUMNS, (
                    run_id, candidate["id"], candidate.get("rank"), rating.get("met"),
                    rating.get("evaluated"), scope,
                    record.get("status", "not_screened") if not complete else (
                        "complete" if row.get("enumeration_complete")
                        and not row.get("alignments_capped") else "reported_counts_only"),
                    record.get("id"), result.get("searched_at") or record.get("created"),
                    key[2], key[3], row.get("hit_count") if complete else None,
                    *(counts.get(str(n)) for n in (19, 18, 17, 16)),
                    row.get("details_retained", sum(
                        h.get("candidate") == candidate["id"] for h in result.get("hits", [])
                    )) if complete else None,
                    row.get("enumeration_complete") if complete else None,
                    row.get("alignments_capped") if complete else None,
                    row.get("on_target_hits") if complete else None,
                    "reference alignments, including isoforms and alternate records",
                    " ".join(EXPORT_NOTES[1:]),
                ), strict=True)))
    return exported


def screening_sheets(
    run_id: str, candidates: Sequence[Mapping[str, Any]], records: Sequence[Mapping[str, Any]],
) -> list[Sheet]:
    """Keep summary counts, retained details and every relevant search record separate."""
    summary = screening_summary(run_id, candidates, records)
    ids = {c["id"] for c in candidates}
    searches, alignments = [], []
    for record in records:
        result = record.get("result") or {}
        searched = set(result.get("candidates_searched", [])) | {
            row["candidate"] for row in result.get("specificity", [])
        }
        if searched and not searched & ids:
            continue
        backend, reference = record_reference(record)
        searches.append((
            record.get("id"), record.get("scope"), record.get("status"), record.get("created"),
            result.get("searched_at"), backend, reference, ", ".join(sorted(searched & ids)),
            record.get("error"), json.dumps(result.get("parameters") or {}, sort_keys=True),
            json.dumps(result.get("reference") or {}, sort_keys=True),
            json.dumps(result.get("software") or {}, sort_keys=True),
            json.dumps(result.get("search_policy") or {}, sort_keys=True),
            " ".join(result.get("notes") or []),
        ))
        if record.get("status") != "complete":
            continue
        placements = {
            (p["hit"]["candidate"], p["hit"].get("accession"), p["hit"].get("subject_start"),
             p["hit"].get("subject_end"), p["hit"].get("subject_strand")): p
            for p in result.get("proximities", []) if p.get("hit")
        }
        for hit in result.get("hits", []):
            if hit.get("candidate") not in ids:
                continue
            place = placements.get((hit["candidate"], hit.get("accession"),
                                    hit.get("subject_start"), hit.get("subject_end"),
                                    hit.get("subject_strand")), {})
            alignments.append((
                record.get("id"), record.get("scope"), backend, reference, hit["candidate"],
                hit.get("accession"), hit.get("title"), hit.get("identities"),
                hit.get("alignment_length"), hit.get("gaps"), hit.get("subject_start"),
                hit.get("subject_end"), hit.get("subject_strand"), hit.get("query_start"),
                hit.get("query_end"), ", ".join(map(str, hit.get("mismatch_positions", []))),
                hit.get("query_sequence"), hit.get("subject_sequence"),
                place.get("nearest_gene"), place.get("distance_to_tss"),
                "Retained detail; may be a subset of the reported alignments",
            ))
    return [
        Sheet("Off-target summary", SUMMARY_COLUMNS,
              tuple(tuple(row[key] for key in SUMMARY_COLUMNS) for row in summary)),
        Sheet("Retained alignments", (
            "search_id", "scope", "backend", "reference", "candidate", "accession", "title",
            "matching_bases", "alignment_length", "gaps", "subject_start", "subject_end",
            "subject_strand", "query_start", "query_end", "mismatch_positions_1_based",
            "aligned_query", "aligned_subject", "nearest_annotated_gene", "distance_to_tss",
            "detail_coverage",
        ), tuple(alignments)),
        Sheet("Screening searches", (
            "search_id", "scope", "status", "created", "searched_at", "backend", "reference",
            "candidates_in_selection", "error", "parameters_json", "reference_json",
            "software_json", "search_policy_json", "notes",
        ), tuple(searches)),
        Sheet("Screening readme", ("interpretation",), tuple((note,) for note in EXPORT_NOTES)),
    ]
