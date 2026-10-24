"""Versioned literature observations, separate from design rules and efficacy predictions."""

import csv
import hashlib
import re
from copy import deepcopy
from functools import lru_cache
from importlib import resources
from io import StringIO
from typing import Any

ARCHIVE = "https://zenodo.org/records/21871232"
VERSION = "1.0.0"


@lru_cache(maxsize=1)
def _database() -> tuple[dict[str, list[dict[str, Any]]], dict[str, str]]:
    data = resources.files("sarna.app").joinpath("data")
    records_bytes = data.joinpath("sarna_records_v1.csv").read_bytes()
    sources_bytes = data.joinpath("sarna_sources_v1.csv").read_bytes()
    sources = {
        row["publication_id"]: row
        for row in csv.DictReader(StringIO(sources_bytes.decode("utf-8-sig")))
    }
    genes: dict[str, list[dict[str, Any]]] = {}
    for row in csv.DictReader(StringIO(records_bytes.decode("utf-8-sig"))):
        if row["species_of_origin"] != "Homo sapiens":
            continue
        source = sources[row["publication_id"]]
        doi = source["doi"].strip()
        pmid = source["pmid"].strip()
        patent = source["patent_number"].strip()
        url = (f"https://doi.org/{doi}" if re.fullmatch(r"10\.\d{4,9}/\S+", doi) else
               f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid.isdigit() else
               f"https://patents.google.com/patent/{patent}"
               if re.fullmatch(r"[A-Z]{2}\d+[A-Z]\d?", patent) else "")
        genes.setdefault(row["ensembl_gene_id"], []).append({
            **row, "source": {**source, "url": url},
        })
    return genes, {
        "records_sha256": hashlib.sha256(records_bytes).hexdigest(),
        "sources_sha256": hashlib.sha256(sources_bytes).hexdigest(),
    }


def for_gene(gene_id: str) -> dict[str, Any]:
    """Retrieve human gene-associated reports without treating controls as active saRNAs."""
    database, hashes = _database()
    records = database.get(gene_id.split(".")[0], [])
    grouped: dict[str, dict[str, Any]] = {}
    for record in records:
        sequence_id = record["sequence_id"]
        group = grouped.setdefault(sequence_id, {
            "sequence_id": sequence_id,
            "sense": record["sense_sequence_normalized"],
            "antisense": record["antisense_sequence_normalized"],
            "observations": [],
        })
        group["observations"].append(deepcopy(record))
    return {
        "version": VERSION, "archive": ARCHIVE, "license": "CC BY 4.0", **hashes,
        "gene_id": gene_id, "record_count": len(records),
        "sequences": list(grouped.values()),
        "scope": "Gene-associated literature observations, including possible controls and "
        "mechanistic reagents. Effects retain their source scale; RNA/protein/function, "
        "statistical support and replication require source review. No effect is inferred "
        "from a numeric value. Sequence identity does not establish matching chemistry, "
        "targeting role or efficacy in your cell. Wang criteria and ranking are unchanged.",
    }


def locate(evidence: dict[str, Any], snapshot: dict[str, Any]) -> dict[str, Any]:
    """Match full reported strands or explicitly annotated, complementary dTdT cores.

    Only plain nucleotide strands of 12-100 bases are searched. A match is
    sequence evidence only. Database sense/antisense assignments, treatment
    equivalence and experimental roles are not independently verified by this operation.
    """
    result = deepcopy(evidence)
    promoter = snapshot["promoter"]["sequence"]
    candidates = {
        (candidate["offset"], candidate["target_site_dna"]): candidate["id"]
        for candidate in snapshot["candidates"]
    }
    for group in result["sequences"]:
        queries: dict[tuple[str, str], list[str]] = {}
        observation_ids = [record["record_id"] for record in group["observations"]]
        for strand in (group["sense"], group["antisense"]):
            if re.fullmatch(r"[ACGTU]{12,100}", strand):
                dna = strand.replace("U", "T")
                complement = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                for query in (dna, complement):
                    queries[query, "Full reported strand"] = observation_ids
        sense, antisense = (group[key].replace("U", "T") for key in ("sense", "antisense"))
        if (re.fullmatch(r"[ACGT]{14,100}TT", sense)
                and re.fullmatch(r"[ACGT]{14,100}TT", antisense)
                and sense[:-2] == antisense[:-2].translate(str.maketrans("ACGT", "TGCA"))[::-1]):
            annotated = [record["record_id"] for record in group["observations"]
                         if re.search(r"(?:^|,)\s*dTdT\s*overhang\s*(?:,|$)",
                                      record["modifications_reported"], re.I)]
            if annotated:
                for query in (sense[:-2], antisense[:-2]):
                    queries[query, "Paired core; source annotates dTdT overhang"] = annotated
        matches = []
        for (query, basis), record_ids in sorted(queries.items()):
            offset = promoter.find(query)
            while offset >= 0:
                matches.append({
                    "offset": offset, "end": offset + len(query), "target": query,
                    "candidate_id": candidates.get((offset, query)),
                    "basis": basis, "record_ids": record_ids,
                })
                offset = promoter.find(query, offset + 1)
        group["matches"] = matches
    return result
