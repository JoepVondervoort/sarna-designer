"""Explicit strand normalization, bounded file imports and promoter matching."""

import csv
import io
import re
from datetime import UTC, datetime
from pathlib import PurePosixPath
from typing import Any
from xml.etree import ElementTree as ET
from zipfile import BadZipFile, ZipFile

from sarna import __version__
from sarna.design.pipeline import DEFAULT_COMPARISON_TABLES
from sarna.reproducibility import json_hash, software_identity
from sarna.sequence.nucleic import NucleicSequence
from sarna.thermo.guide import compare_asymmetry
from sarna.validation.dataset import load_chapter4_dataset, load_chapter4_locations

MAX_RECORDS = 200
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def core(value: str) -> tuple[str, str]:
    """Keep overhang chemistry separate; never silently trim unmarked terminal bases."""
    value = re.sub(r"\s+", "", value).upper()
    overhang = ""
    for suffix in ("[DT][DT]", "[UU]"):
        if value.endswith(suffix):
            value, overhang = value[: -len(suffix)], suffix
            break
    if not re.fullmatch("[ACGTU]{4,100}", value):
        raise ValueError("Use 4-100 A/C/G/T/U bases; mark a 3-prime overhang as [dT][dT] or [UU].")
    if "T" in value and "U" in value:
        raise ValueError("Mixed T/U in a core: specify a DNA or RNA sequence.")
    return value.replace("U", "T"), overhang


def reverse(sequence: str) -> str:
    """Return the DNA reverse complement, written 5-prime to 3-prime."""
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def normalize(row: dict[str, Any]) -> dict[str, Any]:
    """Validate supplied strand roles and a fully paired core without inventing overhangs."""
    role = str(row.get("role") or "unknown").lower()
    if role not in ("unknown", "target", "sense", "antisense"):
        raise ValueError("Role must be unknown, target, sense or antisense.")
    first = str(row.get("sense") or row.get("sequence") or "").strip()
    if row.get("sense"):
        role = "sense"
    if not first and row.get("antisense"):
        first, role = str(row["antisense"]), "antisense"
    direction = str(row.get("direction") or "5to3")
    if direction not in ("5to3", "3to5"):
        raise ValueError("Direction must be 5to3 or 3to5.")
    dna, overhang = core(first)
    if direction == "3to5":
        if overhang:
            raise ValueError("Enter overhang-bearing strands in 5to3 direction.")
        dna = dna[::-1]
    partner = str(row.get("antisense") or "").strip() if role == "sense" else ""
    if partner:
        paired, partner_overhang = core(partner)
        if direction == "3to5":
            paired = paired[::-1]
        if paired != reverse(dna):
            raise ValueError(
                "The supplied cores are not reverse complements. "
                "Check direction, mismatches and overhang notation."
            )
    else:
        partner_overhang = ""
    return {
        "name": str(row.get("name") or "Imported saRNA")[:150],
        "sequence": first,
        "antisense": partner,
        "role": role,
        "direction": direction,
        "dna": dna,
        "overhang": overhang,
        "partner_overhang": partner_overhang,
        "source": str(row.get("source") or "")[:1000],
        "partner_inferred": not bool(partner),
    }


def input_identity(row: dict[str, Any]) -> str:
    """Identify supplied inputs without predictions, matches or browser-only fields."""
    normalized = normalize(row)
    return json_hash(
        {
            key: normalized[key]
            for key in (
                "name",
                "sequence",
                "antisense",
                "role",
                "direction",
                "source",
            )
        }
    )


def file_rows(filename: str, content: bytes) -> list[dict[str, str]]:
    """Read text/CSV or first Excel sheet with strict archive and row limits."""
    if len(content) > 2_000_000:
        raise ValueError("Upload at most 2 MB.")
    if filename.lower().endswith(".xlsx"):
        try:
            with ZipFile(io.BytesIO(content)) as archive:
                if sum(item.file_size for item in archive.infolist()) > 8_000_000:
                    raise ValueError("Expanded workbook exceeds 8 MB.")
                strings = []
                if "xl/sharedStrings.xml" in archive.namelist():
                    strings = [
                        "".join(e.itertext())
                        for e in ET.fromstring(archive.read("xl/sharedStrings.xml"))
                    ]
                workbook = ET.fromstring(archive.read("xl/workbook.xml"))
                sheets = workbook.find(NS + "sheets")
                if sheets is None or not len(sheets):
                    raise ValueError("Workbook has no worksheets.")
                relation_id = sheets[0].get(
                    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
                )
                relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
                relation = next((item for item in relations if item.get("Id") == relation_id), None)
                if relation is None or relation.get("TargetMode") == "External":
                    raise ValueError("The first worksheet must be stored in this workbook.")
                target = relation.get("Target", "")
                sheet_path = str(
                    PurePosixPath(target.lstrip("/"))
                    if target.startswith("/")
                    else PurePosixPath("xl") / target
                )
                if ".." in PurePosixPath(sheet_path).parts:
                    raise ValueError("Invalid worksheet path.")
                sheet = ET.fromstring(archive.read(sheet_path))
                table = []
                for row in sheet.iter(NS + "row"):
                    values: dict[int, str] = {}
                    column = 0
                    for cell in row:
                        reference = cell.get("r")
                        if reference is None:
                            column += 1
                        else:
                            if not re.fullmatch(r"[A-Z]{1,3}[1-9][0-9]{0,6}", reference):
                                raise ValueError("Invalid Excel cell reference.")
                            letters = re.sub("[0-9]", "", reference)
                            column = 0
                            for letter in letters:
                                column = column * 26 + ord(letter) - 64
                        if column > 16384:
                            raise ValueError("Excel column exceeds the worksheet limit.")
                        value = cell.findtext(NS + "v", "")
                        if cell.get("t") == "s":
                            value = strings[int(value)]
                        elif cell.get("t") == "inlineStr":
                            inline = cell.find(NS + "is")
                            value = "".join(inline.itertext()) if inline is not None else ""
                        values[column - 1] = value
                    table.append([values.get(i, "") for i in range(max(values, default=-1) + 1)])
                    if len(table) > MAX_RECORDS + 1:
                        raise ValueError("Import at most 200 sequences at a time.")
        except (BadZipFile, KeyError, ET.ParseError, IndexError) as exc:
            raise ValueError("Could not read this .xlsx workbook.") from exc
    else:
        if not filename.lower().endswith((".txt", ".csv")):
            raise ValueError("Choose .txt, .csv or .xlsx.")
        text = content.decode("utf-8-sig")
        table = list(csv.reader(io.StringIO(text), delimiter="\t" if "\t" in text else ","))
    table = [row for row in table if any(cell.strip() for cell in row)]
    if not table:
        raise ValueError("No sequences found.")
    header = [cell.strip().lower() for cell in table[0]]
    if any(key in header for key in ("sequence", "sense", "antisense")):
        rows = [dict(zip(header, row, strict=False)) for row in table[1:]]
    else:
        rows = [
            {"sequence": row[0]} if len(row) == 1 else {"name": row[0], "sequence": row[1]}
            for row in table
        ]
    if len(rows) > MAX_RECORDS:
        raise ValueError("Import at most 200 sequences at a time.")
    return rows


def references(symbol: str) -> list[dict[str, Any]]:
    """Expose curated cores and contextual observations without changing validation labels."""
    aliases = {"P21": "CDKN1A", "P53": "TP53", "VEGF": "VEGFA"}
    symbol = aliases.get(symbol.upper(), symbol.upper())
    names = {
        "SCN1A": ("SCN1A-10", "SCN1A-7"),
        "CDKN1A": ("P21-322",),
        "VEGFA": ("VEGF",),
        "TP53": ("P53",),
        "VHL": ("VHL-821",),
        "CEBPA": ("CEBPA",),
        "KLF4": ("KLF4",),
    }.get(symbol, ())
    dataset = load_chapter4_dataset()
    rows = []
    locations = load_chapter4_locations()
    for name in names:
        record = dataset.by_name(name)
        evidence = record.evidence
        if name == "SCN1A-7":
            evidence += (
                "; researcher reported downregulation: HeLa and AF22 by qPCR; "
                "AF22 NES also RNA-seq. Dataset classification: not activating."
            )
        rows.append(
            normalize(
                {
                    "name": name,
                    "sense": record.target.bases + ("[dT][dT]" if name.startswith("SCN1A") else ""),
                    "source": evidence,
                    "role": "sense",
                }
            )
            | {
                "reference": True,
                "known_location": {
                    "assembly": str(locations[name].target.assembly),
                    "contig": locations[name].target.contig,
                    "start": locations[name].target.start,
                    "end": locations[name].target.end,
                    "strand": str(locations[name].target.strand),
                    "source": "Repository sequence-mapped reference coordinates",
                    "anchor_id": "scn1a-long-2008" if name.startswith("SCN1A") else None,
                    "nearest_transcript": locations[name].nearest_transcript,
                },
            }
        )
    return rows


def compare(
    snapshot: dict[str, Any],
    rows: list[dict[str, Any]],
    *,
    refresh_models: bool = False,
) -> list[dict[str, Any]]:
    """Locate every exact match in transcript-oriented DNA and reuse saved rule verdicts."""
    promoter = snapshot["promoter"]["sequence"]
    candidates = {(c["offset"], c["target_site_dna"]): c for c in snapshot["candidates"]}
    output = []
    provenance = (
        {
            "schema_version": 1,
            "analyzed_at": datetime.now(UTC).isoformat(),
            "engine_version": __version__,
            "package_sha256": software_identity()["package_sha256"],
            "configuration_sha256": json_hash(snapshot.get("configuration")),
            "scope": "Paired cores; saved model assessments and inputs; overhangs excluded.",
        }
        if rows
        else {}
    )
    for row in rows:
        row = normalize(row) | {
            "reference": bool(row.get("reference")),
            "known_location": row.get("known_location"),
        }
        location = row["known_location"]
        window = snapshot["promoter"]
        outside = bool(
            location
            and window.get("assembly") == location["assembly"]
            and window.get("contig") == location["contig"]
            and (location["end"] <= window["start"] or location["start"] >= window["end"])
        )
        dna = row["dna"]
        targets = [reverse(dna)] if row["role"] == "antisense" else [dna]
        if row["role"] == "unknown" and reverse(dna) != dna:
            targets.append(reverse(dna))
        matches: list[dict[str, Any]] = []
        for target in targets:
            offset = promoter.find(target)
            while offset >= 0:
                if len(matches) >= 500:
                    raise ValueError("More than 500 matches: use a longer core or smaller region.")
                candidate = candidates.get((offset, target))
                matches.append(
                    {
                        "offset": offset,
                        "end": offset + len(target),
                        "input_orientation": "same" if target == dna else "reverse complement",
                        "target": target,
                        "candidate": candidate,
                        "sense": target.replace("T", "U"),
                        "antisense": reverse(target).replace("T", "U"),
                        "thermodynamics": (
                            candidate["guide_prediction"]["assessments"]
                            if candidate and not refresh_models
                            else [
                                a.to_dict()
                                for a in compare_asymmetry(
                                    NucleicSequence.parse(target),
                                    DEFAULT_COMPARISON_TABLES,
                                    end_length=(
                                        candidate["guide_prediction"].get("end_length_nt", 4)
                                        if candidate
                                        else 4
                                    ),
                                ).assessments
                            ]
                        ),
                    }
                )
                matches[-1]["strand_analysis"] = strand_analysis(
                    target,
                    window.get("strand"),
                    matches[-1]["thermodynamics"],
                )
                offset = promoter.find(target, offset + 1)
        analysis = None
        if matches:
            analysis = matches[0]["strand_analysis"]
        elif row["role"] != "unknown":
            target = targets[0]
            analysis = strand_analysis(
                target,
                location.get("strand") if location else None,
                [
                    assessment.to_dict()
                    for assessment in compare_asymmetry(
                        NucleicSequence.parse(target),
                        DEFAULT_COMPARISON_TABLES,
                    ).assessments
                ],
            )
        output.append(
            row
            | {
                "matches": matches,
                "strand_analysis": analysis,
                "analysis_provenance": provenance,
                "status": (
                    "matched"
                    if matches
                    else "Outside current promoter window"
                    if outside
                    else "Not found in current promoter window"
                ),
            }
        )
    return output


def strand_analysis(
    target: str,
    genomic_strand: str | None,
    assessments: list[dict[str, Any]],
) -> dict[str, Any]:
    """Describe base composition and complementarity without asserting physical binding."""
    sense_dna = target
    opposite_dna = reverse(target)
    opposite_strand = {"+": "-", "-": "+"}.get(genomic_strand or "")
    strands = []
    for name, rna, paired_dna, paired_strand in (
        ("sense", sense_dna.replace("T", "U"), opposite_dna, opposite_strand),
        ("antisense", opposite_dna.replace("T", "U"), sense_dna, genomic_strand),
    ):
        strands.append(
            {
                "name": name,
                "rna_5to3": rna,
                "length": len(rna),
                "counts": {base: rna.count(base) for base in "ACGU"},
                "gc_percent": round(100 * (rna.count("G") + rna.count("C")) / len(rna), 1),
                "complementary_dna_5to3": paired_dna,
                "complementary_dna_aligned_3to5": paired_dna[::-1],
                "complementary_genomic_strand": paired_strand,
            }
        )
    return {
        "strands": strands,
        "models": [
            {
                **assessment,
                "predicted_passenger": (
                    {"sense": "antisense", "antisense": "sense"}.get(assessment["predicted_guide"])
                ),
            }
            for assessment in assessments
        ],
        "scope": "Sequence complementarity only; actual DNA or promoter-associated RNA binding "
        "and strand loading are not established by these calculations.",
    }
