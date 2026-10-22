"""Explicit promoter references and exact-coordinate transcript groups for the workbench."""

from typing import Any

from sarna.app.models import CustomReference
from sarna.genome.anchors import PromoterAnchor, ReferenceKind, load_anchors
from sarna.genome.coordinates import GenomicInterval
from sarna.genome.ensembl import Gene, Transcript
from sarna.rules.provenance import Provenance


def anchor_record(anchor: PromoterAnchor) -> dict[str, Any]:
    """Add a readable identifier without changing the source or reference kind."""
    return {**anchor.to_dict(), "label": f"{anchor.gene} · {anchor.id}"}


def transcript_record(transcript: Transcript) -> dict[str, Any]:
    """Expose transcript identity and its annotated reference without arbitrary selection."""
    return {
        "id": transcript.ensembl_id, "biotype": transcript.biotype,
        "mane": transcript.is_mane_select, "canonical": transcript.is_canonical,
        "refseq": transcript.refseq_match, "contig": transcript.contig,
        "strand": str(transcript.strand), "tss": transcript.tss_coordinate,
    }


def reference_priority(transcript: Transcript) -> tuple[int, int, int, str]:
    """Order transcripts so the reference most researchers want is offered first.

    MANE Select leads, because it is the one reference this project either picks or refuses to
    pick (ADR 0012); Ensembl's canonical choice follows, then protein-coding transcripts, then
    the rest by identifier so the order is stable across lookups instead of following whatever
    order the API happened to return.

    This is **presentation order only**. Every annotated transcript stays selectable, and being
    first does not make an annotated coordinate a measured start site: chapter 2's own SCN1A
    hits came from a distal region with no Ensembl TSS within 5.5 kb.
    """
    return (
        0 if transcript.is_mane_select else 1,
        0 if transcript.is_canonical else 1,
        0 if transcript.is_protein_coding else 1,
        transcript.ensembl_id,
    )


def tss_groups(transcripts: tuple[Transcript, ...]) -> list[dict[str, Any]]:
    """Group identical contig/strand/TSS coordinates; nearby starts remain distinct.

    Groups come back in the order their first transcript appears, and each group lists its
    transcripts in the order given, so passing transcripts already sorted by
    :func:`reference_priority` puts the MANE Select group (and the MANE transcript inside it)
    at the top without this function needing to know what a preferred reference is.

    These groups describe annotation coordinates, not experimentally verified promoters.
    Transcript identities remain available for deliberate reference selection.
    """
    groups: dict[str, dict[str, Any]] = {}
    for transcript in transcripts:
        group_id = f"tss:{transcript.contig}:{transcript.strand}:{transcript.tss_coordinate}"
        group = groups.setdefault(group_id, {
            "id": group_id, "contig": transcript.contig, "strand": str(transcript.strand),
            "tss": transcript.tss_coordinate, "transcript_ids": [], "mane": False,
        })
        group["transcript_ids"].append(transcript.ensembl_id)
        group["mane"] = group["mane"] or transcript.is_mane_select
    return list(groups.values())


def gene_overview(gene: Gene, transcripts: tuple[Transcript, ...]) -> dict[str, Any]:
    """Freeze complete reference choices and their source beside every gene result."""
    from sarna.app.sequence_import import references

    mane = [transcript for transcript in transcripts if transcript.is_mane_select]
    ordered = tuple(sorted(transcripts, key=reference_priority))
    return {
        "symbol": gene.symbol, "id": gene.ensembl_id, "description": gene.description,
        "assembly": str(gene.interval.assembly), "strand": str(gene.strand),
        "contig": gene.contig, "start": gene.interval.start, "end": gene.interval.end,
        "transcripts": [transcript_record(transcript) for transcript in ordered],
        "tss_groups": tss_groups(ordered),
        "sarna_references": references(gene.symbol),
        "default_transcript": mane[0].ensembl_id if len(mane) == 1 else None,
        "anchors": [anchor_record(anchor) for anchor in load_anchors().for_gene(gene.symbol)],
        "source": "Ensembl REST (cached responses reused when present)",
    }


def custom_anchor(gene: Gene, reference: CustomReference) -> PromoterAnchor:
    """Convert a supplied 1-based base to an inferred anchor on the confirmed gene's strand.

    A position may lie outside the annotated gene span. Its association with the gene is
    supplied by the researcher, and never inferred from overlap with the annotation.
    """
    coordinate = reference.tss_position - 1
    return PromoterAnchor(
        id=reference.id, gene=gene.symbol,
        interval=GenomicInterval(
            assembly=gene.interval.assembly, contig=gene.contig,
            start=coordinate, end=coordinate + 1, strand=gene.strand,
        ),
        reference_point=coordinate, reference_kind=ReferenceKind.INFERRED,
        provenance=Provenance.citation(
            reference.source,
            f"Researcher-supplied {reference.name!r}, genomic position "
            f"{reference.tss_position} (1-based); source and position not independently verified",
        ),
        note=reference.note,
    )


def custom_anchor_record(anchor: PromoterAnchor, reference: CustomReference) -> dict[str, Any]:
    """Expose supplied metadata and its uncertainty without implying source verification."""
    return {
        **anchor_record(anchor), "label": reference.name, "name": reference.name,
        "source": reference.source, "researcher_supplied": True,
        "tss_position": reference.tss_position,
        "caveat": "Researcher-supplied reference. Its coordinate, cited source and association "
        "with this gene have not been independently verified. Offsets carry that uncertainty.",
    }
