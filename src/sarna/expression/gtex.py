"""GTEx Portal client: per-tissue median expression for every transcript of a gene.

GTEx quantifies bulk RNA-seq from post-mortem donor tissue. For this project it answers one
question the annotation cannot (*which* transcript of a gene a tissue actually makes) and
pointedly fails to answer another: it contains no cell lines, so AF22 and HeLa are absent and
always will be.

**Two identifiers, not one.** GTEx keys everything on a *versioned* GENCODE identifier, and
which version depends on the dataset: ``gtex_v8`` was built on GENCODE v26 and ``gtex_v10`` on
v39, so SCN1A is ``ENSG00000144285.17`` in one and ``ENSG00000144285.23`` in the other.
Querying with the version Ensembl reports today returns an empty result and no error, which
reads exactly like "this gene is not expressed". Every lookup therefore resolves the versioned
identifier through ``reference/gene`` first, against the release the dataset declares.

**An empty result is not a zero.** A gene GTEx has never heard of comes back as HTTP 200 with
an empty ``data`` list. That is raised as :class:`GeneNotInGtexError` rather than returned as a
profile full of zeros.

**The annotations disagree in both directions.** For SCN1A, six transcripts Ensembl serves
today are absent from GENCODE v39, and one transcript GENCODE v39 quantifies
(``ENST00000375405``) has since been retired from Ensembl. Neither direction is an error; both
are recorded rather than smoothed over.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import httpx

from sarna import settings
from sarna.expression.models import (
    ExpressionStatus,
    GeneExpression,
    TissueMedian,
    TranscriptExpression,
)
from sarna.rules.provenance import Provenance
from sarna.sources.json_http import (
    BACKOFF_SECONDS,
    DEFAULT_TIMEOUT_SECONDS,
    JsonHttpCache,
    JsonObject,
    SourceError,
)

BASE_URL = "https://gtexportal.org"

DEFAULT_DATASET = "gtex_v10"
"""The newest GTEx release. Chosen over v8 because v8's GENCODE v26 predates the MANE Select
transcripts this app offers first, so v8 cannot quantify them at all."""

DATASET_ANNOTATION: Mapping[str, str] = {
    "gtex_v8": "v26",
    "gtex_v10": "v39",
}
"""Which GENCODE release each dataset was quantified against.

Stated here because the API will not infer it: ``medianTranscriptExpression`` takes a versioned
identifier and silently returns nothing for a version it does not hold. Verified against the
portal on 8 October 2026, ``reference/gene`` with ``gencodeVersion=v39`` returns
``ENSG00000144285.23``, which ``gtex_v10`` then quantifies.
"""

GENOME_BUILD = "GRCh38/hg38"
"""Both supported datasets are on GRCh38. Sent explicitly rather than relying on a default,
for the same reason :class:`sarna.genome.ensembl.EnsemblClient` names its assembly (ADR 0004)."""

PAGE_SIZE = 250
"""Rows per request. The portal's own default; a gene is 25 transcripts x 54 tissues."""

MAX_PAGES = 40
"""A guard against a paging loop that never terminates. 40 pages is 10,000 rows, far above the
1,350 a gene like SCN1A produces."""

CELL_LINE_CAVEAT = (
    "GTEx contains no cell lines. AF22, HeLa and iPSC-derived neurons are not in this dataset, "
    "so it holds no value for them at all."
)
"""Carried by the profile itself, not only by the page that renders it, so a methods export or
a script reporting these numbers states it too."""

BULK_TISSUE_CAVEAT = (
    "Every value is bulk RNA-seq of whole post-mortem donor tissue, so it is the average of "
    "many cell types rather than a measurement in any one of them. A transcript used by a "
    "small population of cells can read as near zero in the tissue that contains them."
)

MEDIAN_CAVEAT = (
    "Each number is a median across donors and hides the spread. The sample groups also differ "
    "74-fold in how many donors they have, from 11 to 818, so two medians shown side by side "
    "are not equally certain however similar they look."
)

ISOFORM_CAVEAT = (
    "Transcript-level values are estimates apportioned between overlapping isoforms, which is "
    "less certain than a gene-level total, especially for transcripts that differ only at "
    "their ends."
)


class GtexError(SourceError):
    """Base class for failures talking to the GTEx Portal."""


class GeneNotInGtexError(GtexError):
    """The gene has no record in the requested dataset's annotation release."""


class UnknownGtexDatasetError(GtexError):
    """A dataset was requested whose annotation release this module does not know."""


class GtexResponseError(GtexError):
    """A response was reachable and parseable but not shaped as documented."""


class GtexClient:
    """A caching GTEx Portal client, pinned to one dataset.

    Attributes:
        dataset_id: The GTEx dataset every number from this client belongs to.
        annotation_release: The GENCODE release that dataset was quantified against.
    """

    __slots__ = ("_http", "annotation_release", "dataset_id")

    def __init__(
        self,
        dataset_id: str = DEFAULT_DATASET,
        cache_dir: Path | None = None,
        timeout: float = DEFAULT_TIMEOUT_SECONDS,
        client: httpx.Client | None = None,
        backoff_seconds: float = BACKOFF_SECONDS,
        offline: bool | None = None,
    ) -> None:
        """Build a client for one GTEx dataset.

        Args:
            dataset_id: Which release to query, e.g. ``gtex_v10``. Never a default inside a
                request: the dataset decides which GENCODE release the identifiers belong to.
            cache_dir: Directory for cached responses.
            timeout: Per-request timeout in seconds.
            client: An ``httpx.Client`` to use instead of creating one, for testing.
            backoff_seconds: Base delay between retries, doubled each attempt.
            offline: Whether to refuse external requests. ``None`` reads the environment.

        Raises:
            UnknownGtexDatasetError: If the dataset's annotation release is not known here.
                Guessing it would mean querying with identifiers from the wrong release and
                reading the empty result as absence of expression.
        """
        if dataset_id not in DATASET_ANNOTATION:
            known = ", ".join(sorted(DATASET_ANNOTATION))
            raise UnknownGtexDatasetError(
                f"{dataset_id!r} is not a dataset this module knows the annotation release "
                f"for; known datasets are {known}. Add it to DATASET_ANNOTATION with the "
                f"GENCODE release it was quantified against, verified against the portal."
            )
        self.dataset_id = dataset_id
        self.annotation_release = DATASET_ANNOTATION[dataset_id]
        self._http = JsonHttpCache(
            name="gtex", base_url=BASE_URL, cache_dir=cache_dir, timeout=timeout,
            client=client, backoff_seconds=backoff_seconds, offline=offline,
        )

    @property
    def offline(self) -> bool:
        """Return whether this client refuses external requests."""
        return self._http.offline

    def release_gene_id(self, gene_id: str) -> str:
        """Return the versioned GENCODE identifier this dataset uses for a gene.

        Args:
            gene_id: An Ensembl gene identifier, versioned or not. Any version supplied is
                discarded, because the caller's version comes from Ensembl today and this
                dataset's comes from an older release.

        Returns:
            The versioned identifier, e.g. ``ENSG00000144285.23``.

        Raises:
            GeneNotInGtexError: If the release holds no such gene.
            GtexResponseError: If the response is shaped unexpectedly, or answers about a
                different gene or a different annotation release than the one asked for.
        """
        stem = gene_id.split(".")[0].strip()
        if not stem:
            raise GtexResponseError("a gene identifier is required")
        payload = self._http.get(
            f"/api/v2/reference/gene?geneId={stem}"
            f"&gencodeVersion={self.annotation_release}"
            f"&genomeBuild={GENOME_BUILD.replace('/', '%2F')}"
        )
        records = _rows(payload, "reference/gene")
        matching = [
            record for record in records
            if str(record.get("gencodeId", "")).split(".")[0] == stem
        ]
        if not matching:
            raise GeneNotInGtexError(
                f"{stem} has no record in GENCODE {self.annotation_release}, the release "
                f"{self.dataset_id} was quantified against. This is not a measurement of zero: "
                f"nothing about this gene was quantified in this dataset."
            )
        record = matching[0]
        reported = str(record.get("gencodeVersion", ""))
        if reported and reported != self.annotation_release:
            raise GtexResponseError(
                f"asked GTEx for {stem} in GENCODE {self.annotation_release} but it answered "
                f"for {reported}; identifiers from two releases must not be mixed"
            )
        return str(record["gencodeId"])

    def tissues(self) -> dict[str, TissueMedian]:
        """Return one blank :class:`TissueMedian` per sample group, carrying its metadata.

        The expression endpoint reports only a tissue's identifier, so the readable name, the
        UBERON term and the donor count come from here and are joined on.

        Returns:
            Keyed by tissue identifier, e.g. ``Brain_Cortex``. The ``median`` of each is zero
            and meaningless; only the metadata fields are populated.
        """
        payload = self._http.get(
            f"/api/v2/dataset/tissueSiteDetail?datasetId={self.dataset_id}&itemsPerPage=100"
        )
        catalogue: dict[str, TissueMedian] = {}
        for record in _rows(payload, "dataset/tissueSiteDetail"):
            tissue_id = str(record.get("tissueSiteDetailId", ""))
            if not tissue_id:
                continue
            catalogue[tissue_id] = TissueMedian(
                tissue_id=tissue_id,
                tissue_name=str(record.get("tissueSiteDetail") or tissue_id),
                median=0.0,
                unit="",
                ontology_id=str(record.get("ontologyId") or ""),
                sample_count=_sample_count(record),
            )
        if not catalogue:
            raise GtexResponseError(
                f"{self.dataset_id} reported no tissues, so no median could be labelled"
            )
        return catalogue

    def expression_for(
        self, gene_id: str, transcript_ids: Sequence[str] = ()
    ) -> GeneExpression:
        """Return per-tissue medians for every transcript of one gene.

        Args:
            gene_id: An Ensembl gene identifier, versioned or not.
            transcript_ids: Unversioned Ensembl transcript identifiers the caller cares about,
                typically the ones it is offering in a picker. Any that this release never
                quantified are included in the result as
                :attr:`ExpressionStatus.NOT_IN_RELEASE`, so the gap is visible rather than
                silently absent. Passing none returns only what GTEx holds.

        Returns:
            The profile, carrying its dataset, annotation release and caveats.

        Raises:
            GeneNotInGtexError: If the gene is absent from this dataset's release.
            GtexResponseError: If a row is missing a field it is documented to have.
        """
        release_gene_id = self.release_gene_id(gene_id)
        catalogue = self.tissues()
        rows = self._median_rows(release_gene_id)
        measured = _group_by_transcript(rows, catalogue, release_gene_id)
        requested = tuple(dict.fromkeys(entry.split(".")[0] for entry in transcript_ids))
        absent = tuple(
            TranscriptExpression(
                transcript_id=identifier, status=ExpressionStatus.NOT_IN_RELEASE
            )
            for identifier in requested
            if identifier not in {entry.transcript_id for entry in measured}
        )
        return GeneExpression(
            gene_id=gene_id.split(".")[0],
            release_gene_id=release_gene_id,
            dataset_id=self.dataset_id,
            annotation_release=self.annotation_release,
            transcripts=measured + absent,
            provenance=Provenance.citation(
                f"GTEx Portal API v2, dataset {self.dataset_id} "
                f"(GENCODE {self.annotation_release})",
                f"Median transcript-level expression across {len(catalogue)} sample groups, "
                f"retrieved from {BASE_URL}",
            ),
            caveats=_caveats(self.annotation_release, absent),
        )

    def _median_rows(self, release_gene_id: str) -> list[Mapping[str, Any]]:
        """Fetch every page of median rows for one gene.

        Raises:
            GtexResponseError: If paging does not terminate within :data:`MAX_PAGES`.
        """
        collected: list[Mapping[str, Any]] = []
        page = 0
        while page < MAX_PAGES:
            payload = self._http.get(
                f"/api/v2/expression/medianTranscriptExpression"
                f"?gencodeId={release_gene_id}&datasetId={self.dataset_id}"
                f"&page={page}&itemsPerPage={PAGE_SIZE}"
            )
            collected.extend(_rows(payload, "medianTranscriptExpression"))
            paging = payload.get("paging_info")
            total_pages = (
                int(paging.get("numberOfPages", 1)) if isinstance(paging, dict) else 1
            )
            page += 1
            if page >= total_pages:
                return collected
        raise GtexResponseError(
            f"paging through {release_gene_id} did not finish within {MAX_PAGES} pages"
        )


def _rows(payload: JsonObject, what: str) -> list[Mapping[str, Any]]:
    """Return the ``data`` list of a v2 response, checking its shape.

    Raises:
        GtexResponseError: If ``data`` is missing or is not a list of objects.
    """
    data = payload.get("data")
    if not isinstance(data, list):
        raise GtexResponseError(
            f"{what} returned no 'data' list; got keys {sorted(payload)[:6]}"
        )
    for entry in data:
        if not isinstance(entry, dict):
            raise GtexResponseError(
                f"{what} returned a {type(entry).__name__} where an object was expected"
            )
    return [entry for entry in data if isinstance(entry, dict)]


def _sample_count(record: Mapping[str, Any]) -> int:
    """Return the RNA-seq donor count for a tissue, or ``0`` if it was not reported."""
    summary = record.get("rnaSeqSampleSummary")
    if not isinstance(summary, dict):
        return 0
    total = summary.get("totalCount")
    return int(total) if isinstance(total, int | float) else 0


def _group_by_transcript(
    rows: Sequence[Mapping[str, Any]],
    catalogue: Mapping[str, TissueMedian],
    release_gene_id: str,
) -> tuple[TranscriptExpression, ...]:
    """Turn flat median rows into one record per transcript.

    Rows keep the order GTEx returned them in within each transcript, which is its own
    ordering and carries no ranking; transcripts come back in first-seen order for the same
    reason.

    Raises:
        GtexResponseError: If a row lacks a transcript, a tissue or a numeric median, or
            belongs to a different gene than the one requested.
    """
    grouped: dict[str, list[TissueMedian]] = {}
    versions: dict[str, str] = {}
    for row in rows:
        versioned = str(row.get("transcriptId", ""))
        tissue_id = str(row.get("tissueSiteDetailId", ""))
        median = row.get("median")
        if not versioned or not tissue_id or not isinstance(median, int | float):
            raise GtexResponseError(
                f"a median row for {release_gene_id} is missing its transcript, tissue or "
                f"value: {dict(row)!r:.200}"
            )
        if str(row.get("gencodeId", release_gene_id)) != release_gene_id:
            raise GtexResponseError(
                f"asked for {release_gene_id} but a row reports {row.get('gencodeId')!r}"
            )
        stem = versioned.split(".")[0]
        versions.setdefault(stem, versioned)
        known = catalogue.get(tissue_id)
        grouped.setdefault(stem, []).append(
            TissueMedian(
                tissue_id=tissue_id,
                tissue_name=known.tissue_name if known else tissue_id,
                median=float(median),
                unit=str(row.get("unit") or "TPM"),
                ontology_id=(
                    known.ontology_id if known else str(row.get("ontologyId") or "")
                ),
                sample_count=known.sample_count if known else 0,
            )
        )
    return tuple(
        TranscriptExpression(
            transcript_id=stem,
            status=ExpressionStatus.MEASURED,
            medians=tuple(medians),
            release_id=versions[stem],
        )
        for stem, medians in grouped.items()
    )


def _caveats(
    annotation_release: str, absent: Sequence[TranscriptExpression]
) -> tuple[str, ...]:
    """Return the caveats that belong with every profile, plus any gap this one has."""
    caveats = [CELL_LINE_CAVEAT, BULK_TISSUE_CAVEAT, MEDIAN_CAVEAT, ISOFORM_CAVEAT]
    if absent:
        subject = "transcript is" if len(absent) == 1 else "transcripts are"
        caveats.append(
            f"{len(absent)} requested {subject} absent from GENCODE {annotation_release} and "
            f"never quantified here. That is a gap in this dataset's annotation, not a "
            f"measurement of zero."
        )
    return tuple(caveats)


def expression_for_gene(
    gene_id: str,
    transcript_ids: Sequence[str] = (),
    dataset_id: str = DEFAULT_DATASET,
    cache_dir: Path | None = None,
) -> GeneExpression:
    """Fetch one gene's transcript expression, caching under the configured cache root.

    Args:
        gene_id: An Ensembl gene identifier, versioned or not.
        transcript_ids: Transcripts the caller is offering, so gaps become visible.
        dataset_id: Which GTEx release to query.
        cache_dir: Override the cache location; defaults to the ``gtex`` cache directory.

    Returns:
        The profile.
    """
    client = GtexClient(
        dataset_id=dataset_id,
        cache_dir=settings.cache_dir("gtex") if cache_dir is None else cache_dir,
    )
    return client.expression_for(gene_id, transcript_ids)


__all__ = [
    "BULK_TISSUE_CAVEAT",
    "CELL_LINE_CAVEAT",
    "DATASET_ANNOTATION",
    "DEFAULT_DATASET",
    "GeneNotInGtexError",
    "GtexClient",
    "GtexError",
    "GtexResponseError",
    "UnknownGtexDatasetError",
    "expression_for_gene",
]
