"""Ensembl REST client: gene symbol to transcript to TSS to promoter.

This is the adapter that makes Workflow A possible (a researcher types `SCN1A` and gets a
promoter) and it is the boundary where three silent assumptions in the legacy pipeline are
turned into explicit, checked decisions.

**The assembly is named, not inherited.** Ensembl's main endpoint serves GRCh38 and a separate
host serves GRCh37. Defect B15 was combining GRCh38 transcript coordinates with GRCh37 signal
tracks without noticing, so this client takes the assembly as a required argument, derives the
host from it, and **asserts the assembly the response reports back**. A server that changed its
default would fail loudly here rather than shifting every coordinate by a few hundred bases.

**The transcript is chosen, not assumed.** A gene has many transcripts with many TSSs, and a
promoter is a property of a transcript, not of a gene. :meth:`EnsemblClient.select_transcript`
returns the MANE Select transcript and otherwise **refuses**, listing what it found. Picking the
first, the longest or the lowest-numbered is how a tool ends up designing against a promoter the
user never saw.

**Coordinates convert once.** Ensembl REST is 1-based inclusive; this project is 0-based
half-open. The conversion happens here and nowhere else.

The client is offline-capable. Every response is cached on disk by request path, so a promoter
fetched once is reproducible without the network, which is what lets the design rules stay
testable with no external dependency.
"""

import json
import logging
import os
import re
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import httpx

from sarna import settings
from sarna.genome.anchors import DEFAULT_REGION_PADDING_BP as ANCHOR_REGION_PADDING_BP
from sarna.genome.anchors import PromoterAnchor
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.promoter import Promoter
from sarna.reproducibility import source_observation
from sarna.sequence.alphabet import Nucleic
from sarna.sequence.nucleic import NucleicSequence

JsonObject = dict[str, Any]
"""A decoded Ensembl JSON object, as the lookup and sequence endpoints return."""

JsonArray = list[dict[str, Any]]
"""A decoded Ensembl JSON array, as the overlap endpoints return."""

JsonPayload = JsonObject | JsonArray
"""Either shape. Ensembl returns an object from lookups and an array from overlaps."""

HOSTS: Mapping[Assembly, str] = {
    Assembly.GRCH38: "https://rest.ensembl.org",
    Assembly.GRCH37: "https://grch37.rest.ensembl.org",
}
"""Which host serves which assembly. Separate hosts, not a query parameter (ADR 0004)."""

DEFAULT_TIMEOUT_SECONDS = 30.0

DEFAULT_UPSTREAM_BP = 2000
"""Bases upstream of the TSS to fetch by default.

Wider than any rule window, so a user sees context around the region the rules permit. Matches
the -2000 to +500 promoter definition used in thesis chapter 4, which makes promoters fetched
here directly comparable with the designs screened there.
"""

DEFAULT_DOWNSTREAM_BP = 500
"""Bases downstream of the TSS to fetch by default, including the TSS itself."""

MAX_REGION_BP = 10_000_000
"""Ensembl's sequence endpoint refuses larger spans; checked here for a clearer message."""


LOGGER = logging.getLogger(__name__)

MAX_ATTEMPTS = 6
"""Attempts per request.

Ensembl's public endpoint intermittently serves an HTML error page with a 500 status, and
observed behaviour is that the same request succeeds on a later attempt. Four attempts over
seven seconds proved too few in practice; six over about a minute is enough, and a cached
response means the cost is paid once per promoter rather than once per run.
"""

BACKOFF_SECONDS = 2.0
"""Base delay between retries, doubled each attempt: 2, 4, 8, 16, 32 seconds."""

# The User-Agent is built per request from sarna.settings.user_agent(), so setting
# SARNA_CONTACT_EMAIL extends Ensembl the courtesy their usage guidance asks for. It used to be a
# module constant carrying a dangling "https://github.com/" and no address at all.


class EnsemblError(RuntimeError):
    """Base class for failures talking to Ensembl."""


class EnsemblUnavailableError(EnsemblError):
    """Raised when Ensembl cannot be reached or returns an error status."""


class EnsemblCacheMissError(EnsemblUnavailableError):
    """A required response is absent while external requests are disabled."""


class GeneNotFoundError(EnsemblError):
    """Raised when a gene symbol has no Ensembl record."""


class AssemblyMismatchResponseError(EnsemblError):
    """Raised when Ensembl reports a different assembly than the one requested.

    Never ignored, and never corrected by adjusting coordinates: a response from the wrong
    assembly means the host or its default changed, and continuing would reintroduce defect B15.
    """


class NoManeSelectTranscriptError(EnsemblError):
    """Raised when no MANE Select transcript exists, so the choice is the user's.

    Carries the alternatives so a caller can present them. This is a question, not a failure:
    many genes legitimately lack a MANE Select transcript, and for those the promoter depends on
    a judgement the tool must not make silently.
    """

    def __init__(self, symbol: str, alternatives: "tuple[Transcript, ...]") -> None:
        """Record the gene and the transcripts the user must choose between."""
        self.symbol = symbol
        self.alternatives = alternatives
        listing = "\n".join(
            f"  {t.ensembl_id}  {t.biotype:<16} TSS at {t.contig}:{t.tss_coordinate}"
            f"{' (Ensembl canonical)' if t.is_canonical else ''}"
            for t in alternatives
        )
        super().__init__(
            f"{symbol} has no MANE Select transcript, so which promoter to design against is "
            f"a scientific choice this tool will not make for you. "
            f"{len(alternatives)} protein-coding transcript(s) found:\n{listing}\n"
            f"Pass one to promoter_for() explicitly."
        )


@dataclass(frozen=True, slots=True)
class Gene:
    """An Ensembl gene record, with its span on a named assembly."""

    ensembl_id: str
    symbol: str
    interval: GenomicInterval
    description: str

    @property
    def contig(self) -> str:
        """Return the contig the gene lies on."""
        return self.interval.contig

    @property
    def strand(self) -> Strand:
        """Return the strand the gene is transcribed from."""
        return self.interval.strand


@dataclass(frozen=True, slots=True)
class Transcript:
    """One transcript of a gene, and the TSS that follows from its orientation."""

    ensembl_id: str
    biotype: str
    interval: GenomicInterval
    is_mane_select: bool
    is_canonical: bool
    refseq_match: str | None = None
    """The RefSeq accession a MANE Select transcript agrees with, e.g. ``NM_001165963.4``.

    Recorded because MANE Select's value is precisely that two independent annotation projects
    agreed, and the accession is what lets a reader check that agreement.
    """

    @property
    def contig(self) -> str:
        """Return the contig the transcript lies on."""
        return self.interval.contig

    @property
    def strand(self) -> Strand:
        """Return the strand the transcript is transcribed from."""
        return self.interval.strand

    @property
    def tss_coordinate(self) -> int:
        """Return the TSS as a 0-based genomic coordinate.

        The transcript's 5' end **in its own orientation**: the lowest coordinate on the plus
        strand and the highest on the minus strand. Using ``start`` on both is defect B1, and
        for a minus-strand gene it places the promoter at the 3' end of the transcript, which is
        the opposite end of the gene.
        """
        return self.interval.five_prime

    @property
    def is_protein_coding(self) -> bool:
        """Return whether this transcript codes for protein."""
        return self.biotype == "protein_coding"


class EnsemblClient:
    """A caching Ensembl REST client, pinned to one assembly.

    Attributes:
        assembly: The assembly every coordinate from this client belongs to.
        cache_dir: Where responses are stored, or ``None`` for no caching.
    """

    __slots__ = (
        "_backoff",
        "_client",
        "_offline_from_environment",
        "_timeout",
        "assembly",
        "base_url",
        "cache_dir",
        "observations",
        "offline",
    )

    def __init__(
        self,
        assembly: Assembly = Assembly.GRCH38,
        cache_dir: Path | None = None,
        timeout: float = DEFAULT_TIMEOUT_SECONDS,
        client: httpx.Client | None = None,
        backoff_seconds: float = BACKOFF_SECONDS,
        offline: bool | None = None,
    ) -> None:
        """Build a client for one assembly.

        Args:
            assembly: Which assembly to query. Determines the host; never a default.
            cache_dir: Directory for cached responses. Created if absent.
            timeout: Per-request timeout in seconds.
            client: An ``httpx.Client`` to use instead of creating one, for testing.
            backoff_seconds: Base delay between retries, doubled each attempt. Zero in tests,
                so a retry path can be exercised without waiting for it.
            offline: Whether to refuse external lookups and serve only cached responses.
                ``None`` takes the answer from :data:`sarna.settings.OFFLINE_VARIABLE`. Resolved
                once here rather than per request, so one run cannot be half offline.
        """
        self.assembly = assembly
        self.base_url = HOSTS[assembly]
        self.cache_dir = cache_dir
        self.offline = settings.is_offline() if offline is None else offline
        self._offline_from_environment = offline is None and self.offline
        self._timeout = timeout
        self._client = client
        self._backoff = backoff_seconds
        self.observations: list[dict[str, Any]] = []
        if cache_dir is not None:
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError as error:
                # A cache is an optimisation. INSTALL.md offers a shared network folder so a
                # group fetches each gene once, and a laptop off the VPN then has an unreachable
                # path; raising here would turn a missing convenience into a tool that does not
                # start. Caching switches off and every lookup goes to Ensembl as before.
                LOGGER.warning(
                    "caching is off: %s cannot be used as a cache directory (%s)",
                    cache_dir, error,
                )
                self.cache_dir = None

    def _read_cache(self, cached: Path | None) -> JsonPayload | None:
        """Return a cached payload, or ``None`` if there is none or it is unusable.

        A damaged entry is a miss, not a failure. Interrupt the app mid-write -- Ctrl-C, a
        closed laptop, a full disk -- and the file is left truncated; a bare ``json.loads``
        then raised out of the client on every later lookup of that gene, for ever, and nothing
        told a colleague to go and delete it. Refetching replaces it.
        """
        if cached is None or not cached.exists():
            return None
        try:
            payload: JsonPayload = json.loads(cached.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            LOGGER.warning("discarding an unreadable cache entry: %s", cached.name)
            return None
        return payload

    @staticmethod
    def _write_cache(cached: Path | None, payload: JsonPayload) -> None:
        """Write a cache entry so a reader never sees half of it.

        Written beside the target and renamed, which is atomic within a directory. That is what
        makes an interrupted write harmless, and what would make a cache directory shared
        between colleagues safe to read while somebody else is filling it.
        """
        if cached is None:
            return
        temporary: Path | None = None
        try:
            with NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=cached.parent,
                prefix=f"{cached.name}.", suffix=".part", delete=False,
            ) as handle:
                temporary = Path(handle.name)
                json.dump(payload, handle)
            os.replace(temporary, cached)
        except OSError as error:
            # A cache that cannot be written makes the tool slower, not broken.
            LOGGER.warning("could not cache %s: %s", cached.name, error)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def _cache_path(self, path: str) -> Path | None:
        """Return the cache file for ``path``, or ``None`` if caching is off."""
        if self.cache_dir is None:
            return None
        safe = path.replace("/", "_").replace("?", "_").replace(";", "_").replace("=", "-")
        return self.cache_dir / f"{self.assembly}_{safe}.json"

    def _get(self, path: str) -> JsonObject:
        """Fetch ``path``, returning parsed JSON, using the cache when present.

        Retries transport failures, 5xx responses and unparseable bodies. Ensembl's public
        endpoint intermittently serves an HTML error page with a 500 status under load, and
        observed behaviour is that the same request succeeds within a few attempts, so a
        single-shot GET would fail for reasons that have nothing to do with the query.

        A 404 is not retried: a missing gene is an answer, not a transient fault.

        Raises:
            GeneNotFoundError: If Ensembl has no such record.
            EnsemblUnavailableError: If every attempt failed.
        """
        cached = self._cache_path(path)
        from_cache = self._read_cache(cached)
        if isinstance(from_cache, dict):
            self.observations.append(source_observation(f"{self.base_url}{path}", from_cache, True))
            return from_cache
        if self.offline:
            self._refuse_offline(path)

        url = f"{self.base_url}{path}"
        payload = self._get_with_retries(url, path)
        if not isinstance(payload, dict):
            raise EnsemblUnavailableError(
                f"expected a JSON object from {url} but got {type(payload).__name__}"
            )
        self._write_cache(cached, payload)
        self.observations.append(source_observation(url, payload, False))
        return payload

    def _get_with_retries(self, url: str, path: str) -> JsonPayload:
        """Perform the request, retrying what is worth retrying."""
        headers = {"Accept": "application/json", "User-Agent": settings.user_agent()}
        last_problem = "no attempt was made"
        attempts = 0
        for attempt in range(MAX_ATTEMPTS):
            if attempt:
                self._sleep(self._backoff * (2 ** (attempt - 1)))
            attempts += 1
            try:
                response = self._request(url, headers)
            except httpx.HTTPError as error:
                last_problem = f"{type(error).__name__}: {error}"
                continue
            if _is_missing_record(response):
                raise GeneNotFoundError(
                    f"Ensembl has no record at {path}. Check the symbol's spelling and case: "
                    f"Ensembl matches HGNC symbols exactly."
                )
            if response.is_error:
                last_problem = f"HTTP {response.status_code}: {response.text[:160]}"
                if response.status_code < httpx.codes.INTERNAL_SERVER_ERROR:
                    break
                continue
            try:
                decoded: JsonPayload = response.json()
                return decoded
            except ValueError as error:
                # A 200 carrying HTML. Seen from Ensembl's error pages under load.
                last_problem = f"response was not JSON ({error})"
        plural = "attempt" if attempts == 1 else "attempts"
        raise EnsemblUnavailableError(
            f"could not fetch {url} after {attempts} {plural}; last problem was "
            f"{last_problem}. The design engine works offline on a pasted promoter sequence; "
            f"only gene lookup needs the network."
        )

    def _request(self, url: str, headers: Mapping[str, str]) -> httpx.Response:
        """Issue one GET, through the injected client if there is one."""
        if self._client is not None:
            return self._client.get(url, headers=dict(headers), timeout=self._timeout)
        return httpx.get(
            url, headers=dict(headers), timeout=self._timeout, follow_redirects=True
        )

    @staticmethod
    def _sleep(seconds: float) -> None:
        """Wait between retries. Skipped entirely at a zero backoff."""
        if seconds > 0:
            time.sleep(seconds)

    def overlap(self, interval: GenomicInterval, feature: str) -> JsonArray:
        """Return Ensembl features of type ``feature`` overlapping ``interval``, as raw records.

        Deliberately returns unparsed records. Interpreting a regulatory feature is the
        regulatory layer's job, and this client's job is to be the only place that knows how to
        talk to Ensembl; keeping the two apart stops the genome adapter from depending on the
        annotation model it feeds.

        Args:
            interval: The region to query, on this client's assembly.
            feature: An Ensembl feature type, e.g. ``"regulatory"`` or ``"motif"``.

        Returns:
            The records, which may be empty. Empty means "nothing annotated", which is not the
            same as "nothing there".

        Raises:
            ValueError: If the interval is on another assembly.
            EnsemblUnavailableError: If Ensembl cannot be reached.
        """
        if interval.assembly is not self.assembly:
            raise ValueError(
                f"interval is on {interval.assembly} but this client serves {self.assembly}"
            )
        # Ensembl is 1-based inclusive. The overlap endpoint ignores strand, so none is sent.
        region = f"{interval.contig}:{interval.start + 1}-{interval.end}"
        return self._get_array(f"/overlap/region/homo_sapiens/{region}?feature={feature}")

    def _refuse_offline(self, path: str) -> None:
        """Raise rather than reach the network, naming what was missing.

        Raises:
            EnsemblCacheMissError: Always. Only called once a cache miss is established.
        """
        where = (
            "no cache is configured, so nothing can be served"
            if self.cache_dir is None
            else f"it is not in the cache at {self.cache_dir}"
        )
        recovery = (
            f"{settings.OFFLINE_VARIABLE} forces this mode. Unset it and restart the app "
            "to allow online retrieval."
            if self._offline_from_environment else
            "Turn off 'Use cached Ensembl data only' in the app, or set offline=False "
            "in the client, then retry to fetch and cache this response."
        )
        raise EnsemblCacheMissError(
            f"Cache-only mode: {path} was not requested, and {where}. {recovery}"
        )

    def _get_array(self, path: str) -> JsonArray:
        """Fetch ``path``, which returns a JSON array rather than an object."""
        cached = self._cache_path(path)
        from_cache = self._read_cache(cached)
        if isinstance(from_cache, list):
            self.observations.append(source_observation(f"{self.base_url}{path}", from_cache, True))
            return from_cache
        if self.offline:
            self._refuse_offline(path)
        payload = self._get_with_retries(f"{self.base_url}{path}", path)
        if not isinstance(payload, list):
            raise EnsemblUnavailableError(
                f"expected a JSON array from {path} but got {type(payload).__name__}"
            )
        self._write_cache(cached, payload)
        self.observations.append(source_observation(f"{self.base_url}{path}", payload, False))
        return payload

    def _check_assembly(self, reported: str | None, what: str) -> None:
        """Verify Ensembl's reported assembly matches this client's.

        Raises:
            AssemblyMismatchResponseError: If they differ.
        """
        if reported is None:
            return
        if reported != str(self.assembly):
            raise AssemblyMismatchResponseError(
                f"requested {self.assembly} but Ensembl reports {reported!r} for {what}. "
                f"Continuing would mix assemblies, which is defect B15. Check the host for "
                f"{self.assembly} ({self.base_url})."
            )

    def gene(self, symbol: str) -> Gene:
        """Look up a gene by its HGNC symbol or human Ensembl stable gene ID.

        Args:
            symbol: The gene symbol, e.g. ``"SCN1A"``, or ``ENSG`` stable ID. Gene
                symbols are case-sensitive at Ensembl; a gene-ID version suffix is accepted
                but resolves the current record, rather than pinning an annotation release.

        Returns:
            The gene, with its span on this client's assembly.

        Raises:
            GeneNotFoundError: If no gene has that symbol.
            AssemblyMismatchResponseError: If Ensembl reports a different assembly.
            EnsemblUnavailableError: If Ensembl cannot be reached.
        """
        if re.fullmatch(r"ENSG\d{11}(?:\.\d+)?", symbol):
            return self.gene_by_id(symbol)
        payload = self._get(f"/lookup/symbol/homo_sapiens/{symbol}?expand=0")
        self._check_assembly(payload.get("assembly_name"), f"gene {symbol}")
        return Gene(
            ensembl_id=payload["id"],
            symbol=payload.get("display_name", symbol),
            interval=self._interval_from(payload),
            description=payload.get("description") or "",
        )

    def gene_by_id(self, identifier: str) -> Gene:
        """Resolve a human stable gene ID, refusing other species or feature types.

        A version suffix identifies the input but does not freeze the annotation release:
        the unversioned stable ID is looked up on this client's assembly endpoint.
        """
        if not re.fullmatch(r"ENSG\d{11}(?:\.\d+)?", identifier):
            raise ValueError("Enter a human Ensembl gene ID such as ENSG00000144285.")
        stable_id = identifier.split(".", 1)[0]
        payload = self._get(f"/lookup/id/{stable_id}?expand=0")
        if payload.get("species") != "homo_sapiens" or payload.get("object_type") != "Gene":
            raise ValueError("The Ensembl ID must identify a human Gene record.")
        if payload.get("id") != stable_id:
            raise EnsemblUnavailableError(
                f"requested gene {stable_id} but Ensembl returned {payload.get('id')!r}"
            )
        if not payload.get("assembly_name"):
            raise EnsemblUnavailableError(f"Ensembl omitted the assembly for gene {stable_id}.")
        self._check_assembly(payload["assembly_name"], f"gene {stable_id}")
        return Gene(
            ensembl_id=stable_id, symbol=payload.get("display_name") or stable_id,
            interval=self._interval_from(payload), description=payload.get("description") or "",
        )

    def transcripts(self, gene: Gene) -> tuple[Transcript, ...]:
        """Return every transcript of ``gene``, in Ensembl's order.

        Raises:
            EnsemblUnavailableError: If Ensembl cannot be reached.
        """
        # mane=1 is required: without it Ensembl omits MANE annotations entirely and every
        # transcript looks equally unremarkable, which would send select_transcript to its
        # "no MANE Select" branch for genes that have one.
        payload = self._get(f"/lookup/id/{gene.ensembl_id}?expand=1;mane=1")
        self._check_assembly(payload.get("assembly_name"), f"gene {gene.symbol}")
        canonical = payload.get("canonical_transcript") or ""
        return tuple(
            Transcript(
                ensembl_id=entry["id"],
                biotype=entry.get("biotype", "unknown"),
                interval=self._interval_from(entry),
                is_mane_select=_is_mane_select(entry),
                # Ensembl reports the canonical transcript with a version suffix.
                is_canonical=canonical.split(".")[0] == entry["id"],
                refseq_match=_refseq_match(entry),
            )
            for entry in payload.get("Transcript", [])
        )

    def select_transcript(self, gene: Gene) -> Transcript:
        """Return the gene's MANE Select transcript.

        MANE Select is a single transcript per protein-coding gene agreed between Ensembl and
        RefSeq, which makes it the one choice that is both reproducible and externally
        justifiable. Chapter 4 of the thesis used MANE Select for the same reason.

        Returns:
            The MANE Select transcript.

        Raises:
            NoManeSelectTranscriptError: If there is none. The error carries the protein-coding
                alternatives so the caller can ask the user, rather than guessing.
        """
        transcripts = self.transcripts(gene)
        mane = [t for t in transcripts if t.is_mane_select]
        if len(mane) == 1:
            return mane[0]
        if len(mane) > 1:
            raise NoManeSelectTranscriptError(gene.symbol, tuple(mane))
        coding = tuple(t for t in transcripts if t.is_protein_coding)
        raise NoManeSelectTranscriptError(gene.symbol, coding or transcripts)

    def sequence(self, interval: GenomicInterval) -> NucleicSequence:
        """Fetch the reference sequence for ``interval`` as stored, on the **plus strand**.

        Returns the plus strand regardless of the interval's strand, because that is what a
        reference genome holds. :meth:`Promoter.from_reference_sequence` performs the flip for a
        minus-strand feature, in one place.

        Raises:
            ValueError: If the interval is larger than Ensembl will serve.
            EnsemblUnavailableError: If Ensembl cannot be reached.
        """
        if interval.assembly is not self.assembly:
            raise ValueError(
                f"interval is on {interval.assembly} but this client serves {self.assembly}"
            )
        if interval.length > MAX_REGION_BP:
            raise ValueError(
                f"region of {interval.length} bp exceeds Ensembl's {MAX_REGION_BP} bp limit"
            )
        # Ensembl is 1-based inclusive; this project is 0-based half-open. Converted here only.
        region = f"{interval.contig}:{interval.start + 1}..{interval.end}"
        payload = self._get(
            f"/sequence/region/homo_sapiens/{region}?coord_system_version={self.assembly}"
        )
        return NucleicSequence.parse(payload["seq"], Nucleic.DNA)

    def promoter_for(
        self,
        transcript: Transcript,
        upstream_bp: int = DEFAULT_UPSTREAM_BP,
        downstream_bp: int = DEFAULT_DOWNSTREAM_BP,
    ) -> Promoter:
        """Fetch the promoter of ``transcript``, oriented 5' to 3' in its own direction.

        Args:
            transcript: The transcript whose TSS anchors the promoter.
            upstream_bp: Bases upstream of the TSS to include.
            downstream_bp: Bases downstream of the TSS, including the TSS itself.

        Returns:
            The promoter, carrying its genomic interval and the TSS-relative position of its
            first base, so every candidate designed in it can report both.

        Raises:
            ValueError: If the window is empty or runs off the end of the contig.
            EnsemblUnavailableError: If Ensembl cannot be reached.

        Note:
            A thin wrapper over :meth:`promoter_for_anchor`. An annotated TSS is one kind of
            promoter anchor, so routing both through the same code keeps the strand arithmetic in
            one place instead of two that have to be kept in step.
        """
        return self.promoter_for_anchor(
            anchor_for_transcript(transcript), upstream_bp, downstream_bp
        )

    def promoter_for_anchor(
        self,
        anchor: PromoterAnchor,
        upstream_bp: int = DEFAULT_UPSTREAM_BP,
        downstream_bp: int = DEFAULT_DOWNSTREAM_BP,
        region_padding_bp: int = ANCHOR_REGION_PADDING_BP,
    ) -> Promoter:
        """Fetch the promoter a :class:`~sarna.genome.anchors.PromoterAnchor` describes.

        The general form of :meth:`promoter_for`, which is the special case where the anchor is an
        annotated transcript's TSS. An anchor with no reference point yields a promoter with no
        ``tss_start``, so the TSS rules report not-applicable instead of measuring from an invented
        origin.

        Raises:
            ValueError: If the anchor is on another assembly, or its window runs off the contig.
            EnsemblUnavailableError: If Ensembl cannot be reached.
        """
        if anchor.interval.assembly is not self.assembly:
            raise ValueError(
                f"anchor {anchor.id!r} is on {anchor.interval.assembly} but this client serves "
                f"{self.assembly}"
            )
        window = anchor.window(upstream_bp, downstream_bp, region_padding_bp)
        return Promoter.from_reference_sequence(
            reference_sequence=self.sequence(window.interval),
            interval=window.interval,
            tss_start=window.tss_start,
        )

    def _interval_from(self, payload: Mapping[str, Any]) -> GenomicInterval:
        """Build an interval from an Ensembl record, converting to 0-based half-open."""
        return GenomicInterval(
            assembly=self.assembly,
            contig=str(payload["seq_region_name"]),
            start=int(payload["start"]) - 1,
            end=int(payload["end"]),
            strand=Strand.PLUS if int(payload["strand"]) > 0 else Strand.MINUS,
        )


MISSING_RECORD_PHRASES = ("no valid lookup", "not found", "does not exist")
"""Phrases Ensembl uses when the record is absent rather than the service unwell.

Ensembl answers an unknown gene symbol with **HTTP 400**, not 404, and a body reading "No valid
lookup found for symbol X". Treating that as a transport failure would tell a user the service
is down when they have in fact mistyped a gene name, and would retry five more times first.
Found by the live tests; the offline tests cannot see it.
"""


def _is_missing_record(response: httpx.Response) -> bool:
    """Return whether a response means "no such record" rather than "try again"."""
    if response.status_code == httpx.codes.NOT_FOUND:
        return True
    if response.status_code != httpx.codes.BAD_REQUEST:
        return False
    body = response.text.lower()
    return any(phrase in body for phrase in MISSING_RECORD_PHRASES)


def _mane_annotations(entry: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Return the MANE annotations on a transcript record, if any."""
    return list(entry.get("MANE") or [])


def _refseq_match(entry: Mapping[str, Any]) -> str | None:
    """Return the RefSeq accession a MANE Select annotation agrees with, if present."""
    for annotation in _mane_annotations(entry):
        match = annotation.get("refseq_match")
        if match:
            return str(match)
    return None


def anchor_for_transcript(transcript: Transcript) -> PromoterAnchor:
    """Return the promoter anchor an annotated transcript's TSS defines.

    The bridge between the ordinary path and the general one: a transcript's TSS is one kind of
    anchor, and expressing it as one means :meth:`EnsemblClient.promoter_for` and
    :meth:`EnsemblClient.promoter_for_anchor` share their window arithmetic.
    """
    return PromoterAnchor.from_transcript(
        gene=transcript.ensembl_id,
        transcript_id=transcript.ensembl_id,
        tss=transcript.tss_coordinate,
        interval=transcript.interval,
    )


def _is_mane_select(entry: Mapping[str, Any]) -> bool:
    """Return whether an Ensembl transcript record is MANE Select.

    Ensembl has reported this in more than one shape over time (a top-level ``is_mane`` flag and
    an ``MANE`` list of annotations) so both are checked rather than relying on the current one.
    """
    if entry.get("is_mane") in (True, 1, "1"):
        return True
    return any(
        str(annotation.get("type", "")).lower().startswith("mane_select")
        for annotation in _mane_annotations(entry)
    )


def promoter_for_symbol(
    symbol: str,
    assembly: Assembly = Assembly.GRCH38,
    upstream_bp: int = DEFAULT_UPSTREAM_BP,
    downstream_bp: int = DEFAULT_DOWNSTREAM_BP,
    cache_dir: Path | None = None,
) -> tuple[Gene, Transcript, Promoter]:
    """Resolve a gene symbol all the way to a promoter, the Workflow A path.

    Returns the gene and transcript alongside the promoter, because a promoter without the
    transcript it came from cannot be audited.

    Raises:
        GeneNotFoundError: If the symbol is unknown.
        NoManeSelectTranscriptError: If the transcript choice is the user's to make.
        EnsemblUnavailableError: If Ensembl cannot be reached.
    """
    client = EnsemblClient(assembly=assembly, cache_dir=cache_dir)
    gene = client.gene(symbol)
    transcript = client.select_transcript(gene)
    promoter = client.promoter_for(transcript, upstream_bp, downstream_bp)
    return gene, transcript, promoter


__all__ = [
    "ANCHOR_REGION_PADDING_BP",
    "DEFAULT_DOWNSTREAM_BP",
    "DEFAULT_UPSTREAM_BP",
    "HOSTS",
    "AssemblyMismatchResponseError",
    "EnsemblClient",
    "EnsemblError",
    "EnsemblUnavailableError",
    "Gene",
    "GeneNotFoundError",
    "NoManeSelectTranscriptError",
    "Transcript",
    "anchor_for_transcript",
    "promoter_for_symbol",
]
