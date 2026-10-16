"""A cached, retrying, offline-aware JSON fetcher for one external source.

Three behaviours, which every source in this project needs and which are easy to get subtly
wrong separately:

**Caching is atomic.** An entry is written beside its target and renamed, so an interrupted
write leaves either the old file or the new one and never half of one. A damaged entry counts
as a miss rather than an error: a truncated file from a closed laptop used to raise out of the
client on every later lookup of that record, for ever, with nothing telling a colleague to go
and delete it.

**Retrying is selective.** Transport failures, 5xx responses and a 200 carrying HTML are
retried with exponential backoff. A 4xx is not: a rejected parameter or a missing record is an
answer, and repeating the question does not change it.

**Offline is a refusal, not an empty result.** When external requests are disabled and the
cache has no entry, this raises and names what was missing. Returning nothing would be
indistinguishable from the source having no data, which is the error this project exists to
avoid making.
"""

import hashlib
import json
import logging
import os
import time
from collections.abc import Mapping
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import httpx

from sarna import settings

JsonObject = dict[str, Any]
"""A decoded JSON object."""

DEFAULT_TIMEOUT_SECONDS = 30.0

MAX_ATTEMPTS = 4
"""Attempts before giving up, including the first."""

BACKOFF_SECONDS = 0.5
"""Base delay between retries, doubled each attempt. Zero in tests."""

LOGGER = logging.getLogger(__name__)


class SourceError(RuntimeError):
    """Base class for failures talking to an external source."""


class SourceUnavailableError(SourceError):
    """The source could not be reached, or returned an error status."""


class SourceCacheMissError(SourceUnavailableError):
    """A required response is absent while external requests are disabled."""


class JsonHttpCache:
    """Fetches JSON from one base URL, caching every response on disk.

    Attributes:
        name: Short source name, used as the cache subdirectory and in messages.
        base_url: Everything before the request path, with no trailing slash.
        cache_dir: Where responses are stored, or ``None`` when caching is off.
        offline: Whether external requests are refused.
    """

    __slots__ = (
        "_backoff",
        "_client",
        "_offline_from_environment",
        "_timeout",
        "base_url",
        "cache_dir",
        "name",
        "offline",
    )

    def __init__(
        self,
        name: str,
        base_url: str,
        cache_dir: Path | None = None,
        timeout: float = DEFAULT_TIMEOUT_SECONDS,
        client: httpx.Client | None = None,
        backoff_seconds: float = BACKOFF_SECONDS,
        offline: bool | None = None,
    ) -> None:
        """Build a fetcher for one source.

        Args:
            name: Short identifier, e.g. ``"gtex"``. Appears in cache filenames and errors.
            base_url: Scheme and host, e.g. ``https://gtexportal.org``.
            cache_dir: Directory for cached responses. Created if absent; caching switches off
                with a warning if it cannot be.
            timeout: Per-request timeout in seconds.
            client: An ``httpx.Client`` to use instead of the module-level one, for testing.
            backoff_seconds: Base delay between retries, doubled each attempt.
            offline: Whether to refuse external requests and serve only cached responses.
                ``None`` takes the answer from :data:`sarna.settings.OFFLINE_VARIABLE`,
                resolved once here so one run cannot be half offline.
        """
        self.name = name
        self.base_url = base_url.rstrip("/")
        self.cache_dir = cache_dir
        self.offline = settings.is_offline() if offline is None else offline
        self._offline_from_environment = offline is None and self.offline
        self._timeout = timeout
        self._client = client
        self._backoff = backoff_seconds
        if cache_dir is not None:
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError as error:
                # A cache is an optimisation. A shared network folder that is unreachable off
                # the VPN must not turn a missing convenience into a tool that does not start.
                LOGGER.warning(
                    "caching is off: %s cannot be used as a cache directory (%s)",
                    cache_dir, error,
                )
                self.cache_dir = None

    def get(self, path: str) -> JsonObject:
        """Fetch ``path``, returning a parsed JSON object, using the cache when present.

        Args:
            path: Request path including any query string, e.g. ``/api/v2/thing?id=1``.

        Returns:
            The decoded object.

        Raises:
            SourceCacheMissError: If offline and the response is not cached.
            SourceUnavailableError: If every attempt failed, or the body was not an object.
        """
        cached = self._cache_path(path)
        from_cache = self._read_cache(cached)
        if from_cache is not None:
            return from_cache
        if self.offline:
            self._refuse_offline(path)
        payload = self._fetch(f"{self.base_url}{path}")
        self._write_cache(cached, payload)
        return payload

    def _cache_path(self, path: str) -> Path | None:
        """Return the cache file for ``path``, or ``None`` if caching is off."""
        if self.cache_dir is None:
            return None
        identity = json.dumps([self.base_url, path], separators=(",", ":"))
        digest = hashlib.sha256(identity.encode()).hexdigest()
        return self.cache_dir / f"{self.name}_{digest}.json"

    def _read_cache(self, cached: Path | None) -> JsonObject | None:
        """Return a cached object, or ``None`` if there is none or it is unusable."""
        if cached is None or not cached.exists():
            return None
        try:
            payload = json.loads(cached.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            LOGGER.warning("discarding an unreadable cache entry: %s", cached.name)
            return None
        if not isinstance(payload, dict):
            LOGGER.warning("discarding a cache entry that is not an object: %s", cached.name)
            return None
        return payload

    @staticmethod
    def _write_cache(cached: Path | None, payload: JsonObject) -> None:
        """Write a cache entry so a reader never sees half of it."""
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

    def _refuse_offline(self, path: str) -> None:
        """Raise rather than reach the network, naming what was missing.

        Raises:
            SourceCacheMissError: Always. Only called once a cache miss is established.
        """
        where = (
            "no cache is configured, so nothing can be served"
            if self.cache_dir is None
            else f"it is not in the cache at {self.cache_dir}"
        )
        recovery = (
            f"{settings.OFFLINE_VARIABLE} forces this mode. Unset it and restart to allow "
            "online retrieval."
            if self._offline_from_environment
            else "Allow online retrieval for this source, then retry to fetch and cache it."
        )
        raise SourceCacheMissError(
            f"Cache-only mode: {self.name} was not asked for {path} before, and {where}. "
            f"{recovery}"
        )

    def _fetch(self, url: str) -> JsonObject:
        """Perform the request, retrying what is worth retrying.

        Raises:
            SourceUnavailableError: If every attempt failed, or the body was not an object.
        """
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
            if response.is_error:
                last_problem = f"HTTP {response.status_code}: {response.text[:200]}"
                if response.status_code < httpx.codes.INTERNAL_SERVER_ERROR:
                    break
                continue
            try:
                decoded = response.json()
            except ValueError as error:
                # A 200 carrying HTML, which is how an overloaded origin reports a fault.
                last_problem = f"response was not JSON ({error})"
                continue
            if not isinstance(decoded, dict):
                raise SourceUnavailableError(
                    f"expected a JSON object from {url} but got {type(decoded).__name__}"
                )
            return decoded
        plural = "attempt" if attempts == 1 else "attempts"
        raise SourceUnavailableError(
            f"could not fetch {url} after {attempts} {plural}; last problem was {last_problem}"
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


__all__ = [
    "BACKOFF_SECONDS",
    "DEFAULT_TIMEOUT_SECONDS",
    "MAX_ATTEMPTS",
    "JsonHttpCache",
    "JsonObject",
    "SourceCacheMissError",
    "SourceError",
    "SourceUnavailableError",
]
