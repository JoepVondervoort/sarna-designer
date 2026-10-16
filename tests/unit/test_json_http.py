"""Tests for the shared cached JSON fetcher.

Three behaviours, each of which has a wrong version that looks fine until it is not: an
interrupted cache write that leaves a half file, a retry loop that hammers a 4xx, and an
offline mode that returns nothing instead of refusing.
"""

import json
from collections.abc import Mapping
from pathlib import Path

import httpx
import pytest

from sarna.sources.json_http import (
    MAX_ATTEMPTS,
    JsonHttpCache,
    SourceCacheMissError,
    SourceUnavailableError,
)

PAYLOAD = {"data": [{"value": 1}]}


class CountingTransport:
    """Returns a scripted sequence of responses, repeating the last one, counting calls."""

    def __init__(self, *responses: object) -> None:
        """Record the responses to serve in order."""
        self.responses = list(responses)
        self.requested: list[str] = []

    def get(self, url: str, headers: Mapping[str, str], timeout: float) -> httpx.Response:
        """Serve the next scripted response."""
        self.requested.append(url)
        index = min(len(self.requested) - 1, len(self.responses) - 1)
        reply = self.responses[index]
        if isinstance(reply, Exception):
            raise reply
        if isinstance(reply, int):
            return httpx.Response(reply, text="<html>error</html>")
        if isinstance(reply, str):
            return httpx.Response(200, text=reply)
        return httpx.Response(200, json=reply)


def cache_for(
    *responses: object, cache_dir: Path | None = None, offline: bool | None = False
) -> tuple[JsonHttpCache, CountingTransport]:
    """Build a fetcher backed by a scripted transport, with retries made instant."""
    transport = CountingTransport(*responses)
    fetcher = JsonHttpCache(
        name="test", base_url="https://example.invalid",
        client=transport,  # type: ignore[arg-type]
        backoff_seconds=0.0, cache_dir=cache_dir, offline=offline,
    )
    return fetcher, transport


class TestRetryingWhatIsWorthRetrying:
    def test_a_server_error_is_retried_and_can_succeed(self) -> None:
        fetcher, transport = cache_for(500, 500, PAYLOAD)
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 3

    def test_a_client_error_is_not_retried(self) -> None:
        # A rejected parameter or a missing record is an answer. GTEx answers 422 for a
        # dataset it does not have; asking again four times does not help.
        fetcher, transport = cache_for(422)
        with pytest.raises(SourceUnavailableError, match="HTTP 422"):
            fetcher.get("/thing")
        assert len(transport.requested) == 1

    def test_a_transport_failure_is_retried(self) -> None:
        fetcher, transport = cache_for(httpx.ConnectError("no route"), PAYLOAD)
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 2

    def test_a_two_hundred_carrying_html_is_retried(self) -> None:
        # Observed from overloaded origins, Ensembl's included: a 200 whose body is an error
        # page. Treating it as success would hand a caller nothing with no complaint.
        fetcher, transport = cache_for("<html>down</html>", PAYLOAD)
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 2

    def test_it_gives_up_after_the_documented_number_of_attempts(self) -> None:
        fetcher, transport = cache_for(500)
        with pytest.raises(SourceUnavailableError, match=f"after {MAX_ATTEMPTS} attempts"):
            fetcher.get("/thing")
        assert len(transport.requested) == MAX_ATTEMPTS

    def test_a_json_array_is_refused_without_retrying(self) -> None:
        # Not transient, so retrying would only delay the error. Every endpoint this is used
        # for returns an object.
        fetcher, transport = cache_for([1, 2, 3])
        with pytest.raises(SourceUnavailableError, match="got list"):
            fetcher.get("/thing")
        assert len(transport.requested) == 1


class TestCaching:
    def test_a_second_request_is_served_from_disk(self, tmp_path: Path) -> None:
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path)
        assert fetcher.get("/thing") == PAYLOAD
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 1

    def test_different_paths_are_cached_separately(self, tmp_path: Path) -> None:
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/one")
        fetcher.get("/two")
        assert len(transport.requested) == 2
        assert len(list(tmp_path.glob("*.json"))) == 2

    def test_a_query_string_is_part_of_the_cache_key(self, tmp_path: Path) -> None:
        # Two genes differ only in a query parameter; sharing a cache entry would serve one
        # gene's numbers for another.
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/gene?id=A")
        fetcher.get("/gene?id=B")
        assert len(transport.requested) == 2

    def test_a_cache_filename_contains_no_path_separators(self, tmp_path: Path) -> None:
        fetcher, _ = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/api/v2/thing?id=1&x=2")
        written = list(tmp_path.glob("*.json"))
        assert len(written) == 1
        assert written[0].parent == tmp_path

    def test_a_truncated_cache_entry_is_a_miss_rather_than_a_permanent_failure(
        self, tmp_path: Path
    ) -> None:
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/thing")
        cached = next(iter(tmp_path.glob("*.json")))
        cached.write_text('{"data": [{"val', encoding="utf-8")
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 2
        assert json.loads(cached.read_text(encoding="utf-8")) == PAYLOAD

    def test_a_cache_entry_that_is_not_an_object_is_a_miss(self, tmp_path: Path) -> None:
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/thing")
        next(iter(tmp_path.glob("*.json"))).write_text("[1, 2]", encoding="utf-8")
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 2

    def test_no_partial_file_is_left_behind_after_a_successful_write(
        self, tmp_path: Path
    ) -> None:
        fetcher, _ = cache_for(PAYLOAD, cache_dir=tmp_path)
        fetcher.get("/thing")
        assert list(tmp_path.glob("*.part")) == []

    def test_a_failed_fetch_writes_no_cache_entry(self, tmp_path: Path) -> None:
        fetcher, _ = cache_for(500, cache_dir=tmp_path)
        with pytest.raises(SourceUnavailableError):
            fetcher.get("/thing")
        assert list(tmp_path.iterdir()) == []

    def test_an_unusable_cache_directory_switches_caching_off_rather_than_failing(
        self, tmp_path: Path
    ) -> None:
        # A shared network folder named in INSTALL.md is unreachable off the VPN. That must
        # cost speed, not start-up.
        blocker = tmp_path / "file"
        blocker.write_text("not a directory", encoding="utf-8")
        fetcher, transport = cache_for(PAYLOAD, cache_dir=blocker / "cache")
        assert fetcher.cache_dir is None
        assert fetcher.get("/thing") == PAYLOAD
        assert len(transport.requested) == 1


class TestOfflineIsARefusal:
    def test_it_raises_rather_than_returning_nothing(self, tmp_path: Path) -> None:
        fetcher, transport = cache_for(PAYLOAD, cache_dir=tmp_path, offline=True)
        with pytest.raises(SourceCacheMissError, match="Cache-only mode"):
            fetcher.get("/thing")
        assert transport.requested == []

    def test_the_message_names_the_source_the_path_and_the_cache(
        self, tmp_path: Path
    ) -> None:
        fetcher, _ = cache_for(PAYLOAD, cache_dir=tmp_path, offline=True)
        with pytest.raises(SourceCacheMissError) as raised:
            fetcher.get("/thing")
        message = str(raised.value)
        assert "test" in message
        assert "/thing" in message
        assert str(tmp_path) in message

    def test_it_says_so_when_there_is_no_cache_at_all(self) -> None:
        fetcher, _ = cache_for(PAYLOAD, offline=True)
        with pytest.raises(SourceCacheMissError, match="no cache is configured"):
            fetcher.get("/thing")

    def test_a_cached_response_is_still_served(self, tmp_path: Path) -> None:
        warm, _ = cache_for(PAYLOAD, cache_dir=tmp_path)
        warm.get("/thing")
        offline, transport = cache_for(PAYLOAD, cache_dir=tmp_path, offline=True)
        assert offline.get("/thing") == PAYLOAD
        assert transport.requested == []

    def test_the_environment_forced_mode_names_the_variable(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("SARNA_OFFLINE", "1")
        fetcher = JsonHttpCache(
            name="test", base_url="https://example.invalid", cache_dir=tmp_path,
        )
        assert fetcher.offline is True
        with pytest.raises(SourceCacheMissError, match="SARNA_OFFLINE"):
            fetcher.get("/thing")


class TestTheRequestItSends:
    def test_it_identifies_the_tool_in_the_user_agent(self) -> None:
        sent: list[dict[str, str]] = []

        class Recording:
            def get(
                self, url: str, headers: dict[str, str], timeout: float
            ) -> httpx.Response:
                sent.append(headers)
                return httpx.Response(200, json=PAYLOAD)

        fetcher = JsonHttpCache(
            name="test", base_url="https://example.invalid",
            client=Recording(),  # type: ignore[arg-type]
            offline=False,
        )
        fetcher.get("/thing")
        assert sent[0]["Accept"] == "application/json"
        assert "sarna-designer" in sent[0]["User-Agent"]

    def test_the_base_url_is_joined_without_a_doubled_slash(self) -> None:
        fetcher = JsonHttpCache(
            name="test", base_url="https://example.invalid/",
            client=CountingTransport(PAYLOAD),  # type: ignore[arg-type]
            offline=False,
        )
        assert fetcher.base_url == "https://example.invalid"


@pytest.mark.parametrize(('left', 'right'), [
    ('/x?a=b', '/x/a=b'),
    ('/' + 'a' * 200 + '?id=1', '/' + 'a' * 200 + '?id=2'),
])
def test_distinct_requests_never_reuse_a_colliding_cache_entry(
    tmp_path: Path, left: str, right: str,
) -> None:
    fetcher, transport = cache_for({'id': 1}, {'id': 2}, cache_dir=tmp_path)
    assert fetcher.get(left) == {'id': 1}
    assert fetcher.get(right) == {'id': 2}
    assert fetcher.get(left) == {'id': 1}
    assert len(transport.requested) == 2


def test_cache_identity_includes_source_endpoint(tmp_path: Path) -> None:
    first, _ = cache_for({'id': 1}, cache_dir=tmp_path)
    second, _ = cache_for({'id': 2}, cache_dir=tmp_path)
    second.base_url = 'https://another.invalid'
    assert first.get('/same') == {'id': 1}
    assert second.get('/same') == {'id': 2}
