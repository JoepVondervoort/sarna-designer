"""Concurrent source writes must use distinct temporary files and publish complete JSON."""

import json
import os
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from typing import Any

import pytest

from sarna.genome.ensembl import EnsemblClient
from sarna.sources.json_http import JsonHttpCache


@pytest.mark.parametrize('writer', [JsonHttpCache._write_cache, EnsemblClient._write_cache])
def test_threads_do_not_share_a_temporary_cache_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    writer: Callable[[Path | None, Any], None],
) -> None:
    target = tmp_path / 'response.json'
    barrier = Barrier(2)
    replace = os.replace
    temporary_paths: list[Path] = []

    def together(source: Path, destination: Path) -> None:
        temporary_paths.append(source)
        barrier.wait(timeout=5)
        replace(source, destination)

    monkeypatch.setattr(os, 'replace', together)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda value: writer(target, {'value': value}), [1, 2]))
    assert len(set(temporary_paths)) == 2
    assert json.loads(target.read_text()) in ({'value': 1}, {'value': 2})
    assert list(Path(tmp_path).glob('*.part')) == []
