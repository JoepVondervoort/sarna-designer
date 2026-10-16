"""Researcher notes must persist without mutating science or losing concurrent edits."""

import json
from pathlib import Path
from typing import Any

import pytest

from sarna.app.store import RunStore
from sarna.app.web import AppRouter
from tests.unit.test_app import run_request, service_for
from tests.unit.test_workbook import sheets_of


def test_shortlist_persists_separately_and_rejects_stale_revision(tmp_path: Path) -> None:
    service = service_for(tmp_path)
    saved = run_request(service)
    run_id = saved['app']['run_id']
    candidate = saved['candidates'][0]['id']
    router = AppRouter(service)
    try:
        path = f'/api/runs/{run_id}/shortlist'
        body = {'revision': 0, 'entries': [
            {'candidate': candidate, 'reason': 'Review missing context'},
        ]}
        code, _, raw = router.dispatch('POST', path, json.dumps(body).encode())
        assert code == 200
        assert json.loads(raw)['revision'] == 1
        code, _, raw = router.dispatch('POST', path, json.dumps(body).encode())
        assert code == 409
        assert json.loads(raw)['code'] == 'shortlist-conflict'
        assert service.store.snapshot(run_id) == saved
        reopened = RunStore(service.store.path)
        notes = reopened.shortlist(run_id)
        assert notes['entries'] == body['entries']
        code, _, raw = router.dispatch('GET', f'/api/runs/{run_id}/export.json')
        exported = json.loads(raw)
        assert exported.pop('workspace')['shortlist'] == notes
        assert exported == saved
        code, _, raw = router.dispatch('GET', f'/api/runs/{run_id}/export.html')
        assert code == 200
        assert b'Review missing context' in raw
        assert b'<script src="/workflow.js"' not in raw
        code, _, raw = router.dispatch('GET', f'/api/runs/{run_id}/export.xlsx')
        assert code == 200
        assert sheets_of(raw)['Shortlist'][1][0:2] == [candidate, 'Review missing context']
        code, _, raw = router.dispatch('POST', path, b'{"revision":1,"entries":[]}')
        assert code == 200 and json.loads(raw)['entries'] == []
    finally:
        router.close()


@pytest.mark.parametrize('body', [
    {'revision': True, 'entries': []},
    {'revision': -1, 'entries': []},
    {'revision': 0, 'entries': [{'candidate': 'not-in-this-run', 'reason': ''}]},
    {'revision': 0, 'entries': [{'candidate': 'saRNA-001', 'reason': 'x' * 2001}]},
    {'revision': 0, 'entries': [{'candidate': 'saRNA-001'}] * 2},
])
def test_invalid_shortlist_never_changes_saved_notes(tmp_path: Path, body: dict[str, Any]) -> None:
    service = service_for(tmp_path)
    saved = run_request(service)
    run_id = saved['app']['run_id']
    router = AppRouter(service)
    try:
        code, _, _ = router.dispatch(
            'POST', f'/api/runs/{run_id}/shortlist', json.dumps(body).encode(),
        )
        assert code == 422
        assert service.store.shortlist(run_id)['revision'] == 0
    finally:
        router.close()


def test_history_search_reaches_runs_beyond_the_first_hundred(tmp_path: Path) -> None:
    store = RunStore(tmp_path / 'runs.sqlite3')
    ids = []
    for index in range(105):
        run_id = store.create({'symbol': f'GENE{index}', 'transcript_id': f'reference-{index}'})
        store.fail(run_id, 'Fixture')
        ids.append(run_id)
    assert len(store.list_runs()) == 100
    assert [row['id'] for row in store.list_runs(search='reference-0')] == [ids[0]]
    assert [row['id'] for row in store.list_runs(offset=100, limit=50)] == ids[4::-1]
    assert store.list_runs(status='complete') == []
    assert store.list_runs(search="%' OR 1=1 --") == []
    with pytest.raises(ValueError):
        store.list_runs(limit=1000)
