"""Map investigation preserves coordinates, context gaps and concurrent researcher edits."""

import json
from pathlib import Path
from typing import Any

import pytest

from sarna.app.investigation import SequenceSearch
from sarna.app.sequence_import import compare
from sarna.app.web import AppRouter
from tests.unit.test_app import run_request, service_for


def test_search_overlapping_rna_reverse_and_palindromic_matches():
    assert SequenceSearch(sequence="AAAA").find("AAAAA") == [
        {"start": 0, "end": 4, "orientation": "same"},
        {"start": 1, "end": 5, "orientation": "same"},
    ]
    assert SequenceSearch(sequence="UUUU", orientation="reverse complement").find("AAAA") == [
        {"start": 0, "end": 4, "orientation": "reverse complement"},
    ]
    assert len(SequenceSearch(sequence="ACGT").find("ACGT")) == 1
    assert SequenceSearch(sequence="CCCC").find("AAAA") == []
    with pytest.raises(ValueError, match="500"):
        SequenceSearch(sequence="AAAA").find("A" * 1000)
    with pytest.raises(ValueError, match="mixed"):
        SequenceSearch(sequence="AUTG").find("ATTT")


def test_investigation_roundtrip_conflict_validation_and_exports(tmp_path: Path):
    service = service_for(tmp_path)
    snapshot = run_request(service)
    run_id = snapshot['app']['run_id']
    path = f'/api/runs/{run_id}/investigation'
    router = AppRouter(service)
    try:
        draft: dict[str, Any] = {'revision': 0, 'annotations': [
            {'name': 'Literature region', 'start': 0, 'end': 4, 'source': '<citation>',
             'orientation': 'same'},
        ], 'motifs': []}
        code, _, raw = router.dispatch('POST', path, json.dumps(draft).encode())
        assert code == 200
        assert json.loads(raw)['revision'] == 1
        assert router.dispatch('POST', path, json.dumps(draft).encode())[0] == 409
        draft['revision'] = 1
        draft['annotations'][0]['end'] = snapshot['promoter']['length'] + 1
        assert router.dispatch('POST', path, json.dumps(draft).encode())[0] == 422
        assert service.store.snapshot(run_id) == snapshot
        assert service.store.investigation(run_id)['revision'] == 1
        for endpoint in ('export.json', 'export.html'):
            code, _, raw = router.dispatch('GET', f'/api/runs/{run_id}/{endpoint}')
            assert code == 200 and b'Literature region' in raw
        code, _, raw = router.dispatch('GET', f'/api/runs/{run_id}')
        assert json.loads(raw)['result']['workspace']['investigation']['motifs'] == []
        query = {'sequence': snapshot['promoter']['sequence'][:4], 'orientation': 'same'}
        code, _, raw = router.dispatch(
            'POST', f'/api/runs/{run_id}/sequence-search', json.dumps(query).encode(),
        )
        assert code == 200 and json.loads(raw)['matches'][0]['start'] == 0
    finally:
        router.close()


def test_unmatched_external_core_keeps_context_checks_unavailable(tmp_path: Path):
    snapshot = run_request(service_for(tmp_path))
    # A core absent from this fixture must not acquire made-up TSS/flank/CpG context.
    core = 'A' * 19
    assert core not in snapshot['promoter']['sequence']
    record = compare(snapshot, [{'sequence': core, 'role': 'sense'}])[0]
    assert record['matches'] == []
    assessment = record['sequence_only'][0]['assessment']
    assert assessment['rank'] is None
    assert assessment['tss_span'] is None and assessment['genomic'] is None
    assert any(r['outcome'] == 'not-applicable' for r in assessment['rules'])
    assert any(r['outcome'] == 'fail' for r in assessment['rules'])
    assert 'scope' in record['sequence_only'][0]
    unknown = compare(snapshot, [{'sequence': core, 'role': 'unknown'}])[0]
    assert {a['orientation'] for a in unknown['sequence_only']} == {
        'same', 'reverse complement',
    }
    assert compare(snapshot, [{'sequence': 'AAAA', 'role': 'sense'}])[0]['sequence_only'] == []


def test_motif_selection_accepts_scanned_names_only(tmp_path: Path):
    service = service_for(tmp_path)
    snapshot = run_request(service, motifs=True, motif_shuffles=2)
    path = f"/api/runs/{snapshot['app']['run_id']}/investigation"
    router = AppRouter(service)
    try:
        names = [m['motif'] for m in snapshot['motif_scans'][:2]]
        assert len(names) == 2
        draft = {'revision': 0, 'annotations': [], 'motifs': names}
        assert router.dispatch('POST', path, json.dumps(draft).encode())[0] == 200
        for motifs in (['UNKNOWN'], [names[0], names[0]]):
            draft.update(revision=1, motifs=motifs)
            assert router.dispatch('POST', path, json.dumps(draft).encode())[0] == 422
        assert service.store.investigation(snapshot['app']['run_id'])['motifs'] == names
    finally:
        router.close()
