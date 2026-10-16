"""Run identity and offline replay must expose missing provenance and computational drift."""

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import httpx
import pytest

from sarna.app.replay import replay
from sarna.app.web import AppRouter
from sarna.genome.ensembl import EnsemblClient
from sarna.reproducibility import json_hash, software_identity
from tests.unit.test_app import run_request, service_for
from tests.unit.test_ensembl import GENE_PAYLOAD, GENE_ROUTE, client_for
from tests.unit.test_workbook import sheets_of


@pytest.mark.parametrize('roles', [{}, {'tss-window': 'preference', 'tss-preferred-window': 'off'}])
def test_replay_uses_saved_inputs_without_http(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, roles: dict[str, str],
) -> None:
    saved = run_request(service_for(tmp_path), roles=roles, tss_start=-1000)

    def refuse(*args: Any, **kwargs: Any) -> None:  # noqa: ANN401
        raise AssertionError('Replay must not contact external services')

    monkeypatch.setattr(httpx, 'get', refuse)
    monkeypatch.setattr(EnsemblClient, 'gene', refuse)
    report = replay(saved)
    assert report['output_matches']
    assert report['software_matches']
    assert report['environment_matches']
    assert report['changed_candidates'] == []


def test_gene_replay_preserves_minus_strand_coordinates(tmp_path: Path) -> None:
    saved = run_request(service_for(tmp_path), source='gene', symbol='SCN1A')
    assert saved['promoter']['assembly'] == 'GRCh38'
    assert saved['promoter']['strand'] == '-'
    assert saved['reproducibility']['sources']
    assert replay(saved)['output_matches']


@pytest.mark.parametrize('field', ['promoter', 'configuration', 'candidates', 'counts'])
def test_modified_recorded_content_is_rejected(tmp_path: Path, field: str) -> None:
    saved = run_request(service_for(tmp_path))
    saved[field] = {}
    with pytest.raises(ValueError, match='fingerprint'):
        replay(saved)


def test_historical_runs_do_not_acquire_invented_provenance(tmp_path: Path) -> None:
    saved = run_request(service_for(tmp_path))
    del saved['reproducibility']
    with pytest.raises(ValueError, match='Historical provenance'):
        replay(saved)


def test_matching_output_does_not_claim_matching_software_or_environment(tmp_path: Path) -> None:
    saved = run_request(service_for(tmp_path))
    saved['reproducibility']['software']['package_sha256'] = 'older-code'
    saved['reproducibility']['environment']['python'] = 'older-python'
    report = replay(saved)
    assert report['output_matches']
    assert not report['software_matches']
    assert not report['environment_matches']


def test_changed_outputs_are_reported_without_rewriting_the_snapshot(tmp_path: Path) -> None:
    saved = run_request(service_for(tmp_path))
    saved['candidates'][0]['rank'] = 999
    saved['reproducibility']['candidates_sha256'] = json_hash(saved['candidates'])
    before = deepcopy(saved)
    report = replay(saved)
    assert not report['output_matches']
    assert saved['candidates'][0]['id'] in report['changed_candidates']
    assert saved == before


def test_source_identity_marks_legacy_cache_retrieval_date_unknown(tmp_path: Path) -> None:
    client = client_for({GENE_ROUTE: GENE_PAYLOAD})
    client.cache_dir = tmp_path
    client.gene('SCN1A')
    client.gene('SCN1A')
    live, cached = client.observations
    assert live['response_sha256'] == cached['response_sha256']
    assert live['retrieved_at']
    assert not live['cache_hit']
    assert cached['cache_hit']
    assert cached['retrieved_at'] is None
    assert cached['annotation_release'] is None


def test_package_identity_changes_with_content_and_survives_packaging(tmp_path: Path) -> None:
    (tmp_path / 'engine.py').write_text('first')
    first = software_identity(tmp_path)
    assert first['git_revision'] is None
    (tmp_path / '_build_identity.json').write_text(json.dumps({
        **first, 'git_revision': 'a' * 40, 'git_dirty': False,
    }))
    packaged = software_identity(tmp_path)
    assert packaged['package_sha256'] == first['package_sha256']
    assert packaged['matches_packaged_build']
    assert packaged['git_revision'] == 'a' * 40
    (tmp_path / 'engine.py').write_text('second')
    changed = software_identity(tmp_path)
    assert changed['package_sha256'] != first['package_sha256']
    assert not changed['matches_packaged_build']


def test_excel_carries_run_identity_and_reference_observations(tmp_path: Path) -> None:
    service = service_for(tmp_path)
    saved = run_request(service, source='gene', symbol='SCN1A')
    router = AppRouter(service)
    try:
        url = f"/api/runs/{saved['app']['run_id']}/export.xlsx"
        status, _, body = router.dispatch('GET', url, b'')
        assert status == 200
        workbook = sheets_of(body)
        fields = dict(workbook['Run'][1:])
        assert fields['package SHA-256'] == saved['reproducibility']['software']['package_sha256']
        assert len(workbook['Reference sources']) > 1
        assert 'unknown' in workbook['Reference sources'][1]
    finally:
        router.close()
