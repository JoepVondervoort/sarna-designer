"""Content identities and run manifests; no network or scientific evaluation occurs here."""

import hashlib
import json
import platform
import subprocess
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, distributions, version
from pathlib import Path
from typing import Any

PACKAGE_ROOT = Path(__file__).resolve().parent
IDENTITY_FILE = '_build_identity.json'
ASSET_SUFFIXES = {'.py', '.json', '.csv', '.toml', '.html', '.js', '.css', '.svg'}


def json_hash(value: Any) -> str:  # noqa: ANN401
    """Hash canonical UTF-8 JSON with sorted keys and explicit finite numbers."""
    encoded = json.dumps(value, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def package_files(root: Path) -> dict[str, str]:
    """Identify shipped code and reference assets, excluding generated build metadata."""
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob('*'))
        if path.is_file() and path.suffix in ASSET_SUFFIXES
        and path.name != IDENTITY_FILE and '__pycache__' not in path.parts
    }


def git_identity(root: Path) -> dict[str, Any]:
    """Read only the checkout containing this package, never an unrelated working directory."""
    checkout = root.parent.parent
    unknown = {'git_revision': None, 'git_dirty': None}
    if not (checkout / '.git').exists() or checkout / 'src' / 'sarna' != root:
        return unknown
    try:
        def git(*args: str) -> str:
            return subprocess.run(
                ['git', '-C', str(checkout), *args], check=True, capture_output=True,
                text=True, timeout=5,
            ).stdout.strip()
        return {'git_revision': git('rev-parse', 'HEAD'),
                'git_dirty': bool(git('status', '--porcelain', '--untracked-files=normal'))}
    except (OSError, subprocess.SubprocessError):
        return unknown


def software_identity(root: Path = PACKAGE_ROOT) -> dict[str, Any]:
    """Report exact content plus checkout or packaged build provenance when available."""
    files = package_files(root)
    build: dict[str, Any] = {}
    stamp = root / IDENTITY_FILE
    if stamp.exists():
        try:
            decoded = json.loads(stamp.read_text(encoding='utf-8'))
            if isinstance(decoded, dict):
                build = decoded
        except (OSError, ValueError):
            pass
    git = git_identity(root)
    if git['git_revision'] is None:
        git = {key: build.get(key) for key in ('git_revision', 'git_dirty')}
    return {
        **git, 'package_sha256': json_hash(files), 'files': files,
        'build_package_sha256': build.get('package_sha256'),
        'matches_packaged_build': (
            json_hash(files) == build['package_sha256'] if 'package_sha256' in build else None
        ),
    }


def environment_identity() -> dict[str, Any]:
    """Record installed versions without local paths, repository URLs or credentials."""
    installed = sorted(
        (distribution.metadata['Name'] or '', distribution.version)
        for distribution in distributions()
    )
    return {'python': platform.python_version(), 'implementation': platform.python_implementation(),
            'platform': platform.system(), 'machine': platform.machine(),
            'packages': {name: value for name, value in installed if name}}


def source_observation(url: str, payload: object, cached: bool) -> dict[str, Any]:
    """Identify the returned response; a legacy cache has no known original retrieval date."""
    now = datetime.now(UTC).isoformat()
    return {'source': 'Ensembl REST', 'url': url, 'response_sha256': json_hash(payload),
            'observed_at': now, 'cache_hit': cached,
            'retrieved_at': None if cached else now, 'annotation_release': None,
            'note': 'Annotation release not recorded; a stable ID does not pin a release.'}


def run_manifest(snapshot: dict[str, Any], sources: list[dict[str, Any]],
                 software: dict[str, Any]) -> dict[str, Any]:
    """Freeze input identities and exact core-engine output, with explicit replay boundaries."""
    try:
        engine_version = version('sarna-designer')
    except PackageNotFoundError:
        engine_version = '0+unknown'
    return {
        'schema_version': 1, 'engine_version': engine_version, 'software': software,
        'environment': environment_identity(),
        'sequence_sha256': hashlib.sha256(snapshot['promoter']['sequence'].encode()).hexdigest(),
        'promoter_sha256': json_hash(snapshot['promoter']),
        'configuration_sha256': json_hash(snapshot['configuration']),
        'candidates_sha256': json_hash(snapshot['candidates']),
        'counts_sha256': json_hash(snapshot['counts']),
        'request_sha256': json_hash(snapshot['app']['request']),
        'reference_sha256': json_hash(snapshot['app']['anchor']),
        'sources': sources,
        'replay_scope': 'Core candidate generation, rule verdicts, guide models and ranking.',
        'limitations': [
            'Core replay uses saved sequence and resolved rules; it never refreshes references.',
            'Remote annotation releases and original retrieval dates of legacy caches are unknown.',
            'Chromatin, expression, motifs and BLAST are preserved outputs, not replayed.',
            'Recorded package versions identify an environment; they are not a portable lockfile.',
            'Hashes detect changed content; they are not signatures or biological validation.',
        ],
    }
