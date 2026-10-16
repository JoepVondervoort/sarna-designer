"""Replay an exported design's core calculation offline: python -m sarna.app.replay run.json."""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from sarna.app.service import candidate_record
from sarna.design.config import DesignConfig
from sarna.design.pipeline import design_sarnas
from sarna.genome.coordinates import Assembly, GenomicInterval, Strand
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.reproducibility import environment_identity, json_hash, software_identity
from sarna.rules.resolved import ResolvedRule
from sarna.rules.rule import Rule
from sarna.sequence.nucleic import NucleicSequence


def replay(snapshot: dict[str, Any]) -> dict[str, Any]:
    """Validate recorded inputs and compare a network-free replay with the saved candidates."""
    manifest = snapshot.get('reproducibility')
    if not isinstance(manifest, dict) or manifest.get('schema_version') != 1:
        raise ValueError('No supported reproducibility manifest in this record. '
                         'Historical provenance cannot be reconstructed reliably.')
    for name in ('promoter', 'configuration', 'candidates', 'counts'):
        if json_hash(snapshot[name]) != manifest.get(f'{name}_sha256'):
            raise ValueError(f'The saved {name} content does not match its recorded fingerprint.')
    for name, field in (('request', 'request'), ('reference', 'anchor')):
        if json_hash(snapshot['app'][field]) != manifest.get(f'{name}_sha256'):
            raise ValueError(f'The saved {name} does not match its recorded fingerprint.')
    data = snapshot['promoter']
    if hashlib.sha256(data['sequence'].encode()).hexdigest() != manifest.get('sequence_sha256'):
        raise ValueError('The promoter sequence does not match its recorded fingerprint.')
    interval = None
    if data['assembly'] is not None:
        interval = GenomicInterval(Assembly(data['assembly']), data['contig'],
                                   data['start'], data['end'], Strand(data['strand']))
    promoter = Promoter(
        NucleicSequence.parse(data['sequence']),
        TssPosition(data['tss_start']) if data['tss_start'] is not None else None, interval,
    )
    resolved = []
    for row in snapshot['configuration']['rules']:
        # Reconstruct frozen rule declarations and values, never today's preset or catalogue.
        rule = Rule.model_validate(row)
        values = {parameter['name']: parameter['value'] for parameter in row['parameters']}
        resolved.append(ResolvedRule(rule, overrides=values))
    config = DesignConfig(tuple(resolved))
    result = design_sarnas(promoter, config)
    thermo = next((config.rule(key) for key in (
        'thermodynamic-asymmetry', 'asymmetry-margin',
    ) if key in config), None)
    model = thermo.text('nearest-neighbour-table') if thermo else 'rna-xia-1998'
    end_length = thermo.integer('end-length-nt') if thermo else 4
    ranks = {candidate.id: candidate for candidate in result.ranked}
    candidates = [candidate_record(ranks.get(c.id, c), model, end_length)
                  for c in result.generated]
    counts = {'generated': len(result.generated), 'passing': len(result.passing),
              'rejected': len(result.rejected)}
    changed = [old['id'] for old, new in zip(snapshot['candidates'], candidates, strict=False)
               if old != new]
    current = software_identity()
    environment = environment_identity()
    return {
        'run_id': snapshot['app']['run_id'], 'scope': manifest['replay_scope'],
        'input_fingerprints_match': True,
        'output_matches': (json_hash(candidates) == manifest['candidates_sha256']
                           and counts == snapshot['counts']),
        'saved_candidates': len(snapshot['candidates']), 'replayed_candidates': len(candidates),
        'changed_candidates': changed,
        'software_matches': current['package_sha256'] == manifest['software']['package_sha256'],
        'environment_matches': environment == manifest['environment'],
        'current_software': current, 'limitations': manifest['limitations'],
        'note': 'Matching output is computational agreement, not biological validation. '
                'Software and environment agreement are reported separately.',
    }


def main() -> int:
    """Print a comparison report; exit 0 for matching output, 1 for drift, 2 for invalid input."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path, help='Whole-run JSON exported by the workbench')
    args = parser.parse_args()
    try:
        snapshot = json.loads(args.record.read_text(encoding='utf-8'))
        if not isinstance(snapshot, dict):
            raise ValueError('The exported record must be a JSON object.')
        report = replay(snapshot)
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps({'error': str(error)}))
        return 2
    print(json.dumps(report, indent=2))
    return 0 if report['output_matches'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
