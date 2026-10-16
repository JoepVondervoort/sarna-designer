# Reproducibility of a saved design

A saved result can be reopened without recalculating it. A core replay goes further: it repeats
candidate generation, rule evaluation, guide-model calculations and ranking from the saved
promoter and resolved rule definitions. Neither action establishes biological efficacy.

## What new runs record

The whole-run JSON contains a `reproducibility` manifest:

- The package version, Git commit when available, and whether the source checkout had changes.
- SHA-256 identities for shipped package files, including code and bundled reference assets.
  A dirty Git revision alone cannot identify edited code; the content fingerprint distinguishes it.
- Python, operating-system and installed distribution versions, without machine-specific paths.
- Fingerprints of the promoter sequence, promoter coordinates, resolved configuration, request,
  selected reference, candidate records and counts.
- Each successful Ensembl response used in the run: request URL, canonical JSON response hash,
  observation time and cache/live status. Failed lookups remain run notices rather than evidence.

Canonical JSON hashing sorts object keys, preserves array order, uses UTF-8 and compact separators,
and rejects non-finite numbers. Sequence hashing uses the normalized DNA stored in the snapshot.
The hashes identify content; they are not signatures or an independent check of source truth.

The app exposes a reproducibility disclosure beneath the ranking explanation. JSON and offline
HTML retain the manifest. Excel includes identifying fields on its Run sheet and response metadata
on a Reference sources sheet. Candidate CSVs remain sequence/order formats; keep the whole-run
JSON for replay.

## Run a core replay

Export **JSON** from the workbench, then use the same app installation:

```bash
sarna-designer-replay run.json
```

From a source checkout:

```bash
.venv/bin/python -m sarna.app.replay run.json
```

Replay checks the recorded input/output fingerprints, reconstructs the frozen rules and their
values, and runs against the saved sequence. It never looks up the gene or resolves today's
preset. A changed reference service therefore cannot alter the input to this replay.

The printed report distinguishes `output_matches`, `software_matches` and `environment_matches`.
An output match on different software is computational agreement on this input, not an identical
execution environment or validation of either implementation. Exit codes are 0 for matching
output, 1 for output drift, and 2 for a missing/invalid manifest or input that cannot be replayed.
The original record is never modified.

To make a fresh-reference analysis, submit a new gene-based design in the workbench. That is a
new run, not a replay. Normal online mode can reuse cached responses; merely allowing network
access does not guarantee refreshed references. Do not overwrite a saved result to update it.

## Wheels and source archives

The Hatch build hook embeds `_build_identity.json` in wheels and source archives. It captures
Git metadata while the checkout is available. A wheel built from that source archive retains the
captured revision even on a machine without Git. Runtime hashing also reports whether installed
package content matches the packaged build fingerprint. Unknown revision stays unknown when
building from a source tree without Git or an existing build record.

Build as usual:

```bash
.venv/bin/python -m build
```

Keep the distributed wheel or source archive with the analysis. The fingerprint helps identify
code, but cannot recover a dirty working tree whose source files were not preserved. Update the
installed app after rebuilding to obtain the new replay command and manifest support.

## Boundaries and unknowns

The Ensembl adapter does not pin an annotation release. Stable identifiers and genome assembly
names do not freeze an annotation release. The manifest therefore reports an unknown annotation
release rather than inferring it. Cached responses have a recorded observation date, but their
original retrieval date is unknown because older cache files did not store it. A new live response
has a retrieval date recorded for this run; that is not a database publication date.

Core replay does not re-run chromatin, tissue expression, motif annotations, or BLAST. Their
saved results remain inspectable with the provenance available in their respective records.
Exact raw external response bodies and local bigWig/chain files are not bundled into this manifest.
The Ensembl response hashes identify the responses used; they do not make those bodies retrievable.

Installed package versions are an environment inventory, not a cross-platform dependency lock.
Recreating that environment still requires the recorded packages and a compatible Python/platform.
The report explicitly distinguishes output agreement from environment agreement.

Historical runs remain unchanged. A missing manifest is reported as not recorded; today's code
revision and data are never attributed retroactively to an old analysis.
