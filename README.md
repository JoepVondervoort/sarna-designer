# saRNA Designer

**An auditable research workbench for promoter-targeting small activating RNA design.**

saRNA Designer takes a human gene symbol or promoter sequence through candidate generation,
literature-derived design checks, off-target screening and export of RNA duplexes for ordering.
It connects sequence design with the evidence needed to explain a shortlist: transcript choice,
genomic coordinates, rule provenance, strand-loading predictions and screening records.

Developed by **Joep van de Vondervoort**, the project brings saRNA biology, critical evaluation
of design assumptions and reproducible research software into one workflow.

[![Checks](https://github.com/JoepVondervoort/sarna-designer/actions/workflows/checks.yml/badge.svg)](https://github.com/JoepVondervoort/sarna-designer/actions/workflows/checks.yml)
![Python 3.11](https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-087F8C)](LICENSE)

[Installation](INSTALL.md) · [Scientific evidence](docs/scientific/evidence-guide.md) ·
[Empirical evaluation](docs/scientific/empirical-validation.md) ·
[Architecture decisions](docs/adr/) · [Screenshot gallery](docs/assets/screenshots/README.md)

![SCN1A promoter landscape with transcript reference, design windows, candidate positions and exclusion reasons](docs/assets/screenshots/promoter-landscape.png)

*SCN1A, using cached Ensembl promoter data and the default Wang 2015 configuration.
The interface retains excluded targets and explains the available checks. Screenshots show
software outputs, not experimentally confirmed activation.*

## Why this project exists

Promoter-targeting saRNA design involves decisions that a sequence score alone cannot capture.
A gene can have multiple transcription start sites; design recommendations have different
evidential support; strand and assembly mistakes change the interpretation of a target.
A candidate that passes a checklist still needs experimental validation.

This project replaces disconnected design and annotation scripts with a workflow in which
those decisions are explicit, inspectable and saved with the result.

## From target to research shortlist

| Stage | What the workbench provides |
|---|---|
| **Define the promoter** | Human gene lookup through Ensembl, MANE Select or explicit reference selection, assembly-tagged coordinates and transcript-oriented sequence. Pasted DNA is also supported. |
| **Generate and inspect** | Candidate enumeration, GC and sequence criteria, CpG-island checks where context permits, TSS-relative positions and an interactive promoter map. |
| **Explain the design** | Per-rule citations, rationale and verdicts; mandatory gates separated from preferences; editable settings with recorded departures from the baseline. |
| **Inspect the duplex** | Guide/passenger predictions, terminal thermodynamics, antiparallel strand display and 3′ dTdT overhangs for the project's order chemistry. |
| **Review off-target evidence** | Local human genome and transcript screening against pinned references, with optional remote NCBI BLAST cross-checks. Missing and incomplete evidence remain visible. |
| **Export and reproduce** | Order CSV, screening summaries, Excel, self-contained HTML and whole-run JSON; offline replay of the core design calculation. |

![Candidate inspection showing paired RNA strands, terminal energies, guide prediction and individual rule verdicts](docs/assets/screenshots/duplex-inspection.png)

*Candidate-level inspection connects the proposed duplex to the calculations and design rules
behind it. Strand-loading predictions do not establish activation.*

## Scientific approach and current evidence

The default configuration implements the recommendations in
[Wang et al., 2015](https://doi.org/10.14440/jbm.2015.39).
Mandatory criteria determine inclusion; surviving candidates are ordered by the number of
preferences met, with thermodynamic asymmetry as a tie-break. Equal preference weighting is a
project choice, not a published efficacy model. Optional siRNA-derived criteria and exploratory
hypotheses are distinguished from published saRNA guidance.

The bundled evaluation contains **125 assayed sequences, including 14 labelled activators**.
It makes limitations of the design assumptions measurable:

- The sequence-only preference score **has not demonstrated discrimination above chance**.
  This evaluation does not validate the complete production ranking with its thermodynamic tie-break.
- Known activators can fail published gates: the CEBPA example fails the position-19 criterion,
  and the default location window excludes 10 of 14 labelled activators using recovered locations.
- The dataset is small, heterogeneous and restricted by the original design process.
  It can reveal exclusions of known activators, but cannot establish performance on unselected
  candidates or independently validate hypotheses derived from the same data.

**Use rank to prioritise investigation, not to predict activation.** Off-target alignments
also do not establish biological specificity. The next scientific validation step is an
independent, prospectively defined evaluation across relevant targets and cell contexts.

Read the [empirical evaluation](docs/scientific/empirical-validation.md) for methods and label
caveats, and the [evidence guide](docs/scientific/evidence-guide.md) for interpretation.

## Engineering decisions that protect scientific interpretation

- **Coordinate correctness:** explicit sequence, TSS and genomic coordinate conventions,
  with strand-aware conversion and assembly boundaries covered by regression tests.
- **Traceable decisions:** citations and rationale travel with rules; custom roles and parameter
  changes are recorded rather than silently treated as published recommendations.
- **Explicit unknowns:** unavailable context differs from a negative finding. Sequence motif
  predictions, tissue expression and chromatin observations remain distinct evidence axes.
- **Reproducible calculations:** saved inputs, resolved configurations and software fingerprints
  support offline replay. Annotation enrichment and screening records are preserved, not replayed.
- **A reusable Python core:** the design engine is independent of the browser interface and
  network services, so it can also run in scripts or notebooks.

The [architecture decision records](docs/adr/) explain the trade-offs.
The [legacy inventory](docs/scientific/legacy-inventory.md) records defects found in the
original scripts, which remain preserved under [`legacy/`](legacy/).

## Run locally

**Status: local research workbench. Python 3.11.** The repository currently requires GitHub
read access. New gene retrieval needs internet; cached inputs support offline work.

With [uv](https://docs.astral.sh/uv/) and Git installed:

```bash
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
sarna-designer-check
sarna-designer
```

Open **http://localhost:8001/** and keep the terminal running.
For a GitHub installation, update with `uv tool upgrade sarna-designer`.

The [installation guide](INSTALL.md) covers Windows, macOS, Linux and troubleshooting.
The [Windows GitHub walkthrough](docs/install-via-github.md) includes first-time authentication
and local off-target reference setup. Colleagues without repository access can install a supplied
`sarna_designer-0.2.0-py3-none-any.whl`; see the [wheel instructions](INSTALL.md#step-3-install-sarna-designer).

Local screening references are a separate, explicit download of about 1.2 GB and need at least
20 GiB of free disk space during preparation. Screening covers full-length 19-base cores with
0–3 mismatches; gaps, bulges, shorter seed matches and personal variants are outside its scope.
Native Windows screening and fresh Windows/macOS installations still need platform validation.

## Use the Python engine

```python
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.promoter import Promoter
from sarna.sequence.nucleic import NucleicSequence

with open("promoter.txt") as handle:
    promoter = Promoter(sequence=NucleicSequence.parse(handle.read()))

result = design_sarnas(promoter, WANG_2015.resolve())
for candidate in result.ranked[:10]:
    print(candidate.id, candidate.score_summary)
    print(candidate.antisense_strand)
```

A pasted sequence without a TSS or genomic context leaves location-dependent checks unavailable.
The run retains rejected candidates and their reasons. See the
[workbench guide](docs/development/local-workbench.md) and
[reproducibility guide](docs/development/reproducibility.md) for the full workflow.

## Development and verification

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -e ".[app,dev]"
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m mypy src tests scripts
```

GitHub Actions runs the default test suite, lint and strict type checks. Live-network tests are
excluded by default; optional AF22 track tests skip when their data are absent. Browser checks
exercise the app separately, including saved runs, screening displays, exports and mobile layout.
See [screenshot provenance](docs/assets/screenshots/README.md) and the
[development environment](docs/development/environment.md).

| Directory | Purpose |
|---|---|
| [`src/sarna/`](src/sarna/) | Design engine, sequence and coordinate models, evidence adapters and local app |
| [`tests/`](tests/) | Scientific invariants, regression tests, application and export checks |
| [`docs/scientific/`](docs/scientific/) | Evidence interpretation, evaluation methods and findings |
| [`docs/adr/`](docs/adr/) | Scientific and architectural decisions |
| [`scripts/`](scripts/) | Evaluation, browser verification and research utilities |
| [`legacy/`](legacy/) | Original scripts preserved for comparison |

## Licence and attribution

Code: [Apache License 2.0](LICENSE). Copyright 2026 Joep van de Vondervoort.
Bundled assets retain the terms recorded in [NOTICE](NOTICE), including attribution for the
JASPAR CORE subset. Cite the underlying sources when using their data or design recommendations.
