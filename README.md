# saRNA Designer

An internal research application for designing promoter-targeting **small activating RNAs**
(saRNAs), taking a researcher from a gene symbol or a pasted promoter sequence through to
ranked, annotated, exportable duplex candidates.

> **Status: Phase 2 — the design engine works; there is no interface yet.**
> `design_sarnas(promoter, config)` generates, evaluates, gates and ranks candidates from a
> pasted promoter or an assembly-anchored one, on either strand, and explains every verdict.
> **There is no gene lookup, no API and no UI** — nothing resolves a gene symbol to a promoter
> yet, so Workflow A does not exist. The engine has not been validated against any empirical
> design list, because none is available. Treat its output as a prioritisation aid under review,
> not as a validated result.

## What problem this solves

The existing workflow is two unconnected Python scripts: one scans a pasted promoter string for
19-mer candidates, the other locates already-chosen candidates in a genome and checks histone
marks around them. Neither knows about the other. Between them they disagree on the design
rules, the scientific parameters are hardcoded module constants, there is no concept of a
transcript or a genome assembly in the design path, and minus-strand coordinates are computed
incorrectly.

This project turns that into a maintainable tool that a colleague without programming
experience can use, and whose answers a reviewer can audit.

## Scientific scope

Design rules derive from **Wang J, Place RF, Portnoy V, Huang V, Kang MR, Kosaka M, Ho MKC,
Li L-C. _Inducing gene expression by targeting promoter sequences using small activating RNAs._
J Biol Methods, 2015** ([PMC4379447](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379447/)),
whose Table 1 is the canonical criteria set.

Two rule presets ship (see [ADR 0002](docs/adr/0002-wang-2015-default-legacy-preset.md)):

- **`wang-2015`** — the default. Conformant with Table 1.
- **`legacy-scanner`** — reproduces the original script's *filtering* closely, defects included,
  so historical designs can be re-derived. **Not bit-for-bit**, and the two known divergences are
  declared on the preset itself so they appear in exports and the UI: position 1 G/C is scored as
  a preference where the original required it (so this preset admits candidates the original
  rejected), and ranking is by criteria count rather than the original's weighted sum. See
  [ADR 0010](docs/adr/0010-presets-are-selections-with-overrides.md).

Every design rule carries its **source citation, evidence class and rationale**, visible in the
UI ([ADR 0003](docs/adr/0003-rule-provenance-is-first-class.md)). Rules with no traceable
source are labelled `unsourced` rather than quietly presented as established biology.

### How candidates are ranked

Mandatory criteria are a **gate**: fail one and the candidate is excluded, with the failing rule
named. Everything else is a **checklist of preferences, each worth one point, all weighted
equally** — a candidate reads as "meets 5 of 7 preferences", itemised. Ties break on
thermodynamic asymmetry, reported as a tie-break rather than folded into the score. See
[ADR 0006](docs/adr/0006-criteria-count-ranking.md).

This follows the structure of Wang 2015 Table 1, which separates requirements from criteria
marked "preferably". It is **not** a formula taken from the paper: the paper publishes no
weights and no ranking formula at all. Equal weighting is a deliberate choice, because the
evidence gives no basis for ranking the preferences against each other.

### Honest expectations

Wang et al. report that of designed duplexes, **15–20 % significantly activate** the target
gene, **~70 % have no effect**, and **5–10 % suppress** it. A high rank from this tool is a
prioritisation aid, not a prediction of activity.

## Using the engine

```python
from sarna.design.pipeline import design_sarnas
from sarna.design.presets import WANG_2015
from sarna.genome.promoter import Promoter
from sarna.genome.tss import TssPosition
from sarna.sequence.nucleic import NucleicSequence

promoter = Promoter(
    sequence=NucleicSequence.parse(open("promoter.txt").read()),
    tss_start=TssPosition(-1000),          # optional; omit for a bare sequence
)
result = design_sarnas(promoter, WANG_2015.resolve())

for candidate in result.ranked[:10]:
    print(candidate.id, candidate.tss_span, candidate.score_summary)
    print(" guide 5'-", candidate.antisense_strand, "-3'")
    print(" ", candidate.guide_comparison.summary)

rejected = result.rejected[0]
print(rejected.id, rejected.rejection_summary)   # never silently dropped
```

`result.to_dict()` is a JSON-serialisable record of the run, including the resolved
configuration, every candidate and every rule verdict.

## Current limitations

- **No gene lookup.** Workflow A (gene symbol to promoter) is not implemented, so a promoter must
  be supplied as sequence.
- **No API or UI.** The engine is a library only.
- **No empirical validation.** The engine reproduces the rules as published and as scripted; it
  has not been checked against measured activation data, and none is available to us.
- Promoter fixtures must be CpG-depleted or real sequence. Uniform random sequence is called as
  one large CpG island by any correct detector, which rejects nearly every candidate — see
  [`docs/scientific/synthetic-promoter-fixtures.md`](docs/scientific/synthetic-promoter-fixtures.md).
- The legacy scripts have **14 catalogued defects**, three scientifically material. They are
  preserved unfixed by design; see [`legacy/README.md`](legacy/README.md).
- The AF22 histone bigWigs (GSE108506) are on **GRCh37/hg19**, while Ensembl REST serves
  **GRCh38** by default. Verified 2026-10-06 from the bigWig headers. The legacy annotation
  pipeline combined the two without converting, so **every regulatory annotation it produced is
  unreliable** (defect B15). How the new engine reconciles the assemblies is an open
  architecture decision; the design path is unaffected.
- `Designs.xlsx`, the empirical design list that would be the most valuable test fixture
  available, is missing and cannot be supplied.

## Repository layout

```
sarna-designer/
├── GLOSSARY.md              Canonical project vocabulary — read this first
├── legacy/                  Original scripts, preserved verbatim. Never edited.
├── src/sarna/               The engine. No HTTP, no UI, no web framework (ADR 0005).
│   ├── sequence/            Validated sequence primitives and CpG-island detection
│   ├── genome/              Assembly-tagged coordinates, TSS positions, promoters
│   ├── rules/               The rule catalogue, parameters, provenance and results
│   ├── thermo/              Nearest-neighbour tables, duplex ends, guide prediction
│   └── design/              Contexts, evaluators, presets and the design pipeline
├── docs/
│   ├── adr/                 Architecture decision records
│   ├── scientific/          Legacy inventory, rule provenance, validation records
│   ├── architecture/        Component and data-flow documentation
│   └── development/         Environment and workflow notes
├── tests/                   Unit, regression, strand and integration tests
└── pyproject.toml           Package metadata and tooling configuration
```

## Setup (WSL)

Prerequisites, verified versions and the rationale for keeping this repository on the Linux
filesystem rather than `/mnt/c` are in
[`docs/development/environment.md`](docs/development/environment.md).

Python 3.11 in an isolated virtualenv, so the miniconda base environment is never mutated:

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -e ".[dev,genomics]"
```

The `genomics` extra pulls in `pyBigWig`, `httpx` and `pyliftover`. The engine itself installs
without it, which is what keeps the rules testable with no network and no 890 MB of tracks
present.

```bash
.venv/bin/python -m pytest                      # full suite; network tests excluded by default
.venv/bin/python -m pytest -m strand            # the coordinate-correctness gate
.venv/bin/python -m pytest --cov --cov-report=term-missing
.venv/bin/python -m mypy src tests              # strict
.venv/bin/python -m ruff check .
```

Tests that need the AF22 bigWigs skip cleanly when the files are absent; point
`SARNA_AF22_DIR` at them if they live somewhere other than `~/af22_data`.

Build and run commands will be added with the API and frontend layers, which do not exist yet.

## Licence

**Undecided.** Default copyright applies — see [`LICENSE`](LICENSE).
