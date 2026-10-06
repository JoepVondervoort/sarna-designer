# saRNA Designer

An internal research application for designing promoter-targeting **small activating RNAs**
(saRNAs), taking a researcher from a gene symbol or a pasted promoter sequence through to
ranked, annotated, exportable duplex candidates.

> **Status: Phase 1 — foundations.**
> The stack is now decided ([ADR 0005](docs/adr/0005-python-core-with-thin-web-ui.md)) and the
> genome-coordinate layer exists. **No design engine, no rules and no scoring exist yet.**
> Nothing here produces scientific results. Do not use it to design anything.

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
- **`legacy-scanner`** — reproduces the original script bit-for-bit, including its known
  deviations and defects, so historical designs remain reproducible.

Every design rule carries its **source citation, evidence class and rationale**, visible in the
UI ([ADR 0003](docs/adr/0003-rule-provenance-is-first-class.md)). Rules with no traceable
source are labelled `unsourced` rather than quietly presented as established biology.

### Honest expectations

Wang et al. report that of designed duplexes, **15–20 % significantly activate** the target
gene, **~70 % have no effect**, and **5–10 % suppress** it. A high rank from this tool is a
prioritisation aid, not a prediction of activity.

## Current limitations

- Nothing is implemented. See status above.
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
│   └── genome/              Assembly-tagged coordinates and reference contig data
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
