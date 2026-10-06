# saRNA Designer

An internal research application for designing promoter-targeting **small activating RNAs**
(saRNAs), taking a researcher from a gene symbol or a pasted promoter sequence through to
ranked, annotated, exportable duplex candidates.

> **Status: Phase 0 — discovery and requirements.**
> No design engine exists yet. What is in this repository today is the legacy-script inventory,
> the domain glossary, the architecture decisions made so far, and the preserved original
> scripts. The architecture and technology stack are **not yet decided**. Nothing here produces
> scientific results. Do not use it to design anything.

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
- Whether the AF22 histone bigWigs (GSE108506) are on the same assembly as Ensembl's default
  (GRCh38) is **unverified**. An assembly mismatch here would invalidate every regulatory
  annotation and is a blocking open question.
- `Designs.xlsx`, the empirical design list that would be the most valuable test fixture
  available, is missing and cannot be supplied.

## Repository layout

```
sarna-designer/
├── GLOSSARY.md              Canonical project vocabulary — read this first
├── legacy/                  Original scripts, preserved verbatim. Never edited.
├── docs/
│   ├── adr/                 Architecture decision records
│   ├── scientific/          Legacy inventory, rule provenance, validation records
│   ├── architecture/        Component and data-flow documentation
│   └── development/         Environment and workflow notes
├── tests/                   Unit, regression, strand and integration tests
└── pyproject.toml           Tooling config (no package metadata until stack is decided)
```

## Setup (WSL)

Prerequisites, verified versions and the rationale for keeping this repository on the Linux
filesystem rather than `/mnt/c` are in
[`docs/development/environment.md`](docs/development/environment.md).

Installation, development, test, build and run commands will be documented here once the stack
is chosen. Documenting them now would mean inventing them.

## Licence

**Undecided.** Default copyright applies — see [`LICENSE`](LICENSE).
