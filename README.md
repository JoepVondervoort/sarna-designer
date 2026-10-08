# saRNA Designer

An internal research application for designing promoter-targeting **small activating RNAs**
(saRNAs), taking a researcher from a gene symbol or a pasted promoter sequence through to
ranked, annotated, exportable duplex candidates.

> **Status: local research workbench.**
> A gene symbol resolves to a promoter through Ensembl, the engine designs and ranks candidates
> on either strand and explains every verdict, the promoter is annotated with CpG islands,
> regulatory features and TF motif predictions, and the whole run renders as one self-contained
> interactive HTML page. The local app adds gene lookup, browser-editable rules, Python reruns,
> saved runs, an off-target screen against NCBI, and exports including a spreadsheet of every
> candidate with its values.
>
> The engine is now validated against measured activity for the first time — 125 assayed
> sequences, 14 of them activators. The position-19 gate excludes CEBPA, a known activator;
> the default window excludes 10 of 14 when recovered locations are used.
> Read [`docs/scientific/empirical-validation.md`](docs/scientific/empirical-validation.md) for
> what that does and does not establish. Treat a high rank as a prioritisation aid, not a
> prediction of activity.

## Installing it — complete steps

Everything needed is below. [`INSTALL.md`](INSTALL.md) is the same thing with more hand-holding
and a troubleshooting section; this is the whole procedure.

### What you need first

| | Needed? | How |
|---|---|---|
| **PowerShell** | Yes, it is where you type | Already on Windows. Press <kbd>Win</kbd>, type `powershell`, Enter. |
| **`uv`** | Yes — the only install | One line, in step 1 below. |
| **Python** | No | `uv` downloads the right version (3.11) into its own folder. Your system Python is untouched. |
| **WSL** | **No** | Not used. See [the one exception](#optional-af22-histone-tracks) below. |
| **Administrator rights** | No | Everything goes under your own user account. |
| **Docker, Anaconda, Visual Studio** | No | None of them. |
| **Git** | Only for route B | [Git for Windows](https://git-scm.com/download/win) |

### Windows, start to finish

**1. Install `uv`** ([what this is](https://docs.astral.sh/uv/)). In PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

If your organisation blocks that, `winget install --id=astral-sh.uv -e` does the same job.

**2. Close PowerShell and open it again.** Required — the new command is only visible to a
fresh window. Check it worked:

```powershell
uv --version
```

**3. Install the tool.** This repository is **private**, so pick the route that applies:

*Route A — from a `.whl` file you were sent* (needs no GitHub account and no Git):

```powershell
uv tool install --python 3.11 "C:\Users\YourName\Downloads\sarna_designer-0.1.0.dev0-py3-none-any.whl[app]"
```

Tip: type `uv tool install --python 3.11 "`, drag the file from Explorer into the window to
paste its path, then type `[app]"` and Enter. The `[app]` is not optional — it installs the
browser interface.

*Route B — from this repository*, if you have been added to it and have
[Git for Windows](https://git-scm.com/download/win):

```powershell
uv tool install --python 3.11 "sarna-designer[app] @ git+https://github.com/JoepVondervoort/sarna-designer"
```

Either way it ends with `Installed 2 executables: sarna-designer, sarna-designer-check`.

**4. Close PowerShell and open it again** — the same PATH reason as step 2.

**5. Check it works:**

```powershell
sarna-designer-check
```

Ends with `RESULT: PASS`. Every failure prints what to do about it.

**6. Run it:**

```powershell
sarna-designer
```

Open the address it prints (`http://localhost:8001/`) in your browser, and leave the
PowerShell window open while you work. <kbd>Ctrl</kbd>+<kbd>C</kbd> stops it.

### macOS and Linux

Identical except step 1:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Optional: AF22 histone tracks

**This is the only part that needs WSL on Windows, and almost nobody needs it.** The
`[genomics]` extra reads local AF22 bigWig files, and it depends on
[pyBigWig](https://pypi.org/project/pyBigWig/), which publishes **no Windows wheel** and lists
only MacOS, POSIX and Unix as supported — it needs `libcurl` and `zlib` headers to build. The
design tool does not use it; chromatin simply reports "not loaded" without it.

If you do want those tracks on Windows, install
[WSL](https://learn.microsoft.com/windows/wsl/install) from an **administrator** PowerShell:

```powershell
wsl --install
```

Reboot, then follow the Linux steps above inside the Ubuntu window, adding `[app,genomics]`
instead of `[app]`, and set `SARNA_AF22_DIR` to the folder holding the bigWigs.

### Giving it to a colleague

Because the repository is private, send them a file rather than a link:

```bash
python -m build --wheel      # one ~320 kB .whl appears in dist/
```

They need no GitHub account, no Git and no Python — that is route A above.
[`INSTALL.md`](INSTALL.md#sending-it-to-a-colleague) has a message you can copy into the email,
and — if you develop inside WSL — how to reach the file from Windows at all, since a path
inside WSL is invisible to Outlook.

### Running it from a checkout instead

```bash
.venv/bin/python -m sarna.app      # then open http://localhost:8001/
```

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
whose Table 1 supplies the default criteria. These are literature recommendations rather than
universal biological requirements; researchers can investigate alternatives in custom configurations.

Four rule presets ship (see [ADR 0002](docs/adr/0002-wang-2015-default-legacy-preset.md) and
[ADR 0011](docs/adr/0011-optional-rules-are-hypotheses-not-findings.md)):

- **`wang-2015`** — the default. Conformant with Table 1.
- **`strict`** — tightens only the *location* criteria. The range-restricted dataset does not
  establish better efficacy for these settings. Two sequence tightenings were tried, measured and
  rejected: requiring zero CpG excludes 4 of the 14 known activators, and requiring an asymmetry
  margin of 0.5 kcal/mol excludes 2 more.
- **`exploratory`** — adds optional criteria from Reynolds et al. 2004, Ui-Tei et al. 2004 and the
  thesis's own positional analysis. Admits exactly the same candidates as the default, since every
  added rule is a preference; it reorders them. Its score denominator is 13, so its numbers are
  **not comparable** with any other preset's.
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

To use the adjustable app, follow the [local workbench guide](docs/development/local-workbench.md).
It retrieves supported human genes through Ensembl and accepts pasted DNA. New genes/windows need
internet; cached responses support offline runs. The earlier generated pages remain described in
the [preview guide](docs/development/trying-the-browser.md). An editable
[flow diagram](docs/architecture/website-flow.drawio) shows the dataflow and the remaining work;
its [SVG preview](docs/architecture/website-flow.svg) is generated from it by
`scripts/render_flow_diagram.py`.

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

### Custom settings

Wang remains the default. Supported binary scientific rules can become gates, preferences or
be switched off. Required generation inputs and continuous tie-break semantics remain enforced.

```python
from sarna.rules.rule import RuleKind

config = WANG_2015.resolve(
    {"gc-content": {"maximum-gc-fraction": 0.65}},
    role_overrides={
        "tss-window": RuleKind.PREFERENCE,
        "tss-preferred-window": None,   # off: avoid scoring location twice
    },
)
result = design_sarnas(promoter, config)
print(config.customizations)             # changes from the frozen preset baseline
```

The generated page labels custom settings and shows disabled rules, sources, consequences and
alternatives. A matching numerical value is distinguished from a deliberately selected sourced
option. Adjusting a rule changes filtering/ranking; it does not establish improved activity.
See [configuration verification](docs/development/custom-configuration-verification.md).

### From a gene symbol

```python
from sarna.genome.ensembl import promoter_for_symbol

gene, transcript, promoter = promoter_for_symbol("SCN1A")   # MANE Select, or it raises
result = design_sarnas(promoter, WANG_2015.resolve())
```

`select_transcript` returns the MANE Select transcript or **refuses**, listing the alternatives: a
promoter is a property of a transcript, and SCN1A alone has 16 distinct TSSs across 75 kb
([ADR 0012](docs/adr/0012-mane-select-or-refuse.md)).

### Checking a rule against measured activity

```python
from sarna.validation.benchmark import benchmark
from sarna.validation.dataset import load_chapter4_dataset

report = benchmark(load_chapter4_dataset(), WANG_2015.resolve())
print(report.summary)
for lost in report.rejected_activators:        # duplexes that worked and the gate would discard
    print(lost.sequence.name, [r.rule_id for r in lost.candidate.excluded_by])
```

## Current limitations

- **The workbench is an initial local release.** Automated offline tests cover jobs, scientific
  configuration, coordinates, persistence and exports; Chromium checks cover cached SCN1A,
  fresh live-Ensembl CHD2/KCNQ2 designs, custom reruns and explicit cache-miss recovery.
  Windows/macOS installers were not verified in this environment. Packaged installers remain
  to be built; see the
  [continuation plan](docs/development/continuation-plan.md).
- **3' overhangs are not modelled.** Candidates are 19 nt cores; what a researcher orders is a
  21-mer with a 3' dTdT overhang on both strands. Order-ready CSV and chemistry validation
  remain to be implemented.
- **Empirical validation is narrow.** 125 assayed sequences exist, but all were designed to pass
  Wang's gate, so the set cannot show a gate to be too permissive — only too strict. The sequences
  have recovered coordinates and 519 bp context for 123 of 125 designs, enabling window and
  flank checks through `benchmark(..., locations=load_chapter4_locations())`. This context is
  insufficient for the default 500 bp CpG-island definition. Ranking remains unvalidated.
- **No cell-type binding data.** Motif predictions are sequence-only and say nothing about any
  cell type. TF-occupancy adapters remain absent; AF22 histone ChIP-seq is a separate evidence axis.
- ~~Licensing metadata needs correction before distribution.~~ **Resolved 2026-10-07.** The code
  is Apache-2.0 and the bundled JASPAR subset is CC BY 4.0, verified against
  [JASPAR's own about page](https://jaspar.elixir.no/about/). The repository's previous
  CC BY-NC-SA label was wrong and *more* restrictive than the real terms. See [`NOTICE`](NOTICE)
  and [ADR 0018](docs/adr/0018-apache-2-and-bundled-data-licensing.md).
- Promoter fixtures must be CpG-depleted or real sequence. Uniform random sequence is called as
  one large CpG island by any correct detector, which rejects nearly every candidate — see
  [`docs/scientific/synthetic-promoter-fixtures.md`](docs/scientific/synthetic-promoter-fixtures.md).
- The legacy scripts have **18 catalogued defects**, three scientifically material. They are
  preserved unfixed by design; see [`legacy/README.md`](legacy/README.md). The `legacy-scanner`
  preset excludes 5 of the 14 duplexes known to activate their target, four of them on its
  zero-CpG rule.
- The AF22 histone bigWigs (GSE108506) are on **GRCh37/hg19**, while Ensembl REST serves
  **GRCh38** by default. Verified 2026-10-06 from the bigWig headers. The legacy annotation
  pipeline combined the two without converting, so **every regulatory annotation it produced is
  unreliable** (defect B15). The current AF22 reader performs guarded GRCh38-to-hg19 liftover,
  preserving conversion failure reasons; the design path retains assembly-tagged coordinates.
- `Designs.xlsx` is still missing, but the thesis's own chapter 4 supplies a labelled set of 125
  assayed sequences, which is now the project's validation fixture.

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
│   ├── design/              Contexts, evaluators, presets and the design pipeline
│   ├── regulatory/          TF motifs, Ensembl Regulatory Build, the evidence model
│   ├── validation/          Assayed sequences and the benchmark that can refute a rule
│   └── view/                View model and the self-contained HTML browser
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

To design against a real gene and open the result:

```bash
.venv/bin/python scripts/build_promoter_view.py SCN1A --html scn1a.html
```

The first run fetches from Ensembl and caches to `.cache/ensembl`; afterwards it works offline.
The page carries its own data, so it can be mailed or committed beside a manuscript.

### Configuration

Everything is optional and `.env.example` documents all of it — a test asserts that file and the
code name the same set of variables, so neither can drift from the other.

| Variable | Effect |
|---|---|
| `SARNA_CONTACT_EMAIL` | Sent to Ensembl in the `User-Agent`, as their usage guidance asks. Unset sends the tool's name and purpose without an address; the placeholder in `.env.example` counts as unset. |
| `SARNA_CACHE_DIR` | Cache root, default `.cache`. On WSL keep it off `/mnt/c`. |
| `SARNA_OFFLINE` | Refuse every external lookup and serve only cached responses, so a recorded run cannot quietly re-fetch a changed annotation. A cache miss fails naming what was missing. |
| `SARNA_ANCHORS` | Extra promoter anchors — a `.toml` file or a directory of them. |
| `SARNA_CELL_LINES` | Extra cell lines, each declaring evidence along three axes. |
| `SARNA_AF22_DIR` | Where the AF22 bigWigs live, default `~/af22_data`. Absence is normal; chromatin then reports unavailable. |

The two registries load shipped entries first and refuse duplicate ids, so your files add to them
and never silently replace a shipped entry.

## Licence

**Apache License 2.0** — see [`LICENSE`](LICENSE). Copyright 2026 Joep Vondervoort.

Bundled data keep their own terms, listed in [`NOTICE`](NOTICE):

| Bundled asset | Terms |
|---|---|
| JASPAR CORE subset | CC BY 4.0 — redistribution and commercial use permitted with attribution |
| Chapter 4 screening dataset | From a publicly defended thesis; redistributed here under Apache-2.0 |
| Promoter anchors, cell-line registry | Compiled for this project; Apache-2.0 |

Ensembl and NCBI BLAST are queried at runtime and redistribute nothing; cite Ensembl in work
that relies on its annotation. [ADR 0018](docs/adr/0018-apache-2-and-bundled-data-licensing.md)
records why Apache-2.0 rather than MIT, and how the JASPAR mislabelling was found and corrected.
