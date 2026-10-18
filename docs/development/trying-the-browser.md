# Open and test the current browser

The **generated pages** below are viewers over a finished run: settings are changed in
Python and the page regenerated. The **local app** runs new designs from a gene symbol,
supports browser-editable settings and saves runs. Automated tests cover the engine, application
routes and exports; separate Chromium checks exercise browser interactions. For installation
and the current workflow, start with the [workbench guide](local-workbench.md).
Fresh Windows/macOS installations still need platform validation.

## Start with the offline examples

In a terminal at the repository root, with the existing `.venv`:

```bash
.venv/bin/python scripts/build_demo_views.py
.venv/bin/python -m http.server 8000 --bind 127.0.0.1 --directory out
```

Open **http://localhost:8000/** in your browser. Leave the terminal running; press Ctrl+C to
stop serving it. If port 8000 is occupied, use 8001 in both the command and the URL. This server
serves files only; it does not run scientific calculations or save runs.

You can also open `out/index.html` directly in a browser without starting a server. On WSL,
open the output directory in Windows Explorer with:

```bash
explorer.exe "$(wslpath -w "$PWD/out")"
```

Then double-click `index.html`. The generated records carry their own data. Optional Google
Fonts may be fetched for styling; fallback fonts work when offline.

The two demo pages use the same synthetic DNA, without genomic coordinates or a TSS. They are
interface examples, not measured gene designs. The default currently gives 10 passing candidates
out of 582. The custom example gives 37, with a 65% upper GC bound, position-19 as a preference
and the TSS gate explicitly off. Three settings change together; counts alone cannot identify
each setting's effect or demonstrate better efficacy.

## Generate a real gene view

The current workspace also contains `out/scn1a.html`, generated from cached Ensembl records:
GRCh38, minus strand, MANE Select transcript ENST00000674923, 2,500 bp promoter, Wang defaults.
It has 14 passing candidates out of 2,482. AF22 histone tracks were deliberately skipped. That
absence is not a negative assay. This preview uses the MANE anchor; the literature anchors
were not chosen.

To regenerate this view, or replace SCN1A with another supported human gene:

```bash
.venv/bin/python scripts/build_promoter_view.py SCN1A \
  --no-chromatin --html out/scn1a.html --out out/scn1a.json
```

The first run needs Ensembl access and the genomics extra described in README. Cached requests
can be enforced without network access by prefixing the command with `SARNA_OFFLINE=1`;
a missing cached response then fails explicitly. Histone files are not required for this example.
Rerun `build_demo_views.py` after generating SCN1A if you want its link added to the preview index.

To compare existing CLI presets, write a second page rather than overwrite the first:

```bash
.venv/bin/python scripts/build_promoter_view.py SCN1A --preset exploratory \
  --no-chromatin --html out/scn1a-exploratory.html --out out/scn1a-exploratory.json
```

Open that file from the browser's address bar or directly. These are full independent design
records; there is no saved-run comparison screen yet. Missing MANE Select causes the CLI to list
alternatives and stop instead of silently choosing an arbitrary transcript.

## What to try manually

1. Click a candidate in the table or promoter track. Check its target, strand sequences,
   measurements, rule outcomes and reasons.
2. Switch **Passing / All**. Rejected candidates must stay visible and name their rejection rules.
   Large tables display the first 400 entries; use the track to find other positions.
3. Drag the track, zoom and return to the whole promoter. Inspect TSS alternatives and toggle
   motif predictions on the real-gene page.
4. Read **Rules and their sources**. Open thresholds and alternatives. On the custom demo,
   verify the changed roles, baseline departures and disabled rule record.
5. Check missing cell-context evidence. Motif predictions do not establish binding in a cell
   type; skipped histone tracks do not establish closed chromatin.

To try another custom value or role, edit the example configuration in
`scripts/build_demo_views.py`, rerun it and refresh the page. Invalid values should fail before
new records are written; previous output remains on disk. Use separate `--out-dir` directories
to preserve earlier examples. Browser controls currently change the view, not the rule settings.

For automated engine checks:

```bash
.venv/bin/python -m pytest -q tests/unit/test_custom_configuration.py
.venv/bin/python -m pytest -q
```

## The local application

`src/sarna/app/` is a local workbench: a loopback-only Python service and the page it serves.

```bash
.venv/bin/python -m sarna.app --port 8001
# then open http://localhost:8001/
curl -s localhost:8001/api/catalogue | head -c 300   # rules, presets and cell lines
```

Automated tests exercise application routes, saved runs, configuration changes and exports.
The Chromium verification script also checks the live browser workflow, including screening
displays and mobile layout. Saved runs can be reopened for inspection. Software verification
does not establish that a candidate will activate its target; see the
[scientific evidence guide](../scientific/evidence-guide.md).

The first run for a gene needs Ensembl access; afterwards `.cache/ensembl` serves it. A request
for a promoter window that is not in the cache fails with the exact query it would have sent.

Most route behavior can be exercised without a socket; selected integration checks also use
real localhost HTTP. For a direct route check:

```python
from pathlib import Path
from sarna.app.service import DesignService
from sarna.app.web import AppRouter

router = AppRouter(DesignService(Path("/tmp/sarna-app-check")))
status, media, body = router.dispatch("GET", "/api/catalogue")
router.close()
```

Runs are stored in `~/.sarna-designer/runs.sqlite3` unless `--data-dir` says otherwise.

## Flow diagrams

Open [website-flow.drawio](../architecture/website-flow.drawio) in draw.io. It contains two
pages: **Current multi-promoter workbench** and **Remaining continuation**.

Page one is the dataflow as built: input, reference choice, one Python run per reference, saved
local comparison, the opt-in off-target screen, and the order-ready duplex that the exports are
written from. Remote BLAST sends candidate sequences to NCBI; local reference screening
keeps candidate queries on the computer. Gene and annotation retrieval contact their source
services, and reference preparation downloads external data. On page two, green marks completed
work. See the current workbench and screening guides for implemented behavior.

[website-flow.svg](../architecture/website-flow.svg) is a browser-viewable preview, and it is
**generated**: edit the `.drawio` and re-run:

```bash
.venv/bin/python scripts/render_flow_diagram.py
```

A test fails if the two disagree, because when both were maintained by hand they drifted.
Neither diagram claims that planned work is already implemented.
