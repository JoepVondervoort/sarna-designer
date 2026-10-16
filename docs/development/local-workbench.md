# Run the adjustable local workbench

From the repository root, with the existing virtual environment:

```bash
.venv/bin/python -m sarna.app --port 8001
```

Open **http://localhost:8001/**. Leave the terminal running; Ctrl+C stops it. The earlier
`http.server --directory out` on port 8000 serves precomputed previews. This command serves
the adjustable app and runs scientific calculations in Python. If 8001 is occupied, choose
another `--port` and the matching URL.

## Try a real design

1. Enter **SCN1A**, another supported human gene symbol, or a supported Ensembl gene ID and
   click **Look up**. Inspect the **Gene overview** and choose a transcript/start-site group or
   a declared literature reference. Transcripts that start at exactly the same coordinate are
   grouped, with their identities retained. MANE Select is selected only when unique. Neither
   a grouped annotated start site nor MANE establishes cell-specific promoter activity;
   literature anchors retain their source and uncertainty.
2. Choose the upstream/downstream span and experimental cell line. GTEx expression can provide tissue context; it does not establish transcription in the
   selected experimental cell line. AF22 histone tracks are an optional
   local input; other lines do not acquire histone evidence simply by being selected.
3. Keep Wang defaults or open **All rules & evidence**. Edit supported values, choose a sourced
   alternative, or make a binary rule a gate, preference or off. Citations, evidence classes,
   ranges, rationale and consequences sit beside the controls. Required generation inputs and
   continuous tiebreaks retain their semantics.
4. Click **Design saRNAs**. Hover or focus map elements for explanations. Zoom and pan; click
   a candidate to inspect its DNA, both RNA cores and every rule verdict. The table is paginated,
   including all excluded targets.
5. Read the **named model's predicted RNA guide strand**, separately from genomic DNA strand.
   Disagreement between models is retained. **Rating counts preferences met / evaluated**, not
   activation probability. Scores across different configurations are not interchangeable.
6. Change a setting and **Apply changes & design**. Saved results stay visible while edits are
   a draft. A fresh Python run reports admission and rank changes on the same promoter/reference.
   Both snapshots remain in **Saved runs**.
7. Export **Order CSV · 21-mers**, **Core sequences CSV**, **JSON**, or **Offline HTML**.
   - **Order CSV** is what you synthesise: both strands as 19 nt RNA cores plus a 3′ `[dT][dT]`
     overhang, 21 nt each, in 5′→3′ orientation. Core, overhang and joined form are separate
     columns, because vendor notation for a chimeric oligo differs between suppliers. Every row
     carries its own sanity report, including the latest completed screening findings per candidate and scope, or an explicit
     not-screened status. A similarity search does not establish biological specificity. Passing
     candidates only, unless you select specific ones. See
     [ADR 0017](../adr/0017-order-ready-duplex-export.md).
   - **Core sequences CSV** is the engine's view: RNA cores alone, no overhangs. Deliberately not
     order-ready.
   - Both CSVs are written in **rank order**, so the first row is rank 1.
   - **HTML** preserves complete candidate explanations and works without the server, including
     the order download, whose rows Python stored in the record. Full records can make HTML files
     large.

## Design against several promoters

1. Look up the gene and inspect its **Gene overview**. Choose a **Promoter reference** and click
   **Add to comparison**, then repeat for each reference you want to investigate. A whole-gene
   overview shows the available references; it does not guarantee that all biological promoters
   have been measured or annotated. Exact shared-start grouping is not activity evidence.
2. Inspect **Selected references**. Set **Upstream (bp)** and **Downstream (bp)** for each
   reference. References using **Use current criteria** share the current editor settings.
   **Capture current criteria** pins the current criteria to one reference; **Use shared
   criteria** returns it to the shared settings. **Captured settings** also shows the saved
   cell and optional annotation controls. The visible current cache-only mode applies at
   submission. Review these choices before designing so a later shared edit does not
   accidentally change the interpretation of a comparison.
3. Click **Design selected promoters**. Python runs the engine separately for each reference.
   **Promoter comparison** shows each run's state and candidate counts. A failed retrieval is
   reported for that reference; completed results remain saved. Cache-only failures can be
   retried explicitly through **Allow Ensembl access & retry failed** when online access is
   permitted. Different windows need their own cached sequence responses.
4. Click **Open run** to inspect a reference's map, candidate ranks, guide models, rule verdicts
   and exact configuration. Rankings are separate within each promoter/configuration. The
   comparison does not choose the most active promoter, combine ranks, or make a cross-promoter
   efficacy prediction. The same genomic duplex can appear under several references with
   different offsets, verdicts and ranks.
5. Reopen the group through **Saved runs**, or preserve it with **Group JSON** and
   **Group offline HTML**. Group exports retain the member contexts and complete explanations;
   they do not recompute designs. Individual child runs remain available for their existing
   exports and reruns.

For SCN1A, include the relevant literature reference rather than assuming its default MANE TSS
covers the earlier experiments. The Long-associated reference and the Dong-defined references
are separate choices; the inferred Dong distal reference retains its caveat.

## Add your own genomic reference

After looking up the gene, open **Your genomic reference** and enter **Reference name**,
**Genomic +1 coordinate (1-based)**, **Source or derivation**, and **Uncertainty / context note**.
Click **Add researcher reference**, then select it for an individual design or comparison.

The genomic coordinate is the GRCh38 base you are treating as the first transcribed base; the
app takes contig and strand from the confirmed gene. Python converts this input to its internal
0-based genomic reference once. The input is not the −1/+1 relative coordinate used for pasted
DNA. The saved reference is labelled researcher-supplied and inferred, with your source and
note. The position and fetched window must fit the primary chromosome; the position may lie
outside the annotated gene span. Entering a publication does not independently validate the
coordinate, a mapped TSS, or its association with this gene. Check the coordinate, assembly and
experimental interpretation before using it.

Human Ensembl gene IDs may include their version suffix, but lookup removes that suffix and
uses the current GRCh38 record. A versioned ID does not freeze an annotation release. Broad
synonym and ambiguous-alias resolution remain future work; use the current gene symbol or
Ensembl gene ID when a name is not recognized. Completed gene runs also retain the available
reference overview returned by their lookup.

New genes/windows require Ensembl access. **Use cached Ensembl data only** prevents requests
and explains cache misses. This control and the effective connection mode are visible next to
the target setup, including after reopening a saved run. Human symbols are case-insensitive:
`chd2` and `CHD2` use the same lookup and cache entry. New installations have no cached SCN1A
data. Alternatively choose
**Your DNA**, with plain DNA or one FASTA record in transcript orientation. Supply the first
base's TSS-relative position if known; there is no position zero. Upstream-window rules currently
support negative positions only. The thesis-derived alternative explicitly implements −2000…−1,
its upstream portion, rather than the full −2000…+500 promoter definition.

## If another gene fails to load

If the error says **Cache-only mode**, the app did not contact Ensembl. Turn off **Use cached
Ensembl data only**, or click **Allow Ensembl access & retry**. The app retrieves missing
gene/transcript/sequence records and caches them; later requests for that same window can run
offline. A cached gene does not mean every adjustable promoter window is cached. The retry also
handles a missing window during design and preserves previously completed results.

If the mode says it is enforced by **`SARNA_OFFLINE`**, that terminal setting takes precedence
over the checkbox. Stop the app with Ctrl+C, then restart from the repository root:

```bash
unset SARNA_OFFLINE
.venv/bin/python -m sarna.app --port 8001
```

After a code update, restart the app and reload the browser so both use the new version. An
unknown gene and an unavailable Ensembl service have different error messages. Keep online
retrieval allowed for a gene/window's first run. An internet connection is needed for uncached
references; the design engine and saved results run locally.

## Install from a checkout

Use **Python 3.11**, as required by `pyproject.toml`. Linux/macOS/WSL:

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -e ".[app]"
.venv/bin/python -m sarna.app
```

Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[app]"
.\.venv\Scripts\python.exe -m sarna.app
```

These are source-install instructions, not validated installers. The `app` extra adds HTTPX for
Ensembl. No FastAPI, Node build or bigWig package is required for the standard workbench. Optional
AF22 tracks additionally need the existing `genomics` extra and actual local track files. WSL was
exercised; fresh dependency downloads and clean Windows/macOS installs were unavailable to verify.

Runs live in **`~/.sarna-designer/runs.sqlite3`** by default. To isolate test runs:

```bash
.venv/bin/python -m sarna.app --data-dir /tmp/sarna-test-runs
```

Back up that directory with the app stopped. Ensembl caches use `SARNA_CACHE_DIR` (default
`.cache` in the launch directory). The server binds to loopback for one local researcher.
External reference requests go to Ensembl; pasted DNA and saved snapshots are processed locally.

## Verification

```bash
.venv/bin/python -m pytest -q tests/unit/test_app.py
.venv/bin/python -m pytest -q
.venv/bin/python -m mypy src tests scripts
.venv/bin/python -m ruff check .
```

The repeatable browser check uses a supplied Chromium executable, an isolated temporary
database, and no browser Python package:

```bash
.venv/bin/python scripts/verify_workbench_browser.py --chrome /path/to/chromium
```

It defaults to pasted DNA and needs no cached gene. Add `--gene SCN1A` for a cached gene run,
and `--screenshots out/browser-check` to retain viewport captures. It checks draft preservation,
Python reruns/diffs, saved-run reopening, tooltips, mobile layout, and an exported HTML opened
as a local file, including rejected-candidate verdicts.

To verify a new gene with a fresh temporary cache, explicitly enable network access:

```bash
.venv/bin/python scripts/verify_workbench_browser.py --chrome /path/to/chromium --gene kcnq2 --online
```

This first checks that cache-only lookup refuses the uncached gene, then clicks the online retry
and completes a live gene-to-promoter-to-candidates workflow. It also checks retry after changing
to an uncached promoter window. `--online` requires `SARNA_OFFLINE` unset/off. The default browser
checks switch off optional motif/regulatory annotations to isolate lookup and design behaviour.

On 7 October 2026, Chromium exercised the real cached SCN1A MANE reference, 2,000 bp upstream
and 500 downstream: **2,482 targets, 14 passing, top rating 6/7**. Making `tss-window` a preference
gave **57 passing**, 43 newly admitted and 13 shared candidates with changed rank. Draft edits
left the original count visible until rerun. Gene reference selection, pagination, selection,
zoom/pan and the rule editor were exercised, with no uncaught browser exceptions.

The full offline suite passed **1,327 tests** (eight live-network tests deselected); strict mypy
and Ruff passed. Combined statement/branch coverage was 97% in the pytest run, which does not
include the standalone browser checks. HTTP serving and command startup were additionally
exercised in the browser/launcher checks, rather than covered by that pytest measurement.

Fresh-cache browser checks on 7 October 2026 also retrieved **CHD2** and **KCNQ2** from live
Ensembl: each generated 2,482 targets from a 2,500 bp MANE-referenced promoter. Wang defaults gave
**13 passing** for CHD2 and **9 passing** for KCNQ2. Lowercase input, explicit online recovery,
saved reopening, custom reruns and offline HTML inspection passed. The original eight network
tests remain separate from these live browser checks.
An additional fresh-cache CHD2 job completed with both default regulatory retrieval and motif
scanning enabled; regulatory retrieval produced no unavailable-source notice.

These are observations on specific source responses, not efficacy estimates or counts guaranteed
for future Ensembl annotations. Automated backend tests use recorded responses for multiple genes
and hand-derived strand/coordinate and two-base thermodynamic expectations. External
expression/occupancy adapters, order chemistry and off-target screening remain outside this
verification.

The multi-reference extension is being verified separately; the test counts above describe
the preceding single-reference version and must not be read as its verification result.

See [ADR 0015](../adr/0015-local-workbench-and-frozen-ui.md) for the local architecture and
[ADR 0016](../adr/0016-multiple-promoter-designs.md) for reference grouping, custom genomic
inputs, separate ranks and remaining limits.

## Research workspace

Use the stage bar to move from target/reference selection to candidates, screening and export.
**Design setup** stays visible after a run; **Go to design setup** returns to its controls. The candidate inspector sits beside the
table on wide screens. **Add to research shortlist** saves a selection with your reason, separately
from scientific results; table checkboxes remain temporary export/screening selections.

**Run history** supports search, status filters and comparison of two saved runs. Read the
[research workflow guide](research-workflow.md) for shortlist persistence, conflict recovery,
offline notes, exports and the limits of cross-run comparison.
