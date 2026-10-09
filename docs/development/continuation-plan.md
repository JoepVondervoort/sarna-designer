# Continuation plan: 7 October 2026

## Objective and authority

Build a human saRNA design workbench that other researchers can use to choose a promoter,
inspect cell-context evidence, design duplexes and export an auditable experimental record.
Scientific traceability takes priority over a persuasive score or interface.

The user's additional requirement is that researchers can adjust the application to their own
experiments. Defaults are explained starting points. Every supported adjustment must show what
it changes, why the default was chosen, where alternatives come from and what remains uncertain.
The configuration editor and its explanations are core functionality.

The user's supplied continuation discussion establishes 16 decisions. Their latest clarification
keeps Wang as the default and selects local installation; those supersede the earlier proposal
to change the default window role. Other papers inform optional recommendations, not universal
rules. These decisions supersede
the corresponding open questions in `CLAUDE.md`; implementation should not ask those questions
again. This document distinguishes implemented foundations from remaining application work.
Read the existing ADRs, glossary and scientific notes before changing their contracts.

## Verified starting point

- The Python engine, Ensembl adapter, rule catalogue, empirical benchmark and static HTML
  renderer exist. There is no interactive application backend or run store.
- Literature-defined anchors already exist in `sarna.genome.anchors`, including three SCN1A
  regions. The distal reference point is explicitly inferred and uncertain; it must not acquire
  a measured-TSS label simply because another annotation source becomes available.
- Cell lines already have a TOML registry and expression, chromatin and TF-occupancy axes.
  This is availability metadata, not an implemented expression-data adapter.
- The AF22 reader already implements enrichment over same-contig background and guarded
  GRCh38-to-hg19 conversion. Preserve those safeguards while generalising it.
- `target-length` already has a maximum of 19. Check the complete supported-length contract
  before describing the chemistry as strictly fixed at 19 nt; an upper bound alone is not that
  contract. Order-ready dTdT duplex export remains absent.
- `tss-window` remains mandatory in Wang's default, as the latest user instruction requires.
  Custom configurations can investigate a preference or disable it, with an explicit departure.
- An existing uncommitted change in `src/sarna/view/template.html` adds visible rule sources,
  reasoning, thresholds, alternatives and configuration-derived window bands. Preserve it,
  review it and exercise it in a browser before counting the transparency work as complete.
  *Committed since, and loaded in a browser by the user. See the status section below.*
- Local verification at the time this plan was written: 1,279 tests passed, eight live-network
  tests deselected; strict mypy over `src tests scripts` and `ruff check .` passed. This does not
  establish coverage, live-source correctness, browser behavior or availability of the actual
  histone tracks. *Superseded: 1,300 pass, mypy is clean, `ruff check .` now fails with 9 errors
  and coverage is 89%, both because of `src/sarna/app/`. See the status section below.*
- The subsequent custom-configuration slice implements effective roles and frozen baseline/choice
  metadata, with benchmark and browser checks recorded in
  [configuration verification](custom-configuration-verification.md). Full milestone 1 and the
  interactive application remain unfinished.
- The thesis PDF is present at `docs/191591_vandeVondervoort_BNW_def_DPR.pdf`. This review used
  the repository's extracted dataset and scientific notes; it did not independently re-extract
  the full thesis. `pdftotext` is unavailable in the environment.

README and CLAUDE.md now reflect settled dTdT chemistry, recovered coordinates, assembly handling
and the latest default/delivery decisions. Source metadata and scientific notes still need
reconciliation; historical passages sometimes precede later corrections.

The supplied discussion mentions “seven defects” without enumerating them. Do not invent their
identity or equate them with the 18 legacy defects. Establish an explicit current defect list
from code, tests and rendered output during the first milestone.

## Status: 7 October 2026, 13:33

**Read this section's timestamp before trusting its numbers.** The local-app slice was being
written by a second agent *while this review ran*: `src/sarna/app/static/` did not exist at
13:19, appeared at 13:28, and `app.js`, `service.py` and `models.py` changed again at 13:32.
Nine `ruff` errors present at 13:25 were gone by 13:32. Re-run the checks rather than trusting
what follows.

Measured at 13:33: `1300 passed, 8 deselected`; `mypy` strict clean over 108 source files;
`ruff check .` clean; statement coverage **89%**, and **every one of the 399 uncovered statements
is in `src/sarna/app/`**, which is at 0%.

### What is done

- **Milestone 1, transparency half.** The page now renders what the rule records carry. Three
  commits landed the export layer (`Provenance.to_dict`, `RuleParameter.to_dict`,
  `RuleOption.to_dict`, `Rule.to_dict`, `DesignConfig.to_dict` joining by rule id) and the
  template that reads it: per-rule citation, evidence class, rationale, consequence, thresholds
  with units and legal ranges, the sourced alternatives not chosen, and context bands derived
  from the configured window rather than hard-coded to Wang's. A generated page carries 18 rule
  records with `rationale`, `consequence`, structured `provenance` and valued parameters; the
  footer's claim about provenance is now true of the data behind it.
- **Milestone 1, defects half.** A declared `target-length` range that crashed at 20 and 21 is
  capped at 19; chromatin signals carry their own `Provenance` so a HeLa profile can no longer
  cite AF22's accession; the cell-line dispatch asks the reader which line its tracks describe
  instead of testing identity against AF22; `src/sarna/settings.py` owns the environment
  contract and the three documented-but-dead registry variables now work; `--full-rules` reaches
  the HTML as well as the JSON.
- **Milestone 2, engine and export.** `Preset.resolve` takes `role_overrides` and `option_ids`.
  Effective roles reach gating, preference counting, ranking, candidate verdicts and exports
  together; the frozen preset baseline, disabled rules, deliberate option identity and custom
  departures all survive into the record. [ADR 0013](../adr/0013-researchers-can-adjust-rule-roles.md)
  and [ADR 0014](../adr/0014-local-delivery-and-wang-default.md) record the accepted policy, and
  [configuration verification](custom-configuration-verification.md) records the hand-computed
  expectations.
- **Browser check, manual.** `scripts/build_demo_views.py` writes offline default and custom
  previews; `out/scn1a.html` is a real cached SCN1A record. The user served `out/` over
  `http.server` and confirmed all three pages load. No automated browser check exists: Playwright
  has no Chromium installed in this environment, so page behaviour is covered by data-level
  assertions in `tests/unit/test_html_render.py` and by manual inspection only.

### What exists but is unexercised

`src/sarna/app/` is a local workbench with both halves present as of 13:33: a backend and the
static interface it serves.

The backend is real and was verified end-to-end by hand, offline, through `AppRouter.dispatch`
with no socket: `GET /api/genes/SCN1A?offline=true` returns 30 transcripts with
`ENST00000674923` as the unique MANE Select and the three declared SCN1A literature anchors;
a submitted run completes through its progress stages to 2,482 generated and 14 passing; rank 1
(`saRNA-1549`) reports `6/7 evaluated preferences` with `antisense` as the predicted guide under
`rna-xia-1998` and both tables agreeing; `export.csv` writes the 14 ranked rows with both strands
5′→3′; `export.json` and `export.html` both answer. An out-of-cache window fails with an error
that names the exact query, the cache path and the two ways to fix it, which is the behaviour the
offline contract promises.

`RunStore` persists to one SQLite file with progress stages, compresses each snapshot, and marks
runs interrupted by a restart. `DesignRequest` validates input and resolves the configuration
through the same Python path as the library, so an invalid custom configuration is rejected before
a job is queued. The server binds loopback only and checks `Host` and `Origin`.

*Superseded (the paragraph that stood here said the five modules had no tests, that nobody had
opened the interface, and that no Chromium existed to do it with. All three are now false. See
the measured re-check below.*

### Measured re-check) 7 October 2026, 15:40

The claims above were re-measured rather than re-read, and the environment itself had changed.

- **`src/sarna/app/` is tested.** `tests/unit/test_app.py`, `test_app_batches.py` and
  `test_app_references.py` cover it at **82%** (589 statements, 94 uncovered). What remains
  uncovered is the socket layer and the launcher: `web.py` at 75% (its `make_server` request
  handlers, lines 263–302) and `__main__.py` at **0%**, 27 statements. Both are exercised by the
  browser and launcher checks instead, which `pytest --cov` cannot see. Neither is a scientific
  path; `models.py` (90%), `references.py` (100%), `service.py` (89%) and `store.py` (99%) are the
  ones that decide anything.
- **Chromium is installed**, at `~/.cache/ms-playwright/chromium-1229/chrome-linux64/chrome`
  (1200 and 1228 are also present). `scripts/verify_workbench_browser.py` needs no browser Python
  package, only `--chrome`. The plan's premise that this check could only ever be manual was
  wrong.
- **The interface has been opened, twice, and both runs pass.** Cached `--gene SCN1A`: 14 passing,
  draft retained across a Python rerun with its diff, history, hover, offline rejected verdicts and
  the mobile layout. Live `--gene KCNQ2 --online` against an **empty cache**: the offline lookup is
  refused with 409 *without fetching*, the visible retry then retrieves gene, transcripts and
  promoter, a lowercase symbol resolves, an uncached-window retry completes without discarding
  earlier results, and the run yields **9 passing** (the same number
  [local-workbench.md](local-workbench.md) records for KCNQ2. No uncaught browser exceptions in
  either run. **This is the gene-to-candidates path for an arbitrary gene, working end to end.**
- **`ruff check .` had one real failure**) `ANN202` on an untyped `factory` in
  `test_app_batches.py`. Fixed; the gate is clean again. The status section above recorded ruff as
  clean at 13:32, so this arrived with the app commit.

Totals after this slice: **1,383 passed**, 8 deselected; strict mypy clean over 113 files;
`ruff check .` clean.

### The `end-length-nt` wiring: the plan was wrong about what was unpinned

`_compare` reads the selected asymmetry rule's `end-length-nt` rather than hard-coding 4. The
status above said **no test names this pipeline path**. That was not true, and the error was found
the only way it could be (by planting the defect and running the suite:

| Defect planted in `_compare` | Caught without the new tests? |
|---|---|
| `end_length=4` hard-coded again | **Yes**) `test_app.py::test_custom_end_length_controls_displayed_guide_and_all_model_comparisons`, on its hand-calculated +2.49 kcal/mol |
| Tiebreak consulted *before* the gate | **Yes** (same test, incidentally, via the same value |
| No-rule fallback changed from 4 to 6 | **No. All 1,378 tests passed.** |

So the wiring was pinned and the plan's D3 was stale) but the **fallback default was entirely
unpinned**, and the documented default could have been changed from 4 to 6 with a green suite.
[ADR 0015](../adr/0015-local-workbench-and-frozen-ui.md) states that default as part of the
accepted contract, so nothing pinned a decision the ADR records.

`TestTheStrandPredictionHonoursTheConfiguredEndLength` in `tests/unit/test_design_pipeline.py`
now covers all four branches of `_compare`'s selection, at the pipeline level where the existing
test works through the app. Its fixture is the part worth keeping: `GGGTAATTTTGACAGGTCA`, found by
enumerating random 19-mers until one had **all three tables agreeing on antisense at the default
four-base end and all three agreeing on sense at six**. Both verdicts are unanimous, so the test
demonstrates the *consequence* the wiring exists to prevent (the opposite strand reported with
full confidence, which is the wrong duplex to order) rather than only pinning an energy value.
Each of the five tests was confirmed against a planted defect.

### Defect list with dispositions

Milestone 1 requires an explicit current defect list. These are confirmed by running the code,
not inferred.

**Dispositions re-measured 15:40 on 7 October.** Six of the nine are closed; the rows keep their
original wording so what was claimed stays visible beside what was found.

| # | Defect | Evidence | Disposition |
|---|---|---|---|
| D1 | `src/sarna/app/` has no tests, 399 uncovered statements across five modules, plus untested `app.js` | `pytest --cov`: every app module 0% | **Closed for the modules that decide anything.** Now 82%: `references.py` 100%, `store.py` 99%, `models.py` 90%, `service.py` 89%, `web.py` 75%, `__main__.py` 0%. The remainder is the socket layer and the launcher, covered by the browser and launcher checks instead. `app.js` still has no unit tests; it is exercised through the browser check |
| D2 | The interface has never been opened in a browser | No Chromium for Playwright here; no manual check recorded | **Closed.** The premise was wrong: Chromium 1229 is installed. `verify_workbench_browser.py` passes for cached SCN1A (14 passing) and for live uncached KCNQ2 (9 passing), with no uncaught exceptions |
| D3 | `_compare` reading `end-length-nt` from the configuration has no test | `grep end_length tests/` finds thermodynamic coverage only | **Stale as written**, the evidence was a grep that missed `test_app.py`, which does pin the wiring. The real gap was the **no-rule fallback default**, which a green 1,378-test suite tolerated being changed from 4 to 6. Now closed, with all four selection branches covered |
| D4 | The generated viewer and the app's interface are now two separate renderers of the same run | `src/sarna/view/template.html` (1,074 lines) and `src/sarna/app/static/app.js` (38 KB); `export.html` inlines the latter | **Already decided, not open.** [ADR 0015](../adr/0015-local-workbench-and-frozen-ui.md) keeps `sarna.view` for existing scripts and historical reports, explicitly not as a second settings editor, and requires a separate compatibility review before migrating its output. The open question below predates that ADR |
| D5 | The frozen `export.html` is 18.9 MB for SCN1A against 2.1 MB for the same gene from `template.html`; it embeds every rule verdict for all 2,482 candidates with no abbreviation | Measured 7 October | **Accepted cost, not a defect.** Re-measured at **19.25 MB**, of which 19.15 MB is the embedded payload and 18.92 MB is the 2,482 candidates (~8.0 KB each: `rules` 4.2 KB over 18 records, `guide_prediction` 2.2 KB). ADR 0015 accepts this in as many words, "complete rejected-candidate verdicts make HTML large, but preserve every explanation offline". The `_compact` abbreviation in `promoter_view.py` remains available if the cost is later judged too high; applying it would be a **change to the ADR**, so it needs the user, not a fix |
| D6 | `export.html` still references `/favicon.svg` and links to `/`, both of which break once the file is archived away from a running app | Measured 7 October: those are the only two external references left | **Closed.** `web.py:53-55` base64-inlines the favicon and rewrites the brand link to `#` |
| D7 | `app.css` sets `font-family:Inter` with no `@font-face` and no stylesheet link, so Inter never loads and the page silently falls back to the system UI font | `grep -c "@font-face\|fonts.googleapis" app.css` → 0 | **Closed** by taking the second option: a deliberate `ui-sans-serif, system-ui, …` stack. Inter is gone, so nothing promises a font it cannot load, and the frozen export stays dependency-free |
| D8 | `tss-window` offers an option labelled "−2000 to +500 bp" whose effective `window-end-bp` is 0 | Recorded earlier in this plan; still open | **Closed.** The option is now labelled "-2000 to -1 bp (upstream part of thesis chapter 4 window)" and its description states that it implements only the upstream portion and does not reproduce the full thesis window |
| D9 | `GLOSSARY.md` defined the duplex overhang as `UU`, contradicting the settled dTdT chemistry | Glossary versus constraint 4 and `docs/scientific/` | Fixed 7 October alongside this status update |

Remaining open: `app.js` has no unit tests (D1's tail). The user's interface item 4 (rebuilding
the promoter visualisation with hover explanations) turned out to be **already implemented**:
`app.js` draws seven labelled tracks (reference/TSS, design windows, CpG islands, regulatory
build, motif matches, passing and excluded targets), each label carrying a tooltip that says what
the element is *and* what it does not establish, with candidates clickable and zoom/pan working.
A screenshot was inspected. If it still is not what was wanted, that needs saying concretely;
nothing further is obviously missing from the description in the feedback.

### Order-ready duplex export: 7 October 2026, 16:10

The gap `CLAUDE.md` called "largest before bench use" is closed, except for the screening half.

`sarna.design.order` builds the two 21-mers and their sanity report;
[ADR 0017](../adr/0017-order-ready-duplex-export.md) records the decisions and their costs. The
workbench serves them at `GET /api/runs/{id}/export.order.csv` and the browser has an
**Order CSV · 21-mers** button; a frozen archive downloads the same rows, which Python stored in
the record so that JavaScript never re-derives chemistry (accepted constraint 11). A test asserts
the stored rows and the live route's output are byte-identical.

What the slice found while building it:

- **Both CSV exports were in position order, not rank order.** Measured on cached SCN1A,
  `export.csv` emitted ranks `8, 3, 7, 5, 13, 4, 2, 1, …`. A researcher taking the top three rows
  ordered the wrong three molecules. This was pre-existing, not introduced here. Both exports now
  sort; the record keeps position order, which is correct for the map and the sortable table.
- **`off-target-screening` reports `not-performed` on every duplex.** Accepted constraint 15 wants
  off-target findings in this report and no off-target analysis exists, so the row states the
  absence rather than being omitted, which would read as a clean screen.
- The five checks were each confirmed against a planted defect: a ribo-`UU` overhang (caught by
  three tests, after one self-referential assertion was rewritten to pin the literal), an overhang
  on one strand only, and `not-performed` silently downgraded to `ok`.

Two live `--online` browser runs failed during this slice with a 45 s timeout. The cause was
established rather than assumed: a direct repro showed Ensembl returning `ReadTimeout` after six
attempts on `lookup/id`, which is gene lookup, upstream of anything this slice touches. A
clean-tree run passing in between was coincidence, not evidence. The check passes on both trees
when Ensembl responds.

### Open question this status raises: answered by ADR 0015

*The question below was written before [ADR 0015](../adr/0015-local-workbench-and-frozen-ui.md)
was read. That ADR answers it: `sarna.view` stays for existing scripts and historical reports and
is explicitly not a second settings editor, and migrating its output requires a separate
compatibility review. The duplication is accepted and bounded rather than open. What the ADR does
not do is remove the underlying risk the question names (two renderers can still disagree about
one run) so if that is to be guarded, it wants a test comparing the two payloads for the same
run, not a new decision.*

There are now two independent renderers for one run: the self-contained
`src/sarna/view/template.html`, which carries the committed transparency work and has data-level
tests, and `src/sarna/app/static/app.js`, which the workbench serves and which `export.html`
inlines to produce a frozen page. Two renderers means every transparency guarantee has to be
built and kept true twice, and the two can disagree about the same run. Converging on one is a
decision with a real cost either way and belongs in an ADR.

### The user's interface feedback, 7 October

Recorded verbatim in substance, because it sets milestone 6's acceptance criteria and contradicts
nothing already accepted:

1. A researcher chooses **any supported human gene**; the tool looks up its promoter and designs
   duplexes against it. Pre-assigned demonstrations are previews, not the product.
2. The **predicted guide strand** must be visible per candidate, with the model that predicted it.
3. The **rating** must be legible and must say what it counts.
4. The **promoter visualisation needs rebuilding**, with hover explanations saying what each
   element means rather than assuming the reader knows.
5. **Every setting is adjustable from the browser**, not by editing Python.

Items 1–3 and 5 are served by the existing backend and blocked only by D1. Item 4 is new design
work on the promoter track and is the one part with no implementation behind it yet.

## Accepted constraints

1. Wang 2015 remains the default, including its window gate. Other settings and roles are
   optional research choices; literature conformance does not establish activation efficacy.
2. Distance is a sortable column with no points. A custom window preference must explain the
   overlap with the existing preferred-window criterion; do not silently give distance extra points.
3. `pos19-a` remains mandatory in the accepted default; methods must name the known CEBPA
   exception. Researchers may change its role in an advanced custom configuration, as accepted
   in [ADR 0013](../adr/0013-researchers-can-adjust-rule-roles.md).
4. Order chemistry is a 19 nt RNA core plus a 3′ dTdT overhang on both strands.
5. Cell-context assessment precedes candidates, is promoter-level and never changes ranking.
6. The assessment has three outcomes and names its sources. Missing data is not a negative assay.
7. The 5× enrichment threshold is an explicit, provenanced parameter and never blocks design.
8. Human-only support is a real, visible limitation.
9. The application runs locally, with accessible installation and launch paths on supported
   Windows, macOS and Linux systems. The library remains the core; exports are frozen records.
10. Transparency includes a truthful interface, per-run methods and live reranking with a diff.
11. Python is the sole design engine; JavaScript displays results rather than implementing rules.
12. Use the accepted MIT/Apache licensing direction. Audit the actual bundled JASPAR assets:
    official JASPAR documentation states CC BY 4.0, contradicting the repository's non-commercial
    label. Unbundling is therefore not established as a licensing necessity. Record exact asset
    terms and attribution, and the exact code licensing form before replacing LICENSE.
13. Adapter implementation order is archived EPDnew, DepMap, FANTOM5, ENCODE.
14. Saved runs live in one local SQLite file.
15. Order export includes CSV and a sanity report, including off-target findings. Findings report
    rather than exclude; build promoter-level off-target analysis first.
16. A short guided setup leads into a workbench with deliberate anchor and cell-line selection.

## Scientific interpretation to preserve

The AF22 measurements support low enrichment of the assayed active histone marks at the SCN1A
regions in that dataset. They do not directly measure baseline transcription, prove complete
silence or establish why every null result occurred. Keep a chromatin warning separate from an
expression observation and from a proposed biological explanation.

Gene-level expression cannot establish which alternative promoter is active. Display its scope.
HeLa expression remains unverified here; the DepMap adapter must first establish that the exact
line and relevant expression dataset are covered. A missing row is unknown, not zero expression.
Neither additional context nor the accepted rule changes establishes predictive accuracy: the
chapter 4 dataset is range-restricted and ranking has not demonstrated performance above chance.

## Researcher-adjustable settings and their explanations

### Controls to deliver

| Area | Researcher controls | Explanation beside the control |
|---|---|---|
| Starting configuration | Wang default, sourced alternatives and saved custom configurations | Which source supports each part, departures and why they were chosen |
| Sequence criteria | Supported GC bounds, CpG limits, repeat thresholds and related numeric/choice parameters | Definition, units, default origin, alternative evidence and direction of the filtering change |
| Location | Promoter span, anchor/reference, supported TSS window and distance sorting | Coordinate convention, total containment, anchor uncertainty and whether the control filters or scores |
| Rule participation | Set supported scientific rules to mandatory, preference or off in advanced custom configurations | Off means not evaluated; a preference can fail without rejection; a failed gate excludes |
| Chromatin assessment | Supported enrichment threshold, query window, marks and evidence source | Background definition, assay limitations and effects on the warning without effects on ranking |
| Evidence and context | Cell line, available samples/tracks and local evidence | Exact sample identity, assay scope, assembly and missing evidence |
| Results and ordering | Sort columns, candidate selection and saved-run comparison | Viewing/selection changes versus changes requiring a Python design run |

Keep the accepted equal-weight preference count. Adjustable scientific criteria do not imply
adding arbitrary score weights. Keep 19 nt cores with dTdT for the supported order workflow;
different lengths or chemistries require explicit engine/evaluator support, rather than an input
that merely accepts another number. Explain any supported limitation at its control.

Recommended first-release scope is the existing catalogue and local custom configurations;
researcher-authored executable rules and a general plugin system can be considered later.
Describe this as the proposed scope, not a user decision already made.

### What each setting must explain

Each control must expose the current value and role, the selected preset's value and role, and
the source rule's published setting where it exists. These are distinct baselines. Provide:

1. A short plain-language definition, units and the operation it affects.
2. Why this default was selected, with citation and precise table/section/page locator where
   available. Cite project policy as project policy and thesis findings as thesis findings.
3. Named alternatives with their own source and reason for choosing them. A source for the rule
   is not automatically a source for every value a user types.
4. The expected computational consequence of increasing/decreasing the value or changing the
   role, and the biological motivation or uncertainty. Separate these claims.
5. Measured consequences on this promoter after a preview/rerun: candidates admitted or excluded,
   preference denominator changes, rank movement and changed rejection reasons.
6. Known exceptions and limits of validation. A benchmark observation must state its dataset;
   an unmeasured custom setting must not acquire an efficacy claim.

The interface can show a brief explanation beside the setting and expand citations, methods
and caveats. The explanation must be available before applying the change. Unsupported evidence
for a custom value is an explicit answer, not a reason to forbid a computationally valid value.

For example, raising the maximum GC fraction from 0.60 to 0.65 relaxes that gate and allows
candidates above 60% and up to 65% GC to pass it. The additional survivors must be measured
after all other gates: some candidates may still fail elsewhere. The app must not imply that
those extra candidates have demonstrated activity, or attribute the custom value to Wang.

### Editing, comparison and reuse

- Start with an explained preset and offer an advanced editor using the same rule metadata.
  Keep the published Wang settings as the default baseline and clearly name custom departures.
- Label changed values and rules. Provide per-setting reset and whole-preset reset. Show which
  settings are unsupported or inapplicable in the current context and why.
- Keep draft settings separate from the displayed completed run. Apply/preview executes Python,
  then shows the resulting changes. If a calculation is pending or fails, identify which
  configuration the visible results still belong to.
- Show excluded candidates and their reasons; do not add up per-rule rejection counts as if they
  were disjoint. Compare settings with the same sequence, anchor and evidence snapshot where
  possible, and name any additional changes that prevent attributing the difference to one rule.
- Save a custom configuration with a name and version locally; support export/import of the
  configuration for collaborators. An optional research note can explain the user's rationale.
  Changing a saved configuration creates a new version and cannot alter a completed run.
- Record custom departures in run summaries, candidate explanations and methods. Changing the
  preference set changes the denominator; different scores are not directly comparable.

## Configuration foundation required by adjustability

Extend the Python contracts before exposing broad controls. The first implementation slice adds
effective roles, selected option identity, a frozen preset baseline, explicit disabled rules and
departures. It preserves catalogue records and shares validation across library entry points.
Configuration schema/versioning, import, user notes and the interactive editing path remain work
for the persistence and app milestones.

The planned configuration record contains:

- Versioned base preset, base values and base roles.
- Effective selected rules, supported roles and resolved values.
- Explicit choice origin: preset, named sourced option or researcher custom value. Record the
  option id when selected; a numerical match alone does not establish the user's source.
- Disabled rules and changed settings, so omissions are visible in methods as deliberate choices.
- Required explanatory metadata, units, supported ranges and source locators.
- Configuration schema version and optional user note; complete resolved data is sufficient to
  reopen a run without applying a changed preset definition.

Use one Python validation path for the library, import/export and eventual UI. Validate finite
numeric values, integer-only fields and relationships between settings, such as minimum GC not
exceeding maximum GC and window start not exceeding window end. Reject unknown and ineffective
overrides instead of silently ignoring them. Revalidate imported custom configurations.

Distinguish an invalid computation from an exploratory scientific choice. Bounds required for
algorithm correctness should be labelled as supported limits, and source recommendations as
recommendations. Every restriction needs a reason. Do not widen a range until the evaluator,
coordinate handling and methods description support it.

An observed item for milestone 1: `tss-window` offers a thesis option labelled “−2000 to +500 bp”
whose effective `window-end-bp` is 0. Its note admits the restriction, but the option label and
actual value still differ. Display the implemented range accurately and explain the source
range and restriction; allowing downstream targeting is a separate support decision.

## Milestone 1: reconcile the current implementation

Review and preserve the existing template changes. Generate representative pages for both
strands, inferred and absent TSS references, and different presets/windows. Exercise candidate
selection, sorting, rule explanations and missing-data states in a browser. Inventory concrete
defects and fix only confirmed ones, with appropriate regression checks.

Update README and CLAUDE.md to distinguish current functionality, accepted choices and remaining
work. Record accepted scientific/architecture changes in ADRs, including measured costs and
superseded statements. Keep `legacy/` immutable.

**Exit:** page claims agree with visible behavior; the defect list has named dispositions;
documentation no longer contradicts settled chemistry, anchor support or assembly handling.

## Milestone 2: implement the accepted rule policy

Support resolved rule roles without mutating global catalogue entries. Preserve Wang's default
rules and parameters. Label a changed run as custom and show its frozen baseline. A changed role
does not establish a better design recommendation. Preserve total-containment
semantics. Distance sorting must state what
position it uses and carry inferred-reference caveats.

Keep `pos19-a` mandatory in the default and pin 19 nt order chemistry. Rerun the chapter 4 benchmark using
coordinates, report exclusions by rule, candidate counts and preference denominators, and record
the costs. “Removing the window's rejections” must not be reported as “every activator passes
all gates”: CEBPA still fails `pos19-a`.

**Exit:** gates, scoring and exported configuration reflect the accepted policy; published and
historical configurations remain reproducible; independent coordinate expectations cover both
strands and window edges.

**Adjustability deliverable:** a validated custom configuration can select supported rules and
values and change scientific rule roles, preserve the chosen source and preset baseline, and
export every departure. Python validates this contract before any UI exposes it. A role change
must reach evaluation, gating, scoring, candidate explanations and exports consistently.

## Milestone 3: generalise cell-context evidence

Extend the existing cell-line registry with explicit tracks: mark, local path or URL, assembly,
assay/sample identity, source accession and processing details. Separate track identity from a
reader tied to AF22. GSM2902406 is AF22 H3K4me3 and links both GSE108506 and GSE111327;
the accession itself is not shown to be wrong. Audit each track and correct the publication
attribution against its source record. Preserve liftover rejection reasons and local track use.

Make enrichment threshold, background sampling, seed, interval and statistic part of the run
record. Produce a three-valued promoter assessment with source, reason and missing-data details.
Do not reuse TF “bound/not bound” vocabulary for expression or histone measurements. Conflicting
evidence remains visible rather than disappearing into a summary label.

**Exit:** unknown, measured-negative and measured-positive evidence are distinguishable; warnings
appear before candidates, remain non-blocking and cannot affect rank.

## Milestone 4: add provenance-aware sources

Build a source registry with release/version, retrieval date, assembly, identifiers, licence,
cache location and content checksum. Add EPD first, then DepMap, FANTOM5 and ENCODE. Verify each
source's current format, coverage and terms from its primary documentation during implementation.
Recorded fixtures must include absent, ambiguous and incompatible records; add targeted live
checks for external contracts.

EPD was retired; implement the final November 2023 EPDnew archive rather than depending on a
continuously updated API. Record archive DOI, version and checksums. EPDnew incorporates FANTOM5
and ENCODE data, so overlapping records are not independent corroboration.

EPD offers measured anchor alternatives without silently replacing an inferred reference.
DepMap adds gene expression in explicitly matched cell models. FANTOM5 adds supported TSS evidence;
ENCODE adds assay-specific cell-context evidence. Support researcher-supplied local evidence for
lines absent from public datasets.

Resolve JASPAR's incorrect licence metadata before distribution; retain attribution and choose
bundling versus cached retrieval based on verified terms and installation practicality.
Explain external retrieval in setup, including what gene/coordinate queries are sent. Keep
researchers' uploaded sequences and results local. “Local app” does not imply optional evidence
retrieval is network-free. Provide an explicit cached/offline mode.

Audit the September 2026 curated saRNA database and its linked primary experiments as optional
external validation material. Check sequence identity/orientation, endpoints, cell models, negative
results, duplicates and overlapping studies before drawing conclusions. Hold out genes/studies
where feasible. Database recommendations do not replace the Wang default or become truth by
being newer; rank performance and empirical exceptions remain claims requiring measurement.

**Exit:** every observation can be traced to an exact source record, and cached/offline runs do
not silently refresh their evidence.

## Milestone 5: freeze runs before building the workbench

Define a versioned run record containing input sequence, assembly and coordinates, anchor and
uncertainty, cell context, resolved rules and kinds, software version, evidence snapshots and
source hashes, candidates, verdicts and ranks. Store snapshots in SQLite with schema migrations;
cached URLs alone are insufficient archival evidence.

Changed settings produce a derived run linked to its parent. Compare candidates by stable
sequence/location/orientation identity rather than rank or a display ID alone. Report admitted
and excluded candidates, rank movement, preference changes and changed reasons.

**Exit:** a saved run reopens offline without consulting today's preset definitions, and a diff
explains a deliberate parameter change.

**Adjustability deliverable:** saved custom configurations are versioned independently of run
snapshots. Import/export preserves choice origins and supported limitations; reopening a frozen
record requires no access to the current preset catalogue for interpretation. Recomputing it
also requires the recorded engine version or an explicit report of version differences.

## Milestone 6: deliver the local workbench

Deliver a browser workbench backed by a thin local Python service. Bind to the loopback interface;
launch it through an application launcher. Keep science in Python and use HTML/CSS/JavaScript for
controls, display and progress. No hosted service or account is required. Store runs and caches
in a documented per-user data directory, separate from the installation, so updates preserve them.

Guided setup selects gene or pasted sequence, deliberate anchor, cell line and configuration.
The workbench shows promoter context, evidence assessment, passing and rejected candidates,
rule explanations, settings, reranking and saved-run comparisons. A pasted sequence can proceed
with explicit unknown genomic context. Exports remain frozen even when the workbench can rerun.

Use the researcher-settings contract above as the editor's acceptance criteria. Complete one
working path from changing a GC bound, through Python validation and rerun, to an explained diff
and saved configuration before adding all controls. Then add rule participation/roles within
the agreed scope, location controls, chromatin controls and the remaining catalogue fields.

**Exit:** a researcher can complete setup, inspect warnings, change a setting, compare and reopen
a run through the application without writing Python.

### Acceptance criteria from the user's 7 October feedback

The application written on 7 October claims to serve most of this and has no tests or browser
check behind the claim. Each criterion is checked by using the running application, not by reading
its code.

- Typing a human gene symbol returns its transcripts and declared literature anchors, and the
  researcher picks one deliberately. No arbitrary transcript is chosen for them.
- Each candidate shows its predicted guide strand, named with the model that predicted it, and
  says so when the two nearest-neighbour tables disagree. The genomic strand of the locus is
  displayed separately and never conflated with it.
- The rating states what it counts (preferences met out of preferences *evaluated*), names the
  gates that could not be checked, and never reads as a probability of activation.
- The promoter visualisation explains itself on hover: what a band, a track, a mark and a
  candidate are, in a sentence a researcher who has not read this repository can act on.
- Every setting the catalogue declares is adjustable in the browser, with its definition, units,
  legal range, default origin, sourced alternatives and the measured effect of the change on this
  promoter. Applying a change runs Python and shows the resulting diff.
- A pasted sequence with no genomic context still works, with the unavailable checks named.

Also verify a sourced alternative, a custom value, an invalid combination, a reset, an imported
configuration and missing-context behavior. The researcher can distinguish the setting's
motivation, its measured effect on this run and the absence of validated efficacy evidence.

## Milestone 7: produce auditable bench exports

Build promoter-level off-target reporting first, with explicit reference/search scope and
matching assumptions. Do not call it comprehensive genome-wide or transcriptome-wide safety.
Missing search results are labelled unavailable; off-target findings never exclude candidates.

Generate order CSV with both strands written 5′→3′, unambiguous RNA cores and `[dT][dT]`
overhang notation, core/total lengths, candidate identity and run identity. The sanity report
checks reverse complementarity of the cores, chemistry, lengths, duplicate orders, orientation
and available off-target evidence. Do not silently interpret mixed DNA/RNA ordering notation.

Generate per-run methods from the saved record: anchor derivation, cell-context evidence,
citations, resolved thresholds, gates/preferences, ranking/tie-breaks, CEBPA exception, empirical
validation limitations and source versions. Bundle frozen HTML, methods, CSV and sanity report.
Methods and chemistry serialization can be developed against milestone 5 before the UI is done.

**Exit:** an exported selection agrees with the saved run and can be reviewed independently of
the application. Claims about off-target coverage and activity remain bounded by actual evidence.

## Immediate next implementation slice

**Superseded on 7 October by what was actually built.** Milestone 1's transparency work and
milestone 2's engine, export and ADRs are done (see the status section above), and a milestone 6
backend was then written out of order with no interface attached. The next slice is therefore to
finish that application rather than to start milestone 3:

1. Test the router and the store offline through `AppRouter.dispatch`, which takes a method and a
   target and needs no socket. Cover an invalid custom configuration, a gene that does not exist,
   a pasted sequence with no TSS, a run that fails mid-way, a reopened snapshot and an export.
2. Test the `app.js` behaviour that carries a scientific claim, the rating's denominator, the
   predicted guide and its model, rejected candidates staying visible with their reasons, and
   missing evidence reading as unknown rather than negative.
3. Add the pipeline test for the configured asymmetry end length (D3).
4. Use the running application against a real cached gene and read what it shows. Record what is
   wrong with it as defects before changing anything.
5. Decide the one-renderer-or-two question (D4) and record it as an ADR before the two diverge
   further.

Only then return to milestone 3. The earlier sequencing argument still holds for everything after
this slice: resolve accepted policy before building adapters around it.

The historical text of this slice follows, since its milestone-2 steps describe work already
landed and its reasoning is still the reasoning.

Deliver milestone 2 in small steps:

1. Record the default-versus-custom policy and published baseline in ADRs.
2. Extend configuration resolution, choice origins and shared validation; exercise a sourced
   option and a hand-edited GC bound through the Python engine.
3. Preserve the Wang default and allow explicitly labelled custom roles and sourced alternatives.
4. Export the resolved values, effective roles, disabled rules and deviations; measure candidate
   changes and benchmark exclusions, then render example pages that explain them.

Run persistence and the first settings-to-diff workbench path depend on these contracts and the
source/evidence snapshot format. They do not need to wait for every public adapter to be complete.
Keep the EPD → DepMap → FANTOM5 → ENCODE adapter order while delivering the application
incrementally with explicit unavailable evidence.

The user has accepted advanced changes of supported scientific rules between mandatory,
preference and off. [ADR 0013](../adr/0013-researchers-can-adjust-rule-roles.md) records this choice,
its costs and the explanation/export requirements. Named preset defaults remain versioned
baselines. Engine and frozen-page support are implemented; the interactive editor, configuration
versioning and import/export remain pending.

## Local installation and release requirements

Start a packaging feasibility check before substantial app work: native genomics dependencies,
including bigWig access, can constrain supported operating systems and architectures. Test what
can actually be bundled; do not promise universal compatibility from Python alone.

The research-user path should install and launch the workbench without editing Python or requiring
Docker. Provide platform packages where feasible, a documented Python installation for developers
and an optional container route for labs. If a platform initially needs an alternative path,
state the support boundary and concrete steps rather than claiming a native installer exists.

Each released platform needs a clean-machine installation and upgrade check: launch, an offline
example design, an invalid custom setting, a saved run, export, restart and reopen. Verify that an
upgrade preserves the local database and provide backup/restore instructions. Include a concise
quick start, optional-source downloads and sizes, supported versions/architectures, troubleshooting,
and an example dataset. Large optional resources must not be mandatory just to start the tool.

Release packaging follows the first complete settings-to-saved-run app path. The feasibility
check starts earlier so distribution constraints can shape dependency choices in time.

## Primary sources for the corrected source strategy

- [EPD retirement notice and final archive](https://www.expasy.org/archives/epd)
- [EPDnew human data sources](https://epd.expasy.org/epd/human/human_database.php)
- [JASPAR documentation and licence](https://jaspar.elixir.no/about/)
- [AF22 H3K4me3 sample GSM2902406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2902406)
- [2026 curated saRNA database paper](https://doi.org/10.3390/molecules31183249)
- [Authors' database and source records](https://github.com/NIC-SBI/sarna-database)

## Workbench verification update: 7 October 2026, after implementation

This supersedes the 13:33 review's untested-app observations while preserving them as history.
The local workbench is implemented and exercised. Launch and installation instructions are in
[local-workbench.md](local-workbench.md); architectural choices are in
[ADR 0015](../adr/0015-local-workbench-and-frozen-ui.md).

- **D1–D3:** 22 new offline tests exercise jobs, multiple recorded genes, minus-strand
  coordinates, custom roles/values/options, missing evidence, guide-comparison end length,
  snapshots and exports. Full suite: 1,322 passed, eight live tests deselected. Ruff and strict
  mypy pass; pytest's combined statement/branch coverage is 97%. The independently derived
  GC/AA two-base RNA difference pins the adjusted pipeline path at +2.49 kcal/mol.
- **D2, browser evidence:** `scripts/verify_workbench_browser.py` passed with cached SCN1A and
  pasted DNA. Actual Chromium checks cover draft preservation, reruns/diff, saved reopening,
  hover, mobile layout and an offline HTML file's full excluded-candidate verdicts. Desktop
  and mobile screenshots were inspected. No Playwright Python dependency was needed.
- **D4–D5:** The live workbench and its frozen export share the same UI. The older CLI viewer
  remains for compatibility. Complete archival verdicts intentionally cost file size, while
  SQLite compresses snapshots. See ADR 0015 for scope and later renderer migration.
- **D6–D8:** Frozen favicon is inline and the root link is neutralized; fonts use an explicit
  system stack. Upstream rule bounds now end at −1, and the thesis alternative clearly names
  its upstream-only extent rather than claiming +500 or attempting a nonexistent zero.
- **Milestone 6:** Gene/transcript/literature-anchor or pasted-DNA setup, all supported controls,
  a new promoter map, named guide predictions and evaluated preference ratings, local saves,
  previous-run comparison, and frozen HTML/JSON/core CSV are implemented. Cached SCN1A gave
  14 passing; changing the window gate to a preference gave 57, with 43 newly admitted and
  13 shared rank changes. This measures changed criteria, not improved activation.

The app remains an initial local release. Live-source access and clean Windows/macOS installs
were not verified here; no expression adapter, generalized cell-specific histone reader,
off-target report, dTdT order export or methods document is claimed. Source version pinning,
licence/metadata reconciliation and these remaining milestones still matter.

## New-gene lookup correction: 7 October 2026

The user's failing `chd2` and `kcnq2` requests contained `offline=true`: Ensembl was never
contacted. Lowercase input also bypassed existing uppercase cache entries, because normalisation
previously ran only after a remote "gene not found" response, not after an offline cache miss.

- Normalise human symbols before cache/network access and in saved requests.
- Place the cache-only control and effective mode in the target setup. Reopening a saved run
  restores its mode visibly. `SARNA_OFFLINE` remains an enforced server setting and the UI
  explains how to restart without it.
- Distinguish uncached/offline records, unknown genes and service outages. An explicit online
  retry can recover either lookup or a missing promoter window during design. Previous results
  remain saved. Offline requests never silently fetch.
- Version 2 of the local SQLite schema adds a failure code; existing snapshots survive the
  migration. No scientific defaults, scores or strand calculations changed in this correction.
- Fresh-cache Chromium runs fetched live CHD2 and KCNQ2 promoters and generated 2,482 targets
  each. Wang defaults yielded 13 and 9 passing candidates respectively. Lookup recovery,
  lowercase input, reruns, saved reopening and frozen reports passed. KCNQ2 additionally passed
  recovery after changing to an uncached 2,501 bp promoter window.
  A separate fresh-cache CHD2 job also completed with default motif scanning and regulatory
  retrieval enabled, without a regulatory-unavailable notice.
- The offline suite now passes 1,327 tests with eight original network tests deselected;
  strict mypy and Ruff pass. Live browser checks are separate from that offline test count.

This supersedes the earlier statement that new live-source access had not been verified.
See [the troubleshooting instructions](local-workbench.md#if-another-gene-fails-to-load).
New genes/windows still require internet unless the necessary reference responses are cached.

## Multi-promoter workbench extension: 7 October 2026

The user's follow-up makes deliberate design for genes with several promoters part of the
current implementation scope. The prior single-reference workbench remains the scientific
foundation; [ADR 0016](../adr/0016-multiple-promoter-designs.md) records the new decision and
revises ADR 0012's earlier rejection of designing against several TSSs without pooling ranks.

- A whole-gene reference overview groups annotated transcripts only when their contig, strand
  and TSS coordinate match exactly. Member transcript IDs remain visible. This is coordinate
  grouping, not an inference of biological promoter equivalence or cell-specific activity.
- Supported human gene symbols and stable Ensembl gene IDs lead to the same reference-selection
  workflow. An optional gene-ID version suffix is removed for current-record lookup; it does
  not pin a release. Broad synonym and ambiguous-alias resolution is not claimed by this slice.
- Researchers can select several annotated/literature references and enter their own genomic
  reference with a name, source/derivation and uncertainty note. Custom GRCh38 coordinates are
  entered 1-based, converted to 0-based once, and labelled researcher-supplied/inferred. Gene
  contig and strand are retained. Position and window must fit the primary chromosome; an
  outside-gene reference is allowed with its gene association explicitly unverified.
  Annotation and citations do not prove measured activity. Completed gene runs freeze the
  available reference overview returned at lookup.
- Each selected reference starts from the shared rule configuration and can retain pinned
  criteria and its own retrieval span. Defaults, supported roles/values, equal preference
  counting, named guide models and 19-nt default cores remain unchanged.
- Python creates separate saved child runs. A saved comparison group links them and reports
  per-reference progress/failures. Successful results survive another reference's failure.
  Reopening uses saved snapshots; group JSON and offline HTML preserve their distinct contexts.
- The comparison retains each reference and configuration. It neither merges ranks nor declares
  a best promoter. A duplex shared by overlapping windows may have different relative positions
  and verdicts under each reference. Existing same-reference rerun comparisons remain separate
  from cross-reference inspection.

The guide's [multi-promoter instructions](local-workbench.md#design-against-several-promoters)
and editable [flow diagram](../architecture/website-flow.drawio) describe this workflow.
Verification of this extension is pending. The 1,327-test count and live CHD2/KCNQ2 observations
above belong to the preceding version; record the new suite/browser results here only after
they complete. Expression/occupancy adapters, measured-source coverage, arbitrary-cell tracks,
order-ready exports, off-target reporting and clean installers remain continuation work.
