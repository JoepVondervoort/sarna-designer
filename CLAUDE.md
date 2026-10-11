# saRNA Designer: working notes for Claude

Python library that turns a promoter into ranked, explained saRNA candidates. Scientific
traceability outranks clean code: every rule carries a citation, every verdict carries its reason.

**State:** gene symbol → promoter → ranked candidates → interactive page all work.
`design_sarnas(promoter, config)` is the engine; `EnsemblClient` resolves a gene to a promoter;
`sarna.regulatory` annotates it; `sarna.view.html.render` writes a self-contained browser.
1,931 tests passed in 49 s, live-network tests excluded; mypy strict and Ruff clean.
**Verified, not assumed:** the whole suite passes under `unshare -rn`, with no network
namespace at all, so it genuinely reaches no external service.
**The unit suite reaches no external service.** It used to: the screening route built a real
`BlastClient`, so every full run queued a live BLAST job at NCBI and then blocked closing the
router until it came back. `DesignService.blast_factory` and `poll_sleep` are the seams that
stop it, mirroring `client_factory`; a test asserts the route uses them.
The earlier combined statement/branch coverage measurement was **97%**; socket/launcher paths
are additionally
exercised by separate browser checks, which are not part of that coverage measurement.
`src/sarna/app/` sits at 82%, the remainder being `web.py`'s request handlers and `__main__.py`.

**`src/sarna/app/` is now an exercised local workbench.** Gene lookup, adjustable supported
rules, Python jobs, draft/rerun diff, SQLite snapshots and JSON/HTML/core-CSV exports work.
The candidate view opens on a **duplex panel**: the two strands drawn antiparallel with a rung
only where the bases pair, both 3' dTdT overhangs, guide and passenger labelled, and the two
end energies with the stacks they were summed from, the measurement the guide call *is*, which
until now existed only inside one English sentence. The promoter map pans by dragging and zooms
on the wheel; excluded targets are binned into a composition histogram above 280 bp because at
full extent they were three overlapping marks per pixel. Read
`docs/scientific/verification-2026-10-07.md` for what has been checked against its source.
`tests/unit/test_app.py` covers recorded sources, coordinate/strand invariants, configuration,
storage and exports. `scripts/verify_workbench_browser.py` uses a supplied Chromium binary:
pasted DNA and cached SCN1A checks passed, including saved reopening, tooltips, mobile layout
and complete excluded-candidate verdicts in an offline HTML file. Desktop/mobile screenshots
were inspected. Fresh-cache live Ensembl browser runs also passed for lowercase `chd2` and
`kcnq2`, with 13 and 9 candidates passing Wang defaults. All three paths, cached SCN1A, pasted DNA and
`--gene kcnq2 --online` from an empty cache, were re-run on 8 October and pass; Chromium 1229
is installed at `~/.cache/ms-playwright/chromium-1229/chrome-linux64/chrome`, so this check is
**not** manual, pass it to `--chrome`. The browser check now also covers the duplex ladder
(a rung exactly where the bases pair, both directions), the promoter map's drag and wheel zoom
(comparing the translated picture against the full redraw), the binned excluded-target track
(its bins must account for every rejection), and the cell-context panel. Cache-only mode is visible; missing
records offer an explicit online retry, while environment-forced offline remains enforced.
Fresh Windows/macOS installations remain unverified.
Read `docs/development/local-workbench.md` and ADR 0015. Earlier mid-session review observations
in the continuation plan are historical; its final verification update supersedes them.

## Accepted continuation decisions: 2026-10-07

Read `docs/development/continuation-plan.md` for the user's accepted continuation decisions;
they supersede corresponding open questions in these older working notes. The latest decision
includes [ADR 0013](docs/adr/0013-researchers-can-adjust-rule-roles.md): researchers may adjust supported
scientific rules between mandatory, preference and off in advanced custom configurations.
`pos19-a` stays mandatory in the default, while custom configurations may change its role.
Explain defaults, alternatives and consequences, and record all custom departures. This is an
accepted requirement; engine support and the interactive editor are implemented.
[ADR 0014](docs/adr/0014-local-delivery-and-wang-default.md) records the latest clarification:
local installation is the delivery target, Wang remains the default including its window gate,
and other papers inform optional recommendations rather than universal rules. This supersedes
the earlier proposal to change the default window to a preference.

## Commands

Use `.venv/bin/python`; there is a miniconda base environment that must stay untouched.

```bash
.venv/bin/python -m pytest -q                      # 1,931 tests; excludes live network
unshare -rn .venv/bin/python -m pytest -q          # same, with no network at all
.venv/bin/python -m pytest -m network              # 8 live-Ensembl tests, excluded by default
.venv/bin/python -m pytest -m tracks               # needs the AF22 bigWigs; skips if absent
.venv/bin/python -m mypy src tests scripts         # strict
.venv/bin/python -m ruff check .
```

**No em dashes anywhere outside `legacy/`.** The user asked for the GitHub page, the
documentation and the application's own text to carry no em dashes and no obvious AI phrasing.
The whole tree was swept on 8 October 2026: prose em dashes became commas, colons or sentence
breaks, paired ones became parentheses, and the escaped placeholders in tables became `n/a` or
a word. `legacy/` keeps its own, being immutable.

`tests/unit/test_writing_style.py` enforces this, so it is a failing test rather than a
convention: the character, its JavaScript escape and its HTML entities are all banned, along
with a short list of phrases that have no legitimate use here. Text that has to *name* one of
them marks that line with `writing-style: allow`, per line and conspicuously, so an exemption
cannot be granted quietly. An en dash stays allowed, because this project writes real ranges.
Generated output under `out/` is not scanned.

`ruff format` is deliberately not part of the workflow, the source is hand-wrapped near 88
columns while `line-length` is 100, so running it reformats most files. `ruff check` is the gate.

Scripts that measure rather than assert, which is the house pattern:

```bash
scripts/validate_against_chapter4.py       # every preset against 125 assayed sequences
scripts/positional_motif_significance.py   # ~2 min; permutation test over the motif grid
scripts/locate_scn1a_designs.py            # needs network; maps the thesis's designs to GRCh38
scripts/af22_chromatin_at_targets.py       # needs network + AF22 tracks; controls included
scripts/locate_dataset_designs.py          # needs network; already run, its CSV is committed
scripts/build_promoter_view.py SCN1A --html out.html    # needs network the first time
scripts/build_demo_views.py                # offline synthetic default/custom preview pages
```

The local app; inspect explanations and data gaps, not just the rank:

```bash
.venv/bin/python -m sarna.app --port 8001            # then open http://localhost:8001/
curl -s localhost:8001/api/catalogue | head -c 200   # rules, presets and cell lines
```

`AppRouter.dispatch(method, target, body)` exercises every route without a socket, which is how
the app's tests should be written.

## Invariants

Violating one of these produces wrong science that type-checks and passes unrelated tests.

- **Positions inside a sequence are 1-based; genomic intervals are 0-based half-open.** Table 1
  speaks of "nucleotide 18", so `base_at` matches the paper. The two conventions never meet:
  one indexes within a sequence, the other locates a sequence on a contig.
- **TSS positions have no zero.** `+1` is the first transcribed base, `-1` the base before it.
  All arithmetic goes through `TssPosition`, which converts to a continuous offset first.
- **Promoters are held 5'→3' in the transcript's own orientation.** `Promoter.
  from_reference_sequence` flips a minus-strand reference once, at the boundary. After that the
  engine is strand-agnostic. **All strand arithmetic lives in three methods on `Promoter`**,
  `genomic_interval`, `offset_of`, `span_of`, where offset 0 is the *highest* coordinate on the
  minus strand, plus `EnsemblClient.promoter_for`, which mirrors them. Defect B1's home.
- **Target sites are DNA; duplex strands are RNA.** CpG content is a property of the DNA locus
  and `cpg_count` refuses to answer for RNA.
- **Three outcomes, not two.** `not-applicable` means the rule could not be evaluated, and it
  neither excludes nor scores. A rule the preset omits produces no result at all. "Omitted",
  "could not be checked" and "failed" are three different statements about a candidate.
- **A motif prediction may not name a cell type.** `Evidence` raises if one does. Sequence
  matching is cell-type agnostic by construction, so attaching a cell type implies a measurement
  nobody made. `CellTypeStatus.NOT_ASSAYED` is distinct from `ASSAYED_NOT_BOUND`, and rendering
  them alike turns ignorance into a finding.
- **`DesignResult.generated` carries no ranks.** Ranking returns *copies*, so `ranked` has them
  and `generated` does not. Export `generated` directly and a renderer silently falls back to
  position order; use `PromoterView.ranked_candidates`.
- **`legacy/` is immutable.** Preserved verbatim as the scientific ground truth; linting and
  formatting exclude it.
- **Chromatin "presence" is enrichment, never raw signal.** The AF22 tracks have no input
  normalisation, so every mark reads about 1.0 everywhere; presence is fold over a 250-window
  same-contig background. A raw threshold would call every mark present at every locus.
- **Synthetic promoter fixtures must be CpG-depleted or real sequence.** Uniform random sequence
  has CpG obs/exp of 1.0 by construction and is called as one giant island, rejecting nearly
  every candidate. The symptom looks exactly like a broken design rule. See
  `docs/scientific/synthetic-promoter-fixtures.md`.

## Where decisions already live

Read the ADR before reopening its question; each records the alternatives and why they lost.
`docs/adr/`:

| Question | ADR |
|---|---|
| Which script is ground truth | 0001 |
| Why two presets exist | 0002 |
| Why every rule carries provenance | 0003 |
| Assembly handling, GRCh37 bigWigs vs GRCh38 | 0004 |
| Python library + thin UI, no framework in the engine | 0005 |
| Ranking: count criteria, never a weighted sum | 0006 |
| Evaluator contract, and why `not-applicable` exists | 0007 |
| TSS numbering, and total containment in a window | 0008 |
| Guide prediction reports every table, never averages | 0009 |
| What a preset is; legacy's declared divergences | 0010 |
| Optional rules are hypotheses; the default carries none | 0011 |
| MANE Select or refuse to pick a transcript | 0012 |
| Researcher-adjustable scientific rule roles and transparent custom configurations | 0013 |
| Local delivery, tested installation paths and Wang as the default | 0014 |
| Local Python jobs, and one UI for live and frozen results | 0015 |
| Several promoter references per gene, each its own design context | 0016 |
| The order-ready duplex: separate export, Python-computed, caveats in the row | 0017 |
| Apache-2.0, and the bundled data's actual terms | 0018 |
| Off-target screening: opt-in, remote, dual-scope, reported never filtered | 0019 |
| Localhost rather than a hosted service, and what would change that | 0020 |
| Transcript expression is evidence beside the reference, never a reordering | 0021 |
| Rank is shaded in one hue, never green to red | 0022 |
| Wang's own 80% threshold, screening in batches, and the on-target hit | 0023 |
| The human assembly for genomic screening, and collecting a slow search | 0024 |

`docs/scientific/`:

| Finding | File |
|---|---|
| 18 defects (B1–B18) in the original scripts, cite by id | `legacy-inventory.md` |
| What the measured data says about each gate | `empirical-validation.md` |
| Why chapter 4's positional hotspots are not usable | `positional-motif-analysis.md` |
| Where the thesis's SCN1A designs actually sit | `scn1a-design-coordinates.md` |
| What TF annotation can and cannot honestly show | `tf-binding.md` |
| Why random sequence is an invalid promoter fixture | `synthetic-promoter-fixtures.md` |
| AF22 chromatin: targets open, SCN1A closed | `af22-chromatin.md` |
| Table 1, the thermodynamics and the coordinates, checked against their sources | `verification-2026-10-07.md` |
| Guide prediction: 13% of candidates flip between tables | `guide-prediction-divergence.md` |
| SCN1A's MANE Select transcript is 0.00 TPM in all 54 GTEx tissues | `transcript-expression.md` |

`GLOSSARY.md` fixes the vocabulary, use its terms rather than coining synonyms.

## The empirical dataset

`sarna.validation.dataset.load_chapter4_dataset()` holds 125 assayed target sequences from thesis
chapter 4, 14 of which activated their gene. It is the only thing here that can show a rule to be
**wrong** rather than merely unconventional, and it has already refuted one preset.

Every sequence was built to pass Wang's gate, so the set is **range-restricted**: it cannot show a
gate too permissive, only one too strict. `CHAPTER4_CAVEATS` carries the rest.

`load_chapter4_locations()` adds recovered GRCh38 coordinates for 123 of the 125, with 519 bp of
context each, so the TSS windows and 3' flank preferences are evaluable, pass `locations=` to
`benchmark`. Island overlap still is not: 519 bp cannot support a 500 bp island detector.

Headline results:

- GC, CpG-count and homopolymer gates at Wang's values each exclude **none** of the 14 activators.
  The position-19 gate does exclude CEBPA. Do not describe the full sequence gate set as lossless.
- Wang's **window** gate excludes **10 of 14**, and the 4 survivors are exactly the 4 published
  literature duplexes. Taking the most favourable TSS anchor for each design it still excludes 8.
  Kept in the default by the user's latest decision; custom roles can investigate the cost.
- No preset's *ranking* beats chance (AUC 0.41–0.53), with or without coordinates.
- `exploratory` scores highest because two of its rules were fitted to this very data. The tests
  assert that, so the circularity stays visible.

## How to work here

- **Ask before defaulting.** Scientific choices, a threshold, a definition, which error to
  prefer, are the user's to make. Present the recommendation with its cost, then wait.
- **Record every accepted choice as an ADR**, including the cost of the option chosen.
- One question at a time.
- Tests pin scientific defaults on purpose. A change to a published threshold should break a test
  that names the paper; if it does not, the test is missing.
- Verify with an independent derivation, not by re-running the implementation: hand-computed
  expected values, brute force alongside an incremental algorithm, or planting a known defect to
  confirm the test catches it.
- Measure claims rather than asserting them, then write the number down. The permutation test over
  chapter 4's motif grid and the exhaustive guide-prediction enumeration are the pattern.
- **Run the thing and read its output before calling it done.** The rank bug above survived a
  green suite and 100% coverage; it was found by looking at a generated page.
- Live-API tests earn their keep. Two Ensembl bugs, `mane=1` being required, and a 400 rather
  than a 404 for an unknown symbol, were invisible to recorded payloads.

## Implementation gaps and remaining decisions

- **Measured alternatives to existing inferred anchors.** Literature-defined anchors already
  exist in `sarna.genome.anchors`. The source adapter remains to be built:
  chapter 2 designed SCN1A against promoters mapped by Dong 2014 and Long 2008, and the Dong distal
  region (which produced 2 of the 3 SCN1A activators) has no Ensembl TSS within 5.5 kb. A tool
  that can only design relative to annotated TSSs could not reproduce this project's own hits.
  Preserve explicit uncertainty when adding measured alternatives. See
  `docs/scientific/scn1a-design-coordinates.md`.
- **Local application and installers.** `sarna.view.html.render` produces a self-contained page;
  the local Python service, editor, run snapshots and frozen workbench export now exist.
  Tested cross-platform installation paths and packaged installers remain to be built.
- **Position 18/19 when target length ≠ 19.** Evaluators use the literal Table 1 position and
  report `not-applicable` if the target is shorter. For a 21-mer, position 19 stops being the 3'
  terminal base and `pos19-a`'s mechanism no longer holds. Maximum is already 19; the library
  still supports 18 with explicit not-applicable results. Order export is fixed at 19 + dTdT.
- **3' overhangs are modelled, and the order export exists.** `sarna.design.order` builds the two
  21-mers, 19 nt RNA core plus `[dT][dT]` on both strands, per chapter 2's supplementary material,
  and the workbench serves them at `GET /api/runs/{id}/export.order.csv`, in rank order, each row
  carrying its own sanity report. [ADR 0017](docs/adr/0017-order-ready-duplex-export.md).
  The ordered strand is deliberately **not** a `NucleicSequence`: a chimeric RNA/dT 21-mer fits
  neither alphabet, so the core stays RNA and the overhang is a separate deoxy value.
  What is still missing is the *screening*: every duplex reports `off-target-screening`
  as `not-performed`, because no off-target analysis exists. Accepted constraint 15 wants those
  findings here, so the order export is honest but not yet complete.
- **Window alternatives and their costs.** As a gate at
  Wang's −1000…−100 it rejects 64% of candidates on a real promoter and 10 of the 14 known
  activators, including KLF4 at −45 and CEBPA (in human trials) at −1088. As a preference it would
  admit all 14 and lose the strongest filter the engine has. The rule now offers the thesis's own
  -2000…+500 window as a sourced upstream-only option, actually -2000…-1, with the restriction
  stated in its label and note. Zero cannot be a TSS coordinate.
  Wang's gate remains the default; custom configurations can change its role or disable it.
- **`pos19-a` excludes CEBPA**, the duplex in MTL-CEBPA and the first saRNA in human trials.
  The gate stays conformant with Table 1 in the default. Custom roles are now supported and
  retain their departure from the frozen baseline in exports.
- **Off-target screening exists, and the order export cites it.**
  `sarna.offtarget` screens a finished run against RefSeq RNA (could the duplex silence an
  mRNA?) and the genome (could it act on another promoter?), opt-in and asynchronous, reported
  as several sortable measures and never as a filter. Genomic hits report distance to the
  nearest annotated TSS with **no threshold**; the one flag says a hit falls inside the window
  that run itself used. Screenings are separate dated records, schema v4.
  [ADR 0019](docs/adr/0019-off-target-screening.md). `export.order.csv` reads the newest
  complete screening per scope and reports its findings with the database and date; no
  screening still reads `not-performed`, and a failed one contributes nothing.
- **The screening panel was rebuilt on 11 October; read
  [ADR 0023](docs/adr/0023-off-target-reporting-thresholds-and-batches.md) before changing it.**
  The defect that prompted it: `parse_hits` held **one** budget of 500 alignments for a whole
  multi-query report, so the first query consumed it and every later candidate came back with
  zero. A real 22 candidate run showed 1 candidate with 500 and 21 with none, printed like
  genuine clean results. The budget is per candidate now, and a candidate that reaches it is
  recorded as capped: no alignments, some alignments and more than were kept are three
  different statements.
  Going back to Wang 2015 for the seed question settled three things. The paper **names no
  seed region** (one use of the word, in quotation marks), so none is drawn. It states exactly
  one homology threshold, **80% for rejecting a control duplex**, so `wang_homology_hits`
  (16 of 19) is reported beside the exact count and leads the interface; the locally invented
  `NEAR_PERFECT_IDENTITIES = 17` stays for the exports that already carry it. Its own
  candidate search prescribes **megablast**, which returns nothing at 19 nt, so `MEGABLAST=off`
  is a declared divergence rather than a code comment. **Wang step 9, a BLAT check rejecting
  any target containing a SNP, is not implemented anywhere**, which is a real Table 1 adjacent
  gap.
  Genomic screening used to count **each candidate's own promoter** as a perfect off-target
  match, making that column meaningless. The own-locus hit is identified by exact coordinate
  overlap, moved to `on_target_hits` rather than deleted, and counted; a run from pasted
  sequence has no coordinates and says so instead of guessing.
  Screening now runs in batches (1/5/10/15/25/all, default 15, plus the next batch, the rest
  and the ticked selection), all resolving to one explicit id list the route validates. Several
  dated records per scope are expected and the panel merges them.
  The panel itself shows **where** a duplex binds: the gene, the matched bases, and a
  base-resolution footprint whose mismatch positions come from BLAST's `Hsp_midline`, plus a
  signed distance-to-TSS plot for genomic hits with the run's own window as the only band.
  `ScreeningRecord.window` exists because the in-window flag is meaningless without it.
  The panel lives in `src/sarna/app/static/offtarget.js`, a second classic script sharing
  app.js's globals. **It must be inlined in `frozen_page` and listed in the static route**, and
  `showResult` calls into it through a `typeof` guard because a frozen archive runs app.js to
  completion before this file is parsed.
- **The promoter scope was searching the wrong database; read
  [ADR 0024](docs/adr/0024-human-assembly-and-collectable-searches.md).** It sent
  `DATABASE=refseq_genomic`, a name that appears **zero times** in NCBI's BLAST database
  directory and **zero times** in the database list on NCBI's own blastn page. The URL API
  accepted it and returned a request id anyway, so the search ran, against something that is
  not 3.3 Gbp of human genome, and a real run hit the full half-hour budget and returned
  nothing. It is now `GPIPE/9606/current/ref_top_level`, the human GRCh38 reference assembly,
  which is what NCBI's own page submits.
  **Database size was not the reason, and measuring refuted the guess.** Reading NCBI's own
  RTOE at submission on 11 October 2026: the human assembly is estimated at **28,502 s, about
  8 hours**, and so is every other genome database (`human_genome`, `refseq_genomes`,
  `RefSeqGene`, `refseq_select_rna`, all 27,757 to 28,745 s). Word size 11 and `EXPECT=10`
  change nothing (28,747 s and 28,765 s). Meanwhile `refseq_rna` and `core_nt` are both
  estimated at **30 s**, and `core_nt` is 332 times the size of the human genome. NCBI queues
  the genome databases and serves the others interactively; that is the whole explanation, so
  do not try to tune this with parameters. Lowering `EXPECT` would also lose findings: a 16 of
  19 match on a 3 Gbp database has an E-value in the hundreds, so Wang's 80% threshold is only
  visible *because* `EXPECT` is 1000. `core_nt` restricted to human genomic molecules does run
  in 27 to 40 s, but its hits are RefSeqGene records and clones with no chromosome
  coordinates, and restricting to chromosome-sized RefSeq records returns nothing, because
  `core_nt` does not hold the assembly. **There is no fast remote route to a human
  chromosome-level search.**
  So a search whose RTOE exceeds the wait budget gets **one** poll at the rate limit and is
  then handed back with its request id, and the interface says up front that genome searches
  are an errand to collect rather than something to wait for.
  **The one route to an interactive promoter search is a local BLAST+ backend, and ADR 0019
  rejected it on a wrong number:** it says "roughly 10 GB of database", and the human genome
  BLAST database is 1,157,233,426 bytes, about 1.16 GB. `BlastClient` is already the seam.
  Undecided; worth putting to the user.
  The assembly holds one species, so the genomic scope sends **no** `ENTREZ_QUERY`. That
  matters beyond speed: `ENTREZ_QUERY` is not a documented URL API parameter and NCBI's page
  says parameters it does not list "should not be relied upon". The transcript scope still
  depends on it, so the two scopes are human in different ways and
  `SearchParameters.species_restriction` records which; the panel words them differently on
  purpose. Moving the transcript scope to `GPIPE/9606/current/rna` is an open question, not
  decided.
  **Running out of patience is now a fourth status, `waiting`, not an error.** NCBI holds a
  request id for about a day, so `SearchStillAtNcbi` carries it out, `store.pause_screening`
  keeps it, and `collect_screening` finishes the original job rather than queueing a duplicate.
  A `waiting` screening survives an application restart; a `running` one does not, because
  `running` is this process's work and `waiting` is NCBI's.
  Two quieter fixes in the same area: the per-candidate cap is a bounded best-of **by bit
  score**, because BLAST groups a report by database sequence and on chromosomes the first 200
  alignments are all chromosome 1's weak ones; and placement is bounded at the 300 strongest,
  measured at 1,022 Ensembl requests for 3,000 scattered hits.
  **The progress bar may only claim a fraction it has.** NCBI reports started or finished and
  nothing between, so the placing phase shows a real fraction, the waiting phase shows elapsed
  against NCBI's own RTOE capped at nine tenths, and past the estimate it stops filling and
  says so. The browser check calls `progressBar` with crafted states and asserts all of that.
- **The experimental cell line control is nearly inert, and now says so.** Measured: **1 of the
  registry's 24 (cell line x evidence axis) pairs is available**, AF22 chromatin, and that one
  needs the chromatin option switched on and the bigWigs installed. The field never touches a
  candidate, a gate, a preference or a rank. All of that availability data already reached the
  browser and none of it was drawn beside the control, which is why it read as doing nothing.
  The setup panel now shows per-axis status, source and note for whatever is typed, sharing
  `axisState` with the results banner so the two cannot word the same fact differently.
- **Transcript expression is wired in, GTEx only.** `sarna.expression` and
  `sarna.sources.json_http`; `GET /api/genes/{symbol}/expression`, a separate request from
  lookup so a GTEx outage cannot take gene lookup down. It never reorders the picker, and the
  browser check asserts the order is identical before and after the evidence loads. The
  evidence is split across two surfaces: a four-line answer beside the picker and a
  comparison of every transcript in the gene overview card, because a 310px column cannot
  hold both. `sarna.expression.organs` groups the 54 sample groups into 21 organs and the
  page shows the strongest few as symbols; only 8 organs get one, the rest are named, and a
  symbol means a median above zero with the group and value always beside it.
  **Corrected 9 October:** the caveat said GTEx "contains no cell lines", which is wrong.
  `gtex_v10` holds two cultured groups, donor-derived fibroblasts and EBV-transformed
  lymphocytes. AF22 and HeLa are still absent, and tests pin the corrected wording.
- **No factor occupancy anywhere.** AF22 chromatin *accessibility* is measured, which is not the
  same thing. `CELL_TYPE_OBSERVATION` for a transcription factor is still unproduced; ENCODE would
  populate HeLa, and AF22 has no public TF ChIP-seq at all. The asymmetry is permanent-ish: the
  tool will always know more about HeLa than about the disease-relevant line.
- ~~Licence audit.~~ **Done 2026-10-07.** Code is Apache-2.0; the bundled JASPAR subset is
  CC BY 4.0, verified against JASPAR's about page, the repository's CC BY-NC-SA label was wrong
  and more restrictive than the real terms. `NOTICE` lists every bundled asset and ships in the
  wheel. The unused `api` extra (fastapi/uvicorn) is removed.
  [ADR 0018](docs/adr/0018-apache-2-and-bundled-data-licensing.md).

## Delivery and updating

Versions are real now. `0.1.0.dev0` was frozen across every build handed to a colleague, so no
saved run could be traced to the code that produced it; releases carry numbers and tags, and
`tests/unit/test_release_metadata.py` fails if the wheel filename quoted in `README.md` or
`INSTALL.md` drifts from `pyproject.toml`.

**Measured, not assumed** (uv 0.12.23, 8 October 2026, two builds sharing a version string):
`uv tool upgrade` tracks the **git commit**, not the version number, so a git-installed tool
updates even when the version has not changed. For a wheel install it is a no-op that prints
`Nothing to upgrade`, which a colleague reads as "I am up to date" when they are not; the
instructions now say so. `docs/development/releasing.md` holds the table and the checklist.

**The outstanding defect is the branch, not the tooling.** `uv tool install ... @ git+https://...`
clones the **default branch**. `main` is 33 commits behind `publish-for-colleagues`, so
`uv tool upgrade` faithfully reinstalls stale code. Nothing a colleague runs can reach the new
work until the branch is merged. That merge and the tag push are the user's to authorise.

## Next

The first `src/sarna/app/` slice is implemented and exercised: gene-to-promoter lookup,
editable rules, visible guide predictions and ratings, hover explanations, reruns with diff,
local persistence and frozen exports. Repeat its browser check when changing UI paths that
carry scientific claims; use `AppRouter.dispatch` for offline route/store/source tests.

Continue the continuation plan: source/licence metadata reconciliation and
configuration versioning; generalise context readers and add pinned archived EPDnew, then DepMap,
FANTOM5 and ENCODE; methods/order exports and tested installers. The current diff compares a
new run against the previously displayed run; arbitrary saved-run comparison remains an extension.
No further scientific/default confirmation is required for that slice. Keep optional source
retrieval explicit and include cached/offline use.

The user's 7 October interface feedback refines milestone 6's acceptance criteria and is recorded
there: pre-assigned demos are not the deliverable, the promoter visualisation needs hover
explanations of what each element means, the rating and the predicted guide strand must be
visible per candidate, and every setting must be adjustable from the browser.

Low active-mark enrichment in AF22 is a contextual observation, not a direct transcription
measurement or proof of the cause of a null result. HeLa expression remains unverified. Keep
warnings separate from expression observations and do not make chromatin a candidate score.
