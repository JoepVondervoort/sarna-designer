# saRNA Designer — working notes for Claude

Python library that turns a promoter into ranked, explained saRNA candidates. Scientific
traceability outranks clean code: every rule carries a citation, every verdict carries its reason.

**State:** gene symbol → promoter → ranked candidates → interactive page all work.
`design_sarnas(promoter, config)` is the engine; `EnsemblClient` resolves a gene to a promoter;
`sarna.regulatory` annotates it; `sarna.view.html.render` writes a self-contained browser.
1125 tests, 100% statement coverage, mypy strict clean. **No HTTP API** — the library and its
scripts are the only interfaces.

## Commands

Use `.venv/bin/python`; there is a miniconda base environment that must stay untouched.

```bash
.venv/bin/python -m pytest -q                      # 1125 tests, ~15 s
.venv/bin/python -m pytest -m network              # 8 live-Ensembl tests, excluded by default
.venv/bin/python -m pytest -m tracks               # needs the AF22 bigWigs; skips if absent
.venv/bin/python -m mypy src tests scripts         # strict
.venv/bin/python -m ruff check .
```

`ruff format` is deliberately not part of the workflow — the source is hand-wrapped near 88
columns while `line-length` is 100, so running it reformats most files. `ruff check` is the gate.

Scripts that measure rather than assert, which is the house pattern:

```bash
scripts/validate_against_chapter4.py       # every preset against 125 assayed sequences
scripts/positional_motif_significance.py   # ~2 min; permutation test over the motif grid
scripts/locate_scn1a_designs.py            # needs network; maps the thesis's designs to GRCh38
scripts/af22_chromatin_at_targets.py       # needs network + AF22 tracks; controls included
scripts/build_promoter_view.py SCN1A --html out.html    # needs network the first time
```

## Invariants

Violating one of these produces wrong science that type-checks and passes unrelated tests.

- **Positions inside a sequence are 1-based; genomic intervals are 0-based half-open.** Table 1
  speaks of "nucleotide 18", so `base_at` matches the paper. The two conventions never meet:
  one indexes within a sequence, the other locates a sequence on a contig.
- **TSS positions have no zero.** `+1` is the first transcribed base, `-1` the base before it.
  All arithmetic goes through `TssPosition`, which converts to a continuous offset first.
- **Promoters are held 5'→3' in the transcript's own orientation.** `Promoter.
  from_reference_sequence` flips a minus-strand reference once, at the boundary. After that the
  engine is strand-agnostic. **All strand arithmetic lives in three methods on `Promoter`** —
  `genomic_interval`, `offset_of`, `span_of` — where offset 0 is the *highest* coordinate on the
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
| Guide prediction: 13% of candidates flip between tables | `guide-prediction-divergence.md` |

`GLOSSARY.md` fixes the vocabulary — use its terms rather than coining synonyms.

## The empirical dataset

`sarna.validation.dataset.load_chapter4_dataset()` holds 125 assayed target sequences from thesis
chapter 4, 14 of which activated their gene. It is the only thing here that can show a rule to be
**wrong** rather than merely unconventional, and it has already refuted one preset.

Two limits decide what any number from it means, and both are in `CHAPTER4_CAVEATS`:

- Every sequence was built to pass Wang's gate, so the set is **range-restricted**. It cannot show
  a gate too permissive. It *can* show one too strict, and that is the finding to look for.
- The sequences arrive **without their promoters**, so five of seven preferences and all the
  location gates report not-applicable. Recovering the coordinates is the single most valuable
  addition; `scripts/locate_scn1a_designs.py` does it for the 15 SCN1A designs.

Headline results: Wang's published sequence gates exclude **none** of the 14 activators, while
every alternative setting measured excludes at least one. No preset's *ranking* beats chance
(AUC 0.41–0.53) — expected, given the range restriction. `exploratory` scores highest because two
of its rules were fitted to this very data, which the tests assert so the circularity stays visible.

## How to work here

- **Ask before defaulting.** Scientific choices — a threshold, a definition, which error to
  prefer — are the user's to make. Present the recommendation with its cost, then wait.
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
- Live-API tests earn their keep. Two Ensembl bugs — `mane=1` being required, and a 400 rather
  than a 404 for an unknown symbol — were invisible to recorded payloads.

## Open decisions

- **No HTTP API or hosted UI.** `sarna.view.html.render` produces a self-contained page, which
  covers the "show me the promoter" need without a server. Whether a FastAPI layer is still wanted
  (ADR 0005 anticipated one) is undecided.
- **Position 18/19 when target length ≠ 19.** Evaluators use the literal Table 1 position and
  report `not-applicable` if the target is shorter. For a 21-mer, position 19 stops being the 3'
  terminal base and `pos19-a`'s mechanism no longer holds. Undecided; keep length at 19.
- **3' overhangs are not modelled.** Candidates carry 19 nt cores. Wang's geometry adds a 3' UU
  overhang to both strands, and what a researcher orders is the 21-mer. Whether the overhang is
  synthetic (`UU`) or taken from flanking genomic sequence is a design decision for the export
  step, not yet made. This is the largest gap before the tool is usable at the bench.
- **`pos19-a` excludes CEBPA**, the duplex in MTL-CEBPA and the first saRNA in human trials.
  Documented, not worked around: the gate stays conformant with Table 1. Making it a preference
  would be defensible and has not been done speculatively.
- **No factor occupancy anywhere.** AF22 chromatin *accessibility* is measured, which is not the
  same thing. `CELL_TYPE_OBSERVATION` for a transcription factor is still unproduced; ENCODE would
  populate HeLa, and AF22 has no public TF ChIP-seq at all. The asymmetry is permanent-ish: the
  tool will always know more about HeLa than about the disease-relevant line.
- **Licence.** Bundled JASPAR matrices are CC BY-NC-SA 4.0, so shipping them makes the package
  non-commercial. `LICENSE` still says undecided. See `docs/scientific/tf-binding.md`.

## Next

Highest value first: recover promoter coordinates for the rest of the chapter 4 dataset, which
makes the location rules testable against measured activity; then the 3' overhang decision and an
order-ready export; then ENCODE occupancy for HeLa. An HTTP API only if the user wants one.

One open scientific question worth putting to the user: chapter 4's SCN1A screens appear to have
been run in undifferentiated AF22, whose SCN1A promoter carries no active mark. If so, those null
results are predicted by the mechanism rather than informative about the duplexes, and the
comparison belongs in differentiated neurons. See `docs/scientific/af22-chromatin.md`.
