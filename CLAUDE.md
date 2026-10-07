# saRNA Designer — working notes for Claude

Python library that turns a promoter into ranked, explained saRNA candidates. Scientific
traceability outranks clean code: every rule carries a citation, every verdict carries its reason.

**State:** gene symbol → promoter → ranked candidates → interactive page all work.
`design_sarnas(promoter, config)` is the engine; `EnsemblClient` resolves a gene to a promoter;
`sarna.regulatory` annotates it; `sarna.view.html.render` writes a self-contained browser.
1,300 tests passed, eight live-network tests excluded; mypy strict and Ruff clean. Coverage was
not remeasured. **No HTTP API** — the library and its
scripts are the only interfaces.

## Accepted continuation decisions — 2026-10-07

Read `docs/development/continuation-plan.md` for the user's accepted continuation decisions;
they supersede corresponding open questions in these older working notes. The latest decision
includes [ADR 0013](docs/adr/0013-researchers-can-adjust-rule-roles.md): researchers may adjust supported
scientific rules between mandatory, preference and off in advanced custom configurations.
`pos19-a` stays mandatory in the default, while custom configurations may change its role.
Explain defaults, alternatives and consequences, and record all custom departures. This is an
accepted requirement; engine support is implemented, while the interactive editor is pending.
[ADR 0014](docs/adr/0014-local-delivery-and-wang-default.md) records the latest clarification:
local installation is the delivery target, Wang remains the default including its window gate,
and other papers inform optional recommendations rather than universal rules. This supersedes
the earlier proposal to change the default window to a preference.

## Commands

Use `.venv/bin/python`; there is a miniconda base environment that must stay untouched.

```bash
.venv/bin/python -m pytest -q                      # 1300 tests, ~17 s
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
scripts/locate_dataset_designs.py          # needs network; already run, its CSV is committed
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
| Researcher-adjustable scientific rule roles and transparent custom configurations | 0013 |
| Local delivery, tested installation paths and Wang as the default | 0014 |

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

Every sequence was built to pass Wang's gate, so the set is **range-restricted**: it cannot show a
gate too permissive, only one too strict. `CHAPTER4_CAVEATS` carries the rest.

`load_chapter4_locations()` adds recovered GRCh38 coordinates for 123 of the 125, with 519 bp of
context each, so the TSS windows and 3' flank preferences are evaluable — pass `locations=` to
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

## Implementation gaps and remaining decisions

- **Measured alternatives to existing inferred anchors.** Literature-defined anchors already
  exist in `sarna.genome.anchors`. The source adapter remains to be built:
  chapter 2 designed SCN1A against promoters mapped by Dong 2014 and Long 2008, and the Dong distal
  region — which produced 2 of the 3 SCN1A activators — has no Ensembl TSS within 5.5 kb. A tool
  that can only design relative to annotated TSSs could not reproduce this project's own hits.
  Preserve explicit uncertainty when adding measured alternatives. See
  `docs/scientific/scn1a-design-coordinates.md`.
- **Local application and installers.** `sarna.view.html.render` produces a self-contained page;
  a local Python service, interactive editor and tested installation paths remain to be built.
- **Position 18/19 when target length ≠ 19.** Evaluators use the literal Table 1 position and
  report `not-applicable` if the target is shorter. For a 21-mer, position 19 stops being the 3'
  terminal base and `pos19-a`'s mechanism no longer holds. Maximum is already 19; the library
  still supports 18 with explicit not-applicable results. Order export is fixed at 19 + dTdT.
- **3' overhangs are not modelled.** Candidates carry 19 nt cores; what gets ordered is a 21-mer.
  The chemistry question is *settled* — chapter 2's supplementary material shows **dTdT** on both
  strands, not ribo-UU — so what remains is building the order-ready export, not deciding what it
  should say. Largest gap before bench use.
- **Window alternatives and their costs.** As a gate at
  Wang's −1000…−100 it rejects 64% of candidates on a real promoter and 10 of the 14 known
  activators, including KLF4 at −45 and CEBPA (in human trials) at −1088. As a preference it would
  admit all 14 and lose the strongest filter the engine has. The rule now offers the thesis's own
  −2000…+500 window as a sourced option, with implementation restricted to an end of 0.
  Wang's gate remains the default; custom configurations can change its role or disable it.
- **`pos19-a` excludes CEBPA**, the duplex in MTL-CEBPA and the first saRNA in human trials.
  The gate stays conformant with Table 1 in the default. Custom roles are now supported and
  retain their departure from the frozen baseline in exports.
- **No factor occupancy anywhere.** AF22 chromatin *accessibility* is measured, which is not the
  same thing. `CELL_TYPE_OBSERVATION` for a transcription factor is still unproduced; ENCODE would
  populate HeLa, and AF22 has no public TF ChIP-seq at all. The asymmetry is permanent-ish: the
  tool will always know more about HeLa than about the disease-relevant line.
- **Licence audit.** Official JASPAR documentation states CC BY 4.0, contradicting the repository's
  CC BY-NC-SA label. Audit the exact bundled assets and correct metadata/attribution before
  distribution; do not assume unbundling is legally required. The precise MIT/Apache code licence
  form still needs recording; `LICENSE` remains undecided.

## Next

Follow the continuation plan. No further scientific/default confirmation is required for its
next slice. Complete browser checks, source/licence metadata reconciliation and configuration
versioning; generalise context readers and add pinned archived EPDnew, then DepMap, FANTOM5 and
ENCODE. Build frozen local run storage, the settings-to-rerun workbench path and tested installers.
Keep optional source retrieval explicit and include cached/offline use.

Low active-mark enrichment in AF22 is a contextual observation, not a direct transcription
measurement or proof of the cause of a null result. HeLa expression remains unverified. Keep
warnings separate from expression observations and do not make chromatin a candidate score.
