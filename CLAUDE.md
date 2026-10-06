# saRNA Designer — working notes for Claude

Python library that turns a promoter sequence into ranked, explained saRNA candidates. Scientific
traceability outranks clean code: every rule carries a citation, every verdict carries its reason.

**State:** the design engine works and is the only thing that works. `design_sarnas(promoter,
config)` generates, gates, ranks and explains, on either strand, from a pasted or
assembly-anchored promoter. 728 tests, 100% statement coverage, mypy strict clean.
No gene lookup, no API, no UI — Workflow A (gene symbol → promoter) does not exist yet.

## Commands

Use `.venv/bin/python`; there is a miniconda base environment that must stay untouched.

```bash
.venv/bin/python -m pytest -q                      # 728 tests, ~3 s
.venv/bin/python -m mypy src tests scripts         # strict
.venv/bin/python -m ruff check .
```

`ruff format` is deliberately not part of the workflow — the source is hand-wrapped near 88
columns while `line-length` is 100, so running it reformats most files. `ruff check` is the gate.

## Invariants

Violating one of these produces wrong science that type-checks and passes unrelated tests.

- **Positions inside a sequence are 1-based; genomic intervals are 0-based half-open.** Table 1
  speaks of "nucleotide 18", so `base_at` matches the paper. The two conventions never meet:
  one indexes within a sequence, the other locates a sequence on a contig.
- **TSS positions have no zero.** `+1` is the first transcribed base, `-1` the base before it.
  All arithmetic goes through `TssPosition`, which converts to a continuous offset first.
- **Promoters are held 5'→3' in the transcript's own orientation.** `Promoter.
  from_reference_sequence` flips a minus-strand reference once, at the boundary. After that the
  engine is strand-agnostic; the only strand arithmetic left is `Promoter.genomic_interval`,
  where offset 0 is the *highest* coordinate on the minus strand. This is defect B1's home.
- **Target sites are DNA; duplex strands are RNA.** CpG content is a property of the DNA locus
  and `cpg_count` refuses to answer for RNA.
- **Three outcomes, not two.** `not-applicable` means the rule could not be evaluated, and it
  neither excludes nor scores. A rule the preset omits produces no result at all. "Omitted",
  "could not be checked" and "failed" are three different statements about a candidate.
- **`legacy/` is immutable.** Preserved verbatim as the scientific ground truth; linting and
  formatting exclude it.
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
| What a preset is; legacy's two declared divergences | 0010 |

`docs/scientific/legacy-inventory.md` catalogues 17 defects (B1–B17) in the original scripts with
severities. Cite them by id. `GLOSSARY.md` fixes the vocabulary — use its terms rather than
coining synonyms.

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
- Measure claims rather than asserting them. `scripts/guide_prediction_divergence.py` is the
  pattern: enumerate exhaustively, write the number down.
- Run the engine and read its output before calling something done.

## Open decisions

- **No empirical validation exists.** `Designs.xlsx`, the only measured design list, is missing.
  "Correct" currently means "matches the paper and the script". Recovering that file would change
  what can be claimed, so it is worth asking about again.
- **Position 18/19 when target length ≠ 19.** Evaluators use the literal Table 1 position and
  report `not-applicable` if the target is shorter. For a 21-mer, position 19 stops being the 3'
  terminal base and `pos19-a`'s mechanism no longer holds. Undecided; keep length at 19.
- **3' overhangs are not modelled.** Candidates carry 19 nt cores. Wang's geometry adds a 3' UU
  overhang to both strands, and what a researcher orders is the 21-mer. Whether the overhang is
  synthetic (`UU`) or taken from flanking genomic sequence is a design decision for the export
  step, not yet made.
- **Assembly reconciliation is unimplemented.** ADR 0004 put liftover inside the bigWig adapter;
  that adapter does not exist, so no regulatory annotation runs yet (defect B15).

## Next

Phase 3: Ensembl client for gene → transcript → TSS → promoter, which makes Workflow A real and
lets promoters be spot-checked against an authoritative source. Then the API, then the UI.
