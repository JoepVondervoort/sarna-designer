# Legacy Script Inventory

Status: **Phase 0 / Discovery.** Living document, updated during requirements grilling.

## 1. Supplied material

| Artifact | Lines | Role | Notes |
|---|---|---|---|
| `saRNA Promoter Scanner.py` | ~560 | Promoter-sequence → 19-mer candidate design + rank | Attached in chat |
| `saRNA by Promoter Sequence.py` | ~560 | — | **Byte-identical duplicate** of the above. Not a second implementation. |
| `EpiMarker checker .py` | 837 | Gene symbol → Ensembl TSS → ±30 kb → locate *known* saRNAs → histone ChIP enrichment + RandomForest | On disk at `/home/jvon/EpiMarker checker .py` |
| `af22_data/*.bw` | 6 files, ~890 MB | AF22 histone ChIP-seq bigWigs, GSE108506, 50 bp bins | On disk, already downloaded |
| `20250428_saRNA_designs_with_epigenomic_marker_scanning/` | empty | Output dir from a prior run | — |
| `Designs.xlsx` | **MISSING** | Input to EpiMarker: per-gene columns of 19-mers labelled positive/non control | Script aborts without it |

Two scripts, not three: the design engine (Scanner) and the epigenomic annotation/validation
script (EpiMarker) are **separate pipelines with no shared code and conflicting rules**.

## 2. Terminology (to be ratified during grilling)

| Term | Working definition in this project |
|---|---|
| Target | The 19 nt promoter-strand (sense) subsequence the duplex is designed against |
| Passenger / sense | Strand identical to the target; printed 5'→3' by the legacy Scanner |
| Guide / antisense | Complement of the target; intended RISC-loaded strand |
| Duplex | 19 bp core + 3' `UU` overhang on both strands (i.e. delivered as 21-mers) |
| Promoter window | Legacy Scanner: *whatever sequence is pasted* — no TSS-relative window exists yet |
| POS_CTRL / NON_CTRL | EpiMarker labels for empirically active vs inactive known saRNAs |

Note: the legacy Scanner has **no concept of a TSS or of genomic coordinates at all.** It
scans a bare pasted string. All TSS-relative logic lives only in EpiMarker.

## 3. Design rules — Promoter Scanner (the design engine)

Ordering below is the order applied in the legacy main loop.

| # | Rule | Location | Current default | Range | Rationale | Evidence class | Configurable? | Effect on ranking |
|---|---|---|---|---|---|---|---|---|
| R1 | Candidate length = 19 nt | `maketargets` | 19 | 18–21 typical | siRNA/saRNA duplex core | Established | Yes | Defines the window; none directly |
| R2 | Scan start offset `i = 3`, end `len-4` | main loop | 3 / −4 | — | **Unexplained.** `−4` is needed so `checkflanks` can read +22; the leading `3` has no justification and silently drops the first 3 candidates | Artefact / likely bug | Remove | Drops candidates |
| R3 | No homopolymer run | `checkconsecutive` | run ≥ 5 rejected | 3–6 | Runs cause synthesis/off-target/structure problems | Heuristic | Yes | Hard filter |
| R4 | `endstability ≥ 0` | `endstability` | ≥ 0 | — | Thermodynamic asymmetry: 5' end of *target* more stable ⇒ 5' end of *guide* less stable ⇒ guide preferentially loaded | Established (Schwarz 2003, Khvorova 2003) | Yes (threshold) | Hard filter **and** scored (`stab`) |
| R5 | GC fraction | `countGC` | **0.40 – 0.65** | 0.3–0.7 | Duplex stability / RISC loading | Empirical | Yes | Hard filter |
| R6 | Position 19 = `A` | `checkpos19` | required | — | ⇒ guide 5' nt is `U`; U/A at guide 5' favours AGO2 loading | Established-ish (Ui-Tei 2004) | Yes (on/off) | Hard filter |
| R7 | Position 18 ∈ {`A`,`U`} | `checkpos18` | required, A=2 / U=1 | — | A/U-rich 3' target end = low-stability guide 5' end | Empirical | Yes | Hard filter **and** dominant score term |
| R8 | End-base asymmetry | `endbase` | required True | — | Single-base proxy for R4 using half the AA/GG stacking energy | Heuristic, crude | Yes | Hard filter |
| R9 | No overlap with any CpG site | `insite` | required | — | Avoid CpG methylation-sensitive targets | Heuristic, project-specific | Yes | Hard filter |
| R10 | Not inside a CpG island | `inisland` | **COMMENTED OUT in main** | — | Takai & Jones island definition | Heuristic | Yes | Currently inert |
| R11 | Position 7 = `U` | `checkpos7` | scored only | — | Guide position-13 preference / cleavage-site geometry | Empirical | Yes | Score term only |
| R12 | 3'-flanking A/U count | `checkflanks` | 4 nt window, scored | 0–8 | A/U-rich flank aids unwinding | Heuristic | Yes | Score term only |

### Scoring formula (legacy)

```
score = stab*stabilitydiff + ind18*indicator18 + ind7*indicator7 + flank*flankscore
        0.01 * ~[0..10]      1000 * {1,2}        1000 * {0,1}       1 * {0..4}
```

**This is not a calibrated numeric score — it is a lexicographic tier system in disguise.**
Because `indicator18 ∈ {1,2}` (0 is already hard-filtered) and `indicator7 ∈ {0,1}`, the
weights 1000/1000 create four discrete tiers (1000, 2000, 2000, 3000 — note tiers *collide*:
`ind18=2, ind7=0` ties `ind18=1, ind7=1` at 2000), with `flankscore` as the within-tier
discriminator and `stab` as a final sub-unit tiebreak. The 2-decimal-looking output score is
pseudo-precision. Candidate for explicit tiered ranking in the new engine.

### Nearest-neighbour thermodynamics

10 Gibbs values hardcoded (`AAE`…`GGE`), sourced from "Wiki: Nucleic Acid Thermodynamics".
`endstability` sums NN stacking energies over the first and last `end`(=4) nt.

## 4. Design rules — EpiMarker (conflicting second implementation)

| Rule | Location | Value | Conflict with Scanner |
|---|---|---|---|
| Length = 19 | `validate_dna` | 19 | consistent |
| GC % | `validate_dna` | **40 – 60 %** | **CONFLICTS with Scanner's 40–65 %** |
| Score | `compute_score` | `0.7*GCscore + 0.3*distanceScore`, ideal = **−500 bp** from TSS | **Completely different, incompatible scoring scheme** |
| Promoter window | `get_region` | **±30 000 bp** around TSS | Scanner has no window at all |
| TSS | `get_gene_and_tss` | Ensembl canonical transcript only | — |

EpiMarker *does not generate* candidates — it locates user-supplied 19-mers by string search.

## 5. Guide-strand prediction (to be preserved)

Legacy behaviour: **guide = antisense, always.** There is no per-candidate guide/passenger
*choice*. Instead, asymmetry is enforced as an *admission criterion*: R4, R6, R7 and R8
together only admit targets whose thermodynamic/sequence asymmetry implies the antisense
strand will be the loaded guide. Candidates where the sense strand would win are rejected
rather than relabelled.

This is a defensible reading of Schwarz/Khvorova and must be preserved as-is. The new engine
should expose the underlying features (`stabilitydiff`, `E5prime`, `E3prime`, end bases,
position 18/19 identity) rather than only the boolean outcome.

## 6. Biological assumptions (all currently implicit)

1. Guide is always the antisense strand.
2. One gene → one canonical transcript → one TSS → one promoter (EpiMarker).
3. Genome assembly is whatever Ensembl REST currently defaults to (GRCh38) — **never pinned**.
4. AF22 bigWig coordinates are assumed compatible with that assembly — **unverified, and a
   likely GRCh37/GRCh38 mismatch given GSE108506's publication date.** Blocking question.
5. 19 nt + UU overhang is the correct delivery format.
6. A single-bp bigWig sample represents the signal over a 19 nt target.
7. CpG overlap is undesirable for a saRNA target.
8. siRNA-derived asymmetry rules transfer to promoter-targeted saRNA. Widely assumed in the
   saRNA literature but mechanistically weaker than for canonical RNAi.

## 7. Bugs and defects found

Severity: **S1** = wrong scientific output, **S2** = wrong behaviour, **S3** = hygiene.

| ID | Sev | Location | Defect |
|---|---|---|---|
| B1 | **S1** | EpiMarker `get_region` + `find_sarna` | For minus-strand genes Ensembl is asked for `strand=-1`, returning the **reverse complement** of the region. `idx = seq.find(dna)` therefore counts from the opposite end, yet `abspos = start + idx` treats it as a forward offset, and `rel = abspos − TSS` keeps the plus-strand sign convention. **All minus-strand positions and all TSS-relative distances are wrong.** Exactly the strand-coordinate class of bug to test explicitly. |
| B2 | **S1** | Scanner NN energy table | Reverse-complement-equivalent dinucleotides are mis-paired. In NN convention `AC ≡ GT`, `UC ≡ GA`, `AG ≡ CT`. The code maps `AC→CAE`, `UC→CTE`, `AG→GAE` — i.e. three pairs assigned to the wrong parameter. Biases every `endstability` value and therefore filter R4 and the score. **Needs your decision: fix, or preserve for traceability?** |
| B3 | **S1** | EpiMarker `find_sarna` return | `return rel, compute_score(gc, rel) if ok else (rel, 0)` parses as `return rel, (score if ok else (rel, 0))`. When GC validation fails the "score" is a **tuple**, not a number. Silently corrupts downstream scoring/ML. |
| B4 | **S2** | EpiMarker `enrichment_analysis` | Signal sampled at a **single base pair** (`TSS+p` … `+1`) from 50 bp-binned bigWigs, and not averaged over the 19 nt target. Also uses `TSS+p` which compounds B1 on minus strands. |
| B5 | **S2** | Scanner `checkconsecutive` | Docstring claims "22-nucleotide subsequence" and "four or more consecutive identical bases"; the loop actually detects runs of **five** within the 19-mer. Code and documented rule disagree. |
| B6 | **S2** | Scanner main loop | `i = 3` start offset silently discards the first three candidate windows with no stated reason (R2). |
| B7 | **S2** | Scanner `CpGislanddetector` | Relies on **list aliasing** — the in-code comment calls it "a strange feature of Python" and the merge logic depends on the alias. Fragile. Also `CpG > 7` where Takai & Jones specify ≥ 7, and `CpGexp` is recomputed with `str.count` every iteration (O(n·200)) while GC is updated incrementally. |
| B8 | **S2** | Scanner main | `inisland` is **commented out** — R10 does not run. Intentional or accidental? |
| B9 | **S2** | Scanner `endbase` / `endstability` / `complement` | Error branches `print(...)` and then **continue with a stale or undefined `energy`**. `endbase` would raise `UnboundLocalError` on the first call; later calls silently reuse the previous value. Violates fail-fast. |
| B10 | **S3** | EpiMarker `setup_environment` | Hardcoded Windows path `C:\Users\joep-\...` + `os.chdir`. Cannot run in WSL. |
| B11 | **S3** | EpiMarker | `pip install` invoked at runtime from inside the script (`ensure_plotly`). |
| B12 | **S3** | Scanner | `from pylab import *` — unused, and star-imports numpy/matplotlib names over builtins. |
| B13 | **S3** | Scanner `complement` | Named/documented as "reverse complement" in its error message but correctly computes a plain complement. |
| B14 | **S3** | Both | Everything is index-arithmetic `while` loops with manual counters; no functions have type hints; `promotername` is a hardcoded module constant. |
| B15 | **S1** | EpiMarker — Ensembl REST ↔ `af22_data/*.bw` | **Assembly mismatch, verified 2026-10-06 (§12).** TSS coordinates come from Ensembl REST, which defaults to **GRCh38**; the AF22 bigWigs are **GRCh37/hg19**. Coordinates are also named differently (`1` vs `chr1`). Every histone-signal lookup therefore reads the wrong locus. This was the blocking open question — it is now closed, and the answer is that the mismatch is real. |

## 8. Code to preserve vs refactor

**Preserve the science, with regression tests pinning current output:**
`endstability` / NN table (pending B2 decision), `checkpos7/18/19`, `checkflanks`,
`checkconsecutive` semantics, `countGC`, `endbase`, `CpGsitedetector`, the rule *ordering*,
and the asymmetry-as-admission-criterion model of guide prediction.

**Refactor (mechanical, behaviour-preserving):**
loop style → comprehensions/generators; module constants → validated config schema;
`print` → structured logging/errors; duplicated 60-line NN if-chain (appears **twice**
verbatim inside `endstability`) → single lookup table.

**Rewrite deliberately, with documented behaviour change:**
B1 strand/coordinate handling, B3 return bug, B4 signal aggregation, EpiMarker's
`compute_score`, `setup_environment`, `ensure_plotly`.

**Retire:** the duplicate second script; `from pylab import *`.

## 9. Dependencies observed

`requests`, `pandas`, `numpy`, `scipy`, `pyBigWig`, `matplotlib`, `seaborn`, `plotly`,
`kaleido`, `openpyxl`, `scikit-learn` (RandomForest), `pylab`/matplotlib.

External services: Ensembl REST (`xrefs/symbol`, `lookup/id`, `sequence/region`),
NCBI GEO FTP (bigWig download). No assembly pinned, no caching, no rate limiting,
no retries on either.

## 10. Environment (verified 2026-10-06)

Python 3.11.4 · Node v20.20.0 · git 2.34.1 · Linux 6.18.40.1-microsoft-standard-WSL2 ·
working dir `/home/jvon` on the Linux filesystem (good — not `/mnt/c`) · not yet a git repo.

---

## 11. Canonical source: Wang et al. 2015

**Wang J, Place RF, Portnoy V, Huang V, Kang MR, Kosaka M, Ho MKC, Li L-C.**
*Inducing gene expression by targeting promoter sequences using small activating RNAs.*
J Biol Methods, 2015. PMC4379447. Supplementary **File S1** = the Excel design tool.

User-confirmed (2026-10-06) as the guideline the legacy designs were based on.
`saRNA Promoter Scanner.py` is confirmed as the ground-truth implementation to pin.

### Wang 2015 Table 1 vs. the legacy Scanner

| Wang 2015 Table 1 criterion | Scanner implementation | Verdict |
|---|---|---|
| Distance from TSS: −100 to −1000 bp | **absent — no TSS concept at all** | **MISSING** |
| CpG islands: avoid | `inisland` written but **commented out** | **INERT** (confirms B8 is a defect) |
| Target size: 19 nt | 19 nt | conforms |
| Target GC content: **40–60 %** | **40–65 %** | **DIVERGES** |
| Consecutive nucleotides ≤ 5 | rejects runs ≥ 5 | conforms (B5 was a *docstring* error; code is right) |
| Simple repeats (di-/tri-nucleotide): avoid | **absent** | **MISSING** |
| CpGs in target ≤ 1 | **zero tolerance** (`insite`) | **DIVERGES — stricter** |
| 1st and 2nd nucleotides preferably G or C | **absent** | **MISSING** |
| 18th nucleotide: A or T | required (`checkpos18`) | conforms |
| 19th nucleotide: A | required (`checkpos19`) | conforms |
| 20th–23rd preferably A or T | `checkflanks`, soft-scored | conforms |
| Thermodynamic asymmetry (NN model, antisense 5' less stable) | `endstability` + `endbase` | conforms in intent; see B2 |
| — *no such rule in the paper* — | **position 7 = `U`, weight `ind7 = 1000`** | **UNSOURCED** |

### Consequences

1. **`ind7 = 1000` has no provenance.** Position 7 appears nowhere in Wang 2015. It is
   nonetheless one of the two dominant ranking terms. Possibly a half-remembered Reynolds
   (2004) siRNA criterion ("U at position 10"). Must be reclassified as *exploratory* unless
   a source is produced.
2. **Wang 2015 assigns no weights whatsoever.** The `1000 / 1000 / 1 / 0.01` weighting is
   entirely local invention, not paper-derived. The paper's own tool exposes only an
   undisclosed "overall quality score".
3. **On GC content, EpiMarker (40–60 %) was right and the Scanner (40–65 %) is the outlier.**
4. **Three Table 1 rules are absent** from the Scanner: TSS distance window, di-/tri-nucleotide
   repeats, and the G/C preference at positions 1–2.
5. The paper classifies rules as a **mix** of mandatory ("avoid") and preference ("preferably")
   with no stated hierarchy. The new engine must make this hierarchy explicit.

### Honest-expectations statement (for the UI)

Wang 2015 reports that of designed duplexes, **15–20 % significantly activate** the target
gene, **~70 % have no effect**, and **5–10 % suppress** it. The application must surface this
base rate rather than implying a high-scoring candidate is likely to work.

### Product requirement captured

User (2026-10-06): the paper "explained why the score was like that… It should be like that in
our new application as well." ⇒ **per-rule score contribution and its rationale must be
visible in the UI for every candidate**, not just a final number. Feeds the Explainable
Ranking requirement.

---

## 12. Genome assembly of the AF22 bigWigs — RESOLVED 2026-10-06

Previously the blocking open question ("are the AF22 bigWigs on the same assembly as Ensembl's
default?"). **They are not.**

Method: chromosome sizes were read directly from the bigWig headers with `pyBigWig` and
compared against the known GRCh37 and GRCh38 lengths.

| Chromosome | Observed in all 6 `.bw` | GRCh37/hg19 | GRCh38 |
|---|---|---|---|
| `chr1` | 249,250,621 | **249,250,621** ✅ | 248,956,422 |
| `chrX` | 155,270,560 | **155,270,560** ✅ | 156,040,895 |

All six files (GSM2902405–GSM2902410) agree: 25 contigs, UCSC-style `chr`-prefixed names,
**GRCh37/hg19**. This is consistent with GSE108506's 2017 deposition date.

### Two distinct incompatibilities

1. **Assembly.** Ensembl REST (`rest.ensembl.org`) serves **GRCh38** by default. The legacy
   EpiMarker script pins no assembly and so silently received GRCh38 TSS coordinates, then
   used them to index hg19 signal tracks.
2. **Contig naming.** Ensembl returns `1`, `2`, `X`; the bigWigs expect `chr1`, `chr2`, `chrX`.
   A lookup with the Ensembl name raises or returns empty rather than wrong data — so this one
   fails loudly, whereas the assembly mismatch fails **silently**.

### Scientific consequence

GRCh37→GRCh38 coordinate shifts are locus-dependent and range from zero to several megabases.
The pipeline samples 50 bp-binned signal in a ±30 kb window, so any shift beyond a few hundred
base pairs puts the read outside the intended region entirely. **Every regulatory annotation
produced by the legacy EpiMarker script must be treated as unreliable**, independently of B1
(minus-strand coordinates) and B4 (single-bp sampling), both of which compound it.

This does not affect the Scanner: it has no coordinate concept at all (§3), which is why the
design path is unaffected and only the annotation path is implicated.

### Not yet decided

How the new engine reconciles the two — pinning the GRCh37 Ensembl endpoint, lifting the
bigWigs over to GRCh38, or lifting query coordinates per-request — is an open architecture
decision requiring an ADR. Whichever is chosen, the engine must **pin and assert the assembly
explicitly at every genomic boundary** rather than relying on a service default, and must
normalise contig naming at the bigWig adapter.

### Addendum — the mitochondrion is *not* shared between Ensembl GRCh37 and UCSC hg19

Validating the project's contig-length tables against the AF22 headers (2026-10-06) found all
**24 nuclear contigs identical**, confirming GRCh37, with one discrepancy:

| Contig | AF22 bigWigs (UCSC hg19) | Ensembl GRCh37 |
|---|---|---|
| `chrM` / `MT` | 16,571 bp (NC_001807) | 16,569 bp (rCRS) |

These are **different sequences**, not the same sequence measured differently, so mitochondrial
coordinates do not transfer between them under any liftover. The bigWig adapter therefore
rejects `MT` explicitly rather than treating hg19 and GRCh37 as synonymous there.

No scientific cost: mitochondrial DNA is not packaged into nucleosomes, so histone ChIP signal
is not interpretable on it, and saRNA targets are nuclear promoters. The point is recorded
because "GRCh37 ≡ hg19" is a tempting shorthand that is wrong in exactly this one place.
