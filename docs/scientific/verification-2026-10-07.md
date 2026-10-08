# Scientific verification, 2026-10-07

What was checked, how, and what it found. The house rule applies throughout: verify by an
independent derivation, not by re-running the implementation.

## 1. The design rules against the published source

Wang J, Place RF, Portnoy V, et al. *Inducing gene expression by targeting promoter sequences
using small activating RNAs.* J Biol Methods 2015;2(1):e14. doi:10.14440/jbm.2015.39
([PMC4379447](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379447)).

Table 1 was read from the paper and compared with the shipped `wang-2015` preset, criterion by
criterion.

| Table 1 | Shipped default | |
|---|---|---|
| Distance from TSS −100 to −1000 bp | `tss-window` −1000 … −100 | ✔ |
| Target size 19 nt | `target-length` = 19 | ✔ |
| Target GC content 40–60 % | 0.40 / 0.60 | ✔ |
| CpG islands: avoid | `not-in-cpg-island`, mandatory | ✔ |
| Consecutive nucleotides ≤ 5 | `maximum-run-length` = 4, rejecting runs ≥ 5 | see below |
| Simple repeats (di/tri): avoid | `no-simple-repeats`, mandatory | ✔ |
| CpGs in target ≤ 1 | `maximum-cpg-count` = 1 | ✔ |
| 1st and 2nd nucleotides preferably G or C | `start-1-gc`, `start-2-gc`, preferences | ✔ |
| 18th nucleotide A or T | `pos18-at`, mandatory | ✔ |
| 19th nucleotide A | `pos19-a`, mandatory | ✔ |
| 20–23rd nucleotides preferably A or T | `flank-20-au` … `flank-23-au`, preferences | ✔ |
| Asymmetric thermodynamic stability | `thermodynamic-asymmetry`, mandatory | ✔ |

Every numeric default already has a test naming the paper, in
`tests/unit/test_catalog_wang2015.py`, so a change to a published threshold breaks a test that
cites its source.

**The run-length row is the one that needs a sentence.** Read literally, "≤ 5" permits a run of
exactly five; the implementation rejects it. This is not an oversight: the legacy scanner is the
declared ground truth ([ADR 0001](../adr/0001-scanner-is-ground-truth.md)), its `checkconsecutive`
compares five positions and so rejects runs of five or more, and the discrepancy between that code
and its own docstring is already recorded as **defect B5** in
[legacy-inventory.md](legacy-inventory.md), which resolves the reading in the same direction. The
measured cost is zero: at this threshold the rule excludes **none** of the 14 activators in
chapter 4's dataset ([empirical-validation.md](empirical-validation.md)). Loosening it by one
would admit more candidates, and the dataset is range-restricted, so it cannot say whether that
would be better. Left as it is, and flagged here rather than changed quietly.

## 2. The thermodynamics, by hand

The guide-strand prediction for SCN1A rank 1 (`saRNA-1549`, target `GTGAAGCCATAGATCTGTA`) was
recomputed from the published Xia et al. 1998 RNA nearest-neighbour ΔG°₃₇ parameters:

- 5′ end `GUGA`: GU −2.24 + UG −2.11 + GA −2.35 = **−6.70 kcal/mol**
- 3′ end `UGUA`: UG −2.11 + GU −2.24 + UA −1.33 = **−5.68 kcal/mol**
- E(3′) − E(5′) = **+1.02 kcal/mol**, which is what the model table prints.

The individual stack values match the published table (GC −3.42, GG −3.26, CG −2.36, GA −2.35,
GU −2.24, CA/UG −2.11, UA −1.33, AA −0.93). The 3′ end is the less stable; the antisense strand's
5′ end pairs with it; antisense is predicted as the guide. All three parameter sets agree here.

A separate test computes a two-base end by hand (GC −3.42, AA −0.93, difference +2.49) through the
whole stack, from `end_stability` to the JSON the browser receives.

## 3. The duplex, as drawn

The ladder computes a rung from the bases rather than assuming one. Verified in the browser for a
real candidate: 19 of 19 columns pair, both 3′ dTdT overhangs present, guide and passenger labelled
from the model's own verdict.

**This check was found to be vacuous and was fixed.** Planting a wrong pairing table
(`G→G`, `C→C`) still produced PASS, because the assertion only required *drawn* rungs to be real
pairs — and a wrong table draws fewer rungs rather than wrong ones. It now requires a rung exactly
where the bases pair, in both directions, and re-planting the same defect fails with
"column 0 pairs G with C (a real pair) but did not draw a rung".

## 4. Coordinates, against the reference genome

The strongest check available, and the only one that does not ask the engine to agree with itself.
For a finished design, each candidate's reported genomic interval is fetched from Ensembl
separately and required to contain the reported target sequence, reverse-complemented on the minus
strand. `tests/network/test_ensembl_live.py`, run live:

- **SCN1A** (minus strand): every passing candidate sits where it says it does.
- **CHD2** (plus strand): likewise, with the flip correctly *not* applied.
- TSS-relative spans agree with genomic coordinates, which are computed independently of them.

All 11 live tests pass, in 11 minutes. This is the test that could catch defect B1 — a sign error
in `genomic_interval`, `offset_of` or `span_of` produces candidates that are internally consistent,
pass every rule, export cleanly and point at the wrong DNA.

## 5. Defects found while verifying

| | Where | Consequence if unfixed |
|---|---|---|
| Telomere crash | `place_hit` clipped the low end of its query interval and not the high end | One alignment within 10 kb of a chromosome end raised out of `run_screening` and failed the whole screening |
| Cache poisoning | responses written non-atomically, read with a bare `json.loads` | An interrupted write left a truncated file; every later lookup of that gene raised, for ever, with nothing telling a colleague to delete it |
| Pool starvation | screenings shared the two-worker pool with design runs | Starting both scopes froze new designs for up to 30 minutes |
| Unreachable cache directory | `mkdir` raised in the constructor | A shared cache folder off the VPN stopped the tool starting |
| Vacuous pairing check | assertion was one-directional | A wrong pairing table would have drawn 19 crosses on a perfect duplex, silently |

## What this does not establish

Nothing here measures activation. The rules are reproduced faithfully from their source and the
arithmetic is correct; whether a duplex meeting them raises transcription is an experimental
question, and on the 125 assayed sequences in chapter 4 **no preset's ranking beats chance**
(AUC 0.41–0.53). That measurement is displayed in the application and is not softened here.
Sequence similarity found by off-target screening is likewise not measured binding.
