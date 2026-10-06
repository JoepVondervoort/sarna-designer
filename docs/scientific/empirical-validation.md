# Empirical validation against measured activity

**What this establishes.** Wang et al. 2015's sequence gates exclude none of the 14 saRNAs known
to activate their target. The legacy script's gates exclude five of them. No preset's *ranking*
separates activators from non-activators, which is the expected result on this dataset and not
evidence against the ranking.

Produced by `scripts/validate_against_chapter4.py` against
`sarna.validation.dataset.load_chapter4_dataset`. Reproducible; the seed is printed with every
p-value.

## What the dataset can and cannot answer

125 target sequences from thesis chapter 4, 14 with reproducible ≥1.5-fold activation and 111
without. Three limits decide what any number here means, and they are not incidental.

**The gate is largely untestable, in one direction.** Every sequence was designed to satisfy
Wang's structural criteria, so the set contains almost nothing the gate should reject. It cannot
show a gate to be too permissive. It *can* show a gate to be too strict — if a sequence that
demonstrably worked is excluded, that is a real finding, and it is the finding this page leads
with.

**Sequences arrive without their promoters.** The chapter tabulates target sites, not coordinates.
So the TSS windows, CpG-island overlap and the 3' flank preferences report not-applicable and
contribute nothing. Everything measured below is the behaviour of the *sequence-intrinsic* rules
alone. Recovering the coordinates — mapping each target back to the promoter it was designed
against — would lift this limitation and is the single most valuable addition to the dataset.

**Cell type dominated the experiment.** AF22 neural stem cells produced no activators at any
locus, including the P21 control that works reliably in HeLa. No sequence rule can represent that,
so a rule fitted to these labels should not be expected to transfer between cell types.

## Wang's gates are well calibrated; the legacy script's are not

Cost of each gate, measured as known activators excluded:

| Gate | Setting | Activators lost | Others lost |
|---|---|---|---|
| GC content | 40-60 % (Wang) | **0** of 14 | 2 of 111 |
| GC content | 40-65 % (legacy) | 0 of 14 | 0 of 111 |
| GC content | 30-52 % (Reynolds) | 2 of 14 | 28 of 111 |
| GC content | 45-55 % | 8 of 14 | 47 of 111 |
| CpG count | ≤ 1 (Wang) | **0** of 14 | 0 of 111 |
| CpG count | 0 (legacy, and `strict` as first written) | **4** of 14 | 19 of 111 |
| Homopolymer run | ≤ 4 (Wang) | **0** of 14 | 4 of 111 |
| Homopolymer run | ≤ 3 | 3 of 14 | 23 of 111 |
| Position 18 A/T | Wang | **0** of 14 | 2 of 111 |
| Position 19 A | Wang | 1 of 14 | 11 of 111 |

Wang's published settings are the only ones in the table that cost nothing. That is a genuine, if
narrow, positive result for the default preset: on the one dataset available, every one of its
sequence gates admits every duplex known to work.

### The CpG gate points the wrong way

The legacy script tolerates **zero** CpG in a target. That excludes P53, SCN1A-10, CEBPA and KLF4
— four of the fourteen known activators, including two of this project's own SCN1A duplexes.

The direction is what makes it interesting. Activators carry *more* CpG than non-activators, not
less: mean 0.29 against 0.17, and 71 % of activators have no CpG against 83 % of non-activators.
A small difference on 14 sequences, so not a finding on its own — but it removes any empirical
case for tightening the gate, and tightening it has a measured cost of 29 % of known activators.

Mechanistically this is unsurprising in hindsight. saRNA targets sit in promoters; promoters are
CpG-rich; and the rule that matters for methylation-related silencing is CpG *island* overlap,
which is a property of the surrounding region rather than of 19 bases. The island rule already
covers that, which leaves the per-target CpG count doing little except excluding good candidates.

**Acted on.** The `strict` preset was first written with `cpg-count: 0`, on the reasoning that
fewer CpGs must be safer. The harness caught it on the first run: that version excluded 5 of 14
activators. `strict` now tightens only location criteria, which this dataset cannot test, and
leaves every sequence criterion at Wang's value. The rejected tightenings are recorded in the
preset's own `divergences` so the next person does not reinstate them.

### The legacy script's dominant ranking term is anti-correlated with activity

`pos7-u` has no source (finding B16 in the legacy inventory) and the original script weighted it at
1000, making it one of two dominant ranking terms. On this dataset it passes **43 % of activators
against 60 % of non-activators**, a 18-point separation in the wrong direction.

With 14 activators that is six sequences against eight, so it is not proof the rule is harmful.
It is the complete absence of support, in the only direction that could have rescued it. A rule
with no literature, no mechanism and an anti-correlation in the only data available has nothing
left to stand on. Recorded as finding B18.

`start-2-gc` behaves similarly at −10 pp across every preset, and it is Wang's own. Worth noting
and not worth acting on: it is a preference, so it costs a point rather than admission, and 14
sequences cannot overturn a published criterion.

## The most clinically advanced saRNA fails the default gate

**CEBPA is excluded by `pos19-a`** in every preset, because its position 19 is C.

This is not an obscure sequence. It is the duplex in MTL-CEBPA (Reebye et al. 2018), the first
saRNA to enter human clinical trials, for liver cancer. It activates its target, it has been
through Phase I, and Wang et al. 2015 Table 1 would have rejected it at the design stage.

Three readings, and the project does not choose between them on 14 sequences:

1. **Position 19 A is over-stated in Table 1.** It is one criterion among many, derived from a
   limited panel, and the mechanism it serves — an A at target position 19 puts a U at the guide's
   5' end, which AGO2 prefers — is a bias rather than a requirement. Ui-Tei et al. 2004 state the
   same criterion more weakly, as A *or* U at the guide's 5' end.
2. **CEBPA works for a reason position 19 does not capture.** Its 57.9 % GC and GC-rich 3' end
   make it unusual among the activators in several ways at once.
3. **Table 1 is right and CEBPA is an exception.** 13 of 14 activators do carry A at position 19.

What the tool does about it: nothing automatic. `pos19-a` stays a mandatory gate, because that is
what the paper says and conformance with the paper is the default preset's job. But this is now a
documented, named case where the gate excludes a working duplex, which is exactly the kind of
thing a researcher should see before trusting a rejection. A future option making `pos19-a` a
preference would be defensible; it is not added speculatively.

## No ranking separates the classes

| Preset | AUC | p | Gate rejections | Activators lost |
|---|---|---|---|---|
| `wang-2015` | 0.463 | 0.71 | 25 / 125 | 1 |
| `strict` | 0.463 | 0.71 | 25 / 125 | 1 |
| `exploratory` | 0.527 | 0.37 | 25 / 125 | 1 |
| `legacy-scanner` | 0.414 | 0.86 | 41 / 125 | 5 |

AUC is the probability that a randomly chosen activator outranks a randomly chosen non-activator,
ties counted as half, so 0.5 is a coin flip. Every preset is indistinguishable from chance, and
three of the four fall slightly below it.

This is the expected result and it is worth being precise about why, because "the ranking does not
work" is the wrong conclusion. The ranking orders candidates by how many *preferences* they meet.
The dataset was built by applying those same criteria and keeping what passed, so it contains
almost none of the variation the ranking exists to order — and of the seven preferences, four
(the 3' flank positions) and one (the preferred TSS window) cannot be evaluated at all without
promoter coordinates. Two preferences are doing all the work. A near-0.5 AUC under those
conditions tells you about the dataset, not about the ranking.

### `exploratory` scoring highest is circular, not encouraging

`exploratory` has the best AUC, and `pos10-g` and `central-stability` are the two best-separating
rules in the entire catalogue (+25 pp and +18 pp). Both were derived **from this dataset**.

The +25 pp for `pos10-g` is the same +24.8 pp from the motif analysis, on the same 125 sequences.
Measuring it here is not validation; it is reading the training set back. And even with two rules
fitted to the data it is scored against, `exploratory` reaches p = 0.37 — which is the clearest
possible statement of how little signal there is to fit.

## What would make this answerable

- **Recover the promoter coordinates.** Would activate five of the seven preferences and the
  location gates, and is the one improvement that needs no new experiments.
- **Hold out the chapter 4 sequences.** `pos10-g` and `central-stability` cannot be assessed on
  the data that produced them. Six published activators are available as an independent set, but
  six is too few.
- **Test one hypothesis prospectively.** Duplexes differing at position 10 and nowhere else would
  make the existing data decisive, because the grid collapses from 76 cells to one.
- **Record the negatives per cell line.** "Did not activate" currently pools HeLa and AF22, and
  the chapter's own strongest effect is that these differ. Separating them would let a cell-type
  term enter the model at all.
