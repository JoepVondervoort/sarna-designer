# Positional motif enrichment in functional saRNAs

**Question.** Chapter 4 of the thesis reports that functional saRNAs carry particular bases at
particular positions (G at 10, U at 15, C depleted at 9) with differences of 20-25 percentage
points against non-functional designs. Should those positions become design rules?

**Answer.** No. They reproduce exactly, and none of them survives correction for the size of the
grid they were found in. They are carried as *optional, off-by-default* rules so a user can test
them deliberately (ADR 0011), not as criteria that shape the default ranking.

Measured by `scripts/positional_motif_significance.py`. Numbers below are reproducible with the
seed the script prints.

## What reproduces

Every value quoted in the chapter recomputes from the published appendix sequences:

| Chapter's claim | Recomputed |
|---|---|
| pos 10 G, +24.8 pp | +24.8 pp |
| pos 15 U, +23.0 pp | +23.0 pp |
| pos 15 G, −22.6 pp | −22.6 pp |
| pos 16 G, +22.1 pp | +22.1 pp |
| pos 9 C, −21.7 pp | −21.7 pp |
| pos 9 G, +20.3 pp | +20.3 pp |
| UG at pos 15, +24.1 pp | +24.1 pp |
| GG at pos 12, +20.5 pp | +20.5 pp |
| CA at pos 14, +16.0 pp | +16.0 pp |
| AA at pos 17, +18.6 pp | +18.6 pp |

The transcription of the appendix into
`src/sarna/validation/data/vondervoort_2026_chapter4.csv` is therefore sound, which is the first
reason to run this analysis at all.

**One correction.** The chapter's discussion lists "AC (pos 9)" among the strongest enriched
dinucleotides at +19.6 pp. The motif at +19.6 pp is **AG** at position 9; AC is +0.8 pp. AG is
also what the single-nucleotide result predicts, since position 9 shows G enriched and C depleted.
This looks like a typo rather than a different calculation.

## Why the enrichments do not survive

Two properties of the measurement, neither visible in a heatmap:

**The grid is large.** 19 positions × 4 bases is 76 single-nucleotide cells; dinucleotides give
270 observed cells and trinucleotides 772. Reporting the largest difference found across 76 cells
is reporting a maximum. A maximum's null distribution sits far to the right of any single cell's,
so the usual intuition about what a 20-point difference means does not apply.

**One sequence is worth 7.1 points.** With 14 activators, a frequency can only be 0, 7.1, 14.3,
21.4 … percent. Every "hotspot" in the table above is three sequences or fewer. The second decimal
place in "+24.8 pp" describes which of 14 sequences fell where, not a property of saRNA biology.

A max-statistic permutation test addresses both at once: shuffle the 125 labels keeping group
sizes at 14 and 111, recompute every cell, keep the largest absolute difference, repeat 20 000
times. No independence assumption is needed, which matters because neighbouring positions and
overlapping k-mers are strongly dependent.

| Motif width | Cells | Observed max | Null median | Null 95th pct | Family-wise p |
|---|---|---|---|---|---|
| 1-mer | 76 | 24.8 pp | 31.1 pp | 41.8 pp | 0.939 |
| 2-mer | 270 | 24.1 pp | 25.8 pp | 35.6 pp | 0.673 |
| 3-mer | 772 | 19.6 pp | 20.5 pp | 28.6 pp | 0.734 |

The result is stronger than "not significant". At every motif width the observed maximum is
**below the median of the null**: a random 14-against-111 split of these same 125 sequences
typically produces a *larger* positional difference than the real activity labels do. Whatever
distinguishes the 14 activators, this analysis cannot see it in positional base composition.

This agrees with the chapter's own PCA, which found functional and non-functional sequences
overlapping with the first two components explaining 5.8% and 5.2% of variance. Both analyses say
the same thing; only the motif heatmaps look like they say otherwise, because a heatmap has no
way to show the null.

## A strand correction that changes the interpretation

The appendix labels its sequence column only "Sequence". The chapter's motif analysis describes
the positions as lying "along the 19-nucleotide guide strand". The sequences are **target/sense**
strands.

Established from the data rather than assumed: 113 of 125 carry A at position 19 and 123 of 125
carry A or U at position 18, which are Wang Table 1 criteria written about the *target* site. Read
as guide strands the same criteria would appear reverse-complemented, requiring U at position 1,
which only 15 of 125 have. The reverse complement of the stored P21-322 entry is
`UACUUGGAGAAUGAGUUGG`, the published dsP21-322 antisense strand.

This matters mechanistically, not just clerically. **Position 15 of a target is position 5 of its
guide.** Counting on the wrong strand puts a finding at the wrong end of the wrong molecule, and
the chapter's interpretations are anchored accordingly: "weak bases near the 3′ end" of what it
calls the guide are in fact near the 3′ end of the *target*, which is the **5′ end of the guide**.

Re-anchored, the observation becomes more interesting rather than less. A/U-rich sequence at the
guide's 5′ end is exactly what asymmetric RISC loading predicts: a thermodynamically weak 5′ end
is what makes a strand the guide (Schwarz et al. 2003; Khvorova et al. 2003), and AGO2 prefers a
5′ U (Ui-Tei et al. 2004). So the chapter's stability gradient (strong bases centrally, weak
bases toward target position 15-17) is a restatement, on the target strand, of established
strand-selection thermodynamics the engine already models explicitly in `sarna.thermo.guide`.

That is the useful residue of this analysis: not a new positional rule, but corroboration that
the asymmetry criterion already in the default preset is pointing at something real. It also
means a "position 15 U" rule would partly double-count the asymmetry rule rather than add
information.

## What would make this answerable

The limitation is the 14 activators, not the method. Rough guide, assuming the effect sizes are
real as reported: separating a 25-point difference from this null at 80% power needs on the order
of 50-60 activators, which at the field's 15-20% hit rate means screening 300-400 designs. Routes
that need fewer:

- **Test the specific hypotheses prospectively.** Pre-register "G at position 10" as one rule,
  design duplexes that differ at position 10 and nothing else, and the grid collapses from 76
  cells to 1. The same data would then be decisive.
- **Pool with published screens.** The constraint is a shared activity threshold and comparable
  assay conditions, which is why this dataset's own 14 positives are already heterogeneous.
- **Model the locus, not the sequence.** Chapter 4's largest effect by far is that AF22 produced
  no activators at any locus, including the P21 control. A sequence-only model cannot represent
  the variable that mattered most.

## Bearing on the design engine

- The positional preferences are catalogued as `EXPLORATORY`, off by default, with this document
  as their provenance (ADR 0011). Selecting one states a hypothesis; it does not apply a finding.
- The default `wang-2015` preset is unchanged. Nothing here justifies altering it.
- `pos7-u`, already carried as unsourced and off by default, is in the same category. This
  analysis gives no support for it either: position 7 appears nowhere in the top differences.
- The dataset's own caveats (`CHAPTER4_CAVEATS`) travel with any export that uses it, so a
  downstream reader is told what the set cannot show.
