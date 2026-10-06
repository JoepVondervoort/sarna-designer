# AF22 chromatin at the target promoters

**Two findings, pulling in opposite directions.**

All six of chapter 4's target promoters are **open** in AF22, as is the P21 positive control —
H3K4me3 at 17–147× background, with H3K27me3 depleted at every one. So the chapter's proposed
explanation for AF22's failure, that restrictive chromatin limits promoter accessibility, is **not
supported at these loci**. The limitation lies downstream of chromatin.

SCN1A is the exception, and it is closed — no mark enriched at any of its four candidate TSSs. By
the mechanism the chapter itself states, a saRNA cannot activate a locus that is not already
transcribed, so SCN1A saRNAs screened in undifferentiated AF22 are being tested at a promoter that
is not yet active.

Produced by `scripts/af22_chromatin_at_targets.py` from GSE108506, 2026-10-06. This is also the
work that closes defect **B15**: the GRCh37/GRCh38 reconciliation ADR 0004 put behind the bigWig
adapter now exists.

## Why the controls come first

The tracks are GRCh37/hg19; the engine is GRCh38. At SCN1A's MANE Select TSS those differ by
**856,510 bp** — GRCh38 chr2:166,128,019 is hg19 chr2:166,984,529. Reading a GRCh38 coordinate
straight out of these files does not return noisy data, it returns a different gene's data, and the
legacy pipeline did exactly that.

So a null result is only interpretable beside a positive control:

| Gene | H3K4me3 | H3K27ac | H3K27me3 |
|---|---|---|---|
| GAPDH | **87.8×** | 28.4× | 0.21× |
| ACTB | **146.7×** | 60.0× | 0.49× |

Fold over the median of 250 random same-contig windows of the same width, ±500 bp of the MANE
Select TSS. Housekeeping promoters land two orders of magnitude above background, which is what
makes everything below readable.

## The target promoters are open

| Gene | H3K4me3 | H3K27ac | H3K27me3 | Verdict |
|---|---|---|---|---|
| CDKN1A (P21 control) | 41.0× | 8.2× | 0.55× | open |
| CHD2 | 26.6× | 26.2× | 0.53× | open |
| DEPDC5 | 122.0× | 18.2× | 0.10× | open |
| NPRL2 | 44.3× | 8.9× | 0.38× | open |
| NPRL3 | 17.0× | 11.1× | 0.11× | open |
| TSC1 | 146.8× | 29.9× | 0.00× | open |
| TSC2 | 69.5× | 8.8× | 0.15× | open |
| SCN1A | 1.15× | 2.08× | 0.92× | **closed** |

Every target except SCN1A carries strong H3K4me3 and H3K27ac, and **no target is
polycomb-repressed** — H3K27me3 is depleted at all eight, between 0.00× and 0.92×.

### This contradicts the chapter's stated hypothesis

Chapter 4's discussion proposes that AF22 is less permissive to RNAa

> possibly due to lower transfection efficiency, reduced availability of AGO2 or other cofactors,
> or more restrictive chromatin states that limit promoter accessibility

and elaborates that AF22 chromatin is "comparatively compact, with fewer active histone marks
(e.g. H3K4me3, H3K27ac) at many promoters".

At these promoters it is not. The third hypothesis is the one the data can speak to, and the data
goes against it. That is useful rather than merely negative: it leaves the first two.

The P21 result sharpens it. CDKN1A's promoter is open at 41× H3K4me3 in AF22, and P21 saRNA still
produced only modest, variable activation there (1.04–2.80× across chapter 4's plates, against up
to 14.6× in HeLa). A duplex failing at a demonstrably open, actively marked promoter is not failing
because the promoter is shut. The remaining candidates are **delivery** (lipofection into slowly
dividing cells), **AGO2 abundance**, **nuclear import of the loaded complex**, and the **RITA
cofactors** RHA and CTR9 that Portnoy et al. 2016 place at the promoter. Chapter 3 measured AGO2;
chapter 5 took on delivery. Neither cofactor availability nor nuclear import has been measured.

## SCN1A is closed, at every TSS its designs target

| TSS | Designs there | H3K4me3 | H3K27ac | H3K27me3 |
|---|---|---|---|---|
| MANE Select, 166,128,019 | 0 | 1.15× | 2.08× | 0.92× |
| cluster A, 166,149,213 | 5 (incl. SCN1A-10) | 0.11× | 1.54× | 1.25× |
| cluster B, 166,073,669 | 5 | 0.47× | 0.94× | 1.23× |
| cluster C, ~166,079,200 | 5 (incl. SCN1A-12, -14) | 1.20× | 0.95× | 0.81× |

Nothing is enriched anywhere. **Not a coordinate artefact**, checked three ways:

- The lifted locus lands where hg19 puts SCN1A: GRCh38 165,984,640–166,182,806 maps to hg19
  chr2:166,841,150–167,039,316, against a published hg19 SCN1A of about 166,845,670–166,984,523.
- Maximum H3K4me3 anywhere across the whole 198 kb locus: **10.4**.
- Maximum H3K4me3 in the 500 kb immediately upstream, same contig, same track: **164.3**.

The track has strong peaks next door and essentially nothing across SCN1A. The region is readable
and the gene is unmarked.

### What that means for the model

This is biologically unsurprising and methodologically important. AF22 are neuroepithelial **stem**
cells; SCN1A encodes Nav1.1 and is expressed in mature neurons, particularly
fast-spiking interneurons. A mature-neuron gene being closed in an undifferentiated progenitor is
what one would expect.

The chapter states the relevant rule itself:

> saRNAs do not "switch on" completely silent genes; rather, they enhance transcription from loci
> that already exhibit a detectable baseline level of expression.

Taken together: **screening SCN1A saRNAs in undifferentiated AF22 tests them at a promoter that is
not yet active**, so a null result there is predicted by the mechanism and is weak evidence about
the duplexes. The informative experiment is in differentiated AF22 neurons, which chapter 5 does
produce. Worth checking whether the chapter 4 SCN1A screens were in progenitors or in differentiated
cells; the methods section describes progenitor culture.

The same argument does *not* apply to the mTOR genes or CHD2, whose promoters are open in AF22. For
those, the null result is real and needs a different explanation.

## Limits of this measurement

- **No input normalisation.** These are read-depth tracks with no matched input or IgG, so
  "enrichment" is fold over a median random-window background on the same contig, not a called
  peak. The housekeeping controls at 88–147× show the measure discriminates; it should not be
  compared across datasets or treated as a quantitative occupancy.
- **The window is 1 kb.** H3K4me3 peaks sharply at the TSS, so this is the right width for that
  mark and may miss a mark offset from the TSS. A displaced H3K27ac enhancer signal would not be
  caught.
- **One dataset, one line, no replicate comparison.** GSE108506 provides single tracks per mark.
- **Absence of a mark is not silence.** A gene can be transcribed at low level without detectable
  H3K4me3, and a 5× threshold is a convention — though nothing here sits near it, with every open
  promoter above 16× and every SCN1A position below 1.3×.
- **Progenitors are not neurons.** Everything above describes AF22 in the proliferating state.

## Consequences for the tool

- The adapter is in `sarna.regulatory.chromatin` and refuses what it cannot guarantee: a
  mitochondrial query, an unmapped or ambiguous coordinate, a contig change, a strand flip, or an
  interval whose length changes by more than 2% under conversion. Each returns a stated reason, not
  a number. A signal value that might be from the wrong locus is worse than none, because it looks
  like data.
- Chromatin signal is the project's **only** cell-type-specific evidence, and it is a
  `CELL_TYPE_OBSERVATION` in the evidence model. Motif predictions remain cell-type agnostic by
  construction.
- Defect B15 is closed for this path. Query coordinates are converted per request, which is the
  third of ADR 0004's three options; the bigWigs are not rewritten and no GRCh37 endpoint is pinned.
