# Transcription factor binding: what the tool can honestly show

**The short version.** The tool can show where a factor *could* bind (motif match, from sequence),
and where a region is annotated regulatory *in some cell type* (Ensembl Regulatory Build). It cannot
currently show where a factor *is* bound in HeLa or AF22, because no ChIP-seq data for those lines
is wired in. The three are different claims and the model keeps them apart by construction.

## Three claims, not one

`sarna.regulatory.evidence.EvidenceKind` makes the distinction mandatory on every annotation:

| Kind | Claim | Cell type | Abundance |
|---|---|---|---|
| `motif-prediction` | the sequence could bind this factor | **none, by construction** | very high |
| `aggregated-annotation` | this region is regulatory in some cell type | pooled, unspecified | moderate |
| `cell-type-observation` | this factor was observed bound here, in this cell line | named | rare |

One rule is enforced in code: **a motif prediction may not name a cell type.** Constructing an
`Evidence` with `MOTIF_PREDICTION` and a cell type raises, because it would imply a measurement
nobody made. The same applies to aggregated annotations.

A fourth state matters as much as the three claims. `CellTypeStatus` separates `OBSERVED_BOUND`,
`ASSAYED_NOT_BOUND` and `NOT_ASSAYED`. Most factor-by-cell-type combinations are in the third state.
An interface that renders "nobody looked" the same as "not bound" converts ignorance into a negative
finding, which is the specific failure this model exists to prevent.

## Motif scanning, and why it over-reports

Sixteen JASPAR CORE matrices ship, each with a stated reason for being in the set. SP1 because
Li et al. 2006 name SP1 sites in the p21 promoter beside the founding saRNA target; REST because it
represses SCN1A in non-neuronal cells and is the best single candidate explanation for a
HeLa-versus-AF22 difference at that locus; SOX2 and PAX6 because they define the AF22 state; and so
on. Scanning all ~800 CORE matrices would produce several hundred hits over a 2 kb promoter and
nothing a researcher could act on.

Two corrections to the conventional method are applied, both because the naive version misleads.

### The relative-score threshold admits matches less likely than background

JASPAR's convention is a **relative score** (the raw log-odds rescaled onto [0, 1] over the
matrix's achievable range) with a default cutoff of 0.80. That makes thresholds comparable between
matrices, but it depends on how negative each matrix's *minimum* is, and a matrix with one near-zero
count has a very negative minimum. So 0.80 of the range can still be a **negative log-odds**,
meaning the match is less likely under the motif than under uniform background.

Measured over the sixteen curated matrices, four are in that state at relative score 0.80:

| Motif | Raw log-odds at relative 0.80 |
|---|---|
| NRF1 | −10.2 |
| MEF2C | −1.7 |
| YY1 | −1.4 |
| ETS1 | −1.1 |

Before the fix, NRF1 matched **141 times in 408 bp** of GC-rich sequence. The module now requires
the raw log-odds to exceed zero in addition to the relative threshold. That is the weakest claim
worth drawing on a diagram: the sequence is at least more consistent with the motif than with
background. Total hits over that test sequence fell from over 230 to 95, and NRF1's 141 went to 0.

**The rigorous alternative is not implemented.** A per-matrix p-value from the exact score
distribution (as TFMPvalue computes) would replace both thresholds with one interpretable number.
It is the right answer and it is more machinery than a prioritisation aid currently warrants. Until
then, the threshold is a convention and the shuffled expectation below is what makes it readable.

### A hit count means nothing without the count chance would give

Every scan reports `expected_hits`: the mean number of matches over shuffles of the same sequence,
preserving base composition. Promoter sequence is GC-rich and so are the motifs of most
promoter-binding factors, so an impressive hit count is often entirely explained by composition.
`MotifScan.summary` says so in words when the observed count is at or below the expectation.

Shuffling destroys dinucleotide structure, which slightly overstates the expectation for
CpG-dependent motifs. Stated rather than corrected: a dinucleotide-preserving shuffle is more
machinery than this estimate warrants.

## Ensembl Regulatory Build, and reading absence correctly

The Build classifies segments as promoter, enhancer, CTCF binding site, open chromatin and TF
binding site, from epigenomic data **pooled across a panel of cell types**. Two consequences, both
enforced in `sarna.regulatory.features`:

- A feature present means "regulatory in at least one surveyed cell type". It does not say which,
  so the evidence kind forbids naming one.
- A feature absent means "not detected in that panel". The panel includes neither AF22 nor any
  neural progenitor. `RegulatoryAnnotation.no_features_note` returns the caveat as a *value* when
  the region is empty, so an interface has to handle the empty case rather than drawing a clean
  promoter.

One implementation note worth recording: Ensembl reports the kind in each record's `description`
field. The `feature_type` field is always the literal string `regulatory`, so reading it instead,
the obvious thing to do, makes every feature "other".

Measured consequence, from `docs/scientific/scn1a-design-coordinates.md`: the Build annotates a
promoter, an enhancer and a CTCF site at SCN1A's MANE Select TSS, and **nothing** at either of the
two TSSs the thesis designed against. Suggestive, and not evidence that those are not promoters.

## What is missing, and what it would take

**Cell-type-specific binding.** The `CELL_TYPE_OBSERVATION` kind exists and nothing produces one
yet. Two routes:

- **ENCODE ChIP-seq for HeLa** is abundant. HeLa-S3 is one of the most-assayed lines in existence,
  with hundreds of factors. Wiring in ENCODE's API would populate the HeLa column properly, and
  HeLa is the cell line in which chapter 4's activators actually worked.
- **AF22 has no public ChIP-seq.** The project already holds AF22 histone bigWigs (GSE108506), which
  give chromatin state rather than factor occupancy. H3K4me3 and H3K27ac at a promoter is evidence
  it is active in AF22, which is arguably more useful for explaining chapter 4's result than any
  single factor would be. Blocked on the GRCh37/GRCh38 reconciliation (defect B15, ADR 0004).

The asymmetry is worth stating plainly: the tool will be able to say much more about HeLa than about
AF22, and AF22 is the disease-relevant model. That is a property of the available data, not of the
tool, and the interface should say so rather than letting an empty AF22 column read as a negative
result.

## Licensing constraint

JASPAR matrices are distributed under **CC BY-NC-SA 4.0**, non-commercial, share-alike. Bundling
them in `src/sarna/regulatory/data/` therefore constrains this package's own licence, which
[`LICENSE`](../../LICENSE) records as undecided. For an internal research tool this is no obstacle.
It would become one if the tool were ever distributed commercially, and the alternative is to fetch
matrices at runtime rather than ship them. Flagged here so the decision is made deliberately rather
than discovered later.
