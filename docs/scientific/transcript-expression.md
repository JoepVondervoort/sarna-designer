# Transcript expression: what GTEx says about SCN1A, and what it cannot say

Measured 8 October 2026 against the GTEx Portal API v2, dataset `gtex_v10`, annotation GENCODE
v39. Reproduce with:

```bash
.venv/bin/python -c "
from sarna.expression.gtex import GtexClient
p = GtexClient().expression_for('ENSG00000144285')
e = p.for_transcript('ENST00000674923')
print(e.status, e.highest.median, e.highest.tissue_name)
"
```

## The finding

**SCN1A's MANE Select transcript is measured at 0.00 TPM in all 54 GTEx tissues**, including all
13 nervous-system groups. The transcript is `ENST00000674923`, confirmed as `MANE_Select` and
Ensembl canonical from Ensembl's own record for the gene, which is the transcript this
application offers first (ADR 0012).

The brain signal sits elsewhere:

| Transcript | Biotype | Brain, cortex | Highest brain group |
|---|---|---|---|
| `ENST00000674923` | protein coding, **MANE Select** | 0.00 | 0.00, every group |
| `ENST00000637038` | nonsense mediated decay | 1.54 | **3.02**, frontal cortex BA9 |
| `ENST00000375405` | retired from Ensembl | 0.02 | **2.68**, cerebellar hemisphere |
| `ENST00000303395` | protein coding | **1.63** | 2.49, frontal cortex BA9 |
| `ENST00000641575` | protein coding | 0.88 | 0.88, cortex |

25 of the gene's transcripts are quantified in this release. **13 of those 25 are measured at
zero in every one of the 54 groups.**

## What this does and does not mean

**It does not mean the MANE Select transcript is the wrong promoter to design against.** Three
reasons to be careful with it:

1. **Bulk tissue averages over cell types.** A transcript used by one neuronal population reads as
   near zero in whole cortex, which contains far more of everything else. GTEx cannot separate
   them.
2. **The highest brain value belongs to a decay target.** `ENST00000637038` is annotated
   `nonsense_mediated_decay`. Being abundant in a median does not make its promoter a sensible
   thing to activate.
3. **Isoform quantification is an estimate.** Transcript-level values are apportioned between
   overlapping isoforms, which is less certain than a gene-level total, and SCN1A's transcripts
   overlap heavily.

**It does mean the researcher should see it.** Reaching a defensible reference and reaching a
transcribed one are different problems, and this is a case where the answers differ. The
application shows both and resolves neither.

## The two annotations disagree in both directions

Ensembl serves 30 SCN1A transcripts today. GENCODE v39, which `gtex_v10` was built on, holds 25.

- **In Ensembl, not in GENCODE v39:** `ENST00000713692`, `ENST00000714447`, `ENST00000714448`,
  `ENST00000714449`, `ENST00000714450`, `ENST00001144292`. Nothing was ever quantified for them.
  Rendering that as 0.00 would invent six measurements.
- **In GENCODE v39, not in Ensembl today:** `ENST00000375405`. It carries a real 2.68 TPM in
  cerebellar hemisphere. Dropping it because the current annotation no longer lists it would hide
  a measurement.

`ExpressionStatus.NOT_IN_RELEASE` exists for the first case, and is kept distinct from a measured
zero for the same reason `CellTypeStatus.NOT_ASSAYED` is kept distinct from `ASSAYED_NOT_BOUND`.

## The versioning trap

GTEx keys everything on **versioned** GENCODE identifiers, and the version depends on the
dataset:

| Dataset | Annotation | SCN1A |
|---|---|---|
| `gtex_v8` | GENCODE v26 (2017) | `ENSG00000144285.17` |
| `gtex_v10` | GENCODE v39 | `ENSG00000144285.23` |

Ensembl reports a different version again. Querying `medianTranscriptExpression` with the wrong
version returns **HTTP 200 with an empty `data` list**, which is indistinguishable from "this gene
is not expressed" unless you know to look. The same is true of a gene GTEx has never heard of.

Both are now failures rather than empty results: the versioned identifier is resolved through
`reference/gene` against the release the dataset declares, and a gene with no record raises
`GeneNotInGtexError`.

`gtex_v10` is the default because **v8 cannot quantify MANE Select transcripts at all**: GENCODE
v26 predates them.

## What GTEx cannot answer

**GTEx holds nothing for the cells this project uses.** 52 of its 54 `gtex_v10` sample groups
are post-mortem donor tissue; the other two are cultured, donor-derived fibroblasts and
EBV-transformed lymphocytes. Neither is AF22, HeLa or an iPSC-derived neuron, and none of those
will appear, so "no GTEx value for AF22" is a statement about GTEx.

An earlier version of this page, and of the caveat the application shows, said GTEx "contains no
cell lines". That was wrong, and the two cultured groups are the reason. The claim that matters
for this project survives being exact, and a test now pins the corrected wording.

This is the same asymmetry recorded in [`af22-chromatin.md`](af22-chromatin.md) and in the
working notes: the tool knows more about tissues and about HeLa than about the disease-relevant
line, and adding sources tends to widen that gap rather than close it. Making the gap visible is
the point.

## Organs, and which ones get a symbol

GTEx samples at a finer grain than anyone thinks in: 13 of its groups are nervous system, 3 are
artery, 2 are heart. `sarna.expression.organs` maps all 54 onto **21 organs**, written out
rather than inferred from the names, because `Nerve_Tibial` belongs with the nervous system
without starting with `Brain`, and a prefix rule would file `Cells_Cultured_fibroblasts` under
whatever it happened to match.

Eight organs carry a symbol:

| Organ | Symbol | Why this one |
|---|---|---|
| Brain and nerve | 🧠 | anatomical |
| Heart | 🫀 | anatomical |
| Lung | 🫁 | anatomical |
| Blood | 🩸 | anatomical |
| Skeletal muscle | 💪 | the closest honest depiction |
| Thyroid | 🦋 | the butterfly is a shape convention in wide medical use |
| Kidney | 🫘 | the bean is a shape convention in wide medical use |
| Cultured cells | 🧫 | **not** an organ, and the dish says so |

The other 13 are named in words. Most organs have no Unicode emoji, and picking a lookalike
object for a liver or a pancreas would be both ambiguous and undignified in a tool that reports
measurements. An empty label is better than a wrong one.

**A symbol means: GTEx measured a median above zero in at least one of that organ's sample
groups.** It is always shown with the group and the value that produced it, so "brain" can be
seen to have come from frontal cortex rather than from an average over 13 regions that nobody
pooled. Its absence means the median was zero there, which is a measurement and not a gap. A
transcript measured at zero in every group shows no symbols at all, rather than greyed-out
ones, because a greyed symbol reads as "looked at, unclear".

For SCN1A the strip makes the gene's tissue specificity visible at a glance: every expressed
transcript is dominated by brain, with small contributions from lung and reproductive tract and
nothing elsewhere.

## Donor counts

GTEx's groups differ **74-fold** in how many donors stand behind a median. In `gtex_v10` the
smallest is kidney medulla at **11 donors** and the largest is skeletal muscle at **818**. A
median over 11 people and a median over 818 print to the same two decimal places and are not
remotely equally certain, so the count travels with every number.

GTEx also files **cervical spinal cord** among its 13 "Brain" groups. The application reports the
highest single group and states that it is not an average over them, because averaging regions
nobody pooled would be a measurement nobody made.
