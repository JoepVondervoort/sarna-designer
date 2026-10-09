# How much the nearest-neighbour table changes the guide prediction

**Measured 2026-10-06.** Reproduce with `.venv/bin/python scripts/guide_prediction_divergence.py`.

## Why this was measured

[ADR 0009](../adr/0009-guide-prediction-reports-both-tables.md) decided that guide-strand
prediction reports every selected parameter set's verdict rather than one, on the grounds that the
`legacy-scanner` and `wang-2015` presets will sometimes disagree. That decision is only worth its
complexity if the disagreement is common enough to matter. This is the measurement.

## Method

Only the terminal 4 nt at each end of a target enter an end-stability sum
([`ends.py`](../../src/sarna/thermo/ends.py)), so the bases between them cannot affect the
prediction. Holding those constant and enumerating both ends covers **every input the calculation
can distinguish**: 4⁴ × 4⁴ = 65 536 combinations. This is exhaustive rather than a sample from one
promoter, so the rates below are properties of the parameter sets themselves, not of a particular
gene.

Each combination was assessed with all three shipped tables at the default `end_length` of 4 and
the legacy comparison tolerance, and the predicted guide strands compared pairwise.

## Result

| Comparison | Different guide predicted | Rate |
|---|---|---|
| RNA (Xia 1998) vs legacy script table | 8 990 / 65 536 | **13.72 %** |
| RNA (Xia 1998) vs DNA (SantaLucia 1998) | 8 538 / 65 536 | **13.03 %** |
| DNA (SantaLucia 1998) vs legacy script table | 1 050 / 65 536 | **1.60 %** |

The first disagreement found between the RNA parameters and the legacy table is
`AAAAACGGATCCAAGAAAT`, where the RNA parameters predict the **sense** strand and the legacy table
the **antisense** strand (opposite answers to the question the gate is there to decide.

## Interpretation) this corrects an emphasis in the inventory

The [legacy inventory](legacy-inventory.md) lists defect B2, the three mis-assigned stacks, as the
substantive problem with the legacy table. The measurement says the mis-assignment is the
**smaller** of two problems by roughly a factor of eight:

- **Using DNA parameters for an RNA:RNA duplex changes the predicted guide for 13.0 % of
  candidates.** This is the dominant effect. It is not a bug in the script's arithmetic; it is a
  modelling choice that was never stated. An saRNA duplex is RNA:RNA, so RNA parameters are the
  ones that describe it.
- **Defect B2 changes the predicted guide for 1.60 % of candidates**, isolated by comparing the
  legacy table against the correctly-assigned DNA parameters it was derived from. Real, and real
  enough to move candidates across a mandatory gate, but secondary.

Neither number is negligible for a mandatory gate. Together they mean that **roughly one candidate
in seven gets a different guide-strand verdict** under `legacy-scanner` than under `wang-2015`,
which is why ADR 0009 shows both rather than picking one.

## What this does not tell you

- It says nothing about which prediction is **correct**. It measures disagreement between
  parameter sets, not agreement with experiment. Validating the prediction itself needs measured
  loading data, which this project does not have.
- The rates are uniform over sequence space. Real promoters are GC-biased, so the rate on an
  actual promoter will differ. Once a promoter fixture exists, the same comparison should be run
  across its candidate set and reported beside these numbers.
- `end_length` is held at 4 throughout. A different end definition would give different rates.
