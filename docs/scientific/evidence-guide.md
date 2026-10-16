# Reading the evidence in saRNA Designer

A high checklist score is not an activation probability. The application implements design
criteria and records their consequences. It has not demonstrated that its ranking predicts
activation in an independent experiment. Passing available gates can coexist with unchecked
gates; unknown evidence must remain visible.

## Published guidance

[Wang et al. 2015](https://doi.org/10.14440/jbm.2015.39), Table 1 and the target-selection
protocol, provide the default saRNA design recommendations. Their publication supports attribution
of those recommendations. It does not make every threshold a universal biological necessity.
The bundled CEBPA example fails position 19 A, and recovered project coordinates place known
activators outside the default upstream window. These are counterexamples within the recorded
evidence, not estimates of the benefit of removing a gate.

The default remains the Wang preset. Changing a gate to a preference is a recorded researcher
choice. Gate, preference and off are software roles, not evidence grades.

## Evidence transferred from another setting

[Reynolds et al. 2004](https://doi.org/10.1038/nbt936) and
[Ui-Tei et al. 2004](https://doi.org/10.1093/nar/gkh247) investigated siRNA silencing.
Optional criteria attributed to them extrapolate to promoter-targeted activation. Support for
knockdown does not establish activation efficacy. Sharing a proposed loading mechanism is a
reason to investigate transfer, not confirmation that transfer succeeds.

[Schwarz et al. 2003](https://doi.org/10.1016/S0092-8674(03)00759-1) studied asymmetric
assembly of the RNAi complex. [Xia et al. 1998](https://doi.org/10.1021/bi9809425) measured
RNA duplex thermodynamics. Here terminal-energy calculations support a model prediction of strand
preference. They do not measure loading in the experimental cell or establish activation.

## Project choices and unsupported hypotheses

Equal preference weighting, thermodynamic tie-breaking, and their combination into the displayed
order are project choices. Wang does not supply that ranking formula. Counting four flank
positions separately also gives that group more opportunities to contribute; equal points do
not imply equal or independent biological evidence.

Chapter 4 positional hypotheses were derived from the same dataset used in the current benchmark.
Their apparent performance on that dataset is not independent validation. The unsourced position-7
rule remains explicitly unsupported. Neither a source citation nor a manually chosen threshold
turns a hypothesis into a validated predictor.

The historical `established`, `empirical`, and `exploratory` identifiers are retained in saved
records for compatibility. Read their accompanying evidence-scope explanation: `established`
refers to support in the original scientific setting, not established saRNA efficacy.

## What sequence annotation cannot establish

[Takai and Jones 2002](https://doi.org/10.1073/pnas.052410099) provide sequence criteria for
CpG islands. Using those criteria in this application's sliding-window implementation is a
project choice, not evidence of methylation or accessibility. A short sequence or insufficient
flanks cannot exclude an island extending beyond the supplied region. Such negative checks
report unknown; a detected overlap still reports the observed sequence feature.

An annotated transcript start, a predicted TF motif, tissue expression, and a histone signal
are different observations. None alone establishes cell-specific promoter activity or saRNA
response. BLAST results describe returned sequence alignments with recorded parameters. No hits
is not proof of biological specificity, and an alignment is not measured off-target activity.

## What the bundled benchmark measures

See [empirical validation](empirical-validation.md) and the
[dataset definition](../../src/sarna/validation/dataset.py). The 125 labelled sequences include
14 positives drawn from heterogeneous studies and assay contexts. Some project positives were
reproducible without reaching statistical significance. Negative labels apply to the tested
conditions, not all possible conditions.

The app's displayed AUC measures sequence-only preference points. It does not include recovered
coordinates or the production thermodynamic tiebreak. No preset demonstrates discrimination above
chance in that analysis. This is absence of demonstrated predictive benefit, not proof of no
benefit. Range restriction, small positive counts and heterogeneous conditions limit inference.
Independent evaluation is needed before claiming predictive improvement.

## Correctness changes and existing records

New runs record per-rule evidence scope and curated source links in JSON and HTML; Excel carries
those fields on its Rules sheet. Old snapshots remain unchanged. Missing historical evidence
metadata is labelled as not recorded, rather than filled from today's catalogue.

Screening exports use the latest completed finding per candidate and scope, preserving other
batches. Failed or incomplete searches do not replace completed evidence.

The generic JSON cache now hashes endpoint plus full request path. Old filenames are retained on
disk but are not reused because their request identity cannot be verified. Cache-only users may
need an online retrieval to populate the new entries. Saved run snapshots are unaffected.
