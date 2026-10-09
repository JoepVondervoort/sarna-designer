"""Transcript-level expression: which transcripts of a gene are actually used, and where.

The workbench offers transcripts in a preferred order (MANE Select first (ADR 0012)) but
"the reference annotation's representative transcript" and "the transcript this tissue
actually makes" are different claims, and for this project's own flagship gene they disagree.
In GTEx v10, SCN1A's MANE Select transcript ``ENST00000674923`` has a median of **0.00 TPM in
all 54 tissues**, including all 13 brain regions, while ``ENST00000303395`` carries 1.63 TPM in
cortex and ``ENST00000637038`` 1.54. A researcher designing an saRNA for SCN1A in neurons is
choosing a promoter, and which transcript's promoter that is turns out to matter.

This subpackage adds that evidence beside the choice without making the choice. It never
reorders transcripts, never scores a candidate and never filters: ordering stays annotation's
job, and expression is shown as what it is, a median over post-mortem bulk tissue.

**It is almost all tissue, and none of it is the cells this project uses.** GTEx v10 has two
cultured sample groups, donor-derived fibroblasts and EBV-transformed lymphocytes, and 52
post-mortem tissue groups. It has no AF22, no HeLa and no iPSC-derived neuron. "No GTEx value
for AF22" is a statement about GTEx, not about AF22, and :mod:`sarna.expression.models` keeps
the three ways a number can be missing apart so an interface cannot render ignorance as a
measurement.
"""
