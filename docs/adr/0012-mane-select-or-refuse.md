# ADR 0012: The transcript is MANE Select, or the tool refuses to choose

**Status:** accepted
**Date:** 2026-10-06

## Context

A promoter is a property of a **transcript**, not of a gene. Workflow A asks a user for a gene
symbol, so something has to bridge the gap, and the bridge is not neutral.

Measured on the project's own primary target: `SCN1A` has 30 transcripts, 12 of them
protein-coding, with **16 distinct transcription start sites** and five distinct ones among the
coding transcripts. They span chr2:166,073,669 to chr2:166,149,213 — over 75 kb. Choosing between
them moves the promoter window by up to that distance, which is to say the choice decides which
sequence gets designed against, entirely.

The legacy pipeline took the first transcript Ensembl returned. The user's brief lists "never
assume one TSS" among the silent assumptions this project exists to remove.

Candidate policies:

- **First returned.** Depends on Ensembl's ordering, which is not a scientific property.
- **Longest transcript.** A real convention, but length has no bearing on which promoter is active
  in the cell type being transfected.
- **Most upstream TSS.** Deterministic, and wrong about as often as it is right: for SCN1A it
  picks a TSS 21 kb from the MANE Select one.
- **Ensembl canonical.** Reasonable, and usually the same as MANE Select, but Ensembl-internal.
- **MANE Select.** One transcript per protein-coding gene, agreed between Ensembl/GENCODE and
  NCBI RefSeq, with a matched RefSeq accession.
- **Ask the user, always.** Maximally honest, maximally obstructive.

## Decision

**`select_transcript` returns the MANE Select transcript. Where there is none, it raises, carrying
the alternatives.**

MANE Select is chosen for a property the others lack: it is an *agreement between two independent
annotation projects*, so it is both reproducible and externally justifiable to a reviewer. The
matched RefSeq accession is recorded on the transcript (`Transcript.refseq_match`) because the
accession is what lets a reader check the agreement rather than take it on trust. Thesis chapter 4
used MANE Select for its own designs, so promoters fetched here are directly comparable with that
work.

The refusal matters as much as the choice. `NoManeSelectTranscriptError` is **a question, not a
failure**: many genes legitimately have no MANE Select transcript, and for those the promoter
depends on a judgement the tool must not make. The error lists the protein-coding alternatives
with their TSSs so a caller can present them, and says in as many words that the choice is the
user's. More than one MANE Select annotation also raises, rather than being resolved by taking the
first.

Two related decisions fall out:

**The default window is −2000 to +500 bp**, wider than any rule window. Wang's criteria permit
−1000 to −100, so a user sees context either side of the region the rules allow, and the window
matches chapter 4's promoter definition.

**The assembly is a required argument that is checked against the response.** Ensembl serves
GRCh38 from one host and GRCh37 from another; the client derives the host from the assembly and
raises `AssemblyMismatchResponseError` if a response disagrees. This is defect B15 made
unrepeatable.

## Consequences

**Good:**

- The promoter a run used is traceable to a named transcript and a RefSeq accession.
- A gene with no MANE Select transcript produces a question with a list, not a silently wrong
  promoter.
- Verified against the live API: the client reproduces the transcript's TSS from the promoter's own
  interval and `tss_start` by separate arithmetic, which is an independent check on the
  minus-strand case rather than a restatement of it.

**Costs:**

- **Workflow A is not fully automatic.** Some genes need a round trip to the user. That is the
  intended behaviour, and it means the UI must have a transcript-choice screen before it can claim
  to handle arbitrary genes.
- **MANE Select is not always the biologically relevant promoter.** It is the agreed
  representative, not the one active in the cell type under study. For SCN1A specifically, none of
  the thesis's own 15 SCN1A designs targets the MANE Select promoter — they cluster around two
  other TSSs (see `docs/scientific/scn1a-design-coordinates.md`). A tool that only ever offered
  MANE Select would not have been able to reproduce this project's own earlier work.
- **MANE Select covers protein-coding genes only.** A lncRNA target needs a different policy,
  undecided.

## Follow-on, not yet decided

The SCN1A finding argues for an explicit **"design against a different transcript"** affordance
rather than MANE-Select-or-nothing, and for showing *all* TSSs on the promoter view so a user can
see what they are not designing against. Both are UI decisions; the client already returns every
transcript, so no engine change is needed for either.

## Alternatives rejected

**Default to Ensembl canonical when MANE Select is absent.** Tempting, and it would make Workflow A
always succeed. Rejected because it converts the exact question this ADR exists to surface — which
promoter? — back into a silent default, and because "canonical" is an Ensembl-internal designation
a reviewer cannot check against a second source.

**Design against every TSS and merge the results.** Scientifically the most complete answer, and it
may eventually be right. Rejected for now because it multiplies the candidate set by the number of
TSSs with no way to tell the user which results to believe, and because the ranking has no term
for "how likely is this TSS to be the active one in your cell type". Revisit once the regulatory
annotation layer can speak to promoter activity per cell type.
