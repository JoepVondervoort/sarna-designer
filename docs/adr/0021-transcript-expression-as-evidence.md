# ADR 0021: Transcript expression is shown beside the reference choice, and never reorders it

**Status:** accepted
**Date:** 2026-10-08

## Context

[ADR 0012](0012-mane-select-or-refuse.md) settled *which* transcript this tool offers first: MANE
Select, or it refuses to choose. That ADR is about reaching a defensible reference. It says
nothing about whether the reference it reaches is a transcript the researcher's cells actually
make, and the user asked for exactly that:

> it should show the "main" transcript on top and it should also show as much as possible which
> transcripts are transcripted in which cells. So for a researcher it is nice to know if certain
> transcripts are transcript in the brain, or af22 cells for example.

The two questions can disagree, and for this project's own flagship gene they do. Measured
against GTEx v10 on 8 October 2026:

| SCN1A transcript | Brain, cortex | Highest of any of the 54 tissues |
|---|---|---|
| `ENST00000674923` (**MANE Select**, Ensembl canonical) | **0.00 TPM** | **0.00 TPM**, in every tissue |
| `ENST00000637038` (nonsense mediated decay) | 1.54 TPM | 3.02 TPM, frontal cortex BA9 |
| `ENST00000303395` (protein coding) | 1.63 TPM | 2.49 TPM, frontal cortex BA9 |

The transcript this tool puts at the top of the picker has a median of zero in all 54 GTEx
tissues, including all 13 nervous-system groups. The MANE Select designation was independently
confirmed from Ensembl's own record for the gene, so this is not an identifier mix-up.

That is not an argument for dropping ADR 0012. A median of zero in post-mortem bulk tissue is
weak evidence about a transcript, and the second-highest brain value belongs to a transcript
annotated `nonsense_mediated_decay`, whose promoter is a questionable design target. It is an
argument for the researcher seeing both facts at the moment they choose.

Candidate sources considered:

- **GTEx, transcript level.** 54 post-mortem tissue groups, 13 of them nervous system, with
  per-transcript medians. No cell lines at all. **Chosen by the user.**
- **GTEx, gene level.** Cannot answer the question, which is about transcripts.
- **FANTOM5 CAGE.** Measures transcription start usage directly, which is closer to what a
  promoter-targeting tool wants, and covers some cell lines. Harder to join to transcript
  identifiers and not yet built.
- **DepMap / Expression Atlas.** Have cell lines, including HeLa. Still no AF22.
- **Nothing.** Leaves the picker as a list of identifiers with no indication which matter.

## Decision

**Expression is evidence shown beside the reference choice. It never reorders the picker, never
scores a candidate and never filters one.**

1. **The picker's order stays the annotation's.** MANE Select first, then Ensembl canonical, then
   protein coding, then by identifier. Expression annotates the options in place. The browser
   check asserts that the order is byte-identical before and after the evidence loads, because a
   list sorted by expression is a recommendation, and this tool is not in a position to make one
   from bulk tissue medians.
2. **It is a separate request from gene lookup.** Lookup is on the path to designing anything and
   GTEx is not. Folding them together would make every lookup wait on a second service and fail
   when that service is down. A failed retrieval renders as a failed retrieval, carrying no
   transcript records at all, because a column of blanks reads as measured absence.
3. **Three silences stay distinct**, mirroring `CellTypeStatus`: measured and zero; absent from
   the annotation release the dataset was built on; and not covered by the resource at all. The
   first is a finding, the second is a gap, and the third says nothing about the sample.
4. **The dataset names its own annotation release.** GTEx keys on versioned GENCODE identifiers,
   and `gtex_v8` is GENCODE v26 while `gtex_v10` is v39. Querying with the wrong version returns
   an empty result and HTTP 200, which reads exactly like "not expressed". Every lookup resolves
   the versioned identifier through `reference/gene` first and refuses a reply about another
   release.
5. **A gene GTEx has never seen raises** rather than returning a profile full of zeros, for the
   same reason.
6. **`gtex_v10` is the default**, not v8, because v8's GENCODE v26 predates the MANE Select
   transcripts this app offers first and cannot quantify them at all.

## Consequences

- A researcher choosing a promoter for SCN1A now sees, at the moment of choosing, that the
  default reference is not detectably transcribed in any GTEx tissue. They may still choose it.
  Annotation and expression are different kinds of evidence and the tool presents both rather
  than resolving the conflict.
- **AF22 is still not answered, and the panel says so in those words.** GTEx has no cell lines.
  Adding this source makes the gap more visible rather than less, which is the intended effect:
  the tool will keep knowing more about bulk tissue than about the disease-relevant line until a
  source that covers it exists.
- The annotations disagree in both directions. Six SCN1A transcripts Ensembl serves today are
  absent from GENCODE v39, and one that v39 quantifies, `ENST00000375405`, has since been retired
  from Ensembl. Both are reported rather than dropped.
- Medians carry their donor counts, which differ **74-fold** between GTEx groups: 11 donors for
  kidney medulla against 818 for skeletal muscle. Two numbers printed to the same precision are
  not equally certain, and the panel says so.
- GTEx files cervical spinal cord among its 13 "Brain" groups. The panel reports the highest
  single group and states that it is not an average, because averaging regions nobody pooled
  would invent a measurement.
- A second external service is now reachable from the app. It is cached, offline-aware and
  refuses rather than returning nothing, on the same terms as Ensembl, and
  `sarna.sources.json_http` is that behaviour written once.
- **Not decided here:** whether expression should ever influence ranking. It currently cannot,
  and nothing in this ADR opens that door. Doing so would need a source measuring the cells being
  transfected, which is the gap this ADR documents rather than closes.
