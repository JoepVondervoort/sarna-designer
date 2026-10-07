# Glossary

The project's canonical vocabulary. Terms here are used consistently in code, API, UI and docs.
No implementation details — this file defines *what words mean*, nothing else.

Status: in progress. Built during requirements grilling.

---

## Sequence and strand

**Target**
The 19 nt subsequence of the promoter, read 5'→3' on the promoter's **sense strand**, that a
duplex is designed against. Positions within a target are numbered **1–19 from its 5' end**, so
"position 19" is the 3'-terminal nucleotide. All position-specific rules use this numbering.

**Sense strand** / **Passenger strand**
The duplex strand whose sequence is identical to the Target. Not intended to be loaded into
RISC. The legacy script prints this as "Passenger strand". The two terms are synonyms in this
project; prefer **passenger** when discussing strand *fate*, **sense** when discussing
*sequence identity*.

**Antisense strand** / **Guide strand**
The duplex strand complementary to the Target, and the strand intended to be loaded into RISC.
Prefer **guide** when discussing strand *fate*, **antisense** when discussing *sequence
identity*. See *Guide prediction* — in this project guide fate is enforced, not predicted
per-candidate.

**Duplex**
The deliverable: a 19 bp Target/antisense core plus a 2 nt `UU` 3' overhang on **both** strands,
i.e. two 21-mers. The overhang is not part of the Target and is never scanned or scored.

**Complement** vs **reverse complement**
**Complement** maps each base to its partner, preserving index order. **Reverse complement**
additionally reverses the sequence. The legacy script computes a *complement* and renders it
3'→5' alongside the 5'→3' Target, which displays a correctly annealed duplex. Its internal
error message calling this a "reverse complement" is wrong. Both operations exist in this
project and must never be conflated.

---

## Genomic context

**TSS** — Transcription Start Site
The first transcribed base of a given **Transcript**, as annotated by the chosen annotation
source. A TSS belongs to a transcript, never to a gene. Written **+1**; there is no position 0.

**Promoter anchor**
What a **promoter window** is measured from. Usually a **TSS**, but not always: a promoter mapped
in a publication may have no annotated transcript at all. Chapter 2 of the thesis designed SCN1A
against promoters from Dong et al. 2014 and Long et al. 2008, and the Dong distal region — which
produced two of three SCN1A activators — has no Ensembl transcript within 5.5 kb. Every anchor
carries its citation and a **reference kind** saying how well its zero point is established:
`annotated-tss`, `mapped-tss`, `inferred`, or `none`. An anchor of kind `none` has no zero point,
so no offset can be computed from it and the TSS rules report *not-applicable*.

**Promoter window**
A half-open interval defined **relative to a specific promoter anchor**, e.g. `-1000 … -100`.
Negative numbers are upstream, positive downstream. A promoter window is meaningless without the
anchor it is relative to, and an offset is only as trustworthy as its anchor's reference kind.

**Assembly**
A named genome build, e.g. `GRCh38`, `GRCh37`. Every genomic coordinate in this project carries
its assembly. Coordinates from different assemblies are never compared, mixed or silently
converted.

**Offset** vs **genomic coordinate**
An **offset** is a TSS-relative position (`-347`). A **genomic coordinate** is an
assembly-anchored absolute position (`chr2:166,075,123` on `GRCh38`). Both are retained for
every candidate. Converting between them requires the TSS *and* the strand.

**Strand**
`+` or `-`, the DNA strand a transcript is transcribed from. For a `-` strand gene, increasing
TSS-relative offsets run in the direction of *decreasing* genomic coordinates. This inversion is
the single largest source of error in the legacy code and is tested explicitly.

---

## Design rules and scoring

**Design rule**
One named, independently evaluable criterion applied to a candidate, carrying its own
provenance and evidence class. A rule is either a **filter** or a **score term**, never both
implicitly.

**Filter** (mandatory rule)
A rule whose failure **rejects** the candidate. Rejection is always recorded with the rule that
caused it; rejected candidates are retained and inspectable, never discarded silently.

**Score term** (preference rule)
A rule that contributes to ranking but cannot reject a candidate.

**Evidence class**
How well-supported a rule is. Exactly one of:
- **Established** — strong, direct experimental support.
- **Empirical** — a useful, reported design observation without established mechanism.
- **Exploratory** — plausible but unvalidated. Never enabled in a default preset.

**Provenance**
For every rule: the specific source it came from (citation and locator, e.g. "Wang et al. 2015,
Table 1"), or `legacy-script` where the only source is the original implementation, or
`unsourced` where no source has been found. Shown in the UI next to the rule.

**Preset**
A complete, named, immutable set of rule settings. Presets never hide values: the fully resolved
configuration is always inspectable and exportable. See ADR 0002.

**Candidate**
One Target plus everything derived from it: offsets, genomic coordinates, duplex strands,
per-rule results, score contributions, warnings and status.

**Candidate status**
Exactly one of **generated** (window enumerated), **rejected** (failed ≥1 filter), **passing**
(all filters satisfied), **ranked** (passing, with a rank assigned). Statuses are cumulative
stages, not independent labels.

---

## Thermodynamics

**Nearest-neighbour (NN) model**
Estimating duplex stability by summing stacking free energies of adjacent base-pair steps.

**End stability**
The summed NN free energy over the terminal `end` nucleotides (default 4) of a Target's 5' and
3' ends. Energies are negative; **more negative = more stable**.

**Stability difference** (`stabilitydiff`)
`E(3' end) − E(5' end)` of the **Target**. A **positive** value means the Target's 5' end is
more stable, which implies the guide strand's 5' end is *less* stable — the condition that
favours guide loading. Positive is desirable. The sign convention is counter-intuitive and is
preserved from the legacy script deliberately.

**Thermodynamic asymmetry**
The principle that the duplex strand whose 5' end is less stably paired is preferentially loaded
into RISC. The basis of filters R4, R6, R7 and R8.

---

## Regulatory annotation

**Experimental TF binding**
TF occupancy measured on real material, e.g. a ChIP-seq peak. Carries a cell type/line and a
dataset accession.

**Predicted TF motif**
A sequence match to a TF's binding model. Carries no evidence that the TF is bound. **Never
displayed or exported as equivalent to experimental binding.**

**CpG site**
A single ordered `CG` dinucleotide.

**CpG island**
A region satisfying a stated island definition (this project uses a Takai & Jones-style
definition). A CpG island is *not* a CpG site and the two rules are independent.

**Database annotation** vs **locally calculated feature**
An annotation **retrieved** from an external source, versus one **computed** here from sequence.
Kept conceptually and visually separate even when they describe the same thing.

**Cell type / cell line**
The biological context a regulatory dataset was generated in. Every regulatory track names its
cell context and dataset. When no dataset exists for the chosen context, the application says so
rather than substituting another. Cell lines are **data, not code**: they are declared in
`cell_lines.toml` with a citation, and a user adds more through `SARNA_CELL_LINES`.

**Data axis**
One of the three kinds of evidence a cell line may have, each answering a different question:
**expression** (is the gene transcribed here at all?), **chromatin** (is the promoter accessible?)
and **TF occupancy** (which factors are bound?). Expression is asked first, because a saRNA
enhances transcription from a locus that already has some rather than switching a silent gene on.
"Not available" is a first-class answer and stays visible.

---

## Terms deliberately avoided

**"saRNA score"** — too vague; say *score contribution*, *final score* or *rank*.
**"Promoter"** alone, where a window is meant — say **promoter window** with its TSS.
**"Position"** alone — say **target position** (1–19), **offset** (TSS-relative) or
**genomic coordinate**.
**"The promoter sequence"** for a pasted input — that is an **input sequence**; it becomes a
promoter window only once a TSS is known.
