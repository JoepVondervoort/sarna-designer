# Glossary

The project's canonical vocabulary. Terms here are used consistently in code, API, UI and docs.
No implementation details, this file defines *what words mean*, nothing else.

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
identity*. The duplex is *designed* so that the antisense strand is the guide; see
*Guide prediction* for what a thermodynamic model expects for a given candidate, which is
reported alongside that intent and can disagree with it.

**Duplex**
The deliverable: a 19 bp Target/antisense core plus a 2 nt **dTdT** 3' overhang on **both**
strands, i.e. two 21-mers. The overhang is deoxythymidine, not ribo-`UU`, chapter 2's
supplementary material settles that for this project's own duplexes. It is a synthesis choice
and does not change the 19 nt core, which is what every rule is evaluated against. The overhang
is not part of the Target and is never scanned or scored.

**Complement** vs **reverse complement**
**Complement** maps each base to its partner, preserving index order. **Reverse complement**
additionally reverses the sequence. The legacy script computes a *complement* and renders it
3'→5' alongside the 5'→3' Target, which displays a correctly annealed duplex. Its internal
error message calling this a "reverse complement" is wrong. Both operations exist in this
project and must never be conflated.

---

## Genomic context

**TSS**. Transcription Start Site
The first transcribed base of a given **Transcript**, as annotated by the chosen annotation
source. A TSS belongs to a transcript, never to a gene. Written **+1**; there is no position 0.

**Promoter anchor**
What a **promoter window** is measured from. Usually a **TSS**, but not always: a promoter mapped
in a publication may have no annotated transcript at all. Chapter 2 of the thesis designed SCN1A
against promoters from Dong et al. 2014 and Long et al. 2008, and the Dong distal region, which
produced two of three SCN1A activators, has no Ensembl transcript within 5.5 kb. Every anchor
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
- **Established**, strong, direct experimental support.
- **Empirical**, a useful, reported design observation without established mechanism.
- **Exploratory**, plausible but unvalidated. Never enabled in a default preset.

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

**Rating**
What a ranked candidate reports to a reader: **preferences met out of preferences evaluated**,
itemised, plus the gates that could not be checked. The denominator is the number of preferences
that could actually be evaluated for *this* candidate, never the number configured, because a
*not-applicable* preference is not a failed one. A rating is a count of criteria, not a
probability of activation, and two configurations with different preference sets produce ratings
that are not comparable with each other.

---

## Configuration and runs

**Rule role**
Whether a rule acts as a **filter**, a **score term** or not at all, *in one configuration*.
Two roles are distinguished and both travel with a run: the **catalogue role**, which is the
role the rule declares and which no run may change, and the **effective role**, which is the
role this configuration resolved it to. A rule's effective role reaches gating, preference
counting, ranking, candidate verdicts and exports together, or it has not been changed at all.

**Off**
A rule a configuration has deliberately removed. *Off* is not *failed* and not
*not-applicable*: the rule produced no result because it was not evaluated, and the run records
that the omission was a choice. See [the evaluator contract](docs/adr/0007-evaluators-take-a-candidate-context.md).

**Custom configuration**
A configuration that departs from its base preset in any way, a parameter value, a selected
option or an effective role. A custom configuration is labelled as custom wherever it appears,
and it carries its **frozen baseline** so a reader can see what it departed from.

**Frozen baseline**
The complete resolved configuration of the base preset *as it was when the run happened*, stored
with the run. A preset definition can change between versions, so naming the preset is not
enough to reconstruct what was run.

**Departure**
One recorded difference between a custom configuration and its frozen baseline, naming the rule,
what changed and the **choice origin**: preset default, named sourced option, or researcher
custom value. A value that merely happens to equal a published option's value does not acquire
that option's provenance, origin is declared, never inferred from a number.

**Run**
One execution of the engine over one promoter with one resolved configuration, with its inputs,
evidence snapshot, candidates, verdicts and ranks. A run is the unit that gets saved, reopened,
exported and compared.

**Run snapshot**
The frozen record of a completed run. Reopening a snapshot consults no preset definition, no
remote source and no current software default; everything needed to interpret it is inside it.

**Derived run**
A run produced by changing a setting on an earlier run, linked to its parent so the difference
can be explained. Candidates are compared across runs by sequence, location and orientation,
never by rank or display id, which both move.

**Viewer** vs **workbench**
A **viewer** displays one finished run and cannot change it; the self-contained HTML page is a
viewer. A **workbench** is the local application in which a researcher chooses a gene, adjusts
settings and starts new runs. A viewer is what a workbench exports.

---

## Thermodynamics

**Nearest-neighbour (NN) model**
Estimating duplex stability by summing stacking free energies of adjacent base-pair steps.

**End stability**
The summed NN free energy over the terminal `end` nucleotides (default 4) of a Target's 5' and
3' ends. Energies are negative; **more negative = more stable**.

**Stability difference** (`stabilitydiff`)
`E(3' end) − E(5' end)` of the **Target**. A **positive** value means the Target's 5' end is
more stable, which implies the guide strand's 5' end is *less* stable, the condition that
favours guide loading. Positive is desirable. The sign convention is counter-intuitive and is
preserved from the legacy script deliberately.

**Thermodynamic asymmetry**
The principle that the duplex strand whose 5' end is less stably paired is preferentially loaded
into RISC. The basis of filters R4, R6, R7 and R8.

**Nearest-neighbour table**
One published set of stacking free energies, e.g. `rna-xia-1998`. Tables disagree, so every
thermodynamic number names the table it came from.

**Guide prediction**
Which strand of a designed duplex a nearest-neighbour table expects RISC to retain. Distinct
from design *intent*: the duplex is built so that the antisense strand is the guide, and for a
given candidate the prediction can disagree with that intent. Predictions are reported **per
table and never averaged**, because 13% of candidates flip between tables
(`docs/scientific/guide-prediction-divergence.md`). A displayed prediction always names its
table, and disagreement between tables is shown rather than resolved.

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

**"saRNA score"**, too vague; say *score contribution*, *final score* or *rank*.
**"Promoter"** alone, where a window is meant, say **promoter window** with its TSS.
**"Position"** alone, say **target position** (1–19), **offset** (TSS-relative) or
**genomic coordinate**.
**"The promoter sequence"** for a pasted input, that is an **input sequence**; it becomes a
promoter window only once a TSS is known.
**"Target strand"**, ambiguous between three different things: the DNA strand the Target is
read from, the duplex strand complementary to the Target, and the strand RISC retains. Say
**sense strand** or **antisense strand** for sequence identity, **guide strand** for strand
fate, **predicted guide** for a model's expectation, and **genomic strand** for the locus.
**"X% match"**, **"confidence"** or **"likelihood of activation"** for a rating, the rating
counts criteria met out of criteria evaluated. Nothing in this project estimates a probability
of activation, and no ranking here has been shown to beat chance.
