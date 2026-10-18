# ADR 0023: Wang's own homology threshold, screening in batches, and the on-target hit

**Status:** accepted
**Date:** 2026-10-11
**Deciders:** Joep Vondervoort

Extends [ADR 0019](0019-off-target-screening.md), which decided that screening is opt-in,
remote, dual-scope and reported rather than filtering. None of that changes here. What changes
is *how much* is searched at a time, *which* alignments are called strong, and what happens to
the alignment every genomic candidate has against its own locus.

## Context

The first real screening of a 22 candidate run produced a table in which one candidate carried
500 alignments and the other 21 carried zero. None of those zeros was a measurement. The parser
held a single budget of 500 alignments for the whole multi-query job and returned as soon as it
filled, so the first query consumed it and the remaining 21 were never read. Reproduced against
a synthetic three query report of 600, 4 and 3 alignments: 500 parsed, all attributed to the
first query, the other two reported clean.

That defect is fixed separately and needs no decision. Reading the output it produced raised
three questions that do.

Going back to the source settled part of it. **Wang J, Place RF, Portnoy V, Huang V, Kang MR,
Kosaka M, Ho MKC, Li L-C. *Inducing gene expression by targeting promoter sequences using small
activating RNAs.* J Biol Methods, 2015 (PMC4379447)** says more about off-targets than this
project had carried across:

- It states exactly one homology threshold, **80%**, used to reject a control duplex showing
  "significant homology" to anything.
- Its own candidate search, steps 7 and 8, prescribes **megablast** against reference genome
  sequences. At 19 nt megablast returns nothing, so this project already sets `MEGABLAST=off`
  and reproduces `blastn-short` explicitly. That divergence was a code comment and is now
  recorded below.
- It uses the word "seed" once, in quotation marks, in the sentence "Each activating duplex
  likely possesses a different 'seed' sequence and target non-overlapping regions." It defines
  no seed region, gives no positions, and makes no design criterion of it.
- Its step 9, a BLAT check rejecting any target containing a SNP, is **not implemented in this
  project at all**. Recorded here as a known gap, not decided by this ADR.

## Decisions

### Two identity columns, one of them sourced

`near_perfect_hits` at 17 of 19 identities was chosen locally and its own docstring says so. It
is kept, because exports and saved runs carry it, but it is no longer what the interface leads
with. Two measures are now reported side by side:

| Column | Threshold | Where it comes from |
|---|---|---|
| Exact | the whole core, no gaps, no mismatches | the only unambiguous case |
| Wang 80% | `ceil(0.80 x core length)`, so 16 of 19 | Wang 2015's control duplex rejection rule |

**The cost, stated.** Wang applies 80% to *control* duplexes, where the question is whether a
scrambled sequence is clean enough to serve as a negative control. Applying it to candidate
off-targets is an extension by analogy, not something the paper asks for. It is also looser
than 17, so more candidates are flagged and some of those flags will be noise. What it buys is
that the number a researcher reads first has a citation behind it, which 17 never had. The
label in the interface names the paper so the analogy is visible rather than implied.

Rejected: replacing 17 with 16 outright. Saved runs and the xlsx export already carry the 17
column, and silently changing what a stored column means is worse than carrying two.

### Five matches at 80% is called promiscuous, as a display rule

A candidate with five or more alignments at Wang's threshold is summarised as binding in many
places and its alignment list is collapsed behind a count rather than drawn in full.

Five is a project-selected count, without literature validation. It is therefore declared
as what it is: a **display** threshold that changes no rank, excludes no candidate, and alters
no stored measure. The underlying alignments stay in the record and the list expands on demand.

**The cost.** It cannot make a search finish sooner. One BLAST job runs every query regardless
of what the client does with the answer, so the only real saving in wall clock time comes from
screening fewer candidates, which is the next decision. What collapsing buys is a readable page,
a smaller payload and skipped transcript placement calls.

### Screening runs in researcher-chosen batches

The fixed top 50 becomes a choice: 1, 5, 10, 15, 25 or every passing candidate, with 15 the
default, plus "the next 15" and "all the rest" once some are done, and "screen my ticked
selection" reusing the table's existing checkboxes.

Mechanically all of these resolve to one thing. The browser sends an explicit list of candidate
identifiers and the server validates that each belongs to the run and passed its gates. A
request for more than one job's worth runs as several sequential jobs, one at a time, each its
own dated record exactly as ADR 0019 requires, never as parallel submissions to NCBI.

**The cost.** A run's screening is now several records per scope instead of one, so the
interface has to merge them and show which search each candidate's row came from. That is more
state than one record and one table. It is accepted because the previous behaviour, one fixed
cap of 50 with no way to ask for fewer, is what made a genomic search time out after 1800
seconds and return nothing at all.

Unchanged: one job per request, never one per candidate, and never faster than NCBI's requested
15 second poll interval.

### The genomic on-target hit is removed from the findings, and counted

Every candidate matches its own promoter locus perfectly, so an unhandled genomic screening
reports at least one perfect match for every candidate and the column a researcher reads first
means nothing.

The hit is identified by **exact coordinate overlap** with the candidate's own genomic interval,
which the run already records. No threshold is invented and no sequence is guessed at.

It is removed from `hits`, from every count and from the table. It is **not** deleted: it moves
to `on_target_hits` on the record, and the count of what was moved is stated beside the results.
A zero that had something taken out of it must not read like a zero that never did.

Two honest limits, both reported rather than papered over:

- A run designed from pasted sequence has no genomic coordinates, so nothing can be identified.
  Those screenings say so instead of guessing.
- The transcript scope gets no removal at all. A promoter target is usually not transcribed, so
  there is normally nothing to remove, and matching by gene name in a BLAST description is a
  guess this project will not make silently.

**The cost.** A reader of the raw record sees two hit lists where there was one, and a reader of
the table sees a number that has had something removed from it. Both are stated in the record
and on the page.

### No seed region is drawn

An alignment's footprint is drawn against the core, showing which bases matched and where the
mismatches fall. Positions 2 to 8 are **not** marked as a seed.

The convention comes from siRNA off-target tooling, not from Wang 2015, which names no seed
region. Marking one would put a borrowed mechanism on the screen with this tool's authority
behind it, and [ADR 0011](0011-optional-rules-are-hypotheses-not-findings.md) already refuses
that pattern for rules.

**The cost.** Seed position is genuinely how a reader judges whether a 16 of 19 match matters,
and withholding the annotation leaves that judgement entirely to them. The footprint shows the
positions, so a researcher who holds that model can apply it; the tool does not apply it for
them. Revisit with a citation that covers saRNA rather than siRNA.

## Divergences from Wang 2015 recorded, not resolved

| Wang 2015 | This project | Why |
|---|---|---|
| Candidate search uses megablast | `MEGABLAST=off`, `blastn-short` reproduced explicitly | megablast returns nothing at 19 nt, which is indistinguishable from a clean result |
| Control duplexes rejected at 80% homology | 80% reported as a column, nothing rejected | ADR 0019: screening reports, it never filters |
| Step 9 rejects targets containing a SNP (BLAT) | not implemented | a real gap; needs a SNP source and its own decision |
| "Expression of putative off-target genes should be screened" | stated in the interface as the paper's own caveat | it is an experimental instruction, not something software can do |

## Consequences

- `CandidateSpecificity` gains `wang_homology_hits` and the threshold it was measured at.
  `near_perfect_hits` keeps its meaning and its tests.
- `ScreeningRecord` gains `on_target_hits` and the explicit candidate list it searched.
- Alignment parsing is capped per candidate, and a candidate whose list was capped says so, so
  "no alignments", "some alignments" and "more than we kept" stay three different statements.
- The interface shows where each alignment lands: the gene, the matched span against the core,
  the strand, and for genomic hits the contig, position and signed distance to the nearest
  annotated start site.
