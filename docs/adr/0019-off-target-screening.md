# ADR 0019 — Off-target screening: opt-in, remote, dual-scope, reported never filtered

Date: 2026-10-07. Status: accepted. Implements accepted constraint 15 in part, and **amends its
stated ordering**; see the last section.

## Context

Every order export said `off-target-screening: not-performed` ([ADR 0017](0017-order-ready-duplex-export.md)),
which was honest and useless. A researcher given 14 duplexes to synthesise had no way to ask
where else those sequences occur. This ADR records the decisions behind the answer.

## Decisions

**Opt-in and asynchronous.** The design run completes and is fully readable first; a button
starts a search that reports when it finishes. This is the user's own stated requirement and it
is also the right default for a second reason: starting a search sends this run's candidate
sequences to NCBI, and that should be a deliberate act. The queries are 19-mers drawn from the
public human genome, so no novel sequence is disclosed — but a search does tell NCBI which gene
the researcher is working on, and the interface says so beside the button.

**Remote NCBI rather than a local BLAST install.** A local install means BLAST+ binaries plus
roughly 10 GB of database, which contradicts the copy-paste installation accepted for colleagues
with no programming experience ([INSTALL.md](../../INSTALL.md)). Remote costs an internet
dependency and the privacy note above. `sarna.offtarget.ncbi.BlastClient` is a seam: a local
backend can implement the same interface later without touching anything above it.

**One job per run, never one per candidate.** NCBI limits request rates and can block an IP —
which would be the whole institution's, not one researcher's. Every candidate goes in a single
multi-FASTA submission, polling never runs faster than NCBI's requested 15 seconds, and the
client refuses to go faster. Screening is capped at the top 50 ranked candidates: a run can
generate thousands, and only the ones someone might order are worth a shared public service's
time. A truncated screening says how many it left out.

**`blastn-short` parameters, explicitly.** At 19 nt the default `blastn` word size and E-value
return *nothing*, which is indistinguishable from a clean result — the exact misreading this
project exists to prevent. `TASK` is not a documented URL API parameter, so the configuration is
reproduced explicitly: `MEGABLAST=off`, `WORD_SIZE=7`, a generous `EXPECT`, low-complexity
filtering off. Tests assert on what is *sent*, not only on what is parsed, because a
misconfigured search fails silently and in the reassuring direction.

**Both scopes, because saRNA carries two different risks.** A transcript search (`refseq_rna`)
asks whether the duplex could silence an mRNA — the classic siRNA concern, and what "BLAST it"
usually means. A genomic search asks whether it could act on *another promoter*, which is the
risk specific to what this tool designs and which a transcript search cannot see. Shipping only
the first under the label "off-target" would imply a check nobody ran.

**Genomic hits are placed, not classified.** Each one reports the distance to the nearest
annotated TSS and which gene it belongs to. **No threshold decides whether a hit "is" in a
promoter**, because this project has not defended such a definition and inventing one would turn
a convenience into a scientific claim. One flag is set: whether the hit falls inside *the same
window this run designed against*, which is a statement about the run's own configuration. The
Ensembl Regulatory Build was considered and rejected for this: `tf-binding.md` and
`af22-chromatin.md` already record that it pools cell types and misses neuronal promoters, so
its silence would be weak evidence presented as reassurance.

**Several measures, never a score.** Best E-value, near-perfect hits (≥17 of 19), longest
complementary stretch, perfect full-length matches, and for genomic searches the count inside the
run's window — each reported and sortable on its own. A single "specificity score" would have to
weigh one perfect match elsewhere against forty weak partial ones, and for some choice of weights
ranks them identically. Those are different experimental problems and the choice belongs to the
researcher. ADR 0006 already refuses weighted sums for ranking; this is the same reasoning.

**Findings report; they never filter.** No candidate is removed, reordered or marked unusable.
The summary text is tested to contain no instruction to discard anything. A researcher may
knowingly order a duplex with a perfect off-target match, and the tool's job is to make sure they
know.

**A screening is a separate, dated record.** NCBI's databases change week to week, so a result
is only meaningful with its own date, database and parameters. Writing it into the design
snapshot would both break the frozen-record rule of [ADR 0015](0015-local-workbench-and-frozen-ui.md)
and silently overwrite an earlier observation. Screenings live in their own table keyed to the
run; re-screening adds a record rather than replacing one, and both can be compared. Schema
version 4.

**A clean result is not proof of specificity**, and the summary says so in those terms: one
search, of one database, on one day, with parameters that can miss things.

## Amendment to accepted constraint 15

Constraint 15 says "build promoter-level off-target analysis first". Both scopes are built, but
**transcript screening is the one delivered as the primary button**, with genomic beside it. The
reason is that `refseq_rna` is one API call and catches the failure mode most likely to confound
an experiment, while the genomic scope depends on Ensembl annotation being reachable for every
hit. The constraint's intent — that the saRNA-specific promoter risk not be ignored — is met;
its ordering is not. Recorded here rather than quietly departed from.

## Costs and limits

Screening needs internet; the designs do not, and a failed screening leaves the run untouched.
NCBI can be slow or refuse service, and the client gives up after thirty minutes. The 50-candidate
cap is arbitrary in the sense that no measurement chose it.

**The order export still reports `not-performed`.** A screening lives in its own record and the
order CSV is generated from the run alone, so until the two are joined the export's answer stays
the honest one. That join is the obvious next slice.

`NEAR_PERFECT_IDENTITIES = 17` is a reporting threshold with no biological claim behind it. It
changes what is counted, never what is kept, and the raw alignments travel beside it.

Nothing here measures binding. Sequence similarity is sequence similarity; whether a duplex
actually engages any of these loci is unmeasured, and so is whether engagement would change
expression.
