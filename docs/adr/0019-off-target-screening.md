# ADR 0019: Off-target screening: opt-in, remote, dual-scope, reported never filtered

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
public human genome, so no novel sequence is disclosed, but a search does tell NCBI which gene
the researcher is working on, and the interface says so beside the button.

**Remote NCBI rather than a local BLAST install.** A local install means BLAST+ binaries plus
roughly 10 GB of database, which contradicts the copy-paste installation accepted for colleagues
with no programming experience ([INSTALL.md](../../INSTALL.md)). Remote costs an internet
dependency and the privacy note above. `sarna.offtarget.ncbi.BlastClient` is a seam: a local
backend can implement the same interface later without touching anything above it.

**One job per run, never one per candidate.** NCBI limits request rates and can block an IP,
which would be the whole institution's, not one researcher's. Every candidate goes in a single
multi-FASTA submission, polling never runs faster than NCBI's requested 15 seconds, and the
client refuses to go faster. Screening is capped at the top 50 ranked candidates: a run can
generate thousands, and only the ones someone might order are worth a shared public service's
time. A truncated screening says how many it left out.

**`blastn-short` parameters, explicitly.** At 19 nt the default `blastn` word size and E-value
return *nothing*, which is indistinguishable from a clean result, the exact misreading this
project exists to prevent. `TASK` is not a documented URL API parameter, so the configuration is
reproduced explicitly: `MEGABLAST=off`, `WORD_SIZE=7`, a generous `EXPECT`, low-complexity
filtering off. Tests assert on what is *sent*, not only on what is parsed, because a
misconfigured search fails silently and in the reassuring direction.

**Both scopes, because saRNA carries two different risks.** A transcript search (`refseq_rna`)
asks whether the duplex could silence an mRNA, the classic siRNA concern, and what "BLAST it"
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
run's window, each reported and sortable on its own. A single "specificity score" would have to
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
hit. The constraint's intent (that the saRNA-specific promoter risk not be ignored) is met;
its ordering is not. Recorded here rather than quietly departed from.

## What a live search changed

Two defects survived a green suite and a stub-backed browser test, and were found only by
running one real search against NCBI. Both are the kind this project keeps warning about: they
fail quietly, in the reassuring direction.

**The XML guard rejected every genuine reply.** It refused any document declaring a DTD, and
every real NCBI report opens with `<!DOCTYPE BlastOutput PUBLIC ... NCBI_BlastOutput.dtd>`. The
guard now refuses `<!ENTITY` only, which is the actual billion-laughs vector; `xml.etree` never
fetches an external DTD, so a bare DOCTYPE expands nothing. A test now parses a report carrying
NCBI's real doctype.

**The search was not restricted to human.** `refseq_rna` and `refseq_genomic` span every species
RefSeq holds. One 19-mer returned **500 alignments, 230 of them "near-perfect"**, in a sawfish,
a plant, a grass and a frog. Every one a real alignment, every one irrelevant, and together
enough to make a perfectly specific duplex look promiscuous. Searches now send
`ENTREZ_QUERY=txid9606[ORGN]`, recorded as a declared parameter so the restriction is visible in
the screening record rather than buried in a request.

**Live searches are slow.** The first completed somewhere between 7 and 15 minutes; a second, at
a stricter `EXPECT`, had not finished after 4. The bottleneck is NCBI's queue rather than the
parameters. This vindicates the asynchronous design (a researcher reads their candidates while
it runs) and the 30-minute ceiling is set accordingly. It also means screening is not something
to run casually, which is another reason it is a button rather than part of every design.

## Making it faster, 2026-10-07

Three changes, one of which was a defect rather than a tuning opportunity.

**A screening was occupying the capacity the rest of the app needs.** `AppRouter` held one
`ThreadPoolExecutor(max_workers=2)` and submitted design runs, batch children *and* screenings to
it. A design finishes in under a second; a screening waits on NCBI for seven to fifteen minutes
and may wait thirty. Starting both scopes (which is the supported use, and the one the user asked
for) and therefore occupied every worker, so the next design sat in `Queued` until a search came
back. This ADR's own premise is that the researcher reads their candidates while the search runs;
that only holds if the search is not holding the capacity the reading needs. Screenings now have
their own two-worker pool, which also caps this app at two concurrent NCBI jobs. A test submits
two blocking screenings and then a design, and asserts the design leaves the queue; it fails
against the shared pool and passes against the separate one, in 1.3 s rather than hitting its
21 s timeout.

**NCBI's own estimate was parsed and then ignored.** Every submission is answered with `RTOE`, how
long NCBI expects the search to take, and the search cannot be ready before then. The loop slept a
flat `MINIMUM_POLL_SECONDS` from the moment of submission, so a seven-minute search made about
twenty-eight requests, nearly all of them before completion was possible. Every one is a request
against the shared service this ADR notes can block the whole institution's address. The first
wait is now the estimate, floored at the rate limit and capped at half the budget so a wild
estimate cannot turn an early finish into a wait for the deadline. Later polls return to the
floor: honouring the estimate must not become a backoff, because a job that overruns should still
be collected within one poll of finishing. This removes requests; it does not add latency.

**The two scopes could not run at once.** `startScreening` disabled both buttons, so a researcher
who wanted both answers (the dual-scope requirement this ADR implements) waited two NCBI queues
end to end. Each scope now has its own lane, its own status and its own button, and the server
runs them concurrently. The wall-clock for both answers is now the longer of the two searches
rather than their sum.

What none of this changes is the bottleneck. NCBI's queue decides when a search finishes, and
nothing here makes it decide sooner. A local BLAST backend would, at the cost of ~10 GB of
database and an installation step that contradicts [INSTALL.md](../../INSTALL.md);
`BlastClient` remains the seam where that would go.

## Costs and limits

Screening needs internet; the designs do not, and a failed screening leaves the run untouched.
NCBI can be slow or refuse service, and the client gives up after thirty minutes. The 50-candidate
cap is arbitrary in the sense that no measurement chose it.

**The order export now cites the screening.** ~~A screening lives in its own record and the order
CSV is generated from the run alone, so until the two are joined the export's answer stays the
honest one.~~ Closed 2026-10-08. The export reads the newest *complete* screening per scope and
reports its findings in the row, with the database and the date, because a result without its
provenance is unreadable six months later and a row pasted into a vendor form arrives alone. The
run itself is untouched: the screening stays its own dated record and the export reads it.

Three outcomes are preserved, which is the whole reason this took a shape rather than a boolean.
No screening is still `not-performed`, not `ok`, missing data is not a negative assay. A failed
or unfinished screening contributes nothing, because a search that did not complete is not a
search that found nothing. A screening with perfect, near-perfect or in-window hits is a
`concern`, and a concern still ships the sequences.

`NEAR_PERFECT_IDENTITIES = 17` is a reporting threshold with no biological claim behind it. It
changes what is counted, never what is kept, and the raw alignments travel beside it.

Nothing here measures binding. Sequence similarity is sequence similarity; whether a duplex
actually engages any of these loci is unmeasured, and so is whether engagement would change
expression.
