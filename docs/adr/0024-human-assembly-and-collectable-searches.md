# ADR 0024: Search the human assembly, and collect a slow search instead of losing it

**Status:** accepted
**Date:** 2026-10-11
**Deciders:** Joep Vondervoort

Extends [ADR 0019](0019-off-target-screening.md) and [ADR 0023](0023-off-target-reporting-thresholds-and-batches.md).
What a screening reports does not change here. What it searches, and what happens when it runs
long, do.

## Context

The promoter scope took so long it was unusable. A real run reached the full half-hour wait
budget, was recorded as a failure, and returned nothing, while NCBI was still working on it.

Looking into it turned up something simpler than a performance problem.

**`refseq_genomic` is not a database NCBI offers.** It appears zero times in NCBI's BLAST
database directory and zero times in the database list on NCBI's own blastn page. The URL API
accepted it and returned a request id anyway, so the search ran, against something that is
certainly not the human genome.

The sizes, from NCBI's own published metadata, are the whole argument:

| Database | Letters | Sequences |
|---|---|---|
| `GPIPE/9606/current/ref_top_level` (human GRCh38) | 3,272,089,205 | 639 |
| `refseq_rna` | 251,911,205,044 | 84,802,138 |
| `core_nt` | 1,086,004,251,965 | 131,744,797 |

**Database size turned out not to be the reason, and measuring said so.** The first guess was
that the genomic scope was searching far too much sequence. A timed probe refuted it: one
19-mer against the retired name and one against the human assembly were *both* still running
after 25 minutes. Reading NCBI's own RTOE, which comes back instantly at submission, gave the
real answer:

| Database | NCBI's own estimate |
|---|---|
| `GPIPE/9606/current/ref_top_level` (human assembly, 3.3 Gbp) | **28,502 s, about 8 hours** |
| same, word size 11 instead of 7 | 28,747 s |
| same, `EXPECT=10` instead of 1000 | 28,765 s |
| `human_genome`, `refseq_genomes`, `RefSeqGene`, `refseq_select_rna` | 27,757 to 28,745 s |
| `refseq_rna` + Entrez human (252 Gbp) | **30 s** |
| `core_nt` + Entrez human (1,086 Gbp) | **30 s** |

`core_nt` is 332 times the size of the human genome and is estimated at thirty seconds. So this
is which NCBI backend serves a database, not how big it is: the genome and assembly databases
are queued, and no choice of word size or E-value moves them. Lowering `EXPECT` would not have
helped and would have cost findings: a 16 of 19 match on a 3 Gbp database has an E-value in the
hundreds, so Wang's 80% threshold is only visible *because* `EXPECT` is 1000.

Searching `core_nt` restricted to human genomic molecules does run on the fast path, in 27 to
40 seconds, but its hits land on RefSeqGene records, clones and patent sequences rather than on
chromosomes, so they carry no coordinates to place against a start site; restricting further to
chromosome-sized RefSeq records returns nothing at all, because `core_nt` does not hold the
assembly. There is no fast remote route to a human chromosome-level search.

Two further things came out of the same reading. `ENTREZ_QUERY`, which is how this project
confined searches to human, is **not a documented URL API parameter**, and NCBI's own page says
parameters it does not list "should not be relied upon by users of this API". And placement was
measured at 1,022 Ensembl requests for 3,000 scattered hits, a minute or two, real but not the
problem.

## Decisions

### The genomic scope searches the human reference assembly

`GPIPE/9606/current/ref_top_level`, which is what NCBI's own page submits for "Human genomic
plus transcript" minus the transcript half, since this project searches transcripts as its
other scope. No `ENTREZ_QUERY` is sent with it: the assembly holds one species, so the filter
has nothing to do and the undocumented dependency goes away.

**This fixes correctness, not speed.** The old name was retired *and* slow; the new one is
current, human-only, chromosome-level and still slow, because every genome database at NCBI is.
What it buys is that when the search does come back, its hits are on chromosomes and can be
placed against annotated start sites, which is the whole point of the scope.

**The cost.** This is the *reference* assembly, so alternate haplotypes and patch scaffolds are
out of scope, and a duplex matching only a patch sequence will not be reported. It is also a
different search set from whatever the old name resolved to, so a genomic screening saved
before this change is **not comparable** with one saved after it. That is survivable only
because every screening records the database it used and its date, which ADR 0019 already
required for exactly this reason.

Rejected: `refseq_genomes`, estimated at the same 27,757 s and spanning every species.
Rejected: `core_nt` restricted to human genomic, which is fast but returns unplaceable records.

**Open, and the only route to an interactive promoter search: a local BLAST+ backend.**
ADR 0019 rejected a local install because it meant "BLAST+ binaries plus roughly 10 GB of
database". The human genome BLAST database is **1,157,233,426 bytes**, about 1.16 GB, so that
figure was wrong by roughly nine times and the decision deserves revisiting on the real number.
`sarna.offtarget.ncbi.BlastClient` was built as a seam for exactly this. Not decided here.

### Each scope records how it was confined to human, and the interface says which

The two are not the same guarantee. The assembly is human **by construction**. `refseq_rna` is
human **by filter**, through a parameter NCBI does not document. `SearchParameters` carries
`species_restriction`, and the panel states the difference in words rather than labelling both
"human" and leaving a reader to assume the stronger reading.

**Open, not decided here:** whether the transcript scope should move to
`GPIPE/9606/current/rna` for the same reasons. It would be faster and would remove the last
reliance on `ENTREZ_QUERY`, but RefSeq RNA's wider coverage is arguably the point of that
scope, and the question deserves its own measurement.

### A search NCBI estimates beyond the wait budget is handed back after one look

RTOE arrives with the submission. When it exceeds the time this app is willing to wait, polling
for half an hour only confirms what the first reply already said. One poll goes out at the rate
limit, because an estimate is not a measurement and a finished job should be collected rather
than handed back, and after that the request id is returned with the estimate in plain words.

The interface says so before the search is started, so nobody sits watching a bar for an errand
that takes hours.

### Running out of patience is a fourth state, not a failure

NCBI holds a request id for about a day. The wait budget expiring means this client stopped
waiting; it does not mean the search failed, and it certainly does not mean nothing was found.
A screening in that position is stored as `waiting` with its request id, and the interface
offers to go back for it. Resuming collects the original search and never submits a second.

**The cost.** A fourth status that every reader of a screening record has to handle, and a row
that is neither a result nor an error. The alternative is what happened: a job NCBI was still
running became unreachable, and the only way forward was to queue an identical duplicate
against a shared public service.

A `waiting` screening also survives an application restart, where a `running` one does not.
That asymmetry is deliberate: `running` is this process's work, `waiting` is NCBI's.

### The per-candidate cap keeps the strongest, not the first

BLAST groups a report by database sequence. On a chromosome-level database one sequence carries
thousands of alignments, so taking the first N in document order fills the budget with
chromosome 1's weak partial matches and never reaches a strong alignment on chromosome 7. The
cap is now a bounded best-of by bit score.

This was harmless on `refseq_rna`, where each transcript is small, and became a real defect the
moment the genomic scope pointed at chromosomes.

### Placement is bounded at the 300 strongest, and says so

Each genomic hit costs an Ensembl request per uncached megabase block. The strongest are placed
first; the rest are reported with no distance, and the record states how many. They are still
counted and still listed.

**The cost.** A weak alignment sitting exactly on another gene's start site can go unplaced, and
the bound is a round number with no principle behind it. It is a reporting limit, not a finding,
and it is stated in the record rather than left invisible.

### The progress bar claims a fraction only where one exists

NCBI reports a search as started or finished. There is no completion percentage to be had, and
the only forecast is RTOE, its own estimate at submission.

So: the placing phase shows a true fraction, because it has one. The waiting phase shows elapsed
against RTOE, worded as an estimate and **capped at nine tenths**, and once the estimate is
passed it stops filling and says the search is taking longer than NCBI predicted. A bar sitting
at 99% for twenty minutes would be worse than no bar.

**The cost.** The bar is honest rather than smooth, and during a long wait it is an animation
with a number of seconds beside it rather than a sense of how close the end is. That is all the
information there is.

## Consequences

- Genomic screenings saved before this change name `refseq_genomic` and are not comparable with
  later ones. Nothing migrates; the database and date in each record are what tell them apart.
- `run_screening` splits into submission and `collect_screening`, so a search can be collected
  from a request id without being re-run.
- The store gains a `progress` column and a `waiting` status, at schema version 5.
- The unit suite is now verified to reach no external service by running it in a network
  namespace with no network at all, rather than by inspection.
