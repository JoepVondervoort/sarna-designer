# Off-target screening reliability and interpretation

Screening reports sequence matches, not measured binding, specificity, or changes in gene
expression. It does not alter design ranks. Transcript and genomic coverage are independent.
The UI reports candidate coverage and separates waiting, failed, unavailable, and negative
findings. A completed search does not imply that both scopes or every candidate were searched.

## Evidence and reporting

The main threshold is labelled by the matching-base count (≥16/19 for a 19 nt core). Applying
80% to candidate off-targets is a project interpretation; Wang et al. discuss that threshold
for control duplexes. The count is not a validated saRNA risk cutoff. Alignment counts,
distinct transcript accessions, named genes from source descriptions, and distinct genomic
alignment loci have different units. Named genes do not include unnamed hits and are not
counts of affected genes. Genomic placements report the nearest TSS within 10 kb; they do not
establish promoter activity or exhaustively consider every alternative transcript start.

The transcript search is not a comprehensive seed-risk screen. Evidence for seed-mediated
siRNA silencing is relevant context, not validation of a saRNA activation predictor.

Sources:

- [Wang et al., 2015](https://doi.org/10.14440/jbm.2015.39)
- [Jackson et al., 2006](https://doi.org/10.1261/rna.25706)
- [NCBI usage guidelines](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)
- [NCBI Common URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html)

## Search integrity and recovery

The client requests `XML2_S` (single-file XML2). `XML2` is the multiple-file format and
can return a ZIP archive for multi-query searches; passing that archive to an XML parser
caused the observed invalid-token error at line 1, column 2. The parser supports XML2 and
recorded classic XML, including XML2 strand labels. Unreadable downloads are retried once
using the same RID, then kept collectable with a retrieval warning. On restart, previously
failed XML-parse jobs with retained request IDs and candidate lists become collectable again.
No fresh search is submitted to recover them. It rejects unexpected XML,
missing query identifiers, duplicate queries, and incomplete query coverage. Explicit
READY/no-hits status is distinguished from missing status fields. A server sequence limit
reached or a client alignment limit exceeded marks counts as potentially incomplete.

Aligned query and subject strings, including gaps, are retained. The footprint shows the
DNA sequences in BLAST alignment orientation, not an RNA duplex pairing diagram. T corresponds
to U in RNA. Old records without aligned strings show unavailable subject bases; they are not
reconstructed from identity counts. A fresh search is needed for those missing bases.

Genomic chromosomes come from assembled chromosome accessions. Transcript chromosomes load
on expansion using an Ensembl gene-symbol lookup and are saved with their source and date.
This is gene context, not a transcript-to-genome coordinate mapping. Failed or ambiguous
lookups remain unavailable. Source descriptions, including predicted function, are preserved.
Archived pages show only annotations already recorded and do not fetch more data.

The remote request ID, candidate list, and submitted parameters are saved immediately after
acceptance, permitting collection after restart. Identical active selections within one run
and scope share a reservation, including across browser tabs. Resume is claimed atomically.
An unfinished batch stops subsequent submissions; after collecting it, the user can screen
the remaining candidates. Saved completed results are reused in the panel with their original
dates. A repeat search deliberately produces a new dated observation; there is no silent reuse
of a previous result against a database whose current version is unknown.

## Performance and limits

Requests share process-wide scheduling: at least 10 seconds between contacts, at least 60
seconds between status polls for one RID. Multiple independent app processes do not share a
limiter. Set SARNA_NCBI_EMAIL to provide NCBI's requested contact; the tool name is always sent.
Retrieval retries transient failures, with bounded Retry-After handling. Submissions are not
retried automatically because an interrupted response may already have created a job.

The first status check follows the policy interval, not RTOE: an estimate is not a lower bound.
A worker collects for up to three minutes for genomic searches or ten minutes for transcript
searches, then leaves an unfinished remote job collectable rather than occupying a worker for
half an hour. A long estimate can hand the job back after the first status check.
Waiting shows an indeterminate bar and a labelled estimate; only countable annotation progress
uses a percentage. The server controls queue/search time. No promise of faster remote genomic
searches is implied by these changes.

TSS annotation uses cached megabase blocks, at most 300 placements, a 60-second soft budget
checked between lookups, and round-robin allocation of strongest hits per candidate. A single
in-flight lookup can exceed the soft budget. Failed blocks are remembered as failures for that
attempt, and the service stops further remote annotation attempts after three failures.
Unknown/unannotated context never becomes zero window hits. Source response caches are reused;
optional chromosome lookups do not delay the primary search result.

To measure local processing without consuming NCBI capacity:

```bash
.venv/bin/python scripts/benchmark_offtarget.py --queries 25 --alignments 120
```

A development measurement on 16 October 2026 parsed 3,000 synthetic chromosome alignments in
0.0385 seconds and performed 300 cached placements in 0.0011 seconds (one underlying stub
lookup). This is a local processing benchmark, **not** an end-to-end genomic search benchmark.
For a deployment benchmark, record submission, first READY, report retrieval, and annotation
completion separately. Remote queue time must not be attributed to parser performance.
A pinned local reference/backend is a future deployment option for predictable throughput;
it is not installed or silently substituted here.

## Automatic collection while a run is open

An unfinished remote job is checked again automatically after a one-minute pause while its
run remains open. The browser calls the resume endpoint for the existing screening ID; it
never submits the sequences again. The server's shared scheduler still enforces NCBI pacing.
Timers are deduplicated and guard both run and screening IDs, so navigating to a different run
cannot trigger a stale collection. Closing the run pauses browser-driven checking; reopening
restores it. Unreadable reports require an explicit retry instead of an endless automatic loop.
If another tab is already collecting the same search, this tab follows its local status.
The count in each scope explicitly means candidates with completed results, not zero hits.
