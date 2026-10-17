# Screening API contract and local reference setup

## Immediate defects and interpretation

The browser's empty POST to `/api/screenings/{screening_id}/resume` omitted the JSON media
type. The loopback HTTP handler requires `application/json` for POSTs and correctly rejected
it with HTTP 415 before the router ran. Router-only tests and mocked browser requests missed
that boundary. The client now sends `{}` with its media type. API transport errors also
return JSON, and the client reports status codes if an unexpected non-JSON response arrives.
The same-origin protection and JSON requirement remain in force.

In the reported run, each completed transcript search requested one candidate and returned
100 alignments. The genomic search also requested one candidate. The UI's batch selector
was set to one. A single candidate row is not a single off-target: expanding it shows the
reported transcript alignments. The matching-base threshold count, distinct transcript count,
and named-gene count have different units. None is a count of biologically affected genes.

The primary buttons now state the number of candidates the next search covers. Explicit
“Search all”, “Search remaining”, and “Search selected” actions distinguish selection scope
from the batch-size control. Selecting one is respected rather than silently expanded.

## Routes and naming

There is no universal officially mandated URL spelling. This API uses plural resource names,
opaque screening IDs, and HTTP methods with separate read and mutation semantics. The custom
`resume` action is a pragmatic action route; Google's API style instead uses colon-suffixed
custom methods. Either convention needs a consistent contract. Renaming this route would not
fix the failure and is unnecessary for this bug fix.

| Method and path | Meaning | Request | Success |
| --- | --- | --- | --- |
| `POST /api/runs/{run_id}/screenings` | Reserve/submit candidate screening | JSON: scope and candidate IDs, or count | 202, screening ID, reuse flag |
| `GET /api/runs/{run_id}/screenings` | List stored screening jobs | No body | 200, metadata list |
| `GET /api/screenings/{screening_id}` | Read current status and completed result | No body | 200, status record |
| `POST /api/screenings/{screening_id}/resume` | Collect the existing remote search | JSON `{}` | 202, same screening ID |
| `POST /api/screenings/{screening_id}/annotations` | Retrieve/save optional gene context | JSON: accession list | 200, annotated/unavailable entries |

Accepted screening actions return `Location: /api/screenings/{screening_id}` and
`Retry-After: 5`. Status records include `poll_after_seconds`; the browser polls local running
jobs no faster than five seconds. A local GET returning 200 means the **status request**
succeeded. Only the JSON body's `status` establishes whether screening has completed.
These local GETs are not individual BLAST searches or NCBI polls.

`screening_id` identifies a local dated job. `request_id` is the external NCBI RID. They are
not interchangeable. Every result retains scope, sequences' candidate IDs, submitted
parameters, database, time, and reporting limits. A repeated active selection is reserved
atomically; completed results are not silently reused against an unspecified current database.

Transport failures use HTTP status plus a JSON `detail` and machine-readable `code`. Resume
conflicts distinguish already collecting and already complete; the browser follows the
existing job. Invalid media types remain 415, other-origin requests remain 403, invalid
candidate selections are rejected, and missing query results never become zero-hit findings.
A versioned OpenAPI description is appropriate before offering this as a public API.

Sources: [HTTP semantics, RFC 9110](https://www.rfc-editor.org/rfc/rfc9110.html),
[Google resource-oriented API guide](https://docs.cloud.google.com/apis/design),
[Google custom methods](https://google.aip.dev/136).

## Local primary screening backend

Public remote BLAST is useful as an optional external cross-check, but is a poor primary
interactive backend here: it has externally controlled queues, mutable databases, per-RID
polling restrictions, heuristic seeding, and reporting limits. Polling changes cannot remove
those limitations. A local parser benchmark is not evidence of fast remote screening.

In particular, the current word size of seven does not guarantee enumeration of all
full-length 16/19 matches. Three mismatches can split a 19-base alignment into runs shorter
than seven exact bases. A count labelled “≥16/19” counts returned alignments; it does not
establish that every such match was searched exhaustively.

The app now implements a separately provisioned, versioned human reference/index. The design requirements are:

1. Pin the human assembly and annotation releases, reference FASTA checksums, index build
   settings, aligner version, and coordinate conventions. Provide a supported installation
   path on Windows as well as Linux; do not silently download large references on page load.
2. Batch all explicitly requested passing candidates into a short-sequence aligner suitable
   for the defined query length. Enumerate full-length exact/mismatch alignments under an
   explicit computational policy. A mismatch budget defines the search, not a validated
   biological risk threshold. Do not silently drop repetitive queries or return only their
   best alignment. Bowtie 1's documented exhaustive `-v`/`-a` mode is a candidate to benchmark;
   this is not a claim that its defaults or another aligner's heuristic mode are sufficient.
3. Search genome and transcript reference sets separately. Annotate genomic loci against a
   local transcript/TSS index rather than hundreds of online annotation requests. Keep
   full-length mismatch screening separate from gapped similarity and any RNAi-derived
   seed analysis; document what each omits and avoid calling any one a safety score.
4. Remove the candidate's own genomic locus using coordinates. Keep raw alignment counts,
   distinct loci, transcript isoforms, and named genes separate. Declare ambiguous mappings,
   unplaced contigs, truncation, and unknown annotation. Record actual nucleotide alignments.
5. Cache by normalized query sequence, scope, reference/index checksums, tool version, search
   policy, and annotation version. Repeated candidate sequences can reuse computation while
   retaining each candidate's identity. Store timestamps and distinguish cached results.
6. Benchmark a full shortlist, repeats, exact matches, mismatch positions that defeat BLAST
   seeding, reverse strands, repetitive sequences, chromosome boundaries, alternative
   contigs, and on-target removal. Use brute-force enumeration on small synthetic references
   as an independent correctness oracle. Publish genome/transcript end-to-end times only
   after measurements on the supported installation; no new runtime guarantee is implied.

This replacement needs a reference installation/index lifecycle and platform benchmark. It
is not implemented by the current HTTP fix and must not be represented as already available.
Existing remote jobs can still be collected rather than discarded.

Sources: [NCBI usage policy](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html),
[NCBI BLAST seeding explanation](https://www.ncbi.nlm.nih.gov/books/NBK20261/),
[Bowtie manual](https://bowtie-bio.sourceforge.net/manual.html).


## Implemented local setup contract

- `GET /api/references`: local installation status, pinned catalog release, supported platform,
  disk requirement, download progress and installed manifest. No external request or download.
- `POST /api/references/setup`, JSON `{}`: explicitly start/resume the background installer.
  Returns 202; duplicate clicks in the same process reuse the running task. An OS file lock
  prevents simultaneous reference preparation by multiple app processes. Offline mode rejects
  downloads. A ready current release is reused without downloading or rebuilding.
- `POST /api/references/pause`, JSON `{}`: request a safe pause. It does not delete verified files.
- `POST /api/runs/{id}/screenings` accepts `backend: "local"` or `backend: "ncbi"`. The website
  explicitly sends its selection and defaults to local. For compatibility with existing API
  clients, omitted `backend` retains the previous NCBI behavior; new clients should always
  specify it. Missing local data is an error, never an implicit remote fallback.

Local results add `backend`, `reference`, `software`, `search_policy`, `enumeration_complete` and
`details_retained`. Summary counts cover every eligible alignment; `hits` and `proximities`
contain up to 200 strongest detail alignments per candidate. This detail bound is **not**
`alignments_capped`: the search was fully enumerated and counts are complete within its policy.
E-values are null because Bowtie does not calculate BLAST E-values. The bit-score field
is also null: neither BLAST statistic is invented for local alignments.

Policy: Bowtie 1 `-f -v 3 -a --best -p 2`, both strands, all mismatch strata, no gapped or
seed-only alignments, unambiguous 19-base DNA queries. Reference windows containing uncalled
bases are excluded. Distinct transcript accessions and genome positions remain different units.
Alternate/patch contigs are searched and can represent overlapping biological loci. Their
coordinates remain contig coordinates, not falsely translated to chromosome coordinates.

Own-locus removal uses the existing coordinate-overlap rule only in genomic scope. Nearest
TSS is calculated from the midpoint of the aligned reference interval against transcript starts
in the pinned GFF, within 10 kb, with upstream/downstream signed in transcript orientation.
TSS proxies use non-partial mRNA, lnc_RNA, transcript and primary_transcript models only.
Other RNA accessions remain available for chromosome mapping, but mature small-RNA boundaries
are not automatically promoted to TSSs. These are annotated 5-prime boundaries, not measured
cell-specific starts. This eligibility rule is a conservative project choice.
The run-window flag concerns that **nearest** TSS, not every alternative transcript promoter.
Missing transcript-to-genome mappings remain unavailable. Reference description text preserves
source-predicted functions; none is experimentally validated by this search.

Genome indexing uses `bowtie-build-s --packed --bmaxdivn 8 --threads 2`; the smaller transcript reference uses `--bmaxdivn 4 --threads 2` without packed strings. One native index search or build runs at a time per app process to limit memory peaks. Duplicate query
sequences share native work; successful per-sequence compressed alignments are cached by
sequence, scope, reference source checksums, tool checksum and policy version. A changed
reference invalidates reuse. A failed process never creates a successful result/cache. Resource
limits are 30 minutes of native search, 2 GiB raw output, or 250,000 alignments per batch;
exceeding one fails explicitly, without presenting partial counts as complete. These are
engineering limits, not biological thresholds. Results do not change candidate rankings.

Files download to resumable `.part` paths and must pass pinned official checksums before use.
Bowtie archives use pinned SHA256; NCBI files use the published MD5 manifest, with SHA256 also
recorded after verification. Index completion markers are written only after successful builds and bind source checksum, tool checksum, build options and FASTA normalization policy.
Before large downloads, a tiny independent exhaustive oracle checks the native tools on both strands and all four mismatch strata. A known-reference positive control also opens each completed human index before activation. An atomic
active pointer is written only after all inputs, indexes and annotations are prepared. Old
release directories remain in place. App updates ship the reference catalog; status checks do
not discover or install mutable datasets automatically.

Validation includes actual Linux Bowtie index building and searches against small source-shaped
references, compared with an independent exhaustive Hamming-distance oracle on both strands,
including distributed mismatches, duplicate queries, coordinate removal, caching and detail
limits. This does not establish full-human runtime or native Windows compatibility; those
require machine-specific validation before performance promises.

Primary sources: [Bowtie manual](https://bowtie-bio.sourceforge.net/manual.shtml),
[official Bowtie releases](https://github.com/BenLangmead/bowtie/releases),
[pinned NCBI release files and checksums](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/).

Local records also identify the current app source/package hash and mark reuse of cached
alignments per candidate. An explicit candidate list takes precedence over the default batch
size. Native searches batch up to 50 candidates; remote website batches remain limited to 25.

`annotation_complete` means all reported loci were queried against this eligible annotation
index. It does not mean every biological promoter or TSS is known.

Annotation parsing follows [NCBI's GFF3 specification](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/): accession.version identifiers and 1-based feature bounds are preserved. The documented `partial=true` flag can describe internal or terminal incompleteness; excluding all such models from TSS proxies is deliberately conservative. The four eligible feature types are a project policy, not a validated saRNA predictor. Source annotations may be predicted or incomplete. The aligner treats reference records as linear; circular-origin-spanning alignments are not included.

Local requests capture the active reference release when accepted. That release is part of the
active-job deduplication key and is passed to the worker explicitly, so an intervening reference
update cannot silently change the reference of a queued search. Earlier installed releases are
kept available. TSS plots use retained details and are labelled as a subset whenever full counts
exceed those details.

Progress and coverage follow the selected engine/reference; old NCBI waiting jobs cannot block
new local searches or count as completed local coverage. Their saved request IDs are retained,
and can be collected when the remote engine is selected. The results table keeps the latest
saved finding per candidate across methods and labels the method/reference; it explicitly notes
when these differ from the engine selected for a new search. Switching engines stops stale
frontend polling and further queued submissions from that old frontend operation. Local status
reads use a one-second minimum while visible and honor larger server hints; NCBI-provider
limits are unchanged. Complete immutable records already loaded are reused on engine switches.
