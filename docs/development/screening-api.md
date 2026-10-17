# Screening API contract and backend recommendation

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

## Recommended primary screening backend

Public remote BLAST is useful as an optional external cross-check, but is a poor primary
interactive backend here: it has externally controlled queues, mutable databases, per-RID
polling restrictions, heuristic seeding, and reporting limits. Polling changes cannot remove
those limitations. A local parser benchmark is not evidence of fast remote screening.

In particular, the current word size of seven does not guarantee enumeration of all
full-length 16/19 matches. Three mismatches can split a 19-base alignment into runs shorter
than seven exact bases. A count labelled “≥16/19” counts returned alignments; it does not
establish that every such match was searched exhaustively.

Replace the primary path with a separately provisioned, versioned human reference/index:

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
