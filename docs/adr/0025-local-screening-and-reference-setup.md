# ADR 0025: Local short-sequence screening and explicit reference provisioning

**Status:** implemented; native Windows validation pending
**Date:** 2026-10-17

Extends [ADR 0019](0019-off-target-screening.md),
[ADR 0023](0023-off-target-reporting-thresholds-and-batches.md), and
[ADR 0024](0024-human-assembly-and-collectable-searches.md).

## Decision

The website defaults to local screening. Human genome, human RefSeq RNA, GFF annotation and
assembly metadata are pinned to GRCh38.p14 / GCF_000001405.40 / RefSeq release 110. This choice
supplies a reproducible baseline; it does not claim the annotation is current. Official source
checksums and additional SHA256 provenance are retained. Reference releases are supplied with
app versions and installed separately through an explicit in-app action. Normal installation,
page load, status checks and app updates do not download a genome.

Official Bowtie 1 native tools are installed automatically after SHA256 verification:
Linux x64 1.3.1 and Windows x64 1.2.2 (release tag v1.2.2_p1). The Windows archive avoids an
additional Perl/Python wrapper or WSL setup. Its different version is recorded, not hidden.
Other platforms are marked unsupported by automatic setup rather than given a broken download.
Archive contents include upstream license files; no binaries or genome files are committed.

Indexes are built once on the user's computer. Genome preparation uses packed strings and
`--bmaxdivn 8 --threads 2` to limit memory; the smaller transcript reference uses unpacked
strings with `--bmaxdivn 4 --threads 2` to avoid unnecessary packed-string overhead. A hosted, prebuilt index bundle could reduce
preparation time, but requires artifact hosting, checksum publication and platform validation
that this project does not currently provide. This implementation does not invent such a
service. Verified downloaded inputs and completed indexes are reused after interruptions;
partial index builds are restarted, not trusted. A tiny exhaustive both-strand/mismatch self-test checks the tools before large downloads.
All six human index components and a positive-control search must succeed before an atomic active pointer exposes a release. Prior releases remain
in place during updates. Disk preflight conservatively requires 20 GiB free; native memory and
preparation speed depend on the machine.

Search policy is full-length, ungapped 19-base cores, both strands, zero to three mismatches:
`bowtie-align-s -f -v 3 -a --best -p 2`. No `--strata`, `-k`, `-m`, masking or species-filter
heuristic suppresses eligible alignments. Uncalled reference bases are excluded. All mismatch
strata are counted, including distributed mismatches that may have no exact seven-base seed.
No E-value is invented for Bowtie output. Each retained nucleotide footprint is checked against
actual reference bases and oriented into the original query's 5′ to 3′ direction.

Genome and transcript counts stay separate. The established own-locus coordinate-overlap rule
applies only to genome matches. Local GFF transcript starts supply the nearest TSS proxy within
10 kb, signed in transcript orientation; this distance is descriptive. Proxies are annotated
5-prime boundaries of non-partial mRNA, lnc_RNA, transcript and primary_transcript models.
Mature small-RNA starts are retained for accession mapping but are not automatically treated
as promoters. The eligibility rule is a conservative project choice, not empirical evidence
that every model start is a cell-specific TSS. The run-window flag
concerns that nearest TSS, not every alternative promoter. Transcript accession mappings supply
chromosome context without online gene-symbol lookups. Alternate/patch sequences are searched,
but their coordinates remain contig coordinates and are not translated to primary chromosomes.

All enumerated alignment counts are retained; up to 200 detail alignments per candidate are
stored. This display/storage boundary is distinct from a truncated search. Native failure or
resource-budget exhaustion fails the job without returning partial counts or a false zero.
Resource budgets are documented engineering choices, not biological thresholds.

Successful per-sequence alignments are cached by query, scope, verified reference identities,
tool archive hash and policy version. Duplicate candidates share native work. The app serializes
native searches and index building to limit memory peaks. Search/annotation computation is
local. NCBI BLAST remains explicitly selectable with its existing queue, reporting limits and
saved-RID recovery; a missing local reference never silently switches to NCBI.

## Evidence and limits

This policy is computationally defined, not an experimentally validated saRNA risk model.
It excludes gapped/bulged interactions, shorter seed effects, sample-specific variation,
chromatin accessibility, expression context and biological binding. Reference-alternate hits
can represent overlapping biological loci. No candidate ranking changes.

Native Linux integration tests prepare indexes from tiny source-shaped FASTA/GFF/report files
and compare actual Bowtie output with an independent exhaustive Hamming-distance oracle. They
cover both strands, distributed mismatches, duplicates, caching, own-locus removal, detail limits
and exports. Windows execution cannot be inferred from Linux tests and remains a separate check.

Sources: [Bowtie manual](https://bowtie-bio.sourceforge.net/manual.shtml),
[official Bowtie releases](https://github.com/BenLangmead/bowtie/releases),
[NCBI pinned release and checksums](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/).

The reference-pack release ID must change when inputs, native tool archives, derived schema or
preparation rules change. Completion markers bind actual input/tool hashes and build options;
an incompatible staged index is rebuilt. Existing complete records do not count as coverage
for another engine or reference release. Old NCBI queues remain retained and do not block the
local engine. New search controls and saved evidence methods are labelled separately.
