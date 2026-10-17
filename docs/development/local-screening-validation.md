# Local screening validation, 17 October 2026

The implemented search was checked against the full pinned human reference, not just tiny
fixtures or a parser benchmark. This is computational validation of a defined sequence search;
it is not experimental validation of saRNA specificity, binding, activation or silencing.

## Reference and policy

- GRCh38.p14, GCF_000001405.40, NCBI RefSeq annotation release 110.
- Official genomic/RNA FASTA, GFF and assembly report checksums are pinned in
  `sarna.offtarget.references`; SHA256 identities are also recorded after verification.
- Full-length 19-base DNA cores, both strands, zero to three mismatches, no gaps.
- Only callable A/C/G/T reference windows; linear records, without circular-origin wrapping.
- Bowtie 1.3.1 on Linux x64, two search threads, an 8 GB machine. Windows uses the separately
  pinned official 1.2.2 archive and remains untested natively in this environment.

## Independent complete-reference check

[scripts/local_hamming_oracle.cpp](../../scripts/local_hamming_oracle.cpp) implements a separate
rolling, two-bit Hamming counter. It does not use Bowtie, an FM index, screening summaries, or
the normalized FASTA to determine matches. It scans the original verified NCBI FASTA streams,
resets at record boundaries and uncalled bases, and compares each window against both the
query and its reverse complement. A mismatch is counted once per differing nucleotide using
an XOR/bit-pair mask and population count. The compiled counter was first checked against a
simple exhaustive Python oracle on a small known-answer reference.

Three computational examples were then checked against the original whole-genome and
transcript FASTA streams. These are test cases, not recommended experimental reagents.

| Example core, DNA 5′–3′ | Genomic counts at 0 / 1 / 2 / 3 mismatches | Transcript counts at 0 / 1 / 2 / 3 mismatches |
| --- | --- | --- |
| GTGAAGCCATAGACTCGTA | 0 / 1 / 8 / 228 | 0 / 1 / 3 / 50 |
| GCTCTACTCTCAGTTCCAA | 1 / 1 / 64 / 884 | 0 / 0 / 8 / 229 |
| CACGCTCCCAGCCCGAACG | 1 / 0 / 16 / 129 | 9 / 0 / 0 / 95 |

All four mismatch strata from the native cached alignments matched these independent totals
in both scopes. The reference-positive example also matched the original genomic FASTA at
`NC_000017.11:7687410–7687429`, in 0-based half-open coordinates. The application removed
that single own-locus alignment when the candidate carried those coordinates, giving 145
reported genomic alignments. It retained all nine transcript matches, which represent RefSeq
transcript records rather than nine different genes.

The application's complete-summary counts, perfect counts, near-perfect counts, retained
footprints and own-locus handling were checked as well. A newly queried reverse complement
also passed after the completed reference was moved to the app's default data directory.

## Measured screening time

For the three examples together, with existing prepared indexes:

| Scope | First screening in this test | Repeat using cached alignments |
| --- | ---: | ---: |
| Genome, including local annotation | 4.167 s | 0.062 s |
| Transcripts, including local annotation | 0.738 s | 0.029 s |

These timings are observations on this machine, not a guarantee for other machines, colder
filesystem caches, larger batches or repetitive queries. No NCBI BLAST submission or Ensembl
annotation request was needed. Of the 395 distinct transcript accessions in this cohort,
393 had primary-chromosome context and two were explicitly unavailable.

The first preparation is substantially longer: genome indexing took about 39 minutes and
transcript indexing about 28 minutes on this machine, plus downloading and preparation.
This cost is paid once per reference pack; normal app updates reuse it. Interrupted downloads
and completed staged indexes were successfully reused during validation. Keep the app process
running during the initial setup. The installed reference plus verified download cache occupies
about 9 GiB here; preparation conservatively requires 20 GiB free.

## Other checks and limits

The automated checks cover native tool self-tests, independent strand/mismatch oracles,
duplicate queries, saved own-locus alignments, resumed downloads, checksum rejection, staged
index reuse, cached searches, bounded detail storage, exports, and missing-reference failures.
Browser checks exercise setup/pause/resume, unavailable local data, local-versus-NCBI coverage,
null BLAST statistics, complete counts versus retained details, and mixed saved methods.
The built wheel was inspected for its source/resource hashes and contains no genome, index or
native binary payload; setup downloads those explicitly.

Counts remain complete within the declared policy when a job succeeds. At most 200 strongest
detail alignments per candidate are stored, and plots using those details are labelled as a
subset. Resource-budget exhaustion fails explicitly instead of returning partial counts or
false zeroes. A zero within this policy does not establish biological specificity.

Promoter-distance proxies use non-partial mRNA, lnc_RNA, transcript and primary_transcript
models. Other RNA accession mappings are preserved without automatically treating mature
small-RNA boundaries as promoters. The local GFF index contains 188,352 mapped RNA features,
185,486 eligible under this conservative proxy policy. These are annotation-derived 5′
boundaries, not measured, cell-specific TSSs. The nearest eligible start is sought within
10 kb, with transcript-oriented signed distance; alternative starts, cell activity and binding
remain untested. Alternate/patch records can represent overlapping biological loci.

## Repeating the independent count

On Linux with a C++ compiler, build the optional developer-only oracle:

```bash
g++ -O3 -march=native -std=c++17 scripts/local_hamming_oracle.cpp -o /tmp/sarna-oracle
```

Pass the original FASTA stream and each query plus its reverse complement:

```bash
gzip -dc genomic.fna.gz | /tmp/sarna-oracle /dev/stdin \
  GTGAAGCCATAGACTCGTA TACGAGTCTATGGCTTCAC \
  GCTCTACTCTCAGTTCCAA TTGGAACTGAGAGTAGAGC \
  CACGCTCCCAGCCCGAACG CGTTCGGGCTGGGAGCGTG
```

Standard output gives each argument's index followed by counts at 0, 1, 2 and 3 mismatches.
Sum the two orientations for each example. Standard error prints up to 100 exact positional
controls for the third example's argument pair, as used in this validation. This utility is
not an app dependency and does not run during installation.

Sources: [pinned NCBI reference and checksum manifest](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/),
[Bowtie manual](https://bowtie-bio.sourceforge.net/manual.shtml),
[NCBI GFF3 specification](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).
