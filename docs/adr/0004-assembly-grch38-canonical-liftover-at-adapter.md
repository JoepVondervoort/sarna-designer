# ADR 0004 — GRCh38 is canonical; GRCh37 conversion happens inside the bigWig adapter

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The blocking open question from Phase 0 is now closed, and the answer is unfavourable. Verified
2026-10-06 by reading the bigWig headers directly (see
[`legacy-inventory.md` §12](../scientific/legacy-inventory.md)):

| Source | Assembly | Contig naming |
|---|---|---|
| Ensembl REST (default) | **GRCh38** | `1`, `2`, `X` |
| AF22 bigWigs, GSE108506 | **GRCh37/hg19** | `chr1`, `chr2`, `chrX` |

All six bigWigs agree (`chr1` = 249,250,621; `chrX` = 155,270,560), consistent with the 2017
deposition date. The legacy EpiMarker script pins no assembly and therefore indexed hg19 signal
with GRCh38 coordinates — defect **B15**, severity S1. GRCh37→GRCh38 shifts are locus-dependent
and reach megabases, while the pipeline samples 50 bp bins in a ±30 kb window, so the signal
read was frequently from an unrelated locus.

This is not a question we can defer to implementation. Both assemblies have to exist somewhere
in the system; the only choice is **where the boundary sits**.

Three options were considered:

1. **Pin Ensembl to GRCh37** (`grch37.rest.ensembl.org`). No conversion anywhere. Rejected:
   it puts the entire tool, including every coordinate a user reads or exports, on a deprecated
   assembly, to accommodate one annotation input.
2. **Lift the bigWigs to GRCh38 once** with CrossMap. Rejected: ~890 MB of reprocessing,
   liftover of pre-binned signal is lossy in a way that is hard to characterise, and the derived
   files would silently diverge from the published GEO record that cites them.
3. **Lift query coordinates per request.** Chosen.

## Decision

**GRCh38 is the single canonical assembly of this application.** Every coordinate in the
domain model, the API, the UI and every export is GRCh38. There is no second coordinate space
in the design or scoring layers.

GRCh37 exists in exactly one place: **inside the bigWig adapter**, which is the only component
permitted to know that the AF22 tracks are hg19. For each signal query it

1. receives a GRCh38 interval,
2. converts it to GRCh37 via a UCSC chain file,
3. normalises the contig name to the `chr`-prefixed form the files use,
4. reads the signal, and
5. returns values tagged with the assembly they were read on.

Consequent rules, all of which are testable:

- **Every genomic coordinate carries its assembly in the type**, rather than being a bare `int`.
  A coordinate without an assembly is not representable. This is the mechanism that makes B15
  impossible to reintroduce.
- **No external genomics call may rely on a service default.** The Ensembl client pins and then
  *asserts* GRCh38 in the response. A default that changes upstream must break a test, not
  silently alter results.
- **Unmappable intervals are reported, never guessed.** Where the chain file has no GRCh37
  equivalent — or maps an interval to a different contig, or splits it — the adapter returns an
  explicit "no data for this region" result. It never falls back to unconverted coordinates.
- **Contig naming is normalised at the adapter only.** `chr`-prefixing is a file format detail
  of these particular tracks and does not leak into the domain model.

## Consequences

- A chain file (`hg38ToHg19.over.chain.gz`, UCSC) becomes a required data dependency. It is
  version-pinned and checksummed; it is data, not code, so it lives outside the repository
  alongside the bigWigs and is documented in the environment notes.
- Signal lookups gain a conversion step. This is negligible against the bigWig read itself, and
  conversions are cacheable per TSS window.
- Some loci will legitimately return no epigenomic annotation. This is correct behaviour and
  must be surfaced in the UI as "not mappable to the AF22 assembly" — *not* as zero signal.
  Zero signal and absent signal are scientifically different statements and must never render
  the same way.
- The annotation layer must state the provenance of what it shows: cell line AF22, GSE108506,
  GRCh37, lifted. Under [ADR 0003](0003-rule-provenance-is-first-class.md) this is the same
  obligation design rules carry, applied to data.
- Historical outputs from the legacy EpiMarker script cannot be reproduced by the new engine,
  because reproducing them would mean reproducing B15. Any comparison against them is a
  comparison against known-wrong values and must be labelled as such.
- Should other assay tracks be added later on GRCh38 natively, they need no adapter conversion;
  the boundary is per-source, not global.
