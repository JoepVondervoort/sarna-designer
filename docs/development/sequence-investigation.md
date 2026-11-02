# Sequence investigation and design reports

The Candidates workspace provides these workflows through the promoter map’s tool buttons.
See the [interface hierarchy](interface-hierarchy.md) for the page order and disclosure behaviour.

The two related sequence workflows are:

- **Evaluate my saRNA & reference duplexes** imports externally designed duplex cores,
  locates every exact occurrence, and compares matched cores with the saved run's rules.
  An unmatched core of the configured target length receives a separate sequence-only
  assessment using the run's saved rule settings. It has no design rank, genomic placement
  or off-target screen. TSS, promoter flanks and CpG context remain unavailable. Unknown
  strand roles receive separate assessments of both distinct orientations. Other core
  lengths retain sequence/thermodynamic inspection without a rule rating.
- **Search and annotate the promoter** searches 4-10,000 bases in the loaded promoter.
  Choose same orientation, reverse complement or both. RNA U is normalized to DNA T;
  mixed T/U, ambiguity codes and mismatches are not supported. Overlapping occurrences
  are retained, palindromes are deduplicated, and searches exceeding 500 matches require
  a longer query. No match describes this window only.

Click a search result to zoom to it and populate a region's coordinates. Give the region
an informative name and source/citation, then save it. Manual regions use **1-based,
inclusive** first/last bases. Saved annotations and APIs use **0-based, half-open**
coordinates. Annotations appear in purple on the map; current search results appear in
amber. Searches are temporary; named regions persist.

Motif display choices support search by TF name or matrix accession and multiple checked
factors. Show all and hide all are explicit choices. Existing focus in the motif evidence
table highlights one factor within that selection. Display choices never rescore candidates.

The distance tool compares candidates, imported matches, regulatory features, displayed
motifs, named regions and an in-window reference TSS. It reports interval overlap or the
number of bases between interval edges, plus signed start displacement in transcript
orientation. Touching intervals have zero gap and zero overlap. The TSS is a point at the
start of the +1 base; a point inside a region is reported explicitly. These geometric
measurements do not establish regulatory effects.

## Saved state

Annotations and motif choices are stored in a separate SQLite table, with a revision check
that rejects stale writes. They do not change scientific snapshots, ranking or replay
fingerprints. On conflict, reload saved map settings before retrying; pending annotation
form fields remain available. Whole-run and group JSON/HTML exports include this state.
Offline HTML permits inspection, sequence search, measurements and report downloads;
annotation editing requires the running app.

## Downloads

The Export workspace adds:

- **Design summary HTML:** selected candidates, configured versus assessable preference
  counts, unavailable gates, individual verdicts, promoter/reference settings, selected
  candidates' screening counts and coverage, displayed motifs, annotations, citations and
  software/reference provenance. Detailed inputs are expandable. Open the document and
  print to PDF. Whole-run JSON remains the complete machine-readable archive.
- **Promoter alignment HTML:** the full transcript-oriented promoter in 80-base blocks,
  with selected saRNA target sites, displayed predicted TF motifs, imported matches and
  named regions aligned to the same coordinates. Motif rows show promoter DNA at the
  predicted interval, not a TF protein sequence. A coordinate/source table retains
  orientation information.
- **Coordinate CSV:** the same features, with explicit zero-based, half-open promoter and
  genomic intervals, assembly, promoter strand, orientation and sources. Minus-strand
  promoters map offsets using the genomic window's end. Unanchored promoters leave genomic
  coordinates blank.

These downloads follow explicit candidate checkboxes, or all design-passing candidates
when no explicit selection exists. A managed empty screening selection remains empty.
The motif checklist applies to the alignment, coordinate table and summary regulatory
context; single-factor highlighting does not discard other checked motifs.

Ranking still uses the number of preferences met, followed by the saved thermodynamic
tie-break. The displayed assessable denominator can vary when context is missing. The
configured denominator is shown separately; the percentage bar is descriptive and does
not determine rank. No new efficacy model or probability is introduced.
