# 0029: Plot-first screening and explicit analysis exports

The screening workspace opens with paired transcript and genomic bars for each candidate,
stacked by saved full-length 19/19, 18/19, 17/19 and 16/19 identity counts. Both scopes share
one linear count scale. Unknown counts have a separate marker; incomplete searches have
an explicit dashed outline. Hover/focus shows design rank, score, sequence and saved search
context. Clicking opens that candidate's detail row in the selected scope.

One set of reference, cohort and ordering controls drives both the plot and detailed table.
The default is all screened candidates in original design-rank order. Optional cohorts use
design ranks 1–10 or 1–15. Count sorts use the selected scope and preserve rank for ties;
missing and incomplete evidence stays last. Evidence from different references is not pooled.
The detailed table, secondary plots and methods remain available through disclosures.
The promoter landscape is now always visible at the top of Candidates, superseding its
expandable placement in ADR 0028. Biological context and motif predictions remain expandable.

Each export format has a visible content description:

- Excel adds screening summary, retained alignments, dated searches and limitations sheets.
- Screening CSV contains per-candidate/scope/reference counts and search provenance, without
  alignment detail. Unscreened/failed/running counts are blank, distinct from complete zeroes.
- Ordering CSV retains ordered duplexes and brief screening findings; core CSV contains
  sequences/design scores without screening analysis.
- JSON and offline HTML preserve the whole run and every dated screening record. Retained
  alignments remain a subset, not a reconstruction of all reported matches. The offline report
  can reproduce the screening CSV and ordering CSV without contacting the service.
- SVG saves the current paired-bar view, including its legend and ordering metadata.

Selected/default-passing candidates apply to Excel and the three CSV downloads. JSON and
HTML preserve the whole run; SVG follows the plot controls. Original stored design snapshots
and the lightweight interactive API remain unchanged. No plot or export changes scientific
ranking or establishes biological specificity.
