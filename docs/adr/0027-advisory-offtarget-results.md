# 0027: Off-target findings accompany the top designs

Supersedes the automatic exclusions and export restrictions in ADR 0026, following user
review on 18 October 2026. The original design gates and rank still determine the top 10/15.
All design-passing candidates remain selectable, including those with sequence matches or
incomplete screening. Off-target screening is advisory, not a validated safety decision.

The top list shows complete per-identity counts for transcripts and genome at 19, 18 and 17
matching bases. Missing counts are never displayed as zero; incomplete evidence is labelled.
Historical Python selection assessments remain in saved reports and exports for provenance,
but their exclusion reasons are now review notes. They do not filter or reorder the top list,
remove selected candidates, or block export. Rank 1 also carries a compact screening disclosure.

The genomic distance plot remains available in a collapsed disclosure. It states the number
of plotted points divided by all reported alignments, plus the retained-detail count. This
fraction is plot coverage: retained alignment details within 10 kb of an eligible annotated
start site anywhere in the reference. It is neither the fraction of all alignments near a TSS
nor proximity specifically to the analysed gene. Unretained alignments can also be near TSSs.
The 200-detail cap and nearest-start annotation limit remain explicit. Counts and original
screening records are unchanged.


## Candidate comparison plots

Three expandable plots accompany the displayed top 10/15 candidates and preserve design rank.

- Match-count bars use saved 19/18/17 identity counts with one shared linear scale across
  both scopes and all displayed candidates. Exact counts remain visible beside each bar.
  Incomplete evidence is labelled; unavailable counts have no bar and are never zero-filled.
- The criteria heatmap uses saved per-rule preference verdicts, rather than deriving cells
  from the aggregate score. Met, missed and not-evaluated states have distinct symbols and
  colours. Hover/focus reveals the saved explanation. Mandatory gates and tie-breakers are
  outside this preference comparison.
- The mismatch map has candidate and scope selectors. At most 40 strongest usable retained
  alignments are drawn, with shown/retained/reported counts. Only full-length, ungapped hits
  with complete mismatch positions are usable. Columns follow the target DNA 5′ to 3′ on
  either subject strand; they do not imply guide or seed orientation. Missing details,
  partial and gapped alignments are not rendered as matches. This is a selected subset,
  not a biological risk estimate.

Selectors preserve open disclosures and keyboard focus; wide matrices scroll inside their
panels on narrow screens. All three plots work from saved data, including frozen HTML exports,
without extra searches or annotation requests.
