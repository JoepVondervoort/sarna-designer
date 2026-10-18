# 0028: Four persistent workspace tabs

Following user approval of the four-tab layout, Design, Candidates, Off-target screening and
Export replace the previous scroll-to-section workflow. Exactly one tab panel is visible;
controls remain mounted so switching tabs does not reset the draft, selected candidates,
notes, expanded plots, or background screening jobs. Tabs support arrow keys, Home and End,
linked tab/panel labels and managed focus. A compact common header names the saved gene,
reference, run and explicit selection count. Screening activity remains visible in its tab.

- Design owns target/reference selection, gene and expression context, cell context and all
  settings. Numeric rule overrides sit in an Advanced settings disclosure. Screening
  reference preparation is available in its own tab even before a design has been run.
- Candidates opens after a new or reopened run. It owns the top list with compact screening
  counts, the preference heatmap, full candidate inspection and researcher notes. The
  top-ranked duplex and promoter/context tracks are expandable. The ranking-performance
  headline remains visible, with the detailed measurement and limitations expandable.
- Off-target screening owns reference preparation, search controls, coverage/progress,
  complete result counts, retained alignment details, match-count bars, mismatch maps and
  the TSS-distance plot. Candidate screening links activate this tab and open that
  candidate's saved detail row. Searches continue when another workspace is active.
- Export previews the candidates that the sequence exports will include, showing the
  engine-produced ordered duplexes. It owns downloads, group exports and reproducibility
  records. With no explicit selection, the existing all-design-passing export behaviour
  is stated. Download formats and the scientific record are unchanged.

Scientific ranking, screening counts and export selection rules are unchanged by navigation.
Batch comparisons stay accessible while child runs complete, and group downloads remain
available even if there is no completed child. Frozen HTML retains the same tabs and saved
observations; Design shows read-only reference/settings context. Duplex printing hides every
workspace panel and prints only the generated duplex sheet, regardless of the active tab.


## Screening order and design context

The detailed screening table defaults to original design-rank order and shows the same
saved rank and preferences-met/evaluated score as Candidates. An explicit sort control
also offers fewest reported matches, exact alignments and other displayed counts. Count
sorts start ascending; headings and a direction control allow reversing the order. Ties
retain original rank. Missing counts and incomplete/capped evidence do not sort as clean
zeroes. Zero-match candidates remain individual rows in the selected order, with an
explicit statement that absence of reported matches does not establish specificity.
