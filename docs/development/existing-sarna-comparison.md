# Compare existing saRNAs

## Published database observations

Gene lookup now includes **Previously reported sequences & experimental results**, using
the bundled [saRNA Database v1.0.0](https://zenodo.org/records/21871232) by Sušjan-Leite
and colleagues (CC BY 4.0). Human observations are selected by Ensembl gene ID rather
than guessed gene aliases. Search sequences, reported effects, cell models, conditions
and source citations; observations remain grouped by the database's sequence identity.

The gene overview marks database-reported GRCh38 sites on a separate literature lane.
These coordinates have not been independently verified, and their database TSS may differ
from the selected reference. Clicking a marker opens its observations and original sources.
Sites on other assemblies remain visible in the observation details rather than being
plotted on GRCh38.

New gene design runs also locate complete reported nucleotide strands in the loaded
promoter, allowing reverse complements and U/T normalization. A paired core can be matched
after removing terminal TT only for observations explicitly annotated with a dTdT overhang,
and only when the remaining sense/antisense cores are exact reverse complements. Unmarked
terminal bases are never guessed away. Plain strands shorter than 12 or longer than 100
bases and unsupported sequence notation remain browsable without automatic matching.

Exact matches appear on the promoter map and as **Already reported** candidate badges,
including excluded candidates. Badges describe sequence evidence, not proof of identical
experimental treatment or verified targeting role. The observation IDs supporting each
match are listed. The same sequence can have different results in different contexts.

Effects retain their original values and scale. Numeric values are not automatically
classified as activation or failure; qualitative values, controls, unknowns and conflicting
observations remain visible. The release does not consistently provide structured RNA,
protein, functional-rescue, statistical-support or replication fields, so those interpretations
require review of the original source. A source without a resolvable identifier is explicitly
labelled; papers link through DOI/PMID and patents through their publication number.

Wang remains the default. Literature evidence does not change gates, preference points or
ranking. Users can continue to adjust criteria deliberately through existing controls.
This first integration does not infer new sequence-pattern recommendations.

Literature observations, release identity and CSV hashes are frozen with new runs and
preserved in JSON, offline HTML and dedicated Excel sheets. Historical runs are not enriched
on reopening: make a new run to obtain the bundled database evidence. The database requires
no runtime download. Absence of a record or match does not establish that a gene or sequence
has never been tested.

## Researcher imports and existing reference duplexes

After completing a design run, open **Existing saRNAs & reference duplexes** in
Candidates. Add one sequence per line, or `name,sequence`, select its strand role
and direction, and preview before adding. Repeat to add more records. Uploads
accept UTF-8 TXT/CSV and Excel `.xlsx` (first worksheet), with at most 200 records
per run and 2 MB per file. Download the CSV template in the app.

Supported header fields are `name`, `sequence`, `role`, `direction`, `sense`,
`antisense`, and `source`. Roles are `unknown`, `target`, `sense`, `antisense`;
direction is `5to3` (default) or `3to5`. Supplied sense and antisense columns
represent the two duplex strands in the selected direction. Cores must be exact
reverse complements; mismatched duplex models are not supported. File imports
show row-level errors and preview valid records before saving.

U is replaced with T for same-strand DNA matching. Complementation is a separate
operation that reverses direction and pairs A with T, C with G. Unknown strands
are searched in both orientations against the run's transcript-oriented promoter.
Every exact occurrence is listed. A missing match means absent from the loaded
window, not absent from the genome. Mismatch-tolerant matching is not implemented.

Write known terminal overhangs explicitly as `[dT][dT]` or `[UU]`; these are
excluded from core matching and retained as metadata. Unmarked terminal bases
remain part of the core. Other chemical modifications require a separate,
unmodified core for comparison. A missing partner is inferred; missing overhang
chemistry remains unknown. Predictions describe paired cores, not modifications.

SCN1A automatically displays SCN1A-10 and SCN1A-7; CDKN1A/p21 displays dsP21-322.
SCN1A cores and dTdT annotations come from repository records. SCN1A-7 additionally
records researcher-reported downregulation in HeLa and AF22 by qPCR, with RNA-seq
in AF22 NES. This contextual evidence does not replace the benchmark's existing
“not activating” classification. References without a match remain listed.

Click a reference or import to zoom to its match and show both strands, per-model
guide/passenger predictions, and the matching candidate's current-run rule rating
(including exclusions). Select a newly generated candidate for comparison and
base overlap; nearby candidates can be opened from the reference details.
Different-length imports have model predictions but no current-run rule rating.

Imports are stored separately from scientific snapshots in SQLite and are
included in run/group JSON and offline HTML exports. Offline HTML supports
inspection; importing requires the running local Python app. Original run
configuration and candidate results are preserved.

Curated reference coordinates are also shown in the whole-gene overview, so a
reference remains locatable when the chosen promoter window excludes its site.
The table distinguishes **Outside current promoter window** using the recorded
assembly and chromosome. **Choose this promoter** selects the reference in
Design and preserves your current criteria and upstream/downstream spans.
Review the draft and submit when ready. Selecting a gene-map reference marker
also selects a draft without starting a run. The previous run remains saved. A known
coordinate is labelled as sequence-mapped repository evidence; it is not drawn
as an exact match inside a different promoter window.

Before running a design, gene lookup shows **Known saRNAs by promoter region**
in the Design gene overview. The table lists all available promoter references
whose current upstream/downstream window fully contains a curated site. Choose
a row to update the draft promoter without submitting a run. The promoter dropdown
also lists contained saRNA names, and the selection note reports when the chosen
window has no curated sites. Coverage updates as the window size changes; it is
coordinate-based guidance until exact matching is performed on the loaded sequence.

The curated literature positive controls include VEGF (VEGFA), P53 (TP53),
dsP21-322 (CDKN1A), VHL-821 (VHL), CEBPA and KLF4, using the thesis dataset's
recorded cores, citations and sequence-mapped coordinates. Historical labels
VEGF/P53/P21 resolve to the corresponding gene symbols. Unrecorded overhangs
remain unknown.

Design completion keeps the current workspace tab and scroll position. A saved
results notice in Design offers **Review candidates when ready**; the Candidates
tab is also available. Completing a promoter batch does not switch tabs either.
Offline reports initially open Candidates for inspection.

Landscape reference labels use separate callouts connected to the genomic core
bars. Overlapping callouts are placed on additional rows rather than drawn over
one another. Purple “↑ Activating” marks the curated positive references; rose
“↓ Downregulating” marks SCN1A-7 with the researcher-reported HeLa/AF22 context.
Imported sequences have teal “Imported” labels with no assigned effect. The same
colors and text labels appear in the reference table and inspection panel;
activity labels describe the recorded experimental setting, not predictions.

Generated candidates are labelled **Same core as [control]** when their complete
paired core equals a curated reference core (including an exchanged sense/
antisense orientation). Overlap, a shifted window, or a different length does not
qualify. These badges appear in the candidate workbench, inspection, top duplex,
and preference heatmap; matching landscape bars get the reference color, and
reference callouts name the corresponding generated candidate. Core identity
is not a claim of identical overhang chemistry or experimental efficacy.

**Analyse both strands** expands each reference's nucleotide composition,
numbered RNA bases, GC percentage, predicted guide/passenger sequences by model,
and antiparallel complementary DNA sequence. Genomic +/− is assigned from the
promoter orientation, or the recorded reference locus when outside the window.
Sense RNA complements the opposite genomic strand; antisense RNA complements
the gene-oriented genomic strand. Without a genomic orientation, the strand
label stays unknown. A single supplied strand has an explicitly inferred partner.
Sequence complementarity does not establish physical binding, substrate identity
(DNA versus promoter-associated RNA), or measured strand loading. Sources linked
in the panel include https://pmc.ncbi.nlm.nih.gov/articles/PMC4797303/ and
https://pmc.ncbi.nlm.nih.gov/articles/PMC4783471/.

Existing duplexes now show the candidate-style antiparallel ladder, guide and
passenger sequence cards, and end-stability evidence automatically. The selected
run model is used where assessed, with disagreement displayed. This reference
view shows RNA cores only so that unrecorded overhangs are never invented.
Detailed base composition and per-model tables remain available below it.

Matching control names and effects are shown beside candidate IDs in the top
saRNA/off-target table without inspection. The workbench opens when controls
are present and has an always-visible summary of exact matches, including those
excluded by current gates; rank, passing filters, export selection and off-target
findings remain independent of control identity. Control names can be searched
in the workbench.


Comparison analysis is frozen when calculated. New runs retain control results
inside the original snapshot; imported analyses are saved separately with model
assessments, engine version, software fingerprint, calculation timestamp, and
configuration fingerprint. Reopening a run and exporting it do not recalculate.
**Recalculate comparisons with current software** explicitly creates a new
revision; previous comparison revisions remain stored. Legacy analyses are
preserved and labelled when their original software provenance is unavailable.
Historical runs without frozen controls can add them through explicit recalculation.

Sequence writes carry only supplied input fields and a revision number, rather
than nested candidate results. A stale window receives a conflict and can
**Reload saved sequences** before retrying; pending input remains visible.
Adding a new input retains existing analyses under their original model versions.
The HTTP import route accepts the base64 envelope for a 2 MB file without
increasing the transport limit on unrelated routes. Excel worksheet order is
resolved using workbook relationships, including renamed or deleted sheets.
Map panning, zooming and resizing update only the map and preserve reference
panel disclosures and focus.

Browser regressions are available with:

```bash
.venv/bin/python scripts/verify_workbench_browser.py --chrome /path/to/chrome
```

They exercise explicit review navigation, real file uploads, repeated sequence
additions, saved comparison reads, stale-window recovery, persistent strand
panels, screening, exports, and offline reports without public screening calls.
