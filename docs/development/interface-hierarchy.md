# Interface hierarchy

Version 0.3.1 reorganises the workbench around the current task. The four workspace tabs
remain Design, Candidates, Off-target screening and Export. Scientific calculations,
ranking and screening policies are unchanged.

## Design

Target and promoter reference inputs come first, followed by the starting configuration
and primary design action. Experimental context, expression evidence, promoter comparisons,
custom coordinates, parameter overrides and source options open on demand. The gene map
stays visible; alternative start sites, literature observations and expression details
have separate disclosures. The top bar has one History entry, with quick-open and searchable
history inside it.

## Candidates

The initial view shows a short result summary, the promoter map and an immediately usable
candidate list with ten rows per page. The initial candidate's data can be preloaded without
opening its inspector. Clicking a candidate opens sequence information and concise notices;
paired ordering diagrams, model comparisons, individual rules with sources, and biological
context are separate inspector sections. Closing inspection returns to the list/map.

Map tools (Evaluate saRNA, Find sequence, Regions, Motifs and Measure distance) open
individually without clearing their inputs when closed. Existing reference-site lanes are
shown when evaluating a sequence; literature lanes appear when the reported-evidence panel
is open. Exclusion explanations and map interpretation are available below the map.

Saved shortlist notes, checklist/screening comparisons, literature observations and ranking
validation remain below the primary candidate workflow. Validation limits are summarised
next to the results, with the complete measurement under Ranking, rules and evidence.

## Screening

A compact reference-status row remains accessible before a run exists. Download size and
required temporary disk space are shown before installation; file breakdown, installation
methods and source versions are in Reference details. Installation progress and pause/retry
controls remain accessible in the compact row.

With a run, one primary action screens the active candidates against transcripts and genome.
Per-scope actions, engine selection and batch controls are under Search options. Choosing
NCBI reveals a data-sharing notice beside the primary submission action. Coverage, running
jobs, failed jobs and incomplete evidence stay visible. Detailed alignments, comparison
plots and methods open separately.

## Export

The active selection is declared before the three principal choices: Order sequences,
Share a summary and Analyse results. Additional formats live under More downloads. A separate,
initially closed sequence/evidence preview follows the choices. Opening it preserves its
state and keyboard focus when screening updates change the evidence. Whole-run formats
continue to include every candidate; sequence exports follow the active selection.

## Selection and state

Candidate checkboxes form the shared active selection for screening and sequence exports.
When no rows are checked, all design-passing candidates are the explicit default, except
for an intentionally empty managed screening selection. The header, screening action and
export summary name this scope. Screening still omits candidates rejected by design gates;
explicit ordering selections can include them, with their failed gates retained.

The saved shortlist is a separate collection of research notes. “Use saved shortlist” copies
its candidate IDs into the active selection; saving a note does not silently change downloads.
Shortlist saves and reloads retain map-investigation metadata.

Routine saved-state messages are quiet. Missing checks, failed gates, incomplete screening,
custom configuration, failed saves and stale-draft notices remain visible where relevant.
Rule sources and methodological explanations are attached to the corresponding disclosure.
Background screening updates preserve opened evidence sections and identifiable focused
controls; switching workspaces retains draft fields, selection and tool input.

## Verification

The unit suite checks workspace ownership, default disclosure states, download ordering,
reference setup availability and scientific invariants. The Chromium script checks the
initial list, optional tools, inspector state, shared screening action, export preview,
keyboard focus, offline reports and mobile layouts. Its screenshots use isolated fixtures;
they are UI verification artifacts, not biological screening results.
