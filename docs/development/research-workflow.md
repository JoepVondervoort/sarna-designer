# Research workflow

The workbench keeps design calculations in Python and researcher decisions in a separate local
workspace. A shortlist is a selection for review, not a model prediction or an experimental result.

## Navigate the work

The stage bar links Target, Reference & context, Candidates, Screening and Export. Candidates and
later stages become available after a design result exists. Design setup stays visible alongside
results on desktop and above them on small screens. **Go to design setup** or the first two stages
navigate back to its controls without hiding the panel. Settings remain drafts until submitted
as a new design; saved scientific results are
not changed by workspace navigation.

On a wide desktop, candidate explanations sit beside the table and remain visible while moving
between candidates. On narrower screens the inspector follows the table in one column. The same
sequence, coordinate, guide-model and rule explanations remain available.

## Build a shortlist with reasons

Inspect a candidate and choose **Add to research shortlist**. Record why you selected it, any
missing evidence, or a question to investigate. Notes save automatically to the local SQLite store.
The status reports saved, saving, or failure; it never treats a failed save as successful.

A shortlist belongs to one run. Candidate identifiers from another run are rejected by the API.
Selections do not automatically transfer to a rerun, even when candidate names look the same.
An excluded candidate can be retained for review, with its failed-gate status visible.

Table checkboxes remain a temporary export/screening selection. **Select shortlist for export /
screening** deliberately copies the shortlist into that selection. Screening searches passing
candidates only. Explicit order selections can include excluded candidates under the existing
export behavior; shortlisting does not remove or override their failed-gate reports.

Notes are stored separately from the frozen scientific snapshot and reproducibility fingerprints.
Saving them cannot change rules, candidate ranks, sequence identities or core replay results.
JSON and offline HTML attach a `workspace.shortlist` record. Excel adds a Shortlist sheet with all
saved selections and reasons for the run, independently of any candidate-row export filter.

Exports from the app wait for pending shortlist saves. If notes cannot be saved, resolve the error
before exporting the workspace or use **Download shortlist notes** to preserve local edits in a
separate JSON file. Offline HTML carries read-only archived notes, while its temporary candidate
selection remains usable.

Each save includes the last known revision. A stale browser receives a conflict rather than
replacing another window's work. Its local notes remain in that window. Download them if needed
before choosing **Reload saved notes**, which replaces local edits with the saved version.

## Find and compare saved runs

**Run history** searches gene/label, reference metadata or run ID, with a status filter and pages
of 50 runs. Searching operates in SQLite and can reach records older than the first 100 results.
Recent promoter groups remain available from the same dialog. The desktop quick-open dropdown
also remains available.

Select two completed runs and choose **Compare two selected runs**. Selection order defines the
first and second columns. Comparison shows reference, assembly, cell context, configuration,
software identity and passing counts, followed by changes to effective rule roles, values and
named options.

Candidate admission/rank differences are counted only when the saved promoter and reference match.
Different input or reference produces an explicit independent-context comparison. Neither higher
passing counts nor changed ranks establish improved biological activity. A matching cell name
also does not establish equivalent cell-specific evidence; inspect each run's source records.

Comparison reads existing snapshots; it does not rerun the engine or update either result.
