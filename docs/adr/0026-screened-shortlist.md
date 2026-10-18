# 0026: A separate shortlist after complete local screening

> Automatic exclusions and export restrictions are superseded by [ADR 0027](0027-advisory-offtarget-results.md).

The user requests a top 10 or 15 that considers off-target findings and excludes close matches.
This supersedes the earlier decision to make screening display-only **for automatic selection**;
the frozen design scores, original ranking and manual research notes remain unchanged.

The Python policy `conservative-full-length-v1` excludes a candidate when either scope contains
any full-length 18/19 or 19/19 alignment, or at least five 17/19 alignments. These cutoffs are
**project choices without experimental validation as saRNA safety thresholds**. They are not
Wang's 80% control-duplex rule and do not estimate binding, efficacy or probability of harm.
Reference alignments include isoforms, alternate and patch records; the count is deliberately
conservative and must not be described as a count of independent genes or biological targets.
Transcript matches to the intended gene remain relevant to potential RNA silencing.

Both RNA and genomic searches must be complete on the same installed local reference. The
intended genomic locus must be known for coordinate-based own-hit removal. Missing evidence,
old reports without selection assessments, unknown own-locus coordinates, and capped reports
cannot qualify. A new result with missing evidence supersedes older evidence for that candidate
and scope. Refreshing older local reports reuses native alignment caches.

Complete per-identity counts are computed from all alignments before detail retention. The
200-detail display cap cannot change selection. Each local row carries the policy version,
thresholds, evidence label, outcome and reasons. The UI filters using these Python outcomes,
then retains original design rank: no arbitrary efficacy/specificity weighted score is added.
Top 10/15 means *up to* that many; it never fills a vacancy with an excluded or pending candidate.
The order remains provisional while other passing candidates lack screening decisions.

Choosing the screened list replaces the temporary export selection explicitly. Subsequent
failures remove managed candidates from that selection. Once this selection is in use,
sequence/selected spreadsheet exports reject empty, pending or excluded selections instead of
falling back to exporting all designs. Whole-run JSON/HTML remain accessible as scientific
records. A separate selection JSON records every decision, original rank, score, source search
ID, reference release, timestamp, identity counts and policy. Manual notes remain separate.

The screen links to Wang et al. (2015, PMC4379447), Meng et al. (2016, PMC4797303), and
Jackson et al. (2006, PMC1484447). The first two support RNAa sequence-design and recognition
context; the latter demonstrates siRNA seed-mediated off-target silencing. They do **not**
validate this project's numeric selection cutoffs. Seed effects, gaps, bulges, variants and
cell-specific activity remain outside the current local screening policy.

The reference setup status API also supplies a per-asset byte plan from pinned download
metadata. First-time setup opens automatically and shows download sizes before any download
is requested. Downloads, verification/index building and reuse are explained separately.
The installed view states no download is needed; normal application updates preserve references.
