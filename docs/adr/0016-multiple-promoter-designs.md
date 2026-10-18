# ADR 0016: Select multiple promoter references and preserve separate design contexts

Date: 2026-10-07. Status: accepted.

The workbench supports adjustable design for human genes with multiple promoters.

## Context and decision

A gene can have several annotated transcript start sites and additional literature-defined
promoter references. SCN1A already demonstrates why a unique MANE Select transcript is not
enough: its thesis designs lie near references other than MANE. A transcript's annotated start
site does not establish that its promoter is active in the researcher's experimental cell type.

The local workbench exposes a whole-gene reference overview and allows the researcher to select
one or several references. Transcript records with exactly the same contig, strand and start
coordinate are grouped for selection, while member transcript identities remain available.
The grouping uses exact coordinates; it does not infer a shared biological promoter from
proximity and supplies no cell-specific activity evidence. Declared literature anchors remain
distinct and retain their source and reference uncertainty.

Unique MANE Select remains the annotated default for an initial single-reference design.
When it is absent or ambiguous, the app requires an explicit choice. Selecting several
references is a deliberate research action, not an automatic search for the most active TSS.

Each selected reference produces an ordinary Python design run. Shared settings provide the
starting configuration, with supported per-reference adjustments. The group records its
member requests and saved run identities in local SQLite. Each child retains its own sequence,
coordinates, anchor provenance, effective configuration, complete verdicts and ranks. Group
views compare these contexts without pooling candidate ranks or presenting a single cross-
promoter activation score. A failure for one reference remains visible and does not erase the
successful runs for the others.

This revises ADR 0012's earlier rejection of designing against several TSSs: separate ranked
contexts can be inspected without the promoter-activity term that a pooled ranking would
require. The MANE policy, coordinate invariants and evidence limitations still apply.

## Researcher-supplied reference

The researcher may enter a custom **GRCh38 genomic position in 1-based coordinates** with a
name, source and explanatory note. The gene lookup supplies contig and transcription strand.
The Python adapter converts the supplied position to a 0-based reference coordinate once,
then uses the existing anchor/window and promoter-orientation methods. This input is distinct
from a pasted sequence's TSS-relative first-base position, which uses −1/+1 and has no zero.
The reference position and retrieval window must fit the primary chromosome. A position can
lie outside the annotated gene span, allowing an alternative promoter reference, but its
association with that gene remains the researcher's claim rather than a verified finding.

A custom reference is explicitly researcher-supplied and inferred. The supplied source is
recorded, not independently validated as a mapped TSS. Offsets and targeting-window verdicts
retain that uncertainty. The app must not change the reference label to measured merely
because a citation or cell name was entered.

## Settings and interpretation

Wang remains the default, including its window gate. Supported values and gate/preference/off
roles remain adjustable with their rationale, sources and consequences visible. Equal-weight
preference counting remains unchanged. Generation limits and chemistry support remain engine
constraints rather than unrestricted controls; the workbench's default target cores remain
19 nt. A core-sequence CSV does not add dTdT chemistry or off-target checks.

Every candidate stays tied to its reference and configuration. Ratings count evaluated
preferences, not activation probability, and cannot be treated as comparable probabilities
across promoters or changed rule sets. The named thermodynamic guide prediction remains
separate from genomic strand and does not establish the strand or identity of a target RNA.

Reopening a saved group retrieves its records without recalculating. Frozen exports preserve
the separate contexts and explanations; running a changed design continues to require the
local Python engine. Candidate-relative positions alone are not valid identifiers for matching
results from different anchors. The same genomic duplex may appear in more than one child run
with different TSS-relative offsets, verdicts and ranks.

## Costs and remaining limits

Multiple references increase sequence retrieval, computation and storage roughly with the
number selected. Visible group progress and per-reference failures are required. Reusing cached
gene records does not imply that every chosen promoter window is cached. Offline mode remains
explicit and never silently permits remote requests.

Human GRCh38 primary contigs remain the supported genomic scope. Human Ensembl gene IDs are
accepted, including an optional version suffix. The suffix is removed for lookup of the current
assembly record; entering it does not pin an annotation release. Stable gene-ID input reduces
dependence on a displayed gene symbol; this does not promise arbitrary synonyms or ambiguous
alias resolution. Complete available gene-reference overviews are frozen with completed gene
runs, preserving the options returned at that time. The interface is a reference overview rather
than an exon-aware genome browser.
It does not establish exhaustive biological promoter coverage.

EPDnew/FANTOM5 measured promoter adapters, expression and occupancy sources, generalized local
cell tracks, order-ready chemistry, off-target reporting, methods documents and installers
remain continuation work. This slice adds deliberate reference selection and provenance rather
than evidence of promoter activity or saRNA efficacy.

See [the workbench guide](../development/local-workbench.md) for usage and measured verification.
