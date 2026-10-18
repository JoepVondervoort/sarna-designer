# ADR 0015: Local Python jobs and one UI for live and frozen workbench results

Date: 2026-10-07. Status: accepted.

## Context and decision

Serving `out/` with `http.server` cannot retrieve a new gene or apply browser edits. The user
requires gene-to-promoter design, adjustable supported rules, visible ratings/guide predictions
and an explained interactive promoter map.

`python -m sarna.app` serves the workbench on loopback. A stdlib HTTP adapter dispatches validated
requests into the existing engine, using two background workers and local SQLite snapshots.
HTTPX is the sole additional workbench dependency. JavaScript controls drafts and displays
returned results; scientific evaluation and ranking remain Python. This adapter can be replaced
without changing the engine.

Any symbol accepted by the human GRCh38 Ensembl adapter can be retrieved. The user can choose
returned transcripts or declared literature anchors, or paste DNA. Only unique MANE Select is
chosen automatically. Missing mandatory source responses fail explicitly; optional annotation
failure does not alter ranks. Cell selection does not imply an implemented expression/occupancy
adapter or arbitrary-cell histone reader.

Saved snapshots retain input sequence/coordinates, effective configuration, preset baseline,
named-option identity, available rule declarations, full verdicts and engine version. Reopening
uses the record rather than current defaults. A Python rerun compares candidate sequence and
location only when promoter sequence/reference are unchanged.

The live workbench and its offline HTML export use **one UI implementation** with inline
CSS/JS/data. Setup is disabled in archives. Complete rejected-candidate verdicts make HTML large,
but preserve every explanation offline; SQLite compresses the underlying record. The older
`sarna.view` CLI viewer remains for existing scripts/historical reports, not as another settings
editor. Migrating its output later requires a separate compatibility review.

CSV is explicitly a **core-sequence export**, not order-ready chemistry. *Superseded in part by
[ADR 0017](0017-order-ready-duplex-export.md), which adds a separate `export.order.csv` carrying
the dTdT 21-mers; this file keeps its core-sequence meaning unchanged, but both exports are now
written in rank order rather than position order.* A named model's RNA
guide prediction is shown separately from other model predictions and genomic DNA strand.
This does not predict the identity of a promoter-associated target RNA. Ratings count criteria,
not activation probability.

## Correctness required by adjustment

`end-length-nt` now also controls side-by-side guide comparisons. The asymmetry gate supplies
the length, then the tiebreak if the gate is off, otherwise the existing four-base default.
Default calculations are unchanged. An independent two-base GC/AA calculation pins a
+2.49 kcal/mol RNA difference and the pipeline wiring.

Upstream window bounds previously allowed zero, which `TssPosition` rightly rejects. Their
maximum is now −1. The thesis-derived alternative is explicitly **the upstream portion** and
applies −2000…−1, rather than labelling an unrepresentable −2000…+500 range as implemented.
Wang's −1000…−100 default is unchanged.

## Costs and limits

The stdlib server reduces installation dependencies in the current environment, where dependency
downloads failed. It checks same origin/host, bounds JSON bodies and serves static files through
explicit routes. It is for a local single-user application, not public/multi-user hosting.
Windows/macOS installers and fresh installs remain unverified.

Ensembl annotation is not a measured-promoter source. Full upstream release/source-response
pinning, EPDnew, DepMap, FANTOM5, ENCODE, generalized local tracks, methods documents, order
chemistry and off-target reports remain continuation work. Coordinates, sequence and source
declarations are frozen today; comprehensive upstream dataset version provenance is incomplete.

See [the launch/testing guide](../development/local-workbench.md) for measured browser results.
