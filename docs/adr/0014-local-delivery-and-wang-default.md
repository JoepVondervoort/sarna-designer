# ADR 0014: Local delivery with Wang as the default

- **Status**: Accepted; local application and packaging pending
- **Date**: 2026-10-07
- **Deciders**: Joep Vondervoort

## Context and decision

The user selected local hosting and accessible installation so researchers can use the tool in
their own environments. They clarified that papers are possible recommendations, not universal
truth, and that Wang's paper remains the default.

Keep the current Wang 2015 default, including the mandatory TSS window and position-19 criterion.
This supersedes the earlier continuation proposal to make a window preference the default.
Researchers may change supported roles and thresholds in clearly labelled custom configurations
under [ADR 0013](0013-researchers-can-adjust-rule-roles.md). Explain published support, project
interpretations, known exceptions and custom departures without implying validated efficacy.
Other literature and curated datasets can support optional alternatives and external validation;
they do not replace the default solely because they are newer.

Keep Python as the sole scientific engine. Deliver a local browser workbench through a thin
loopback service and an accessible launcher. HTML/CSS/JavaScript can implement the interface;
scientific evaluation and ranking execute in Python. Use local SQLite for run snapshots, as
already accepted, and preserve frozen exports as independent records.

Provide tested installation paths for supported Windows, macOS and Linux platforms. Aim for
research-user packages without requiring Python programming or Docker. Establish a real support
matrix through a packaging feasibility check; optional native dependencies may require alternate
installation routes. Containers and developer Python installation remain additional paths.
Separate installation files from per-user runs/caches and document updates and backup/restore.

## Costs and limits

Keeping the Wang window gate preserves a published baseline but retains known exclusions:
the chapter 4 benchmark with recovered locations reports 10 of 14 measured activators excluded
by that window. Custom alternatives can examine this cost; fewer exclusions do not establish
better rank prediction. This range-restricted dataset cannot validate general efficacy.

Local delivery avoids maintaining a public service, but shifts responsibility to packaging,
platform testing and update support. Python does not by itself guarantee native dependencies
work on every machine. Optional source retrieval still sends documented gene/coordinate queries
outside the lab; cached/offline use must be explicit. Uploaded sequences and results stay local.

## Verification

Pin the unchanged default independently of custom-role tests. Verify effective roles consistently
across evaluation, gates, ranks and archived configuration, without modifying the catalogue.
Before each platform release, verify clean installation, launch, an offline example, custom
validation, save/export, restart/reopen and upgrade with preserved data.
