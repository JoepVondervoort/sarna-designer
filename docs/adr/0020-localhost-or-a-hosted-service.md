# ADR 0020: Localhost, or a hosted service?

Date: 2026-10-07. Status: accepted. Revisits the delivery decision in
[ADR 0014](0014-local-delivery-and-wang-default.md) now that colleagues are about to use the tool.

## Context

The workbench runs as a local application: `sarna-designer` starts a Python HTTP server bound to
`127.0.0.1`, opens a browser at it, and keeps everything in SQLite under the user's own directory.
Each colleague installs it ([INSTALL.md](../../INSTALL.md)) and runs their own copy.

The obvious question, asked directly: is that the best way? A hosted web app would mean no
installation, one shared cache, shared run history and one place to update. Those are real
advantages and this ADR does not pretend otherwise.

## Decision

**Stay local for now.** Three reasons, none of them inertia.

**1. The gene symbol is the research question.** Running locally means the only parties who learn
which gene is being worked on are Ensembl, and NCBI if the researcher opts into a screening, and
the interface says so beside that button. A hosted service adds its operator, anyone who can read
its database, and its backups. For unpublished work on a disease gene, that is a disclosure
decision, not an implementation detail, and it should not be made silently by a convenience.

**2. Centralising inverts the usual API argument.** Normally one shared server means fewer calls
to an upstream service, which is better for everyone. Here, [ADR 0019](0019-off-target-screening.md)
records that NCBI rate-limits and can block an **IP address**, and that a block would hit the whole
institution rather than one researcher. A hosted service funnels every screening through one
address, concentrating exactly the risk the current design spreads. Ensembl is different, it is
happy to be cached, and the shared-cache option below addresses that without a server.

**3. What runs today is not a server, and hosting it would be a rewrite.** `http.server` is
documented as not for production. There is no authentication, no TLS, no per-user separation; the
single SQLite file and the two-worker pool assume one user; the request handler's defence is that
it only accepts `Host: localhost`. Deploying it as-is would be irresponsible, and deploying it
properly means an ASGI application, a session store, an identity provider and somebody accountable
for uptime. The `fastapi`/`uvicorn` extra was deliberately removed in
[ADR 0018](0018-apache-2-and-bundled-data-licensing.md) because nothing used it; re-adding it is a
project, not a flag.

## What the local path costs, and what is being done about it

**Every colleague fetches the same reference data.** Five people designing against SCN1A fetch the
same gene, transcripts and sequence five times. `SARNA_CACHE_DIR` already points the response cache
anywhere, including a shared network folder, so a lab can warm one cache and read it. That was not
safe before today: entries were written with a plain `write_text`, so a reader could see a
half-written file. Writes are now atomic, which is the precondition for sharing a directory.

**Runs are not shareable.** A colleague cannot send a link to a result. They can send the offline
HTML export, which is self-contained and carries the full record, that is the honest substitute,
and it has the advantage that the recipient gets a frozen artefact rather than a live page whose
answer may change.

**Updating means re-installing.** `uv tool upgrade sarna-designer` is the whole procedure, and
`sarna-designer-check` verifies the result. This is a documentation problem, not an architecture
one.

**Someone has to install Python.** `uv` fetches its own interpreter, so the install is two
commands and does not touch the system Python. Fresh Windows and macOS runs remain unverified by a
human, which is the real open risk in the local path and is not fixed by hosting. It is fixed by
somebody running `sarna-designer-check` on those machines.

## What would change this decision

Any one of these, and the answer flips:

- Colleagues need to **share runs** with each other rather than files.
- More than a handful of people use it, so per-machine installation becomes the bottleneck.
- Someone who cannot install software needs access.
- An institutional host with SSO becomes available, which removes the authentication objection and
  most of the operational one.

If that happens, the path is not to expose the current server. It is an ASGI application behind
institutional identity, with per-user run storage, and with screening still submitted under
credentials that make the NCBI rate limit attributable to a person rather than to a building. The
engine itself (`design_sarnas`, the rules, the exports) is already a library with no web
framework in it ([ADR 0005](0005-python-core-with-thin-web-ui.md)), so that rewrite touches `sarna.app` and
nothing below it. That is the property worth protecting, and it is protected.

## A middle option, considered and not taken

Packaging a double-clickable desktop application (PyInstaller, Briefcase) removes the install step
without a server and keeps every privacy property. It was not taken now because it needs a build
and a code-signing identity for each of Windows and macOS, and neither platform's install path has
yet been verified by a human even in its simple form. Verifying the simple path first is the
cheaper order.

## Costs of this decision

Colleagues each carry an installation, each warm their own cache unless the shared directory is
set up, and cannot hand each other a URL. Nobody is accountable for the tool being up, which is
also the point. And the decision is being made on reasoning rather than on observed multi-user
demand, because there is not yet any multi-user use to observe.
