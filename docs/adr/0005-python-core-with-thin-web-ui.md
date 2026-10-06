# ADR 0005 — A pure-Python design engine as a library, with a thin web UI over it

- **Status**: Accepted
- **Date**: 2026-10-06
- **Deciders**: Joep Vondervoort

## Context

The stack was deliberately left undecided through Phase 0 so that it would follow the
requirements rather than precede them. Those requirements are now known, and three of them
constrain the choice directly:

- **The science must be auditable and regression-tested.** The `legacy-scanner` preset has to
  reproduce `legacy/sarna_promoter_scanner.py` bit-for-bit ([ADR 0002](0002-wang-2015-default-legacy-preset.md)),
  which means the engine must be callable and assertable without a UI, a server or a browser.
- **The target user is a colleague without programming experience.** A library alone does not
  meet the brief.
- **Ranking must be explainable per candidate, per rule** ([ADR 0003](0003-rule-provenance-is-first-class.md)).
  The UI is not a table of scores; it has to render provenance, evidence class and per-term
  contributions, which is real interface work.

The ecosystem also constrains it. The legacy code is Python; the bigWig reader (`pyBigWig`),
sequence handling and the coordinate liftover required by [ADR 0004](0004-assembly-grch38-canonical-liftover-at-adapter.md)
are all Python-first with no credible JavaScript equivalent. Re-implementing the science in
another language would mean re-deriving it, which is precisely the failure mode this project
exists to end.

## Decision

Three layers, with a hard dependency direction.

**1. `sarna` — the engine, a plain importable Python package.** Candidate generation, the rule
registry, filtering, scoring and explanation. It has **no knowledge of HTTP, of the UI, or of
any web framework**, and no global state. Rules are the data records ADR 0003 specifies.
Genomics and epigenomics access sits behind interfaces, so the engine is testable with no
network and no 890 MB of bigWigs present. This layer is where the regression suite points, and
it is the deliverable that outlives any UI choice.

**2. A FastAPI service** that adapts the engine to HTTP: request validation, error mapping,
caching of external calls, export generation. It contains **no design rules**. Pydantic models
at the boundary satisfy the project's validate-at-boundaries requirement, and the schema is
generated from them rather than maintained by hand.

**3. A thin frontend** for the researcher-facing workflow: enter a gene symbol or paste a
promoter, choose a preset, see ranked candidates with their rule breakdown, export.
"Thin" is a constraint, not a description — it holds no scientific logic and computes no
scores, so a wrong number can only ever come from the engine.

Streamlit was considered and rejected for the UI: it is markedly faster to a usable tool, but
the per-rule provenance and contribution rendering ADR 0003 mandates is exactly the kind of
layout it makes awkward, and that rendering is a core requirement rather than a nicety. A full
React application was considered and judged correct in shape but premature in cost — the
decision here does not preclude it, because the frontend is the only layer it would replace.

## Consequences

- The engine can be exercised end-to-end from pytest before any server or page exists, so
  Phase 2 can start on the science rather than on scaffolding.
- Because the engine is UI-agnostic, a CLI is nearly free and will be added for batch and
  reproducibility work.
- `pyproject.toml` gains real package metadata, replacing the tooling-only configuration that
  deliberately deferred this decision. The repository stays a **single package with a single
  virtualenv**, not a monorepo; there is one deployable and one test suite, and a monorepo would
  add structure for a second consumer that does not exist.
- Strict `mypy` applies to the engine from the first commit. The typed coordinates ADR 0004
  requires are only load-bearing if type checking is enforced.
- The frontend's thinness must be actively defended in review. The moment it computes a score,
  the audit guarantee this ADR exists to protect is gone.
- Python 3.11 is pinned to match the verified environment, in an isolated `.venv/` so the
  miniconda base environment is never mutated.
